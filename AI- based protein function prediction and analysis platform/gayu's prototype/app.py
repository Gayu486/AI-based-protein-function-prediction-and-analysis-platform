"""
Flask Application Entry Point
AI-BASED PROTEIN FUNCTION PREDICTION AND ANALYSIS PLATFORM

Provides REST endpoints for:
- Live protein sequence validation (/api/validate)
- Multi-tier biological analysis & ML prediction (/api/analyze)
- Report exports in HTML, JSON, and FASTA formats (/api/report/*)
- Curated reference sample retrieval (/api/samples)
"""

import os
import json
import logging
from datetime import datetime, timezone
from flask import Flask, render_template, request, jsonify, make_response, send_from_directory

from bio_analysis.validator import validate_protein_sequence
from bio_analysis.physicochemical import analyze_physicochemical_properties
from bio_analysis.uniprot_service import fetch_uniprot_by_accession, search_uniprot_by_sequence
from bio_analysis.predictor import predict_protein_function
from bio_analysis.report_generator import generate_html_report, generate_json_report, generate_fasta_export

# Configure logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(name)s: %(message)s")
logger = logging.getLogger(__name__)

app = Flask(__name__, static_folder="static", template_folder="templates")
app.config['MAX_CONTENT_LENGTH'] = 2 * 1024 * 1024  # 2 MB max upload limit


@app.route('/')
def index():
    """Serves the main application dashboard."""
    return render_template('index.html')


@app.route('/health', methods=['GET'])
def health():
    """Service health check endpoint."""
    return jsonify({
        'status': 'healthy',
        'service': 'AI-Based Protein Function Prediction Platform',
        'timestamp': datetime.now(timezone.utc).isoformat()
    })


@app.route('/api/samples', methods=['GET'])
def get_sample_proteins():
    """Returns curated preset sample proteins."""
    samples_path = os.path.join(app.static_folder, 'data', 'sample_proteins.json')
    try:
        if os.path.exists(samples_path):
            with open(samples_path, 'r', encoding='utf-8') as f:
                samples = json.load(f)
            return jsonify({'success': True, 'samples': samples})
    except Exception as e:
        logger.error(f"Error loading sample proteins: {e}")
    return jsonify({'success': False, 'samples': []}), 500


@app.route('/api/validate', methods=['POST'])
def validate_sequence_endpoint():
    """
    Lightweight endpoint for live sequence validation as user types or pastes.
    """
    data = request.get_json(silent=True) or {}
    sequence_input = data.get('sequence', '')

    validation_result = validate_protein_sequence(sequence_input)
    return jsonify(validation_result)


@app.route('/api/analyze', methods=['POST'])
def analyze_protein():
    """
    Primary analysis endpoint.
    Accepts sequence, accession, or FASTA file upload.
    Runs validation, physicochemical analysis, biological database lookup,
    and machine learning function prediction.
    """
    sequence_text = ""
    accession_query = ""
    header_name = ""

    # Check for multipart file upload
    if 'file' in request.files:
        file = request.files['file']
        if file and file.filename != '':
            try:
                content = file.read().decode('utf-8', errors='replace')
                sequence_text = content
            except Exception as e:
                return jsonify({'success': False, 'error': f"Failed to read uploaded file: {str(e)}"}), 400

    # If no file, check JSON body or form data
    if not sequence_text:
        if request.is_json:
            req_data = request.get_json(silent=True) or {}
            sequence_text = req_data.get('sequence', '').strip()
            accession_query = req_data.get('accession', '').strip().upper()
        else:
            sequence_text = request.form.get('sequence', '').strip()
            accession_query = request.form.get('accession', '').strip().upper()

    # If user provided a UniProt accession directly and no sequence
    if accession_query and not sequence_text:
        record = fetch_uniprot_by_accession(accession_query)
        if not record:
            return jsonify({
                'success': False,
                'error': f"The UniProt accession '{accession_query}' was not found in UniProtKB or the curated reference database. Please check the accession ID or submit a raw amino acid sequence."
            }), 404

        sequence_text = record.get('sequence', '')
        header_name = f"UniProtKB:{record.get('accession')} {record.get('protein_name')}"

    # Validate sequence
    validation = validate_protein_sequence(sequence_text)
    if not validation['is_valid']:
        return jsonify({
            'success': False,
            'validation': validation,
            'error': "Sequence validation failed. " + " ".join(validation['errors'])
        }), 400

    clean_seq = validation['cleaned_sequence']
    if validation['header']:
        header_name = validation['header']

    try:
        # 1. Physicochemical analysis
        physicochemical = analyze_physicochemical_properties(clean_seq)

        # 2. Biological database identification & alignment
        if accession_query:
            # Direct accession was queried
            direct_record = fetch_uniprot_by_accession(accession_query)
            if direct_record:
                from bio_analysis.alignment import align_sequences
                aln = align_sequences(clean_seq, direct_record.get('sequence', clean_seq))
                identification = {
                    'matched': True,
                    'status': 'Identified Match by Accession',
                    'record': direct_record,
                    'alignment': aln,
                    'database_source': direct_record.get('data_source', 'UniProtKB')
                }
            else:
                identification = search_uniprot_by_sequence(clean_seq)
        else:
            identification = search_uniprot_by_sequence(clean_seq)

        # 3. Machine Learning Function Predictions
        predictions = predict_protein_function(clean_seq)

        # 4. Compile analysis result payload
        timestamp = datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M:%S UTC')
        result_payload = {
            'success': True,
            'timestamp': timestamp,
            'header': header_name,
            'sequence': clean_seq,
            'validation': validation,
            'physicochemical': physicochemical,
            'identification': identification,
            'predictions': predictions
        }

        return jsonify(result_payload)

    except Exception as e:
        logger.exception("Internal error during protein analysis:")
        return jsonify({
            'success': False,
            'error': f"An error occurred during biological analysis: {str(e)}. Please check the input sequence and retry."
        }), 500


@app.route('/api/report/html', methods=['POST'])
def export_html_report():
    """Generates and serves a downloadable standalone HTML report."""
    data = request.get_json(silent=True)
    if not data or not data.get('sequence'):
        return make_response("Invalid or missing analysis payload for report generation.", 400)

    html_content = generate_html_report(data)
    acc = data.get('identification', {}).get('record', {}).get('accession', 'analysis') if data.get('identification') else 'analysis'
    filename = f"Protein_Analysis_Report_{acc}_{datetime.now(timezone.utc).strftime('%Y%m%d_%H%M%S')}.html"

    response = make_response(html_content)
    response.headers['Content-Type'] = 'text/html; charset=utf-8'
    response.headers['Content-Disposition'] = f'attachment; filename="{filename}"'
    return response


@app.route('/api/report/json', methods=['POST'])
def export_json_report():
    """Generates and serves a downloadable JSON export of the analysis payload."""
    data = request.get_json(silent=True)
    if not data or not data.get('sequence'):
        return make_response("Invalid or missing analysis payload for JSON export.", 400)

    json_str = generate_json_report(data)
    acc = data.get('identification', {}).get('record', {}).get('accession', 'analysis') if data.get('identification') else 'analysis'
    filename = f"Protein_Analysis_Data_{acc}_{datetime.now(timezone.utc).strftime('%Y%m%d_%H%M%S')}.json"

    response = make_response(json_str)
    response.headers['Content-Type'] = 'application/json; charset=utf-8'
    response.headers['Content-Disposition'] = f'attachment; filename="{filename}"'
    return response


@app.route('/api/report/fasta', methods=['POST'])
def export_fasta_file():
    """Generates and serves a downloadable FASTA file."""
    data = request.get_json(silent=True)
    if not data or not data.get('sequence'):
        return make_response("Invalid or missing analysis payload for FASTA export.", 400)

    fasta_str = generate_fasta_export(data)
    acc = data.get('identification', {}).get('record', {}).get('accession', 'sequence') if data.get('identification') else 'sequence'
    filename = f"{acc}.fasta"

    response = make_response(fasta_str)
    response.headers['Content-Type'] = 'text/plain; charset=utf-8'
    response.headers['Content-Disposition'] = f'attachment; filename="{filename}"'
    return response


if __name__ == '__main__':
    port = int(os.environ.get('PORT', 5000))
    print(f"\n=======================================================")
    print(f" AI-BASED PROTEIN FUNCTION PREDICTION AND ANALYSIS PLATFORM")
    print(f" Academic Bioinformatics & Machine Learning Prototype")
    print(f" Server running at: http://127.0.0.1:{port}")
    print(f"=======================================================\n")
    app.run(host='0.0.0.0', port=port, debug=True)
