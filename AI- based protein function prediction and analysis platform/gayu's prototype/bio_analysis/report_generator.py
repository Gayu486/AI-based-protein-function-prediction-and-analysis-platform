"""
Report Generator Module
Generates comprehensive, publication-ready analysis reports in HTML, JSON, and FASTA formats.
Includes:
- Analysis metadata, timestamps, and sequence metrics
- Full physicochemical properties profile
- Database identification & sequence similarity statistics (Identity %, Coverage %, E-value)
- Gene Ontology functional annotations & evidence grades
- Protein families, domains & functional sites
- Computational AI predictions (Enzyme class, Subcellular localization, TM helices)
- External database links and scientific limitation disclosures
"""

import json
from datetime import datetime, timezone
from typing import Dict, Any


def generate_json_report(analysis_result: Dict[str, Any]) -> str:
    """Serializes analysis results into formatted JSON string."""
    return json.dumps(analysis_result, indent=2, default=str)


def generate_fasta_export(analysis_result: Dict[str, Any]) -> str:
    """Generates a standard FASTA file with descriptive header."""
    meta = analysis_result.get('identification', {})
    rec = meta.get('record', {}) if meta else {}
    name = rec.get('protein_name', 'Analyzed_Protein') if rec else 'Analyzed_Protein'
    acc = rec.get('accession', 'UNKNOWN') if rec else 'UNKNOWN'
    seq = analysis_result.get('sequence', '')
    length = len(seq)

    header = f">sp|{acc}|{name} Length={length} AnalyzedBy=AI-Based-Protein-Platform"
    
    # Format sequence in standard 60-character lines
    lines = [header]
    for i in range(0, length, 60):
        lines.append(seq[i:i + 60])
        
    return "\n".join(lines)


def generate_html_report(analysis_result: Dict[str, Any]) -> str:
    """
    Renders a standalone, self-contained, publication-grade HTML report
    with embedded clean CSS for both screen viewing and printing (PDF export).
    """
    timestamp = analysis_result.get('timestamp', datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M:%S UTC'))
    seq = analysis_result.get('sequence', '')
    length = len(seq)
    phys = analysis_result.get('physicochemical', {})
    ident = analysis_result.get('identification', {})
    record = ident.get('record') if ident else None
    alignment = ident.get('alignment') if ident else None
    pred = analysis_result.get('predictions', {})
    validation = analysis_result.get('validation', {})

    # HTML builder
    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>Protein Analysis Report - {record.get('protein_name', 'Sequence Analysis') if record else 'Sequence Analysis'}</title>
  <style>
    :root {{
      --primary: #1e3a8a;
      --primary-light: #3b82f6;
      --secondary: #0d9488;
      --dark: #0f172a;
      --gray-light: #f8fafc;
      --border: #e2e8f0;
      --text: #334155;
      --accent: #6366f1;
    }}
    * {{ box-sizing: border-box; margin: 0; padding: 0; }}
    body {{
      font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif;
      color: var(--text);
      background: #f1f5f9;
      line-height: 1.6;
      padding: 30px 15px;
    }}
    .report-card {{
      max-width: 960px;
      margin: 0 auto;
      background: #ffffff;
      border-radius: 12px;
      box-shadow: 0 10px 25px -5px rgba(0,0,0,0.08);
      overflow: hidden;
      border: 1px solid var(--border);
    }}
    .report-header {{
      background: linear-gradient(135deg, #1e293b 0%, #0f172a 100%);
      color: #ffffff;
      padding: 35px 40px;
      position: relative;
    }}
    .report-title {{
      font-size: 22px;
      font-weight: 700;
      letter-spacing: -0.5px;
      color: #38bdf8;
      margin-bottom: 5px;
      text-transform: uppercase;
    }}
    .report-subtitle {{
      font-size: 15px;
      color: #94a3b8;
    }}
    .report-meta-grid {{
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
      gap: 15px;
      margin-top: 25px;
      padding-top: 20px;
      border-top: 1px solid #334155;
      font-size: 13px;
    }}
    .meta-item strong {{ color: #e2e8f0; display: block; margin-bottom: 2px; }}
    .report-body {{
      padding: 40px;
    }}
    .section-title {{
      font-size: 17px;
      font-weight: 700;
      color: var(--dark);
      margin-top: 35px;
      margin-bottom: 15px;
      padding-bottom: 8px;
      border-bottom: 2px solid #e2e8f0;
      display: flex;
      align-items: center;
      gap: 10px;
    }}
    .section-title:first-child {{ margin-top: 0; }}
    .badge {{
      display: inline-block;
      padding: 3px 10px;
      border-radius: 20px;
      font-size: 11px;
      font-weight: 600;
      text-transform: uppercase;
      letter-spacing: 0.5px;
    }}
    .badge-reviewed {{ background: #dcfce7; color: #166534; border: 1px solid #bbf7d0; }}
    .badge-unreviewed {{ background: #fef3c7; color: #92400e; border: 1px solid #fde68a; }}
    .badge-source {{ background: #e0e7ff; color: #3730a3; border: 1px solid #c7d2fe; }}
    .badge-calc {{ background: #f3e8ff; color: #6b21a8; border: 1px solid #e9d5ff; }}
    .badge-pred {{ background: #fee2e2; color: #991b1b; border: 1px solid #fecaca; }}

    .grid-2 {{
      display: grid;
      grid-template-columns: 1fr 1fr;
      gap: 20px;
    }}
    .grid-3 {{
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(240px, 1fr));
      gap: 15px;
    }}
    .stat-box {{
      background: var(--gray-light);
      padding: 16px 20px;
      border-radius: 8px;
      border: 1px solid var(--border);
    }}
    .stat-label {{ font-size: 12px; font-weight: 600; text-transform: uppercase; color: #64748b; margin-bottom: 4px; }}
    .stat-val {{ font-size: 20px; font-weight: 700; color: var(--dark); }}
    .stat-desc {{ font-size: 12px; color: #64748b; margin-top: 4px; }}

    table {{
      width: 100%;
      border-collapse: collapse;
      font-size: 13px;
      margin: 15px 0;
    }}
    th, td {{
      padding: 10px 14px;
      text-align: left;
      border-bottom: 1px solid var(--border);
    }}
    th {{ background: #f8fafc; font-weight: 600; color: #475569; }}
    tr:hover td {{ background: #f8fafc; }}

    .sequence-block {{
      font-family: ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas, monospace;
      font-size: 12px;
      background: #0f172a;
      color: #38bdf8;
      padding: 20px;
      border-radius: 8px;
      overflow-x: auto;
      line-height: 1.8;
      letter-spacing: 1px;
      word-break: break-all;
    }}

    .alignment-box {{
      font-family: ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas, monospace;
      font-size: 11px;
      background: #f8fafc;
      border: 1px solid #cbd5e1;
      padding: 15px;
      border-radius: 6px;
      overflow-x: auto;
      white-space: pre;
      line-height: 1.4;
      color: #1e293b;
    }}

    .callout {{
      background: #eff6ff;
      border-left: 4px solid var(--primary-light);
      padding: 16px 20px;
      border-radius: 0 8px 8px 0;
      margin: 15px 0;
      font-size: 13px;
    }}
    .callout-warning {{
      background: #fffbeb;
      border-left-color: #f59e0b;
      color: #92400e;
    }}
    .btn-print {{
      display: inline-block;
      margin-top: 20px;
      padding: 10px 24px;
      background: #2563eb;
      color: white;
      text-decoration: none;
      border-radius: 6px;
      font-weight: 600;
      font-size: 13px;
      cursor: pointer;
      border: none;
    }}
    .btn-print:hover {{ background: #1d4ed8; }}

    @media print {{
      body {{ background: #ffffff; padding: 0; }}
      .report-card {{ box-shadow: none; border: none; }}
      .no-print {{ display: none !important; }}
    }}
  </style>
</head>
<body>
  <div class="report-card">
    <div class="report-header">
      <div class="report-title">AI-Based Protein Function Prediction & Analysis Platform</div>
      <div class="report-subtitle">Comprehensive Biological Sequence & Computational Characterization Report</div>
      <div class="report-meta-grid">
        <div class="meta-item">
          <strong>Analysis Timestamp</strong>
          <span>{timestamp}</span>
        </div>
        <div class="meta-item">
          <strong>Sequence Length</strong>
          <span>{length} Amino Acids</span>
        </div>
        <div class="meta-item">
          <strong>Identified Record</strong>
          <span>{record.get('protein_name', 'Novel / Unannotated Sequence') if record else 'Unknown / Query Sequence'}</span>
        </div>
        <div class="meta-item">
          <strong>Accession ID</strong>
          <span>{record.get('accession', 'N/A') if record else 'N/A'}</span>
        </div>
      </div>
    </div>

    <div class="report-body">
      <div class="no-print" style="margin-bottom: 25px; text-align: right;">
        <button class="btn-print" onclick="window.print()">Print or Save as PDF</button>
      </div>

      <!-- SECTION 1: PROTEIN IDENTIFICATION -->
      <h2 class="section-title">
        1. Biological Identification & Sequence Alignment
        <span class="badge badge-source">UniProtKB Integration</span>
      </h2>
"""

    if record:
        reviewed_badge = '<span class="badge badge-reviewed">Reviewed (Swiss-Prot)</span>' if record.get('reviewed') else '<span class="badge badge-unreviewed">Unreviewed (TrEMBL)</span>'
        html += f"""
      <table style="margin-bottom: 20px;">
        <tr><th style="width: 25%;">Protein Name</th><td><strong>{record.get('protein_name', 'N/A')}</strong> ({record.get('entry_name', 'N/A')}) {reviewed_badge}</td></tr>
        <tr><th>Gene Symbol</th><td>{record.get('gene_name', 'Not annotated')}</td></tr>
        <tr><th>Source Organism</th><td>{record.get('organism', 'N/A')} (NCBI Taxon ID: {record.get('taxon_id', 'N/A')})</td></tr>
        <tr><th>UniProt Accession</th><td><a href="https://www.uniprot.org/uniprotkb/{record.get('accession')}" target="_blank">{record.get('accession')}</a></td></tr>
        <tr><th>Database Source</th><td>{ident.get('database_source', 'UniProtKB')}</td></tr>
      </table>
"""
        if alignment:
            html += f"""
      <div class="grid-3" style="margin-bottom: 15px;">
        <div class="stat-box">
          <div class="stat-label">Sequence Identity</div>
          <div class="stat-val" style="color: #2563eb;">{alignment.get('identity_pct')}%</div>
          <div class="stat-desc">{alignment.get('identical_count')} / {alignment.get('alignment_length')} aligned residues</div>
        </div>
        <div class="stat-box">
          <div class="stat-label">Alignment Coverage</div>
          <div class="stat-val" style="color: #059669;">{alignment.get('coverage_pct')}%</div>
          <div class="stat-desc">Of query protein length</div>
        </div>
        <div class="stat-box">
          <div class="stat-label">Karlin-Altschul E-Value</div>
          <div class="stat-val" style="color: #7c3aed;">{alignment.get('e_value_str')}</div>
          <div class="stat-desc">Raw score: {alignment.get('alignment_score')}</div>
        </div>
      </div>
"""
            # Visual alignment block snippet
            blocks = alignment.get('blocks', [])
            if blocks:
                html += """<div class="stat-label" style="margin-top: 15px; margin-bottom: 5px;">Pairwise Sequence Alignment Snippet:</div><div class="alignment-box">"""
                for blk in blocks[:4]:
                    html += f"Query   {blk['query_start']:4d}  {blk['query_str']}  {blk['query_end']}\n"
                    html += f"              {blk['match_str']}\n"
                    html += f"Sbjct   {blk['subject_start']:4d}  {blk['subject_str']}  {blk['subject_end']}\n\n"
                html += "</div>"
    else:
        html += """
      <div class="callout callout-warning">
        <strong>No High-Identity Database Match:</strong> The submitted sequence did not yield a close homologue (> 25% identity) in Swiss-Prot. The characteristics and computational predictions below reflect sequence-derived properties.
      </div>
"""

    # SECTION 2: FUNCTIONAL ANNOTATIONS
    html += f"""
      <h2 class="section-title">
        2. Functional Annotations & Gene Ontology (GO)
        <span class="badge badge-source">Biological Evidence</span>
      </h2>
"""
    if record and record.get('function'):
        html += f"""
      <div class="callout">
        <strong>Annotated Biological Function:</strong><br>
        {record.get('function')}
      </div>
"""
        go_data = record.get('go_terms', {})
        html += """
      <div class="grid-3" style="margin-top: 15px;">
        <div class="stat-box">
          <div class="stat-label" style="color: #2563eb;">Molecular Function (GO)</div>
          <ul style="padding-left: 18px; font-size: 12px; margin-top: 8px;">
"""
        for item in go_data.get('molecular_function', [])[:5]:
            html += f"<li>{item['name']} <small style='color:#64748b;'>({item['evidence']})</small></li>"
        if not go_data.get('molecular_function'):
            html += "<li>None annotated</li>"
        html += """
          </ul>
        </div>
        <div class="stat-box">
          <div class="stat-label" style="color: #059669;">Biological Process (GO)</div>
          <ul style="padding-left: 18px; font-size: 12px; margin-top: 8px;">
"""
        for item in go_data.get('biological_process', [])[:5]:
            html += f"<li>{item['name']} <small style='color:#64748b;'>({item['evidence']})</small></li>"
        if not go_data.get('biological_process'):
            html += "<li>None annotated</li>"
        html += """
          </ul>
        </div>
        <div class="stat-box">
          <div class="stat-label" style="color: #d97706;">Cellular Component (GO)</div>
          <ul style="padding-left: 18px; font-size: 12px; margin-top: 8px;">
"""
        for item in go_data.get('cellular_component', [])[:5]:
            html += f"<li>{item['name']} <small style='color:#64748b;'>({item['evidence']})</small></li>"
        if not go_data.get('cellular_component'):
            html += "<li>None annotated</li>"
        html += """
          </ul>
        </div>
      </div>
"""
    else:
        html += "<p style='font-size: 13px; color: #64748b;'>No external functional database annotations available for this specific query sequence.</p>"

    # SECTION 3: CALCULATED PHYSICOCHEMICAL PROPERTIES
    html += f"""
      <h2 class="section-title">
        3. Calculated Sequence Characteristics
        <span class="badge badge-calc">Physicochemical Engine</span>
      </h2>
      <div class="grid-3" style="margin-bottom: 20px;">
        <div class="stat-box">
          <div class="stat-label">Molecular Weight</div>
          <div class="stat-val">{phys.get('molecular_weight_kda', 'N/A')} <span style="font-size:13px; font-weight:normal;">kDa</span></div>
          <div class="stat-desc">{phys.get('molecular_weight_da', 'N/A')} Daltons</div>
        </div>
        <div class="stat-box">
          <div class="stat-label">Theoretical pI</div>
          <div class="stat-val">{phys.get('isoelectric_point', 'N/A')}</div>
          <div class="stat-desc">Net charge at pH 7: {phys.get('net_charge_ph7', 'N/A')}</div>
        </div>
        <div class="stat-box">
          <div class="stat-label">GRAVY (Hydropathicity)</div>
          <div class="stat-val">{phys.get('gravy', 'N/A')}</div>
          <div class="stat-desc">{phys.get('hydropathy_nature', 'N/A')} profile</div>
        </div>
        <div class="stat-box">
          <div class="stat-label">Instability Index</div>
          <div class="stat-val">{phys.get('instability_index', 'N/A')}</div>
          <div class="stat-desc">{phys.get('stability_class', 'N/A')}</div>
        </div>
        <div class="stat-box">
          <div class="stat-label">Aliphatic Index</div>
          <div class="stat-val">{phys.get('aliphatic_index', 'N/A')}</div>
          <div class="stat-desc">Ikai (1980) thermo-stability metric</div>
        </div>
        <div class="stat-box">
          <div class="stat-label">Extinction Coeff (280nm)</div>
          <div class="stat-val">{phys.get('extinction', {}).get('extinction_reduced', 'N/A')}</div>
          <div class="stat-desc">M<sup>-1</sup> cm<sup>-1</sup> (Reduced state)</div>
        </div>
      </div>
"""

    # SECTION 4: COMPUTATIONAL AI PREDICTIONS
    html += f"""
      <h2 class="section-title">
        4. Computational Function & Localization Predictions
        <span class="badge badge-pred">Machine Learning Module</span>
      </h2>
"""
    enz_p = pred.get('enzyme_prediction', {})
    loc_p = pred.get('subcellular_localization', {})
    struct_p = pred.get('structural_motifs', {})

    html += f"""
      <div class="grid-2" style="margin-bottom: 20px;">
        <div class="stat-box">
          <div class="stat-label">Enzyme Classification Prediction</div>
          <div class="stat-val" style="font-size: 16px; margin: 6px 0; color: #1e3a8a;">{enz_p.get('label', 'N/A')}</div>
          <div class="stat-desc">Model Confidence: <strong>{enz_p.get('confidence', 'N/A')}%</strong></div>
          <div style="font-size: 12px; margin-top: 8px; color: #475569;">
            Top EC Class: <strong>{enz_p.get('top_ec_class', {}).get('name', 'N/A') if enz_p.get('top_ec_class') else 'N/A (Non-enzyme)'}</strong>
          </div>
        </div>
        <div class="stat-box">
          <div class="stat-label">Subcellular Compartment Prediction</div>
          <div class="stat-val" style="font-size: 16px; margin: 6px 0; color: #0f766e;">{loc_p.get('top_compartment', {}).get('name', 'N/A') if loc_p.get('top_compartment') else 'N/A'}</div>
          <div class="stat-desc">Localization Probability: <strong>{loc_p.get('top_compartment', {}).get('probability', 'N/A') if loc_p.get('top_compartment') else 'N/A'}%</strong></div>
          <div style="font-size: 12px; margin-top: 8px; color: #475569;">
            Transmembrane Helices: <strong>{struct_p.get('tm_count', 0)} detected</strong>
          </div>
        </div>
      </div>
"""

    # SECTION 5: SUBMITTED SEQUENCE
    html += f"""
      <h2 class="section-title">5. Submitted Protein Sequence</h2>
      <div class="sequence-block">{seq}</div>

      <!-- SECTION 6: SCIENTIFIC LIMITATIONS -->
      <h2 class="section-title" style="color: #b45309;">6. Scientific Limitations & Methodology Disclosures</h2>
      <div class="callout callout-warning">
        <ul style="padding-left: 20px; font-size: 12px; line-height: 1.7;">
          <li><strong>Distinction of Data Tiers:</strong> Physicochemical characteristics (MW, pI, GRAVY, instability index) are mathematically computed directly from the input sequence. Database annotations are retrieved from official UniProtKB REST endpoints. Pairwise alignments are computed using BLOSUM62 dynamic programming.</li>
          <li><strong>Model Scope:</strong> Computational function classifications and subcellular localization predictions are inferred using a 35-descriptor Random Forest model trained on Swiss-Prot benchmark reference sets. They provide exploratory hypotheses and should not be used as standalone clinical or medical diagnostics.</li>
          <li><strong>Experimental Corroboration:</strong> In silico predictions must be corroborated by wet-lab biological assays (e.g., enzyme kinetics, crystallographic structure determination, mass spectrometry).</li>
        </ul>
      </div>

      <div style="margin-top: 35px; padding-top: 20px; border-top: 1px solid #e2e8f0; font-size: 11px; color: #94a3b8; text-align: center;">
        Generated by AI-Based Protein Function Prediction and Analysis Platform &bull; Academic Bioinformatics Project &bull; UniProtKB, NCBI, AlphaFold, RCSB PDB
      </div>
    </div>
  </div>
</body>
</html>"""
    return html
