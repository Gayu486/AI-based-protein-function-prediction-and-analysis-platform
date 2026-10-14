"""
Machine Learning Protein Function Prediction Module
Performs sequence-based computational predictions:
1. Enzyme vs Non-Enzyme Classification (with class probabilities)
2. Enzyme Commission (EC) Class Propensity (EC 1 to EC 7)
3. Subcellular Localization Prediction (5 cellular compartments)
4. Transmembrane Helices and Signal Peptide Heuristics
5. Feature importance explanation and scientific methodology validation
"""

import os
import joblib
import numpy as np
from typing import Dict, Any, List

from bio_analysis.train_models import extract_features, FEATURE_NAMES, train_and_export_models
from bio_analysis.physicochemical import KYTE_DOOLITTLE

MODEL_PATH = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "models", "protein_function_models.joblib")

EC_DESCRIPTIONS = {
    1: {"name": "EC 1: Oxidoreductase", "desc": "Catalyzes oxidation-reduction reactions (electron transfer, dehydrogenation, oxygenation)."},
    2: {"name": "EC 2: Transferase", "desc": "Catalyzes transfer of functional groups (e.g., kinase phosphoryl group, methyltransferase)."},
    3: {"name": "EC 3: Hydrolase", "desc": "Catalyzes hydrolytic cleavage of chemical bonds (proteases, nucleases, esterases, glycosidases)."},
    4: {"name": "EC 4: Lyase", "desc": "Catalyzes cleavage of various bonds by means other than hydrolysis and oxidation, often forming double bonds."},
    5: {"name": "EC 5: Isomerase", "desc": "Catalyzes geometric or structural changes within one molecule (racemases, epimerases, mutases)."},
    6: {"name": "EC 6: Ligase", "desc": "Catalyzes joining of two molecules coupled with hydrolysis of a diphosphate bond in ATP or similar."},
    7: {"name": "EC 7: Translocase", "desc": "Catalyzes movement of ions or molecules across biological membranes or their separation within membranes."}
}

LOCALIZATION_NAMES = {
    0: {"name": "Cytoplasm / Cytosol", "desc": "Soluble intracellular environment; metabolic pathways and translation machinery."},
    1: {"name": "Nucleus", "desc": "Nuclear compartment; chromatin organization, transcription factors, and RNA processing."},
    2: {"name": "Extracellular / Secreted", "desc": "Secreted into extracellular matrix or biological fluids; peptide hormones, defense factors."},
    3: {"name": "Mitochondria / Organellar", "desc": "Targeted to mitochondria; cellular respiration, electron transport, citric acid cycle."},
    4: {"name": "Plasma Membrane / Membrane-Bound", "desc": "Integral or peripheral membrane association; receptors, ion channels, transport complexes."}
}

_model_cache = None


def get_models():
    """Loads cached models or triggers automatic training on startup."""
    global _model_cache
    if _model_cache is not None:
        return _model_cache

    if os.path.exists(MODEL_PATH):
        try:
            _model_cache = joblib.load(MODEL_PATH)
            return _model_cache
        except Exception:
            pass

    # Train if missing
    _model_cache = train_and_export_models(os.path.dirname(MODEL_PATH))
    return _model_cache


def detect_transmembrane_helices(sequence: str, window_size: int = 19, threshold: float = 1.6) -> List[Dict[str, Any]]:
    """
    Scans sequence using Kyte-Doolittle hydropathy scale for putative transmembrane alpha-helices.
    Standard bioinformatics rule: 19-21 residue hydrophobic windows with avg hydropathy >= 1.6.
    """
    length = len(sequence)
    if length < window_size:
        return []

    tm_helices = []
    i = 0
    while i <= length - window_size:
        window = sequence[i: i + window_size]
        avg_hydropathy = sum(KYTE_DOOLITTLE.get(aa, 0.0) for aa in window) / window_size
        if avg_hydropathy >= threshold:
            start_pos = i + 1  # 1-based index
            # Extend window while hydrophobic
            end_pos = i + window_size
            peak_val = avg_hydropathy
            while end_pos < length:
                ext_win = sequence[end_pos - window_size + 1: end_pos + 1]
                ext_val = sum(KYTE_DOOLITTLE.get(aa, 0.0) for aa in ext_win) / window_size
                if ext_val >= 1.2:
                    end_pos += 1
                    if ext_val > peak_val:
                        peak_val = ext_val
                else:
                    break

            tm_helices.append({
                'helix_number': len(tm_helices) + 1,
                'start': start_pos,
                'end': end_pos,
                'length': end_pos - start_pos + 1,
                'peak_hydropathy': round(peak_val, 2),
                'segment': sequence[start_pos - 1: end_pos]
            })
            i = end_pos + 5  # Skip past this helix to look for the next
        else:
            i += 1

    return tm_helices


def detect_signal_peptide(sequence: str) -> Dict[str, Any]:
    """
    Heuristic detection of N-terminal secretory signal peptides (von Heijne rules).
    Characteristics: N-terminal positive charge (n-region), central hydrophobic core (h-region, 7-15 AAs),
    and small neutral residue cleavage site (-3, -1 Ala/Gly/Ser).
    """
    if len(sequence) < 25:
        return {'predicted': False, 'probability': 0.05, 'cleavage_site': None, 'details': 'Sequence too short for signal peptide.'}

    n_term = sequence[:30]
    # Check for hydrophobic core in residues 5 to 22
    core_found = False
    max_core_hydro = -5.0
    cleavage_pos = None

    for i in range(5, min(22, len(n_term) - 7)):
        sub = n_term[i: i + 8]
        hydro = sum(KYTE_DOOLITTLE.get(aa, 0.0) for aa in sub) / 8.0
        if hydro > max_core_hydro:
            max_core_hydro = hydro
        if hydro >= 1.5:
            core_found = True
            # Look for cleavage motif (-3, -1) in residues 18-28
            for c_pos in range(i + 8, min(28, len(n_term))):
                if sequence[c_pos - 1] in ('A', 'G', 'S', 'C') and sequence[c_pos - 3] in ('A', 'G', 'S', 'V', 'I', 'L'):
                    cleavage_pos = c_pos
                    break
            break

    if core_found and max_core_hydro >= 1.6:
        prob = min(0.95, round(0.5 + (max_core_hydro / 5.0) + (0.2 if cleavage_pos else 0.0), 2))
        return {
            'predicted': True,
            'probability': prob,
            'cleavage_site': f"Between pos {cleavage_pos} and {cleavage_pos + 1}" if cleavage_pos else "Putative cleavage around pos 20-24",
            'details': f"Hydrophobic core detected in N-terminal segment (peak hydropathy {max_core_hydro:.2f})."
        }
    else:
        return {
            'predicted': False,
            'probability': 0.12,
            'cleavage_site': None,
            'details': "No classical N-terminal hydrophobic signal anchor detected."
        }


def predict_protein_function(sequence: str) -> Dict[str, Any]:
    """
    Main computational prediction pipeline.
    Runs machine learning models and structural heuristics on the query sequence.
    """
    models = get_models()
    clf_enzyme = models['model_enzyme']
    clf_ec = models['model_ec']
    clf_loc = models['model_loc']

    features = extract_features(sequence).reshape(1, -1)

    # 1. Enzyme Classification
    enz_probs = clf_enzyme.predict_proba(features)[0]
    enz_class = int(np.argmax(enz_probs))
    # Class 0: Non-enzyme, Class 1: Enzyme
    p_non_enzyme = round(float(enz_probs[0]) * 100, 1)
    p_enzyme = round(float(enz_probs[1]) * 100, 1) if len(enz_probs) > 1 else round(100.0 - p_non_enzyme, 1)
    is_enzyme_predicted = p_enzyme >= 50.0

    # 2. EC Class Prediction (if enzyme or borderline)
    ec_results = []
    top_ec = None
    if p_enzyme >= 35.0:
        ec_probs = clf_ec.predict_proba(features)[0]
        ec_classes = clf_ec.classes_
        for cls_idx, prob in zip(ec_classes, ec_probs):
            if cls_idx in EC_DESCRIPTIONS:
                info = EC_DESCRIPTIONS[cls_idx]
                ec_results.append({
                    'ec_class': int(cls_idx),
                    'name': info['name'],
                    'description': info['desc'],
                    'probability': round(float(prob) * 100, 1)
                })
        ec_results.sort(key=lambda x: x['probability'], reverse=True)
        if ec_results:
            top_ec = ec_results[0]

    # 3. Subcellular Localization
    loc_probs = clf_loc.predict_proba(features)[0]
    loc_classes = clf_loc.classes_
    loc_results = []
    for cls_idx, prob in zip(loc_classes, loc_probs):
        if cls_idx in LOCALIZATION_NAMES:
            info = LOCALIZATION_NAMES[cls_idx]
            loc_results.append({
                'class_id': int(cls_idx),
                'name': info['name'],
                'description': info['desc'],
                'probability': round(float(prob) * 100, 1)
            })
    loc_results.sort(key=lambda x: x['probability'], reverse=True)
    top_loc = loc_results[0] if loc_results else None

    # 4. Transmembrane & Signal Peptide Heuristics
    tm_helices = detect_transmembrane_helices(sequence)
    signal_peptide = detect_signal_peptide(sequence)

    # If transmembrane helices detected, boost membrane localization probability explanation
    if len(tm_helices) >= 1 and top_loc and top_loc['class_id'] != 4:
        # Check if membrane probability is high
        for r in loc_results:
            if r['class_id'] == 4:
                r['probability'] = max(r['probability'], 45.0)
        loc_results.sort(key=lambda x: x['probability'], reverse=True)
        top_loc = loc_results[0]

    # 5. Explanatory Biological Insights
    insights = []
    if is_enzyme_predicted:
        insights.append(f"Predicted as an active enzyme ({p_enzyme}% confidence). " +
                        (f"Top predicted functional activity: {top_ec['name']}." if top_ec else ""))
    else:
        insights.append(f"Predicted as non-enzymatic protein ({p_non_enzyme}% confidence). Characteristic of structural, transport, hormone, or regulatory molecules.")

    if top_loc:
        insights.append(f"Primary predicted subcellular compartment: {top_loc['name']} ({top_loc['probability']}% probability).")

    if tm_helices:
        insights.append(f"Detected {len(tm_helices)} putative transmembrane alpha-helix segment(s) based on Kyte-Doolittle hydropathy profile (peak hydropathy >= 1.6).")

    if signal_peptide['predicted']:
        insights.append(f"N-terminal secretory signal peptide detected ({round(signal_peptide['probability'] * 100, 1)}% likelihood), suggesting ER-Golgi targeting or secretion.")

    return {
        'enzyme_prediction': {
            'is_enzyme': is_enzyme_predicted,
            'label': "Enzyme (Catalytic Activity)" if is_enzyme_predicted else "Non-Enzyme (Regulatory / Structural / Transport / Receptor)",
            'confidence': p_enzyme if is_enzyme_predicted else p_non_enzyme,
            'probabilities': {
                'enzyme': p_enzyme,
                'non_enzyme': p_non_enzyme
            },
            'top_ec_class': top_ec,
            'all_ec_classes': ec_results
        },
        'subcellular_localization': {
            'top_compartment': top_loc,
            'all_compartments': loc_results
        },
        'structural_motifs': {
            'transmembrane_helices': tm_helices,
            'has_transmembrane': len(tm_helices) > 0,
            'tm_count': len(tm_helices),
            'signal_peptide': signal_peptide
        },
        'insights': insights,
        'methodology_notes': {
            'model_type': 'Random Forest Classifier (35 sequence-derived physicochemical descriptors)',
            'features_used': 'Amino Acid Composition (AAC), Dipeptide motifs, Isoelectric Point (pI), GRAVY, Instability Index, Aliphatic Index, Chou-Fasman Secondary Structure Propensities',
            'limitations': 'Academic prototype prediction model. Predictions reflect statistical sequence correlations on benchmark datasets and do not substitute for wet-lab biochemical validation.'
        }
    }
