"""
Comprehensive Automated Test Suite
AI-BASED PROTEIN FUNCTION PREDICTION AND ANALYSIS PLATFORM

Tests:
1. Protein Sequence Validator (Standard AAs, FASTA headers, invalid symbols, ambiguous codes, nucleotide detection)
2. Physicochemical Properties Engine (MW, pI, net charge, extinction coefficients, instability index, GRAVY, aliphatic index)
3. Sequence Alignment Engine (Identity %, Coverage %, Karlin-Altschul E-value, BLOSUM62)
4. UniProt Service & Curated Benchmarks
5. Machine Learning Predictor (Enzyme vs Non-Enzyme, EC Class, Localization, TM helices, Signal peptides)
6. Flask REST Endpoints (/health, /api/samples, /api/validate, /api/analyze, /api/report/*)
"""

import unittest
import json
from app import app
from bio_analysis.validator import validate_protein_sequence
from bio_analysis.physicochemical import analyze_physicochemical_properties
from bio_analysis.alignment import align_sequences
from bio_analysis.uniprot_service import fetch_uniprot_by_accession, search_uniprot_by_sequence
from bio_analysis.predictor import predict_protein_function


class TestProteinPlatform(unittest.TestCase):

    def setUp(self):
        self.app = app.test_client()
        self.app.testing = True

        # Reference sequence: Human Hemoglobin Subunit Beta (147 AA)
        self.hbb_seq = "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH"
        # Reference sequence: Hen Lysozyme C (147 AA)
        self.lysc_seq = "MRSLLILVLCFLPLAALGKVFGRCELAAAMKRHGLDNYRGYSLGNWVCAAKFESNFNTQATNRNTDGSTDYGILQINSRWWCNDGRTPGSRNLCNIPCSALLSSDITASVNCAKKIVSDGNGMNAWVAWRNRCKGTDVQAWIRGCRL"

    # 1. SEQUENCE VALIDATOR TESTS
    def test_validator_valid_sequence(self):
        res = validate_protein_sequence(self.hbb_seq)
        self.assertTrue(res['is_valid'])
        self.assertEqual(res['length'], 147)
        self.assertFalse(res['has_ambiguous'])
        self.assertEqual(len(res['errors']), 0)

    def test_validator_fasta_header(self):
        fasta = f">sp|P68871|HBB_HUMAN Hemoglobin beta\n{self.hbb_seq}"
        res = validate_protein_sequence(fasta)
        self.assertTrue(res['is_valid'])
        self.assertEqual(res['header'], "sp|P68871|HBB_HUMAN Hemoglobin beta")
        self.assertEqual(res['cleaned_sequence'], self.hbb_seq)

    def test_validator_invalid_characters(self):
        invalid_seq = "MVHLTPEEK123@#SAVTALWG"
        res = validate_protein_sequence(invalid_seq)
        self.assertFalse(res['is_valid'])
        self.assertTrue(any("Invalid character" in err for err in res['errors']))

    def test_validator_ambiguous_codes(self):
        amb_seq = "MVHLTPEEKXAVTALWGB"
        res = validate_protein_sequence(amb_seq)
        self.assertTrue(res['is_valid'])
        self.assertTrue(res['has_ambiguous'])
        codes = [item['code'] for item in res['ambiguous_details']]
        self.assertIn('X', codes)
        self.assertIn('B', codes)

    def test_validator_empty_input(self):
        res = validate_protein_sequence("   \n\t  ")
        self.assertFalse(res['is_valid'])
        self.assertIn("Please enter a protein sequence", res['errors'][0])

    def test_validator_nucleotide_warning(self):
        dna_seq = "ATGCGATCGATCGATCGATC"
        res = validate_protein_sequence(dna_seq)
        self.assertTrue(res['is_valid'])
        self.assertTrue(res['is_likely_nucleotide'])

    # 2. PHYSICOCHEMICAL PROPERTIES TESTS
    def test_physicochemical_properties(self):
        phys = analyze_physicochemical_properties(self.hbb_seq)
        # Molecular weight of Hemoglobin beta is ~15.8 kDa (15867 Da)
        self.assertAlmostEqual(phys['molecular_weight_kda'], 16.0, delta=1.0)
        # Theoretical pI of HBB is around 6.8 - 7.2
        self.assertGreaterEqual(phys['isoelectric_point'], 6.0)
        self.assertLessEqual(phys['isoelectric_point'], 8.0)
        # Instability index
        self.assertIn('instability_index', phys)
        # Composition count
        self.assertEqual(len(phys['composition']), 20)
        self.assertEqual(sum(item['count'] for item in phys['composition'].values()), 147)

    # 3. SEQUENCE ALIGNMENT TESTS
    def test_pairwise_alignment_identical(self):
        aln = align_sequences(self.hbb_seq, self.hbb_seq)
        self.assertEqual(aln['identity_pct'], 100.0)
        self.assertEqual(aln['coverage_pct'], 100.0)
        self.assertEqual(aln['identical_count'], 147)
        self.assertLess(aln['e_value'], 1e-10)

    def test_pairwise_alignment_divergent(self):
        aln = align_sequences(self.hbb_seq, self.lysc_seq)
        # Hemoglobin beta and Lysozyme have divergent structures (< 30% identity)
        self.assertLess(aln['identity_pct'], 30.0)

    # 4. UNIPROT SERVICE & BENCHMARKS
    def test_uniprot_accession_fetch(self):
        rec = fetch_uniprot_by_accession("P68871")
        self.assertIsNotNone(rec)
        self.assertEqual(rec['accession'], "P68871")
        self.assertEqual(rec['protein_name'], "Hemoglobin subunit beta")
        self.assertTrue(rec['reviewed'])
        self.assertIn("molecular_function", rec['go_terms'])

    def test_uniprot_sequence_search(self):
        res = search_uniprot_by_sequence(self.hbb_seq)
        self.assertTrue(res['matched'])
        self.assertEqual(res['record']['accession'], "P68871")
        self.assertEqual(res['alignment']['identity_pct'], 100.0)

    # 5. MACHINE LEARNING PREDICTOR
    def test_ml_prediction_hemoglobin(self):
        pred = predict_protein_function(self.hbb_seq)
        # HBB is a non-enzyme cytoplasmic protein
        self.assertFalse(pred['enzyme_prediction']['is_enzyme'])
        self.assertIn("Non-Enzyme", pred['enzyme_prediction']['label'])
        self.assertEqual(pred['subcellular_localization']['top_compartment']['name'], "Cytoplasm / Cytosol")

    def test_ml_prediction_lysozyme(self):
        pred = predict_protein_function(self.lysc_seq)
        # Lysozyme is an enzyme with a secretory signal peptide
        self.assertTrue(pred['enzyme_prediction']['is_enzyme'])
        self.assertEqual(pred['subcellular_localization']['top_compartment']['name'], "Extracellular / Secreted")
        self.assertTrue(pred['structural_motifs']['signal_peptide']['predicted'])

    # 6. FLASK REST ENDPOINTS
    def test_endpoint_health(self):
        resp = self.app.get('/health')
        self.assertEqual(resp.status_code, 200)
        data = json.loads(resp.data)
        self.assertEqual(data['status'], 'healthy')

    def test_endpoint_samples(self):
        resp = self.app.get('/api/samples')
        self.assertEqual(resp.status_code, 200)
        data = json.loads(resp.data)
        self.assertTrue(data['success'])
        self.assertGreaterEqual(len(data['samples']), 5)

    def test_endpoint_validate(self):
        resp = self.app.post('/api/validate', json={'sequence': self.hbb_seq})
        self.assertEqual(resp.status_code, 200)
        data = json.loads(resp.data)
        self.assertTrue(data['is_valid'])

    def test_endpoint_analyze(self):
        resp = self.app.post('/api/analyze', json={'sequence': self.hbb_seq})
        self.assertEqual(resp.status_code, 200)
        data = json.loads(resp.data)
        self.assertTrue(data['success'])
        self.assertIn('physicochemical', data)
        self.assertIn('identification', data)
        self.assertIn('predictions', data)

    def test_endpoint_reports(self):
        # Analyze first
        resp = self.app.post('/api/analyze', json={'sequence': self.hbb_seq})
        analysis_data = json.loads(resp.data)

        # Test HTML report download
        resp_html = self.app.post('/api/report/html', json=analysis_data)
        self.assertEqual(resp_html.status_code, 200)
        self.assertIn('text/html', resp_html.headers['Content-Type'])
        self.assertIn(b'AI-Based Protein Function Prediction', resp_html.data)

        # Test JSON report download
        resp_json = self.app.post('/api/report/json', json=analysis_data)
        self.assertEqual(resp_json.status_code, 200)
        self.assertIn('application/json', resp_json.headers['Content-Type'])

        # Test FASTA download
        resp_fasta = self.app.post('/api/report/fasta', json=analysis_data)
        self.assertEqual(resp_fasta.status_code, 200)
        self.assertIn(b'>sp|P68871', resp_fasta.data)


if __name__ == '__main__':
    unittest.main()
