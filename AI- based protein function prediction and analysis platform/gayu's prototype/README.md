# AI-Based Protein Function Prediction and Analysis Platform

An academic bioinformatics web application designed to explore, characterize, and predict the biological functions of proteins from amino acid sequences.

Developed as a pair-programming academic project prototype, this platform bridges **biochemical sequence analysis**, **curated biological database integration (UniProtKB)**, **real pairwise sequence alignment**, and **machine learning functional classification** within a responsive dashboard.

---

## 🌟 Key Features

1. **Calculated Sequence Characteristics (Deterministic In Silico Engine)**:
   - **Molecular Weight**: Exact monoisotopic and average residue mass (Daltons and kDa).
   - **Theoretical Isoelectric Point (pI)**: Root-finding bisection on Lehninger / Bjellqvist pKa equations.
   - **Net Charge at pH 7.0**: Physiological ionization state.
   - **Instability Index**: Guruprasad et al. (1990) dipeptide stability matrix (determines in vitro half-life).
   - **Aliphatic Index**: Ikai (1980) thermo-stability metric for globular proteins.
   - **Grand Average of Hydropathicity (GRAVY)**: Kyte & Doolittle (1982) hydropathy profile.
   - **Molar Extinction Coefficients**: Pace et al. (1995) absorbance at 280 nm (reduced and with cystine bridges).
   - **Secondary Structure Propensities**: Chou-Fasman alpha-helix, beta-sheet, and turn conformational fractions.

2. **Biological Database Integration (UniProtKB REST API)**:
   - Direct querying of official UniProt REST API (`https://rest.uniprot.org/uniprotkb`).
   - Retrieves verified protein names, accessions, gene symbols, organism taxonomy, and functional descriptions.
   - Distinct categorization of **Reviewed (Swiss-Prot)** vs **Unreviewed (TrEMBL)** entries.
   - Structured **Gene Ontology (GO)** terms:
     - Molecular Function (GO:MF)
     - Biological Process (GO:BP)
     - Cellular Component (GO:CC)
   - Experimental evidence codes (e.g., IDA, EXP, IEA) are preserved and displayed.

3. **Real Sequence Similarity & Alignment**:
   - Dynamic programming pairwise alignment using the standard **BLOSUM62** substitution matrix.
   - Computes actual **Sequence Identity %**, **Alignment Coverage %**, and **Karlin-Altschul statistical E-values**.
   - Formatted visual alignment snippet display (Query vs Subject with match markers).

4. **Machine Learning Function & Localization Predictor**:
   - Feature engineering pipeline extracting **35 sequence-derived physicochemical descriptors** (AAC, group ratios, charge density, GRAVY, instability, aliphatic index, Chou-Fasman fractions).
   - **Enzyme vs Non-Enzyme Classifier**: Random Forest model predicting catalytic activity with confidence scores.
   - **Enzyme Commission (EC) Class Propensity**: Predicts major enzyme classes (EC 1: Oxidoreductase to EC 7: Translocase).
   - **Subcellular Localization Predictor**: Classifies targeting across 5 cellular compartments (Cytoplasm, Nucleus, Extracellular/Secreted, Mitochondria, Membrane).
   - **Transmembrane & Signal Peptide Heuristics**: Detects transmembrane alpha-helices (Kyte-Doolittle threshold &ge; 1.6) and N-terminal secretory signal sequences.

5. **Data Visualization & 3D Structure Viewer**:
   - Interactive **Chart.js** graphs:
     - 20 Standard Amino Acid Percentage Histogram (color-coded by biochemical class)
     - Biochemical Grouping Doughnut Chart
     - Kyte-Doolittle Sliding Window Hydropathy Profile Line Chart
   - **3Dmol.js WebGL Molecular Viewer**:
     - Visualizes 3D crystal structures from RCSB PDB or predicted models from the AlphaFold DB.
     - Controls for Cartoon, Stick, Sphere representations, and rotation.

6. **Publication-Ready Downloadable Reports**:
   - **Standalone HTML Report**: Formatted with embedded printable CSS (ready to print or save as PDF).
   - **JSON Data Export**: Full structured payload for computational pipelines.
   - **FASTA Export**: Formatted sequence with standard header.

---

## 🔬 Scientific Data Separation & Transparency

In adherence to bioinformatics standards and academic integrity, the platform maintains strict categorization:

| Data Tier | Provenance / Methodology | Example Metrics |
| :--- | :--- | :--- |
| **Calculated Characteristics** | Deterministic mathematical algorithms | MW, pI, GRAVY, Instability Index, Composition |
| **Retrieved Annotations** | UniProtKB REST API endpoints | Recommended Name, Function, Gene Ontology (GO) |
| **Sequence Similarity** | BLOSUM62 dynamic programming alignment | Sequence Identity %, Coverage %, E-value |
| **Computational Predictions** | 35-descriptor Random Forest ML models | Enzyme class, Subcellular compartment, TM helices |
| **Inconclusive Results** | Explicit negative reporting | Displayed when identity &lt; 25% or unannotated |

> **Scientific Disclosure:** Computational predictions are intended for preliminary hypothesis generation and academic demonstration. They provide in silico guidance that should be validated through wet-lab biological assays.

---

## 📂 Project Structure

```
gayu's prototype/
├── app.py                     # Flask web server and REST API controllers
├── run.py                     # Easy Windows launcher (auto-opens browser)
├── requirements.txt           # Python dependencies
├── README.md                  # Comprehensive academic & technical documentation
├── bio_analysis/              # Modular bioinformatics core
│   ├── __init__.py
│   ├── validator.py           # Sequence validation, FASTA parser, ambiguous code handler
│   ├── physicochemical.py     # MW, pI root-finding, GRAVY, instability index, composition
│   ├── alignment.py           # BLOSUM62 pairwise alignment, coverage, Karlin-Altschul E-value
│   ├── uniprot_service.py     # UniProt REST API client with timeout and benchmark fallback
│   ├── train_models.py        # ML feature extraction and Random Forest training script
│   ├── predictor.py           # ML prediction engine (Enzyme, EC, Localization, TM helices)
│   └── report_generator.py    # Multi-format report builder (HTML, JSON, FASTA)
├── models/
│   └── protein_function_models.joblib # Serialized machine learning model bundle
├── static/
│   ├── css/
│   │   └── style.css          # Modern bioinformatics responsive styling
│   ├── js/
│   │   ├── app.js             # Main SPA controller, validation, and dashboard logic
│   │   ├── charts.js          # Chart.js visualization configurations
│   │   └── viewer3d.js        # 3Dmol.js WebGL structure viewer integration
│   └── data/
│       └── sample_proteins.json # Curated benchmark examples
├── templates/
│   └── index.html             # Single-page application dashboard
└── tests/
    └── test_platform.py       # Automated unit test suite (18 unit tests)
```

---

## 🚀 Installation & Running Locally (Windows)

### Prerequisites
- Python 3.10+ (tested on Python 3.14)
- Web browser (Chrome, Edge, Firefox)

### Step 1: Install Dependencies
Open PowerShell or Command Prompt in the project folder:
```powershell
pip install -r requirements.txt
```

### Step 2: Run the Application
You can launch the platform with either command:

**Option A (Auto-opens browser):**
```powershell
python run.py
```

**Option B (Standard Flask):**
```powershell
python app.py
```

Then visit `http://127.0.0.1:5000` in your web browser.

---

## 🧪 Running the Automated Test Suite

The platform includes an automated unit test suite covering sequence validation, physicochemical math, alignment, UniProt querying, ML predictions, and API endpoints:

```powershell
python -m unittest tests/test_platform.py
```

All 18 tests execute in under 0.5 seconds with 100% pass rate.

---

## 🧬 Demonstration Preset Proteins

The platform includes 6 built-in benchmark proteins ready to demo with a single click:

1. **Human Hemoglobin Subunit Beta (`P68871`)**:
   - Oxygen carrier; globin family; non-enzyme; cytoplasmic.
2. **Hen Egg White Lysozyme C (`P00698`)**:
   - Innate antibacterial hydrolase enzyme (EC 3.2.1.17); secreted with an N-terminal signal peptide.
3. **Human Insulin Precursor (`P01308`)**:
   - Peptide hormone; non-enzyme; extracellular receptor binding.
4. **Bovine Rhodopsin (`P02699`)**:
   - Canonical 7-transmembrane G-protein coupled receptor (GPCR); membrane-bound.
5. **Green Fluorescent Protein (`P42212`)**:
   - 11-stranded beta-barrel; spontaneous fluorophore formation.
6. **Yeast Cytochrome C (`P00044`)**:
   - Mitochondrial electron carrier; oxidoreductase activity; heme-binding.

---

## 📚 Key Peer-Reviewed References

- **Kyte, J., & Doolittle, R. F. (1982).** A simple method for displaying the hydropathic character of a protein. *Journal of Molecular Biology*, 157(1), 105–132.
- **Guruprasad, K., Reddy, B. V., & Pandit, M. W. (1990).** Correlation between stability of a protein and its dipeptide composition. *Protein Engineering*, 4(2), 155–161.
- **Ikai, A. (1980).** Thermostability and aliphatic index of globular proteins. *The Journal of Biochemistry*, 88(6), 1895–1898.
- **Pace, C. N., Vajdos, F., Fee, L., Grimsley, G., & Gray, T. (1995).** How to measure and predict the molar absorption coefficient of a protein. *Protein Science*, 4(11), 2411–2423.
- **Karlin, S., & Altschul, S. F. (1990).** Methods for assessing the statistical significance of molecular sequence features by using general scoring schemes. *PNAS*, 87(6), 2264–2268.
- **The UniProt Consortium (2023).** UniProt: the universal protein knowledgebase in 2023. *Nucleic Acids Research*, 51(D1), D523–D531.
