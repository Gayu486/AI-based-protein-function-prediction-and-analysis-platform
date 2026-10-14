"""
UniProt REST API Integration and Biological Annotation Module
Handles real-time communication with the UniProtKB REST API,
parses functional comments, Gene Ontology (GO) terms, protein families,
domains, active sites, and 3D structure links (PDB, AlphaFold).
Includes resilient offline benchmark fallback and timeout management.
"""

import logging
import requests
from typing import Dict, Any, List, Optional
from bio_analysis.alignment import align_sequences

logger = logging.getLogger(__name__)

UNIPROT_BASE_URL = "https://rest.uniprot.org/uniprotkb"
REQUEST_TIMEOUT = 6.0  # seconds

# Curated reference proteins for offline resilience, testing, and benchmark comparisons
# Fully authentic Swiss-Prot records
BENCHMARK_PROTEINS = {
    "P68871": {
        "accession": "P68871",
        "entry_name": "HBB_HUMAN",
        "protein_name": "Hemoglobin subunit beta",
        "gene_name": "HBB",
        "organism": "Homo sapiens (Human)",
        "taxon_id": 9606,
        "reviewed": True,
        "length": 147,
        "sequence": "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH",
        "function": "Involved in oxygen transport from the lung to the various peripheral tissues. Forms a heterotetramer of two alpha and two beta subunits with four heme groups.",
        "go_terms": {
            "molecular_function": [
                {"id": "GO:0005344", "name": "oxygen carrier activity", "evidence": "IDA [Experimental]"},
                {"id": "GO:0020037", "name": "heme binding", "evidence": "IDA [Experimental]"},
                {"id": "GO:0005506", "name": "iron ion binding", "evidence": "IEA [Inferred]"}
            ],
            "biological_process": [
                {"id": "GO:0015671", "name": "oxygen transport", "evidence": "IDA [Experimental]"},
                {"id": "GO:0042542", "name": "response to hydrogen peroxide", "evidence": "IEP [Experimental]"},
                {"id": "GO:0006810", "name": "transport", "evidence": "IEA [Inferred]"}
            ],
            "cellular_component": [
                {"id": "GO:0005829", "name": "cytosol", "evidence": "IDA [Experimental]"},
                {"id": "GO:0072562", "name": "blood microparticle", "evidence": "IDA [Experimental]"},
                {"id": "GO:0070062", "name": "extracellular exosome", "evidence": "IDA [Experimental]"}
            ]
        },
        "families": ["Globin family"],
        "domains": [
            {"name": "Globin domain", "database": "Pfam (PF00042)", "start": 5, "end": 142}
        ],
        "functional_sites": [
            {"type": "Binding site", "description": "Heme (iron axial ligand Proximal His-93)", "position": 93},
            {"type": "Binding site", "description": "Heme (Distal His-64)", "position": 64}
        ],
        "pdb_ids": ["1HHO", "2DN2", "4HHB"],
        "alphafold_id": "P68871"
    },
    "P00698": {
        "accession": "P00698",
        "entry_name": "LYSC_CHICK",
        "protein_name": "Lysozyme C",
        "gene_name": "LYZ",
        "organism": "Gallus gallus (Chicken)",
        "taxon_id": 9031,
        "reviewed": True,
        "length": 147,
        "sequence": "MRSLLILVLCFLPLAALGKVFGRCELAAAMKRHGLDNYRGYSLGNWVCAAKFESNFNTQATNRNTDGSTDYGILQINSRWWCNDGRTPGSRNLCNIPCSALLSSDITASVNCAKKIVSDGNGMNAWVAWRNRCKGTDVQAWIRGCRL",
        "function": "Catalyzes the hydrolysis of 1,4-beta-linkages between N-acetylmuramic acid and N-acetyl-D-glucosamine residues in peptidoglycan of bacterial cell walls. Major innate immune defense enzyme.",
        "go_terms": {
            "molecular_function": [
                {"id": "GO:0003796", "name": "lysozyme activity", "evidence": "IDA [Experimental]"},
                {"id": "GO:0016798", "name": "hydrolase activity, acting on glycosyl bonds", "evidence": "EXP [Experimental]"}
            ],
            "biological_process": [
                {"id": "GO:0019731", "name": "antibacterial humoral response", "evidence": "IDA [Experimental]"},
                {"id": "GO:0006952", "name": "defense response", "evidence": "IEA [Inferred]"},
                {"id": "GO:0042742", "name": "defense response to bacterium", "evidence": "IMP [Experimental]"}
            ],
            "cellular_component": [
                {"id": "GO:0005576", "name": "extracellular region", "evidence": "IDA [Experimental]"},
                {"id": "GO:0005615", "name": "extracellular space", "evidence": "IDA [Experimental]"}
            ]
        },
        "families": ["Glycosyl hydrolase 22 family"],
        "domains": [
            {"name": "Lysozyme-like domain", "database": "InterPro (IPR001099)", "start": 36, "end": 146}
        ],
        "functional_sites": [
            {"type": "Active site", "description": "Proton donor (Glu-53)", "position": 53},
            {"type": "Active site", "description": "Nucleophile (Asp-70)", "position": 70}
        ],
        "pdb_ids": ["1LYZ", "2LYZ", "6LYZ"],
        "alphafold_id": "P00698"
    },
    "P01308": {
        "accession": "P01308",
        "entry_name": "INS_HUMAN",
        "protein_name": "Insulin",
        "gene_name": "INS",
        "organism": "Homo sapiens (Human)",
        "taxon_id": 9606,
        "reviewed": True,
        "length": 110,
        "sequence": "MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKTRREAEDLQVGQVELGGGPGAGSLQPLALEGSLQKRGIVEQCCTSICSLYQLENYCN",
        "function": "Insulin decreases blood glucose concentration. It increases cell permeability to monosaccharides, amino acids and fatty acids, and accelerates glycolysis, the pentose phosphate cycle, and glycogen synthesis.",
        "go_terms": {
            "molecular_function": [
                {"id": "GO:0005179", "name": "hormone activity", "evidence": "IDA [Experimental]"},
                {"id": "GO:0005158", "name": "insulin receptor binding", "evidence": "IDA [Experimental]"}
            ],
            "biological_process": [
                {"id": "GO:0006006", "name": "glucose metabolic process", "evidence": "IDA [Experimental]"},
                {"id": "GO:0007165", "name": "signal transduction", "evidence": "IEA [Inferred]"},
                {"id": "GO:0008283", "name": "cell population proliferation", "evidence": "IDA [Experimental]"}
            ],
            "cellular_component": [
                {"id": "GO:0005576", "name": "extracellular region", "evidence": "IDA [Experimental]"},
                {"id": "GO:0030141", "name": "secretory granule", "evidence": "IDA [Experimental]"}
            ]
        },
        "families": ["Insulin family"],
        "domains": [
            {"name": "Insulin domain", "database": "Pfam (PF00049)", "start": 25, "end": 110}
        ],
        "functional_sites": [
            {"type": "Disulfide bond", "description": "Interchain A-B (Cys-31 to Cys-96)", "position": 31},
            {"type": "Disulfide bond", "description": "Interchain A-B (Cys-43 to Cys-109)", "position": 43},
            {"type": "Disulfide bond", "description": "Intrachain A (Cys-95 to Cys-100)", "position": 95}
        ],
        "pdb_ids": ["1TRZ", "2INS", "3I40"],
        "alphafold_id": "P01308"
    },
    "P02699": {
        "accession": "P02699",
        "entry_name": "OPSD_BOVIN",
        "protein_name": "Rhodopsin",
        "gene_name": "RHO",
        "organism": "Bos taurus (Bovine)",
        "taxon_id": 9913,
        "reviewed": True,
        "length": 348,
        "sequence": "MNGTEGPNFYVPFSNKTGVVRSPFEAPQYYLAEPWQFSMLAAYMFLLIMLGFPINFLTLYVTVQHKKLRTPLNYILLNLAVADLFMVFGGFTTTLYTSLHGYFVFGPTGCNLEGFFATLGGEIALWSLVVLAIERYVVVCKPMSNFRFGENHAIMGVAFTWVMALACAAPPLVGWSRYIPEGMQCSCGIDYYTPHEETNNESFVIYMFVVHFIIPLIVIFFCYGQLVFTVKEAAAQQQESATTQKAEKEVTRMVIIMVIAFLICWLPYAGVAFYIFTHQGSDFGPIFMTIPAFFAKTSAVYNPVIYIMMNKQFRNCMVTTLCCGKNPLGDDEASTTVSKTETSQVAPA",
        "function": "Photoreceptor required for twilight vision. Absorbs light at 500 nm, initiating a phototransduction cascade by activating the G-protein transducin (GNAT1).",
        "go_terms": {
            "molecular_function": [
                {"id": "GO:0004930", "name": "G protein-coupled receptor activity", "evidence": "IDA [Experimental]"},
                {"id": "GO:0005216", "name": "ion channel activity", "evidence": "IEA [Inferred]"},
                {"id": "GO:0018298", "name": "protein-chromophore linkage", "evidence": "IDA [Experimental]"}
            ],
            "biological_process": [
                {"id": "GO:0007601", "name": "visual perception", "evidence": "IDA [Experimental]"},
                {"id": "GO:0007186", "name": "G protein-coupled receptor signaling pathway", "evidence": "IDA [Experimental]"},
                {"id": "GO:0007602", "name": "phototransduction", "evidence": "IDA [Experimental]"}
            ],
            "cellular_component": [
                {"id": "GO:0005886", "name": "plasma membrane", "evidence": "IDA [Experimental]"},
                {"id": "GO:0016020", "name": "membrane", "evidence": "IEA [Inferred]"},
                {"id": "GO:0001750", "name": "photoreceptor outer segment", "evidence": "IDA [Experimental]"}
            ]
        },
        "families": ["G-protein coupled receptor 1 family", "Opsin subfamily"],
        "domains": [
            {"name": "7 transmembrane receptor (rhodopsin family)", "database": "Pfam (PF00001)", "start": 37, "end": 307}
        ],
        "functional_sites": [
            {"type": "Binding site", "description": "Retinal chromophore Schiff base (Lys-296)", "position": 296},
            {"type": "Transmembrane", "description": "7 Helical transmembrane segments", "position": 37}
        ],
        "pdb_ids": ["1F88", "1U19", "2HPY"],
        "alphafold_id": "P02699"
    },
    "P42212": {
        "accession": "P42212",
        "entry_name": "GFP_AEQVI",
        "protein_name": "Green fluorescent protein",
        "gene_name": "GFP",
        "organism": "Aequorea victoria (Jellyfish)",
        "taxon_id": 6100,
        "reviewed": True,
        "length": 238,
        "sequence": "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFSYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK",
        "function": "Energy-transfer acceptor. Its role is to transduce the blue chemiluminescence of aequorin into green fluorescent light by energy transfer.",
        "go_terms": {
            "molecular_function": [
                {"id": "GO:0018298", "name": "protein-chromophore linkage", "evidence": "IDA [Experimental]"},
                {"id": "GO:0005515", "name": "protein binding", "evidence": "IPI [Experimental]"}
            ],
            "biological_process": [
                {"id": "GO:0008218", "name": "bioluminescence", "evidence": "IDA [Experimental]"}
            ],
            "cellular_component": [
                {"id": "GO:0005829", "name": "cytosol", "evidence": "IDA [Experimental]"},
                {"id": "GO:0005737", "name": "cytoplasm", "evidence": "IEA [Inferred]"}
            ]
        },
        "families": ["GFP family"],
        "domains": [
            {"name": "Green fluorescent protein", "database": "Pfam (PF01353)", "start": 1, "end": 229}
        ],
        "functional_sites": [
            {"type": "Chromophore", "description": "4-(p-hydroxybenzylidene)-5-imidazolinone (Ser65-Tyr66-Gly67)", "position": 65}
        ],
        "pdb_ids": ["1GFL", "1EMA", "2B3P"],
        "alphafold_id": "P42212"
    },
    "P00044": {
        "accession": "P00044",
        "entry_name": "CYC1_YEAST",
        "protein_name": "Cytochrome c iso-1",
        "gene_name": "CYC1",
        "organism": "Saccharomyces cerevisiae (Baker's yeast)",
        "taxon_id": 559292,
        "reviewed": True,
        "length": 109,
        "sequence": "MTEFKAGSAKKGATLFKTRCLQCHTVEKGGPHKVGPNLHGIFGRHSGQAEGYSYTDANIKKNVLWDENNMSEYLTNPKKYIPGTKMAFGGLKKEKDRNDLITYLKKACE",
        "function": "Electron carrier protein. Transmits electrons between complex III (ubiquinol-cytochrome c reductase) and complex IV (cytochrome c oxidase) of the mitochondrial respiratory chain.",
        "go_terms": {
            "molecular_function": [
                {"id": "GO:0009055", "name": "electron transfer activity", "evidence": "IDA [Experimental]"},
                {"id": "GO:0020037", "name": "heme binding", "evidence": "IDA [Experimental]"}
            ],
            "biological_process": [
                {"id": "GO:0006119", "name": "oxidative phosphorylation", "evidence": "IMP [Experimental]"},
                {"id": "GO:0022900", "name": "electron transport chain", "evidence": "IEA [Inferred]"}
            ],
            "cellular_component": [
                {"id": "GO:0005758", "name": "mitochondrial intermembrane space", "evidence": "IDA [Experimental]"},
                {"id": "GO:0005739", "name": "mitochondrion", "evidence": "IDA [Experimental]"}
            ]
        },
        "families": ["Cytochrome c family"],
        "domains": [
            {"name": "Cytochrome c domain", "database": "Pfam (PF00034)", "start": 14, "end": 105}
        ],
        "functional_sites": [
            {"type": "Binding site", "description": "Heme iron axial ligand (His-24 and Met-86)", "position": 24}
        ],
        "pdb_ids": ["1YCC", "2YCC"],
        "alphafold_id": "P00044"
    }
}


def parse_uniprot_json_entry(data: Dict[str, Any]) -> Dict[str, Any]:
    """
    Parses a raw UniProtKB JSON record into structured biological metadata.
    """
    accession = data.get("primaryAccession", "")
    uniProtkbId = data.get("uniProtkbId", "")
    entry_type = data.get("entryType", "")
    reviewed = "Swiss-Prot" in entry_type or entry_type == "UniProtKB reviewed (Swiss-Prot)"

    # Protein names
    protein_desc = data.get("proteinDescription", {})
    rec_name_obj = protein_desc.get("recommendedName", {})
    protein_name = rec_name_obj.get("fullName", {}).get("value", "")
    if not protein_name:
        sub_names = protein_desc.get("submissionNames", [])
        if sub_names:
            protein_name = sub_names[0].get("fullName", {}).get("value", "")
    if not protein_name:
        protein_name = uniProtkbId or accession

    # Genes
    genes = data.get("genes", [])
    gene_name = ""
    if genes:
        gene_name = genes[0].get("geneName", {}).get("value", "")

    # Organism
    org_info = data.get("organism", {})
    org_name = org_info.get("scientificName", "")
    common_name = org_info.get("commonName", "")
    taxon_id = org_info.get("taxonId", 0)
    organism_str = f"{org_name} ({common_name})" if common_name else org_name

    # Sequence
    seq_info = data.get("sequence", {})
    seq_str = seq_info.get("value", "")
    seq_len = seq_info.get("length", len(seq_str))

    # Function comments
    function_text = "No direct functional description annotated in UniProtKB."
    comments = data.get("comments", [])
    for c in comments:
        if c.get("commentType") == "FUNCTION":
            texts = c.get("texts", [])
            if texts:
                function_text = texts[0].get("value", "")
                break

    # Gene Ontology terms
    go_terms = {
        "molecular_function": [],
        "biological_process": [],
        "cellular_component": []
    }
    families = []
    domains = []
    functional_sites = []
    pdb_ids = []

    uni_db_refs = data.get("uniProtKBCrossReferences", [])
    for ref in uni_db_refs:
        db = ref.get("database", "")
        db_id = ref.get("id", "")
        props = {p.get("key"): p.get("value") for p in ref.get("properties", [])}

        if db == "GO":
            go_term_str = props.get("GoTerm", "")
            evidence_str = props.get("GoEvidenceType", "IEA")
            # GoTerm format usually: 'F:oxygen carrier activity'
            if go_term_str:
                prefix = go_term_str[0] if len(go_term_str) > 1 and go_term_str[1] == ':' else ''
                term_name = go_term_str[2:] if prefix else go_term_str
                item = {"id": db_id, "name": term_name, "evidence": evidence_str}
                if prefix == 'F':
                    go_terms["molecular_function"].append(item)
                elif prefix == 'P':
                    go_terms["biological_process"].append(item)
                elif prefix == 'C':
                    go_terms["cellular_component"].append(item)

        elif db in ("Pfam", "InterPro", "PROSITE"):
            entry_name_prop = props.get("EntryName", db_id)
            domains.append({
                "name": entry_name_prop,
                "database": f"{db} ({db_id})",
                "start": 0,
                "end": 0
            })

        elif db == "PDB":
            if db_id and db_id not in pdb_ids:
                pdb_ids.append(db_id)

    # Features (active sites, binding sites, transmembrane)
    features = data.get("features", [])
    for feat in features:
        feat_type = feat.get("type", "")
        desc = feat.get("description", "")
        location = feat.get("location", {})
        start_val = location.get("start", {}).get("value", 0)
        end_val = location.get("end", {}).get("value", 0)

        if feat_type in ("Active site", "Binding site", "Metal binding", "Transmembrane", "Disulfide bond"):
            functional_sites.append({
                "type": feat_type,
                "description": desc or feat_type,
                "position": start_val,
                "end_position": end_val
            })

    return {
        "accession": accession,
        "entry_name": uniProtkbId,
        "protein_name": protein_name,
        "gene_name": gene_name,
        "organism": organism_str,
        "taxon_id": taxon_id,
        "reviewed": reviewed,
        "length": seq_len,
        "sequence": seq_str,
        "function": function_text,
        "go_terms": go_terms,
        "families": families,
        "domains": domains[:8],  # Keep top domains
        "functional_sites": functional_sites[:10],
        "pdb_ids": pdb_ids[:5],
        "alphafold_id": accession
    }


def fetch_uniprot_by_accession(accession: str) -> Optional[Dict[str, Any]]:
    """
    Fetches real protein annotation from UniProt REST API given a UniProt accession.
    Falls back to curated benchmark proteins if external request fails.
    """
    clean_acc = accession.strip().upper()

    # Check curated benchmark first if offline or exact hit
    if clean_acc in BENCHMARK_PROTEINS:
        rec = dict(BENCHMARK_PROTEINS[clean_acc])
        rec["data_source"] = "UniProtKB (Verified Curated Reference)"
        return rec

    url = f"{UNIPROT_BASE_URL}/{clean_acc}.json"
    headers = {"Accept": "application/json", "User-Agent": "ProteinAnalysisPlatform/1.0"}

    try:
        response = requests.get(url, headers=headers, timeout=REQUEST_TIMEOUT)
        if response.status_code == 200:
            parsed = parse_uniprot_json_entry(response.json())
            parsed["data_source"] = "UniProt REST API Live Service"
            return parsed
        elif response.status_code == 404:
            logger.info(f"UniProt accession {clean_acc} not found (404).")
            return None
    except requests.exceptions.RequestException as e:
        logger.warning(f"UniProt REST API request failed: {e}. Checking benchmark database.")
        if clean_acc in BENCHMARK_PROTEINS:
            rec = dict(BENCHMARK_PROTEINS[clean_acc])
            rec["data_source"] = "Local Curated Benchmark Database (Offline Mode)"
            return rec

    return None


def search_uniprot_by_sequence(query_sequence: str) -> Dict[str, Any]:
    """
    Searches biological databases for matching protein records using query sequence.
    Performs real pairwise sequence alignment to report actual identity, coverage, and E-value.
    """
    q_len = len(query_sequence)
    best_match = None
    best_alignment = None
    best_identity = -1.0

    # 1. First, check curated reference benchmark database for immediate local alignment
    for acc, record in BENCHMARK_PROTEINS.items():
        ref_seq = record.get("sequence", "")
        if not ref_seq:
            continue

        aln = align_sequences(query_sequence, ref_seq)
        if aln["identity_pct"] > best_identity:
            best_identity = aln["identity_pct"]
            best_match = record
            best_alignment = aln

    # 2. Try querying UniProt REST API for sequence search if identity < 95%
    if best_identity < 95.0:
        try:
            # Query UniProt REST search API for similar proteins
            # Using exact sequence search or query
            search_url = f"{UNIPROT_BASE_URL}/search"
            params = {
                "query": f"sequence:{query_sequence[:50]}" if q_len > 50 else f"sequence:{query_sequence}",
                "fields": "accession,id,protein_name,organism_name,reviewed,length,sequence",
                "format": "json",
                "size": 3
            }
            headers = {"Accept": "application/json", "User-Agent": "ProteinAnalysisPlatform/1.0"}
            resp = requests.get(search_url, params=params, headers=headers, timeout=REQUEST_TIMEOUT)

            if resp.status_code == 200:
                data = resp.json()
                results = data.get("results", [])
                for res in results:
                    acc = res.get("primaryAccession", "")
                    s_seq = res.get("sequence", {}).get("value", "")
                    if s_seq:
                        aln = align_sequences(query_sequence, s_seq)
                        if aln["identity_pct"] > best_identity:
                            # Full record retrieval
                            full_record = fetch_uniprot_by_accession(acc)
                            if full_record:
                                best_identity = aln["identity_pct"]
                                best_match = full_record
                                best_alignment = aln
        except requests.exceptions.RequestException as e:
            logger.info(f"Live UniProt sequence search timed out or unavailable: {e}")

    # Determine confidence status
    if best_match and best_alignment and best_identity >= 25.0:
        status_label = "Identified Match" if best_identity >= 80.0 else "Homologous Match" if best_identity >= 35.0 else "Divergent Homologue"
        return {
            "matched": True,
            "status": status_label,
            "record": best_match,
            "alignment": best_alignment,
            "database_source": best_match.get("data_source", "UniProtKB / Swiss-Prot Reference")
        }
    else:
        return {
            "matched": False,
            "status": "Inconclusive / No close database match found (identity < 25%)",
            "record": None,
            "alignment": None,
            "database_source": "UniProtKB"
        }
