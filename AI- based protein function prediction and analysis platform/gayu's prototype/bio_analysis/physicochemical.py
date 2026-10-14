"""
Physicochemical Properties Calculator Module
Implements accurate, peer-reviewed bioinformatics algorithms for protein sequence analysis:
- Molecular weight and atomic composition
- Theoretical Isoelectric Point (pI) via iterative root-finding
- Net charge at physiological pH (pH 7.0)
- Molar extinction coefficients (reduced and with cystine bridges)
- Instability index (Guruprasad et al., 1990)
- Aliphatic index (Ikai, 1980)
- Grand Average of Hydropathicity (GRAVY; Kyte & Doolittle, 1982)
- Kyte-Doolittle sliding window hydropathy profile
- Secondary structure conformational propensities (Chou-Fasman)
- Amino acid distribution and biochemical grouping
"""

from typing import Dict, Any, List, Tuple
from bio_analysis.validator import STANDARD_AMINO_ACIDS, AA_NAMES

# Monoisotopic residue molecular weights in Daltons (amino acid minus water H2O)
RESIDUE_WEIGHTS = {
    'A': 71.0788,
    'R': 156.1875,
    'N': 114.1038,
    'D': 115.0886,
    'C': 103.1388,
    'E': 129.1155,
    'Q': 128.1307,
    'G': 57.0519,
    'H': 137.1411,
    'I': 113.1594,
    'L': 113.1594,
    'K': 128.1741,
    'M': 131.1926,
    'F': 147.1766,
    'P': 97.1167,
    'S': 87.0782,
    'T': 101.1051,
    'W': 186.2132,
    'Y': 163.1760,
    'V': 99.1326,
    # Ambiguous approximations:
    'B': 114.5962,  # avg of D and N
    'Z': 128.6231,  # avg of E and Q
    'X': 110.0,     # average residue mass
    'J': 113.1594,  # L or I
    'U': 150.0388,  # Selenocysteine
    'O': 237.3000   # Pyrrolysine
}
WATER_WEIGHT = 18.01528  # Free amino and carboxyl termini

# Kyte-Doolittle Hydropathy scale (Kyte & Doolittle, 1982)
KYTE_DOOLITTLE = {
    'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5,
    'Q': -3.5, 'E': -3.5, 'G': -0.4, 'H': -3.2, 'I': 4.5,
    'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8, 'P': -1.6,
    'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2,
    'B': -3.5, 'Z': -3.5, 'X': 0.0, 'J': 4.15, 'U': 0.0, 'O': 0.0
}

# Chou-Fasman secondary structure conformational parameters
# P(a): alpha-helix, P(b): beta-sheet, P(turn): beta-turn
CHOU_FASMAN = {
    'A': {'a': 1.42, 'b': 0.83, 'turn': 0.66},
    'R': {'a': 0.98, 'b': 0.93, 'turn': 0.95},
    'N': {'a': 0.67, 'b': 0.89, 'turn': 1.56},
    'D': {'a': 1.01, 'b': 0.54, 'turn': 1.46},
    'C': {'a': 0.70, 'b': 1.19, 'turn': 1.19},
    'Q': {'a': 1.11, 'b': 1.10, 'turn': 0.98},
    'E': {'a': 1.51, 'b': 0.37, 'turn': 0.74},
    'G': {'a': 0.57, 'b': 0.75, 'turn': 1.56},
    'H': {'a': 1.00, 'b': 0.87, 'turn': 0.95},
    'I': {'a': 1.08, 'b': 1.60, 'turn': 0.47},
    'L': {'a': 1.21, 'b': 1.30, 'turn': 0.59},
    'K': {'a': 1.16, 'b': 0.74, 'turn': 1.01},
    'M': {'a': 1.45, 'b': 1.05, 'turn': 0.60},
    'F': {'a': 1.13, 'b': 1.38, 'turn': 0.60},
    'P': {'a': 0.57, 'b': 0.55, 'turn': 1.52},
    'S': {'a': 0.77, 'b': 0.75, 'turn': 1.43},
    'T': {'a': 0.83, 'b': 1.19, 'turn': 0.96},
    'W': {'a': 1.08, 'b': 1.37, 'turn': 0.96},
    'Y': {'a': 0.69, 'b': 1.47, 'turn': 1.14},
    'V': {'a': 1.06, 'b': 1.70, 'turn': 0.50},
    'B': {'a': 0.84, 'b': 0.71, 'turn': 1.51},
    'Z': {'a': 1.31, 'b': 0.73, 'turn': 0.86},
    'X': {'a': 1.00, 'b': 1.00, 'turn': 1.00},
    'J': {'a': 1.14, 'b': 1.45, 'turn': 0.54},
    'U': {'a': 0.70, 'b': 1.19, 'turn': 1.19},
    'O': {'a': 1.16, 'b': 0.74, 'turn': 1.01}
}

# Dipeptide Instability Weight Values (DIWV) from Guruprasad et al., 1990
# 400 dipeptide stability matrix condensed dictionary
DIWV_DATA = {
    'W': {'W': 1.0, 'C': 1.0, 'M': 24.68, 'H': 24.68, 'K': 1.0, 'R': 1.0, 'Q': 1.0, 'E': 1.0, 'T': -14.03, 'S': 1.0, 'N': 13.34, 'A': -14.03, 'D': 1.0, 'F': 1.0, 'Y': 1.0, 'V': -7.49, 'I': 1.0, 'P': 1.0, 'L': 13.34, 'G': -7.49},
    'C': {'W': 24.68, 'C': 1.0, 'M': 33.6, 'H': 33.6, 'K': 1.0, 'R': 1.0, 'Q': -6.54, 'E': 1.0, 'T': 33.6, 'S': 1.0, 'N': 1.0, 'A': -1.4, 'D': 20.26, 'F': 1.0, 'Y': 1.0, 'V': -6.54, 'I': 1.0, 'P': 20.26, 'L': 20.26, 'G': 1.0},
    'M': {'W': 1.0, 'C': 1.0, 'M': -1.4, 'H': 58.28, 'K': 1.0, 'R': -6.54, 'Q': -6.54, 'E': 1.0, 'T': -1.4, 'S': 44.94, 'N': 1.0, 'A': 13.34, 'D': 1.0, 'F': 1.0, 'Y': 24.68, 'V': 1.0, 'I': 1.0, 'P': 44.94, 'L': 1.0, 'G': 1.0},
    'H': {'W': -1.4, 'C': 1.0, 'M': 1.0, 'H': 1.0, 'K': 24.68, 'R': 1.0, 'Q': 1.0, 'E': 1.0, 'T': -6.54, 'S': 1.0, 'N': 24.68, 'A': 1.0, 'D': 1.0, 'F': -9.37, 'Y': 44.94, 'V': 1.0, 'I': 44.94, 'P': -1.4, 'L': 1.0, 'G': -9.37},
    'K': {'W': 1.0, 'C': 1.0, 'M': 33.6, 'H': 1.0, 'K': 1.0, 'R': 33.6, 'Q': 24.68, 'E': 1.0, 'T': 1.0, 'S': 1.0, 'N': 1.0, 'A': 1.0, 'D': 1.0, 'F': 1.0, 'Y': 1.0, 'V': -7.49, 'I': -7.49, 'P': -6.54, 'L': -7.49, 'G': -7.49},
    'R': {'W': 1.0, 'C': 1.0, 'M': 1.0, 'H': 20.26, 'K': 1.0, 'R': 58.28, 'Q': 20.26, 'E': 1.0, 'T': 1.0, 'S': 1.0, 'N': 1.0, 'A': 1.0, 'D': -14.03, 'F': 1.0, 'Y': -6.54, 'V': 1.0, 'I': 1.0, 'P': 20.26, 'L': 1.0, 'G': -7.49},
    'Q': {'W': 1.0, 'C': -1.4, 'M': 1.0, 'H': 1.0, 'K': 1.0, 'R': 1.0, 'Q': 20.26, 'E': 20.26, 'T': 1.0, 'S': 44.94, 'N': 1.0, 'A': 1.0, 'D': 20.26, 'F': -6.54, 'Y': -6.54, 'V': -6.54, 'I': 1.0, 'P': 20.26, 'L': 1.0, 'G': 1.0},
    'E': {'W': -14.03, 'C': 44.94, 'M': 1.0, 'H': -6.54, 'K': 1.0, 'R': 1.0, 'Q': 20.26, 'E': 33.6, 'T': 1.0, 'S': 20.26, 'N': 1.0, 'A': 1.0, 'D': 20.26, 'F': 1.0, 'Y': 1.0, 'V': 1.0, 'I': 20.26, 'P': 20.26, 'L': 1.0, 'G': 1.0},
    'T': {'W': -14.03, 'C': 1.0, 'M': 1.0, 'H': 1.0, 'K': 1.0, 'R': 1.0, 'Q': -6.54, 'E': 20.26, 'T': 1.0, 'S': 1.0, 'N': -14.03, 'A': 1.0, 'D': 1.0, 'F': 13.34, 'Y': 1.0, 'V': 1.0, 'I': 1.0, 'P': 1.0, 'L': 1.0, 'G': -7.49},
    'S': {'W': 1.0, 'C': 33.6, 'M': 28.57, 'H': 1.0, 'K': 1.0, 'R': 20.26, 'Q': 20.26, 'E': 1.0, 'T': 1.0, 'S': 20.26, 'N': 1.0, 'A': 1.0, 'D': 1.0, 'F': 1.0, 'Y': 1.0, 'V': 1.0, 'I': 1.0, 'P': 44.94, 'L': 1.0, 'G': 13.34},
    'N': {'W': -9.37, 'C': -1.4, 'M': 1.0, 'H': 1.0, 'K': 24.68, 'R': 1.0, 'Q': -6.54, 'E': 1.0, 'T': 1.0, 'S': 1.0, 'N': 1.0, 'A': 1.0, 'D': 1.0, 'F': -14.03, 'Y': 1.0, 'V': 1.0, 'I': 44.94, 'P': -1.4, 'L': 1.0, 'G': -6.54},
    'A': {'W': 1.0, 'C': 44.94, 'M': 1.0, 'H': 1.0, 'K': 1.0, 'R': 1.0, 'Q': 1.0, 'E': 1.0, 'T': 1.0, 'S': 1.0, 'N': 1.0, 'A': 1.0, 'D': -7.49, 'F': 1.0, 'Y': 1.0, 'V': 1.0, 'I': 1.0, 'P': 20.26, 'L': 1.0, 'G': 1.0},
    'D': {'W': 1.0, 'C': 1.0, 'M': 1.0, 'H': 1.0, 'K': -7.49, 'R': -6.54, 'Q': 1.0, 'E': 1.0, 'T': -14.03, 'S': 20.26, 'N': 1.0, 'A': 1.0, 'D': 1.0, 'F': -6.54, 'Y': 1.0, 'V': 1.0, 'I': 1.0, 'P': 1.0, 'L': 1.0, 'G': 1.0},
    'F': {'W': 1.0, 'C': 1.0, 'M': 1.0, 'H': 1.0, 'K': 1.0, 'R': 1.0, 'Q': 1.0, 'E': 1.0, 'T': 1.0, 'S': 1.0, 'N': 1.0, 'A': 1.0, 'D': 13.34, 'F': 1.0, 'Y': 33.6, 'V': 1.0, 'I': 1.0, 'P': 20.26, 'L': 1.0, 'G': 1.0},
    'Y': {'W': -14.03, 'C': 1.0, 'M': 44.94, 'H': 13.34, 'K': 1.0, 'R': -14.03, 'Q': 1.0, 'E': -6.54, 'T': -7.49, 'S': 1.0, 'N': 1.0, 'A': 24.68, 'D': 24.68, 'F': 1.0, 'Y': 13.34, 'V': 1.0, 'I': 1.0, 'P': 13.34, 'L': 1.0, 'G': -7.49},
    'V': {'W': 1.0, 'C': 1.0, 'M': 1.0, 'H': 1.0, 'K': -7.49, 'R': 1.0, 'Q': 1.0, 'E': 1.0, 'T': -7.49, 'S': 1.0, 'N': 1.0, 'A': 1.0, 'D': -14.03, 'F': 1.0, 'Y': -6.54, 'V': 1.0, 'I': 1.0, 'P': 20.26, 'L': 1.0, 'G': -7.49},
    'I': {'W': 1.0, 'C': 1.0, 'M': 1.0, 'H': 13.34, 'K': -7.49, 'R': 1.0, 'Q': 1.0, 'E': 44.94, 'T': 1.0, 'S': 1.0, 'N': 1.0, 'A': -7.49, 'D': 1.0, 'F': 1.0, 'Y': 1.0, 'V': -24.33, 'I': 1.0, 'P': -1.4, 'L': 20.26, 'G': -7.49},
    'P': {'W': -1.4, 'C': -6.54, 'M': -6.54, 'H': 1.0, 'K': 1.0, 'R': -6.54, 'Q': 20.26, 'E': 18.38, 'T': 1.0, 'S': 20.26, 'N': 1.0, 'A': 20.26, 'D': -6.54, 'F': 20.26, 'Y': 1.0, 'V': 20.26, 'I': 1.0, 'P': 20.26, 'L': 1.0, 'G': -7.49},
    'L': {'W': 24.68, 'C': 1.0, 'M': 1.0, 'H': 1.0, 'K': -7.49, 'R': 20.26, 'Q': 33.6, 'E': 1.0, 'T': 1.0, 'S': 1.0, 'N': 1.0, 'A': -7.49, 'D': 1.0, 'F': 1.0, 'Y': 1.0, 'V': 1.0, 'I': 1.0, 'P': 20.26, 'L': 1.0, 'G': -7.49},
    'G': {'W': 13.34, 'C': 1.0, 'M': 1.0, 'H': 1.0, 'K': -7.49, 'R': 1.0, 'Q': 1.0, 'E': -6.54, 'T': -7.49, 'S': 1.0, 'N': -7.49, 'A': -7.49, 'D': 1.0, 'F': 1.0, 'Y': -21.64, 'V': 1.0, 'I': -7.49, 'P': 1.0, 'L': 1.0, 'G': 13.34}
}

# Standard Lehninger / Bjellqvist pKa values for Isoelectric point determination
PKA_VALUES = {
    'N_TERM': 9.69,
    'C_TERM': 2.34,
    'K': 10.53,  # Lysine
    'R': 12.48,  # Arginine
    'H': 6.00,   # Histidine
    'D': 3.86,   # Aspartic acid
    'E': 4.25,   # Glutamic acid
    'C': 8.33,   # Cysteine
    'Y': 10.07   # Tyrosine
}


def calculate_net_charge(sequence: str, pH: float) -> float:
    """Calculates net electrical charge of a protein sequence at a given pH."""
    # N-terminal amino group (positive when protonated)
    charge = 1.0 / (1.0 + 10.0 ** (pH - PKA_VALUES['N_TERM']))
    # C-terminal carboxyl group (negative when deprotonated)
    charge -= 1.0 / (1.0 + 10.0 ** (PKA_VALUES['C_TERM'] - pH))

    # Basic residues (positive when protonated)
    k_count = sequence.count('K')
    r_count = sequence.count('R')
    h_count = sequence.count('H')
    charge += k_count * (1.0 / (1.0 + 10.0 ** (pH - PKA_VALUES['K'])))
    charge += r_count * (1.0 / (1.0 + 10.0 ** (pH - PKA_VALUES['R'])))
    charge += h_count * (1.0 / (1.0 + 10.0 ** (pH - PKA_VALUES['H'])))

    # Acidic and neutral-ionizable residues (negative when deprotonated)
    d_count = sequence.count('D')
    e_count = sequence.count('E')
    c_count = sequence.count('C')
    y_count = sequence.count('Y')
    charge -= d_count * (1.0 / (1.0 + 10.0 ** (PKA_VALUES['D'] - pH)))
    charge -= e_count * (1.0 / (1.0 + 10.0 ** (PKA_VALUES['E'] - pH)))
    charge -= c_count * (1.0 / (1.0 + 10.0 ** (PKA_VALUES['C'] - pH)))
    charge -= y_count * (1.0 / (1.0 + 10.0 ** (PKA_VALUES['Y'] - pH)))

    return charge


def compute_isoelectric_point(sequence: str) -> float:
    """
    Computes theoretical isoelectric point (pI) using bisection root-finding method.
    Accuracy: 0.01 pH units. Range: 0.0 to 14.0.
    """
    low_ph = 0.0
    high_ph = 14.0
    precision = 0.001

    while (high_ph - low_ph) > precision:
        mid_ph = (low_ph + high_ph) / 2.0
        charge = calculate_net_charge(sequence, mid_ph)
        if charge > 0:
            low_ph = mid_ph
        else:
            high_ph = mid_ph

    return round((low_ph + high_ph) / 2.0, 2)


def compute_instability_index(sequence: str) -> Tuple[float, str]:
    """
    Computes the instability index (Guruprasad et al., 1990).
    II = (10 / L) * sum_{i=1}^{L-1} DIWV[x_i, x_{i+1}]
    Returns (index_value, classification_string).
    """
    length = len(sequence)
    if length <= 1:
        return 0.0, "Stable"

    score_sum = 0.0
    for i in range(length - 1):
        r1 = sequence[i]
        r2 = sequence[i + 1]
        weight = DIWV_DATA.get(r1, {}).get(r2, 1.0)
        score_sum += weight

    index = round((10.0 / length) * score_sum, 2)
    classification = "Stable (Predicted in vitro half-life > 16 hours)" if index < 40.0 else "Unstable (Predicted in vitro half-life < 5 hours)"
    return index, classification


def compute_aliphatic_index(sequence: str) -> float:
    """
    Computes aliphatic index (Ikai, 1980).
    AI = X(Ala) + a * X(Val) + b * (X(Ile) + X(Leu))
    where a = 2.9, b = 3.9 (relative volume of aliphatic side chains).
    """
    length = len(sequence)
    if length == 0:
        return 0.0

    n_ala = sequence.count('A')
    n_val = sequence.count('V')
    n_ile = sequence.count('I')
    n_leu = sequence.count('L')

    x_ala = (n_ala / length) * 100.0
    x_val = (n_val / length) * 100.0
    x_ile = (n_ile / length) * 100.0
    x_leu = (n_leu / length) * 100.0

    ai = x_ala + (2.9 * x_val) + (3.9 * (x_ile + x_leu))
    return round(ai, 2)


def compute_gravy(sequence: str) -> float:
    """
    Computes Grand Average of Hydropathicity (GRAVY; Kyte & Doolittle, 1982).
    Positive score: overall hydrophobic. Negative score: overall hydrophilic.
    """
    length = len(sequence)
    if length == 0:
        return 0.0

    total_hydropathy = sum(KYTE_DOOLITTLE.get(aa, 0.0) for aa in sequence)
    return round(total_hydropathy / length, 3)


def compute_hydropathy_profile(sequence: str, window_size: int = 9) -> Dict[str, Any]:
    """
    Computes Kyte-Doolittle sliding window hydropathy profile for charting.
    Returns window positions and average hydropathy values.
    """
    length = len(sequence)
    if length < window_size:
        # Fall back to smaller window if sequence is short
        window_size = max(3, length if length % 2 == 1 else length - 1)

    half_window = window_size // 2
    positions: List[int] = []
    scores: List[float] = []

    for i in range(half_window, length - half_window):
        window_seq = sequence[i - half_window: i + half_window + 1]
        score = sum(KYTE_DOOLITTLE.get(aa, 0.0) for aa in window_seq) / window_size
        positions.append(i + 1)  # 1-based index
        scores.append(round(score, 2))

    return {
        'window_size': window_size,
        'positions': positions,
        'scores': scores
    }


def compute_extinction_coefficients(sequence: str, molecular_weight: float) -> Dict[str, Any]:
    """
    Computes molar extinction coefficients at 280 nm (M^-1 cm^-1)
    and Absorbance 0.1% (1 g/L) using Pace et al., 1995:
    - Reduced (no cystines): E = (N_W * 5500) + (N_Y * 1490)
    - With cystine bonds: E = (N_W * 5500) + (N_Y * 1490) + (N_Cys // 2) * 125
    """
    n_trp = sequence.count('W')
    n_tyr = sequence.count('Y')
    n_cys = sequence.count('C')

    coeff_reduced = (n_trp * 5500) + (n_tyr * 1490)
    coeff_cystine = coeff_reduced + ((n_cys // 2) * 125)

    abs_reduced = round(coeff_reduced / molecular_weight, 3) if molecular_weight > 0 else 0.0
    abs_cystine = round(coeff_cystine / molecular_weight, 3) if molecular_weight > 0 else 0.0

    return {
        'extinction_reduced': coeff_reduced,
        'extinction_cystine': coeff_cystine,
        'absorbance_0_1_reduced': abs_reduced,
        'absorbance_0_1_cystine': abs_cystine,
        'tryptophan_count': n_trp,
        'tyrosine_count': n_tyr,
        'cysteine_count': n_cys,
        'disulfide_pairs_possible': n_cys // 2
    }


def compute_secondary_structure_propensities(sequence: str) -> Dict[str, float]:
    """
    Computes average secondary structure conformational propensities using Chou-Fasman parameters.
    Returns normalized fractional propensities for alpha-helix, beta-sheet, and turn/coil.
    """
    length = len(sequence)
    if length == 0:
        return {'helix_propensity': 0.0, 'sheet_propensity': 0.0, 'turn_propensity': 0.0}

    total_helix = sum(CHOU_FASMAN.get(aa, {'a': 1.0})['a'] for aa in sequence)
    total_sheet = sum(CHOU_FASMAN.get(aa, {'b': 1.0})['b'] for aa in sequence)
    total_turn = sum(CHOU_FASMAN.get(aa, {'turn': 1.0})['turn'] for aa in sequence)

    avg_h = total_helix / length
    avg_s = total_sheet / length
    avg_t = total_turn / length

    total_sum = avg_h + avg_s + avg_t
    return {
        'helix_fraction': round((avg_h / total_sum) * 100, 1),
        'sheet_fraction': round((avg_s / total_sum) * 100, 1),
        'turn_fraction': round((avg_t / total_sum) * 100, 1),
        'raw_helix': round(avg_h, 3),
        'raw_sheet': round(avg_s, 3),
        'raw_turn': round(avg_t, 3)
    }


def analyze_physicochemical_properties(sequence: str) -> Dict[str, Any]:
    """
    Comprehensive physicochemical analysis orchestrator.
    Returns full dictionary of calculated sequence metrics.
    """
    length = len(sequence)
    if length == 0:
        return {}

    # 1. Molecular Weight
    mw = sum(RESIDUE_WEIGHTS.get(aa, 110.0) for aa in sequence) + WATER_WEIGHT
    mw_da = round(mw, 2)
    mw_kda = round(mw / 1000.0, 3)

    # 2. Isoelectric Point and Net Charge
    pI = compute_isoelectric_point(sequence)
    charge_ph7 = round(calculate_net_charge(sequence, 7.0), 2)

    # 3. Extinction Coefficients
    extinction = compute_extinction_coefficients(sequence, mw)

    # 4. Stability & Indices
    instability_index, stability_class = compute_instability_index(sequence)
    aliphatic_index = compute_aliphatic_index(sequence)
    gravy = compute_gravy(sequence)

    # 5. Sliding window profile
    hydropathy_profile = compute_hydropathy_profile(sequence)

    # 6. Secondary structure propensities
    sec_struct = compute_secondary_structure_propensities(sequence)

    # 7. Amino Acid Composition Breakdown
    composition = {}
    for aa in sorted(STANDARD_AMINO_ACIDS):
        count = sequence.count(aa)
        pct = round((count / length) * 100.0, 2)
        name, code3 = AA_NAMES[aa]
        composition[aa] = {
            'code1': aa,
            'code3': code3,
            'name': name,
            'count': count,
            'percentage': pct
        }

    # 8. Biochemical Grouping
    aliphatic_aa = sum(sequence.count(aa) for aa in 'AVILM')
    aromatic_aa = sum(sequence.count(aa) for aa in 'FWY')
    polar_uncharged_aa = sum(sequence.count(aa) for aa in 'STNQ')
    positive_charged_aa = sum(sequence.count(aa) for aa in 'KRH')
    negative_charged_aa = sum(sequence.count(aa) for aa in 'DE')
    special_aa = sum(sequence.count(aa) for aa in 'GPC')

    biochemical_groups = {
        'Aliphatic / Hydrophobic (A, V, I, L, M)': {
            'count': aliphatic_aa,
            'percentage': round((aliphatic_aa / length) * 100, 2),
            'color': '#3B82F6'
        },
        'Aromatic (F, W, Y)': {
            'count': aromatic_aa,
            'percentage': round((aromatic_aa / length) * 100, 2),
            'color': '#8B5CF6'
        },
        'Polar Uncharged (S, T, N, Q)': {
            'count': polar_uncharged_aa,
            'percentage': round((polar_uncharged_aa / length) * 100, 2),
            'color': '#10B981'
        },
        'Basic / Positively Charged (K, R, H)': {
            'count': positive_charged_aa,
            'percentage': round((positive_charged_aa / length) * 100, 2),
            'color': '#EF4444'
        },
        'Acidic / Negatively Charged (D, E)': {
            'count': negative_charged_aa,
            'percentage': round((negative_charged_aa / length) * 100, 2),
            'color': '#F59E0B'
        },
        'Special Conformation (G, P, C)': {
            'count': special_aa,
            'percentage': round((special_aa / length) * 100, 2),
            'color': '#6B7280'
        }
    }

    return {
        'sequence_length': length,
        'molecular_weight_da': mw_da,
        'molecular_weight_kda': mw_kda,
        'isoelectric_point': pI,
        'net_charge_ph7': charge_ph7,
        'instability_index': instability_index,
        'stability_class': stability_class,
        'aliphatic_index': aliphatic_index,
        'gravy': gravy,
        'hydropathy_nature': 'Hydrophobic' if gravy > 0 else 'Hydrophilic',
        'extinction': extinction,
        'secondary_structure': sec_struct,
        'composition': composition,
        'biochemical_groups': biochemical_groups,
        'hydropathy_profile': hydropathy_profile
    }
