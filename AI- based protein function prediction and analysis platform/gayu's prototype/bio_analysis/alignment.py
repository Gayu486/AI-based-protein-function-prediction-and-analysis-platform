"""
Sequence Alignment and Similarity Scoring Module
Computes real pairwise alignment metrics between query protein sequence
and candidate database sequences.
Calculates:
- Exact sequence identity percentage
- Conservative similarity percentage (BLOSUM62)
- Alignment coverage percentage
- Alignment score and Karlin-Altschul statistical E-value approximation
- Formatted alignment visualization string with 1-based coordinates
"""

import math
from typing import Dict, Any, Tuple, List

# Standard BLOSUM62 Substitution Matrix
BLOSUM62 = {
    ('A', 'A'): 4, ('A', 'R'): -1, ('A', 'N'): -2, ('A', 'D'): -2, ('A', 'C'): 0, ('A', 'Q'): -1, ('A', 'E'): -1, ('A', 'G'): 0, ('A', 'H'): -2, ('A', 'I'): -1, ('A', 'L'): -1, ('A', 'K'): -1, ('A', 'M'): -1, ('A', 'F'): -2, ('A', 'P'): -1, ('A', 'S'): 1, ('A', 'T'): 0, ('A', 'W'): -3, ('A', 'Y'): -2, ('A', 'V'): 0,
    ('R', 'R'): 5, ('R', 'N'): 0, ('R', 'D'): -2, ('R', 'C'): -3, ('R', 'Q'): 1, ('R', 'E'): 0, ('R', 'G'): -2, ('R', 'H'): 0, ('R', 'I'): -3, ('R', 'L'): -2, ('R', 'K'): 2, ('R', 'M'): -1, ('R', 'F'): -3, ('R', 'P'): -2, ('R', 'S'): -1, ('R', 'T'): -1, ('R', 'W'): -3, ('R', 'Y'): -2, ('R', 'V'): -3,
    ('N', 'N'): 6, ('N', 'D'): 1, ('N', 'C'): -3, ('N', 'Q'): 0, ('N', 'E'): 0, ('N', 'G'): 0, ('N', 'H'): 1, ('N', 'I'): -3, ('N', 'L'): -3, ('N', 'K'): 0, ('N', 'M'): -2, ('N', 'F'): -3, ('N', 'P'): -2, ('N', 'S'): 1, ('N', 'T'): 0, ('N', 'W'): -4, ('N', 'Y'): -2, ('N', 'V'): -3,
    ('D', 'D'): 6, ('D', 'C'): -3, ('D', 'Q'): 0, ('D', 'E'): 2, ('D', 'G'): -1, ('D', 'H'): -1, ('D', 'I'): -3, ('D', 'L'): -4, ('D', 'K'): -1, ('D', 'M'): -3, ('D', 'F'): -3, ('D', 'P'): -1, ('D', 'S'): 0, ('D', 'T'): -1, ('D', 'W'): -4, ('D', 'Y'): -3, ('D', 'V'): -3,
    ('C', 'C'): 9, ('C', 'Q'): -3, ('C', 'E'): -4, ('C', 'G'): -3, ('C', 'H'): -3, ('C', 'I'): -1, ('C', 'L'): -1, ('C', 'K'): -3, ('C', 'M'): -1, ('C', 'F'): -2, ('C', 'P'): -3, ('C', 'S'): -1, ('C', 'T'): -1, ('C', 'W'): -2, ('C', 'Y'): -2, ('C', 'V'): -1,
    ('Q', 'Q'): 5, ('Q', 'E'): 2, ('Q', 'G'): -2, ('Q', 'H'): 0, ('Q', 'I'): -3, ('Q', 'L'): -2, ('Q', 'K'): 1, ('Q', 'M'): 0, ('Q', 'F'): -3, ('Q', 'P'): -1, ('Q', 'S'): 0, ('Q', 'T'): -1, ('Q', 'W'): -2, ('Q', 'Y'): -1, ('Q', 'V'): -2,
    ('E', 'E'): 5, ('E', 'G'): -2, ('E', 'H'): 0, ('E', 'I'): -3, ('E', 'L'): -3, ('E', 'K'): 1, ('E', 'M'): -2, ('E', 'F'): -3, ('E', 'P'): -1, ('E', 'S'): 0, ('E', 'T'): -1, ('E', 'W'): -3, ('E', 'Y'): -2, ('E', 'V'): -2,
    ('G', 'G'): 6, ('G', 'H'): -2, ('G', 'I'): -4, ('G', 'L'): -4, ('G', 'K'): -2, ('G', 'M'): -3, ('G', 'F'): -3, ('G', 'P'): -2, ('G', 'S'): 0, ('G', 'T'): -2, ('G', 'W'): -2, ('G', 'Y'): -3, ('G', 'V'): -3,
    ('H', 'H'): 8, ('H', 'I'): -3, ('H', 'L'): -3, ('H', 'K'): -1, ('H', 'M'): -2, ('H', 'F'): -1, ('H', 'P'): -2, ('H', 'S'): -1, ('H', 'T'): -2, ('H', 'W'): -2, ('H', 'Y'): 2, ('H', 'V'): -3,
    ('I', 'I'): 4, ('I', 'L'): 2, ('I', 'K'): -3, ('I', 'M'): 1, ('I', 'F'): 0, ('I', 'P'): -3, ('I', 'S'): -2, ('I', 'T'): -1, ('I', 'W'): -3, ('I', 'Y'): -1, ('I', 'V'): 3,
    ('L', 'L'): 4, ('L', 'K'): -2, ('L', 'M'): 2, ('L', 'F'): 0, ('L', 'P'): -3, ('L', 'S'): -2, ('L', 'T'): -1, ('L', 'W'): -2, ('L', 'Y'): -1, ('L', 'V'): 1,
    ('K', 'K'): 5, ('K', 'M'): -1, ('K', 'F'): -3, ('K', 'P'): -1, ('K', 'S'): 0, ('K', 'T'): -1, ('K', 'W'): -3, ('K', 'Y'): -2, ('K', 'V'): -2,
    ('M', 'M'): 5, ('M', 'F'): 0, ('M', 'P'): -2, ('M', 'S'): -1, ('M', 'T'): -1, ('M', 'W'): -1, ('M', 'Y'): -1, ('M', 'V'): 1,
    ('F', 'F'): 6, ('F', 'P'): -4, ('F', 'S'): -2, ('F', 'T'): -2, ('F', 'W'): 1, ('F', 'Y'): 3, ('F', 'V'): -1,
    ('P', 'P'): 7, ('P', 'S'): -1, ('P', 'T'): -1, ('P', 'W'): -4, ('P', 'Y'): -3, ('P', 'V'): -2,
    ('S', 'S'): 4, ('S', 'T'): 1, ('S', 'W'): -3, ('S', 'Y'): -2, ('S', 'V'): -2,
    ('T', 'T'): 5, ('T', 'W'): -2, ('T', 'Y'): -2, ('T', 'V'): 0,
    ('W', 'W'): 11, ('W', 'Y'): 2, ('W', 'V'): -3,
    ('Y', 'Y'): 7, ('Y', 'V'): -1,
    ('V', 'V'): 4
}


def get_blosum_score(a: str, b: str) -> int:
    """Returns BLOSUM62 score for a pair of amino acids."""
    pair = (a, b)
    if pair in BLOSUM62:
        return BLOSUM62[pair]
    rev_pair = (b, a)
    if rev_pair in BLOSUM62:
        return BLOSUM62[rev_pair]
    if a == b:
        return 4
    return -1


def compute_karlin_altschul_evalue(raw_score: float, query_len: int, db_size_letters: int = 210000000) -> float:
    """
    Computes statistical E-value approximation using Karlin-Altschul statistics.
    E = K * m * n * exp(-lambda * S)
    For BLOSUM62: lambda ≈ 0.267, K ≈ 0.041.
    """
    if raw_score <= 0:
        return 10.0

    lambda_param = 0.267
    k_param = 0.041
    m = max(1, query_len)
    n = db_size_letters

    exponent = -lambda_param * raw_score
    if exponent > 700:
        return 10.0
    if exponent < -700:
        return 0.0

    e_val = k_param * m * n * math.exp(exponent)
    return max(0.0, e_val)


def align_sequences(query: str, subject: str) -> Dict[str, Any]:
    """
    Performs dynamic programming pairwise global/semi-global alignment between query and subject.
    Optimized for fast responsiveness in web application.
    """
    q_len = len(query)
    s_len = len(subject)

    if q_len == 0 or s_len == 0:
        return {
            'identical_count': 0,
            'similar_count': 0,
            'identity_pct': 0.0,
            'similarity_pct': 0.0,
            'coverage_pct': 0.0,
            'e_value': 10.0,
            'e_value_str': '> 10',
            'alignment_score': 0,
            'alignment_length': 0,
            'aligned_query': '',
            'match_line': '',
            'aligned_subject': '',
            'blocks': []
        }

    # If sequences are identical, compute directly
    if query == subject:
        raw_score = sum(get_blosum_score(aa, aa) for aa in query)
        e_val = compute_karlin_altschul_evalue(raw_score, q_len)
        e_val_str = "0.0" if e_val < 1e-180 else f"{e_val:.2e}" if e_val < 0.01 else f"{e_val:.3f}"
        
        match_line = "|" * q_len
        return {
            'identical_count': q_len,
            'similar_count': q_len,
            'identity_pct': 100.0,
            'similarity_pct': 100.0,
            'coverage_pct': 100.0,
            'e_value': e_val,
            'e_value_str': e_val_str,
            'alignment_score': raw_score,
            'alignment_length': q_len,
            'aligned_query': query,
            'match_line': match_line,
            'aligned_subject': subject,
            'blocks': format_alignment_blocks(query, match_line, subject)
        }

    # Dynamic programming alignment with affine gap penalty
    gap_open = -10
    gap_extend = -1

    # Matrix initialization
    dp = [[0] * (s_len + 1) for _ in range(q_len + 1)]
    trace = [[0] * (s_len + 1) for _ in range(q_len + 1)]  # 1: diag, 2: up (gap in subject), 3: left (gap in query)

    for i in range(1, q_len + 1):
        dp[i][0] = gap_open + (i - 1) * gap_extend
        trace[i][0] = 2
    for j in range(1, s_len + 1):
        dp[0][j] = gap_open + (j - 1) * gap_extend
        trace[0][j] = 3

    for i in range(1, q_len + 1):
        q_aa = query[i - 1]
        for j in range(1, s_len + 1):
            s_aa = subject[j - 1]
            match_score = get_blosum_score(q_aa, s_aa)

            diag = dp[i - 1][j - 1] + match_score
            up = dp[i - 1][j] + (gap_extend if trace[i - 1][j] == 2 else gap_open)
            left = dp[i][j - 1] + (gap_extend if trace[i][j - 1] == 3 else gap_open)

            best = max(diag, up, left)
            dp[i][j] = best

            if best == diag:
                trace[i][j] = 1
            elif best == up:
                trace[i][j] = 2
            else:
                trace[i][j] = 3

    # Traceback
    aligned_q = []
    aligned_s = []
    i, j = q_len, s_len

    while i > 0 or j > 0:
        if i > 0 and j > 0 and trace[i][j] == 1:
            aligned_q.append(query[i - 1])
            aligned_s.append(subject[j - 1])
            i -= 1
            j -= 1
        elif i > 0 and (j == 0 or trace[i][j] == 2):
            aligned_q.append(query[i - 1])
            aligned_s.append('-')
            i -= 1
        else:
            aligned_q.append('-')
            aligned_s.append(subject[j - 1])
            j -= 1

    aligned_q.reverse()
    aligned_s.reverse()

    aligned_q_str = "".join(aligned_q)
    aligned_s_str = "".join(aligned_s)
    total_aligned = len(aligned_q_str)

    identical = 0
    similar = 0
    match_chars = []

    for a, b in zip(aligned_q_str, aligned_s_str):
        if a == b and a != '-':
            identical += 1
            similar += 1
            match_chars.append('|')
        elif a != '-' and b != '-':
            blosum = get_blosum_score(a, b)
            if blosum > 0:
                similar += 1
                match_chars.append('+')
            else:
                match_chars.append(' ')
        else:
            match_chars.append(' ')

    match_line_str = "".join(match_chars)
    raw_score = dp[q_len][s_len]

    # Metrics
    identity_pct = round((identical / max(1, total_aligned)) * 100.0, 1)
    similarity_pct = round((similar / max(1, total_aligned)) * 100.0, 1)
    coverage_pct = round((min(q_len, s_len - aligned_q_str.count('-')) / max(1, q_len)) * 100.0, 1)

    e_val = compute_karlin_altschul_evalue(raw_score, q_len)
    if e_val < 1e-180:
        e_val_str = "0.0"
    elif e_val < 0.01:
        e_val_str = f"{e_val:.2e}"
    else:
        e_val_str = f"{e_val:.3f}"

    return {
        'identical_count': identical,
        'similar_count': similar,
        'identity_pct': identity_pct,
        'similarity_pct': similarity_pct,
        'coverage_pct': coverage_pct,
        'e_value': e_val,
        'e_value_str': e_val_str,
        'alignment_score': raw_score,
        'alignment_length': total_aligned,
        'aligned_query': aligned_q_str,
        'match_line': match_line_str,
        'aligned_subject': aligned_s_str,
        'blocks': format_alignment_blocks(aligned_q_str, match_line_str, aligned_s_str)
    }


def format_alignment_blocks(aligned_q: str, match_line: str, aligned_s: str, block_size: int = 60) -> List[Dict[str, Any]]:
    """
    Formats alignment into 60-residue visual blocks with start/end residue numbers.
    """
    blocks = []
    total_len = len(aligned_q)

    q_pos = 1
    s_pos = 1

    for start in range(0, total_len, block_size):
        end = min(start + block_size, total_len)

        chunk_q = aligned_q[start:end]
        chunk_m = match_line[start:end]
        chunk_s = aligned_s[start:end]

        # Calculate coordinates
        q_residues_in_chunk = sum(1 for c in chunk_q if c != '-')
        s_residues_in_chunk = sum(1 for c in chunk_s if c != '-')

        q_start = q_pos if q_residues_in_chunk > 0 else q_pos
        q_end = q_pos + q_residues_in_chunk - 1 if q_residues_in_chunk > 0 else q_pos

        s_start = s_pos if s_residues_in_chunk > 0 else s_pos
        s_end = s_pos + s_residues_in_chunk - 1 if s_residues_in_chunk > 0 else s_pos

        blocks.append({
            'query_str': chunk_q,
            'match_str': chunk_m,
            'subject_str': chunk_s,
            'query_start': q_start,
            'query_end': q_end,
            'subject_start': s_start,
            'subject_end': s_end
        })

        q_pos += q_residues_in_chunk
        s_pos += s_residues_in_chunk

    return blocks
