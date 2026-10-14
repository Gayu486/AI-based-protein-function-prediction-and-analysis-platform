"""
Protein Sequence Validator Module
Validates raw protein sequences and FASTA format input.
Detects invalid characters, ambiguous amino acids, potential nucleotide inputs,
and handles whitespace and formatting.
"""

import re
from typing import Dict, Any, List, Tuple

# Standard 20 proteinogenic amino acids
STANDARD_AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")

# Ambiguous and non-standard amino acid codes (IUPAC)
AMBIGUOUS_AMINO_ACIDS = {
    'B': 'Aspartic acid or Asparagine (Asx)',
    'Z': 'Glutamic acid or Glutamine (Glx)',
    'J': 'Leucine or Isoleucine (Xle)',
    'X': 'Unspecified / unknown amino acid',
    'U': 'Selenocysteine (Sec)',
    'O': 'Pyrrolysine (Pyl)'
}

# Amino acid full names and single-letter codes
AA_NAMES = {
    'A': ('Alanine', 'Ala'),
    'C': ('Cysteine', 'Cys'),
    'D': ('Aspartic acid', 'Asp'),
    'E': ('Glutamic acid', 'Glu'),
    'F': ('Phenylalanine', 'Phe'),
    'G': ('Glycine', 'Gly'),
    'H': ('Histidine', 'His'),
    'I': ('Isoleucine', 'Ile'),
    'K': ('Lysine', 'Lys'),
    'L': ('Leucine', 'Leu'),
    'M': ('Methionine', 'Met'),
    'N': ('Asparagine', 'Asn'),
    'P': ('Proline', 'Pro'),
    'Q': ('Glutamine', 'Gln'),
    'R': ('Arginine', 'Arg'),
    'S': ('Serine', 'Ser'),
    'T': ('Threonine', 'Thr'),
    'V': ('Valine', 'Val'),
    'W': ('Tryptophan', 'Trp'),
    'Y': ('Tyrosine', 'Tyr'),
}

MIN_SEQUENCE_LENGTH = 5
MAX_SEQUENCE_LENGTH = 10000


def parse_fasta(raw_input: str) -> Tuple[str, str]:
    """
    Parses input text which may contain FASTA header lines.
    Returns (header, sequence_str).
    """
    lines = raw_input.strip().splitlines()
    header = ""
    sequence_parts = []

    for line in lines:
        stripped = line.strip()
        if not stripped:
            continue
        if stripped.startswith(">"):
            if not header:
                header = stripped[1:].strip()
        else:
            sequence_parts.append(stripped)

    return header, "".join(sequence_parts)


def validate_protein_sequence(raw_input: str) -> Dict[str, Any]:
    """
    Validates protein sequence string.
    
    Returns dictionary with:
    - is_valid: bool
    - cleaned_sequence: str (uppercase, no whitespace)
    - header: str (if FASTA header present)
    - length: int
    - has_ambiguous: bool
    - ambiguous_details: list of dicts
    - errors: list of str
    - warnings: list of str
    - is_likely_nucleotide: bool
    """
    errors: List[str] = []
    warnings: List[str] = []

    if not raw_input or not raw_input.strip():
        return {
            'is_valid': False,
            'cleaned_sequence': '',
            'header': '',
            'length': 0,
            'has_ambiguous': False,
            'ambiguous_details': [],
            'errors': ['Please enter a protein sequence or upload a FASTA file.'],
            'warnings': [],
            'is_likely_nucleotide': False
        }

    header, raw_seq = parse_fasta(raw_input)

    # Remove all whitespace, linebreaks, digits and common delimiters
    # First check for invalid characters and their exact positions in the raw sequence
    cleaned_seq_chars = []
    invalid_chars_found = []

    # Track 1-based index in the sequence (excluding whitespace and digits)
    char_pos = 0
    for idx, char in enumerate(raw_seq):
        if char.isspace():
            continue
        if char.isdigit():
            # Many sequence copies have residue numbers (e.g. 10 20 30)
            continue
        if char in ('-', '*', '.'):
            # Alignment gap or stop codon
            warnings.append(f"Ignored gap/stop symbol '{char}' at raw character position {idx + 1}.")
            continue

        char_upper = char.upper()
        char_pos += 1

        if char_upper in STANDARD_AMINO_ACIDS or char_upper in AMBIGUOUS_AMINO_ACIDS:
            cleaned_seq_chars.append(char_upper)
        else:
            invalid_chars_found.append((char, char_pos))

    cleaned_sequence = "".join(cleaned_seq_chars)
    length = len(cleaned_sequence)

    # Check for invalid characters
    if invalid_chars_found:
        unique_invalid = sorted(list(set(c for c, _ in invalid_chars_found)))
        first_few = [f"'{c}' at position {pos}" for c, pos in invalid_chars_found[:5]]
        err_msg = f"Invalid character(s) found in sequence: {', '.join(unique_invalid)}. "
        err_msg += f"Examples: {', '.join(first_few)}"
        if len(invalid_chars_found) > 5:
            err_msg += f" (and {len(invalid_chars_found) - 5} more)."
        err_msg += " Only standard amino acids (A, C, D, E, F, G, H, I, K, L, M, N, P, Q, R, S, T, V, W, Y) are accepted."
        errors.append(err_msg)

    # Check length
    if length == 0 and not errors:
        errors.append("No valid amino acid residues found in the provided input.")
    elif length < MIN_SEQUENCE_LENGTH and not errors:
        errors.append(f"Sequence length is {length} AA. Minimum required length is {MIN_SEQUENCE_LENGTH} amino acids for meaningful biological analysis.")
    elif length > MAX_SEQUENCE_LENGTH:
        errors.append(f"Sequence length ({length} AA) exceeds maximum permitted limit of {MAX_SEQUENCE_LENGTH} amino acids.")

    # Check if input is likely DNA or RNA nucleotides (composed only of A, C, G, T/U, N)
    is_likely_nucleotide = False
    if length >= 15:
        nucleotide_chars = set("ACGTUN")
        seq_chars = set(cleaned_sequence)
        if seq_chars.issubset(nucleotide_chars):
            is_likely_nucleotide = True
            warnings.append(
                "Notice: Sequence contains exclusively nucleotide characters (A, C, G, T/U). "
                "If this is a nucleotide (DNA/RNA) sequence, please translate it into an amino acid sequence before analysis."
            )

    # Check for ambiguous amino acids
    ambiguous_details = []
    has_ambiguous = False
    for code, desc in AMBIGUOUS_AMINO_ACIDS.items():
        count = cleaned_sequence.count(code)
        if count > 0:
            has_ambiguous = True
            ambiguous_details.append({
                'code': code,
                'description': desc,
                'count': count,
                'percentage': round((count / length) * 100, 2) if length > 0 else 0
            })
            warnings.append(f"Sequence contains {count} ambiguous residue(s) '{code}' ({desc}). Standard calculations will handle these with standard non-hydrophobic approximations.")

    is_valid = len(errors) == 0

    return {
        'is_valid': is_valid,
        'cleaned_sequence': cleaned_sequence,
        'header': header,
        'length': length,
        'has_ambiguous': has_ambiguous,
        'ambiguous_details': ambiguous_details,
        'errors': errors,
        'warnings': warnings,
        'is_likely_nucleotide': is_likely_nucleotide
    }
