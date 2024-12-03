#busa_utils/biocompatible.py
"""
This file contains the functions used to validate the sequences generated for the DNA storage process, based on
the organism where the sequences will be stored.
"""

def check_nullomers(sequence, observed_sequences):
    """Check if a sequence is a nullomer."""
    return sequence not in observed_sequences
