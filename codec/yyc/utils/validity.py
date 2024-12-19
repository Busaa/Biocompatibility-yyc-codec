"""
Name: Validity

Coder: HaoLing ZHANG (BGI-Research)[V1]

Function(s):
Check the validity of requested DNA sequence.
The validity describes the friendly index for DNA synthesis, sequencing and related operations
"""

import re
import subprocess
import math


def check(sequence, max_homopolymer=math.inf, max_content=1, min_free_energy=None):
    """
    Check the validity of requested DNA sequence.

    :param sequence: requested DNA sequence.
    :param max_homopolymer: maximum length of homopolymer.
    :param max_content: maximum content of C and G, which means GC content is in [1 - max_content, max_content].
    :param min_free_energy: the free energy of DNA sequence is lower than required min free energy.

    :return: whether the DNA sequence can be considered as valid for DNA synthesis and sequencing.
    """
    if not homopolymer(sequence, max_homopolymer):
        return False
    if not cg_content(sequence, max_content):
        return False
    if  is_palindrome(sequence): #check if the sequence is a palindrome modified by Busaa 2024
        return False
    if is_nullomer(sequence, "busa_inputs/genome.fasta")==False: #check if the sequence is a nullomer modified by Busaa 2024
        return False

    return True


def homopolymer(sequence, max_homopolymer):
    """
    Check the max homopolymer of requested DNA sequence.

    :param sequence: DNA sequence needs detecting.
    :param max_homopolymer: maximum length of homopolymer.

    :return: whether the DNA sequence can be considered as valid for DNA synthesis and sequencing.
    """
    if max_homopolymer > len(sequence):
        return True

    missing_segments = ["A" * (1 + max_homopolymer), "C" * (1 + max_homopolymer), "G" * (1 + max_homopolymer),
                        "T" * (1 + max_homopolymer)] # creates a list of all possible homopolymers in DNA sequence

    for missing_segment in missing_segments:
        if missing_segment in "".join(sequence): # checks if any of the homopolymers are in the DNA sequence 
            return False
    return True


def cg_content(motif, max_content):
    """
    Check the C and G content of requested DNA sequence.

    :param motif: requested DNA sequence.
    :param max_content: maximum content of C and G, which means GC content is in [1 - max_content, max_content].

    :return: whether the DNA sequence can be considered as valid for DNA synthesis and sequencing.
    """
    return (1 - max_content) <= float(motif.count("C") + motif.count("G")) / float(len(motif)) <= max_content


def fold(motif, min_free_energy):
    """
    Call RNAfold to calculate hairpin MFE of a motif

    :param motif: requested DNA sequence.
    :param min_free_energy: min free energy.

    :return: whether the free energy of DNA sequence is lower than required min free energy.
    """
    if min_free_energy is None:
        return True

    process = subprocess.Popen('echo "%s" | RNAfold --noPS --noGU --noconv -T 59.1' % motif,
                               stdout=subprocess.PIPE, shell=True)
    process.wait()
    if process.returncode == 0:
        line = process.stdout.read().decode().split('\n')[1]
        m = re.search("(\S+)\s+\(\s*(\S+)\)", line)
        if m:
            if min_free_energy > float(m.group(2)):
                return True

    return False

## Modifications made to the original code:
### by Busaa 2024
def generate_complement_strand(sequence):
	""" 
	Generate the complementary DNA strand for a given sequence. 
	:param sequence: DNA sequence being codified by the YYC (5' -> 3'). 
	:return: Complementary sequence (3' -> 5'). 
	""" 

	complement_map = {
        "A": "T", 
        "T": "A", 
        "C": "G", 
        "G": "C"
    } 
	return "".join(complement_map[base] for base in sequence)

def is_palindrome(sequence):
    """ 
    Check if a given DNA sequence is a palindrome. 
    :param sequence: DNA sequence being codified by the YYC (5' -> 3'). 
    :return: Whether the DNA sequence is a palindrome. 
    """ 

    return sequence == generate_complement_strand(sequence)[::-1] #invert the complementary sequence and check if it is equal to the original sequence

def is_nullomer(sequence_codec, sequence_ref_path):
    """
    Check if the sequence is a nullomer, i.e., not found in the reference FASTA sequence.

    :param sequence_codec: DNA sequence being codified by the YYC (5' -> 3').
    :param sequence_ref_path: Path to the reference FASTA file.
    :return: True if the sequence is a nullomer (not present in the reference), False otherwise.
    """
    with open(sequence_ref_path, "r") as ref_file:
        ref_lines = ref_file.readlines()

    # Concatenate all sequence lines while ignoring headers (lines starting with '>')
    ref_sequence = "".join(line.strip() for line in ref_lines if not line.startswith(">"))

    # Check if the sequence is absent in the reference
    return sequence_codec not in ref_sequence
