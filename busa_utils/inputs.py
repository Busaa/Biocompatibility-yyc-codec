# busa_utils/inputs.py
"""
This file contains the functions to manage and setup the inputs of the Codec process.
"""

import sys
sys.path.insert(0, "/home/busa/PROJETOS/DNA-STORAGE/codec")
sys.path.insert(0, "/home/busa/PROJETOS/DNA-STORAGE")
import os
import math
from yyc import scheme 


def create_directories(directories):
    """Create directories if they do not exist."""
    for directory in directories:
        if not os.path.exists(directory):
            os.makedirs(directory)

def check_input_files(file_paths, logger):
    """Check if input files exist."""
    for file_path in file_paths:
        if not os.path.exists(file_path):
            logger.error("Missing input file: {}".format(file_path))
            return False
    return True

def setup_YYC(
    base_reference=[0, 1, 0, 1],
    current_code_matrix=[
        [1, 1, 0, 0],
        [1, 0, 0, 1],
        [1, 1, 0, 0],
        [1, 1, 0, 0]
    ],
    support_bases="A",
    support_spacing=0,
    max_ratio=0.8,
    search_count=100,
    max_homopolymer=4,
    max_content=1,
    min_free_energy=None
):
    """
    Initialize the codec, support base, Ying and Yang rules, and return the YYC instance.

    Args:
        base_reference (list): Initial Yang rule (default: [0, 1, 0, 1]).
        current_code_matrix (list): Initial Ying rules (default: predefined matrix).
        support_bases (str): Initial support base (default: "A").
        support_spacing (int): Support spacing (default: 0).
        max_ratio (float): Maximum allowed ratio (default: 0.8).
        search_count (int): Number of searches allowed (default: 100).
        max_homopolymer (int or float): Maximum allowed homopolymer length (default: math.inf).
        max_content (float): Maximum allowed GC content (default: 1).
        min_free_energy (float or None): Minimum free energy (default: None).

    Returns:
        tool (YYC): Initialized YYC codec instance.
    """
    ### Instância do codec YYC
    tool = scheme.YYC(
        base_reference=base_reference,
        current_code_matrix=current_code_matrix,
        support_bases=support_bases,
        support_spacing=support_spacing,
        max_ratio=max_ratio,
        search_count=search_count,
        max_homopolymer=max_homopolymer,
        max_content=max_content,
        min_free_energy=min_free_energy
    )


    return tool


def archive_binary_parting(file_path, parts_path, part_size, logger):
    """
    Splits an input file into binary parts of a fixed size and saves them in a specified directory.

    Args:
        file_path (str): Path to the input file.
        parts_path (str): Directory where the parts will be saved.
        part_size (int): Size of each part in bytes.
        logger (Logger): Logger instance for logging messages.

    Returns:
        None
    """
    # Ensure the input file exists
    if not os.path.exists(file_path):
        logger.error("Input file not found: {}".format(file_path))
        return

    # Create the directory for saving parts if it does not exist
    if not os.path.exists(parts_path):
        os.makedirs(parts_path)
        logger.info("Created directory for parts: {}".format(parts_path))

    try:
        # Open the input file for binary reading
        with open(file_path, 'rb') as file:
            part_number = 0
            while True:
                # Read a fixed-size chunk
                data = file.read(part_size)
                if not data:  # End of file
                    break

                # Generate the part file name
                part_file_path = os.path.join(parts_path, "part_{}.bin".format(part_number))

                # Save the part in the specified directory
                with open(part_file_path, 'wb') as part_file:
                    part_file.write(data)

                logger.info("Saved part {} to {}".format(part_number, part_file_path))
                part_number += 1

        logger.info("Successfully split the file into {} parts.".format(part_number))

    except Exception as e:
        logger.error("An error occurred during file splitting: {}".format(e))