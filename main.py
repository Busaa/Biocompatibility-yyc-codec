"""
Python 3.5.6
Ying-Yang-Codec adaptation to be able to make better outputs of sequences for DNA-storage in-vivo.
Arthur Busanello UFRGS 2024
"""
# Configure paths
import sys
sys.path.insert(0, "/home/busa/PROJETOS/DNA-STORAGE/codec")
sys.path.insert(0, "/home/busa/PROJETOS/DNA-STORAGE")

# Calling the libraries and functions that we will be using
import os
import pandas as pd
import logging
import math
from yyc import pipeline
from yyc import scheme
from busa_utils.inputs import create_directories, check_input_files, setup_YYC, archive_binary_parting
from busa_utils.outputs import save_outputs
from busa_utils.log import setup_logging

# Defining the variables, parameters, and paths
## Inputs archives
input_file = "busa_inputs/texto.txt"  # Example input text file
parts_files_path = "busa_inputs/binary_parts" # Example parts of the input file after binary division
genome_path = "busa_inputs/Bsub-Cohn-genome.fasta"  # Genome file for nullomer analysis
## Encoding variables
output_directory = "busa_outputs/"
dna_parts_path = os.path.join(output_directory, "dna_parts")
segment_length = 120 # Length of the segments of DNA for encoding (incorporation process)
gc_content_range = (30, 70)  # Minimum and maximum GC content for sequences
## Decoding variables
decoded_output_file = os.path.join(output_directory, "decoded_output.txt")
## Output logs variables
log_file = "log.md"
log_level = logging.INFO # Set to logging.DEBUG for more detailed logs, logging.INFO for less

# CODE

# Create necessary directories
create_directories([output_directory, os.path.join(output_directory, "logs"),parts_files_path, dna_parts_path])

#  Setup logging
logger = setup_logging(log_level=log_level, log_file=output_directory+"logs"+log_file)
logger.info("0.0 Initializations ...")

#  Check input files
if not check_input_files([input_file, genome_path], logger):
    logger.error("Missing input files. Exiting...")
    exit(1)

# Parting the input file
archive_binary_parting(input_file, parts_path=parts_files_path, part_size=1200, logger=logger)

# Creating the encoding scheme
logger.info(" Initializing the YYC scheme ...")
tool = setup_YYC(
    base_reference=[0, 1, 0, 1], # Yang Rule. 2 bases must be 1 and other 2 must be 0 
    current_code_matrix=[
        [1, 1, 0, 0],
        [1, 0, 0, 1],
        [1, 1, 0, 0],
        [1, 1, 0, 0]
    ], # Ying Rule. Y (4) = previous base, X (2) = current base pair possibility acording to Yang Rule 
    support_bases="A", # base before oficial data? 
    support_spacing=0, #sapce between support base and current base
    max_ratio=0.8, # max ratio of 0:1 in the binary file, decide good and bad data
    search_count=100, #search count for best enconding possibility in the incporporation process
    max_homopolymer=4, # max homopolymer length
    max_content=0.6, # max GC content
    min_free_energy=None # min free energy of the DNA sequence in encoding
)


# Perform encoding
for part in sorted(os.listdir(parts_files_path)):  # Sorting to ensure parts are processed in order
    input_codec = os.path.join(parts_files_path, part)

    # Check if the item is a file (not a directory)
    if os.path.isfile(input_codec):
        logger.info("Starting the encoding process for part: {}...".format(part))

        # Construct output path for the encoded file
        encoded_output_path = os.path.join(dna_parts_path, "{}_encoded".format(part))

        # Perform encoding
        encoded_sequences = pipeline.encode(
            method=tool, 
            input_path=input_codec, 
            output_path=encoded_output_path,
            model_path=None, 
            verify=None, 
            need_index=True, 
            segment_length=120, 
            need_log=False
        )

        logger.info("Encoding completed successfully for part: {}.".format(part))
    else:
        logger.warning("Skipping non-file item in parts directory:{}".format(part))

"""
# 5. Perform decoding (if applicable)
logger.info("Starting the decoding process...")
pipeline.decode(
    method=None, 
    model_path=None, 
    input_path=None, 
    output_path=None,
    verify=None, 
    has_index=True, 
    need_log=False
)
logger.info("Decoding completed successfully.")

"""