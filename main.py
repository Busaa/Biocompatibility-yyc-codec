"""
Python 3.5.6
Ying-Yang-Codec adaptation to be able to make better outputs of sequences for DNA-storage in-vivo.
Arthur Busanello UFRGS 2024
"""
# Configure paths
import sys
sys.path.insert(0, "/home/busa/projetos/Biocompatibility-yyc-codec-main/codec")
sys.path.insert(0, "/home/busa/projetos/Biocompatibility-yyc-codec-main")

# Calling the libraries and functions that we will be using
import os
import pandas as pd
import logging
import math
from yyc import pipeline
from yyc import scheme
from busa_utils.inputs import create_directories, check_input_files, setup_YYC, archive_binary_parting
from busa_utils.outputs import encode_files_in_directory, decode_files_in_directory
from busa_utils.log import setup_logging


# Defining the variables, parameters, and paths
## Inputs archives
input_file = "busa_inputs/texto.txt"  # Example input text file
parts_files_path = "busa_inputs/binary_parts" # Example parts of the input file after binary division
genome_path = "busa_inputs/Bsub-Cohn-genome.fasta"  # Genome file for nullomer analysis
parts_size = 1000  # Size of the parts to be generated from the input file
## Encoding variables
output_directory = "busa_outputs"
dna_parts_dir = "busa_outputs/dna_parts"  # Path to save the DNA parts generated in the encoding process
models_dir = "busa_outputs/models"  # Path to save the models used in the encoding process
segment_length = 120 # Length of the segments of DNA for encoding (incorporation process)
gc_content_range = (30, 70)  # Minimum and maximum GC content for sequences
primer_5 = "GCTTCTGCTTGTCGTCG"  # 5' primer for the encoding process
primer_3 = "GCTTCTGCTTGTCGTCG"  # 3' primer for the encoding process
## Decoding variables
decoded_output_dir = "busa_outputs/decoded_output"  # Path to save the decoded output files
## Output logs variables
log_file = "log.md"
log_dir = "busa_outputs/logs"  # Directory to save logs
log_level = logging.DEBUG # Set to logging.DEBUG for more detailed logs, logging.INFO for less

# CODE

# Create necessary directories
create_directories([output_directory,parts_files_path, dna_parts_dir, log_dir, models_dir, decoded_output_dir])

#  Setup logging
logger = setup_logging(log_level=log_level, log_file=os.path.join(log_dir, log_file))
logger.info("0.0 Initializations ...")

#  Check input files
if not check_input_files([input_file, genome_path], logger):
    logger.error("Missing input files. Exiting...")
    exit(1)

# Parting the input file
archive_binary_parting(input_file, parts_path=parts_files_path, part_size=parts_size, logger=logger)

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
logger.info("Starting the encoding process...")
encode_files_in_directory(
    input_dir=parts_files_path, 
    output_dir=dna_parts_dir, 
    model_dir=models_dir, 
    pipeline=pipeline, 
    tool=tool, 
    segment_length=120, 
    need_index=True, 
    need_log=False,
    verify=None, 
    logger=logger
)
logger.info("Encoding completed successfully.")


# Perfoming test for decoding
logger.info("Starting the decoding process...")
decode_files_in_directory(
    dna_parts_dir=dna_parts_dir, 
    model_dir=models_dir, 
    output_dir=decoded_output_dir, 
    pipeline=pipeline, 
    tool=tool, 
    logger=logger
)
logger.info("Decoding completed successfully.")

# Uniting the decoded files in one and comparing to the original input file
logger.info("Comparing the original input file with the decoded output file...")
decoded_files = os.listdir(decoded_output_dir)
decoded_files.sort()
decoded_data = ""
for file in decoded_files:
    decoded_data += open(os.path.join(decoded_output_dir, file), "r").read()
original_data = open(input_file, "r").read()
if decoded_data == original_data:
    logger.info("The original input file and the decoded output file are the same.")
    print("The original input file and the decoded output file are the same.")