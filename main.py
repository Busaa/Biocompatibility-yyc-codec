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
parts_size = 5000  # Size of the parts to be generated from the input file
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
log_level = logging.INFO # Set to logging.DEBUG for more detailed logs, logging.INFO for less

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

# Unindo os arquivos decodificados em um único arquivo sem duplicação
logger.info("Unindo os arquivos decodificados e salvando no arquivo 'decoded_info.txt'...")

# Caminho para o arquivo combinado
decoded_file_path = os.path.join(output_directory, "decoded_info.txt")

# Lê e combina os arquivos decodificados em ordem e salva em 'decoded_info.txt'
with open(decoded_file_path, "w") as decoded_file:
    for file in sorted(os.listdir(decoded_output_dir)):
        file_path = os.path.join(decoded_output_dir, file)
        if os.path.isfile(file_path):  # Garante que apenas arquivos sejam processados
            with open(file_path, "r") as f:
                data = f.read()
                decoded_file.write(data + "\n")  # Adiciona quebra de linha para separar conteúdos

logger.info("Arquivos decodificados foram combinados e salvos em '{}'.".format(decoded_file_path))
print("Arquivos decodificados foram combinados e salvos em '{}'.".format(decoded_file_path))

# Comparando o arquivo combinado com o original
with open(input_file, "r") as f:
    original_data = f.read()

with open(decoded_file_path, "r") as f:
    decoded_data = f.read()

# Comparando os dados originais com os decodificados
similarity = sum(1 for a, b in zip(decoded_data, original_data) if a == b) / max(len(decoded_data), len(original_data)) * 100

if similarity == 100:
    logger.info("Os arquivos original e decodificado são idênticos (100% de similaridade).")
    print("Os arquivos original e decodificado são idênticos (100% de similaridade).")
else:
    logger.info("Os arquivos original e decodificado têm {:.2f}% de similaridade.".format(similarity))
    print("Os arquivos original e decodificado têm {:.2f}% de similaridade.".format(similarity))


