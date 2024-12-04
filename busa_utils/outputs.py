#busa_utils/outputs.py
"""
This file contains the functions to save the outputs of the encoding and decoding processes.
"""
import os
import pandas as pd
from busa_utils.inputs import create_directories

## Endoncing outputs

def save_encoding_log_and_table(
    input_path, dna_output_path, log_path, output_table_path, 
    dna_dir, segment_length, primer_5, primer_3, logger
):
    """
    Saves a detailed log and a CSV table summarizing the DNA encoding process.

    Args:
        input_path (str): Path to the original input file.
        dna_output_path (str): Path to the encoded DNA output file.
        log_path (str): Path to save the consolidated log file.
        output_table_path (str): Path to save the final CSV table.
        dna_dir (str): Directory to save the final DNA file.
        segment_length (int): Length of each DNA segment.
        primer_5 (str): 5' primer sequence.
        primer_3 (str): 3' primer sequence.
        logger (Logger): Logger instance for logging messages.

    Returns:
        None
    """
    try:
        # Initialize storage for table data
        table_data = []

        # Read original and encoded DNA content
        with open(input_path, "r") as original_file, open(dna_output_path, "r") as dna_file:
            original_content = original_file.read()
            dna_content = dna_file.read()

        # Convert original content to binary
        original_bits = "".join(format(ord(c), "08b") for c in original_content)

        # Write consolidated log
        with open(log_path, "w") as log_file:
            log_file.write("## Consolidated Log for DNA Encoding\n\n")
            log_file.write("### General Information\n")
            log_file.write("- Original file size (bits): {}\n".format(len(original_bits)))
            log_file.write("- Encoded DNA size (nucleotides): {}\n".format(len(dna_content)))
            log_file.write("- Bits per nucleotide: {:.2f}\n\n".format(len(original_bits) / len(dna_content)))
            log_file.write("### Detailed Information\n")

            # Process segments
            segment_count = len(dna_content) // segment_length
            for segment_idx in range(segment_count):
                dna_segment = dna_content[segment_idx * segment_length:(segment_idx + 1) * segment_length]
                bit_start = segment_idx * (len(original_bits) // segment_count)
                bit_end = (segment_idx + 1) * (len(original_bits) // segment_count)
                original_segment_bits = original_bits[bit_start:bit_end]

                # Append to log
                log_file.write("#### Segment {}\n".format(segment_idx + 1))
                log_file.write("- DNA: {}\n".format(dna_segment))
                log_file.write("- Original bits: {}\n".format(original_segment_bits))

                # Process segment details
                index_seq = dna_segment[:16]  # Example: First 16 nucleotides for index
                dna_data_seq = dna_segment[16:-16]  # Exclude primers and index
                full_sequence = "{}{}{}{}".format(primer_5, dna_data_seq, index_seq, primer_3)

                table_data.append({
                    "Index": segment_idx,
                    "Part": "Segment_{}".format(segment_idx + 1),
                    "Nucleotides": len(dna_segment),
                    "Bits": len(original_segment_bits),
                    "Bits/Nucleotide": len(original_segment_bits) / len(dna_segment),
                    "Primer_5": primer_5,
                    "DNA_Data": dna_data_seq,
                    "Index": index_seq,
                    "Primer_3": primer_3,
                    "Full_Sequence": full_sequence,
                    "Original_Bytes": original_content[bit_start // 8:bit_end // 8],
                    "Original_Bits": original_segment_bits
                })

            log_file.write("\n---\n\n")
        

        # Save table to CSV
        df_table = pd.DataFrame(table_data)
        df_table.to_csv(output_table_path, index=False)

        logger.info("Encoding summary saved as a CSV file at {}".format(output_table_path))

    except Exception as e:
        logger.error("Error during log and table generation: {}".format(e))



def encode_files_in_directory(
    input_dir, output_dir, model_dir, pipeline, tool, segment_length=120, need_index=True, need_log=False, verify=None, logger=None
):
    """
    Encodes all files in the input directory and saves the encoded files and models.

    Parameters:
        input_dir (str): Directory containing the input files to be encoded.
        output_dir (str): Directory to save the encoded DNA files.
        model_dir (str): Directory to save the encoding models.
        pipeline: The encoding pipeline module.
        tool: The encoding tool (e.g., YYC instance).
        segment_length (int): Length of each DNA segment.
        need_index (bool): Whether to include indices in the encoding.
        need_log (bool): Whether to generate detailed logs during encoding.
        logger: Logger instance for logging messages.

    Returns:
        List[str]: List of paths to the encoded files.
    """
    if logger:
        logger.info("Starting encoding process for all files in '{}'.".format(input_dir))

    encoded_files = []

    for file_name in sorted(os.listdir(input_dir)):  # Ensure files are processed in order
        input_file_path = os.path.join(input_dir, file_name)

        if os.path.isfile(input_file_path):
            # Define output paths
            encoded_file_path = os.path.join(output_dir, "{}_encoded.dna".format(file_name))
            model_file_path = os.path.join(model_dir, "{}.pkl".format(file_name))

            if logger:
                logger.info("Encoding file: {}".format(file_name))
                logger.debug(
                    "Encoding parameters: input={}, output={}, model={}, segment_length={}, need_index={}, need_log={}".format(
                        input_file_path, encoded_file_path, model_file_path, segment_length, need_index, need_log
                    )
                )

            try:
                # Perform encoding
                pipeline.encode(
                    method=tool,
                    input_path=input_file_path,
                    output_path=encoded_file_path,
                    model_path=model_file_path,
                    need_index=need_index,
                    segment_length=segment_length,
                    verify=verify,
                    need_log=need_log
                )
                encoded_files.append(encoded_file_path)

                if logger:
                    logger.info("Successfully encoded file: {}".format(file_name))
            except Exception as e:
                if logger:
                    logger.error("Error encoding file {}: {}".format(file_name, e))
        else:
            if logger:
                logger.warning("Skipping non-file item: {}".format(file_name))

    return encoded_files



def decode_files_in_directory(
    dna_parts_dir, model_dir, output_dir, pipeline, tool, logger=None
):
    """
    Decodes all DNA files in the input directory using corresponding models.

    Parameters:
        dna_parts_dir (str): Directory containing the encoded DNA files.
        model_dir (str): Directory containing the decoding models.
        output_dir (str): Directory to save the decoded binary files.
        pipeline: The decoding pipeline module.
        tool: The decoding tool (e.g., YYC instance).
        logger: Logger instance for logging messages.

    Returns:
        List[str]: List of paths to the decoded files.
    """
    if logger:
        logger.info("Starting decoding process for all files in '{}'.".format(dna_parts_dir))

    decoded_files = []

    for file_name in sorted(os.listdir(dna_parts_dir)):  # Ensure files are processed in order
        dna_file_path = os.path.join(dna_parts_dir, file_name)
        model_file_path = os.path.join(model_dir, "{}.pkl".format(file_name.replace("_encoded.dna", "")))
        decoded_file_path = os.path.join(output_dir, "{}_decoded.bin".format(file_name.replace("_encoded.dna", "")))

        if os.path.isfile(dna_file_path):
            if logger:
                logger.info("Decoding file: {}".format(file_name))
                logger.debug(
                    "Decoding parameters: dna_file={}, model_file={}, decoded_output={}".format(
                        dna_file_path, model_file_path, decoded_file_path
                    )
                )

            try:
                # Perform decoding
                pipeline.decode(
                    method=tool,
                    input_path=dna_file_path,
                    model_path=model_file_path,
                    output_path=decoded_file_path
                )
                decoded_files.append(decoded_file_path)

                if logger:
                    logger.info("Successfully decoded file: {}".format(file_name))
            except Exception as e:
                if logger:
                    logger.error("Error decoding file {}: {}".format(file_name, e))
        else:
            if logger:
                logger.warning("Skipping non-file item: {}".format(file_name))

    return decoded_files


