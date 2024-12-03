#busa_utils/outputs.py
"""
This file contains the functions to save the outputs of the encoding and decoding processes.
"""

def save_outputs(output_path, data):
    """Save data to the output file."""
    with open(output_path, 'w') as f:
        f.write(data)


