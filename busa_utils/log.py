# busa_utils/log.py
"""
This file contains the functions to configure and manage logs for the Codec process with Markdown output support.
"""

import logging
import os

def setup_logging(log_level=logging.INFO, log_file=None):
    """
    Set up logging for the project with optional Markdown log file output. 
    The looger will keep writing in the log file if is the same, it wont overwrite it.

    Args:
        log_level (int): The logging level (e.g., logging.DEBUG, logging.INFO, logging.WARNING).
                         Sets the minimal level of log messages to be shown.
        log_file (str, optional): If provided, logs will be saved to this Markdown file.

    Returns:
        logger: Configured logger object for the project.
    """
    # Create a logger for the project
    logger = logging.getLogger("codec_logger")
    logger.setLevel(log_level)

    # Define the log message format in Markdown
    markdown_formatter = logging.Formatter(
        "**[%(asctime)s]** - `%(name)s` - **%(levelname)s**\n> %(message)s\n"
    )

    if not log_file:
        # Add console handler if no log file is provided
        console_handler = logging.StreamHandler()
        console_handler.setLevel(log_level)
        console_handler.setFormatter(markdown_formatter)
        logger.addHandler(console_handler)
    else:
        # Ensure the file exists or create it
        if not os.path.exists(log_file):
            with open(log_file, 'w') as f:
                f.write("# Log File\n\n")  # Add a Markdown title to the log file

        # Add a file handler to write logs to the Markdown file
        file_handler = logging.FileHandler(log_file, mode='a')
        file_handler.setLevel(log_level)
        file_handler.setFormatter(markdown_formatter)
        logger.addHandler(file_handler)

    return logger
