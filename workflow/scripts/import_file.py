from pathlib import Path
import time
import shutil
import gzip
import bz2


def magic_decompress(filename):
    """
    Decompresses a file if it is compressed in gzip or bzip format, or simply opens it if it is a plaintext file.

    Args:
        filename (str): The path to the file to be opened.
    Returns:
        An open binary file handle.
    """

    # Define a dictionary of magic bytes for identifying compressed file formats
    magic_dict = {
        b"\x1f\x8b\x08": (gzip.open, 'rb'),
        b"\x42\x5a\x68": (bz2.BZ2File, 'rb'),
    }

    # Determine the maximum length of any magic byte sequence
    max_len = max(len(x) for x in magic_dict)

    # Read the first `max_len` bytes of the file
    with open(filename, 'rb') as f:
        file_start = f.read(max_len)

    # Check if the file starts with any of the known magic byte sequences
    for magic, (fn, flag) in magic_dict.items():
        if file_start.startswith(magic):
            # If the file is compressed, use the appropriate decompression function
            return fn(filename, flag)

    # If the file is not compressed, simply open it as a plaintext file
    return open(filename, 'rb')

def copycat(out_filepath, input_filepaths):
    """
    Concatenates compressed and plaintext files into a new file.

    Args:
        out_filepath (str): The path to the output file.
        input_filepaths (list): A list of paths to the input files.

    Returns:
        None
    """
    # Get the total number of input files and initialize a counter for processed files
    file_count = len(input_filepaths)
    actual_count = 0

    # Display a progress message indicating the number of input files and the start time
    print(f'{time.strftime("%H:%M:%S", time.localtime())} Creating {out_filepath} from {file_count} source files. This may take some time.')

    # Create a gzip file handle for the output file
    with gzip.open(out_filepath, 'wb') as out_file:
        # Iterate through the input files
        for input_filepath in input_filepaths:
            # Increment the file count and display a progress message
            actual_count += 1
            print(f'\x1b[1K\r{time.strftime("%H:%M:%S", time.localtime())} {actual_count}/{file_count} Loaded file {input_filepath}',end='')  # end='\x1b[1K\r'overwrite instead of new line
            # Open the input file using the `magic_decompress` function to handle various compression formats
            with magic_decompress(input_filepath) as input_file:
                # Copy the contents of the input file to the output file
                shutil.copyfileobj(input_file, out_file)
    # Display a message indicating the completion of the file concatenation
    print(f'{time.strftime("%H:%M:%S", time.localtime())} {out_filepath} was imported.')

if __name__ == "__main__":
    copycat(snakemake.output[0], list(snakemake.input))