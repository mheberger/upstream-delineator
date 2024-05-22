"""
Learning how to write a Python script which takes command line arguments.

Run it from the command line like this.
$ python test_args.py lower.txt upper.txt

Alternatively, the main function can be run from another python script:

>> from test_args import process_data
>> process_data(input_file, output_file)

This way, the user has two options for how to run the file, either
from the command line or via another Python script.

"""

import argparse


def process_data(input_file, output_file):
    # Example function to demonstrate processing
    with open(input_file, 'r') as infile:
        data = infile.read()

    # Perform some data processing (this is just an example)
    processed_data = data.upper()  # For instance, converting the text to uppercase

    with open(output_file, 'w') as outfile:
        outfile.write(processed_data)


def main():
    # Create the parser
    parser = argparse.ArgumentParser(description="Process some data from input file and write to output file.")

    # Add the arguments
    parser.add_argument('input', help="Input filename")
    parser.add_argument('output', help="Output filename")

    # Parse the arguments
    args = parser.parse_args()

    # Call the process_data function with the provided filenames
    process_data(args.input, args.output)


if __name__ == "__main__":
    main()

