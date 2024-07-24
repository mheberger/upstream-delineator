r"""
Delineation of watershed subbasins, using data from
MERIT-Basins and MERIT-Hydro.
Created by Matthew Heberger, May 2024.

See README for more detailed instructions.

Quick Start:

First, set parameters in the file config.py.

Run this script from the command line with two required arguments:
$ python subbasins.py outlets.csv testrun

or with a full file path as follows on Windows or Linux:
$ python subbasins.py C:\Users\matt\Desktop\outlets.csv test
$ python subbasins.py /home/files/outlets.csv test

or in Python as follows:
>> from subbasins import delineate
>> delineate('outlets.csv', 'testrun')

"""

# Standard Python libraries. See requirements.txt for recommended versions.
import argparse
import sys
sys.path.insert(0, ".")

# My stuff
from upstream_delineator.delineator_utils.delineate import delineate


def _run_from_terminal():
    """
    Routine which is run when you call subbasins.py from the command line.
    should be run like:
    python subbasins.py outlets.csv run_name
    """
    # Create the parser for command line inputs.
    description = "Delineate subbasins using data in and input CSV file. Writes " \
        "a set of output files beginning with the output prefix string."

    parser = argparse.ArgumentParser(description=description)

    # Add the arguments
    parser.add_argument('input_csv', help="Input CSV filename, for example 'gages.csv'")
    parser.add_argument('output_prefix', help="Output prefix, a string. The output files will start with this string")

    # Parse the arguments
    args = parser.parse_args()

    # Call the main function, passing the command line arguments
    delineate(args.input_csv, args.output_prefix)


def main():
    # Run directly, for convenience or during development and debugging
    input_csv = 'outlets.csv'
    out_prefix = 'iceland'
    delineate(input_csv, out_prefix)


if __name__ == "__main__":
    if len(sys.argv) > 1:
        # Run with command-line arguments
        _run_from_terminal()
    else:
        main()
