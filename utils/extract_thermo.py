"""
	Data processing script: scans a LAMMPS log file and extracts the thermo output.
"""

import re
import argparse


# Configuring the argument parser
parser = argparse.ArgumentParser(description="Processing a LAMMPS log file to extract the thermodynamics quantities.")
parser.add_argument("log_file", help="The name of the log file.")
parser.add_argument("res_file", help="The name of the results file.")


if __name__ == '__main__':
	# Parsing the arguments passed to the script
	args = parser.parse_args()

	# Compiling the regular expressions
	pre_header_re = re.compile("^Per MPI rank memory allocation")
	header_re = re.compile("^[a-zA-Z]+\s")
	tail_re = re.compile("^Loop time of")

	# Opening the files
	with open(args.log_file, "r") as log, open(args.res_file, "w") as res:
		recording = False

		# Reading the log file
		for line in log:
			if recording:
				if tail_re.match(line):
					recording = False		# Canceling the recording
				elif header_re.match(line):
					res.write("# " + line)	# Writing the header as a comment
				else:
					res.write(line)			# Writing the record
			elif pre_header_re.match(line):
				recording = True			# Starting recording
