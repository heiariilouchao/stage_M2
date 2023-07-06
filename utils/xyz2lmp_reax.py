"""XYZ2LMP
Converts an XYZ file to a LAMMPS data file.

Arguments
---------
An input file which contains the informations relative to the LAMMPS system, that is:
	the configuration file name: blabla.xyz
	the output file name: blabla.txt
	the atom types and labels (according to the XYZ file): A B C
		their masses: 1 2 3
		their charges: 4 5 6
	the space boundaries: -5 -5 -5 5 5 5
"""

import argparse
import pandas as pd
import numpy as np


def read_input(file_name: str) -> tuple:
	"""
	Reads the input file, finding the informations using regular expressions.

	Parameters
	----------
	file_name: str
			The name of the input file.

	Returns
	-------
	tuple
			The data extracted from the input file, that is the configuration file name, the output file name, the types of atoms, their masses and charges, the number of types of bonds and angles, the thresholds for the bonds and the space boundaries.
	"""

	import re

	configuration_filename = ""
	output_filename = ""
	types = []
	masses = []
	charges = []
	boundaries = []

	with open(file_name, "r") as file:
		for line in file:
			if configuration_filename == "" and re.match("^configuration_file:", line):
				configuration_filename = line.split()[1]
				continue

			if output_filename == "" and re.match("^output_file:", line):
				output_filename = line.split()[1]
				continue

			if len(types) == 0 and re.match("^atom_types:", line):
				types = line.split()[1:]
				continue

			if len(masses) == 0 and re.match("^\s*masses:", line):
				masses = [float(m) for m in line.split()[1:]]
				continue

			if len(charges) == 0 and re.match("^\s*charges:", line):
				charges = [float(c) for c in line.split()[1:]]
				continue

			if len(boundaries) == 0 and re.match("^space_boundaries:", line):
				boundaries = [float(b) for b in line.split()[1:]]
				continue

	return configuration_filename, output_filename, types, masses, charges, boundaries


def read_xyz(file_name: str) -> pd.DataFrame:
	"""
	Reads the configuration file and converts it to a Pandas Dataframe.
	"""

	return pd.read_csv(file_name, skiprows=2, delim_whitespace=True, names=["Type", "x", "y", "z"],
		    dtype={"Type": "category", "x": "float64", "y": "float64", "z": "float64"})


# Using argparse for its convenience:
# for the usage and help messages, and to raise an error if the argument is missing
my_parser = argparse.ArgumentParser(description="Converts an XYZ file to a LAMMPS data file.")
my_parser.add_argument("input_file", help="The input file, containing informations relative to the simulation.")


if __name__ == '__main__':
	arguments = my_parser.parse_args()

	# Reading the input file
	(conf_filename, output_filename, atom_types, masses,
  charges, boundaries) = read_input(arguments.input_file)
	
	# Reading the configuration file
	data = read_xyz(conf_filename)

	# Adding the missing columns to the Dataframe, that is:
	# the charges and the atoms labels according to their types
	for i in range(len(atom_types)):
		data.loc[data["Type"] == atom_types[i], ["Charge", "Label"]] = [charges[i], i+1]
	data["Label"] = data["Label"].convert_dtypes()


	# Writing everything to the LAMMPS data file
	with open(output_filename, "w") as f:
		f.write("LAMMPS Description\n\n")
		f.write(f"  {len(data)} atoms\n")
		f.write(f"  {len(atom_types)} atom types\n")
		f.write(f"  {boundaries[0]} {boundaries[3]} xlo xhi\n")
		f.write(f"  {boundaries[1]} {boundaries[4]} ylo yhi\n")
		f.write(f"  {boundaries[2]} {boundaries[5]} zlo zhi\n\n\n")
		f.write("Masses\n\n")
		for i in range(len(atom_types)):
			f.write(f"{i+1} {masses[i]} # {atom_types[i]}\n")
		f.write("\n")
		f.write("Atoms\n\n")
		for d in data.itertuples():
			f.write(
				f"{d[0]+1} {d[6]} {d[5]} {d[2]} {d[3]} {d[4]}\n")
