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
	the number of bond types: 2
		their thresholds (inf sup): 1 2 3 4
	the number of angle types: 1
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
	N_tbonds = 0
	thresholds = []
	N_tangles = 0
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

			if N_tbonds == 0 and re.match("^bond_types:", line):
				N_tbonds = int(line.split()[1])
				continue

			if len(thresholds) == 0 and re.match("^\s*thresholds:", line):
				thresholds = [float(t) for t in line.split()[1:]]
				continue

			if N_tangles == 0 and re.match("^angle_types:", line):
				N_tangles = int(line.split()[1])
				continue

			if len(boundaries) == 0 and re.match("^space_boundaries:", line):
				boundaries = [float(b) for b in line.split()[1:]]
				continue

	return configuration_filename, output_filename, types, masses, charges, N_tbonds, thresholds, N_tangles, boundaries


def read_xyz(file_name: str) -> pd.DataFrame:
	"""
	Reads the configuration file and converts it to a Pandas Dataframe.
	"""

	return pd.read_csv(file_name, skiprows=2, delim_whitespace=True, names=["Type", "x", "y", "z"], dtype={"Type": "category", "x": "float64", "y": "float64", "z": "float64"})


def compute_bonds(data: pd.DataFrame, N_tbonds: int, threshes: list) -> np.ndarray:
	"""
	Computes a 3D-array of atoms indices which are bonded to each other, based on the inter-atomic distance and the thresholds given by the user.

	Parameters
	----------
	data: pandas.Dataframe
			The Dataframe generated by read_xyz.
	N_tbonds: int
			The number of types of bonds, given by the user in the input file.
	threshes: list
			The thresholds (inferior and superior limits) corresponding to the types of bonds.

	Returns
	-------
	numpy.ndarray:
			A 3D-array with the first axis for the bond types, the second for the different bonds and the third for the atoms concerned in the bond.
	"""

	positions = np.asarray(
		[data["x"], data["y"], data["z"]], dtype="float32").T
	bonds = [[] for _ in range(N_tbonds)]
	for i in range(len(positions)):
		for j in range(i+1, len(positions)):
			for k in range(0, N_tbonds):
				if threshes[2*k] <= np.sqrt(np.sum((positions[i] - positions[j])**2)) <= threshes[2*k+1]:
					bonds[k].append([i, j])
	return np.array(bonds, dtype="uint16", ndmin=3)


def angles_condition(data: pd.DataFrame, bonded_to: list) -> bool:
	return np.all(data["Type"][bonded_to] == 'H')


def compute_angles(data: pd.DataFrame, N_tangles: int, threshes: list) -> np.ndarray:
	"""
	Computes a 3D-array of atoms indices which are bonded to each other two by two and whose angles must be constrained.

	Parameters
	----------
	data: pandas.Dataframe
			The Dataframe generated by read_xyz.
	N_tangles: int
			The number of angles, given by the user in the input file.
	threshes: list
			The thresholds (inferior and superior limits) corresponding th the types of bonds.

	"""

	positions = np.asarray(
		[data["x"], data["y"], data["z"]], dtype="float64").T
	angles = [[] for _ in range(N_tangles)]
	for i in range(len(positions)):
		bonded_to = []
		for j in range(len(positions)):
			for k in range(0, 2*N_tangles, 2):
				if threshes[k] <= np.sqrt(np.sum((positions[i] - positions[j])**2)) <= threshes[k+1]: # No condition on the bond type...
					bonded_to.append(j)
					if len(bonded_to) == 2 and angles_condition(data, bonded_to):
						bonded_to.insert(1, i)
						angles[k].append(bonded_to)
						break

	return np.array(angles, dtype="uint16", ndmin=3)


'''
Old implementation: unreadable, too complicated, too much calls to numpy functions, probably inefficient
def compute_molecules(bonds: np.ndarray) -> np.ndarray:
	"""
	
	"""
	molecules = []
	indices = dict()
	for k, t in enumerate(bonds):
		for i, b in enumerate(t):
			if len(molecules) > 0 and np.any(np.isin(b, np.ravel(molecules))):
				for a in b:
					if a in np.ravel(molecules):
						indices[k, i] = list(np.argwhere(a == molecules).flatten())[0]
						break
			else:
				molecules.append(list(b))
	for k, v in indices.items():
		molecules[v].extend(list(bonds[k[0]][k[1]]))
		molecules[v] = list(np.unique(molecules[v]))
	return np.array(molecules, dtype="uint8")
'''


def compute_molecules(bonds: np.ndarray) -> np.ndarray:
	"""
	Computes a 2D-array of atoms indices which are forming molecules.

	Parameters
	----------
	bonds: numpy.ndarray
			The array generated by compute_bonds.

	Returns
	-------
	numpy.ndarray
			A 2D-array with the first axis corresponding to molecules and the second axis to the atoms indices.
	"""

	molecules = []
	other = []
	bonds = np.concatenate(bonds, axis=0).tolist()
	for bond in bonds:
		flag = False
		for atom in bond:
			for molecule in molecules:
				if atom in molecule:
					other[:] = bond[:]
					other.remove(atom)
					molecule.extend(other)
					flag = True
					break
			if flag:
				break
		else:
			molecules.append(bond)
	return np.array(molecules, dtype="uint16")


# Using argparse for its convenience:
# for the usage and help messages, and to raise an error if the argument is missing
my_parser = argparse.ArgumentParser(description="Converts an XYZ file to a LAMMPS data file.")
my_parser.add_argument("input_file", help="The input file, containing informations relative to the simulation.")


if __name__ == '__main__':
	arguments = my_parser.parse_args()

	# Reading the input file
	(conf_filename, output_filename, atom_types, masses,
  charges, N_tbonds, thresholds, N_tangles, boundaries) = read_input(arguments.input_file)
	
	# Reading the configuration file
	data = read_xyz(conf_filename)

	# Adding the missing columns to the Dataframe, that is:
	# the charges and the atoms labels according to their types
	for i in range(len(atom_types)):
		data.loc[data["Type"] == atom_types[i], [
			"Charge", "Label"]] = [charges[i], i+1]
	data["Label"] = data["Label"].convert_dtypes()
	
	# Computing the bonds and angles
	bonds = compute_bonds(data, N_tbonds, thresholds)
	# angles = compute_angles(data, N_tangles, thresholds)
	
	# Computing the molecules
	molecules = compute_molecules(bonds)
	mol_ID = np.zeros((len(data)), dtype="uint16")
	for i in range(len(data)):
		for j, molecule in enumerate(molecules):
			if i in molecule:
				mol_ID[i] = j
	data["Molecule"] = pd.Series(mol_ID)

	# Writing everything to the LAMMPS data file
	with open(output_filename, "w") as f:
		f.write("LAMMPS Description\n\n")
		f.write(f"  {len(data)} atoms\n")
		# f.write(f"  {len(np.concatenate(bonds, axis=0))} bonds\n")
		# f.write(f"  {len(np.concatenate(angles, axis=0))} angles\n\n")
		f.write(f"  {len(atom_types)} atom types\n")
		# f.write(f"  {N_tbonds} bond types\n")
		# f.write(f"  {N_tangles} angle types\n\n")
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
				f"{d[0]+1} {d[7]+1} {d[6]} {d[5]} {d[2]} {d[3]} {d[4]}\n")
		'''
		f.write("\n")
		f.write("Bonds\n\n")
		for i in range(N_tbonds):
			for j in range(len(bonds[i])):
				f.write(f"{j+1} {i+1} ")
				for k in bonds[i][j]:
					f.write(f"{k+1} ")
				f.write("\n")
		f.write("\n")
		f.write("Angles\n\n")
		for i in range(N_tangles):
			for j in range(len(angles[i])):
				f.write(f"{j+1} {i+1} ")
				for k in angles[i][j]:
					f.write(f"{k+1} ")
				f.write("\n")
		'''
