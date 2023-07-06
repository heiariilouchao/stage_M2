"""

"""

import argparse
import numpy as np
import io
import pandas as pd


def read_parameters(filename: str) -> tuple:
	"""
	Reads the parameters from the configuration file.

	Parameters
	----------
	filename: str
			The name of the configuration file
	
	Returns
	-------
	N_atoms: int
			The number of atoms.

	bounds: numpy.ndarray
			The bounds of the simulation box.
	"""

	text = ""

	# Opening the file and reading only the first 256 characters (~ the first 9 lines)
	with open(filename, "r") as f:
		text = f.read(256)
	
	data = text.splitlines()
	
	if len(data) < 8: # We need to read the first 8 lines
		raise RuntimeError("Could not read the parameters from the configuration file.")

	# Extracting the important fields
	return int(data[3]), np.array([data[i].split() for i in range(5, 8)], dtype="float64")


def read_configuration(file: io.TextIOWrapper) -> pd.DataFrame:
	"""
	Reads a single configuration from a file.

	Parameters
	----------
	file: io.TextIOWrapper
			The configuration file.
	
	Returns
	-------
	data: pandas.DataFrame
			The DataFrame containing the current positions and unwrapped positions of the atoms.
	"""

	indices = []
	types = []
	positions = []
	upositions = []

	for l in range(9 + N_atoms):
		if l < 9: # The positions are not yet reached
			if file.readline() == '': # If the read line is empty the end of the file was reached
				return None
		else:
			line = file.readline()
			if line == '':
				return None

			# Splitting the different fields of the line and extracting the important ones
			fields = line.split()
			indices.append(fields[0])
			types.append(fields[1])
			positions.append([x for x in fields[2:5]])
			upositions.append([x for x in fields[5:]])
	
	# Converting the data into numpy arrays
	indices = np.array(indices, dtype="uint32")
	types = np.array(types)
	positions = np.array(positions, dtype="float64")
	upositions = np.array(upositions, dtype="float64")
	
	# Constructing and returning the data as a pandas DataFrame
	return pd.DataFrame({"indices": indices,
						 "types": pd.Series(types, dtype="category"),
						 "x": positions[:, 0],
						 "y": positions[:, 1],
						 "z": positions[:, 2],
						 "ux": upositions[:, 0],
						 "uy": upositions[:, 1],
						 "uz": upositions[:, 2]})


def compute_distance(r1: np.ndarray, r2: np.ndarray, bounds: np.ndarray, hbounds: np.ndarray) -> float:
	"""
	Computes the distance between two positions taking into account the periodic boundary conditions.

	Parameters
	----------
	r1: numpy.ndarray
			The position of the first atom.
	
	r2: numpy.ndarray
			The position of the second atom.
	
	bounds: numpy.ndarray
			The lengths of the simulation box.

	hbounds: numpy.ndarray
			The half-lengths of the simulation box.

	Returns
	-------
	r12: float
			The distance between r1 and r2 (r2 - r1).
	"""

	# The element-wise subtraction
	distance = r2 - r1

	'''
	# Vectorized way: somehow, it is slower than the loop...
	distance = np.where(distance < -hbounds, distance + 2 * hbounds, distance)
	distance = np.where(distance > hbounds, distance - 2 * hbounds, distance)
	'''
	# Applying the periodic boundary conditions over the differences
	for i in range(len(distance)):
		if distance[i] < -hbounds[i]: # If the distance is lower than the negative half-box
			distance[i] += bounds[i]
		elif distance[i] > hbounds[i]: # Else if the distance is higher than the half-box
			distance[i] -= bounds[i]

	# Returns the square root of the dot product of the vector with itself
	return np.sqrt(np.dot(distance, distance))


def increment_histogram(data: pd.DataFrame, types: list, bounds: np.ndarray, hbounds: np.ndarray, histogram: np.ndarray, delta: float) -> np.ndarray:
	"""
	Increments the histogram.

	Parameters
	----------
	data: pandas.DataFrame
			The DataFrame returned by read_configuration, containing the positions of the atoms.
	
	types: list
			The list of types labels of the pair.

	bounds: numpy.ndarray
			The lengths of the simulation box.

	hbounds: numpy.ndarray
			The half-lengths of the simulation box.

	histogram: numpy.ndarray
			The histogram to increment.

	delta: float
			The space step (cutoff distance / number of bins).
	
	Returns
	-------
	histogram: numpy.ndarray
			The incremented histogram.
	"""

	r1 = data.loc[data["types"] == types[0], ["x", "y", "z"]].to_numpy()

	# Copying the histogram to increment it
	new_histogram = histogram.copy()

	if types[0] != types[1]:
		r2 = data.loc[data["types"] == types[1], ["x", "y", "z"]].to_numpy()
		# Looping over all pairs of atoms
		for r_i in r1:
			for r_j in r2:
				dist = compute_distance(r_i, r_j, bounds, hbounds)
				index = int(dist / delta)
				if index < len(histogram):
					new_histogram[index] += 1
	else:
		N_1 = len(r1)
		# Looping over the same array for faster computation
		for i in range(N_1):
			for j in range(i + 1, N_1): # Looping only over the 'other' atoms
				dist = compute_distance(r1[i], r1[j], bounds, hbounds)
				index = int(dist / delta)
				if index < len(histogram):
					new_histogram[index] += 2   # Incrementing the entry by 2 because the pair counts twice

	# Returning the incremented copy of the histogram
	return new_histogram



def compute_RDF(histogram: np.ndarray, r: np.ndarray, delta: float, N_steps: int, N_pair: list, volume: float) -> np.ndarray:
	"""
	Computes the Radial Distribution Function.

	Parameters
	----------
	histogram: numpy.ndarray
			The histogram fully incremented.
	
	r: numpy.ndarray
			The distance range.

	delta: float
			The space step (cutoff distance / number of bins).

	N_steps: int
			The number of configurations the histogram got incremented over.
	
	N_pair: list
			The number of atoms for the pair.
	
	volume: float
			The volume of the simulation box.
	
	Returns
	-------
	RDF: numpy.ndarray
			The Radial Distribution Function.
	"""

	# Computing the number of pairs
	N = N_pair[0]
	if len(N_pair) == 2:
		N *= N_pair[1]
	else:
		N *= (N - 1)
	
	# Computing the mean number of atoms
	n = histogram / N_steps
	n_ideal = 4. / 3. * np.pi * N / volume * ((r + delta)**3 - r**3)
	
	# Returning the RDF
	return n / n_ideal


def compute_MSD(initial: np.ndarray, current: pd.DataFrame, N_atoms: int) -> float:
	"""
	Computes the Mean Squared Displacement.

	Parameters
	----------
	initial: numpy.ndarray
			The initial unwrapped positions.
	
	current: pandas.DataFrame
			The DataFrame returned by read_configuration, containing the current unwrapped positions.
	
	N_atoms: int
			The number of atoms.
	
	Returns
	-------
	MSD: float
			The current Mean Squared Displacement.
	"""

	# The difference between the current and initial positions
	# To ensure that the positions match to the initial array we sort the values by the indices
	displacement = current.sort_values("indices")[["ux", "uy", "uz"]].to_numpy() - initial

	# The tensor dot product is used for efficient computation
	return np.tensordot(displacement, displacement) / N_atoms


# The argument parser
my_parser = argparse.ArgumentParser(description="Processes the output data generated by LAMMPS for configurations.")
# The script's arguments and options
my_parser.add_argument("configuration_file", type=str, help="The dump file from LAMMPS containing the configurations.")
my_parser.add_argument("-nbins", default=50, type=int, help="The number of bins for the PDF. Default is 50.")
my_parser.add_argument("-maxsteps", default=50, type=int, help="The maximum number of configurations to extract from the file. Default is 50.")
my_parser.add_argument("-display", default="true", type=str, help="To display the plots or not. Default is true.")


if __name__ == '__main__':
	# Parsing the arguments and options
	args = my_parser.parse_args()
	conf_file = args.configuration_file
	N_bins = args.nbins
	N_maxsteps = args.maxsteps
	display = args.display.lower() in ["true", "t", "1"]

	# Reading the configuration file
	print("Reading the parameters...")
	N_atoms, limits = read_parameters(conf_file)

	# Extracting/Computing the simulation parameters
	bounds = np.subtract(limits[:, 1], limits[:, 0])
	hbounds = bounds / 2.
	cutoff = np.min(hbounds)
	delta = cutoff / N_bins
	V = np.prod(bounds)
	
	# Initializing the arrays
	OH_list = ["O", "H"]
	OH_histogram = np.zeros((N_bins), dtype="uint32")
	OO_list = ["O", "O"]
	OO_histogram = np.zeros((N_bins), dtype="uint32")
	HH_list = ["H", "H"]
	HH_histogram = np.zeros((N_bins), dtype="uint32")
	MSD = np.zeros((N_maxsteps), dtype="float64")

	# Reading the file
	print("Reading the file...")
	with open(conf_file, 'r') as file:
		data = read_configuration(file)

		# Compute some paramaters once before reading the whole file
		N_steps = 0
		N_OH = [data[data["types"] == t].count()["types"] for t in np.unique(OH_list)]
		N_OO = [data[data["types"] == t].count()["types"] for t in np.unique(OO_list)]
		N_HH = [data[data["types"] == t].count()["types"] for t in np.unique(HH_list)]

		# Extracting the initial positions for the MSD
		initial_pos = data.sort_values("indices")[["ux", "uy", "uz"]].to_numpy()

		# Reading the file configuration by configuration
		print("Reading the configurations...")
		while data is not None and N_steps < N_maxsteps:
			print(f"\rStep: {N_steps + 1} / {N_maxsteps}", end='')   # For the impatient
			N_steps += 1

			# Incrementing the histograms
			OH_histogram = increment_histogram(data, OH_list, bounds, hbounds, OH_histogram, delta)
			OO_histogram = increment_histogram(data, OO_list, bounds, hbounds, OO_histogram, delta)
			HH_histogram = increment_histogram(data, HH_list, bounds, hbounds, HH_histogram, delta)

			# Computing the MSD
			MSD[N_steps - 1] = compute_MSD(initial_pos, data, N_atoms)
			RMSD = np.sqrt(MSD[:N_steps])

			# Reading the next configuration
			data = read_configuration(file)

	# Computing the PDFs
	print("\nAveraging the RDFs...")
	r = np.linspace(delta, cutoff, N_bins, endpoint=False)
	OH_PDF = compute_RDF(OH_histogram, r, delta, N_steps, N_OH, V)
	OO_PDF = compute_RDF(OO_histogram, r, delta, N_steps, N_OO, V)
	HH_PDF = compute_RDF(HH_histogram, r, delta, N_steps, N_HH, V)
	# Shifting the distance range for the PDF (placing it at the center of the bins)
	r += .5 * delta
	t = np.array(range(N_steps))

	# Plotting
	if display:
		import matplotlib.pyplot as plt

		plt.figure()
		plt.title("Radial Distribution Function")
		plt.xlabel("r [Å]")
		plt.ylabel("g")
		plt.plot(r, OH_PDF, label="OH")
		plt.plot(r, OO_PDF, label="OO")
		plt.plot(r, HH_PDF, label="HH")
		plt.legend()

		plt.figure()
		plt.title("Mean Squared Displacement")
		plt.xlabel("t")
		plt.ylabel("MSD [Å²]")
		plt.plot(t, MSD[:N_steps])
		
		plt.figure()
		plt.title("Root Mean Squared Displacement")
		plt.xlabel("t")
		plt.ylabel("RMSD [Å]")
		plt.plot(t, RMSD)

		plt.show()
	
	# Saving the results
	np.savetxt("../data/rdf.dat",
	    		np.transpose([r, OH_PDF, OO_PDF, HH_PDF]),
				fmt="%5f",
				header="{:5s} {:5s} {:5s} {:5s}".format("r", "g_OH", "g_OO", "g_HH"))
	np.savetxt("../data/msd.dat",
	    		np.transpose([t, MSD[:N_steps], RMSD]),
				fmt="%5f",
				header="{:5s} {:5s} {:5s}".format("t", "MSD", "RMSD"))
