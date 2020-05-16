# Python code for translating SNSPH SDF data into input files for PRISM
# A necessary step if we want to use PRISM for postprocessing SNSPH particles

# Last modified 16 May 2020 by Greg Vance

import os.path
import glob
import numpy as np
import sdfpy

# Assumed standard file names for files supporting the SDF data
STD_ABUN_FILE_NAME = "abun.dat"
STD_ZTPI_FILE_NAME = "zonetopartid.dat"

# The mass and radius of the sun in CGS units according to SNSPH
MSUN = 1.9889E+33 # g
RSUN = 6.955E+10 # cm

# Helpful array algorithms

def test_if_sorted(array):
	"""Return True if the given Numpy array is sorted in ascending order."""
	
	i = 0
	while i < array.size - 1:
		if array[i] > array[i + 1]:
			return False
		i += 1
	return True

def test_if_unique(array):
	"""Return True if the given Numpy array contains only unique elements."""
	
	return np.unique(array).size == array.size

def find_one_index(array, value, assume_sorted=False):
	"""Return index such that array[index] == value. Elements of the array are
	assumed to be unique and IndexError is raised if value cannot be found.
	The finding process can be sped up using a binary search if the array is
	already known to be sorted.
	"""
	
	# When the array is sorted, we can use a binary search to find the value
	if assume_sorted:
		
		low_index = 0
		high_index = array.size - 1
		while low_index <= high_index:
			middle_index = low_index + (high_index - low_index) / 2
			if value < array[middle_index]:
				high_index = middle_index - 1
			elif value == array[middle_index]:
				return middle_index
			else: # value > array[middle_index]:
				low_index = middle_index + 1
	
	# Without the sorted assumption, default to a basic linear search
	else:
		
		index = 0
		while index < array.size:
			if array[index] == value:
				return index
			index += 1
	
	# Raise an error in case of search failure
	raise IndexError("array does not contain value: %s" % (repr(value)))

def find_all_indices(array, value, assume_sorted=False):
	"""Return an array of indices indicating every place where array == value.
	Elements of array are not assumed to be unique and an empty array will be
	returned if value is not found. The finding process can be sped up using a
	modified binary search if the array is already known to be sorted.
	"""
	
	# When the array is sorted, we can use a modified binary search algorithm
	# There may be multiple copies of value, but they'll be in one big clump
	if assume_sorted:
		
		n = array.size
		low_index1 = 0
		high_index1 = n - 1
		low_index2 = 0
		high_index2 = n - 1
		index1 = None
		index2 = None
		
		# Search for index1, the first index in the clump
		while low_index1 <= high_index1:
			middle_index1 = low_index1 + (high_index1 - low_index1) / 2
			if value < array[middle_index1]:
				high_index1 = middle_index1 - 1
				high_index2 = middle_index1 - 1
			elif value == array[middle_index1]:
				low_index2 = middle_index1
				if middle_index1 == 0 or array[middle_index1 - 1] < value:
					index1 = middle_index1
					break
				else:
					high_index1 = middle_index1 - 1
			else: # value > array[middle_index1]
				low_index1 = middle_index1 + 1
				low_index2 = middle_index1 + 1
		
		# If we couldn't find a start for the clump, then there is no clump
		if index1 is None:
			return np.array([], 'int')
		
		# Search for index 2, the last index in the clump
		while low_index2 <= high_index2:
			middle_index2 = low_index2 + (high_index2 - low_index2) / 2
			if value < array[middle_index2]:
				high_index2 = middle_index2 - 1
			elif value == array[middle_index2]:
				if middle_index2 == n - 1 or array[middle_index2 + 1] > value:
					index2 = middle_index2
					break
				else:
					low_index2 = middle_index2 + 1
			else: # value > array[middle_index2]
				low_index2 = middle_index2 + 1
		
		# This should never happen if we already found a start of the clump
		# It probably indicates that the array wasn't actually sorted
		# Either that, or something is wrong with our search algorithm
		if index2 is None:
			raise IndexError("unexplained failure in find_all_indices")
		
		# Use the clump's limits to produce the output array
		return np.arange(index1, index2 + 1, dtype='int')
	
	# If unsorted, then just use Numpy to brute-force check every element
	else:
		
		# The nonzero function returns a single-element tuple here,
		# so index into it to produce the 1d array that I'm expecting
		return np.nonzero(array == value)[0]

# File translation classes

class Abundances:
	"""Class to handle reading of the abun.dat file format and provide access
	methods for the data.
	"""
	
	def __init__(self, file_name):
		
		# Sanity-check the file name
		self.file_name = file_name
		if not os.path.exists(self.file_name):
			raise ValueError("abun file does not exist: %s" \
				% (repr(self.file_name)))
		
		# The first five lines are the header info, so start by reading those
		# Convert the header data to Numpy arrays where sequences are expected
		with open(self.file_name, "r") as abun_file:
			header = [abun_file.readline().strip() for h in xrange(5)]
		self.network_size = int(header[0])
		self.n_zones = int(header[1])
		self.protons = np.array(header[2].split(), "int32")
		self.neutrons = np.array(header[3].split(), "int32")
		self.isotopes = np.array([iso for iso in header[4].split() \
			if iso != "r"])
		
		# Use Numpy to process all of the float data following the header
		# The numbers in these files always seem to be single-precision,
		# so don't waste RAM by using the default float64 for this task
		file_array = np.loadtxt(self.file_name, "float32", skiprows=5)
		
		# Extract the first column, which is the list of zone radii
		# The rest of the columns should all contain abundance data
		self.radii = np.copy(file_array[:,0])
		self.abun_data = np.copy(file_array[:,1:])
		del file_array
		
		# Before proclaiming success, go ahead and sanity-check *everything*
		# The first two and last two tests check the file for self-consistency
		# Tests 3 and 4 ensure that the proton and neutron numbers make sense
		# Tests 5 and 6 ensure that the list of isotopes has no repeats
		assert self.protons.size == self.network_size
		assert self.neutrons.size == self.network_size
		assert np.all(self.protons >= 0) and np.all(self.protons < 1000)
		assert np.all(self.neutrons >= 0) and np.all(self.neutrons < 1000)
		assert test_if_unique(self.protons * 1000 + self.neutrons)
		assert test_if_unique(self.isotopes)
		assert self.radii.shape == (self.n_zones,)
		assert self.abun_data.shape == (self.n_zones, self.network_size)
		
		# Produce an array containing the atomic mass numbers A
		self.atomic_masses = self.protons + self.neutrons
	
	# Simple access methods
	
	def get_network_size(self):
		return self.network_size
	
	def get_n_zones(self):
		return self.n_zones
	
	def get_protons(self):
		return np.copy(self.protons)
	
	def get_neutrons(self):
		return np.copy(self.neutrons)
	
	def get_isotopes(self):
		return np.copy(self.isotopes)
	
	def get_radii(self):
		return np.copy(self.radii)
	
	def get_abun_data(self):
		return np.copy(self.abun_data)
	
	def get_atomic_masses(self):
		return np.copy(self.atomic_masses)
	
	# More complicated access methods
	
	def find_isotope_index(self, iso_or_Z, N=None):
		"""Find the index of a given isotope in the network list. With one
		argument, it expects a string that is the name of an isotope. With two
		arguments, it expects proton number, then neutron number as ints.
		"""
		
		# I don't trust myself to sort these, so just use a linear search
		if N is None:
			return find_one_index(self.isotopes, iso_or_Z)
		else:
			proton_match = find_all_indices(self.protons, iso_or_Z)
			neutron_match = find_one_index(self.neutrons[proton_match], N)
			return proton_match[neutron_match]
		
	def extract_abun_column(self, iso_or_Z, N=None):
		"""Return a copy of the abundance column with the abundance of a given
		isotope across all zones.
		"""
		
		index = self.find_isotope_index(iso_or_Z, N)
		return np.copy(self.abun_data[:,index])
	
	def construct_zone_abun(self, zone, tuple_format):
		"""Return a list of all the abundance data for a given zone. The
		tuple_format parameter specifies the desired format of each list
		entry. For example, the tuple ('iso', 'Z', 'N', 'A', 'X', 'Y').
		"""
		
		# Build up a list of all the columns requested by the format tuple
		# These are all just *references* to the appropriate arrays of data
		# This ought to run fast since we aren't actually copying anything yet
		# One exception: the values for Y aren't stored anywhere in this class
		# If Y is requested, we'll have to calculate those values on the fly
		columns = list()
		for item in tuple_format:
			if item == "iso":
				columns.append(self.isotopes)
			elif item == "Z":
				columns.append(self.protons)
			elif item == "N":
				columns.append(self.neutrons)
			elif item == "A":
				columns.append(self.atomic_masses)
			elif item == "X":
				columns.append(self.abun_data[zone,:])
			elif item == "Y":
				# Note [section 3.3.3 of PRISM manual]: X(Z,A) := Y(Z,A) * A
				columns.append(self.abun_data[zone,:] \
					/ np.array(self.atomic_masses, "float32"))
			else:
				raise ValueError("tuple_format contains unknown item: %s" \
					% (repr(item)))
		
		# Once the references have been organized, zip them all together
		return zip(*columns)

class ZoneToPartId:
	"""Class to handle reading of the zonetopartid.dat file format and provide
	access methods for the data.
	"""
	
	def __init__(self, file_name):
		
		# Sanity-check the file name
		self.file_name = file_name
		if not os.path.exists(self.file_name):
			raise ValueError("ztpi file does not exist: %s" \
				% (repr(self.file_name)))
		
		# This is a simple file containing two columns of integers
		self.particle_ids, self.zones = np.loadtxt(self.file_name, "int32",
			unpack=True)
		
		# Store the total length of the file
		self.length = self.particle_ids.size
		
		# Run a few sanity checks on the file contents
		# The first checks that the particle ids are unique
		# The other two ensure that there are no negative values
		assert test_if_unique(self.particle_ids)
		assert np.all(self.zones >= 0)
		assert np.all(self.particle_ids >= 0)
		
		# For fast searching, make sure the arrays are sorted by particle id
		sort_order = np.argsort(self.particle_ids)
		self.particle_ids = np.copy(self.particle_ids[sort_order])
		self.zones = np.copy(self.zones[sort_order])
	
	# Simple access methods
	
	def get_particle_ids(self):
		return np.copy(self.particle_ids)
	
	def get_zones(self):
		return np.copy(self.zones)
	
	def get_length(self):
		return self.length
	
	# More complicated access methods
	
	def find_particle_id_index(self, particle_id):
		"""Find the index of a given particle id in the particle ids list."""
		
		return find_one_index(self.particle_ids, particle_id,
			assume_sorted=True)
	
	def find_zone(self, particle_id):
		"""Return the zone associated with a particular particle id."""
		
		index = self.find_particle_id_index(particle_id)
		return self.zones[index]
	
	def find_particle_ids(self, zone):
		"""Return an array containing all the particle ids associated with a
		given zone.
		"""
		
		indices = find_all_indices(self.zones, zone)
		return np.copy(self.particle_ids[indices])

class SdfToPrismTranslator:
	"""Primary translator class to handle translating a directory of SNSPH SDF
	data to create inputs for PRISM in order to facilitate the postprocessing
	of yields.
	"""
	
	def __init__(self, sdf_dir):
		
		# Sanity-check the directory path
		self.sdf_dir = sdf_dir
		if not os.path.exists(self.sdf_dir):
			raise ValueError("sdf_dir path does not exist: %s" \
				% (repr(self.sdf_dir)))
		
		# Make a sorted list of all SDF file names in the directory
		# These are the files with extensions made up entirely of digits
		self.sdf_file_names = list()
		for file_name in glob.glob(os.path.join(self.sdf_dir, "*.*")):
			dot_ext = os.path.splitext(file_name)[1]
			if dot_ext == "":
				continue
			if dot_ext[0] == '.' and dot_ext[1:].isdigit():
				self.sdf_file_names.append(file_name)
		self.sdf_file_names.sort()
		self.n_sdfs = len(self.sdf_file_names)
		
		# Assume that the other necessary files have standard names
		self.abun_file_name = os.path.join(self.sdf_dir, STD_ABUN_FILE_NAME)
		self.ztpi_file_name = os.path.join(self.sdf_dir, STD_ZTPI_FILE_NAME)
		
		# Open all the SDF files simultaneously and keep them in a list
		# I might need to re-think this strategy if it requires too much RAM
		self.sdf_files = [sdfpy.SDFRead(sdf_file_name) for sdf_file_name \
			in self.sdf_file_names]
		self.sdf_files.sort(key=lambda sdf: sdf.parameters["tpos"])
		
		# Construct an array containing all unique particle ids from the SDFs
		self.particle_ids = reduce(np.union1d, [sdf["ident"] for sdf \
			in self.sdf_files])
		self.n_particles = self.particle_ids.size
		
		# Sanity-check some of the contents of the SDFs
		# First, make sure there are no repeated particle ids
		# Then, ensure self-consistency with the number of particles
		for sdf in self.sdf_files:
			assert test_if_unique(sdf["ident"])
			assert sdf.parameters["npart"] == sdf["ident"].size
			assert sdf.parameters["npart"] <= self.n_particles
		
		# To hopefully decrease search times, check if each SDF is sorted
		self.sdf_is_sorted = np.array([test_if_sorted(sdf["ident"]) for sdf \
			in self.sdf_files], "bool")
		
		# Read the other necesary files using the classes for them
		self.abun = Abundances(self.abun_file_name)
		self.ztpi = ZoneToPartId(self.ztpi_file_name)
		
		# Run every remaining sanity check I can think of
		# The first test ensures that the files agree on the number of zones
		# The second test ensures that they agree on the set of particle ids
		assert np.all(self.ztpi.get_zones() <= self.abun.get_n_zones() - 1)
		assert np.all(np.in1d(self.particle_ids, self.ztpi.particle_ids,
			assume_unique=True))
	
	# Simple access methods
	
	def get_n_sdfs(self):
		return self.n_sdfs
	
	def get_n_particles(self):
		return self.n_particles
	
	def get_particle_ids(self):
		return np.copy(self.particle_ids)
	
	# More complicated access methods
	
	def find_particle_id_indices(self, particle_id):
		"""Return an array containing the index of a given particle id in all
		SDFs ordered by the SDF tpos values.
		"""
		
		indices = np.full(self.n_sdfs, -1, dtype="int")
		n_fails = 0
		
		for i, sdf in enumerate(self.sdf_files):
			
			# If the index from the previous SDF works, then use that
			if i > 0 and sdf["ident"][indices[i - 1]] == particle_id:
				indices[i] = indices[i - 1]
			
			# Otherwise, we'll have to search for the particle id
			else:
				try:
					indices[i] = find_one_index(sdf["ident"], particle_id,
						assume_sorted=self.sdf_is_sorted[i])
				except IndexError:
					n_fails += 1
		
		if n_fails > 0:
			if n_fails == self.n_sdfs:
				raise IndexError("unable to find particle id: %s" \
					% (repr(particle_id)))
			print "Warning: SDF data for particle id %d is incomplete" \
				% (particle_id)
		
		return indices
	
	# Method for imposing a temperature cutoff on the particles
	
	def impose_temperature_cutoff(self, temp_cutoff):
		"""Return a boolean array the same size as the array of particle ids
		indicating those particles whose trajectories reach or exceed the
		given temperature cutoff value (in kelvin). For high enough
		temperature cutoff values (e.g., 1e8 K), this method can be used to
		filter out particles that never get hot enough for nucleosynthesis to
		occur. Only the particle ids which are marked as True by this method
		actually need to be postprocessed by PRISM.
		"""
		
		return reduce(np.logical_or, [(sdf["temp"] >= temp_cutoff) for sdf \
			in self.sdf_files])
	
	# Translation methods for producing PRISM inputs
	
	def write_initial_composition_file(self, particle_id, output_file_name):
		"""Given a particle id and an ouput file name, write a PRISM "initial
		composition" file detailing the initial composition of that particle
		by matching the particle id with the correct zone abundances from the
		1d progenitor model.
		"""
		
		# First, match the particle id with the correct 1d zone index
		zone_index = self.ztpi.find_zone(particle_id)
		
		# Pull the data for the initial composition of that zone
		# PRISM wants this information in the format (Z, A, X)
		prism_format = ("Z", "A", "X")
		composition_tuples = self.abun.construct_zone_abun(zone_index,
			prism_format)
		
		# Convert the tuples to simple space-delimited strings and join them
		composition_strings = list()
		for Z_A_X in composition_tuples:
			string = "%d %d %.6E" % Z_A_X
			composition_strings.append(string)
		output = "\n".join(composition_strings)
		
		# Write the final string of data to the output file
		with open(output_file_name, "w") as output_file:
			output_file.write(output)
	
	def write_trajectory_file(self, particle_id, output_file_name,
		HEADER=None):
		"""Given a particle ID and an output file name, write a PRISM
		"trajectory" file detailing the thermodynamic evolution of that
		particle by extracting the temperature and density of that particle
		from each SDF along with that SDF's simulation time. PRISM trajectory
		files are allowed to begin with a single HEADER line of arbitrary
		length that will be ignored by PRISM.
		"""
		
		# Default HEADER if none is provided
		if HEADER is None:
			HEADER = "PRISM Trajectory (t, T9, RHO) for Particle ID %d" \
				% (particle_id)
		
		# Ensure that the HEADER is exactly one single line
		if HEADER != "" and HEADER[-1] == '\n':
			HEADER = HEADER[:-1]
		if HEADER.count('\n') != 0:
			raise ValueError("PRISM trajectory file HEADER must be one line")
		
		# Pull the tpos value, particle temp, and particle rho from each SDF
		sdf_trajectory_tuples = list()
		particle_indices = self.find_particle_id_indices(particle_id)
		for sdf, particle_index in zip(self.sdf_files, particle_indices):
			if particle_index == -1:
				continue
			tpos = sdf.parameters["tpos"]
			temp = sdf["temp"][particle_index]
			rho = sdf["rho"][particle_index]
			sdf_trajectory_tuples.append((tpos, temp, rho))
		
		# The times, temps, and rhos from the SDFs are all in SNSPH code units
		# These numbers need to be converted to the units that PRISM expects
		# PRISM wants to see seconds, gigakelvin, and grams/centimeter**3
		trajectory_tuples = list()
		for tpos, temp, rho in sdf_trajectory_tuples:
			t = tpos * 100.0  # s per SNSPH_TIME
			T9 = temp * 1E-9  # GK per SNSPH_TEMP
			RHO = rho * (1E-6 * MSUN * RSUN ** -3)  # g/cm3 per SNSPH_DENSITY
			trajectory_tuples.append((t, T9, RHO))
		
		# Convert the tuples to simple space-delimited strings and join them
		trajectory_strings = [HEADER]
		for t_T9_RHO in trajectory_tuples:
			string = "%.8E %.8E %.8E" % t_T9_RHO
			trajectory_strings.append(string)
		output = "\n".join(trajectory_strings)
		
		# Write the final string of data to the output file
		with open(output_file_name, "w") as output_file:
			output_file.write(output)

# If this file is being executed instead of imported, run a few tests
if __name__ == "__main__":
	
	TEST_DIR = "/home/gsvance/vconv_3d_asym_data/vconv_snsph_data/"
	t = SdfToPrismTranslator(TEST_DIR)
	
	t.abun.get_network_size()
	t.abun.get_n_zones()
	t.abun.get_protons()
	t.abun.get_neutrons()
	t.abun.get_isotopes()
	t.abun.get_radii()
	t.abun.get_abun_data()
	t.abun.get_atomic_masses()
	
	t.abun.find_isotope_index('fe61')
	t.abun.find_isotope_index(26, 61 - 26)
	t.abun.extract_abun_column('fe61')
	t.abun.extract_abun_column(26, 61 - 26)
	t.abun.construct_zone_abun(1234, ('X', 'Y', 'Z', 'A', 'iso', 'N'))
	
	t.ztpi.get_particle_ids()
	t.ztpi.get_zones()
	t.ztpi.get_length()
	
	t.ztpi.find_particle_id_index(2468)
	t.ztpi.find_zone(2468)
	t.ztpi.find_particle_ids(1234)
	
	t.get_n_sdfs()
	t.get_n_particles()
	t.get_particle_ids()
	
	t.find_particle_id_indices(13579)
	t.impose_temperature_cutoff(1E+8)
	t.write_initial_composition_file(5500, "test/test_icomp.txt")
	t.write_trajectory_file(5500, "test/test_traj.txt")
	
	print "All tests ran without raising exceptions"

