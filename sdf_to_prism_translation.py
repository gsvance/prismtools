# Python code for translating SNSPH SDF data into input files for PRISM
# A necessary step if we want to use PRISM for postprocessing SNSPH particles

# Last modified 15 May 2020 by Greg Vance

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
		
		# The first five lines are the header info, so start with those
		# Convert all of this data to Numpy arrays where possible
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
		# The middle two tests ensure that the array of isotopes is unique
		assert self.protons.size == self.network_size
		assert self.neutrons.size == self.network_size
		assert np.unique(np.column_stack([self.protons, self.neutrons]), \
			axis=0).shape == (self.network_size, 2)
		assert np.unique(self.isotopes).size == self.network_size
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
		
		if N is None:
			index_array = np.nonzero(self.isotopes == iso_or_Z)
		else:
			protons_match = (self.protons == iso_or_Z)
			neutrons_match = (self.neutrons == N)
			index_array = np.nonzero(protons_match & neutrons_match)
		
		# The constructor method checks that all the isotopes are unique
		# This should only be triggered if the isotope does not exist
		if index_array.size != 1:
			raise IndexError("failed to find isotope index: (%s, %s)" \
				% (repr(iso_or_Z), repr(N)))
		
		return index_array[0]
	
	def find_abun_column(self, iso_or_Z, N=None):
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
		assert np.unique(self.particle_ids).size == self.length
		assert np.all(self.zones >= 0)
		assert np.all(self.particle_ids >= 0)
	
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
		
		index_array = np.nonzero(self.particle_ids == particle_id)
		
		# The constructor method checks that all the particle ids are unique
		# This should only be triggered if the particle id does not exist
		if index_array.size != 1:
			raise IndexError("failed to find particle id: %s" \
				% (repr(particle_id)))
		
		return index_array[0]
	
	def find_zone(self, particle_id):
		"""Return the zone associated with a particular particle id."""
		
		index = self.find_particle_id_index(particle_id)
		return self.zones[index]
	
	def find_particle_ids(self, zone):
		"""Return an array containing all the particle ids associated with a
		given zone.
		"""
		
		return np.copy(self.particle_ids[self.zones == zone])

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
		
		# Read the other necesary files using the classes for them
		self.abun = Abundances(self.abun_file_name)
		self.ztpi = ZoneToPartId(self.ztpi_file_name)
		
		# Run every remaining sanity check I can think of
		# The second test ensures that the files agree on the number of zones
		# The third test ensures that they agree on the set of particle ids
		assert np.all(self.ztpi.get_zones() <= self.abun.get_n_zones() - 1)
		assert np.all(np.in1d(self.particle_ids, self.ztpi.particle_ids,
			assume_unique=True))
	
	# Simple access methods
	
	def get_n_sdfs(self):
		return self.n_sdfs
	
	def get_n_particles(self):
		return self.n_particles
	
	def get_particle_ids(self):
		particle_ids = list(self.particle_id_set)
		particle_ids.sort()
		return particle_ids
	
	# More complicated access methods
	
	def find_particle_id_indices(self, particle_id):
		"""
		"""
		
		raise NotImplementedError
	
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
			in sdf_files])
	
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
		particle_indices = self.find_particle_indices(particle_id)
		for sdf, particle_index in zip(self.sdf_files, particle_indices):
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

