# Python code for translating SNSPH SDF data into formats expected by PRISM

# Last modified 11 May 2020 by Greg Vance

import os.path
import glob
import numpy as np
import sdfpy

# The mass and radius of the sun in CGS units according to SNSPH
MSUN = 1.9889E33
RSUN = 6.955E10

class Abundances:
	"""Class to handle reading of the abun.dat file format and provide access
	methods for the data."""
	
	def __init__(self, file_name):
		
		# Sanity-check the file name
		self.file_name = file_name
		if not os.path.exists(self.file_name):
			raise ValueError("abun file does not exist: %s" \
				% (self.file_name))
		
		# The first five lines are the header info, so start with those
		with open(self.file_name, "r") as abun:
			header = [abun.readline().strip() for i in xrange(5)]
		self.network_size = int(header[0])
		self.n_zones = int(header[1])
		self.protons = [int(Z) for Z in header[2].split()]
		self.neutrons = [int(N) for N in header[3].split()]
		self.isotopes = [iso for iso in header[4].split() if iso != "r"]
		
		# Enlist numpy to process all of the float data after the header
		# The numbers in these files always seem to be single-precision
		# Don't waste RAM using dtype="float64" for this task
		file_array = np.loadtxt(self.file_name, "float32", skiprows=5)
		
		# Extract the first column, which is the list of zone radii
		# The rest of the columns should contain abundance data
		self.radii = np.copy(file_array[:,0])
		self.abun_data = np.copy(file_array[:,1:])
		del file_array
		
		# Before proclaiming success, go ahead and sanity-check *everything*
		assert len(self.protons) == self.network_size
		assert len(self.neutrons) == self.network_size
		assert len(self.isotopes) == self.network_size
		assert self.radii.shape == (self.n_zones,)
		assert self.abun_data.shape == (self.n_zones, self.network_size)
		
		# Produce a list with tuples containing both protons and neutrons
		self.protons_neutrons = list(zip(self.protons, self.neutrons))
		
		# Produce a list of atomic mass numbers A
		self.atomic_masses = [sum(pair) for pair in self.protons_neutrons]
	
	# Simple access methods
	
	def get_network_size(self):
		return self.network_size
	
	def get_n_zones(self):
		return self.n_zones
	
	def get_protons(self):
		return list(self.protons)
	
	def get_neutrons(self):
		return list(self.neutrons)
	
	def get_isotopes(self):
		return list(self.isotopes)
	
	def get_radii(self):
		return np.copy(self.radii)
	
	def get_abun_data(self):
		return np.copy(self.abun_data)
	
	def get_protons_neutrons(self):
		return list(self.protons_neutrons)
	
	def get_atomic_masses(self):
		return list(self.atomic_masses)
	
	# More complicated access methods
	
	def find_isotope_index(self, iso_or_Z, N=None):
		"""Find the index of a given isotope in the network list. With one
		argument, it expects a string that is the name of an isotope. With two
		arguments, it expects proton number, then neutron number as ints."""
		
		if N is None:
			return self.isotopes.index(iso_or_Z)
		else:
			return self.protons_neutrons.index((iso_or_Z, N))
	
	def find_abun_column(self, iso_or_Z, N=None):
		"""Return a copy of the abundance column with the abundance of a given
		isotope across all zones."""
		
		index = self.find_isotope_index(iso_or_Z, N)
		return np.copy(self.abun_data[:,index])
	
	def construct_zone_abun(self, zone, tuple_format):
		"""Return a list of all the abundance data for a given zone. The
		tuple_format parameter specifies the desired format of each list
		entry. For example, the tuple ('iso', 'Z', 'N', 'A', 'X', 'Y')."""
		
		zone_abun = list()
		for i in xrange(self.network_size):
			
			next_entry = list()
			for string in tuple_format:
				if string == "iso":
					next_entry.append(self.isotopes[i])
				elif string == "Z":
					next_entry.append(self.protons[i])
				elif string == "N":
					next_entry.append(self.neutrons[i])
				elif string == "A":
					next_entry.append(self.atomic_masses[i])
				elif string == "X":
					next_entry.append(self.abun_data[zone,i])
				elif string == "Y":
					next_entry.append(self.abun_data[zone,i] \
						/ self.atomic_masses[i])
				else:
					raise ValueError("tuple_format has unidentified" \
						" quantity \"%s\"" % (string))
				
			zone_abun.append(tuple(next_entry))
		
		return zone_abun

class ZoneToPartId:
	"""Class to handle reading of the zonetopartid.dat file format and provide
	access methods for the data."""
	
	def __init__(self, file_name):
		
		# Sanity-check the file name
		self.file_name = file_name
		if not os.path.exists(self.file_name):
			raise ValueError("ztpi file does not exist: %s" \
				% (self.file_name))
		
		# This is a simple file that just contains two coluns of ints
		self.particle_ids, self.zones = list(), list()
		with open(self.file_name, "r") as ztpi:
			for line in ztpi:
				stripped = line.strip()
				if stripped == "":
					continue
				else:
					pair = [int(string) for string in stripped.split()]
					assert len(pair) == 2
					particle_id, zone = tuple(pair)
					self.particle_ids.append(particle_id)
					self.zones.append(zone)
		
		# Store the length of the file
		self.length = len(self.particle_ids)
		
		# As a sanity check, ensure there are no repeated particle ids
		assert len(set(self.particle_ids)) == self.length
		
		# Also check that all the zone indices are non-negative
		assert min(self.zones) >= 0
		
		# And the particle ids should all be unsigned ints
		assert min(self.particle_ids) >= 0
	
	# Simple access methods
	
	def get_particle_ids(self):
		return list(self.particle_ids)
	
	def get_zones(self):
		return list(self.zones)
	
	def get_length(self):
		return self.length
	
	# More complicated access methods
	
	def find_particle_id_index(self, particle_id):
		"""Find the index of a given particle id in the particle ids list."""
		
		return self.particle_ids.index(particle_id)
	
	def find_zone(self, particle_id):
		"""Return the zone associated with a particular particle id."""
		
		index = self.find_particle_id_index(particle_id)
		return self.zones[index]
	
	def find_particle_ids(self, zone):
		"""Return a list of all the particle ids associated with a particular
		zone."""
		
		zone_particle_ids = list()
		for i in xrange(self.length):
			if self.zones[i] == zone:
				zone_particle_ids.append(self.particle_ids[i])
		
		return zone_particle_ids

class SdfToPrismTranslator:
	"""Class to handle translating a directory of SNSPH SDF data to inputs for
	PRISM in order to facilitate the postprocessing of yields."""
	
	def __init__(self, sdf_dir):
		
		# Sanity-check the directory path
		self.sdf_dir = sdf_dir
		if not os.path.exists(self.sdf_dir):
			raise ValueError("sdf_dir path does not exist: %s" \
				% (self.sdf_dir))
		
		# Make a sorted list of all SDF file names in the directory
		# These are the files with extensions made up entirely of digits
		self.sdf_name_list = list()
		for file_name in glob.glob(os.path.join(self.sdf_dir, "*")):
			dot_ext = os.path.splitext(file_name)[1]
			if dot_ext == "":
				continue
			if dot_ext[0] == '.' and dot_ext[1:].isdigit():
				self.sdf_name_list.append(file_name)
		self.sdf_name_list.sort()
		self.n_sdfs = len(self.sdf_name_list)
		
		# Assume that the other necessary files have standard names
		self.abun_file_name = os.path.join(self.sdf_dir, "abun.dat")
		self.ztpi_file_name = os.path.join(self.sdf_dir, "zonetopartid.dat")
		
		# Open all the SDF files simultaneously and keep them in a list
		# I might need to re-think this strategy if it requires too much RAM
		self.sdf_list = [sdfpy.SDFRead(sdf_name) for sdf_name in \
			self.sdf_name_list]
		self.sdf_list.sort(key=lambda x: x.parameters["tpos"])
		
		# Construct a set object containing all particle ids from the SDFs
		self.particle_id_set = set()
		for sdf in self.sdf_list:
			self.particle_id_set.update(sdf["ident"])
		self.n_particles = len(self.particle_id_set)
		
		# Read the other necesary files using the classes for them
		self.abun = Abundances(self.abun_file_name)
		self.ztpi = ZoneToPartId(self.ztpi_file_name)
		
		# Run every remaining sanity check I can think of
		assert max(self.ztpi.get_zones()) <= self.abun.get_n_zones() - 1
		assert self.particle_id_set.issubset(set(self.ztpi.particle_ids))
	
	# Simple access methods
	
	def get_n_sdfs(self):
		return self.n_sdfs
	
	def get_n_particles(self):
		return self.n_particles
	
	def get_particle_ids(self):
		particle_ids = list(self.particle_id_set)
		particle_ids.sort()
		return particle_ids
	
	# Method for imposing a temperature cutoff
	
	def get_particle_ids_with_temp_cutoff(self, temp_cutoff):
		"""Return a list of all particle ids whose trajectories reach or
		exceed the given temperature cutoff value (in kelvin). For high enough
		temperature cutoff values (e.g., 1e8 K), this method can be used to
		filter out particles that never get hot enough for nucleosynthesis to
		occur. Only the particle ids returned by this method actually need to
		be postprocessed by PRISM."""
		
		# For each particle id present, check its temperature in every SDF
		# If we find a single SDF with high enough temperature, save the id
		# Once we find one SDF where this is true, we can check the next id
		hot_particle_ids_list = list()
		for particle_id in self.particle_id_set:
			for sdf in self.sdf_list:
				index_array = np.argwhere(sdf["ident"] == particle_id)
				assert index_array.size == 1
				particle_index = index_array[0]
				if sdf["temp"][particle_index] >= temp_cutoff:
					hot_particle_ids_list.append(particle_id)
					break
		
		hot_particle_ids_list.sort()
		return hot_particle_ids_list
	
	# Translation methods for producing PRISM inputs
	
	def write_initial_composition_file(self, particle_id, output_file_name):
		"""Given a particle id and an ouput file name, write a PRISM "initial
		composition" file detailing the initial composition of that particle
		by matching the particle id with the correct zone abundances from the
		1d progenitor model."""
		
		# First, match the particle id with the correct 1d zone index
		zone_index = self.ztpi.find_zone(particle_id)
		
		# Pull the data for the initial composition of that zone
		# PRISM wants this information in the format (Z, A, X)
		prism_format = ("Z", "A", "X")
		composition_tuples = self.abun.construct_zone_abun(zone_index,
			prism_format)
		
		# Convert the tuples to simple space-delimited strings and join them
		composition_strings = list()
		for Z, A, X in composition_tuples:
			string = "%d %d %.6E" % (Z, A, X)
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
		length that will be ignored by PRISM."""
		
		global MSUN, RSUN
		
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
		for sdf in self.sdf_list:
			tpos = sdf.parameters["tpos"]
			index_array = np.argwhere(sdf["ident"] == particle_id)
			assert index_array.size == 1
			particle_index = index_array[0]
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
		for t, T9, RHO in trajectory_tuples:
			string = "%.8E %.8E %.8E" % (t, T9, RHO)
			trajectory_strings.append(string)
		output = "\n".join(trajectory_strings)
		
		# Write the final string of data to the output file
		with open(output_file_name, "w") as output_file:
			output_file.write(output)

