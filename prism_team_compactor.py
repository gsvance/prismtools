#!/usr/bin/env python
# Python wrapper program for parallelization of PRISM (part 3 of 3)
# Goal 1: produce post-processed yields for every particle from an SNSPH run
# Goal 2: avoid having to submit 950k sbatch jobs--that way lies madness
# Goal 3: organize all the output from PRISM after the processing is done
# This script exists to clean up all of the PRISM outputs after the fact
# Probably best to run this script as an interactive job

# Last modified 1 Aug 2020 by Greg Vance

import sys
import os
import time as tm

import json
import glob
import h5py

import numpy as np

import sdf_to_prism_translation as s2p

# Prefixes for input and output files
TEAMMATE_DATA_PREFIX = "teammate"
HDF5_DATA_PREFIX = "compdata"

def main():
	
	info = get_compactor_info()
	
	# Clean up the inputs and store them in more convenient variables
	# Note that the "output dir" here is the output dir from the PRISM jobs
	# Most of this script's output will actually go to files in the hdf5 dir
	sdf_dir = os.path.abspath(info["sdf dir"])
	output_dir = os.path.abspath(info["output dir"])
	hdf5_dir = os.path.abspath(info["hdf5 dir"])
	team_size = int(info["team size"])
	hdf5_count = int(info["hdf5 count"])
	del info
	
	# We will need some of the SDF data for the HDF5 files, so load it now
	print "Loading SDF data into a translator..."
	translator = s2p.SdfToPrismTranslator(sdf_dir)
	print "SDF data loaded"
	print
	
	# Assemble a list of the PRISM data files from the teammate scripts
	# Next, work out the set of particle ids and the full nuclear network
	print "Organizing input data from teammate scripts..."
	data_file_names = get_data_file_names(output_dir, team_size)
	particle_ids, isotopes = get_particle_ids_and_isotopes(translator,
		data_file_names)
	print "Data files found, particle id list ready, isotope list ready"
	print
	
	# Set up object-oriented approaches to our input and output files
	print "Starting up objects for input and output file management..."
	hdf5_files = HDF5Group(hdf5_dir, hdf5_count, particle_ids, isotopes)
	data_files = DataIter(data_file_names, translator, isotopes)
	print "Manager objects are ready"
	print
	
	# Pull data blocks from the input files particle-by-particle
	# Write them to the output HDF5 files until we run out of data
	print "Transferring data to HDF5 files now..."
	n_part = 0
	t0 = tm.time()
	data_block = data_files.get_next_data_block()
	while data_block is not None:
		hdf5_files.write_cycle(*data_block)
		n_part += 1
		if n_part % 5000 == 0:
			print "  %d/%d particles transferred (%.0f%%)" % (n_part, \
				len(particle_ids), n_part * 100.0 / len(particle_ids))
			tn = tm.time()
			speed = float(n_part) / (tn - t0)
			t_left = (len(particle_ids) - n_part) / speed
			print "  estimated %.0f minutes remaining" % (t_left / 60.0)
		data_block = data_files.get_next_data_block()
	print "All data transferred"
	print
	
	print "Starting clean up phase..."
	
	# Close all of the output HDF5 files
	hdf5_files.close_all()
	
	# Delete all the leftover files that aren't needed any longer
	clean_up_files(output_dir)
	
	print "Clean up done"
	print

def get_compactor_info():
	"""Parse this script's single command line argument, read in the compactor
	info from the specified JSON file, run a few quick checks, and then return
	the dict of compactor info.
	"""
	
	# Check for the one command line argument that we want to see
	# We are only expecting to get the name of a small input JSON file
	if len(sys.argv) != 2:
		raise TypeError("one command line argument is required: " \
			+ "name of JSON compactor info file")
	json_file_name = sys.argv[1]
	
	# Read the little JSON file and load the dict of compactor info
	with open(json_file_name, "r") as json_file:
		info = json.load(json_file)
	
	# Pretty-print the input values we read from the file
	print "Read compactor info from input file:\n%s" % (json_file_name)
	print json.dumps(info, indent=2, separators=(",", ": "))
	print
	
	# Sanity-check a few of the input values before returning them
	if not os.path.exists(info["sdf dir"]):
		raise ValueError("SDF directory does not exist: %s" \
			% (info["sdf dir"]))
	if not os.path.exists(info["output dir"]):
		raise ValueError("output directory does not exist: %s" \
			% (info["output dir"]))
	if not os.path.exists(info["hdf5 dir"]):
		raise ValueError("hdf5 directory does not exist: %s" \
			% (info["hdf5 dir"]))
	if info["team size"] <= 0:
		raise ValueError("team size must be a positive integer")
	if info["hdf5 count"] <= 0:
		raise ValueError("number of hdf5 files must be a positive integer")
	
	return info

def get_data_file_names(output_dir, team_size):
	"""Produce a sorted list of all the teammate script output data files
	containing the final compositions from PRISM that this script needs to
	read and organize. We are looking for a set of team_size files, so raise
	an error if we don't have that number of files.
	"""
	
	# Use glob to match all of the file names, then sort the list
	match_string = os.path.join(output_dir, "%s_*.dat" \
		% (TEAMMATE_DATA_PREFIX))
	names_list = sorted(glob.glob(match_string))
	
	# Check the number of files and return the list if no problems
	if len(names_list) != team_size:
		raise ValueError("expected to find %d teammate files, but found %d" \
			% (team_size, len(names_list)))
	else:
		return names_list

def get_particle_ids_and_isotopes(translator, data_file_names):
	"""Rifle through both the SDF data and the PRISM data files to get
	comprehensive lists of the set of particle ids and the netwrok isotopes
	that need to be handled in this script.
	"""
	
	# Start by pulling all the available data from the original SDF files
	particle_ids = set([int(pid) for pid in translator.get_particle_ids()])
	nz_array = translator.abun.get_protons()
	nn_array = translator.abun.get_neutrons()
	isotopes = set()
	for i in xrange(translator.abun.get_network_size()):
		isotopes.add((int(nz_array[i]), int(nn_array[i])))
	
	# Read the teammate script data files one-by-one and line-by-line
	# Make sure we don't miss any more particle ids or network isotopes
	for data_file_name in data_file_names:
		with open(data_file_name, "r") as data_file:
			for line in data_file:
				
				line_parts = line.strip().split()
				
				# Lines with 2 parts have an id and whether PRISM was used
				if len(line_parts) == 2:
					particle_id, prism_or_abun = line_parts
					particle_ids.add(int(particle_id))
				
				# Lines with 3 parts have isotope Z, A, and mass fraction X
				elif len(line_parts) == 3:
					Z, A, X = line_parts
					isotopes.add((int(Z), int(A) - int(Z)))
	
	# Make a sorted list out of the set of unique particle ids
	particle_ids = sorted(particle_ids)
	
	# Make a sorted list out of the set of unique isotope tuples
	# Move the free neutrons, protons, and alphas to the very end
	isotopes = sorted(isotopes)
	isotopes.remove((0, 1))  # free neutrons
	isotopes.remove((1, 0))  # protons (hydrogen-1)
	isotopes.remove((2, 2))  # alphas (helium-4)
	isotopes.extend([(0, 1), (1, 0), (2, 2)])
	
	return particle_ids, isotopes

class HDF5Group:
	"""Object-oriented approach to managing the output HDF5 files."""
	
	def __init__(self, hdf5_dir, hdf5_count, particle_ids, isotopes):
		"""Set up all the HDF5 files and prepare them for data writing."""
		
		self.hdf5_dir = hdf5_dir
		self.hdf5_count = hdf5_count
		self.particle_ids = list(particle_ids)
		self.isotopes = list(isotopes)
		
		# Calculate the number of particle ids to store in each HDF5 file
		# The last file will have slightly fewer unless this divides evenly
		self.ids_per_file = len(self.particle_ids) / self.hdf5_count
		remainder = len(self.particle_ids) % self.hdf5_count
		if remainder != 0:
			self.ids_per_file += 1
		
		# Open all of the new HDF5 files in a big list
		self.files = list()
		for i in xrange(self.hdf5_count):
			file_name = os.path.join(self.hdf5_dir, "%s_%03d.h5" \
				% (HDF5_DATA_PREFIX, i + 1))
			self.files.append(h5py.File(file_name, "w"))
		
		# Add the nz and nn datasets to each of the new HDF5 files
		# Also set each file's "version" attributes (might be important)
		sh = (len(self.isotopes),)
		dt = np.dtype([("data", "int32")])
		nz_array = np.empty(sh, dtype=dt)
		nn_array = np.empty(sh, dtype=dt)
		nz_array["data"][:] = [isotope[0] for isotope in self.isotopes]
		nn_array["data"][:] = [isotope[1] for isotope in self.isotopes]
		for hdf5 in self.files:
			hdf5.create_dataset("nz", shape=sh, dtype=dt, data=nz_array)
			hdf5.create_dataset("nn", shape=sh, dtype=dt, data=nn_array)
			hdf5.attrs["HDF5_version"] = np.array(["1.8.2"], dtype="S6")
			hdf5.attrs["SE_version"] = np.array(["1.2"], dtype="S4")
			hdf5.flush()
	
	def find_file_index_by_particle_id(self, particle_id):
		"""Given a particle id, figure out which of the HDF5 files that
		particle should be stored in.
		"""
		
		particle_index = self.particle_ids.index(particle_id)
		return particle_index / self.ids_per_file
	
	def write_cycle(self, particle_id, fmass, time, rho, temp, mass):
		"""Given all the data for a single particle, write one "cycle" to an
		appropriate output HDF5 file using that blob of data."""
		
		# Figure out which file the new cycle needs to be put into
		file_index = self.find_file_index_by_particle_id(particle_id)
		hdf5 = self.files[file_index]
		
		# Create the new HDF5 group with the requisite "cycle#####" name
		cycle_name = "cycle%010d" % (particle_id)
		cycle = hdf5.create_group(cycle_name)
		
		# Create the cycle's fmass dataset and fill it with numbers
		s1 = (len(self.isotopes),)
		d1 = np.dtype([("data", "float64")])
		fmass_array = np.empty(s1, dtype=d1)
		fmass_array["data"][:] = fmass
		cycle.create_dataset("fmass", shape=s1, dtype=d1, data=fmass_array)
		
		# Create the cycle's SE_DATASET and fill it with trajectory data
		s2 = (len(time),)
		d2 = np.dtype([("time", "float64"), ("rho", "float64"),
			("temp", "float64")])
		se_array = np.empty(s2, dtype=d2)
		se_array["time"][:] = time
		se_array["rho"][:] = rho
		se_array["temp"][:] = temp
		cycle.create_dataset("SE_DATASET", shape=s2, dtype=d2, data=se_array)
		
		# Set the attributes... I don't know how many of these are necessary
		# At the very least, the mass value will be needed by burn_query
		cycle.attrs["tfinal"] = np.array([99.9], dtype="float64")
		cycle.attrs["looperror"] = np.array([999], dtype="int32")
		cycle.attrs["mass"] = np.array([mass], dtype="float64")
		
		# Flush any extra data from the buffer and write it to the file
		# Is this actually necessary?
		#hdf5.flush()
	
	def close_all(self):
		"""Close the full set of HDF5 files."""
		
		while len(self.files) > 0:
			hdf5 = self.files.pop()
			hdf5.close()

class DataIter:
	"""Object-oriented approach to managing the input data files containing
	the composition data from PRISM. Also manages the SDf translator object as
	needed to obtain masses and trajectory data."""
	
	def __init__(self, data_file_names, translator, isotopes):
		"""Set up the data manager object and open its first file."""
		
		self.data_file_names = list(data_file_names)
		self.translator = translator  # Big object, do not make a copy
		self.isotopes = list(isotopes)
		
		# Set up the isotope index records for quick finding of isotopes
		self.isotope_index_records = dict()
		for i in xrange(len(self.isotopes)):
			isotope = self.isotopes[i]
			self.isotope_index_records[isotope] = i
		
		# Start reading the set of data files
		self.file_index = 0
		file_name = self.data_file_names[self.file_index]
		self.current_file = open(file_name, "r")
		self.next_particle = self.get_line()
	
	def get_line(self):
		"""Internal method for more-or-less extracting one line at a time from
		the set of files. If we reach EOF, open the next file seamlessly. When
		there are no more files, return None."""
		
		# Read lines from the current file until we find one that isn't blank
		line = self.current_file.readline()
		while line == "\n":
			line = self.current_file.readline()
		
		# Empty string means EOF, so close this file and move on
		if line == "":
			self.current_file.close()
			self.file_index += 1
			
			# If there are more files, then open the next one
			if self.file_index < len(self.data_file_names):
				file_name = self.data_file_names[self.file_index]
				self.current_file = open(file_name, "r")
				return self.get_line()
			
			# If we're out of files, then signal it with None
			else:
				self.current_file = None
				return None
		
		# When we do get a file line, return it as a list of split parts
		return line.strip().split()
	
	def find_isotope_index(self, isotope):
		"""Find and return the list index of a given isotope tuple."""
		
		return self.isotope_index_records[isotope]
	
	def get_next_data_block(self):
		"""Return the next full data block to be written to the output HDF5
		files. Reads and combines data from several places to accomplish this.
		"""
		
		# Check if we've run out of data files
		if self.next_particle is None:
			return None
		
		# Pull out the particle info from the next line of the file
		particle_id, prism_or_abun = self.next_particle
		particle_id = int(particle_id)
		
		# Start an array of fmass numbers and fill it with composition data
		fmass_array = np.zeros(len(self.isotopes), dtype="float64")
		line = self.get_line()
		while line is not None and len(line) == 3:
			Z, A, X = line
			isotope = (int(Z), int(A) - int(Z))
			isotope_index = self.find_isotope_index(isotope)
			fmass_array[isotope_index] = float(X)
			line = self.get_line()
		
		# Save the first line that doesn't have three elements for later on
		# This is probably the two-element line we need for the next particle
		self.next_particle = line
		
		# Grab the particle's trajectory and mass from the SDFs
		particle_traj = self.translator.write_trajectory_file(particle_id)
		particle_mass = self.translator.find_particle_mass(particle_id)
		
		# Return the data block as a big tuple for the HDF5Group class
		# Convert temp array from T9 (GK) to temps in kelvin
		return (particle_id, fmass_array, particle_traj["t"],
			particle_traj["RHO"], 1E+9 * particle_traj["T9"], particle_mass)

def clean_up_files(output_dir):
	"""Clean up subroutine. For now, do nothing."""
	
	pass
	return

main()

