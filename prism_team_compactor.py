#!/usr/bin/env python
# Python wrapper program for parallelization of PRISM (part 3 of 3)
# Goal 1: produce post-processed yields for every particle from an SNSPH run
# Goal 2: avoid having to submit 950k sbatch jobs--that way lies madness
# Goal 3: organize all the output from PRISM after the processing is done
# This script exists to clean up all of the PRISM outputs after the fact
# Probably best to run this script as an interactive job

# Last modified 29 Jul 2020 by Greg Vance

import sys
import os

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
	
	# 
	hdf5_files = HDF5Group(hdf5_dir, hdf5_count, particle_ids, isotopes)
	
	
	
	
	
	
	
	
	# 
	hdf5_files.close_all()
	
	
	
	
	
	
	
	
	print "copying fmass data now"
	for data_file_name in data_file_names:
		particle_id, source, cycle = None, None, None
		Z, A, X = None, None, None
		with open(data_file_name, "r") as data_file:
			for line in data_file:
				line_bits = line.strip().split()
				if len(line_bits) == 2:
					particle_id, source = line_bits
					cycle = "cycle%010d" % (int(particle_id))
					i_hdf5 = simple_index(particle_ids, int(particle_id)) \
						/ ids_per_file
					#print "i_hdf5", i_hdf5
					#hdf5_files[i_hdf5][cycle].attrs.create("source", source)
				elif len(line_bits) == 3:
					Z, A, X = line_bits
					ZA = np.array([(int(Z), int(A) - int(Z))],
						dtype=isotopes.dtype)
					index = simple_index(isotopes, ZA)
					#print "index", index
					hdf5_files[i_hdf5][cycle]["fmass"]["data"][index] \
						= float(X)
	
	
	
	
	
	# Delete all the leftover teammate script files that aren't needed
	clean_up_files(output_dir)

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
		isotopes.add((int(nz[i]), int(nn[i])))
	
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
	"""
	"""
	
	def __init__(self, hdf5_dir, hdf5_count, particle_ids, isotopes):
		
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
			file_name = os.path.join(self.hdf5_dir, "%s_03d.h5" \
				% (HDF5_DATA_PREFIX, i + 1))
			self.files.append(h5py.File(file_name, "w"))
		
		# Add the nz and nn datasets to each of the new HDF5 files
		sh = (len(self.isotopes),)
		dt = np.dtype([("data", "int32")])
		nz_array = np.empty(sh, dtype=dt)
		nn_array = np.empty(sh, dtype=dt)
		nz_array["data"][:] = [isotope[0] for isotope in self.isotopes]
		nn_array["data"][:] = [isotope[1] for isotope in self.isotopes]
		for hdf5 in self.files:
			hdf5.create_dataset("nz", shape=sh, dtype=dt, data=nz_array)
			hdf5.create_dataset("nn", shape=sh, dtype=dt, data=nn_array)
	
	def find_file_index_by_particle_id(self, particle_id):
		particle_index = self.particle_ids.index(particle_id)
		return particle_index / self.ids_per_file
	
	def write_cycle(self, particle_id, fmass, time, rho, temp, mass):
		
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
		# At the very least, the mass value is going to be needed
		cycle.attrs["tfinal"] = np.array(np.nan, dtype="float64")
		cycle.attrs["looperror"] = np.array(9999, dtype="int32")
		cycle.attrs["mass"] = np.array(mass, dtype="float64")
	
	def close_all(self):
		while len(self.files) > 0:
			hdf5 = self.files.pop()
			hdf5.close()





















def simple_index(array, target):
	"""
	"""
	
	for i in xrange(array.size):
		if array[i] == target:
			return i
	
	raise IndexError

def clean_up_files(output_dir):
	"""
	"""
	
	return

main()

