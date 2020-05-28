#!/usr/bin/env python
# Python wrapper program for parallelization of PRISM (part 3 of 3)
# Goal 1: produce post-processed yields for every particle from an SNSPH run
# Goal 2: avoid having to submit 950k sbatch jobs--that way lies madness
# Goal 3: organize all the output from PRISM after the processing is done
# This script exists to clean up all of the PRISM outputs after the fact
# Probably best to run this script as a big interactive job

# Last modified 28 May 2020 by Greg Vance

import sys
import os

import json
import glob
import h5py

import numpy as np

#import sdf_to_prism_translation as s2p

# 
TEAMMATE_DATA_PREFIX = "teammate"

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
	
	# Read the sdf data into Python once again
	# We need this so we can put the trajectory data into the hdf5 files
	#translator = s2p.SdfToPrismTranslator(sdf_dir)
	#particle_ids = translator.get_particle_ids()
	
	# mostly need to consolidate output into a single big abun hdf5 file
	# this should try to match the expectations of the data mining pipeline
	
	# Get all the isotopes from all the files
	isotopes_array = get_complete_isotopes(output_dir)
	
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
			% (info(["hdf5 dir"]))
	if info["team size"] <= 0:
		raise ValueError("team size must be a positive integer")
	if info["hdf5 count"] <= 0:
		raise ValueError("number of hdf5 files must be a positive integer")
	
	return info

def get_complete_isotopes(output_dir):
	"""
	"""
	
	complete_isotopes = set()
	
	data_file_names = sorted(glob.glob(os.path.join(output_dir, "%s_*.dat" \
		% (TEAMMATE_DATA_PREFIX))))
	for data_file_name in data_file_names:
		with open(data_file_name, "r") as data_file:
			
			for line in data_file:
				bits = line.split()
				if len(bits) == 3:
					Z, A, X = bits
					isotope = (int(Z), int(A) - int(Z))
					complete_isotopes.add(isotope)
	
	nz_nn = {"names": ("nz", "nn"), "formats": ("int32", "int32")}
	return np.array(sorted(list(complete_isotopes)), dtype=nz_nn)

def clean_up_files(output_dir):
	"""
	"""
	
	return

main()

