#!/usr/bin/env python
# Python wrapper program for parallelization of PRISM (part 2 of 3)
# Goal 1: produce post-processed yields for every particle from an SNSPH run
# Goal 2: avoid having to submit 950k sbatch jobs--that way lies madness
# Goal 3: organize all the output from PRISM after the processing is done
# This script is tasked with running all of the PRISM jobs assigned to it
# Run this script as a non-interactive cluster job!

# Last modified 22 May 2020 by Greg Vance

import sys
import os
import time

import json
import subprocess

import sdf_to_prism_translation as s2p

# Standard prefixes for each type of file that needs creating
ICOMPOSITION_PREFIX = "icomp"
TRAJECTORY_PREFIX = "traj"
CONTROL_PREFIX = "ctrl"
FCOMPOSITION_PREFIX = "fcomp"
PRISM_STDOUT_PREFIX = "prism"
PRISM_STDERR_PREFIX = "prism"

# Path to my PRISM install directory on Agave
PRISM_DIR = "/home/gsvance/prism/prism-1.5.0/"
# Path to the default control.json file that came with PRISM
DEFAULT_CONTROL = os.path.join(PRISM_DIR, "input/control.json")

# Flags to set for debug testing
PARTICLES_LIMIT = 4  # stop after processing some number of particles
DELETE_FILES = False  # whether to clean up files as we go or not

def main():
	
	# Start timer and print welcome messages
	t1 = time.time()
	print "<< PRISM Parallelization Teammate Script >>"
	print "Beginning setup procedure..."
	print
	
	parameters = read_input_parameters()
	
	# Extract the input parameters into shorter variable names
	sdf_dir = os.path.abspath(parameters["sdf dir"])
	output_dir = os.path.abspath(parameters["output dir"])
	team_size = int(parameters["team size"])
	team_rank = int(parameters["team rank"])
	temp_cut = float(parameters["temp cut"])
	del parameters
	
	# Print our team teammate status as we understand it
	print "This job is assigned to run team rank %d of team size %d." \
		% (team_rank, team_size)
	print
	
	# Set up the data translator object... this can take a while to do
	# Only need to do this once, so it should be a small part of the runtime
	print "Initializing SDF-to-PRISM data translator object..."
	translator = s2p.SdfToPrismTranslator(sdf_dir)
	print "Translator object is constructed and ready to proceed."
	
	# Get the full array of particle ids from the translator
	particle_ids = translator.get_particle_ids()
	print "SDF directory has data for %d particles." % (particle_ids.size)
	
	# Impose a temperature cutoff to limit the number of necessary PRISM runs
	print "Imposing temperature cutoff at %.3E kelvin..." % (temp_cut)
	particle_is_hot = translator.impose_temperature_cutoff(temp_cut)
	print "Found %d particles that will require processing by PRISM." \
		% (np.sum(particle_is_hot))
	print
	
	# Conclude the setup and report the amount of time it consumed
	print "Setup procedure is complete."
	t2 = time.time()
	t_setup = t2 - t1
	print "The setup took %.3f minutes to finish." % (t_setup / 60.)
	print
	
	################################
	##### SETUP CODE ENDS HERE #####
	################################
	
	# Switch over to PRISM's local install directory before we try running it
	# There's a lot of deafault files that it looks for in the subdirectories
	# We *could* run it elsewhere, but it's a lot easier to run it from here
	os.chdir(PRISM_DIR)
	
	# Establish the name of this teammate's big data output file
	# We will consolidate all final composition data here as we go along
	# Open the file in write mode for a moment to clear its contents
	data_file_name = os.path.join(output_dir, "teammate_%04d.dat" \
		% (team_rank))
	open(data_file_name, "w").close()
	print "Consolidation data file name: %s" % (data_file_name)
	
	# Loop over all the particles that were assigned to this teammate
	# THE RULE HERE: We are responsible for every index in the particle id
	# array that is congruent to our team_rank (modulo team_size).
	# Example: team_rank = 6, team_size = 10
	# We would take the particles at indices 6, 16, 26, 36, 46, ...
	print "Starting main loop over assigned particle ids..."
	print
	n_loops = 0
	for particle_index in xrange(team_rank, particle_ids.size, team_size):
		
		# If a debug limit was set on the number of particles, then check that
		if PARTICLES_LIMIT is not None and n_loops >= PARTICLES_LIMIT:
			print "Main loop reached PARTICLES_LIMIT debug value."
			print "Breaking out of main loop now."
			print
			break
		
		particle_id = particle_ids[particle_index]
		print "Now considering particle id %08d at index %d in array." \
			% (particle_id, particle_index)
		
		# If this particle needs post-processing, then go ahead and run PRISM
		if particle_is_hot[particle_index]:
			print "This particle will require processing by PRISM."
			print "Preparing files needed by PRISM..."
			file_names = prepare_prism_files(particle_id, translator,
				output_dir)
			print "Files are ready. Starting up PRISM subprocess..."
			run_prism_once(particle_id, file_names)
			print "PRISM subprocess has finished running."
			fcomposition_file_name = file_names["final composition"]
		
		# If the particle does NOT need post-processing, then don't run PRISM
		else:
			print "This particle will NOT require any processing by PRISM."
			print "Copying initial composition directly to output file..."
			fcomposition_file_name = write_fake_output_file(particle_id,
				translator, output_dir)
		
		# Consolidate the latest output data into our big output data file
		print "Consolidating final composition data into big file..."
		consolidate(particle_id, fcomposition_file_name, data_file_name)
		print "Consolidation completed."
		
		n_loops += 1
		print
	
	# The loop is done, so just finish up
	print "Main loop is now complete."
	print "A total of %d particles were considered in the loop." % (n_loops)
	t3 = time.time()
	t_loop = t3 - t2
	t_total = t3 - t1
	print "The main loop took %.3f minutes to run." % (t_loop / 60.)
	print "The entire program with setup took %.3f minutes to run." \
		% (t_total / 60.)
	print
	
	print "Script execution is now complete. Bye!"

def read_input_parameters():
	"""Parse the script's command line argument, read the input parameters
	from the specified JSON file, run some quick checks, and return the dict
	of parameters.
	"""
	
	# Check the number of command line arguments that were received
	# We are only expecting the name of a small JSON file with inputs
	if len(sys.argv) != 2:
		raise TypeError("one command line argument is required: " \
			+ "name of JSON input file")
	json_file_name = sys.argv[1]
	
	# Read the little JSON file and load the dict of input parameters
	with open(json_file_name, "r") as json_file:
		parameters = json.load(json_file)
	
	# Pretty-print the input values we read from the file
	print "Read parameters from input file:\n%s" % (json_file_name)
	print json.dumps(parameters, indent=2, separators=(",", ": "))
	print
	
	# Sanity-check a few of the input values before returning to main
	# Existence of SDF directory will be checked by the translator object
	if not os.path.exists(parameters["output dir"]):
		raise ValueError("output directory does not exist: %s" \
			% (parameters["output dir"]))
	if parameters["team rank"] not in range(parameters["team size"]):
		raise ValueError("team rank/size of %d/%d is not allowable" \
			% (parameters["team rank"], parameters["team size"]))
	if parameters["temp cut"] < 0.0:
		raise ValueError("bad temperature cutoff at %g K" \
			% (parameters["temp cut"]))
	
	return parameters

def prepare_prism_files(particle_id, translator, output_dir):
	"""Given a particle id, prepare all of the necessary input files that
	PRISM needs in order to postprocess yields for that particle. Return a
	dict of paths to the newly created files. Produce the following files:
	  - ASCII initial composition file
	  - ASCII trajectory file
	  - JSON control file
	Details about their formats can be found in the PRISM manual PDF. In
	addition to the particle id, this function will need access to the data
	translator and the path to the output directory where all the files should
	be saved. New files will be written to the output directory and PRISM will
	be instructed to save its output to the same output directory.
	"""
	
	# Create this particle's initial composition ASCII file for PRISM
	icomposition_file_name = os.path.join(output_dir, "%s_%08d.dat" \
		% (ICOMPOSITION_PREFIX, particle_id))
	translator.write_initial_composition_file(particle_id,
		icomposition_file_name)
	
	# Create this particle's trajectory ASCII file for PRISM
	trajectory_file_name = os.path.join(output_dir, "%s_%08d.dat" \
		% (TRAJECTORY_PREFIX, particle_id))
	translator.write_trajectory_file(particle_id, trajectory_file_name)
	
	# Create the PRISM JSON control file for handling this particle
	control_file_name = os.path.join(output_dir, "%s_%08d.json" \
		% (CONTROL_PREFIX, particle_id))
	fcomposition_file_name = os.path.join(output_dir, "%s_%08d.dat" \
		% (FCOMPOSITION_PREFIX, particle_id))
	write_control_file(control_file_name, icomposition_file_name,
		trajectory_file_name, fcomposition_file_name)
	
	# Set two file names to collect PRISM's stdout and stderr
	prism_stdout_file_name = os.path.join(output_dir, "%s_%08d.out" \
		% (PRISM_STDOUT_PREFIX, particle_id))
	prism_stderr_file_name = os.path.join(output_dir, "%s_%08d.err" \
		% (PRISM_STDERR_PREFIX, particle_id))
	
	# Return all the file names so they don't need to be reconstructed later
	file_names = dict()
	file_names["initial composition"] = icomposition_file_name
	file_names["trajectory"] = trajectory_file_name
	file_names["control"] = control_file_name
	file_names["final composition"] = fcomposition_file_name
	file_names["prism stdout"] = prism_stdout_file_name
	file_names["prism stderr"] = prism_stderr_file_name
	return file_names

def write_control_file(control_file_name, icomposition_file_name,
	trajectory_file_name, fcomposition_file_name, stop_temp=None):
	"""Write the JSON control file for a single particle using the requested
	file names and PRISM stop temperature (in kelvin). This mostly involves
	reading the deafult control file, making a few modifications, and then
	writing the modified data to a new JSON file. If no stop temperature is
	provided the control file will tell PRISM to stop when it reaches the end
	of the trajectory file instead.
	"""
	
	# Read the JSON data from the default control file that came with PRISM
	with open(DEFAULT_CONTROL, "r") as default_control_file:
		control = json.load(default_control_file)
	
	# For now, use the standard nuclear data files that came with PRISM
	#control["nuclear"]["datasets"] = [DEFAULT]
	
	# Change the conditions object to direct PRISM to the correct input files
	control["conditions"]["initial_composition"]["path"] \
		= icomposition_file_name
	control["conditions"]["trajectory"]["path"] = trajectory_file_name
	
	# Alter the network object to set PRISM's start and stop triggers
	control["network"]["start"] = {"comment":
		"Left blank to start calculation at beginning of trajectory file"}
	if stop_temp is not None:
		stop_T9 = stop_temp * 1E-9  # Convert stop temperature to GK
		control["network"]["stop"] = {"T9": stop_T9, "comment":
			"Stops calculation when trajectory reaches %.3E GK" % (stop_T9)}
	else:
		control["network"]["stop"] = {"comment":
			"Left blank to stop calculation at end of trajectory file"}
	#control["network"]["extent"] = [DEFAULT]
	
	# Use the output object to strongly limit PRISM's set of output files
	for output_key in control["output"].keys():
		control["output"][output_key]["active"] = False
	control["output"]["x"]["active"] = True
	control["output"]["x"]["path"] = fcomposition_file_name
	
	# Write the modified JSON data out to create the new control file
	# Use pretty printing when creating this file to make life easy for me
	with open(control_file_name, "w") as new_control_file:
		json.dump(control, new_control_file, indent=2, separators=(",", ": "))

def run_prism_once(particle_id, file_names):
	"""Given a particle id and the dict of related file names, start a
	subprocess and run PRISM once. If all goes well, clean up any unneeded
	files to avoid overwhelming the file system. If the exit code from PRISM
	indicates a problem, print an alert and don't delete the files--they could
	be important for manual debugging later.
	"""
	
	# Open two file streams to catch the stdout and stderr from PRISM
	stdout_file = open(file_names["prism stdout"], "w")
	stderr_file = open(file_names["prism stderr"], "w")
	
	# Run a PRISM subprocess using files that were prepared earlier
	# Remember that Python has already chdir'ed to PRISM's install directory
	command = ["./prism", "-c", file_names["control"]]
	prism_exit_code = subprocess.call(command, stdout=stdout_file,
		stderr=stderr_file)
	
	# Close the file streams collecting stdout and stderr from PRISM
	stdout_file.close()
	stderr_file.close()
	
	# Draw attention to any unhappy exit codes from PRISM
	if prism_exit_code != 0:
		print "ALERT: PRISM RETURNED BAD EXIT CODE %d" % (prism_exit_code)
		print "ALERT: THE OFFENDING PARTICLE ID IS %08d" % (particle_id)
	
	# Immediately delete any files that are no longer needed
	# This isn't just to keep the file system tidy---if we don't do this,
	# we could wind up dealing with literally *millions* of leftovers
	# If the exit code was unhappy, then we might want to keep the files
	if prism_exit_code == 0 and DELETE_FILES:
		os.remove(file_names["initial composition"])
		os.remove(file_names["trajectory"])
		os.remove(file_names["control"])
		os.remove(file_names["prism stdout"])
		os.remove(file_names["prism stderr"])

def write_fake_output_file(particle_id, translator, output_dir):
	"""To avoid running PRISM on cold particles, copy the initial composition
	of the particle straight to a fake PRISM output file. Return the name of
	the new composition file that we created.
	"""
	
	# The PRISM input and output files both have the same format
	# Several ASCII lines, each listing (Z, A, X) for one species
	fcomposition_file_name = os.path.join(output_dir, "%s_%08d.dat" \
		% (FCOMPOSITION_PREFIX, particle_id))
	translator.write_initial_composition_file(particle_id, \
		fcomposition_file_name)
	
	return fcomposition_file_name

def consolidate(particle_id, fcomposition_file_name, data_file_name):
	"""Consolidate the latest final composition file from PRISM into a big
	data file that consolidates everything we've done so far. Open the data
	file in append mode, then append the particle id and all the data from the
	final composition file.
	"""
	
	# Open final composition file for reading and data file for appending
	fcomposition_file = open(fcomposition_file_name, "r")
	data_file = open(data_file_name, "a")
	
	# Write the particle id to the data file, then write the composition lines
	data_file.write("%s\n" % (particle_id))
	for line in fcomposition_file:
		stripped = line.strip()
		if stripped != "":
			data_file.write(" ".join(stripped.split()) + "\n")
	
	# Close the files and delete the unneeded final composition file
	fcomposition_file.close()
	data_file.close()
	if DELETE_FILES:
		os.remove(fcomposition_file_name)

# This file is organized in the style of C--define functions, then call main
main()

