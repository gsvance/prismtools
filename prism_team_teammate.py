#!/usr/bin/env python
# 

# Last modified 16 May 2020 by Greg Vance

import sys
import os.path
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
DELETE_FILES = False  # whether to clean up the input files or not

def main():
	
	# Start timer and print welcome messages
	t1 = time.time()
	print "<< PRISM Parallelization Teammate Script >>"
	print "Beginning setup procedure...\n"
	
	parameters = read_input_parameters()
	
	# Extract the input parameters into shorter variable names
	sdf_dir = parameters["sdf dir"]
	output_dir = parameters["output dir"]
	team_size = parameters["team size"]
	team_rank = parameters["team rank"]
	temp_cut = parameters["temp cut"]
	del parameters
	
	# Print our team status as we understand it
	print "This job is assigned to run team rank %d of team size %d\n." \
		% (team_rank, team_size)
	
	# Set up the data translator object... this can take a while
	print "Initializing SDF-to_PRISM data translator object..."
	translator = s2p.SdfToPrismTranslator(sdf_dir)
	print "Translator object is constructed and ready to proceed."
	
	# Read the array of particle ids from the translator
	particle_ids = translator.get_particle_ids()
	print "SDF directory contains data for %d particles." (particle_ids.size)
	
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
	
	# Loop over all the particles that were assigned to this teammate
	# THE RULE: We are responsible for every index in the particle id array
	# that is congruent to our team_rank (modulo team_size).
	# Example: team_rank = 6, team_size = 10
	# We are responsible for the particles at indices 6, 16, 26, 36, 46, ...
	print "Starting main loop over assigned particle ids...\n"
	for particle_index in xrange(team_rank, particle_ids.size, team_size):
		
		particle_id = particle_ids[particle_index]
		print "Now considering particle id %08d at index %d." \
			% (particle_id, particle_index)
		
		# If this particle needs post-processing, then go ahead and run PRISM
		if particle_is_hot[particle_index]:
			print "This particle will require processing by PRISM."
			print "Preparing files needed by PRISM..."
			file_names = prepare_prism_files(particle_id, translator, \
				output_dir)
			print "Files ready. Starting up PRISM subprocess..."
			run_prism_once(particle_id, file_names)
			print "PRISM has finished running."
		
		# If the particle does NOT need post-processing, then don't run PRISM
		else:
			print "This particle will NOT require any processing by PRISM."
			print "Copying initial composition directly to output file..."
			write_fake_output_file()
		
		# consolidate the latest output file into a big output file
		consolidate()
		
		print
	
	# The loop is done, so just finish up
	print "Main loop is complete."
	t3 = time.time()
	t_loop = t3 - t2
	t_total = t3 - t1
	print "The main loop took %.3f minutes to run." % (t_loop / 60.)
	print "The entire program with setup took %.3f mniutes to run." \
		% (t_total / 60.)
	print
	
	print "Execution is now complete. Bye!"

def read_input_parameters():
	"""Parse the script's command line argument, read the input parameters
	from the specified JSON file, run some quick checks, and return the dict
	of parameters.
	"""
	
	# Check the number of command line arguments that were received
	# We're only expecting the name of a small input JSON file
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
	
	# Sanity-check a few of the input values before starting
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
	Details about their formats can be found in the PRISM manual PDF."""
	
	# Create this particle's initial composition ASCII file for PRISM
	icomposition_file_name = os.path.join(output_dir, "%s_%08d.dat" \
		% (ICOMPOSITION_PREFIX, particle_id))
	translator.write_initial_composition_file(particle_id,
		icomposition_file_name)
	
	# Create this particle's trajectory ASCII file for PRISM
	trajectory_file_name = os.path.join(output_dir, "%s_%08d.dat" \
		% (TRAJECTORY_PREFIX, particle_id))
	translator.write_trajectory_file(particle_id, trajectory_file_name)
	
	# Create the PRISM control file for processing this particle
	control_file_name = os.path.join(output_dir, "%s_%08d.json" \
		% (CONTROL_PREFIX, particle_id))
	fcomposition_file_name = os.path.join(output_dir, "%s_%08d.dat" \
		% (FCOMPOSITION_PREFIX, particle_id))
	write_control_file(control_file_name, icomposition_file_name,
		trajectory_file_name, fcomposition_file_name)
	
	# Set file names to collect PRISM's stdout and stderr
	prism_stdout_file_name = os.path.join(output_dir, "%s_%08d.out" \
		% (PRISM_STDOUT_PREFIX, particle_id))
	prism_stderr_file_name = os.path.join(output_dir, "%s_%08d.err" \
		% (PRISM_STDERR_PREFIX, particle_id))
	
	# Return all the file names so they don't need to be reconstructed later
	file_names = dict()
	file_names["initial composition"] = icomposition_file_name
	file_names["trajectory"] = trajectory_file_name
	file_names["control"] = control_file_name
	file_names["prism stdout"] = prism_stdout_file_name
	file_names["prism stderr"] = prism_stderr_file_name
	return file_names

def write_control_file(control_file_name, icomposition_file_name,
	trajectory_file_name, fcomposition_file_name, stop_temp=None):
	"""Write the JSON control file for a single particle using the requested
	file names and PRISM stop temperature (in kelvin). This mostly involves
	reading the deafult control file, making a few modifications, and then
	writing the modified data to a new JSON file."""
	
	# Read the JSON data from the default control file that came with PRISM
	with open(DEFAULT_CONTROL, "r") as default_control_file:
		control = json.load(default_control_file)
	
	# For now, use the standard nuclear data files
	#control["nuclear"]["datasets"] = :DEFAULT:
	
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
	#control["network"]["extent"] = :DEFAULT:
	
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
	"""
	"""
	
	# assemble a prism command
	# run prism with subprocess
	# check that prism was happy
	# clean up unneeded files

def	prism_minion_main():

	# Relocate to PRISM's install directory---it likes to run from there
	# Not whether each process can have a different working directory...
	os.chdir(PRISM_DIR)
	
	# Send your rank to the administrator to request your first assignment
	comm.send(rank, dest=admin_rank, tag=REQUEST_TAG)
	order = comm.recv(source=admin_rank, tag=ASSIGN_TAG)
	
	# Loop until the administrator sends the terminate signal
	while order != TERMINATE_SIGNAL:
		
		# The order will be a dict of file names if not the terminate_signal
		file_names = order
		
		# Open two file streams to catch the reporting from PRISM
		stdout_file = open(file_names["prism stdout"], "w")
		stderr_file = open(file_names["prism stderr"], "w")
		
		# Run a PRISM subprocess using the files provided to you
		command = ["./prism", "-c", file_names["control"]]
		prism_exit_code = subprocess.call(command, stdout=stdout_file,
			stderr=stderr_file, shell=True)
		
		# Close the PRISM output files
		stdout_file.close()
		stderr_file.close()
		
		# Draw attention to any unhappy exit codes from PRISM
		if prism_exit_code != 0:
			print "ALERT (PROCESS %d): PRISM RETURNED BAD EXIT CODE %d" \
				% (rank, prism_exit_code)
			print "ALERT (PROCESS %d): OFFENDING PARTICLE ID: %08d" \
				% (rank, file_names["particle id"])
		
		# Immediately delete any input files that are no longer needed
		# This isn't just to keep the file system tidy---if we don't do this,
		# we could wind up dealing with literally *millions* of leftover files
		# If the exit code was unhappy, then we might want to keep the files
		if prism_exit_code == 0 and DELETE_FILES:
			os.remove(file_names["initial composition"])
			os.remove(file_names["trajectory"])
			os.remove(file_names["control"])
			os.remove(file_names["prism stdout"])
			os.remove(file_names["prism stderr"])
		
		# Done! Request another assignment from the administrator
		comm.send(rank, dest=admin_rank, tag=REQUEST_TAG)
		order = comm.recv(source=admin_rank, tag=ASSIGN_TAG)
	
	# Print a final sign-off message when termination occurs
	print "  Process %d is terminating." % (rank)

# This file is organized C-style: define a bunch of functions, then call main
main()

