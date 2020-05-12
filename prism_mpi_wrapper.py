# Python MPI wrapper for administering a big bunch of parallel PRISM threads
# Goal 1: produce postprocessed yields for every particle from an SNSPH run
# Goal 2: avoid submitting 950k sbatch jobs---that way lies madness
# Example usage:
#     module purge
#     module load python/2.7.9
#     mpiexec -n 50 python prism_mpi_wrapper.py sdf/dir/path/ output/dir/path/

# ESSENTIAL INFORMATON: THIS SCRIPT WILL ONLY RUN PROPERLY WITH THE CORRECT
# PYTHON MODULE VERSION LOADED: module load python/2.7.9
# THIS HAS SOMETHING TO DO WITH THE INSTALL DIRECTORY OF mpi4py ON AGAVE...

# Last modified 12 May 2020 by Greg Vance

from mpi4py import MPI  # Provides Python bindings for MPI

import sys
import os
import time

import json
import subprocess

import sdf_to_prism_translation as s2p

# Rank of the MPI process that will administrate all other processes
ADMINISTRATOR_RANK = 0
# Tags for MPI communications (TOTALLY RANDOMLY CHOSEN NUMBERS)
REPORT_TAG = 15
REQUEST_TAG = 40
ASSIGN_TAG = 22
# Signal code for the minion processes to terminate operations
TERMINATE_SIGNAL = -9999

# Temperature cutoff to use when deciding what particles to postprocess
# Also useful as a stop condition when running PRISM
TEMPERATURE_CUTOFF = 1E8 # K

# Standard output file prefixes for each type of file that needs naming
ICOMPOSITION_PREFIX = "icomp"
TRAJECTORY_PREFIX = "traj"
CONTROL_PREFIX = "ctrl"
FCOMPOSITION_PREFIX = "fcomp"
PRISM_STDOUT_PREFIX = "prism"
PRISM_STDERR_PREFIX = "prism"

# Path to the PRISM install directory
PRISM_DIR = "/home/gsvance/prism/prism-1.5.0/"
# Path to the deafult control.json file that came with PRISM
DEFAULT_CONTROL_JSON = os.path.join(PRISM_DIR, "input/control.json")

def main():
	"""Ultra-simple main function to check the rank of each MPI process and
	then direct the process to run one of the other two "main" functions."""
	
	# Find the rank of the MPI process
	comm = MPI.COMM_WORLD
	rank = comm.Get_rank()
	
	# One process is the administrator---it will run the show
	# It manages the input data and farms out PRISM jobs to other processes
	if rank == ADMINISTRATOR_RANK:
		administrator_main()
	
	# All the processes are the legions of PRISM job minions
	# Their task is to run PRISM jobs as directed by the administrator
	else:
		prism_minion_main()

def administrator_main():
	"""Main function for the administrator process to run. The administrator's
	responsibilities are as follows:
	  - Print informative messages to stdout for the user
	  - Record runtime statistics for this script
	  - Read and store all the particle data from SNSPH
	  - Use a translator object to produce PRISM input files
	  - Send particle ids to minion processes who will actually run PRISM"""
	
	# Set the basic MPI variables
	comm = MPI.COMM_WORLD
	rank = comm.Get_rank()
	size = comm.Get_size()
	
	print "Administrator process %d was initialized successfully." % (rank)
	print "<< Python Wrapper Script for PRISM Parallelization >>"
	t1 = time.time()
	print "Beginning administrator setup..."
	print
	
	# Check that we have a least one other process to boss around
	if size < 2:
		raise TypeError("two or more processes are required")
	print "Detected %d other processes for running PRISM jobs." % (size - 1)
	
	# Parse and sanity-check the two required command line arguments
	if len(sys.argv) != 3:
		raise TypeError("two command line arguments required, see .py file")
	sdf_dir, output_dir = tuple(sys.argv[1:])
	if not os.path.exists(sdf_dir):
		raise ValueError("SDF directory does not exist")
	if not os.path.exists(output_dir):
		raise ValueError("output directory does not exist")
	print "Parsed input SDF directory path: %s" % (sdf_dir)
	print "Parsed output directory path: %s" % (output_dir)
	
	# Get the translator set up in the SDF directory... this can take a while
	print "Initializing SDF-to-PRISM data translator object..."
	translator = s2p.SdfToPrismTranslator(sdf_dir)
	print "Translator object is prepared and ready to proceed."
	
	# Read the full list of particle ids from the translator
	particle_ids = translator.get_particle_ids()
	n_particles = len(particle_ids)
	print "SDF directory contains data for %d particles." % (n_particles)
	
	# Impose temperature cutoff to limit the number of PRISM runs
	# This filtering can also take a while to happen....
	print "Imposing temperature cut at %.3E kelvin..." % (TEMPERATURE_CUTOFF)
	hot_particle_ids \
		= translator.get_particle_ids_with_temp_cutoff(TEMPERATURE_CUTOFF)
	n_hot_particles = len(hot_particle_ids)
	print "Found %d particles that will need to be processed by PRISM." \
		% (n_hot_particles)
	
	# Signal all other processes to report in by sending them our rank
	print "Signaling all other processes to report in..."
	for other_rank in xrange(size):
		if other_rank != rank:
			comm.send(rank, dest=other_rank, tag=REPORT_TAG)
	
	# Report total the amount of time that the setup consumed
	print
	print "Administrator setup complete."
	t2 = time.time()
	t_setup = t2 - t1
	print "Setup took %.3f minutes to finish." % (t_setup / 60.)
	
	##############################################
	##### ADMINISTRATOR SETUP CODE ENDS HERE #####
	##############################################
	
	# Run the main loop over every particle id that needs processing
	print 
	print "Starting main loop over particle ids..."
	for hot_particle_id in hot_particle_ids:
		
		# Prepare the PRISM input files for this particle id
		print
		print "Preparing PRISM input files for particle id %08d..." \
			% (hot_particle_id)
		file_names = prepare_prism_input_files(hot_particle_id, translator,
			sdf_dir, output_dir)
		
		# Wait for another process to request a job and then give it to them
		print "Files prepared. Waiting for an available process..."
		assign_rank = comm.recv(source=MPI.ANY_SOURCE, tag=REQUEST_TAG)
		comm.send(file_names, dest=assign_rank, tag=ASSIGN_TAG)
		print "Particle id %08d was assigned to process %d." \
			% (hot_particle_id, assign_rank)
	
	# Wait for the remaining processes to finish, then terminate them
	print
	print "Main loop complete. Waiting for processes to finish..."
	processes_remaining = size - 1
	while processes_remaining > 0:
		process_rank = comm.recv(source=MPI.ANY_SOURCE, tag=REQUEST_TAG)
		comm.send(TERMINATE_SIGNAL, dest=process_rank, tag=ASSIGN_TAG)
		processes_remaining -= 1
	
	print
	print "All non-administrator processes have terminated."
	t3 = time.time()
	t_loop = t3 - t2
	t_total = t3 - t1
	print "PRISM processing took %.3f minutes to run." % (t_loop / 60.)
	print "Setup and processing took %.3f minutes total." % (t_total / 60.)
	
	print
	print "Execution is complete. Bye!"
	print "Process %d terminated." % (rank)

def prepare_prism_input_files(particle_id, translator, sdf_dir, output_dir):
	"""Given a particle id, prepare all of the necessary input files that
	PRISM needs in order to postprocess yields for that particle. Return a
	dict of paths to the newly created files. Produces the following files:
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
		% (FCOMPOSITION_PREFIX, particle_id)
	write_control_file(control_file_name, icomposition_file_name,
		trajectory_file_name, fcomposition_file_name, TEMPERATURE_CUTOFF)
	
	# Set file names to collect PRISM's stdout and stderr
	prism_stdout_file_name = os.path.join(output_dir, "%s_%08d.out" \
		% (PRISM_STDOUT_PREFIX, particle_id))
	prism_stderr_file_name = os.path.join(output_dir, "%s_%08d.err" \
		% (PRISM_STDERR_PREFIX, particle_id))
	
	# Return all the file names so they don't need to be reconstructed later
	# The administrator can send them straight along to another process
	file_names = dict()
	file_names["initial composition"] = icomposition_file_name
	file_names["trajectory"] = trajectory_file_name
	file_names["control"] = control_file_name
	file_names["prism stdout"] = prism_stdout_file_name
	file_names["prism stderr"] = prism_stderr_file_name
	return file_names

def write_control_file(control_file_name, icomposition_file_name,
	trajectory_file_name, fcomposition_file_name, stop_temp,):
	"""Write the JSON control file for a single particle using the requested
	file names and PRISM stop temperature (in kelvin). This mostly involves
	reading the deafult control file, making a few modifications, and then
	writing the modified data to a new JSON file."""
	
	# Read the JSON data from the default control file that came with PRISM
	with open(DEFAULT_CONTROL_JSON, "r") as default_control_file:
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
	stop_T9 = stop_temp * 1E-9  # Convert stop temperature to GK
	control["network"]["stop"] = {"T9": stop_T9, "comment":
		"Stops calculation when trajectory reaches %.3E GK" % (stop_T9)}
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

def	prism_minion_main():
	"""Main function for all the prism minion processes to run. The minions'
	responsibilities are as follows:
	  - Message the administrator to request assignments
	  - Run the actual PRISM subprocesses
	  - Clean up any leftover files that might overwhelm the file system"""
	
	# Set the basic MPI variables
	comm = MPI.COMM_WORLD
	rank = comm.Get_rank()
	size = comm.Get_size()
	
	# Wait to receive the administrator's initial order to report in
	# This might take some time since the administrator has to set up first
	admin_rank = comm.recv(source=MPI.ANY_SOURCE, tag=REPORT_TAG)
	print "  Process %d, standing by." % (rank)
	
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
			stderr=stderr_file)
		
		# Close the PRISM output files
		stdout_file.close()
		stderr_file.close()
		
		# Draw attention to any unhappy exit codes from PRISM
		if prism_exit_code != 0:
			print "ALERT (PROCESS %d): PRISM RETURNED BAD EXIT CODE %d" \
				% (rank, prism_exit_code)
			print "ALERT (PROCESS %d): OFFENDING FILES DICT: %s" \
				% (rank, file_names)
		
		# If the exit code was happy, then we don't need to keep the reports
		if prism_exit_code == 0:
			os.remove(file_names["prism stdout"])
			os.remove(file_names["prism stderr"])
		
		# Immediately delete any input files that are no longer needed
		# This isn't just to keep the file system tidy---if we don't do this,
		# we could wind up dealing with literally *millions* of leftover files
		os.remove(file_names["initial composition"])
		os.remove(file_names["trajectory"])
		os.remove(file_names["control"])
		
		# Done! Request another assignment from the administrator
		comm.send(rank, dest=admin_rank, tag=REQUEST_TAG)
		order = comm.recv(source=admin_rank, tag=ASSIGN_TAG)
	
	# Print a final sign-off message when termination occurs
	print "  Process %d is terminating." % (rank)

# This file is organized C-style: define a bunch of functions, then call main
main()

