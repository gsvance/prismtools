#!/usr/bin/env python
# Python master script for administering a huge number of PRISM trajectories
# Goal: produce postprocessed yields for every particle from an SNSPH run
# Usage:
#     ./prism_master.py /path/to/sdf/dir /path/to/output/dir

# Last modified 6 May 2020 by Greg Vance

import sys
import os
import time
import json
import subprocess
import sdf_to_prism_translation as s2p

# Limit the number of jobs this script is allowed to administer simultaneously
MAX_JOBS = 500
# Amount of time to wait in between checking the number of jobs
JOB_CHECK_WAIT = 60. # seconds

# Standard output file prefixes for each type of file that needs naming
ICOMPOSITION_PREFIX = "icomp"
TRAJECTORY_PREFIX = "traj"
CONTROL_PREFIX = "ctrl"
PRISM_OUTPUT_PREFIX = "fcomp"
JOB_SCRIPT_PREFIX = "prism"
JOB_STDOUT_PREFIX = "part"
JOB_STDERR_PREFIX = "part"
# Path to the deafult control.json file that came with PRISM
DEFAULT_CONTROL_JSON = "/home/gsvance/prism/prism-1.5.0/input/control.json"

def main():
	
	print "Python Master Script for PRISM Parallelization"
	t1 = time.time()
	print "Beginning setup..."
	print
	
	# Parse and sanity-check the two required command line arguments
	print "Checking command line arguments..."
	if len(sys.argv) != 3:
		raise TypeError("%s takes two command line arguments, see py file" \
			% (sys.argv[0]))
	sdf_dir, output_dir = tuple(sys.argv[1:])
	if not os.path.exists(sdf_dir):
		raise ValueError("SDF directory does not exist")
	if not os.path.exists(output_dir):
		raise ValueError("output directory does not exist")
	print "Parsed input SDF directory path: %s" % (sdf_dir)
	print "Parsed output directory path: %s" % (output_dir)
	print
	
	# Get the translator set up in the SDF directory... this can take a while
	print "Setting up SDF-to-PRISM translator object..."
	translator = s2p.SdfToPrismTranslator(sdf_dir)
	particle_id_list = translator.get_particle_ids()
	print "Translator is ready to proceed."
	print
	
	# Report the amount of time that the setup consumed
	print "Setup complete."
	t2 = time.time()
	t_setup = t2 - t1
	print "Setup took %.3f minutes to run." % (t_setup / 60.)
	print
	
	# Run the main loop over every particle id present in the SDF directory
	print "Starting main loop over particle ids in SDF directory..."
	print
	for particle_id in particle_id_list:
		
		# Check the number of jobs I have running or queued right now
		# If there are too many, then wait until there are fewer
		ta = time.time()
		while check_number_of_jobs() >= MAX_JOBS:
			time.sleep(JOB_CHECK_WAIT)
		tb = time.time()
		t_wait = tb - ta
		print "Waited %.0f minutes due to number of jobs" % (t_wait / 60.)
		print
		
		print "Now setting up for particle id %08d..." % (particle_id)
		
		# Create the initial composition file for this particle id
		icomposition_file_name = os.path.join(output_dir, "%s_%08d.dat" \
			% (ICOMPOSITION_PREFIX, particle_id))
		translator.write_initial_composition_file(particle_id,
			icomposition_file_name)
		print "  Wrote icomp file: %s" % (icomposition_file_name)
		
		# Create the trajectory file for this particle id
		trajectory_file_name = os.path.join(output_dir, "%s_%08d.dat" \
			% (TRAJECTORY_PREFIX, particle_id))
		translator.write_trajectory_file(particle_id, trajectory_file_name)
		print "Wrote traj file: %s" % (trajectory_file_name)
		
		# Create the PRISM control file for this particle id
		control_file_name = os.path.join(output_dir, "%s_%08d.json" \
			% (CONTROL_PREFIX, particle_id)
		write_control_file(particle_id, control_file_name)
		print "Wrote ctrl file: %s" % (control_file_name)
		
		# Assemble the command to run PRISM for this particle
		prism_command_list = list()
		prism_command_list.append("prism")
		prism_command_list.extend(["-c", control_file_name])
		prism_command = " ".join(prism_command_list)
		
		# Assemble a job script around the PRISM command
		######
		
		# submit the prism command job script
		# ensure that the submission happened correctly
		
		print 
	
	print "Main loop complete."
	t3 = time.time()
	t_loop = t3 - t2
	t_total = t3 - t1
	print "Main loop took %.3f minutes to run." % (t_loop / 60.)
	print "Setup and main loop took %.3f minutes total." % (r_total / 60.)
	print
	
	print "Execution complete. Bye!"

def check_number_of_jobs():
	"""Use a subprocess command line call to find out how many jobs I have
	running or queued on the cluster right now. Return the number of jobs."""
	
	raise NotImplementedError

def write_control_file(particle_id, file_name):
	"""Write the JSON control file for a given particle id using the requested
	file name. Mostly involves reading the deafult control file, making some
	modifications to that, and then writing a new JSON file."""
	
	raise NotImplementedError

main()
