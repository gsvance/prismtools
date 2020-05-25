#!/usr/bin/env python
# Python wrapper program for parallelization of PRISM (part 1 of 3)
# Goal 1: produce post-processed yields for every particle from an SNSPH run
# Goal 2: avoid having to submit 950k sbatch jobs--that way lies madness
# Goal 3: organize all the output from PRISM after the processing is done
# This script is the setup program for submitting a team of teammate jobs
# Run this script interactively!

# Last modified 25 May 2020 by Greg Vance

import sys
import os.path

import json
import subprocess

# Full path to the Python PRISM team teammate script (part 2 of 3)
# This is the file that needs to be executed in all the jobs we will submit
PRISM_TEAM_TEAMMATE = "/home/gsvance/prismtools/prism_team_teammate.py"

def main():
	
	info = get_dispatch_info()
	
	# Process the inputs a little bit and store them in new variables
	sdf_dir = os.path.abspath(info["sdf dir"])
	output_dir = os.path.abspath(info["output dir"])
	team_size = int(info["team size"])
	job_time = str(info["job time"])
	temp_cut = float(info["temp cut"])
	del info
	
	# Prepare the files for each teammate job that needs dispatching
	job_script_file_names = list()
	for team_rank in xrange(team_size):
		
		# Write little JSON control files to keep each teammate job organized
		json_contents = {
			"sdf dir": sdf_dir,
			"output dir": output_dir,
			"team size": team_size,
			"team rank": team_rank,
			"temp cut": temp_cut
		}
		json_file_name = os.path.join(output_dir, "teammate_%04d.json" \
			% (team_rank))
		with open(json_file_name, "w") as json_file:
			json.dump(json_contents, json_file, indent=2,
				separators=(",", ": "))
		
		# Write an sbatch job script for submitting each teammate job
		job_script_file_name = write_job_script(output_dir, json_file_name,
			team_rank, job_time)
		job_script_file_names.append(job_script_file_name)
		
		print "Wrote JSON and job script files for team rank %d of %d." \
			% (team_rank, team_size)
	
	# This file should be run interactively, not as a cluster job
	# Ask the user whether we want to actually submit the job scripts
	if not ask_user("Submit sbatch job scripts to the cluster?"):
		print "Aborting..."
		return
	
	# Submit all the teammate job scripts using sbatch subprocesses
	# Keep track of the exit codes to make sure this goes smoothly
	for job_script_file_name in job_script_file_names:
		exit_code = sbatch_submit(job_script_file_name)
		while exit_code != 0:
			print "Job submit failed: exit code was %d" % (exit_code)
			if not ask_user("Try to submit again?"):
				print "Aborting..."
				return
			exit_code = sbatch_submit(job_script_file_name)

def get_dispatch_info():
	"""Parse this script's single command line argument, read in the dispatch
	info from the specified JSON file, run a few quick checks, and then return
	the dict of dispatch info.
	"""
	
	# Check for the one command line argument that we want to see
	# We are only expecting to get the name of a small input JSON file
	if len(sys.argv) != 2:
		raise TypeError("one command line argument is required: " \
			+ "name of JSON dispatch info file")
	json_file_name = sys.argv[1]
	
	# Read the little JSON file and load the dict of dispatch info
	with open(json_file_name, "r") as json_file:
		info = json.load(json_file)
	
	# Pretty-print the input values we read from the file
	print "Read dispatch info from input file:\n%s" % (json_file_name)
	print json.dumps(info, indent=2, separators=(",", ": "))
	print
	
	# Sanity-check a few of the input values before returning them
	if not os.path.exists(info["sdf dir"]):
		raise ValueError("SDF directory does not exist: %s" \
			% (info["sdf dir"]))
	if not os.path.exists(info["output dir"]):
		raise ValueError("output directory does not exist: %s" \
			% (info["output dir"]))
	if info["team size"] <= 0:
		raise ValueError("team size must be a positive integer")
	if not set(info["job time"]).issubset(set("0123456789-:")):
		raise ValueError("invalid slurm job time: %s" % (info["job time"]))
	if info["temp cut"] < 0.0:
		raise ValueError("bad temperature cutoff at %g K" \
			% (info["temp cut"]))
	
	return info

def write_job_script(output_dir, json_file_name, team_rank, job_time):
	"""Write an sbatch job script for submitting a single teammate job and
	return the path to the newly created file. The job's stdout and stderr
	files will be saved to the provided output directory. The teammate Python
	script will be invoked with the name of the given JSON file as input when
	the job runs. For file naming, please provide the rank of this teammate
	job in the team. The job script will request the specified amount of wall
	time from slurm, which should be specified as a slurm-compatible string,
	e.g., "2-04:32" indicates 2 days, 4 hours, and 32 minutes.
	"""
	
	# Put the job's stdout and stderr files in the output directory
	# Name them based on the job's team rank and slurm job id
	stdout_file_name = os.path.join(output_dir, "teammate_%04d.%%j.out" \
		% (team_rank))
	stderr_file_name = os.path.join(output_dir, "teammate_%04d.%%j.err" \
		% (team_rank))
	
	# Contruct the contents of the script file using a list of strings
	script_contents = [
		"#!/bin/bash",
		"",
		"#SBATCH -n 1",
		"#SBATCH -t %s" % (job_time),
		"#SBATCH -o %s" % (stdout_file_name),
		"#SBATCH -e %s" % (stderr_file_name),
		"#SBATCH --mail-type=ALL",
		"#SBATCH --mail-user=gsvance@asu.edu",
		"",
		"module purge",
		"module load numpy/python-2x",  # for running Numpy code
		"module load intel/2017x",  # for running PRISM executable
		"",
		"python %s %s" % (PRISM_TEAM_TEAMMATE, json_file_name),
		""
	]
	
	# Join the script contents together and write them to a file
	job_script_file_name = os.path.join(output_dir, "teammate_%04d.sh" \
		% (team_rank))
	with open(job_script_file_name, "w") as job_script_file:
		job_script_file.write("\n".join(script_contents))
	
	return job_script_file_name

def ask_user(question):
	"""Prompt the user with a printed question, read their response, and
	return whether the reply was "yes."
	"""
	
	# Print the question, prompt for input, and then check the input
	response = raw_input(question + " [y/n] ").lower()
	return (response == "y" or response == "yes")

def sbatch_submit(job_script_file_name, print_command=True):
	"""Use a subprocess to submit a job script to the cluster with sbatch.
	Return the exit code from the sbatch command. As usual with exit codes
	from shell commands, 0 indicates success and nonzero indicates problems.
	"""
	
	# Assemble the sbatch command that will submit the job script
	sbatch_command = ["sbatch", job_script_file_name]
	
	# Print the command right before running it if the flag was set
	if print_command:
		print "$ " + " ".join(sbatch_command)
	
	# Run the command and return its exit code
	exit_code = subprocess.call(sbatch_command)
	return exit_code

main()

