#!/usr/bin/env python
# Python wrapper program for parallelization of PRISM (part 1 of 3)
# Goal 1: produce post-processed yields for every particle from an SNSPH run
# Goal 2: avoid having to submit 950k sbatch jobs--that way lies madness
# Goal 3: organize all the output from PRISM after the processing is done
# This script is the setup program for submitting a team of teammate jobs

# Last modified 16 May 2020 by Greg Vance

import sys
import os.path

import json
import subprocess

def main():
	
	info = get_dispatch_info()
	
	# Process the inputs a little bit and store them in new variables
	sdf_dir = os.path.abspath(info["sdf dir"])
	output_dir = os.path.abspath(info["output_dir"])
	team_size = int(info["team size"])
	temp_cut = float(info["temp cut"])
	del info
	
	# Prepare files for each teammate job that needs dispatching
	job_script_file_names = list()
	for team_rank in xrange(team_size):
		
		# Write a little JSON control file for each teammate job
		inputs = {
			"sdf dir": sdf_dir,
			"output dir": output_dir,
			"team size": team_size,
			"team rank": team_rank,
			"temp cut": temp_cut}
		json_file_name = os.path.join(output_dir, "teammate_%04d.json" \
			% (team_rank))
		with open(json_file_name, "w") as json_file:
			json_file.dump(inputs, json_file)
		
		# Write an sbatch job script for each teammate job
		job_script_file_name = write_job_script(output_dir, json_file_name,
			team_rank)
		job_script_file_names.append(job_script_file_name)
	
	# Submit all the teammate job scripts using sbatch subprocesses
	for job_script_file_name in job_script_file_names:
		sbatch_submit(job_script_file_name)

def get_dispatch_info():
	"""Parse this script's single command line argument, read in the dispatch
	info from the specified JSON file, run a few quick checks, and then return
	the dict of dispatch info.
	"""
	
	# Check for the one command line argument that we want to see
	# We're only expecting the name of a small input JSON file
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
		raise ValueError("team size must be a positive int")
	if info["temp cut"] < 0.0:
		raise ValueError("bad temperature cutoff at %g K" \
			% (info["temp cut"]))
	
	return info

def write_job_script(output_dir, json_file_name, team_rank):
	"""
	"""
	
	script = [
		"#!/bin/bash",
		"",
		"#SBATCH -n 1",
		""
		]
	
	job_script_file_name = os.path.join(output_dir, "teammate_%04d.sh" \
		% (team_rank))
	with open(job_script_file_name, "w") as job_script_file:
		job_script_file.write("\n".join(script))
	
	return job_script_file_name

def sbatch_submit(job_script_file_name, print_command=True):
	"""Use a subprocess to submit a job script to the cluster with sbatch.
	Return the exit code from the sbatch command
	"""
	
	# Assemble the sbatch command that will submit the job script
	sbatch_command = ["sbatch", job_script_file_name]
	
	# Print the command right before running it if the flag is set
	if print_command:
		print "$ " + " ".join(sbatch_command)
	
	# Run the command and return its exit code
	exit_code = subprocess.call(sbatch_command)
	return exit_code

main()

