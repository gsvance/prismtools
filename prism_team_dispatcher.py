#!/usr/bin/env python
# 

# Last modified 15 May 2020 by Greg Vance

import os
import json

# Constants here


def main():
	# get inputs: sdf_dir, ouput_dir, team_size
	# maybe from prompting user and not command line?
	
	# Sanity check all of the inputs
	if not os.path.exists(sdf_dir):
		raise ValueError("")
	if not os.path.exists(output_dir):
		raise ValueError("")
	if team_size <= 0:
		raise ValueError("")
	
	# Process the inputs
	sdf_dir = os.path.abspath(sdf_dir)
	output_dir = os.path.abspath(output_dir)
	
	# write little json control files
	for team_rank in xrange(team_size):
		
		lil_json = {"sdf dir": sdf_dir, "output dir": output_dir,
			"team size": team_size, "team rank": team_rank, "temp cut": 1E+8}
		
		json_filename = ""
		with open(json_filename, "w") as json_f:
			json.dump(lil_json, json_f)
	
	# write little job scripts
	
	
	# sbatch submit job scripts
	

def sbatch_submit(job_script_file_name, print_command=True):
	"""Use a subprocess to submit a job script to the cluster with sbatch.
	Return the exit code from the sbatch command."""
	
	# Assemble the sbatch command that will submit the job script
	sbatch_command = ["sbatch", job_script_file_name]
	
	# Print the command right before running it if the flag is set
	if print_command:
		print "$ " + " ".join(sbatch_command)
	
	# Run the command and return its exit code
	exit_code = subprocess.call(sbatch_command)
	return exit_code

main()

