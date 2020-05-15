#!/usr/bin/env python


command line gives little json file

read little json file with inputs

check inputs for sanity

set up translator object

get ready with the particle ids and the temperature cut

report time?

loop over all assigned particles
	
	if hot
		
		write all prism input files
		
		run prism with subprocess
		
		check that prism was happy
		
		clean up unneeded files
	
	if not hot
		
		ignore or write a fake output file with icomp
	
	consolidate the output file into a big output file

clean up (and report times?)

