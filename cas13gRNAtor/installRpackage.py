#!/usr/bin/env python3

import os
from cas13gRNAtor.utils import *

def main():
	software = [
	    "misc/bowtie",
	    "misc/bowtie-build",
	    "misc/RNAfold",
	    "misc/RNAhybrid",
	    "misc/RNAplfold",
	    "misc/RNAhyb.sh",
	]

	#for s in software:
	#	os.system(f"chmod +x {PYTHON_FILE}/{s}")

	check_Rscripts()
	check_Rscripts_tools()

	print ("Installing Rpackages!")
	MOD_NAME = "cas13gRNAtor"
	FASTA = f"{PYTHON_FILE}/misc/data/test.fa"
	os.system(" ".join(['Rscript', f'{PYTHON_FILE}/misc/RfxCas13d_GuideScoring.R', FASTA, f"{PYTHON_FILE}/misc/data/Cas13designGuidePredictorInput.csv", "true"]))
	print ("Finish Installing Rpackages!")

if __name__ == "__main__":
	main()