#!/usr/bin/env python3

import os
import argparse
import time
import logging
import logging.handlers
from tqdm import tqdm
from cas13gRNAtor.utils import *
from cas13gRNAtor.cas13gRNAtor import *

logger = logging.getLogger(__name__)

DESCRIPTION = '''
Cas13gRNAvalidate Validates Previous gRNAs used against an MSA or non-MSA
Scores Shannon Entropy and Conservation Scores for each gRNA
Input: Aligned Sequences.
Output: Conservation and Entropy score in csv as supplementary data'
'''

def get_args():
	parser = argparse.ArgumentParser(
		prog='cas13gRNAvalidate',
		description=DESCRIPTION
	)
	parser.add_argument('-i', '--input', required=True, dest='gRNA',
						help="Takes in a list of gRNAs in a fasta file or a text file, list seperated by a newline")
	parser.add_argument('-m','--MSA', dest = "MSA", metavar="fasta", default = None, help="Sequences to align")
	parser.add_argument('-bc', '--bowtie-conservation', dest = "bowtie", default = None, help="Path to fasta file for Bowtie conservation calculation")
	parser.add_argument('-bo', '--bowtie-offtarget', dest = "offtarget", help="Path to fasta file for Bowtie conservation calculation")
	
	parser.add_argument('--mismatch', default=2, type=int,  choices=[0,1,2,3,4,5,6,7,8,9,10],
						help='Number of mismatches allowed for Bowtie (Default: %(default)s)')

	parser.add_argument('-t', '--threads', default=4, type=int, dest = 'CPU', help="Number of Core to use (Default: %(default)s)")
	parser.add_argument('-p', '--prefix', default = 'RfxCas13d_gRNA' , help='file name prefix to your file names (Default: %(default)s)')
	#parser.add_argument('-o', '--output',  default = '.', help='Default: Save files to current Directory')
	parser.add_argument('-d', '--debug',
						help='Print lots of debugging statements',
						action="store_const",dest="loglevel",const=logging.DEBUG,
						default=logging.INFO)
	'''
	parser.add_argument('--conservation-index', dest="conservation_index", default=None, help='Input an index for conservation calculation if you have!')
	parser.add_argument('--offtarget-index', dest="offtarget_index", default=None, help='Input an index for offtargets if you have!')
	'''
	parser.add_argument('--keep-tmp', action="store_false", dest='temp', help='Keep Temporary files (Default: Does not Keep)')
	parser.add_argument('--no-plot', action="store_false",dest="plot", help='Does not Plot graph')

	args = parser.parse_args()
	if not args.MSA:
		assert args.bowtie, ASSERT_MESSAGE_ALIGNMENT
	elif not args.bowtie:
		assert args.MSA, ASSERT_MESSAGE_ALIGNMENT
	
	
	return args

def open_gRNA_file(args):
	## Check gRNA file -> Text or Fasta 
	gRNA_file = args.gRNA
	crRNA_handle = args.prefix + '_4thQuartile_gRNA.fasta'
	crRNA_file = open(crRNA_handle, 'w')
	crRNA_count = 0
	gRNA_class_list = []
	with open(gRNA_file, 'r') as f: 
		for line in f:
			if line.startswith('>'): continue
			seq = line.strip('\n')
			seq = seq.upper()
			if len(seq) != LENGTH_OF_GRNA: logger.info(f"{seq} length != 23, Will Ignore!")
			for c in break_apart(seq):
				if not c in ALPHABET_NORMAL:
					logger.info(f"{seq} contains invalid Nucleotide Alphabet/s, Will Ignore!")
					continue
			else:
				crRNA_count += 1
				gRNA = Cas13gRNA(f"gRNA_{crRNA_count}", seq)
				gRNA_class_list.append(gRNA)
				crRNA_file.write(f">gRNA_{crRNA_count}\n{seq}\n")

	crRNA_file.close()
	return crRNA_handle, gRNA_class_list

def main():
	args = get_args()
	start_time = time.time()
	logging.basicConfig(level=args.loglevel, format=FORMAT)
	logger.info("cas13gRNAvalidate starting!!")
	temp_files_to_remove = []
	
	crRNA_handle, gRNA_class_list = open_gRNA_file(args) 

	offtarget_summary, temp_files_to_append = bowtie_main(args, crRNA_handle, mode = 'offtarget')
	gRNA_class_list = update_class(gRNA_class_list, offtarget_summary = offtarget_summary)
	temp_files_to_remove += temp_files_to_append

	gRNA_class_list, conservation_summary, temp_files_to_append, scoring_method, consensus_length = get_scores(args, gRNA_class_list, crRNA_handle = crRNA_handle, get_weighted_score = False)
	temp_files_to_remove += temp_files_to_append
	
	write_supplementary_data(args, gRNA_class_list)
	gRNA_freq = {k:0 for k in range(consensus_length)}
	gRNA_freq = write_all(args, gRNA_class_list, gRNA_freq = gRNA_freq, scoring_method = scoring_method, plot = args.plot, weighted = False)
	if scoring_method[0]:
		main_plot(args.prefix, gRNA_freq)
		if args.plot:
			logger.info("Plotting individual gRNAs!")
			if not os.path.isdir(args.prefix + '_best_gRNAs'): os.mkdir(args.prefix + '_best_gRNAs') 
			for g in tqdm(gRNA_class_list):
				if g.found: plot_highly_conserved_gRNA(g.id, g.c_score_list, g.e_score_list, g.seq, g.pos, args.prefix)

	if args.temp:
		for fname in temp_files_to_remove:
			os.remove(fname)
		logger.info("Temporary files removed!")
		
	end_time = time.time()
	total_time = round(end_time - start_time,3)
	logger.info(f"cas13gRNAvalidate has ended in {total_time} seconds!")


if __name__ == "__main__":
	main()
