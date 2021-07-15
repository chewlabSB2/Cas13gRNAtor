#!/usr/bin/env python3

import argparse
import time
import logging
import logging.handlers
from cas13gRNAtor.utils import *

logger = logging.getLogger(__name__)

DESCRIPTION = """
cas13gRNAscore Multiple Sequence Alignment Entropy and Conservation scoring
Input: Aligned Sequences
Output: Conservation and Entropy score in csv as supplementary data
"""

def get_args():
	parser = argparse.ArgumentParser(
		prog='cas13gRNAscore',
		description=DESCRIPTION
	)
	parser.add_argument('-m','--MSA', dest = "MSA", required = True, metavar="fasta", help = "Multiple Sequence Alignment")
	parser.add_argument('-p', '--prefix', default='Supp_Data', metavar="file name",
						help='file name prefix to your file names (Default: %(default)s)')
	parser.add_argument('-d', '--debug',
						help='Print lots of debugging statements',
						action="store_const",dest="loglevel",const=logging.DEBUG,
						default=logging.INFO)
	parser.add_argument('--no-plot', action='store_false', dest="plot",
						help='Do not Plot graph (Default: %(default)s)')
	args = parser.parse_args()
	return args

def main():
	args = get_args()
	start_time = time.time()
	
	logging.basicConfig(level=args.loglevel, format=FORMAT)
	logger.info("cas13gRNAscore starting!!")
	
	aln_array, MSA = read_alignment(True, args.MSA)
	if not MSA: raise AlignmentError(f"{args.MSA} is not aligned or formatted wrongly, please check!")
	if len(aln_array) < 100: 
		logger.debug("Number of Sequences insufficient!!!")
		logger.debug("Conservation and Entropy Score might not be accurate")
	_, conservation, entropy = main_calc(aln_array, args.prefix)
	
	for i in [0, 10, 50, 100, 250]:
		ma_entropy = get_moving_average(entropy, i)
		ma_conservation = get_moving_average(conservation, i)
		plot_everything(ma_entropy, ma_conservation, args.prefix + '_' + str(i) + '_movingAvg')
	
	end_time = time.time()
	total_time = round(end_time - start_time,3)
	logger.info(f"cas13gRNAscore has ended in {total_time} seconds!")


if __name__ == "__main__":
	main()