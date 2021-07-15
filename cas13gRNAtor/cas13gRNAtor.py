#!/usr/bin/env python3

import os
import logging
import logging.handlers
import time
from cas13gRNAtor.utils import *

logger = logging.getLogger(__name__)

DESCRIPTION = '''
RfxCas13d guide RNA Efficacy Prediction, Conservation and Entropy Scoring and Selection
Based on Guide Efficacy, Shannon Entropy and Conservation Scores of gRNA 
Removes gRNA with offtargets and homopolymer
Plot The Shannon Entropy, Conservation and gRNA freq of the Genome
Plot The Shannon Entropy and Conservation Score for the best 100 gRNAs 
'''

#To-do
def get_args():
	parser = argparse.ArgumentParser(
		prog = 'cas13gRNAtor',
		description = DESCRIPTION
	)
	parser.add_argument('-r', '--reference', dest = 'reference', metavar="fasta", help="Reference for Guide Score Calculation")
	parser.add_argument('--guide-score', dest = 'scorefile', metavar="csv", help="Guide Score From RfxCas13d_GuideScoring.R")

	parser.add_argument('-m','--MSA', dest = "MSA", metavar="fasta", help="Sequences to align")
	parser.add_argument('-bc', '--bowtie-conservation', dest = "bowtie", help="Path to fasta file for Bowtie conservation calculation")
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
	parser.add_argument('--keep-homopolymer', action="store_false", dest="homopolymer", help='Keep gRNAs with homopolymer UUUU (Default: False)')
	parser.add_argument('--keep-tmp', action="store_false", dest='temp', help='Keep Temporary files (Default: Does not Keep)')
	parser.add_argument('--no-plot', action="store_false",dest="plot", help='Does not Plot graph')

	args = parser.parse_args()
	if not args.reference:
		assert args.scorefile, ASSERT_MESSAGE_SCORE
	elif not args.scorefile:
		assert args.reference, ASSERT_MESSAGE_SCORE

	if not args.MSA:
		assert args.bowtie, ASSERT_MESSAGE_ALIGNMENT
	elif not args.bowtie:
		assert args.MSA, ASSERT_MESSAGE_ALIGNMENT
			
	return args

def get_gRNAs(args):
	temp_files_to_remove = []
	if not args.scorefile:
		check_Rscripts()
		check_Rscripts_tools()

		logger.info("gRNA scoring and prediction starting!")
		#gRNA scoring from consensus sequence
		reference = list(SeqIO.parse(args.reference, "fasta"))[0]
		reference_len = len(reference)
		logger.info(f"Reference Length is {reference_len}!")
		cmd_gRNA = generate_gRNA_scores(args.reference)
		success_gRNA = run_shell_command(cmd_gRNA)
		if not success_gRNA:
			raise AlignmentError("Error during gRNA scoring")
		logger.info("gRNA scoring and prediction done!")
		cwd = os.path.abspath(os.getcwd())
		csv_file = get_csv_file(cwd, 'CasRxguides.csv')

		Rfiles = [f'{reference.id}_CasRxguides.fa', f'{reference.id}_CasRxguides.csv'] #, 'Consensus_Sequence_CasRxguides.pdf'
		temp_files_to_remove += Rfiles
	else:
		csv_file = args.scorefile

	crRNA_handle, gRNA_classes, max_guide_score = filter_gRNA(args, csv_file)
	temp_files_to_remove.append(crRNA_handle)

	return crRNA_handle, gRNA_classes, temp_files_to_remove, max_guide_score

def bowtie_main(args, crRNA_fasta, fasta_dict = {}, mode = 'offtarget'):
	if mode == 'offtarget' and not args.offtarget:
		return {}, []

	temp_files_to_remove = []
	bowtie_path, bowtie_build_path = check_bowtie(args)
	if mode == 'offtarget':
		_prefix = args.prefix  + '_offtarget' 
		bowtie_reference = args.offtarget
		samfile = args.prefix + '_offtarget.sam'
	elif mode == 'conservation':
		_prefix = args.prefix + '_conservation'
		bowtie_reference = args.bowtie
		samfile = args.prefix + '_conservation.sam'

	logger.info(f"Building Bowtie Index for {mode}! This might take a while")
	bowtie_build_cmdline, bowtie_index_prefix = bowtie_build_cmd(bowtie_reference, _prefix, bowtie_build_path, args.CPU)
	success_bowtie_build = run_shell_command(bowtie_build_cmdline)
	if not success_bowtie_build:
		raise AlignmentError("Error during bowtie-build")

	temp_files_to_remove.append(_prefix + '_output_build.txt')
	temp_files_to_remove.append(_prefix + '_error_build.txt')

	logger.info(f"Bowtie Starting for {mode}!")

	#crRNA_fasta = Reference Fasta File
	mismatch = args.mismatch if args.mismatch <= 3 else 3
	cmdline_bowtie = bowtie_cmd(bowtie_index_prefix, crRNA_fasta, _prefix, mismatch, bowtie_path, args.CPU)
	success_bowtie = run_shell_command(cmdline_bowtie)
	if not success_bowtie:
		raise AlignmentError("Error during offtarget filtering")

	bowtie_summary = bowtie_to_summary(samfile, fasta_dict, mode = mode)
	logger.info(f"Bowtie Done for {mode}!")
	temp_files_to_remove.append(samfile)

	onlyfiles = [f for f in os.listdir() if f.startswith(f"{_prefix}_alignment.index")]
	temp_files_to_remove += onlyfiles

	return bowtie_summary, temp_files_to_remove

def get_scores(args, gRNA_class_list, crRNA_handle = None, max_guide_score = 0, get_weighted_score = True):
	conservation_summary, temp_files_to_append, edit_gRNAs = {}, [], []
	do_bowtie = False
	consensus_length = 0
	MSA_ed, Bowtie_ed = False, False
	if args.MSA:
		aln_array, MSA = read_alignment(True, args.MSA)
		logger.debug(f"There are {aln_array} aligned sequences")
		if MSA:
			if len(aln_array) < 100: 
				logger.debug("Number of Sequences insufficient!!!")
				logger.debug("Conservation and Entropy Score might not be accurate")
			logger.info("Will Calculate Conservation and Entropy Score")
			consensus, conservation, entropy = main_calc(aln_array, args.prefix)

			logger.info("Will find reference gRNAs in the consensus sequence")
			gRNA_class_list = multi_search(args, consensus, gRNA_class_list, conservation, entropy, max_guide_score, get_weighted_score = get_weighted_score)
			consensus_length = len(consensus)
			MSA_ed = True
		else:
			logger.info("Will Bowtie to get Conservation Score") 
			do_bowtie = True

	if args.bowtie or do_bowtie:
		_, MSA = read_alignment(False, args.bowtie)
		if MSA:
			logger.debug("Genome sequences is aligned, Bowtie might be inaccurate!")
		fasta_dict = {k.id:0 for k in gRNA_class_list}
		conservation_summary, temp_files_to_append = bowtie_main(args, crRNA_handle, fasta_dict = fasta_dict, mode = 'conservation')
		gRNA_class_list = update_class(gRNA_class_list, conservation_summary = conservation_summary)
		Bowtie_ed = True
	
	return gRNA_class_list, conservation_summary, temp_files_to_append, (MSA_ed, Bowtie_ed), consensus_length

def main_plot(prefix, gRNA_freq):
	csv_file = prefix + '_position_score.csv'
	entropy, conservation = [], []
	with open(csv_file, 'r') as f:
		next(f) ## skip header
		for line in f:
			line = line.strip('\n')
			line = line.split('\t')
			entropy.append(float(line[1]))
			conservation.append(float(line[2]))
	
	for i in [0, 10, 50, 100, 250]:
		ma_entropy = get_moving_average(entropy, i)
		ma_conservation = get_moving_average(conservation, i)
		plot_everything(ma_entropy, ma_conservation, prefix + '_' + str(i) + '_movingAvg', gRNA_list = gRNA_freq)
	
def main():
	args = get_args()
	start_time = time.time()
	logging.basicConfig(level=args.loglevel, format=FORMAT)
	logger.info("cas13gRNAtor starting!!")
	temp_files_to_remove = []
	crRNA_handle, gRNA_class_list, temp_files_to_append, max_guide_score = get_gRNAs(args)
	temp_files_to_remove += temp_files_to_append

	offtarget_summary, temp_files_to_append = bowtie_main(args, crRNA_handle, mode = 'offtarget')
	gRNA_class_list = update_class(gRNA_class_list, offtarget_summary = offtarget_summary)
	temp_files_to_remove += temp_files_to_append

	gRNA_class_list, conservation_summary, temp_files_to_append, scoring_method, consensus_length = get_scores(args, gRNA_class_list, crRNA_handle = crRNA_handle, max_guide_score =  max_guide_score)
	temp_files_to_remove += temp_files_to_append

	write_supplementary_data(args, gRNA_class_list)
	gRNA_freq = {k:0 for k in range(consensus_length)}
	if scoring_method[1] and not scoring_method[0]:
		gRNA_class_list.sort(key=lambda x: x.c_score_bowtie, reverse=True)
	gRNA_freq = write_all(args, gRNA_class_list, gRNA_freq = gRNA_freq, scoring_method = scoring_method, plot = args.plot)
	if scoring_method[0]:
		write_best(args, gRNA_class_list, plot = args.plot)
		main_plot(args.prefix, gRNA_freq)

	if args.temp:
		for fname in temp_files_to_remove:
			os.remove(fname)      
		logger.info("Temporary files removed!")
		
	end_time = time.time()
	total_time = round(end_time - start_time,3)
	logger.info(f"cas13gRNAtor has Ended in {total_time} seconds!")


if __name__ == '__main__':
	main()

