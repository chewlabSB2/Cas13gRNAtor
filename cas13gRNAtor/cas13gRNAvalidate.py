#!/usr/bin/env python3

import os
import argparse
import time
import logging
import logging.handlers
from tqdm import tqdm
from cas13gRNAtor.utils import *
from cas13gRNAtor.cas13gRNAtor import *
from cas13gRNAtor.simpleSearch import *

logger = logging.getLogger(__name__)

HEADER_REFERENCE = [
	'Found',
	'Sequence_Found',
	'Mismatch',
	'Cigar',
	'Intolerant',
	] 

DESCRIPTION = '''
Cas13gRNAvalidate Validates Previous gRNAs used against an MSA or non-MSA
Scores Shannon Entropy and Conservation Scores for each gRNA
Input: Aligned Sequences.
Output: Conservation and Entropy score in csv as supplementary data'
'''

##include validate score 

def get_args():
	parser = argparse.ArgumentParser(
		prog='cas13gRNAvalidate',
		description=DESCRIPTION
	)
	parser.add_argument('-i', '--input', required=True, dest='gRNA',
						help="Takes in a list of gRNAs in a fasta file or a text file, list seperated by a newline")
	parser.add_argument('-m','--MSA', dest = "MSA", metavar="fasta", default = None, 
						help="Sequences to align")
	parser.add_argument('-bc', '--bowtie-conservation', dest = "bowtie", default = None, 
						help="Path to fasta file for Bowtie conservation calculation")
	parser.add_argument('-bo', '--bowtie-offtarget', dest = "offtarget", 
						help="Path to fasta file for Bowtie conservation calculation")
	parser.add_argument('--mismatch', default=2, type=int,  choices=range(0,10),
						help='Number of mismatches allowed for Bowtie (Default: %(default)s)')
	parser.add_argument('-t', '--threads', default=4, type=int, dest = 'CPU', 
						help="Number of Core to use (Default: %(default)s)")
	parser.add_argument('-p', '--prefix', default = 'RfxCas13d_gRNA' , 
						help='file name prefix to your file names (Default: %(default)s)')
	parser.add_argument('-r', '--reference', dest = 'reference', metavar="fasta", 
						help="Reference for Guide Score Calculation")
	parser.add_argument('--guide-score', dest = 'scorefile', metavar="csv", 
						help="Guide Score From RfxCas13d_GuideScoring.R")
	parser.add_argument('-d', '--debug',
						help='Print lots of debugging statements',
						action="store_const",dest="loglevel",const=logging.DEBUG,
						default=logging.INFO)
	parser.add_argument('--conservation-index', dest="conservation_index", default=None, 
						help='Input an index for conservation calculation if you have!')
	parser.add_argument('--offtarget-index', dest="offtarget_index", default=None, 
						help='Input an index for offtargets if you have!')
	parser.add_argument('--keep-tmp', action="store_false", dest='temp', help='Keep Temporary files (Default: Does not Keep)')
	parser.add_argument('--no-plot', action="store_false",dest="plot", help='Does not Plot graph')

	args = parser.parse_args()
	if not args.MSA:
		assert args.bowtie or args.conservation_index, ASSERT_MESSAGE_ALIGNMENT
	elif not args.bowtie or not args.conservation_index:
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
		curr_id = 'gRNA'
		for line in f:
			if line.startswith('>'): 
				curr_id = line.strip('\n').replace('>','')
				continue
			seq = line.strip('\n')
			seq = seq.upper()
			if len(seq) < LENGTH_OF_GRNA: 
				logger.info(f"{seq} length < 23, Will Ignore!")
				continue

			list_of_seq = []
			if len(seq) > LENGTH_OF_GRNA:
				for s in range(len(seq) - LENGTH_OF_GRNA):
					temp_seq = seq[s:s+LENGTH_OF_GRNA]
					list_of_seq.append(temp_seq)
			else:
				list_of_seq.append(seq)

			for c in break_apart(seq):
				if not c in ALPHABET_NORMAL:
					logger.info(f"{seq} contains invalid Nucleotide Alphabet/s, Will Ignore!")
					continue
			else:
				crRNA_count += 1
				temp_count = 0
				for s in list_of_seq:
					temp_count += 1
					gRNA = Cas13gRNA(f"{curr_id}_{crRNA_count}_{temp_count}", s)
					gRNA_class_list.append(gRNA)
					crRNA_file.write(f">{curr_id}_{crRNA_count}_{temp_count}\n{s}\n")

	crRNA_file.close()
	return crRNA_handle, gRNA_class_list

def open_scorefile(scorefile):
	reference = ''
	with open(scorefile, 'r') as f:
		next(f)
		for i, line in enumerate(f): 
			line = line.strip('\n')
			columns = line.split(',')
			sequence = columns[1].upper()
			if i == 0:
				reference += sequence
			else:
				reference += sequence[-1]
	return reference			

def _get_gRNAs(args, random):
	temp_files_to_remove = []
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
	
	Rfiles = [f'{reference.id}_CasRxguides.fa', f'{reference.id}_CasRxguides.csv'] #, 'Consensus_Sequence_CasRxguides.pdf'
	temp_files_to_remove += Rfiles

	if args.temp:
		for fname in temp_files_to_remove:
			os.remove(fname)
		logger.info("Temporary files removed!")

def search_gRNA(wavelet, gRNA_list, queue):
	consensus = wavelet.ReconstructSequence()
	modified_gRNA = []
	for gRNA in gRNA_list:
		sequence = gRNA.seq 
		matches = bitapSearch(consensus, reverseC(sequence), best = True)
		if matches:
			best = matches.pop(0)
			gRNA.pos_reference = best.pos[0]
			gRNA.cigar_reference = best.cigar 
			gRNA.Intolerant_reference = best.intolerant
			gRNA.mismatch_reference = best.editDistance
			gRNA.gRNA_query_reference = reverseC(best.query)
			gRNA.found_reference = True
		else:
			gRNA.found_reference = False

		modified_gRNA.append(gRNA)
	queue.put(modified_gRNA)

def queue_sam(q, edit_gRNAs):
	while True:
		item = q.get()
		if item == 'Done':
			return
		else:
			for i in item:
				edit_gRNAs.append(i)

def write_gRNAs(args, edit_gRNAs):
	with open(args.prefix + "_validate_guideScores.csv", 'w') as write_all:
		write_all.write(",".join(HEADER_MAIN_1))
		write_all.write("," + ",".join(HEADER_MAIN_2))
		write_all.write("," + ",".join(HEADER_REFERENCE))
		write_all.write('\n')

		for g in edit_gRNAs:
			write_all.write(f"{g.id},{g.seq},{g.pos_reference[0]}")
			write_all.write(f",{g.g_score},{g.rank},{g.SGS},{g.quartile}")
			write_all.write(f",{g.found_reference},{g.gRNA_query_reference},{g.mismatch_reference},{g.cigar_reference},{g.Intolerant_reference}\n")


def update_score(csv_file, edit_gRNAs):

	modified_gRNA = []
	temp_gRNA = [(g.pos_reference[0], g) for g in edit_gRNAs]
	temp_gRNA_pos = [p for p, g in temp_gRNA]
	with open(csv_file, 'r') as f:
		next(f)
		for line in f:
			line = line.strip('\n')
			columns = line.split(',')
			pos = int(columns[2]) - LENGTH_OF_GRNA - 1
			if pos not in temp_gRNA_pos: continue

			guidescore = float(columns[3])
			rank = float(columns[4])
			SGS = float(columns[5])
			quartile = int(columns[6])
			gRNAs = [g for p, g in temp_gRNA if p == pos]

			for gRNA in gRNAs:
				gRNA.g_score = guidescore
				gRNA.rank = rank
				gRNA.SGS = SGS
				gRNA.quartile = quartile
				modified_gRNA.append(gRNA)

	return modified_gRNA


def multi_gRNA_score(args, gRNA_class_list, reference = None):
	from multiprocessing import Process, Queue, Manager
	
	q = Queue()
	temp_files_to_remove = []
	current_processes = []
	number_processors = args.CPU

	if args.reference:
		number_processors -= 1
		p = Process(target=_get_gRNAs, args = (args, None))
		current_processes.append(p)
		wavelet = fasta_2_Wavelet(args.reference)

	manager = Manager()
	edit_gRNAs = manager.list()
	
	del reference

	for i in ranges(gRNA_class_list, number_processors):
		current_gRNAs = gRNA_class_list[i[0]:i[1]]
		if not current_gRNAs: continue
		pp = Process(target = search_gRNA, args=(wavelet, current_gRNAs, q))
		current_processes.append(pp)

	for pp in current_processes:
		pp.start()

	sink = Process(target=queue_sam, args=(q, edit_gRNAs))
	sink.start()

	for pp in current_processes:
		pp.join()

	q.put('Done')
	sink.join()

	cwd = os.path.abspath(os.getcwd())
	edit_gRNAs = list(edit_gRNAs)
	csv_file = args.scorefile if args.scorefile else get_csv_file(cwd, 'CasRxguides.csv') 
	edit_gRNAs = update_score(csv_file, edit_gRNAs)

	write_gRNAs(args, edit_gRNAs)
	return edit_gRNAs

def main():
	args = get_args()
	start_time = time.time()
	logging.basicConfig(level=args.loglevel, format=FORMAT)
	logger.info("cas13gRNAvalidate starting!!")
	temp_files_to_remove = []

	crRNA_handle, gRNA_class_list = open_gRNA_file(args) 
	
	reference = None
	if args.scorefile:
		reference = open_scorefile(args.scorefile)

	if args.scorefile or args.reference:
		gRNA_class_list = multi_gRNA_score(args, gRNA_class_list, reference)

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
