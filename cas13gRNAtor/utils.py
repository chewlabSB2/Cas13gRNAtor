#!/usr/bin/env python3

import shlex
import os
import re
import sys
import subprocess
from textwrap import dedent
import argparse
import logging
import logging.handlers
import numpy as np
from Bio import AlignIO, SeqIO
logging.getLogger('matplotlib').setLevel(logging.WARNING)
import matplotlib.pyplot as plt
from tqdm import tqdm
from cas13gRNAtor.simpleBWA import *

PYTHON_FILE = os.path.dirname(os.path.abspath(__file__))

LENGTH_OF_GRNA = 23
RANGE_START, RANGE_END = 2, 7 #Sensitive region
LENGTH_OF_INTOLERANT = abs(RANGE_END - RANGE_START)

FORMAT = '%(asctime)s - %(module)-16s - %(levelname)s - %(message)s'

ALPHABET = ['A', 'C', 'G', 'T', 'N']
ALPHABET_NORMAL = ['A', 'C', 'G', 'T']

HEADER_MAIN_1 = [
	'GuideName',
	'GuideSeq',
	'MatchPos',
]

HEADER_MAIN_2 = ['GuideScores',
	'Rank',
	'standardizedGuideScores',
	'quartiles',
]

HEADER_MSA = [
	'Found',
	'Entropy Mean score',
	'Entropy SD score',
	'Conservation Mean Score',
	'Conservation SD Score',
	'Intolerant region Entropy Mean',
	'Intolerant region Entropy SD',
	'Intolerant region Conservation Mean',
	'Intolerant region Conservation SD',
]

HEADER_BEST_gRNAs = [
	'GuideName',
	'GuideSeq',
	'gRNA_Conservation Mean',
	'gRNA_Intolerant_region_Conservation Mean',
	'GuideScore',
	'Weighted_Score'
]

HEADER_WEIGHTED_SCORE = 'Weighted Score'
HEADER_BOWTIE_OFFTARGETS = 'offtargets'
HEADER_BOWTIE_CONSERVATION = 'Conservation %'

ASSERT_MESSAGE_ALIGNMENT = "Either a Multiple Sequence Alignment (MSA) Fasta File or Non-MSA Fasta File required"
ASSERT_MESSAGE_SCORE = "Guide Scores from RfxCas13d_GuideScoring.R or Reference File Needed"


class AlignmentError(Exception):
	pass

class Cas13gRNA():
	sensitive_pos = (RANGE_START, RANGE_END) #16-21 => Reverse Pos

	def __init__(self, crRNAid, seq):
		self.id = crRNAid
		self.seq = seq
		self.g_score = None
		self._pos = (0,0)
		self.rank = None 
		self.SGS = None
		self.quartile = None

		self.mode = None
		self.checked = False
		self.found = False

		self.c_mean = 0
		self.e_mean = 0
		self.c_mean_sensitive = 0
		self.e_mean_sensitive = 0

		self.c_SD = 0
		self.e_SD = 0
		self.c_SD_sensitive = 0
		self.e_SD_sensitive = 0

		self.c_score_list = []
		self.e_score_list = []
		self.mismatch = 0
		self.weighted_score = 0

		self.c_score_bowtie = None
		self.offtargets = []
		self.gRNA_len = LENGTH_OF_GRNA

	def scoring(self, cap_guidescore, get_weighted_score = True):
		self.c_mean = round(np.mean(self.c_score_list), 5)
		sample_variance = np.sum(((self.c_score_list - self.c_mean)**2) / (LENGTH_OF_GRNA - 1))
		self.c_SD = round(np.sqrt(sample_variance),5)
		
		self.e_mean = round(np.mean(self.e_score_list), 5)
		sample_variance = np.sum(((self.e_score_list - self.e_mean)**2) / (LENGTH_OF_GRNA - 1))
		self.e_SD = round(np.sqrt(sample_variance),5)
		
		self.c_mean_sensitive = round(np.mean(self.c_score_list[RANGE_START:RANGE_END]), 5)
		sample_variance = np.sum(((self.c_score_list[RANGE_START:RANGE_END] - self.c_mean_sensitive)**2) / (LENGTH_OF_INTOLERANT - 1))
		self.c_SD_sensitive = round(np.sqrt(sample_variance),5)
		
		self.e_mean_sensitive = round(np.mean(self.e_score_list[RANGE_START:RANGE_END]), 5)
		sample_variance = np.sum(((self.e_score_list[RANGE_START:RANGE_END] - self.e_mean_sensitive)**2) / (LENGTH_OF_INTOLERANT - 1))
		self.e_SD_sensitive = round(np.sqrt(sample_variance),5)
		
		self.checked = True
		if get_weighted_score:
			self.weighted_score = (0.5*(self.g_score/cap_guidescore) + 0.35*self.c_mean_sensitive  + 0.15*self.c_mean)
			self.weighted_score = round(self.weighted_score, 5)
		

	@property
	def pos(self):
		return self._pos

	@pos.setter
	def pos(self, start):
		end = start + LENGTH_OF_GRNA
		self._pos = (start, end) #Relative to Strand

def reverseC(seq):
	RT = {'A':'T','C':'G','G':'C','T':'A', 'N':'N'}
	reverseComplement = ''
	for i in seq:
		nt = RT.get(i)
		reverseComplement += nt
	return reverseComplement[::-1]

## Conservation And Entropy Score Calculation

def read_alignment(MSA, handle):
	## Return MSA if alignment works else return Null
	if MSA:
		try:
			return AlignIO.read(handle, 'fasta'), True
		except:
			try:
				test = SeqIO.parse(handle, 'fasta')
				return None , False
			except Exception as error:
				raise AlignmentError("\nERROR: Problem reading in {}: {}".format(handle, str(error)))
	else:
		try:
			test = AlignIO.read(handle, 'fasta')
			return None, True
		except:
			try:
				test = SeqIO.parse(handle, 'fasta')
				return None, False
			except Exception as error:
				raise AlignmentError("\nERROR: Problem reading in {}: {}".format(handle, str(error)))

def prettify_gaps(MSA):
	## Converts all bases to uppercase
	## Replace all gaps by 'N' in all sequences in the alignment post-alignment	
	from Bio.Seq import Seq
	for seq_ in MSA:
		_seq = str(seq_.seq)
		_seq = _seq.upper()
		_seq = _seq.replace('-', 'N')
		seq_.seq = Seq(_seq)

def calc_frequencies(MSA):
	## Returns a matrix of nucleotide frequecy for each position for entropy and conservation score
	aln_array = np.array(MSA) 
	mat = np.zeros((len(ALPHABET), aln_array.shape[1]))
	for ni, nuc in enumerate(ALPHABET):
		mat[ni] = np.mean(aln_array == nuc, axis=0)
	mat /= np.sum(mat, axis=0) 
	return mat

def calc_entropy(mat):
	## Shannon's Entropy for each position in MSA
	return [round(x, 5) for x in (-1) * np.sum(mat * np.log2(mat + 1e-15), axis=0)]

def consensus_seq(mat):
	## Get Consensus Sequence from matrix for wessels gRNA scoring
	consensus = []
	## max value from each array
	consensus_value = np.max(mat, axis=0) 
	## position/index of max value in array
	a = np.argmax(mat, axis=0) 
	
	for i in a:
		consensus.append(ALPHABET[i])

	return ''.join(consensus), consensus_value

def convert_to_fasta(consensus_sequence, prefix):
	from Bio.Seq import Seq
	from Bio.SeqRecord import SeqRecord

	sequence_to_save = SeqRecord(seq = Seq(consensus_sequence), id = 'Consensus', description = 'Consensus_Sequence')
	handle = open(prefix + "_consensus_sequence.fasta","w")
	SeqIO.write(sequence_to_save,handle,"fasta")
	handle.close()

def write_score(entropy_score_table, conservation_score_table, consensus_sequence, prefix):
	csv_file = prefix + '_position_score.csv'
	score_csv = open(csv_file, 'w')
	score_csv.write('Position\tEntropy\tConservation\tBase Pair\n')
	for i in range(len(entropy_score_table)):
		score_csv.write(str(i+1) + '\t' + str(entropy_score_table[i]) + '\t' + str(conservation_score_table[i]) +  '\t' + consensus_sequence[i] + '\n')
	score_csv.close()

def main_calc(MSA, prefix):
	prettify_gaps(MSA)
	matrix = calc_frequencies(MSA)
	entropy_list = calc_entropy(matrix)
	consensus, conservation_list = consensus_seq(matrix)
	conservation_list = [round(i, 5) for i in conservation_list]
	write_score(entropy_list, conservation_list, consensus, prefix)
	convert_to_fasta(consensus, prefix)
	return consensus, conservation_list, entropy_list 

## Bowtie Offtarget and Conservation Scoring
def is_tool(name):
	## Check whether `name` is on PATH and marked as executable.
	from shutil import which
	return which(name) is not None

def check_bowtie(args):
	if not is_tool('bowtie') or not is_tool('bowtie-build'):
		if not is_tool(PYTHON_FILE + '/misc/bowtie') or not is_tool(PYTHON_FILE + '/misc/bowtie-build'): 
			raise AlignmentError("Bowtie is not installed")

	if args.bowtie or args.bowtie_index or args.gRNA_mode[0].lower() == 'low': #Edit
		if not is_tool('bowtie') or not is_tool('bowtie-build'):
			if not is_tool(PYTHON_FILE + '/misc/bowtie') or not is_tool(PYTHON_FILE + '/misc/bowtie-build'): 
				raise AlignmentError("Bowtie is not installed")
			else:
				bowtie_path = PYTHON_FILE + '/misc/bowtie'
				bowtie_build_path = PYTHON_FILE + '/misc/bowtie-build'
		else:           
			bowtie_path = 'bowtie'
			bowtie_build_path = 'bowtie-build'

	return bowtie_path, bowtie_build_path

		
def bowtie_build_cmd(bowtie_reference, prefix, bowtie_path = 'bowtie-build', thread = 4):
	cmd = '%s --threads %d %s %s_alignment.index 1>>%s_output_build.txt 2>>%s_error_build.txt'%(bowtie_path, thread, bowtie_reference, prefix, prefix, prefix)
	return cmd, prefix + '_alignment.index'

def bowtie_cmd(bow, reference_fasta, prefix, mismatch, bowtie_path = 'bowtie', thread = 4):
	cmd = '%s -a -S %s -f %s %s.sam -v %d --threads %d'%(bowtie_path, bow, reference_fasta, prefix, mismatch, thread)
	return cmd

def bowtie_to_summary(file, fasta_dict = {}, mode = "offtarget"):
	file = open(file, '+r')
	filtered_sgrna = dict()
	seq_count = 0
	for line in file:
		line = line.strip('\n')
		if line.startswith('@SQ'):
			match = re.findall(r"\tSN:(\S+)", line)
			name = match[0] if len(match) else "unknown_reference"
			seq_count += 1
			continue
		if line.startswith('@'): continue
		cols = line.split("\t")
		sgrna_name = cols[0]
		mapping_status = int(cols[1])
		genome_name = cols[2]
		if mapping_status == 16:
			if mode == "conservation":
				fasta_dict[sgrna_name] += 1
		
			if mode == "offtarget":
				if not sgrna_name in list(fasta_dict.keys()):
					fasta_dict[sgrna_name] = [genome_name]
				else:
					fasta_dict[sgrna_name].append(genome_name)

			
				## Consider Mutated BP in sensitive region?

	if mode == "conservation":
		fasta_dict = {k:round(v/seq_count * 100, 5) for k,v in fasta_dict.items()}

	return fasta_dict

def update_class(gRNA_class_list, conservation_summary = {}, offtarget_summary = {}):
	updated_gRNA = []
	for gRNA in gRNA_class_list:
		if conservation_summary:
			conservation_perc = conservation_summary[gRNA.id]
			gRNA.c_score_bowtie = conservation_perc
		if offtarget_summary:
			if gRNA.id in list(offtarget_summary.keys()):
				gRNA.offtargets = offtarget_summary[gRNA.id]
		updated_gRNA.append(gRNA)

	return updated_gRNA

## Rscript
def check_Rscripts():
	if not is_tool('Rscript'):
		raise AlignmentError("Rscript is not installed")
	if not os.path.exists(PYTHON_FILE + '/misc/RfxCas13d_GuideScoring.R'):
		raise AlignmentError("RfxCas13d_GuideScoring.R in scripts directory is not found")
	if not os.path.exists(PYTHON_FILE + '/misc/data/Cas13designGuidePredictorInput.csv'):
		raise AlignmentError('Files for gRNA scoring and prediction in data directory is missing. Make sure Cas13designGuidePredictorInput.csv is present')
	if not os.path.exists(PYTHON_FILE + '/misc/data/LocalNTdensityCorrelations.txt'):
		raise AlignmentError('Files for gRNA scoring and prediction in data directory is missing. Make sure LocalNTdensityCorrelations.txt is present')

def check_Rscripts_tools():
	RNAfold = PYTHON_FILE + '/misc/RNAfold'
	RNAplfold = PYTHON_FILE + '/misc/RNAplfold'
	RNAhybrid = PYTHON_FILE + '/misc/RNAhybrid'
	if not is_tool(RNAhybrid):
		raise AlignmentError('RNAhybrid is not installed properly')
	if not is_tool(RNAfold):
		raise AlignmentError('RNAfold is not installed properly')
	if not is_tool(RNAplfold):
		raise AlignmentError('RNAplfold is not installed properly')
	

def generate_gRNA_scores(fasta):
	cmd = 'Rscript ' + PYTHON_FILE + '/misc/RfxCas13d_GuideScoring.R  ' + fasta + ' ' + PYTHON_FILE + '/misc/data/Cas13designGuidePredictorInput.csv true'
	return cmd

def get_csv_file(directory, suffix):
	files = os.listdir(directory)
	temp_list = list()
	for f in files:
		part = f.split('_')
		if part[-1] == suffix:
			temp_list.append(f)
	latest_file = max(temp_list, key=os.path.getmtime)
	return latest_file

def filter_gRNA(args, csv_file):
	#global max_guide_score
	removeUUUU = args.homopolymer
	max_guide_score = 0
	csv_file = open(csv_file, 'r')
	gRNA_classes = [] 
	crRNA_handle = args.prefix + '_4thQuartile_gRNA.fasta'
	filtered_crRNA_fasta = open(crRNA_handle, 'w+')

	next(csv_file)
	for line in csv_file:
		line = line.strip('\n')
		columns = line.split(',')
		sequence = columns[1].upper()
		if 'N' in sequence:
			continue
		if removeUUUU and 'TTTT' in sequence:
			continue
		if int(columns[-1]) >= 3:
			guidescore = float(columns[3])
			if max_guide_score < guidescore:
				max_guide_score = guidescore

			rank = float(columns[4])
			SGS = float(columns[5])
			quartile = int(columns[6])
			gRNA = Cas13gRNA(columns[0], sequence)
			gRNA.g_score = guidescore
			gRNA.rank = rank
			gRNA.SGS = SGS
			gRNA.quartile = quartile
			gRNA.pos = int(columns[2]) - LENGTH_OF_GRNA - 1#0-based
			filtered_crRNA_fasta.write('>' + columns[0] + '\n' + sequence + '\n')
			gRNA_classes.append(gRNA)

	if not args.MSA: filtered_crRNA_fasta.close()
	return crRNA_handle, gRNA_classes, max_guide_score

def check_format_gRNAscore():
	pass


## Plotting

def get_moving_average(score, window_size = 10):
	#score = list(score.values())
	position = [x for x in range(len(score))]
	i = 0

	moving_averages = []
	if window_size == 0: window_size = 1
	while i < len(score) - window_size + 1:
		this_window = score[i : i + window_size]
		window_average = sum(this_window) / window_size
		moving_averages.append(window_average)
		i += 1
	
	#moving_average = dict(zip(position, moving_averages))
	return moving_averages

def plot_everything(entropy, conservation, prefix, gRNA_list = []):
	plt.switch_backend('agg')
	ax1 = plt.subplot(311) if gRNA_list else plt.subplot(211)
	plt.plot([i for i in range(len(entropy))], entropy, linewidth=0.5)
	plt.title('Entropy/Conservation/gRNA scoring plot')
	plt.ylabel('Entropy Score', fontsize=6) 
	plt.rc('ytick',labelsize=4)
	plt.setp(ax1.get_xticklabels(), fontsize=6)

	ax2 = plt.subplot(312, sharex=ax1) if gRNA_list else plt.subplot(212, sharex=ax1)
	plt.plot([i for i in range(len(conservation))], conservation, linewidth=0.5)
	plt.ylabel('Conservation Score', fontsize=6)
	plt.rc('ytick',labelsize=4) 
	plt.setp(ax2.get_xticklabels(), fontsize=6)

	if gRNA_list:
		ax3 = plt.subplot(313, sharex=ax1)
		plt.bar(list(gRNA_list.keys()), gRNA_list.values())
		plt.ylabel('gRNA score frequency', fontsize=6)
		plt.rc('ytick',labelsize=4)
		plt.setp(ax3.get_xticklabels(), fontsize=6)
		plt.xlabel('Nucleotide position')
	
	figure = plt.gcf() # get current figure
	figure.set_size_inches(12, 4)
	plt.savefig(prefix + '_combined_genome_plot.png', dpi = 360)

def break_apart(word):
	return [char for char in word]

def plot_highly_conserved_gRNA(crRNA_id, conservation_score, entropy_score, sequence, pos, prefix):
	os.chdir(prefix + '_best_gRNAs')
	curr_range = [i + 1 for i in range(pos[0], pos[1])]

	tick_spacing = 1    
	import matplotlib.ticker as ticker
	
	plt.switch_backend('agg')
	#plt.title('Entropy/Conservation scoring plot for {}'.format(crRNA_id))
	
	ax1 = plt.subplot(211)
	ax1.xaxis.set_major_locator(ticker.MultipleLocator(tick_spacing))
	plt.plot(curr_range, entropy_score)
	#plt.xticks(np.arange(len(v[3])), v[3])
	plt.title(f'Entropy & Conservation Score Plot for {sequence}')
	plt.axvspan(RANGE_START + pos[0] + 1, RANGE_END + pos[0] + 1, color='r', alpha=0.5)
	plt.ylabel('Entropy Score', fontsize=8)
	plt.rc('ytick',labelsize=6)
	ax1.set_xticks(curr_range)
	ax1.set_xticklabels(break_apart(reverseC(sequence)))
	plt.setp(ax1.get_xticklabels(), fontsize=6)

	#ax2 = plt.subplot(212, sharex=ax1)
	ax2 = plt.subplot(212)
	plt.plot(curr_range, conservation_score)
	#if plot_y: plt.axhline(y=conservation_score, linewidth=0.75, color='k', alpha=0.75)
	plt.axvspan(RANGE_START + pos[0] + 1, RANGE_END + pos[0] + 1, color='r', alpha=0.5)
	plt.ylabel('Conservation Score', fontsize=8)
	
	plt.rc('ytick',labelsize=6)
	plt.setp(ax2.get_xticklabels(), fontsize=6)
	plt.xlabel('Nucleotide position')

	figure = plt.gcf()  # get current figure
	figure.set_size_inches(12, 4)
	crRNA_id = crRNA_id.split(':')[0]
	
	plt.savefig(prefix + '_' + str(crRNA_id) + '_plot.png', dpi=360)
	os.chdir("..")

## MultiProcessing 

def partition(lst, n):
	division = len(lst) / n
	return [lst[round(division * i):round(division * (i + 1))] for i in range(n)]

def ranges(N, nb):
	temp_list = [(r.start, r.stop) for r in partition(range(len(N)), nb)]
	for i in temp_list: yield (i[0],i[1])

def process_gRNA(consensus, gRNA_list, queue, conservation_list, entropy_list, max_guide_score, mismatch = 0, get_weighted_score = True):
	bwa = simpleBWA(consensus)
	modified_gRNA = []
	for gRNA in gRNA_list:
		sequence = gRNA.seq 
		matches = bwa.find_match(reverseC(sequence), mismatch)
		if len(matches) >= 1:
			gRNA.mode = "MSA"
			gRNA.pos = matches[0] 
			#start, end = gRNA.pos, gRNA.pos + LENGTH_OF_GRNA
			start, end = gRNA.pos[0], gRNA.pos[1]
			gRNA.e_score_list = entropy_list[start:end]
			gRNA.c_score_list = conservation_list[start:end]
			gRNA.scoring(max_guide_score, get_weighted_score = get_weighted_score)
			gRNA.found = True
			modified_gRNA.append(gRNA)
		else:
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

def write_all(args, gRNA_class_list, gRNA_freq = {}, scoring_method = (True, False), plot = True, weighted = True):
	with open(args.prefix + "_combined_results.csv", 'w') as write_all:
		write_all.write(",".join(HEADER_MAIN_1))
		if weighted: write_all.write("," + ",".join(HEADER_MAIN_2))
		if scoring_method[0]: ## MSA
			write_all.write("," + ",".join(HEADER_MSA))
			if weighted: write_all.write("," + HEADER_WEIGHTED_SCORE)
		if scoring_method[1]: ## Conservation
			write_all.write("," + HEADER_BOWTIE_CONSERVATION)
		if args.offtarget: ## Offtarget
			write_all.write("," + HEADER_BOWTIE_OFFTARGETS)

		write_all.write('\n')
		for g in gRNA_class_list:
			write_all.write(f"{g.id},{g.seq},{g.pos[0]}")
			if weighted: write_all.write(f",{g.g_score},{g.rank},{g.SGS},{g.quartile}")
			if scoring_method[0]: ## MSA
				if g.found:
					for i in range(g.pos[0], g.pos[1]):
						gRNA_freq[i] += 1
				write_all.write(f",{g.found},{g.e_mean},{g.e_SD},{g.c_mean},{g.c_SD},{g.e_mean_sensitive},{g.e_SD_sensitive},{g.c_mean_sensitive},{g.c_SD_sensitive}")
				if weighted: write_all.write(f",{g.weighted_score}")
			if scoring_method[1]:
				write_all.write(f",{g.c_score_bowtie}")
			if args.offtarget:
				if g.offtargets:
					for genome in g.offtargets:
						write_all.write(f",{genome}")

			write_all.write('\n')

	return gRNA_freq

def write_best(args, gRNA_class_list, top = 100, plot = True):
	if plot:
		if not os.path.isdir(args.prefix + '_best_gRNAs'): os.mkdir(args.prefix + '_best_gRNAs')  
	
	write_bowtie_score = False
	with open(args.prefix + "_bestgRNAs.csv", 'w') as write_best:
		write_best.write(",".join(HEADER_BEST_gRNAs))
		if gRNA_class_list[0].c_score_bowtie > 0: 
			write_best.write("," + HEADER_BOWTIE_CONSERVATION)
			write_bowtie_score = True
		write_best.write("\n")
		count = 0
		for g in gRNA_class_list:
			if count > top: break
			if g.offtargets: continue
			write_best.write(f"{g.id},{g.seq},{g.c_mean},{g.c_mean_sensitive},{g.g_score},{g.weighted_score}")
			if write_bowtie_score: 
				write_best.write(f",{g.c_score_bowtie}")
			write_best.write("\n")
			count += 1
			if plot:
				plot_highly_conserved_gRNA(g.id, g.c_score_list, g.e_score_list, g.seq, g.pos, args.prefix)

def write_supplementary_data(args, gRNA_class_list):
	with open(args.prefix + "_conservation.csv", 'w') as write_conservation:
		with open(args.prefix + "_entropy.csv", 'w') as write_entropy:
			for g in gRNA_class_list:
				entropy_list =  [str(round(x, 5)) for x in g.e_score_list]
				conservation_list =  [str(round(x, 5)) for x in g.c_score_list]
				write_entropy.write(f"{g.id},{g.seq},{g.pos[0]}," + ",".join(entropy_list) + "\n")
				write_conservation.write(f"{g.id},{g.seq},{g.pos[0]}," + ",".join(conservation_list) + "\n")

def multi_search(args, consensus, gRNA_list, conservation_list, entropy_list, max_guide_score, get_weighted_score = True):
	from multiprocessing import Process, Queue, Manager, RawArray
	import ctypes

	number_processors = args.CPU
	current_processes = []
	q = Queue()
	manager = Manager()
	edit_gRNAs = manager.list()
	conservation_list = RawArray(ctypes.c_float, conservation_list)
	entropy_list = RawArray(ctypes.c_float, entropy_list)

	for i in ranges(gRNA_list, number_processors):
		current_gRNAs = gRNA_list[i[0]:i[1]]
		pp = Process(target = process_gRNA, args=(consensus, current_gRNAs, q, conservation_list, entropy_list, max_guide_score, args.mismatch, get_weighted_score))
		current_processes.append(pp)

	for pp in current_processes:
		pp.start()

	sink = Process(target=queue_sam, args=(q, edit_gRNAs))
	sink.start()

	for pp in current_processes:
		pp.join()

	q.put('Done')
	sink.join()

	edit_gRNAs = list(edit_gRNAs)
	if edit_gRNAs[0].weighted_score: edit_gRNAs.sort(key=lambda x: x.weighted_score, reverse=True)
	return edit_gRNAs

## Other Functions
def print_error(message, **kwargs):
	"""
	Formats *message* with *kwargs* using :meth:`str.format` and
	:func:`textwrap.dedent` and uses it to print an error message to
	``sys.stderr``.
	"""
	print("\nERROR: " + dedent(message.format(**kwargs)).lstrip("\n")+"\n", file = sys.stderr) 

def run_shell_command(cmd, raise_errors = False, extra_env = None):
	"""
	Run the given command string via Bash with error checking.
	Returns True if the command exits normally.  Returns False if the command
	exits with failure and "raise_errors" is False (the default).  When
	"raise_errors" is True, exceptions are rethrown.

	If an *extra_env* mapping is passed, the provided keys and values are
	overlayed onto the default subprocess environment.
	"""
	env = os.environ.copy()

	if extra_env:
		env.update(extra_env)

	shargs = ['-c', "set -euo pipefail; " + cmd]

	if os.name == 'posix':
		shellexec = ['/bin/bash']
	else:
		# We try best effort on other systems. For now that means nt/java.
		shellexec = ['env', 'bash']

	try:
		# Use check_call() instead of run() since the latter was added only in Python 3.5.
		subprocess.check_output(
			shellexec + shargs,
			shell = False,
			stderr = subprocess.STDOUT,
			env = env)

	except subprocess.CalledProcessError as error:
		print_error(
			"{out}\nshell exited {rc} when running: {cmd}{extra}",
			out = error.output,
			rc  = error.returncode,
			cmd = cmd,
			extra = "\nAre you sure this program is installed?" if error.returncode==127 else "",
		)
		if raise_errors:
			raise
		else:
			return False

	except FileNotFoundError as error:
		print_error(
			"""
			Unable to run shell commands using {shell}!
			Module requires {shell} to be installed.
			""",
			shell = ' and '.join(shellexec)
		)
		if raise_errors:
			raise
		else:
			return False

	else:
		return True  


