#!/usr/bin/env python3

import os
import logging
import logging.handlers
import time
from functools import reduce
from matplotlib.lines import Line2D
from cas13gRNAtor.utils import *
from cas13gRNAtor.cas13gRNAtor import get_gRNAs

logger = logging.getLogger(__name__)
	
PROG = 'Cas13gRNApanvirus'

DESCRIPTION = f'''
{PROG} allows users to check if Pan-virus gRNAs exist
using either a Bowtie or MSA mode
Bowtie sacrifices on-target efficacy and potential regions of mutations
MSA mode sacrifices sensitivity but considers potential regions of mutations via entropy score 
'''

HEADER_PAN_BOWTIE = [
	"GuideName",
	"Sequence",
	"GuideScore",
	"Percentage",
	"Mismatch",
	]

def get_args():
	parser = argparse.ArgumentParser(
		prog = PROG,
		description = DESCRIPTION
	)

	parser.add_argument('-r', '--reference', dest = 'main', required = True, 
						help="Guide Score for Main Reference")
	parser.add_argument('-t', '--threads', default=4, type=int, dest = 'CPU', 
						help="Number of Core to use (Default: %(default)s)")
	parser.add_argument('-p', '--prefix', default = 'Cas13gRNAevaluate', 
						help='file name prefix to your file names (Default: %(default)s)')	
	parser.add_argument('-d', '--debug',
						help='Print lots of debugging statements',
						action="store_const",dest="loglevel",const=logging.DEBUG,
						default=logging.INFO)
	parser.add_argument('-m', '--mismatch', default=5, type=int,  choices=range(0,10),
						help='Number of mismatches allowed for Bowtie (Default: %(default)s)')
	parser.add_argument('--keep-tmp', action="store_false", dest='temp', 
						help='Keep Temporary files (Default: Does not Keep)')
	parser.add_argument('--no-plot', action="store_false",dest="plot", 
						help='Does not Plot graph')

	subparsers = parser.add_subparsers(help='Select mode', dest='mode')

	msa_parser = subparsers.add_parser('MSA', help='Heuristically Calculate Conservation Score and On-Target Efficacy Score')
	msa_parser.add_argument('-f', '--file',  dest='position', required=True, nargs = '+',
						metavar="txt", help="accepts score files from Cas13gRNAscore module")
	
	bowtie_parser = subparsers.add_parser('bowtie', help='Calculate Conservation Score based using bowtie (Sacrifices On-target Efficacy and disallow possible mutations)')
	bowtie_parser.add_argument('--best', dest = "best_count", type=int, default=10,
						help="Minimun Number of gRNAs to Select")
	bowtie_parser.add_argument('-bc', '--bowtie-conservation', dest = "bowtie", nargs = '+',
						help="Path to fasta file for Bowtie conservation calculation")
	#bowtie_parser.add_argument('--conservation-index', dest="conservation_index", default=None, nargs = '+',
	#					help='Input an index for conservation calculation if you have!')
	bowtie_parser.add_argument('--ignore-ontarget', action="store_false", dest='ignore_efficacy', 
						help='Keep Temporary files (Default: Does not Keep)')
	
	args = parser.parse_args() 
	#if args.conservation_index:
	#	args.temp = False

	return args

def str2bool(v):
	return v.lower() in ("true", "1")

class PostCas13gRNA_counter():
	bowtie_count = 0
	bowtie_all_new = []
	max_mismatch = 3

class PostCas13gRNA():
	#bowtie_count = 0
	#bowtie_all_new = []
	#max_mismatch = 3

	## Main Score File
	def __init__(self, GuideName, GuideSeq, MatchPos, 
				 Reference, GuideScore, ScoreRanking, standardizedGuideScores, quartiles, 
				 Sequence, Consensus_Pos, Mismatch, Cigar, Intolerant, 
				 e_sum, c_mean, c_SD, 
				 e_sum_sensitive, c_mean_sensitive, c_SD_sensitive, Weighted_Score, Weighted_Rank, 
				 Rank = float('inf'), Conservation_Bowtie = None, mainFile = True):
		
		self.id               = GuideName
		self.seq              = GuideSeq.upper()
		self.pos              = int(MatchPos)
		self.pos_consensus    = int(Consensus_Pos)
		self.g_score          = round(float(GuideScore), 5)  
		self.rankscore        = round(float(ScoreRanking), 5) 
		self.SGS              = round(float(standardizedGuideScores), 5)
		self.quartile         = int(quartiles)
		self.reference        = Reference
		self.rank             = Rank

		## Relative to Consensus Sequence
		self.mode             = 'MSA'
		self.found            = True
		self.consensus_query  = Sequence.upper()
		self.cigar            = Cigar
		self.Intolerant       = str2bool(Intolerant)
		self.mismatch         = int(Mismatch)

		self.c_mean     	  = float(c_mean)
		self.e_sum     	      = float(e_sum)
		self.c_mean_sensitive = float(c_mean_sensitive)
		self.e_sum_sensitive  = float(e_sum_sensitive)
		self.c_SD 			  = float(c_SD)
		self.c_SD_sensitive   = float(c_SD_sensitive)
		self.Weighted_rank    = float(Weighted_Rank)
		self.weighted_score   = float(Weighted_Score)

		self.c_score_bowtie   = float(Conservation_Bowtie) if Conservation_Bowtie else None
		self.gRNA_len         = LENGTH_OF_GRNA

		## For Reference Sequence Search
		self.mainFile         = mainFile
		self.pan_virus        = []
		self.ref_checked      = 0 
		self.ref_found        = 0
		self.final_rank       = -1

		## Bowtie Mode
		self.count            = [0 for i in range(PostCas13gRNA_counter.max_mismatch + 1)]
		self.gsXbw_rank       = 0
		self.guidescore_rank  = 0
		self.bowtie_rank      = 0
		self.bowtie_percent   = [0 for i in range(PostCas13gRNA_counter.max_mismatch + 1)]
		self.bowtie_temp_ct   = 0
		self.bowtie_list      = []
		self.bowtie_list_mm   = []

	def update_bowtie_percent(self):
		#count_list = {k:self.bowtie_list_mm.count(k) for k in range(PostCas13gRNA_counter.max_mismatch + 1)}
		current_count = 0
		for m in range(PostCas13gRNA_counter.max_mismatch + 1):
			current_count += self.count[m]
			self.bowtie_percent[m] = round(current_count/PostCas13gRNA_counter.bowtie_count * 100, 5)

	def update_bowtie_rank(self, rank):
		self.bowtie_rank = rank
		self.gsXbw_rank = self.bowtie_rank * self.guidescore_rank
		self.gsXbw_rank = round(math.log10(self.gsXbw_rank),5)

	@staticmethod
	def update_bowtie_count(value):
		PostCas13gRNA_counter.bowtie_count += value

	@staticmethod
	def update_bowtie_list(genome_list):
		PostCas13gRNA_counter.bowtie_all_new += genome_list

	def _get_final_rank(self, pan_virus_quantity):
		if self.ref_found < pan_virus_quantity or self.ref_checked < pan_virus_quantity:
			return

		rank = []
		for g in self.pan_virus:
			rank.append(g.RANK)

		self.final_rank = reduce((lambda x, y: x * y), rank) * self.rank  
		self.final_rank = round(math.log10(self.final_rank),5)

def _open_scorefile(scoreFile, mainFile = True, mismatch_score = 5):
	col_len = len(HEADER_MAIN_1) + len(HEADER_MAIN_2) + len(HEADER_MSA) + 2
	offtarget_mode = conservation_bowtie_mode = False
	rank_count = 0
	gRNA_list = []
	with open(scoreFile, 'r') as f:
		try:
			for line in f:
				if line.startswith('GuideName'):
					col = line.strip('\n').split(',')
					offtarget_mode = HEADER_BOWTIE_OFFTARGETS in col
					conservation_bowtie_mode = HEADER_BOWTIE_CONSERVATION in col
					continue

				line = line.strip('\n')
				col = line.split(',')
				## Ignore Offtargets
				if offtarget_mode:
					if conservation_bowtie_mode and len(col) > col_len + 1: continue
					if not conservation_bowtie_mode and len(col) > col_len: continue
				

				mismatch = int(col[10])
				if mismatch > mismatch_score:
					continue
				rank_count += 1
				conservation_bowtie = col[col_len] if conservation_bowtie_mode else None
				col = col[:col_len]
				gRNA = PostCas13gRNA(*col, Rank = rank_count, Conservation_Bowtie = conservation_bowtie, mainFile = mainFile)
				gRNA_list.append(gRNA)
		except:
			logger.info(f'Error in {scoreFile}. Ignoring File, Format Incorrect!')
			logger.info(f'Please Ensure {scoreFile} is derived from Cas13gRNAtor using Multiple Sequence Alignment')
			return [], None

	return gRNA_list, conservation_bowtie_mode

## MSA Mode

class Virus_Score():

	def __init__(self, name, entropy, conservation, sequence):
		self.name = name
		self.entropy = entropy
		self.conservation = conservation 
		self.sequence = Wavelet_Tree(sequence)

	def write_to_fasta(self):
		from Bio.Seq import Seq
		from Bio.SeqRecord import SeqRecord
		handle = f"{self.name}_consensus"
		sequence = self.sequence.ReconstructSequence()
		record = SeqRecord(Seq(sequence), id=self.name, name=handle)
		with open(f"{handle}.fasta", "w") as output_handle:
			SeqIO.write(record, output_handle, "fasta")
		return f"{handle}.fasta"

def _search_gRNA(pan_virus, mismatch, main_current_gRNAs, gRNA_pos, queue):
	consensus = pan_virus.sequence.ReconstructSequence()
	modified_gRNA = []
	name = pan_virus.name
	del pan_virus
	for gRNA in main_current_gRNAs:
		sequence = gRNA.seq 

		## Best : Lowest Mismatch score
		matches = bitapSearch(consensus, reverseC(sequence), mismatch, best = True)
		gRNA.ref_checked += 1
		
		if matches:
			best = matches.pop(0)

			## Cross-Match gRNAs to pan-virus genome
			if int(best.pos[0]) - 1 in list(gRNA_pos.keys()):
				g = gRNA_pos[best.pos[0] - 1]
				g.mismatch = best.editDistance
				g.Intolerant = best.intolerant
				g.cigar = best.cigar 
				g.name = name
				g.pos_consensus = best.pos[0]
				g.consensus_query = reverseC(best.query) 
				gRNA.pan_virus.append(g)
				gRNA.ref_found += 1
		
		modified_gRNA.append(gRNA)

	queue.put(modified_gRNA)

def _open_position_score(position_score):
	entropy_score = []
	conservation_score = []
	sequence = ''
	name = position_score.replace('_position_score.csv','')

	with open(position_score, 'r') as f:
		next(f)
		for line in f:
			line = line.strip('\n').split('\t')
			entropy = float(line[1])
			conservation = float(line[2])
			bp = line[3]
			entropy_score.append(entropy)
			conservation_score.append(conservation)
			sequence += bp 

	return Virus_Score(name, entropy_score, conservation_score, sequence.upper())

def _update_scores(pan_virus, gRNA_class_list, max_guide_score):
	for g in gRNA_class_list:
		start = g.pos[0]
		end = g.pos[1]
		g.e_score_list = pan_virus.entropy[start:end]
		g.c_score_list = pan_virus.conservation[start:end]
		g.scoring(max_guide_score, get_weighted_score = True)

	gRNA_class_list.sort(key=lambda x: x.weighted_score, reverse=True)
	for i, g in enumerate(gRNA_class_list):
		g.RANK = i+1

	return gRNA_class_list

def _multi_secondary_search(main_args, main_gRNA_class):
	from multiprocessing import Process, Queue, Manager
	mismatch = main_args.mismatch
	number_processors = main_args.CPU
	temp_files_to_remove = []
	manager = Manager()

	for position_file in main_args.position:
		logger.info(f'Processing (Getting On-Target Scores) for {position_file}')
		pan_virus = _open_position_score(position_file)
		handle = pan_virus.write_to_fasta()
		_, gRNA_class_list, temp_files_to_append, max_guide_score, _ = get_gRNAs(main_args, reference = handle)
		
		logger.info(f'Processing (Getting Conservation Scores) for {position_file}')
		gRNA_class_list = _update_scores(pan_virus, gRNA_class_list, max_guide_score)
		temp_files_to_remove += temp_files_to_append
		
		current_processes = []
		q = Queue()
		updated_gRNAs = manager.list()

		gRNA_pos = {int(g.pos[0]):g for g in gRNA_class_list}
		del gRNA_class_list
		logger.info(f'Processing (Finding main gRNAs on supplemented sequence) for {position_file}')
		for i in ranges(main_gRNA_class, number_processors):
			current_gRNAs = main_gRNA_class[i[0]:i[1]]
			if not current_gRNAs: continue
			#_search_gRNA(wavelet, mismatch, main_current_gRNAs, gRNA_pos, queue)
			pp = Process(target = _search_gRNA, args=(pan_virus, mismatch, current_gRNAs, gRNA_pos, q))
			current_processes.append(pp)

		for pp in current_processes:
			pp.start()

		sink = Process(target=queue_sam, args=(q, updated_gRNAs))
		sink.start()

		for pp in current_processes:
			pp.join()

		q.put('Done')
		sink.join()

		updated_gRNAs = list(updated_gRNAs)
		main_gRNA_class = updated_gRNAs
		del updated_gRNAs

	delete_files(temp_files_to_remove, main_args.temp)

	pan_viruses = pan_virus_quantity(main_args.position)

	temp = []
	for g in main_gRNA_class:
		g._get_final_rank(pan_viruses)
		temp.append(g)

	main_gRNA_class = temp
	del temp

	logger.info(f"Found {len(main_gRNA_class)} pan-viruses gRNAs")
	main_gRNA_class = [g for g in main_gRNA_class if g.final_rank > 0]
	main_gRNA_class.sort(key=lambda x: x.final_rank)

	return main_gRNA_class

def pan_virus_quantity(pan_virus):
	return len(pan_virus)	

def _write_gRNAs(args, main_gRNA_class, conservation_bowtie_mode = False):
	PAN_VIRUS = "_pan_virus_results.csv"
	i = 0
	with open(args.prefix + PAN_VIRUS, 'w') as write_all:
		write_all.write(",".join(HEADER_MAIN_1))
		write_all.write("," + ",".join(HEADER_MAIN_2))
		write_all.write("," + ",".join(HEADER_MSA))
		write_all.write("," + HEADER_WEIGHTED_SCORE)
		write_all.write("," + HEADER_RANK_SCORE)
		#if conservation_bowtie_mode:
		#	write_all.write("," + HEADER_BOWTIE_CONSERVATION) 
		write_all.write('\n\n')

		for g in main_gRNA_class:
			if g.final_rank < 0: continue
			i += 1

			write_all.write(f"Rank {i}! Final Rank Score: {g.final_rank}\n")
			write_all.write(f"{'*'*25}\n")
			write_all.write(f"{g.id},{g.seq},{g.pos + 2}")
			write_all.write(f",{g.reference},{g.g_score},{g.rankscore},{g.SGS},{g.quartile}")
			write_all.write(f",{g.consensus_query},{g.pos_consensus},{g.mismatch},{g.cigar}")
			write_all.write(f",{g.Intolerant},{g.e_sum},{g.c_mean},{g.c_SD},{g.e_sum_sensitive},{g.c_mean_sensitive},{g.c_SD_sensitive}")
			write_all.write(f",{g.weighted_score}")
			write_all.write(f",{g.rank}")
			#if conservation_bowtie_mode:
			#	write_all.write(f",{g.c_score_bowtie}")
			write_all.write('\n')

			for p in g.pan_virus:
				write_all.write(f"{p.id},{p.seq},{p.pos[0] + 2}")
				write_all.write(f",{p.reference},{round(p.g_score,5)},{round(p.rank,5)},{round(p.SGS,5)},{p.quartile}")
				write_all.write(f",{p.consensus_query},{p.pos_consensus[0] + 1},{p.mismatch},{p.cigar}")
				write_all.write(f",{p.Intolerant},{p.e_sum},{p.c_mean},{p.c_SD},{p.e_sum_sensitive},{p.c_mean_sensitive},{p.c_SD_sensitive}")
				write_all.write(f",{p.weighted_score}")
				write_all.write(f",{p.RANK}")
				write_all.write('\n')
			write_all.write(f"{'*'*25}\n\n")

def MSA_main(args, main_gRNA_class, conservation_bowtie_mode):
	if not args.position:
		raise AlignmentError('No Score Files are added. Please use Cas13gRNAscore to generate entropy and conservation score of virus')
	
	main_gRNA_class = _multi_secondary_search(args, main_gRNA_class)
	_write_gRNAs(args, main_gRNA_class, conservation_bowtie_mode)
	return

## Bowtie Mode

def _generate_gRNA_handle(main_gRNA_class, prefix):
	handle = prefix + '_pan-virus_gRNAs.fasta'
	with open(handle, "w") as f:
		for g in main_gRNA_class:
			f.write(f'>{g.id}\n{g.seq}\n')

	return handle

def _update_guidescore_rank(main_gRNA_class):
	## Update guidescore rank 
	main_gRNA_class.sort(key=lambda x: x.SGS, reverse=True)
	for i, g in enumerate(main_gRNA_class):
		g.guidescore_rank = (i+1)
	return main_gRNA_class

def _update_bowtie_class(samfile, main_gRNA_class_dict, sensitive_range = (2,7)):
	pattern = r"(\d+)([ACGT])"
	seq_count = 0
	sequences = []
	with open(samfile, 'r') as f:
		for line in f: 
			line = line.strip('\n')
			if line.startswith('@SQ'):
				match = re.findall(r"\tSN:(\S+)", line)
				name = match[0] if len(match) else "unknown_reference"
				sequences.append(name)
				seq_count += 1
				continue

			if line.startswith('@'): continue

			cols = line.split("\t")
			mapping_status = int(cols[1])
			if mapping_status == 16:
				sgrna_name = cols[0]
				genome_name = cols[2]
				MD_tag = cols[12].replace("MD:Z:", "")
				sensitive_count = 0
				intolerant = False
				if MD_tag != "0":
					for count, letter in re.findall(pattern, MD_tag):
						sensitive_count += int(count) 
						if sensitive_count in sensitive_range:
							intolerant = True
							break

				if intolerant: continue
				NM_tag = cols[13].replace("NM:i:", "")
				main_gRNA_class_dict[sgrna_name].count[int(NM_tag)] += 1
				main_gRNA_class_dict[sgrna_name].bowtie_list.append(genome_name)
				main_gRNA_class_dict[sgrna_name].bowtie_list_mm.append(int(NM_tag))

	PostCas13gRNA.update_bowtie_list(sequences)
	PostCas13gRNA.update_bowtie_count(seq_count)

	return main_gRNA_class_dict

def _bowtie_pan_virus(args, bowtie_reference, iteration, crRNA_fasta):
	bowtie_path, bowtie_build_path = check_bowtie(args)
	temp_files_to_remove = []

	_prefix = f"{args.prefix}-{iteration}.tmp"
	samfile = f"{_prefix}.sam"
	logger.info(f"Building Bowtie Index for {bowtie_reference}! This might take a while")
	bowtie_build_cmdline, bowtie_index_prefix = bowtie_build_cmd(bowtie_reference, _prefix, bowtie_build_path, args.CPU)
	success_bowtie_build = run_shell_command(bowtie_build_cmdline)
	if not success_bowtie_build:
		raise AlignmentError("Error during bowtie-build")

	temp_files_to_remove.append(_prefix + '_output_build.txt')
	temp_files_to_remove.append(_prefix + '_error_build.txt')

	logger.info(f"Bowtie Starting for {bowtie_reference}!")
	mismatch = args.mismatch if args.mismatch <= 3 else 3
	cmdline_bowtie = bowtie_cmd(bowtie_index_prefix, crRNA_fasta, _prefix, mismatch, bowtie_path, args.CPU)
	success_bowtie = run_shell_command(cmdline_bowtie)
	if not success_bowtie:
		raise AlignmentError("Error during Conservation Scoring Alignment")

	logger.info(f"Bowtie Done for {bowtie_reference}!")
	temp_files_to_remove.append(samfile)

	onlyfiles = [f for f in os.listdir() if f.startswith(f"{_prefix}_alignment.index")]
	temp_files_to_remove += onlyfiles

	return samfile, temp_files_to_remove

def _update_bowtie_rank(main_gRNA_class, mismatch = 0, update = True, ignore_efficacy = False):
	if update:
		for g in main_gRNA_class:
			g.update_bowtie_percent()

	main_gRNA_class.sort(key=lambda x: x.g_score, reverse=True)
	main_gRNA_class.sort(key=lambda x: x.bowtie_percent[mismatch], reverse=True)

	if ignore_efficacy:
		return main_gRNA_class

	if update:
		rank = 0 
		previous_percent = 0
		for g in main_gRNA_class:
			if g.bowtie_percent[0] != previous_percent:
				rank += 1

			g.update_bowtie_rank(rank)
			previous_percent = g.bowtie_percent[0]

	main_gRNA_class.sort(key=lambda x: x.gsXbw_rank, reverse=False)

	return main_gRNA_class

def _plot_bowtie_gRNA(args, best_gRNAs_dict):
	steps = range(1,args.best_count + 1)

	legend_elements = [Line2D([0], [0], color='red', marker='o', label="0 Mismatch", markersize = 10),
					   Line2D([0], [0], color="black", marker="^", label="1 Mismatch", markersize = 10),
					   Line2D([0], [0], color="blue", marker="h", label="2 Mismatch", markersize = 10),
					   Line2D([0], [0], color="green", marker="s", label="3 Mismatch", markersize = 10)]

	plot_dict = {k:[i[1] for i in v] for k,v in best_gRNAs_dict.items()}
	for k, v in plot_dict.items():
		if len(v) < args.best_count:
			diff = args.best_count - len(v)
			last_score = v[-1]
			for i in range(diff):
				v.append(last_score)
		plot_dict[k] = v

	mismatch_0_pan = plot_dict[0]
	mismatch_1_pan = plot_dict[1]
	mismatch_2_pan = plot_dict[2]
	mismatch_3_pan = plot_dict[3]

	plt.switch_backend('agg')
	fig, ax = plt.subplots()
	if mismatch_2_pan != mismatch_3_pan:
		ax.plot(steps, mismatch_3_pan, linewidth = 1, color="green", marker="s", markersize = 10)
	if mismatch_2_pan != mismatch_1_pan:
		ax.plot(steps, mismatch_2_pan, linewidth = 1, color="blue", marker="h", markersize = 10)
	if mismatch_0_pan != mismatch_1_pan:
		ax.plot(steps, mismatch_1_pan, linewidth = 1, color="black", marker="^", markersize = 10)
	ax.plot(steps, mismatch_0_pan, linewidth = 1, color="red", marker="o", markersize = 10)
	
	ax.set_ylabel("Seqeunces Targeted %", fontsize=18)
	ax.set_xlabel("gRNAs", fontsize=18)
	ax.set_xticks(steps)
	ax.tick_params(labelsize=14)

	ax.legend(handles=legend_elements, loc='center right', fontsize = 12)

	figure = plt.gcf() # get current figure
	figure.set_size_inches(8,6)
	plt.savefig(f'{args.prefix}_occurrences_plot.png', dpi = 360)

def _search_bowtie_gRNA(args, main_gRNA_class, ignore_efficacy = False):
	attempts = 0
	start = 0
	end = min(args.best_count, len(main_gRNA_class) - 1)
	best_gRNAs_dict = {}

	for m in range(PostCas13gRNA_counter.max_mismatch + 1):
		for g in main_gRNA_class:
			g.bowtie_temp_ct = 0
		
		main_gRNA_class = _update_bowtie_rank(main_gRNA_class, mismatch = m, update = False, ignore_efficacy = ignore_efficacy)
		
		remainder_count = PostCas13gRNA_counter.bowtie_count
		remainder_list = PostCas13gRNA_counter.bowtie_all_new
		best_gRNAs = []
		first_sequence_count = 0
		attempt = 0
		g_current = main_gRNA_class[0]
		while True:
			eliminate = [i[1] for i in list(zip(g_current.bowtie_list_mm, g_current.bowtie_list)) if i[0] <= m]
			remainder_list = [s for s in remainder_list if s not in eliminate]
			difference = remainder_count - len(remainder_list)
			first_sequence_count += difference
			first_percentage = round(first_sequence_count/PostCas13gRNA_counter.bowtie_count * 100, 10)
			best_gRNAs.append((g_current.id, round(first_percentage,5), g_current.seq, round(g_current.g_score,5)))
			remainder_count = len(remainder_list)
			if remainder_count == 0 and first_percentage == 100: 
				break

			attempt += 1
			if attempt == (end-1):
				break

			remainder_set = set(remainder_list)
			for g in main_gRNA_class:
				og_set = [i[1] for i in list(zip(g.bowtie_list_mm, g.bowtie_list)) if i[0] <= m]
				og_set = set(og_set)
				count = len(list(og_set & remainder_set))
				g.bowtie_temp_ct = count 

			main_gRNA_class.sort(key=lambda x: x.g_score, reverse=True)
			main_gRNA_class.sort(key=lambda x: x.bowtie_temp_ct, reverse=True)
			g_current = main_gRNA_class[0]
			if not g_current.bowtie_temp_ct:
				break

		best_gRNAs_dict[m] = best_gRNAs
		if remainder_list:
			for r in remainder_list:
				logger.info(f"{m} Mismatches: Failed to find gRNAs in Sequence ID: {r}")

	with open(f"{args.prefix}_best_bowtie_gRNAs.csv", "w") as f:
		f.write(",".join(HEADER_PAN_BOWTIE) + "\n")
		for k, v in best_gRNAs_dict.items():
			for i in v:
				f.write(f"{i[0]},{i[2]},{i[3]},{i[1]},{k}\n")

	if args.plot:
		_plot_bowtie_gRNA(args, best_gRNAs_dict)

def Bowtie_main(args, main_gRNA_class):
	temp_files_to_remove = []
	main_gRNA_class = _update_guidescore_rank(main_gRNA_class)
	gRNA_handle = _generate_gRNA_handle(main_gRNA_class, args.prefix)

	main_gRNA_class_dict = {g.id:g for g in main_gRNA_class}

	del main_gRNA_class
	for i, b in enumerate(args.bowtie):
		samfile, temp_files_to_append = _bowtie_pan_virus(args, b, i, gRNA_handle)
		main_gRNA_class_dict = _update_bowtie_class(samfile, main_gRNA_class_dict)
		temp_files_to_remove += temp_files_to_append

	delete_files(temp_files_to_remove, args.temp)
	main_gRNA_class = list(main_gRNA_class_dict.values())
	main_gRNA_class = _update_bowtie_rank(main_gRNA_class, ignore_efficacy=args.ignore_efficacy)
	_search_bowtie_gRNA(args, main_gRNA_class, ignore_efficacy=args.ignore_efficacy)
	return

def main():
	args = get_args()
	start_time = time.time()
	logging.basicConfig(level=args.loglevel, format=FORMAT)
	logger.info(f"{PROG} Starting!!")
	mode = args.mode

	main_gRNA_class, conservation_bowtie_mode = _open_scorefile(args.main, mainFile = True, mismatch_score = args.mismatch)
	if not main_gRNA_class: raise AlignmentError('Cannot Open Main Score File')
	
	if mode == 'MSA':
		MSA_main(args, main_gRNA_class, conservation_bowtie_mode)
	elif mode == 'bowtie':
		Bowtie_main(args, main_gRNA_class)

	end_time = time.time()
	total_time = round(end_time - start_time,3)
	logger.info(f"{PROG} has Ended in {total_time} seconds!")
	
	
if __name__ == '__main__':
	main()
