#!/usr/bin/env python3

import os
import logging
import logging.handlers
import time
import math
from functools import reduce
from cas13gRNAtor.utils import *
from cas13gRNAtor.cas13gRNAtor import get_gRNAs

logger = logging.getLogger(__name__)
	
PROG = 'Cas13gRNAevaluate'

DESCRIPTION = f'''
{PROG} allows users to check if Pan-virus gRNAs exist
'''

def get_args():
	parser = argparse.ArgumentParser(
		prog = PROG,
		description = DESCRIPTION
	)

	parser.add_argument('-r', '--reference', dest = 'main', required = True, 
						help="Guide Score for Main Reference")
	parser.add_argument('-f', '--file',  dest='position', required=True, nargs = '+',
						metavar="txt", help="accepts score files from Cas13gRNAscore module")
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

	args = parser.parse_args() 
	return args

def str2bool(v):
	return v.lower() in ("true", "1")

class PostCas13gRNA():

	## Main Score File
	def __init__(self, GuideName, GuideSeq, MatchPos, 
				 Reference, GuideScore, ScoreRanking, standardizedGuideScores, quartiles, 
				 Found, Sequence, Mismatch, Cigar, Intolerant, 
				 e_mean, e_SD, c_mean, c_SD, 
				 e_mean_sensitive, e_SD_sensitive, c_mean_sensitive, c_SD_sensitive, Weighted_Score,
				 Rank = float('inf'), Conservation_Bowtie = None, mainFile = True):
		
		self.id               = GuideName
		self.seq              = GuideSeq.upper()
		self.pos              = int(MatchPos)
		self.g_score          = float(GuideScore)
		self.rankscore        = float(ScoreRanking)
		self.SGS              = float(standardizedGuideScores)
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
		self.e_mean     	  = float(e_mean)
		self.c_mean_sensitive = float(c_mean_sensitive)
		self.e_mean_sensitive = float(e_mean_sensitive)
		self.c_SD 			  = float(c_SD)
		self.e_SD 			  = float(e_SD)
		self.c_SD_sensitive   = float(c_SD_sensitive)
		self.e_SD_sensitive   = float(e_SD_sensitive)
		self.weighted_score   = float(Weighted_Score)

		self.c_score_bowtie   = float(Conservation_Bowtie) if Conservation_Bowtie else None
		self.gRNA_len         = LENGTH_OF_GRNA

		## For Reference Sequence Search
		self.mainFile         = mainFile
		self.pan_virus        = []
		self.ref_checked      = 0 
		self.ref_found        = 0
		self.final_rank       = -1

	def _get_final_rank(self, pan_virus_quantity):
		if self.ref_found < pan_virus_quantity or self.ref_checked < pan_virus_quantity:
			return

		rank = []
		for g in self.pan_virus:
			rank.append(g.RANK)

		self.final_rank = reduce((lambda x, y: x * y), rank) * self.rank  
		self.final_rank = round(math.log10(self.final_rank),5)

def _open_scorefile(scoreFile, mainFile = True, mismatch_score = 5):
	col_len = len(HEADER_MAIN_1)  + len(HEADER_MAIN_2) + len(HEADER_MSA) + 1
	offtarget_mode = conservation_bowtie_mode = False
	rank_count = 0
	gRNA_list = []
	with open(scoreFile, 'r') as f:
		try:
			for line in f:
				if line.startswith('GuideName'):
					col = line.strip('\n').split(',')
					offtarget_mode = HEADER_BOWTIE_OFFTARGETS in col
					conservation_bowtie_mode = True if HEADER_BOWTIE_CONSERVATION in col else False
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
		g.scoring(max_guide_score, get_weighted_score = True, store = False)

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
		_, gRNA_class_list, temp_files_to_append, max_guide_score = get_gRNAs(main_args, reference = handle)
		
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
		if conservation_bowtie_mode:
			write_all.write("," + HEADER_BOWTIE_CONSERVATION) 
		write_all.write('\n\n')

		for g in main_gRNA_class:
			if g.final_rank < 0: continue
			i += 1

			write_all.write(f"Rank {i}! Final Rank Score: {g.final_rank}\n")
			write_all.write(f"{'*'*25}\n")
			write_all.write(f"{g.id},{g.seq},{g.pos}")
			write_all.write(f",{g.reference},{g.g_score},{g.rankscore},{g.SGS},{g.quartile}")
			write_all.write(f",{g.found},{g.consensus_query},{g.mismatch},{g.cigar}")
			write_all.write(f",{g.Intolerant},{g.e_mean},{g.e_SD},{g.c_mean},{g.c_SD},{g.e_mean_sensitive},{g.e_SD_sensitive},{g.c_mean_sensitive},{g.c_SD_sensitive}")
			write_all.write(f",{g.weighted_score}")
			write_all.write(f",{g.rank}")
			if conservation_bowtie_mode:
				write_all.write(f",{g.c_score_bowtie}")
			write_all.write('\n')

			for p in g.pan_virus:
				write_all.write(f"{p.id},{p.seq},{p.pos[0] + 1}")
				write_all.write(f",{p.reference},{p.g_score},{p.rank},{p.SGS},{p.quartile}")
				write_all.write(f",{p.found},{p.consensus_query},{p.mismatch},{p.cigar}")
				write_all.write(f",{p.Intolerant},{p.e_mean},{p.e_SD},{p.c_mean},{p.c_SD},{p.e_mean_sensitive},{p.e_SD_sensitive},{p.c_mean_sensitive},{p.c_SD_sensitive}")
				write_all.write(f",{p.weighted_score}")
				write_all.write(f",{p.RANK}")
				write_all.write('\n')
			write_all.write(f"{'*'*25}\n\n")

def main():
	args = get_args()
	start_time = time.time()
	logging.basicConfig(level=args.loglevel, format=FORMAT)
	logger.info(f"{PROG} Starting!!")

	main_gRNA_class, conservation_bowtie_mode = _open_scorefile(args.main, mainFile = True, mismatch_score = args.mismatch)
	if not main_gRNA_class: raise AlignmentError('Cannot Open Main Score File')
	
	main_gRNA_class = _multi_secondary_search(args, main_gRNA_class)
	_write_gRNAs(args, main_gRNA_class, conservation_bowtie_mode)

	end_time = time.time()
	total_time = round(end_time - start_time,3)
	logger.info(f"{PROG} has Ended in {total_time} seconds!")
	
if __name__ == '__main__':
	main()
