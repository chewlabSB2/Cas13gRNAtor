#!/usr/bin/env python3

"""Burrows-Wheeler Aligner - String Search"""

class simpleBWA:
	"""
	A Burrows-Wheeler Aligner class. 
	Contains the 4 core datastructures used in the algorithm, SA, BWT, C and Occ which are created in the constructor which is passed the reference string as an argument
	"""

	alphabet = ["A", "C", "G", "T", "N"]
	# Initialize the Datastructure 
	def __init__(self, reference):
		reference = reference.upper()
		reverse_reference = reference[::-1]
		rotation_list, rotation_list_reverse, suffix_array, bwt = [list() for i in range(4)]
		C, Occ, Occ_reverse = [dict() for i in range(3)]
	
		#initialize 2 auxillary datastructures
		for char in self.alphabet:
			C[char] = 0
			Occ[char] = list()# in Occ, each character has an associated list of integer values (for each index along the reference)
			Occ_reverse[char] = list()
	
		#append the ending character to the reference string
		reference = "%s$" % reference
		reverse_reference = "%s$" % reverse_reference

		#create all the rotation/suffix combinations of the reference and reverse reference, and their starting index positions
		for i in range(len(reference)):
			new_rotation = "%s%s" % (reference[i:],reference[0:i])
			struct = Suffix(new_rotation,i)
			rotation_list.append(struct)
			
			new_rotation_reverse = "%s%s" % (reverse_reference[i:],reverse_reference[0:i])
			struct_rev = Suffix(new_rotation_reverse,i)
			rotation_list_reverse.append(struct_rev)
		
			#create the C datastructure. C(a) = the number of characters 'a' in the Reference that are lexographically smaller than 'a'
			#NOTE, the C datastructure is not required for the reverse reference
			if reference[i]!='$':
				for char in self.alphabet:
					if reference[i] < char:
						C[char] = C[char] + 1   
		
		#sort the rotations/suffixes using the suffix/rotation text as the key
		rotation_list.sort(key=textKey)
		rotation_list_reverse.sort(key=textKey)
	
		#now record the results into 2 seperate lists, the suffix (or S) array and the BWT (or B) array
		#also calculate the auxilliary datastructure Occ (or O)
		for i in rotation_list:
			suffix_array.append(i.pos)#the position of the reordered suffixes forms the Suffix Array elements
			bwt.append(i.text[-1:])#the last character in each rotation (in the new order) forms the BWT string elements
		
			#now construct the Occ (or C) datastructure
			for char in self.alphabet:
				if len(Occ[char]) == 0:
					prev = 0
				else:
					prev = Occ[char][-1]
				if i.text[-1:] == char:
					Occ[char].append(prev+1)
				else:
					Occ[char].append(prev)
					
		#now record the results into 2 seperate lists, the suffix (or S) array and the BWT (or B) array
		#also calculate the auxilliary datastructures, C and Occ (or O)
		for i in rotation_list_reverse:
			#construct the Occ (or C) datastructure
			for char in self.alphabet:
				if len(Occ_reverse[char]) == 0:
					prev = 0
				else:
					prev = Occ_reverse[char][-1]
				if i.text[-1:] == char:
					Occ_reverse[char].append(prev+1)
				else:
					Occ_reverse[char].append(prev)                  
					
		## Store all the useful datastructures as class variables for easy future access
		self.SA = suffix_array
		self.BWT = bwt
		self.C = C
		self.Occ = Occ
		## Occ datastructure for the reverse reference, using to construct the D array (the lower bound on the number of differences allowed), to speed up alignments 
		self.Occ_reverse = Occ_reverse 
		self.n = len(reference)
		## empty list for later use
		self.D = list()

	#get the position(s) of the query in the reference
	def find_match(self,query,num_differences):
		query = query.upper()
		if num_differences == 0:
			return self.exact_match(query)
		else:
			return self.inexact_match(query,num_differences)

	#exact matching - no indels or mismatches allowed
	def exact_match(self, query):
		query = query.upper()
		i = 0
		j = self.n - 1
		
		for x in range(len(query)):
			newChar = query[-x-1]
			newI = self.C[newChar] + self.OCC(newChar,i-1) + 1
			newJ = self.C[newChar] + self.OCC(newChar,j)
			i = newI
			j = newJ
		matches = self.SA[i:j+1]
		return matches

	## inexact matching, z is the max threshold for allowed edits
	def inexact_match(self,query,z):
		self.calculate_d(query)
		SA_indices = self.inexact_recursion(query, len(query)-1, z, 0, self.n-1)
		## Return the values in the SA
		return [self.SA[x] for x in SA_indices]

	## recursion function that effectively "walks" through the suffix tree using the SA, BWT, Occ and C datastructures
	def inexact_recursion(self,query,i,z,k,l):
		tempset = set()
			
		## 2 stop conditions, one when too many differences have been encountered, another when the entire query has been matched, terminating in success
		## Reached the limit of differences at this stage, terminate this path traversal
		if (z < self.get_D(i) and use_lower_bound_tree_pruning) or (z < 0 and not use_lower_bound_tree_pruning):
			## Return empty set   
			return set()
		## Empty query string, entire query has been matched, return SA indexes k:l
		if i < 0:
			for m in range(k,l+1):
				tempset.add(m)
			return tempset
			
		result = set()
		if indels_allowed: result = result.union(self.inexact_recursion(query,i-1,z-insertion_penalty,k,l))#without finding a match or altering k or l, move on down the query string. Insertion
		## For each character in the alphabet
		for char in self.alphabet:
			## Find the SA interval for the char
			newK = self.C[char] + self.OCC(char,k-1) + 1 
			newL = self.C[char] + self.OCC(char,l)
			## If the substring was found
			if newK <= newL:
				if indels_allowed: result = result.union(self.inexact_recursion(query,i,z-deletion_penalty,newK,newL))# Deletion
				
				## If the char was correctly aligned, then continue without decrementing z (differences)
				if char == query[i]: result = result.union(self.inexact_recursion(query,i-1,z,newK,newL))
				## Continue but decrement z, to indicate that this was a difference/unalignment
				else: result = result.union(self.inexact_recursion(query,i-1,z-mismatch_penalty,newK,newL))
		return result

	#calculates the D array for a query, used to prune the tree walk and increase speed for inexact searching
	def calculate_d(self,query):
		k = 0
		l = self.n-1
		z = 0

		## Empty the D array
		self.D = list()
		for i in range(len(query)):
			k = self.C[query[i]] + self.OCC(query[i],k-1,reverse=True) + 1
			l = self.C[query[i]] + self.OCC(query[i],l,reverse=True)
			## If this character has NOT been found
			if k > l:
				k = 0
				l = self.n - 1
				z = z + 1
			self.D.append(z)

	## returns normal Occ value, otherwise returns the reverse Occ if explicitly passed as an argument
	## NOTE Occ('a',-1) = 0 for all 'a'
	def OCC(self,char,index,reverse=False):
		if index < 0:
			return 0
		else:
			if reverse:
				return self.Occ_reverse[char][index]
			else:
				return self.Occ[char][index]
	
	## gets values from the D array
	## NOTE D(-1) = 0
	def get_D(self,index):
		if index < 0:
			return 0
		else:
			return self.D[index]

class Suffix:
	"""
	A simple class with 2 variables, used for sorting and calculating the Suffix Array and BWT array
	Each instance holds the position of the suffix and the suffix (text) itself
	"""
	def __init__(self, text, position):
		self.text = text
		self.pos = position

#this is used to sort the Suffix objects, according to their text key
def textKey( a ): return a.text

#environment variables
debug = True
show_data_structures = True
use_lower_bound_tree_pruning = True #set this to false (in conjunction with debug=True) to see the full search through the suffix trie
#search parameters
indels_allowed = False # turn off for mismatches only, no insertion or deletions allowed
#difference_threshold = 0
insertion_penalty = 1
deletion_penalty = 1
mismatch_penalty = 1
