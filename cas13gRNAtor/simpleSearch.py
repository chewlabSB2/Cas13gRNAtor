#!/usr/bin/env python3

from bitarray import bitarray
import numpy as np

def reverseC(seq):
	RT = {'A':'T','C':'G','G':'C','T':'A', 'N':'N'}
	reverseComplement = ''
	for i in seq:
		nt = RT.get(i)
		reverseComplement += nt
	return reverseComplement[::-1]


class Node:
	def __init__(self, value = []):
		self.left  = None #0 
		self.right = None #1
		self.value = bitarray(''.join([str(i) for i in value]))

class Wavelet_Tree():
	
	def __init__(self, sequence):
		self.ALPHABET = list(set(sequence)) 
		self.root = None
		self.length = len(sequence)
		self.CharBit = {k:None for k in self.ALPHABET}
		x = {k:sequence.count(k) for k in self.ALPHABET}
		self.nt_count = {k:v for k,v in x.items() if v} 
		self.ConstructNode(sequence, alphabet = self.ALPHABET)
		#self.printTree()

	def ReconstructSequence(self):
		sequence = ''
		for i in range(self.length):
			sequence += self.Access(i)
		return sequence

	def ConstructNode(self, string, node = None, alphabet = None, charbit = '', root = True):
		if len(alphabet) == 1:
			self.CharBit[alphabet[0]] = charbit
			return 

		charbitL = charbitR = charbit
		charbitL += '0'
		charbitR += '1'
		Bitmap = []
		x = {}

		for char in alphabet:
			if char not in list(self.nt_count.keys()): continue
			x[char] = self.nt_count[char]

		alphabet.sort(reverse = True)
		keys = alphabet
		mid = len(keys)//2
		SplitCharL = keys[:mid]
		SplitCharR = keys[mid:]
		Sleft = Sright = ''

		for c in list(string):
			if c in SplitCharL:
				Sleft += c
				Bitmap.append(0)
			else:
				Sright += c
				Bitmap.append(1)

		newNode = Node(Bitmap)
		if not self.root:
			self.root = newNode

		if charbit: 
			current_charbit = charbit[-1]
			if current_charbit == '0':
				node.left = newNode
			elif current_charbit == '1':
				node.right = newNode

		self.ConstructNode(Sleft, newNode, SplitCharL, charbitL, root = False)
		self.ConstructNode(Sright, newNode, SplitCharR, charbitR, root = False)
		if not root: return
		
		self.CharBit = {k:v for k,v in self.CharBit.items() if v}
		self.CharBit_inv = {v:k for k, v in self.CharBit.items()}

	def printTree(self):
		if self.root is not None:
			self._printTree(self.root)

	def _printTree(self, node, left = True):
		if node is not None:
			self._printTree(node.left, left)
			print (left, node.value.__sizeof__(), str(node.value) + ' ')
			self._printTree(node.right, False)

	def _BinaryRank(self, c, p, node):
		return node.value[:p].count(c)

	def Rank(self, c, p):
		'''	
		c = Character
		p = position/ offset
		'''
		if p < 0:
			return 0
		charbit = list(self.CharBit[c])
		max_counter = len(charbit)
		i = 0
		node = None
		while True: 
			if not node:
				node = self.root
			else:
				if current_charbit == 0:
					node = node.left
				elif current_charbit == 1:
					node = node.right
			
			current_charbit = int(charbit[i])
			p = self._BinaryRank(current_charbit, p, node)
			i += 1
			if i == len(charbit):
				break
				
		return p

	def _BinaryAccess(self, node, o):
		return node.value[o]

	def Access(self, offset):
		charbit = ''
		node = None
		while True:
			if not node:
				node = self.root
			o = int(self._BinaryAccess(node, offset))
			charbit += str(o)

			if charbit in list(self.CharBit_inv.keys()):
				return self.CharBit_inv[charbit]

			r = int(self._BinaryRank(o, offset, node))
			offset = r 
			if o == 0:
				node = node.left
			elif o == 1:
				node = node.right
	
	def Select(self):
		pass

class Cas13gRNA_temp():
	def __init__(self, reference, query, pos, mismatch):
		self.name = None
		self.query = query.upper()
		self.reference = reference
		self.pos = pos
		self.mismatch = mismatch
		self.cigar = None
		self.intolerant = None
		self.editDistance = None 
		self._edit_distance(reference, self.query)

	def _edit_distance(self, string_x=None, string_y=None, get_backtrace=True):
		string_x = string_x.upper()
		string_y = string_y.upper()

		string_x = reverseC(string_x)
		string_y = reverseC(string_y)
		
		matrix_a = np.zeros([2, len(string_y) + 1], dtype = int)
		matrix_a[0, :] = np.arange(0,len(string_y)+1,1)
		col = np.arange(0,len(string_x)+1,1)

		'''
		default       : 000 #match
		deletion      : 100 
		substituition : 010 #mismatch
		insertion     : 001 
		'''
		B = None
		if get_backtrace:
			#B = [['000' if x > 0 else '100' for x in range(len(string_x) + 1)] for y in range(len(string_y) + 1) ]
			#B[0] = ['000' if x == 0 else '001' for x in range(len(string_x) + 1)]
			B = [[bitarray('000') if x > 0 else bitarray('100') for x in range(len(string_x) + 1)] for y in range(len(string_y) + 1)]
			B[0] = [bitarray('000') if x == 0 else bitarray('001') for x in range(len(string_x) + 1)]
		
		for i in range(1,len(col)):
			matrix_a[1, 0] = col[i]
			for j in range(1,len(matrix_a[1,:])):
				deletion = matrix_a[0, j] + 1
				insertion = matrix_a[1, j - 1] + 1
				if (string_x[i - 1] == string_y[j - 1]):
					substitution = matrix_a[0, j - 1]
				else:
					substitution = matrix_a[0, j - 1] + 1
				score = min(([deletion, insertion, substitution]))	
				matrix_a[1, j] = score
				if get_backtrace:
					temp = [score == deletion, score == substitution, score == insertion]
					B[i][j] = bitarray(temp) #''.join(['1' if x else '0' for x in temp])
			matrix_a[0, :] = matrix_a[1, :]

		self.editDistance = matrix_a[-1, -1]
		
		if get_backtrace:
			bt = self._naive_backtrace(B)
			self.cigar, self.intolerant = self._align_gRNA(string_x, string_y, bt)

		#return editDistance

	def _naive_backtrace(self, B_matrix):
		i, j = len(B_matrix[0])-1, len(B_matrix)-1
		backtrace_idxs = [(i, j)]

		while (i, j) != (0, 0):
			if B_matrix[i][j][1]:
				i, j = i-1, j-1
			elif B_matrix[i][j][0]:
				i, j = i-1, j
			elif B_matrix[i][j][2]:
				i, j = i, j-1
			backtrace_idxs.append((i,j))

		return backtrace_idxs

	def _align_gRNA(self, string_x, string_y, bt, sensitive_region = range(15,22)):
		'''
		string_x - reference gRNA
		string_y - gRNA
		return cigar/ get position of Indels
		'''
		op = []
		mutated_pos = [] #relative to reference
		backtrace = bt[::-1]
		
		for k in range(len(backtrace)-1): 
			x0, y0 = backtrace[k]
			x1, y1 = backtrace[k+1]

			if x1 > x0 and y1 > y0:
				if string_x[x0] == string_y[y0]:
					op.append('M')
				else:
					mutated_pos.append(x0)
					op.append(string_x[x0])
			elif x0 == x1:
				mutated_pos.append(y0)
				op.append('I')
			elif y0 == y1:
				mutated_pos.append(x0)
				op.append('D')

		operations = ['M', 'I', 'D']
		current_op = None
		count = 0 
		cigar = ''
		for char in op:
			if char in operations:
				if char != current_op:
					if count >= 1: cigar += f'{count}{current_op}'
					current_op = char
					count = 1
				else:
					count += 1
			else:
				if current_op:
					if count >= 1: cigar += f'{count}{current_op}'
					count = 0 
					current_op = None
				cigar += f'{char}'
		else:
			if count >= 1 and current_op: cigar += f'{count}{current_op}'

		intolerant = False
		for pos in mutated_pos:
			if pos in sensitive_region:
				intolerant = True
				break

		return cigar, intolerant

def _generateAlphabet(reference, query):
	alphabet = list(set(reference))
	bitap_dict = {}
	for letter in alphabet:
		letterPositionInQuery = 0
		for symbol in query:
			letterPositionInQuery = letterPositionInQuery << 1
			letterPositionInQuery |= int(letter != symbol)
		bitap_dict[letter] = letterPositionInQuery
	return bitap_dict

def bitapSearch(reference, query, mismatch = 10, best = False):
	referenceLen = len(reference)
	queryLen = len(query)

	alphabet = _generateAlphabet(reference, query)

	matrix = [] 
	emptyColumn = (2 << (queryLen - 1)) - 1
	underground = [emptyColumn for i in range(referenceLen + 1)]
	matrix.append(underground)
	gRNAs = []
	skip = []

	for k in range(1, mismatch + 2):
		matrix.append([emptyColumn])

		for columnNum in range(1, referenceLen + 1):
			prevColumn = (matrix[k][columnNum - 1]) >> 1
			letterPattern = alphabet[reference[columnNum - 1]]
			curColumn = prevColumn | letterPattern

			if k > 1:
				#insertColumn = curColumn & (matrix[k - 1][columnNum - 1])
				#deleteColumn = curColumn & (matrix[k - 1][columnNum] >> 1)
				#replaceColumn = curColumn & (matrix[k - 1][columnNum - 1] >> 1)
				#curColumn = insertColumn & deleteColumn & replaceColumn
				
				curColumn = curColumn & (matrix[k - 1][columnNum - 1] >> 1)
				
			matrix[k].append(curColumn)

			if (curColumn & 0x1) == 0:
				startPos = max(0, columnNum - queryLen) # taking in account Replace operation
				if startPos in skip: continue
				endPos = min(columnNum, referenceLen) # taking in account Replace operation
				place = reference[startPos:endPos]
				temp = Cas13gRNA_temp(query, place, (startPos, endPos), k - 1)
				gRNAs.append(temp)
				if best: return gRNAs
				skip.append(startPos)
			
	return gRNAs