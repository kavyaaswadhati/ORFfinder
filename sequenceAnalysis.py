#!/usr/bin/env python3
# Name: Kavya Aswadhati (kaswadha)
# Group Members: None
'''
Read sequences from fasta files, compile them into a genome, calculate composition of genome in nucleotides,
amino acids, codons
'''

from collections import defaultdict

class NucParams:
	'''
	Class that contains the methods that will compile and analyze inputted genome.
	'''
	rnaCodonTable = {
		# RNA codon table
		# U
		'UUU': 'F', 'UCU': 'S', 'UAU': 'Y', 'UGU': 'C',  # UxU
		'UUC': 'F', 'UCC': 'S', 'UAC': 'Y', 'UGC': 'C',  # UxC
		'UUA': 'L', 'UCA': 'S', 'UAA': '-', 'UGA': '-',  # UxA
		'UUG': 'L', 'UCG': 'S', 'UAG': '-', 'UGG': 'W',  # UxG
		# C
		'CUU': 'L', 'CCU': 'P', 'CAU': 'H', 'CGU': 'R',  # CxU
		'CUC': 'L', 'CCC': 'P', 'CAC': 'H', 'CGC': 'R',  # CxC
		'CUA': 'L', 'CCA': 'P', 'CAA': 'Q', 'CGA': 'R',  # CxA
		'CUG': 'L', 'CCG': 'P', 'CAG': 'Q', 'CGG': 'R',  # CxG
		# A
		'AUU': 'I', 'ACU': 'T', 'AAU': 'N', 'AGU': 'S',  # AxU
		'AUC': 'I', 'ACC': 'T', 'AAC': 'N', 'AGC': 'S',  # AxC
		'AUA': 'I', 'ACA': 'T', 'AAA': 'K', 'AGA': 'R',  # AxA
		'AUG': 'M', 'ACG': 'T', 'AAG': 'K', 'AGG': 'R',  # AxG
		# G
		'GUU': 'V', 'GCU': 'A', 'GAU': 'D', 'GGU': 'G',  # GxU
		'GUC': 'V', 'GCC': 'A', 'GAC': 'D', 'GGC': 'G',  # GxC
		'GUA': 'V', 'GCA': 'A', 'GAA': 'E', 'GGA': 'G',  # GxA
		'GUG': 'V', 'GCG': 'A', 'GAG': 'E', 'GGG': 'G'  # GxG
	}


	def __init__(self, inString=''):
		'''
		:param inString: optional fasta file input
		initialize dictionaries for nucleotide composition and amino acid composition
		call methods that will analyze and populate dictionaries
		'''
		# initialize
		self.nucComp = {
			'A': 0, 'T': 0, 'C': 0, 'G': 0, 'U': 0, 'N': 0
		}
		self.aaComp = {
			'A': 0, 'C': 0, 'D': 0, 'E': 0, 'F': 0, 'G': 0, 'H': 0, 'I': 0, 'L': 0,
			'K': 0, 'M': 0, 'N': 0, 'P': 0, 'Q': 0, 'R': 0, 'S': 0, 'T': 0, 'V': 0, 'Y': 0, 'W': 0, '-': 0
		}
		self.codonComp = {}
		self.seqCodons = {}

		#populate
		self.addSequence(inString)
		self.codonComposition()
		self.nucComposition()

	def addSequence(self, inSeq):
		'''
		initializes inputted sequences by populating dictionaries of codon counts and and nucleotide counts
		'''
		# convert the sequence that was taken in to upper case for consistency
		toProcess = inSeq.upper()
		#append nucleotide count dict with each additional nuc in inputted seq
		#only allows for the addition of acceptable nuc values (ATGCUN)
		for char in toProcess:
			self.nucComp[char] += 1
		#print('nucComp:', self.nucComp)
		# split the sequence into a list of codons (sets of 3 chars)
		codons = [toProcess[i:i + 3] for i in range(0, len(toProcess), 3)]
		#print (codons)
		# add each codon to a dictionary with key(codon) value(# of instances)
		for codon in codons:
			if codon in self.seqCodons:
				#print('new codon')
				self.seqCodons[codon]+=1
			else:
				self.seqCodons[codon]=1
		return

	def aaComposition(self):
		'''
		Utilize rnaCodonTable to translate codons to amino acids
		Return dictionary of translated amino acids
		'''
		for codon in self.codonComp:
			if codon in NucParams.rnaCodonTable:
				self.aaComp[NucParams.rnaCodonTable.get(codon)]+= self.codonComp[codon]
			else:
				print("ERROR: codon not found",codon)
		return self.aaComp

	def nucComposition(self):
		'''
		Dictionary of nucleotide composition populated as sequences are added.
		Simply returns dictionary of nucleotide composition of genome
		'''
		# print('ncount:',self.nucComp)
		# self.gc = ((self.nucComp['G']+self.nucComp['C'])/self.nucCount())*100
		return self.nucComp


	def codonComposition(self):
		'''
		"cleans" codons that were added as sequences are added to genome
		Converts DNA to RNA codons, ensures codons are valid (exist in rnaCodonTable)
		'''
		for codon in self.seqCodons:
			newCodon = codon.replace('T','U')
			if newCodon in NucParams.rnaCodonTable:
				self.codonComp[newCodon]=self.seqCodons[codon]
		return

	def nucCount(self):
		'''
		Returns total number of nucleotides in genome.
		'''
		return sum(self.nucComp.values())


class FastAReader:
	def __init__(self, fname=''):
		'''contructor: saves attribute fname '''
		self.fname = fname

	def doOpen(self):
		''' Handle file opens, allowing STDIN.'''
		if self.fname is '':
			return sys.stdin
		else:
			return open(self.fname)

	def readFasta(self):
		''' Read an entire FastA record and return the sequence header/sequence'''
		header = ''
		sequence = ''

		with self.doOpen() as fileH:
			header = ''
			sequence = ''

			# skip to first fasta header
			line = fileH.readline()
			while not line.startswith('>'):
				line = fileH.readline()
			header = line[1:].rstrip()

			for line in fileH:
				if line.startswith('>'):
					yield header, sequence
					header = line[1:].rstrip()
					sequence = ''
				else:
					sequence += ''.join(line.rstrip().split()).upper()
		yield header, sequence


