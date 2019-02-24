#!/usr/bin/env python3
# Name: Kavya Aswadhati (kaswadha)
# Group Members: None
'''
Read sequences from fasta files, compile them into a genome, calculate composition of genome in nucleotides,
amino acids, codons.

'''

import sys
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
		Initialize dictionaries for nucleotide composition and amino acid composition
		call methods that will analyze and populate dictionaries.

		Args:
		     inString: optional fasta file input
		'''
		# dictionaries to be populated
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
		Initialize inputted sequences by populating dictionaries of codon counts and and nucleotide counts.

		Args:
		    inSeq: Uncleaned string of nucleotides

		Returns:
		    None
		'''
		# convert the sequence that was taken in to upper case for consistency
		toProcess = inSeq.upper()

		#append nucleotide count dict with each additional nuc in inputted seq
		#only allows for the addition of acceptable nuc values (ATGCUN)
		for char in toProcess:
			self.nucComp[char] += 1

		# split the sequence into a list of codons (sets of 3 chars)
		codons = [toProcess[i:i + 3] for i in range(0, len(toProcess), 3)]

		# add each codon to a dictionary with key(codon) value(# of instances)
		for codon in codons:
			if codon in self.seqCodons:
				self.seqCodons[codon]+=1
			else:
				self.seqCodons[codon]=1
		return

	def aaComposition(self):
		'''
		Utilize rnaCodonTable to translate codons to amino acids.

		Args:
		    self
		Return:
		     self.aaComp: dictionary of translated amino acids.
		'''
		for codon in self.codonComp:
			if codon in NucParams.rnaCodonTable:
				self.aaComp[NucParams.rnaCodonTable.get(codon)]+= self.codonComp[codon]
			else:
				print("ERROR: codon not found",codon)
		return self.aaComp

	def nucComposition(self):
		'''
		Pass dictionary of nucleotide composition populated as sequences are added in addSequence.
		Args:
		    self
		Return:
		    self.nucComp: dictionary of nucleotide composition of genome
		'''
		return self.nucComp


	def codonComposition(self):
		'''
		'Clean' codons: convert DNA to RNA codons, ensure codons are valid
		(exist in rnaCodonTable) before they are added to the genome.

		Args:
		    self

		Return:
		   None
		'''
		for codon in self.seqCodons:
			# convert to RNA codons
			newCodon = codon.replace('T','U')
			if newCodon in NucParams.rnaCodonTable:
				self.codonComp[newCodon]=self.seqCodons[codon]
		return

	def nucCount(self):
		'''
		Return: total number of nucleotides in genome.
		'''
		return sum(self.nucComp.values())


class FastAReader:
	'''
	Open FASTA files and segregate genome into genes with associated headers.

	input: .fa file or STDIN
	output: str header, str sequence
	'''
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


class ProteinParam:
	# These tables are for calculating:
	#     molecular weight (aa2mw), along with the mol. weight of H2O (mwH2O)
	#     absorbance at 280 nm (aa2abs280)
	#     pKa of positively charged Amino Acids (aa2chargePos)
	#     pKa of negatively charged Amino acids (aa2chargeNeg)
	#     and the constants aaNterm and aaCterm for pKa of the respective termini


	aa2mw = {
		'A': 89.093, 'G': 75.067, 'M': 149.211, 'S': 105.093, 'C': 121.158,
		'H': 155.155, 'N': 132.118, 'T': 119.119, 'D': 133.103, 'I': 131.173,
		'P': 115.131, 'V': 117.146, 'E': 147.129, 'K': 146.188, 'Q': 146.145,
		'W': 204.225, 'F': 165.189, 'L': 131.173, 'R': 174.201, 'Y': 181.189
	}

	# allAA = {'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'L', 'K', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'Y', 'W'}
	mwH2O = 18.015
	aa2abs280 = {'Y': 1490, 'W': 5500, 'C': 125}

	aa2chargePos = {'K': 10.5, 'R': 12.4, 'H': 6}
	aa2chargeNeg = {'D': 3.86, 'E': 4.25, 'C': 8.33, 'Y': 10}
	aaNterm = 9.69
	aaCterm = 2.34

	aaComp = {
		'A': 0, 'C': 0, 'D': 0, 'E': 0, 'F': 0, 'G': 0, 'H': 0, 'I': 0, 'L': 0,
		'K': 0, 'M': 0, 'N': 0, 'P': 0, 'Q': 0, 'R': 0, 'S': 0, 'T': 0, 'V': 0, 'Y': 0, 'W': 0
	}

	def __init__(self, protein):
		'''
		Initialize a protein object, as represented by constitutive aa chars,
		in a protein string into a dictionary aaComp.
		'''
		# make sure all the chars are uppercase
		sequence = protein.upper()
		# for each char in the sequence, if it is a valid aa char, increment the value assigned to the char key
		for acid in sequence:
			if acid in ProteinParam.aaComp:
				ProteinParam.aaComp[acid] += 1

	def aaCount(self):
		'''
		Count the number of valid amino acids in the protein.

		Return:
		    numAcids: number of amino acids in the the protein
		'''
		# sum the values of each key in aaComp
		numAcids = 0
		for acid in ProteinParam.aaComp:
			numAcids += ProteinParam.aaComp[acid]
		return numAcids

	def pI(self):
		'''
		Binary search for the pI (isoelectric point, charge = 0) of the protein.

		Return:
		     pI (accuracy (+/-).01)
		'''
		mid = None
		first = 0.0
		last = 14.0

		# As pH increases, so does charge. Finds the charge of a midpoint and depending on whether
		# it is below, above, or close to 0, either raise the floor or lower the roof (respectively)
		# and probe for charge again.
		while (last-first)>=.01:
			mid = (first + last) / 2
			thisCharge = self._charge_(mid)
			if thisCharge < 0:
				last = mid
			else:
				first = mid
		return float(mid)

	def aaComposition(self):
		'''Just pass the dictionary created in __init__'''
		return ProteinParam.aaComp

	def _charge_(self, pH):
		'''
		Find the charge of the protein at given pH (parameter)
		Args:
		    pH: int to calculate at

		Return:
		    netCharge: float of the net charge of the protein at the given pH
		'''
		# begin with the charge of each terminus
		netPosCharge = (((math.pow(10, ProteinParam.aaNterm)) / (math.pow(10, ProteinParam.aaNterm) + math.pow(10, pH))))
		netNegCharge = (((math.pow(10, pH)) / (math.pow(10, ProteinParam.aaCterm) + math.pow(10, pH))))

		# add the charges of negatively charged aas and positively charged aas seperately.
		for acid in ProteinParam.aaComp:
			if acid in ProteinParam.aa2chargeNeg:
				netNegCharge += (ProteinParam.aaComp[acid] * ((math.pow(10, pH)) / (math.pow(10, ProteinParam.aa2chargeNeg[acid]) + math.pow(10, pH))))
			elif acid in ProteinParam.aa2chargePos:
				netPosCharge += (ProteinParam.aaComp[acid] * ((math.pow(10, ProteinParam.aa2chargePos[acid])) / (math.pow(10, ProteinParam.aa2chargePos[acid]) + math.pow(10, pH))))

		# combine charges for netCharge
		netCharge = (netPosCharge - netNegCharge)
		return netCharge

	def molarExtinction(self, Cysteine=True):

		'''
		Find the molar extinction coefficient for the protein sequence.
		Args:
		    Optional parameter Cysteine for reducing conditions.

		Return:
		    extinction: molar extinction coefficient (measurement of how much light the protein absorbs at
		    a given wavelength, an intrinsic property)
		'''
		extinction = 0
		for acid in ProteinParam.aa2abs280:
			extinction += ((ProteinParam.aa2abs280[acid]) * ProteinParam.aaComp[acid])
		if not Cysteine:
			extinction -= ((ProteinParam.aa2abs280['C']) * ProteinParam.aaComp['C'])
		return extinction

	def massExtinction(self, Cysteine=True):
		'''
		Use the molar extinction coefficient to find the mass extinction coefficient.
		'''
		myMW = self.molecularWeight()
		return self.molarExtinction() / myMW if myMW else 0.0

	def molecularWeight(self):
		'''
		Add the individual molecular weights of the aas and find the
		mass of the protein keeping the loss of water in account.

		Return:
		    mw: molecular weight of the protein.
		'''
		mw = ProteinParam.mwH2O
		for acid in ProteinParam.aaComp:
			mw += (ProteinParam.aaComp[acid] * (ProteinParam.aa2mw[acid] - ProteinParam.mwH2O))
		return mw

