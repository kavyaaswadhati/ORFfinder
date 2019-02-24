#!/usr/bin/env python3
# Name: Kavya Aswadhati (kaswadha)
# Group Members: None
############################################################################################################
#   the sequence will be traversed one nucleotide at a time, the position along the sequence will belie
#   the reading frame we are operating within.
#
#   as we traverse the sequence we will check for start codons and populate a dictionary of lists
#   { rf: [{start:num, stop:num},...]
#   as we encounter stops we will populate with consideration for position in sequence, potential preceeding
#   starts/stops
#
#   this process will be repeated with the complementary strand
#
#   once this data structure is populated, we will find putative genes by comparing each stop with its associated
#   start.
#
#
############################################################################################################
from sequenceAnalysis import FastAReader
from Bio.Seq import Seq
import sys
from operator import attrgetter


class ORFFinder:
	'''
    find Open Reading Frames in a given sequence
    input: sequence
    output: file of putative genes sorted by length
    '''

	def __init__(self, args, inHeader="", inSequence=""):
		'''
        Initialize ORF finder by saving user input parameters and saving the sequence to process.
        Call the methods within ORFFinder that will do the task of parsing the sequence for ORFs
        '''
		self.start = args.start
		self.stop = args.stop

		self.minGene = args.minGene
		# creating boolean longestGene from string longestGene given by input
		if args.longestGene == True:
			self.longestGene = True
		else:
			self.longestGene = False
		self.seqName = inHeader
		self.seqLength = len(inSequence)
		self.numGenes = 0

		# bool to keep track of which strand we are parsing
		self.isComplement = False

		# negative indicates complementary strand
		self.foundCodons = {
			1: [],
			2: [],
			3: [],
			-1: [],
			-2: [],
			-3: []
		}

		# list of dictionaries that contains genes found
		# each gene will be defined by:
		# {[reading frame]: start position, stop position, length}
		self.foundGenes = []
		self.findGenes(inSequence)

		# once strand has been processed, make sure it is known that this is now the complement strand
		# and re-run
		self.isComplement = True
		reverseComp = self.reverse(inSequence)
		self.findGenes(reverseComp)

		# convert foundCodons to foundGenes
		self.putativeGenes()

		# output
		self.printORFs()

	def reverse(self, s):
		'''Generate the reverse complement of the orginal sequence'''
		r = ''
		for i in s:
			r = i + r
		return r

	def findGenes(self, sequence):
		'''
        Takes input of sequence and populates foundCodons by calling
        addFoundCodon(index) and addStop(index)
        '''
		# entries will be ie
		# 1:[{'starts':[0,3,6]
		#      'stop':9},
		#     {'starts':[12,15,18]
		#       'stop':21}, ... ]
		# where 1 as the key indicates that we are in reading frame 1, and the value associated with this
		# is a list of dictionaries defined by one stop codon and all its associated start codons

		index = 0
		while index < self.seqLength - 3:
			if sequence[index: index + 3] in self.start:
				# add a dict of {reading frame (1,2,or 3), position, start codon}
				self.addFoundCodons(index + 1)
			elif sequence[index: index + 3] in self.stop:
				self.enterStop(index + 1)
			index += 1

	def addFoundCodons(self, position):
		'''
        Takes input of position and adds to appropriate part of foundCodons(position).
        Keeps in account reading frame and strand complementarity
        '''

		# establish which reading frame we are operating in (- for  strand)
		if self.isComplement:
			readingFrame = -(((position - 1) % 3) + 1)
		else:
			readingFrame = (((position - 1) % 3) + 1)

		item = self.foundCodons[readingFrame]
		# print ('position',position)
		# print('item', item)
		# case : not empty list for reading frame
		if len(item):
			# checking the last dictionary keyed to the reading frame
			assocCodons = self.foundCodons[readingFrame][-1]
			if 'stop' in assocCodons:
				# if it already has a stop codon, we need to start a new dictonary
				newAssoc = {'starts': [position]}
				self.foundCodons[readingFrame].append(newAssoc)

			else:
				# if there is no stop codon, append position to last start codon list
				assocCodons['starts'].append(position)
		# preliminary case, this is the first start encountered in reading frame
		else:
			firstAssoc = {'starts': [position]}
			self.foundCodons[readingFrame].append(firstAssoc)


	def enterStop(self, position):
		'''
        Takes input of position to place stop codon in appropriate positon in
        foundCodons. Takes into consideration previous starts, reading frame, strand complementarity
        '''
		# ensure that new codon is added to negative rf if from complementary strand
		if self.isComplement:
			readingFrame = -(((position - 1) % 3) + 1)
		else:
			readingFrame = ((position - 1) % 3) + 1
		item = self.foundCodons[readingFrame]

		# case 1: there are start codons in the list for associated reading frame
		if len(item):
			# assocCodons is the dictionary of starts for item
			assocCodons = self.foundCodons[readingFrame][-1]
			# print('item',item)
			if 'stop' in assocCodons:
				newAssoc = {'stop': position + 3}
				self.foundCodons[readingFrame].append(newAssoc)
			else:
				self.foundCodons[readingFrame][-1]['stop'] = position
		# case 2: edge case, stop with no starts preceding it
		else:
			firstAssoc = {'stop': position + 3}
			# print('first Assoc', firstAssoc)
			self.foundCodons[readingFrame].append(firstAssoc)
		return

	def putativeGenes(self):
		'''
        Parses through foundCodons to extract valid putative genes and populates them in foundGenes
        '''
		for rf, readingFrames in self.foundCodons.items():
			# accounting for which strand we are in (will have to adjust position saving)
			complement = False
			if rf < 0: complement = True

			for index, assoc in enumerate(readingFrames):
				# front end fragment case: first codon is a stop
				if index == 0 and 'starts' not in assoc:
					if (assoc['stop']+2)>= self.minGene:
						newGene = {'rf': rf,
								   'startPos': 1,
								   'stopPos': assoc['stop']+2,
								   'length': assoc['stop']+2
								   }
						self.numGenes += 1
						self.foundGenes.append(newGene)
						continue
				# for internal stops with no associated starts
				if 'starts' not in assoc:
					continue
				if 'stop' not in assoc:
					stopPos = self.seqLength
				else:
					# typical case, start and stops exist, comparing stop w/ first start
					stopPos = assoc['stop']+2

				# if not only looking for the longest gene per stop, save a length comparison for
				# stop and all its associated starts, check if they satisfy the min length, and save them
				if not self.longestGene:
					# print('lg false')
					lens = []
					for start in assoc['starts']:
						lens.append([start, (stopPos - start)])
					for len in lens:
						print ('len =', len)
						if len[1] > self.minGene:
							newGene = {'rf': rf,
									   'startPos': len[0],
									   'stopPos': stopPos,
									   'length': len[1]
									   }
							self.numGenes += 1
							self.foundGenes.append(newGene)
				# only saving the longest gene
				elif self.longestGene:
					# print ('lg true')
					firstStartPos = assoc['starts'][0]
					length = stopPos - firstStartPos
					if (length) > self.minGene:
						newGene = {'rf': rf,
								   'startPos': firstStartPos,
								   'stopPos': stopPos,
								   'length': length
								   }
						self.numGenes += 1
						self.foundGenes.append(newGene)
		return


	def printORFs(self):
		'''
        Takes information from foundGenes and formats/outputs it to the outFile
        '''
		# used when file is specified in cl
		# print(self.seqName, file=self.file)
		print(self.seqName)
		for gene in sorted(self.foundGenes, key=lambda gene: gene['length'], reverse=True):
			printStr = '{:+d} {:>5d}..{:<5d} {:<5d}'.format(
				gene['rf'],
				gene['startPos'],
				gene['stopPos'],
				gene['length'])
			print(printStr)

			# cl file name print implementation
			# print(printStr, file=self.file)


class CommandLine():
	'''
    Class for handling user input in the command line provided by Prof. Bernick
    '''

	def __init__(self, inOpts=None):
		import argparse
		self.parser = argparse.ArgumentParser(description='Program prolog')

		# removed because infile and outfile are now specified by stdin, stdout
		# self.parser.add_argument('inFile', action='store', help='input file name')
		# self.parser.add_argument('outFile', action='store', help='output file name')

		self.parser.add_argument('-lG', '--longestGene', action='store', nargs='?', const=True, default=False,
			help='longest Gene in an ORF')
		self.parser.add_argument('-mG', '--minGene', type=int, choices=(9, 100, 200, 300, 500, 1000), default=100,
			action='store', help='minimum Gene length')
		self.parser.add_argument('-s', '--start', action='append', default=['ATG'], nargs='?',
			help='start Codon')  # allows multiple list options
		self.parser.add_argument('-t', '--stop', action='append', default=[' TAG', 'TGA', 'TAA'], nargs='?',
			help='stop Codon')  # allows multiple list options
		self.parser.add_argument('-v', '--version', action='version', version='% (prog)s 0.1')
		if inOpts is None:
			self.args = self.parser.parse_args()
		else:
			self.args = self.parser.parse_args(inOpts)


def main(inCL=None):
	if inCL is None:
		myCMD = CommandLine()
	else:
		myCMD = CommandLine(inCL)
	print(myCMD.args)
	myReader = FastAReader()
	for head, seq in myReader.readFasta():
		newORFs = ORFFinder(myCMD.args, head, seq)
	# file.close()


if __name__ == "__main__":
	main()
	# main(['tass2.fa','tass2ORFdata-ATG-100.txt','--longestGene'])
