#!/usr/bin/env python
"""
A Python script to reverse complement a genomic sequence, find a 1Mb gap,
and extract the distal bit before the 1Mb gap.
"""
import sys
import argparse
import pysam

def init_argparse() -> argparse.ArgumentParser:
	"""
	A function to parse the arguments.
	"""
	parser = argparse.ArgumentParser(
		usage="%(prog)s [OPTION] [FILE]...",
		description="Check for a 1Mb gap representing the rDNA array"
	)
	parser.add_argument(
		"-v", "--version", action="version",
		version = f"{parser.prog} version 1.0.0"
	)
	parser.add_argument('-f', '--fasta', type=str, required=True, help='fasta file for assembly reference')
	parser.add_argument('-a', '--haplotype', type=str, required=True, help='name of haplotype in fasta file')
	#parser.add_argument('-o', '--output', type=str, required=True, help='output filename and path')
	return parser

def find_gap(args) -> None:
	"""
	A function to find the gap of a haplotype in a fasta file.
	"""
	fastaFile = args.fasta
	fastaObject = pysam.FastaFile(fastaFile)
	#extract just the given haplotype from the fasta file
	#region = fastaObject.fetch(region=args.haplotype)
	region=fastaObject.fetch(fastaObject.references[0])

	currentGap = False
	currentGapStart = 0
	firstGapStart = 0
	lastGapEnd = 0
	#loop through each char in the string, looking for N to start or continue a gap
	for i in range(len(region)):
		char = region[i]
		if(currentGap==False and char=='N'):
			#print("gap starts at position", i)
			currentGap = True
			currentGapStart = i
			if(firstGapStart == 0):
				firstGapStart = i
		#any character other than N ends a gap
		if(currentGap==True and char!='N'):
			#print("gap ends before position", i)
			currentGap = False
			gapDistance = i - currentGapStart
			lastGapEnd = i
			#print(gapDistance)
			#check for 1Mb gap
			#if(gapDistance==1000000 or gapDistance==1000001):
				#print("end of distal bit is", currentGapStart)
				#with open(args.output, "w") as file1:
				#	file1.write(region[0:currentGapStart] + '\n')
				#return currentGapStart
	#print("Error: did not find 1Mb gap for " + args.haplotype)
	print(str(firstGapStart) + "," + str(lastGapEnd))

def main() -> None:
	"""
	A main function to call the other functions to parse arguments and find
	distal bits of a genomic sequence.
	"""
	parser = init_argparse()
	args = parser.parse_args()

	find_gap(args)

if __name__ == "__main__":
	main()
