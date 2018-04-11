# -*- coding: utf-8 -*-
"""
Downsampling script for paired end data that retains information about non-linear alignments.
If the specified value for downsampling is less than the number of reads in the original bam, 
then it will simply copy the old bam to the new path while omitting any pairs that have no
alignment whatsoever. This makes it so that regardless if downsampling was actually performed,
the output is in the same "style" (i.e. includes no pairs that are complete unmapped).

Author: Kevin Wu, Jan. 2017
"""

import argparse
import os
import random

import pybloomfilter
import pysam
import time
import multiprocessing
import shutil


BF_PROB = 0.0001
MAX_ITERS = 10 # Maximum times we'll loop through a bam

def partition(lst, n):
    """Given a list, partition it into n chunks of approximately equal size"""
    # https://stackoverflow.com/questions/2659900/python-slicing-a-list-into-n-nearly-equal-length-partitions
    return [lst[i::n] for i in xrange(n)]


class allReadNamesProcess(multiprocessing.Process):
    """
    Outputs all read names from given chromosomes and writes to a centralized queue
    This is meant for the non-downsampling branch of this script.
    """
    def __init__(self, bamfile, chromosome_list, outqueue):
        multiprocessing.Process.__init__(self)
        self.out_queue = outqueue
        self.chromosomes = chromosome_list
        self.bam = bamfile

    def run(self):
        """
        Walks through each chromosome, outputting every single readname
        """
        for chromosome in self.chromosomes:
            with pysam.AlignmentFile(self.bam, 'rb') as alignment:
                for read in alignment.fetch(chromosome):
                    self.out_queue.put(read.query_name)


class sampleReadNamesProcess(multiprocessing.Process):
    """
    Sample read names from the given chromosomes and writes them to a centralized queue
    """
    def __init__(self, bamfile, perChromSingleEndReads, perChromUniqueReadnamesTarget, outQueue):
        multiprocessing.Process.__init__(self)
        self.bam = bamfile
        # The dictionary keys in the perChrom things should be the same
        assert isinstance(perChromSingleEndReads, dict)
        assert isinstance(perChromUniqueReadnamesTarget, dict)
        for key in perChromSingleEndReads.keys():
            assert key in perChromUniqueReadnamesTarget
        
        self.chromosomes = perChromSingleEndReads.keys()
        self.perChromSingleEndReads = perChromSingleEndReads
        self.perChromReadnamesTarget = perChromUniqueReadnamesTarget

        # assert isinstance(outQueue, multiprocessing.Queue)
        self.outQueue = outQueue

    def run(self):
        """
        Walk through this thread's assigned chromosomes, getting readnames for each.
        """
        for chromosome in self.chromosomes:
            # Sanity checks
            assert chromosome in self.perChromSingleEndReads
            assert chromosome in self.perChromReadnamesTarget

            # Initialize values for this chromosome's sampling
            targetReadnamesCount = self.perChromReadnamesTarget[chromosome]
            numSingleEndReads = self.perChromSingleEndReads[chromosome]
            acceptedProbability = float(targetReadnamesCount) / float(numSingleEndReads)
            acceptedReadnames = pybloomfilter.BloomFilter(targetReadnamesCount * 2, BF_PROB)

            numIters = 0
            keepGoing = True
            with pysam.AlignmentFile(self.bam, 'rb') as alignment:
                while keepGoing and numIters < MAX_ITERS:
                    for read in alignment.fetch(region=chromosome):
                        read_name, flag = read.query_name, int(read.flag)
                        # If read is supplementary/seconadary, skip it
                        if (flag & 2048) or (flag & 256):
                            continue
                        if random.random() <= acceptedProbability:
                            # False output from add indicates that item was not already in bloom filter
                            # So we tehn add it to the output queue
                            if not acceptedReadnames.add(read.query_name):
                                self.outQueue.put(read_name)
                                if len(acceptedReadnames) >= targetReadnamesCount:
                                    keepGoing = False
                                    break
                    numIters += 1
            # Write to stdout
            print("Processed %s %s: (readnames selected/readnames desired) (%i/%i)" % (os.path.basename(self.bam), chromosome, len(acceptedReadnames), targetReadnamesCount))


class readnameConsumer(multiprocessing.Process):
    """
    Consumes readnames and writes the output bam.
    """
    def __init__(self, totalTargetReadnameCount, bamfileIn, bamfileOut, readnameQueue):
        multiprocessing.Process.__init__(self)
        # assert isinstance(readnameQueue, multiprocessing.Queue)
        assert isinstance(totalTargetReadnameCount, int)
        self.bamIn = bamfileIn
        self.bamOut = bamfileOut
        self.queue = readnameQueue
        self.readnameCount = totalTargetReadnameCount
    
    def run(self):
        # Combine the readnames from each of the child processes into one single bloom filter
        masterBF = pybloomfilter.BloomFilter(self.readnameCount * 2, BF_PROB)
        counter = 0
        while True:
            readname = self.queue.get()
            if readname is None:
                break
            masterBF.add(readname)
            counter += 1

        print("Finished sampling readnames. Writing output bam...")
        # Finished getting all the reads, now write the output bam
        with pysam.AlignmentFile(self.bamIn, 'rb') as bamSource:
            with pysam.AlignmentFile(self.bamOut, 'wb', template=bamSource) as bamOutput:
                for read in bamSource.fetch():
                    if read.query_name in masterBF:
                        bamOutput.write(read)
        print("Finished downsampling %s to %s. Outputted (selected/desired) %i/%i readnames." % (self.bamIn, self.bamOut, counter, self.readnameCount))

def getTotalReads(bam):
    totalReads = 0
    perChromCount = {}
    stats = pysam.idxstats(bam)
    for line in stats.split('\n'):
        tokenized = line.split()
        if len(tokenized) == 0 or tokenized[0] == "*": continue
        c = int(tokenized[2]) + int(tokenized[3]) # mapped + unmapped reads
        perChromCount[tokenized[0]] = c
        totalReads += c
    return totalReads, perChromCount


def downsample(inputBam, outputBam, targetNumUniqueReadnames, numThreads=4):
    # Sanity checks
    if not os.path.isfile(inputBam):
        raise RuntimeError("Cannot find %s" % inputBam)
    indexFile = inputBam + ".bai"
    if not os.path.isfile(indexFile):
        raise RuntimeError("Input bam %s must be sorted" % inputBam)
    
    totalSingleEndReads, perChromSingleEndReads = getTotalReads(inputBam)
    assert totalSingleEndReads != -1

    # Skip downsampling if we already have fewer than downsampled amount
    if targetNumUniqueReadnames * 2 >= totalSingleEndReads:
        print("%s does not need downsampling. Copying old bam to new path while omitting completely unaligned pairs" % inputBam)
        # Copy the file over
        needs_downsampling = False
    else:
        print("Starting downsampling for %s" % inputBam)
        needs_downsampling = True

    # Get the chromosomes included in this bam
    with pysam.AlignmentFile(inputBam, 'rb') as bIn:
        chromosomes = bIn.references

    # Split the chromosomes into per-thread lists of chromosomes
    chromosomeChunks = partition(chromosomes, numThreads)
    readnamesQueue = multiprocessing.Queue()
    # Initialize the readnames consumer
    consumerProcess = readnameConsumer(targetNumUniqueReadnames, inputBam, outputBam, readnamesQueue)
    consumerProcess.start()
    producerProcesses = []
    for chunk in chromosomeChunks:
        if needs_downsampling:
            perThreadSingleEndReads = {}
            for chrom in chunk:
                perThreadSingleEndReads[chrom] = perChromSingleEndReads[chrom]
            perThreadTargetReadnames = {}
            for chrom in chunk:
                perThreadTargetReadnames[chrom] = int(float(perChromSingleEndReads[chrom]) / float(totalSingleEndReads) * float(targetNumUniqueReadnames))
            chunk_process = sampleReadNamesProcess(inputBam, perThreadSingleEndReads, perThreadTargetReadnames, readnamesQueue)
        else:
            chunk_process = allReadNamesProcess(inputBam, chunk, readnamesQueue)

        chunk_process.start()
        producerProcesses.append(chunk_process)

    # Wait for all the per-chrom threads to terminate
    for proc in producerProcesses:
        proc.join()

    # Signal to the writer that we're done
    readnamesQueue.put(None)

    # Wait for the writer to finish
    consumerProcess.join()


def buildParser():
    parser = argparse.ArgumentParser(
        description = __doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--bam", type=str, required=True,
                        help="Bam file to downsample")
    parser.add_argument("--downsampled_pairs", type=int, required=True,
                        help="Number of unique pairs to downsample to")
    parser.add_argument("--output", type=str, required=True,
                        help="Path to output bam file")
    parser.add_argument("--threads", type=int, default=4,
                        help="Number of threads to use")
    return parser


if __name__ == "__main__":
    # main()
    parser = buildParser()
    args = parser.parse_args()

    downsample(args.bam, args.output, args.downsampled_pairs, args.threads)
