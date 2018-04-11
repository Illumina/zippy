import cProfile
import os
import gzip
import csv
import time
from collections import defaultdict
from itertools import combinations
from pybloomfilter import BloomFilter
#from pybloom import ScalableBloomFilter, BloomFilter #pybloom used cryptographic hashes in a bloom filter.  This is a bad idea.
import numpy

class BloomKmerFinder():
    """
    Finds all kmers that show up more than a certain number of times.  Can choose to ignore dimerized reads
    or do only dimerized reads.  Useful for finding common kmers in unmapped reads.  We use a bloom filter
    to do this, so it is very fast, but requires a few GB of ram to keep the filter in memory.
    """
    def __init__(self, params, k, exclude_monomers=False, exclude_dimers=False):
        self.params = params
        self.k = k #kmer for global matching
        self.primer_k = 15  #kmer for primer matching
        self.exclude_monomers = exclude_monomers
        self.exclude_dimers = exclude_dimers
        self.run_dimer_detector = False
        self.bloom = BloomFilter(1e9, 0.01, None) 
        self.count_map = defaultdict(int)
        if exclude_monomers or exclude_dimers:
            self.run_dimer_detector = True
        self.reads_read = 0 # how many lines we've looked at
        self.dimer_reads_read = 0 # how many lines we've looked at with at least 2 primers
        self.monomer_reads_read = 0 # how many lines we've looked at with exactly 1 primer
        self.out_stats_file = open(os.path.join(self.params.output_dir,'count_stats'), 'w')

    def reverse_complement(self, seq):
        rev_map = {'A':'T','C':'G','G':'C','T':'A', 'N':'N'}
        new_seq = ''.join([rev_map[x] for x in seq]) #sub letters
        new_seq = new_seq[::-1] # reverse it
        return new_seq

    def parse_probe_list(self):
        """
        Creates probe map, which is a map of probe names to sequence.
        """
        probe_map = {}
        with open(self.params.probe_list, 'r') as f:
            c = csv.reader(f, delimiter="\t")
            for line in c:
                probe_map[line[0]] = line[1].upper()
        return probe_map

    def build_kmer_map(self, probe_map, k):
        """
        Builds a map from kmer to probenames that have this kmer.
        Also does reverse complements.
        """
        kmer_map = defaultdict(set)
        for (probe, seq) in probe_map.iteritems():
            seq_rc = self.reverse_complement(seq)
            for i in range(0, len(seq)-k):
                kmer_map[seq[i:i+k]].add(probe)
                kmer_map[seq_rc[i:i+k]].add(probe+"rc")
        return kmer_map

    def run_matcher(self, input_file, kmer_map):
        """
        Goes through a fastq, and registers all the kmers in it.
        """
        if input_file is None: # we don't require the input files... one of them can be undefined
            return
        debug = 0
        with open(input_file, 'r') as f:
            counter = 0 # 0: header 1: sequence 2: junk 3: quality
            for line in f:
                if counter == 0:
                    read_name = line.strip()
                if counter == 1:
                    line = line.upper()
                    if self.run_dimer_detector:
                        probe_matches = self.find_matches(line.strip(), kmer_map)
                        if len(probe_matches) > 1 and not self.exclude_dimers:
                            self.dimer_reads_read += 1
                            self.register_kmers(line.strip())
                        elif len(probe_matches) == 1 and not self.exclude_monomers:
                            self.monomer_reads_read += 1
                            self.register_kmers(line.strip())
                        elif len(probe_matches) == 0:
                            debug += 1
                            self.register_kmers(line.strip())                            
                    else:
                        self.register_kmers(line.strip())
                    self.reads_read += 1
                counter += 1
                counter = counter % 4
        print '{} dimer: {}'.format(input_file, self.dimer_reads_read)
        print '{} monomer: {}'.format(input_file, self.monomer_reads_read)
        print '{} none: {}'.format(input_file, debug)
        print '{} total: {}'.format(input_file, self.reads_read)

    def register_kmers(self, read):
        """
        Adds the read and its reverse complement to our bloom filter, and if we have seen it before,
        adds it to the count map.  The idea is that the bloom filter can approximately determine
        if we've seen something before or not, and to the count map are added all kmers that the bloom
        filter reports that we've seen before.
        """
        for i in range(0, len(read)-self.k):
            seq = read[i:i+self.k]
            seq_rc = self.reverse_complement(seq)
            if self.bloom.add(seq):
                self.count_map[seq]+=1
            if self.bloom.add(seq_rc):
                self.count_map[seq_rc]+=1


    def find_matches(self, line, kmer_map):
        """
        For a single read, reports all found primers
        """
        in_primer = None
        matches = []
        for i in range(0, len(line)-self.primer_k):
            sub_seq = line[i:i+self.primer_k]
            if in_primer is None: #we are not currently in a primer.    
                if len(kmer_map[sub_seq]) == 1: # If we see a uniquely mappable kmer, we enter a primer.
                    (in_primer,) = kmer_map[sub_seq]
                    matches.append(in_primer)
                else: # Otherwise, we continue
                    continue
            else: # we are in the middle of seeing a primer sequence
                if in_primer in kmer_map[sub_seq]: # we see this primer again, and are thus still reading it.  Continue.
                    continue
                elif len(kmer_map[sub_seq]) == 1: # We no longer see our current primer, but this sequence is mappable to another primer.  We are now in a different primer.
                    (in_primer,) = kmer_map[sub_seq]
                    matches.append(in_primer)
                else: # We aren't in our current primer, and aren't uniquely in a different primer.
                    in_primer = None
        return matches


    def output_stats(self, kmer_map):
        """
        We print the top-two unique maximal strings, and then a sorted list of all kmers that appear
        at least twice in our reads.
        """
        sorted_map = sorted(self.count_map.items(), key=lambda x: -x[1])
        first_string = self.extend(sorted_map[0][0], self.count_map)
        for (kmer, count) in sorted_map[1:]:
            if kmer not in first_string and self.reverse_complement(kmer) not in first_string:
                second_string = self.extend(kmer, self.count_map)
                second_score = count
                break
        self.out_stats_file.write("{}\t{}\n".format(sorted_map[0][1], first_string))
        self.out_stats_file.write("{}\t{}\n".format(second_score, second_string))
        for (kmer, count) in sorted_map:
            probe_matches = self.find_matches(kmer, kmer_map)            
            if len(probe_matches) == 0:
                self.out_stats_file.write("{}\t{}\n".format(kmer, count))
            else:
                self.out_stats_file.write("{}\t{}\t{}\n".format(probe_matches, kmer, count))

    def extend(self, seed, kmer_map):
        """
        Given a kmer, we greedily extend it in both directions by looking for kmers that differ by 1 on either side.  We add
        the new kmer if its count is at least half of our peak kmer.
        """
        final_string = [seed]
        value = kmer_map[seed]
        forward_extend = True
        current_seed = seed
        while forward_extend:
            extender = current_seed[1:]
            new_kmers = [extender+x for x in 'ACGT']
            new_scores = [kmer_map[x] for x in new_kmers]
            if numpy.max(new_scores)>value*0.5: #we extend
                new_kmer = new_kmers[numpy.argmax(new_scores)]
                if new_kmer == current_seed: #we hit a pathological (recursive) read
                    forward_extend = False
                final_string.append(new_kmer[-1])
                current_seed = new_kmer
            else:
                forward_extend = False
        reverse_extend = True
        current_seed = seed
        while reverse_extend:
            extender = current_seed[:-1]
            new_kmers = [x+extender for x in 'ACGT']
            new_scores = [kmer_map[x] for x in new_kmers]
            if numpy.max(new_scores)>value*0.5: #we extend
                new_kmer = new_kmers[numpy.argmax(new_scores)]
                if new_kmer == current_seed: #we hit a pathological read
                    reverse_extend = False
                final_string = [new_kmer[0]]+final_string
                current_seed = new_kmer
            else:
                reverse_extend = False
        return ''.join(final_string)


    def run(self):
        """
        Main execution function for kmer finder
        """
        if self.run_dimer_detector:
            probe_map = self.parse_probe_list()
            kmer_map = self.build_kmer_map(probe_map, self.primer_k) #kmer-map for dimer detection
        else:
            kmer_map = defaultdict(set)
        for input_file in [self.params.input_file1, self.params.input_file2]:
            self.run_matcher(input_file, kmer_map)
        self.output_stats(kmer_map)
        self.out_stats_file.close()

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--probe_list', type=str, help="Needed if you want to filter by dimers.  tsv with 2 columns: (probe_name, sequence).  If you are looking for adapters or other short sequences, they should be added to the probe list.")
    parser.add_argument('--kmer', type=int, default=30, help="How big a fragment size to count")
    parser.add_argument('--input_file1', help="A fastq file with reads to analyze")
    parser.add_argument('--input_file2', help="Another fastq file (Optional)")
    parser.add_argument('--output_dir')
    parser.add_argument('--exclude_monomers', dest='exclude_monomers', action='store_true', help="Whether we exclude primer monomers from kmer counting")
    parser.set_defaults(exclude_monomers=False)
    parser.add_argument('--exclude_dimers', dest='exclude_dimers', action='store_true', help="Whether we exclude primer dimers from kmer counting")
    parser.set_defaults(exclude_dimers=False)
    params = parser.parse_args()
    bloomy = BloomKmerFinder(params, params.kmer, params.exclude_monomers, params.exclude_dimers)
    start = time.time()
    bloomy.run()
    print time.time()-start