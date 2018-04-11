import os
import gzip
import csv
from collections import defaultdict
from itertools import combinations

class PrimerDimerMiner():
    """
    A primer dimer detector.  Looks for kmers that can be uniquely mapped to a primer or its reverse complement,
    and reports all reads that have kmers belonging to at least 2 such primers.  Note that a primer and its reverse
    complement are counted separately.  We also report autodimers, defined as non-contiguous kmers that map to the same
    primer.
    #TODO: we ignore whether reads are paired end; not sure if there's anything smart we want to do there.
    #TODO: all matches must be exact.  if there's a read error somewhere, we won't catch it, or we might call it an autodimer
    # if it results in 15bp matches on each end of the error.
    """
    def __init__(self, params):
        self.params = params
        self.k = self.params.kmer
        self.reads_read = 0 # how many lines we've looked at
        self.dimer_reads_read = 0 # how many lines we've looked at with at least 2 primers
        self.monomer_reads_read = 0 # how many lines we've looked at with exactly 1 primer
        self.pairs_map = defaultdict(int) # map from tuple of 2 primers (alphabetical) to dimer count
        self.out_reads_file = open(os.path.join(self.params.output_dir,'matched_reads'), 'w')
        self.out_pairs_file = open(os.path.join(self.params.output_dir,'pair_stats'), 'w')
        self.out_stats_file = open(os.path.join(self.params.output_dir,'run_stats'), 'w')

    def reverse_complement(self, seq):
        rev_map = {'A':'T','C':'G','G':'C','T':'A'}
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

    def run_matcher_helper(self, f, kmer_map):
        counter = 0 # 0: header 1: sequence 2: junk 3: quality
        for line in f:
            if counter == 0:
                read_name = line.strip()
            if counter == 1:
                line = line.upper()
                probe_matches = self.find_matches(line.strip(), kmer_map)
                if len(probe_matches) > 1:
                    self.output_matches(read_name, probe_matches, line.strip())
                    self.register_pairs(probe_matches)
                    self.dimer_reads_read += 1
                elif len(probe_matches) == 1:
                    self.monomer_reads_read += 1
                self.reads_read += 1
            counter += 1
            counter = counter % 4

    def run_matcher(self, input_file, kmer_map):
        """
        Goes through a fastq, and finds all dimerized reads
        """
        if input_file is None: # we don't require the input files... one of them can be undefined
            return
        try:
            with gzip.open(input_file, 'r') as f:
                self.run_matcher_helper(f, kmer_map)
        except IOError: #not a gzip
            with open(input_file, 'r') as f:
                self.run_matcher_helper(f, kmer_map)


    def find_matches(self, line, kmer_map):
        """
        For a single read, reports all found primers
        """
        in_primer = None
        matches = []
        for i in range(0, len(line)-self.k):
            sub_seq = line[i:i+self.k]
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

    def output_matches(self, read_name, probe_matches, line):
        self.out_reads_file.write(read_name+"\n")
        self.out_reads_file.write(line+"\n")
        self.out_reads_file.write("\t".join(probe_matches)+"\n")

    def register_pairs(self, probe_matches):
        """
        Increments the occurrence count of all unique dimer pairs found in a line.
        """
        seen_pairs = set()
        probe_matches = sorted(probe_matches)
        for pair in combinations(probe_matches, 2):
            if pair not in seen_pairs:
                self.pairs_map[pair] += 1
                seen_pairs.add(pair)

    def output_stats(self):
        self.out_pairs_file.write("\t".join(["Primer 1", "Primer 2", "# reads containing both", "% reads containing both"])+"\n")
        for (pair, count) in sorted(self.pairs_map.items(), key=lambda x: -x[1]):
            self.out_pairs_file.write("\t".join([pair[0], pair[1], str(count), "{:.4%}".format(float(count)/float(self.dimer_reads_read))])+"\n")
        self.out_stats_file.write("Number of reads seen: {}".format(self.reads_read)+"\n")
        self.out_stats_file.write("Number of reads with dimers: {}, ({:.4%})\n".format(self.dimer_reads_read, float(self.dimer_reads_read)/float(self.reads_read)))
        self.out_stats_file.write("Number of reads with exactly 1 primer: {}, ({:.4%})\n".format(self.monomer_reads_read, float(self.monomer_reads_read)/float(self.reads_read)))




    def run(self):
        """
        Main execution function for PDM
        """
        probe_map = self.parse_probe_list()
        kmer_map = self.build_kmer_map(probe_map, self.k)
        for input_file in [self.params.input_file1, self.params.input_file2]:
            self.run_matcher(input_file, kmer_map)
        self.output_stats()
        self.out_reads_file.close()
        self.out_stats_file.close()
        self.out_pairs_file.close()


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--probe_list', type=str, help="tsv with 2 columns: (probe_name, sequence).  If you are looking for adapters or other short sequences, they should be added to the probe list.")
    parser.add_argument('--kmer', type=int, default=15, help="How big a fragment size to require for a match")
    parser.add_argument('--input_file1', help="A fastq file with reads to analyze")
    parser.add_argument('--input_file2', help="Another fastq file (Optional)")
    parser.add_argument('--output_dir')
    params = parser.parse_args()
    if not os.path.exists(params.output_dir):
        os.makedirs(params.output_dir)
    pdm = PrimerDimerMiner(params)
    pdm.run()    