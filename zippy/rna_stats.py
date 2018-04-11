"""
Combines star_stats, rseqc insert size and rrna count into a single easy to read file.
"""
import subprocess
import numpy
import csv
import os
import pysam

def star_keep_line(line):
    if 'Started job on' in line:
        return False
    if 'Started mapping on' in line:
        return False
    if 'Finished on' in line:
        return False
    if 'Mapping speed' in line:
        return False
    if '|' in line:
        return True

def parse_starlog(fname):
    results = []
    number_reads = -1
    unique_reads = -1
    with open(fname, 'r') as f:
        for line in f:
            line = line.strip()
            if star_keep_line(line):
                split_line = line.split('|')
                split_line = [x.strip() for x in split_line]
                if 'Number of input reads' in line:
                    number_reads = int(split_line[1])
                if 'Uniquely mapped reads number' in line:
                    unique_reads = int(split_line[1])
                results.append("\t".join(split_line)) #we swap the | divider with : 
    return (results, number_reads, unique_reads)

def parse_insert_size(fname):
    data = []
    results = []
    with open(fname, 'rb') as f:
        c = csv.reader(f, delimiter='\t')
        for line in c:
            if 'sameTranscript=No' in line or 'sameChrom=No' in line or 'dist=genomic' in line: #we remove non-traditional transcripts
                continue
            if int(line[1]) > 1000: #we remove the extreme large insert sizes as not 'real'
                continue
            data.append(int(line[1]))
    results.append('Average insert size\t{}'.format(numpy.mean(data)))
    results.append('Standard deviation insert size\t{}'.format(numpy.std(data)))
    return results

def get_intron_read_count(params):
    subprocess.call("bedtools intersect -a {bam_path} -b {intron_bed} -u -wa -f 0.1 -sorted > {temp_path}.temp.bam".format(bam_path=params.bam_path, intron_bed=params.intron_bed, temp_path=params.temp_path), shell=True)
    read_names = set()
    samfile = pysam.AlignmentFile("{}.temp.bam".format(params.temp_path), "rb")
    for line in samfile:
        read_names.add(line.query_name)
    samfile.close()
    #count = subprocess.check_output("samtools view -c {temp_path}.temp".format(temp_path=params.temp_path), shell=True)
    #subprocess.call("rm {temp_path}.temp".format(temp_path=params.temp_path), shell=True)
    return len(read_names)

def get_ribosome_read_count(params):
    ''' deprecated '''
    subprocess.call("samtools view {bam_path} -b -L {ribosome_bed} -o {temp_path}.bam".format(bam_path=params.bam_path, ribosome_bed=params.ribosome_bed, temp_path=params.temp_path), shell=True)
    read_names = set()
    samfile = pysam.AlignmentFile("{}.bam".format(params.temp_path), "rb")
    for line in samfile:
        read_names.add(line.query_name)
    samfile.close()
    print "{temp_path}.bam".format(temp_path=params.temp_path)
    #subprocess.call("rm {temp_path}.bam".format(temp_path=params.temp_path), shell=True)
    return len(read_names)

def get_unique_read_count(bam_path, bed_path, temp_bam_file):
    exit_code = subprocess.call("samtools view {bam_path} -b -L {bed_path} -o {temp_bam_file}".format(bam_path=bam_path, bed_path=bed_path, temp_bam_file=temp_bam_file), shell=True)
    if exit_code != 0:
        raise RuntimeError('Samtools exited with nonzero exit code.')
    read_names = set()
    samfile = pysam.AlignmentFile(temp_bam_file, "rb")
    for line in samfile:
        read_names.add(line.query_name)
    samfile.close()
    subprocess.call("rm {}".format(temp_bam_file), shell=True)
    return len(read_names)

def get_dup_count(bam_path):
    dup_count = subprocess.check_output("samtools view -c -f 1024 -q 255 {}".format(bam_path), shell=True)
    return int(dup_count.strip())/2 #per paired-read

def analyze_rna(params):
    results = []
    sample_name = os.path.basename(params.output_path)
    results.append('Sample\t{}'.format(sample_name))
    number_reads = -1
    if params.starlog_path is not None:
        (star_results, number_reads, unique_reads) = parse_starlog(params.starlog_path)
        results.extend(star_results)
    if params.ribosome_bed is not None:
        #ribosome_read_count = get_ribosome_read_count(params)
        ribosome_read_count = get_unique_read_count(params.bam_path, params.ribosome_bed, params.temp_path+'.ribosome.bam')
        results.append('Ribosomal reads\t{}'.format(ribosome_read_count))
        if number_reads > 0:
            results.append('Ribosomal read percentage\t{:.2%}'.format(float(ribosome_read_count)/float(number_reads))) 
        print results
    #if params.intron_bed is not None:
    #    intron_read_count = get_intron_read_count(params)
    #    results.append('Intron-overlapping reads\t{}'.format(intron_read_count))
    #    if number_reads > 0:
    #        results.append('Intron-overlapping read percentage\t{:.2%}'.format(float(intron_read_count)/float(number_reads))) #2*x because the read count is paired end reads, but our ribosomal count is not.
    #    print results        
    if params.manifest_bed is not None:
        manifest_read_count = get_unique_read_count(params.bam_path, params.manifest_bed, params.temp_path+'.manifest.bam')        
        #manifest_read_count = subprocess.check_output("samtools view {bam_path} -c -L {manifest_bed}".format(bam_path=params.bam_path, manifest_bed=params.manifest_bed), shell=True)
        results.append('Manifest reads\t{}'.format(manifest_read_count))
        if number_reads > 0:
            results.append('Manifest read percentage\t{:.2%}'.format(float(manifest_read_count)/float(number_reads))) #2*x because the read count is paired end reads, but our ribosomal count is not.
        print results        
    if params.rseqc_inner_distance is not None:
        results.extend(parse_insert_size(params.rseqc_inner_distance))
    if params.dup_stats:
        dup_count = get_dup_count(params.bam_path)
        results.append('Duplicate uniquely mapping reads\t{}'.format(dup_count))
        if unique_reads > 0:
            results.append('Duplicate read percentage (of all uniquely mapped reads)\t{:.2%}'.format(float(dup_count)/float(number_reads))) #2*x because the read count is paired end reads, but our ribosomal count is not.
    with open(params.output_path, 'w') as fout:
        print results
        fout.write('\n'.join(results))

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(usage='', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('bam_path', help='The bam to analyze.')
    parser.add_argument('output_path', help='Where to place the stats file.')
    parser.add_argument('--ribosome_bed')
    parser.add_argument('--intron_bed')
    parser.add_argument('--manifest_bed')
    parser.add_argument('--temp_path')
    parser.add_argument('--rseqc_inner_distance')
    parser.add_argument('--starlog_path', help="Star's Final.log.out file.")
    parser.add_argument('--dup_stats', dest='dup_stats', action='store_true', help="Whether we generate duplicate statistics")
    parser.set_defaults(dup_stats=False)
    input_params = parser.parse_args()
    analyze_rna(input_params)