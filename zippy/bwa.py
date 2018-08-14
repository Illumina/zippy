import os.path
import math
from pyflow import WorkflowRunner

class BWAWorkflow(WorkflowRunner):
    def __init__(self, output_dir,
                 bwa_exec, samtools_exec, genome_fa, cores, mem,
                 fastq, sample='', args=''):
        self.output_dir = output_dir
        self.bwa_exec = bwa_exec
        self.samtools_exec = samtools_exec
        self.genome_fa = genome_fa
        self.cores = cores
        self.mem = mem
        self.fastq = fastq
        self.sample = sample
        self.args = args

    def workflow(self):
        # Create the output directory and create a dummy samplesheet
        cmd = "mkdir -p {}".format(self.output_dir)
        self.addTask(label="make_out_dir", command=cmd, isForceLocal=True)
        out_bam = os.path.join(self.output_dir, "out.bam")

    
        # from fastq:
        if len(self.fastq) == 2:
            fastq = " ".join(self.fastq)
        elif len(self.fastq) == 1:
            fastq = self.fastq
        else:
            raise("More than two FASTQs passed to bwamem!")
        if self.args != "":
            self.args=" "+self.args
        else:
            self.args = " -M -R \'@RG\\tID:1\\tLB:{0}\\tPL:ILLUMINA\\tSM:{0}\'".format(self.sample)
        # Figure out the number of cores we can use for alignment and for bam compression
        total_threads = self.cores * 2  # At least 2
        addtl_compression_threads = max(int(0.1 * total_threads), 1) # At a minimum, allocate one extra thread for bam compression
        bwa_threads = total_threads  # Since BWA's output is thread-dependent, we don't decrement here in order to avoid surprises
        assert bwa_threads >= 1
        cmd = "%s mem" % self.bwa_exec \
              + " -t %i" % (bwa_threads) \
              + self.args \
              + " %s %s" % (self.genome_fa, fastq) \
              + " | %s view -@ %i -1 -o %s -" % (self.samtools_exec, addtl_compression_threads, out_bam)
        self.flowLog(cmd)
        self.addTask(label="bwamem", command=cmd, nCores=self.cores,
                     memMb=self.mem, dependencies="make_out_dir")


        # Sort BAM
        out_sorted_bam = os.path.join(self.output_dir, "out.sorted.bam")
        out_temp = os.path.join(self.output_dir, "tmp")
        # Calculate resources for sort
        sort_threads = self.cores * 2
        mem_per_thread = int(math.floor(float(self.mem) / sort_threads * 0.75))  # Per thread, underallocate to allow some overhead
        cmd = self.samtools_exec \
              + " sort %s" % out_bam \
              + " -O bam" \
              + " -o " + out_sorted_bam \
              + " -T " + out_temp \
              + " -@ %i" % sort_threads
        self.addTask(label="sort_bam", command=cmd, nCores=self.cores, memMb=self.mem, dependencies="bwamem")

        # Clean up the unsorted BAM
        cmd = "rm {}".format(out_bam)
        self.addTask(label="del_unsorted_bam", command=cmd, dependencies="sort_bam")


if __name__ == "__main__":
    pass