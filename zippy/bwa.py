import os.path
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
        cmd = "%s mem" % self.bwa_exec \
              + " -t %i" % (self.cores*2) \
              + self.args \
              + " %s %s" % (self.genome_fa, fastq) \
              + " | %s view -b -o %s -" % (self.samtools_exec, out_bam)
        self.flowLog(cmd)
        self.addTask(label="bwamem", command=cmd, nCores=self.cores,
                     memMb=self.mem, dependencies="make_out_dir")


        # Sort BAM
        out_sorted_bam = os.path.join(self.output_dir, "out.sorted.bam")
        out_temp = os.path.join(self.output_dir, "tmp")
        cmd = self.samtools_exec \
              + " sort %s" % out_bam \
              + " -O bam" \
              + " -o " + out_sorted_bam \
              + " -T " + out_temp \
              + " -@ %i" % self.cores
        self.addTask(label="sort_bam", command=cmd, nCores=self.cores, memMb=min(1024 * 32, self.mem), dependencies="bwamem")

        # Clean up the unsorted BAM
        cmd = "rm {}".format(out_bam)
        self.addTask(label="del_unsorted_bam", command=cmd, dependencies="sort_bam")


if __name__ == "__main__":
    pass