from pyflow import WorkflowRunner

class Bcl2Fastq(WorkflowRunner):
    def __init__(self, bcl2fastq, runfolder_dir, output_dir, sample_sheet, args='', max_job_cores=16):
        optional_argument_defaults = {  #'demultiplexing-threads': 10,
                                        'processing-threads': 16}
        self.max_job_cores = max_job_cores
        command = "{} --runfolder-dir {} --output-dir {} --sample-sheet {}".format(
            bcl2fastq, runfolder_dir, output_dir, sample_sheet)
        if args != '':
            command+=" {}".format(args)
        default_additions = []
        for x in optional_argument_defaults.keys():
            if x not in command:
                default_additions.append("--{} {}".format(x, optional_argument_defaults[x]))
        if len(default_additions) > 0:
            command = command+" "+" ".join(default_additions)
        self.command = command
        
    def workflow(self):
        self.flowLog('bcl2fastq command: %s' % self.command)
        self.addTask('bcl2fastq', self.command, nCores=self.max_job_cores)
        
if __name__ == "__main__":
    pass