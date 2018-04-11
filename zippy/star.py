import os
import re
import random
import string
import subprocess
from collections import defaultdict
from pyflow import WorkflowRunner

char_set = string.ascii_uppercase + string.ascii_lowercase + string.digits

def find_fastq_files(fastq_dir, rm_undetermined=False, no_crawl=False):
    '''
    @author: jtoung
    '''
    fastq_files = {}
    if no_crawl:
        find = subprocess.Popen(['find', '-L', fastq_dir, '-maxdepth', '1', '-name', '*.fastq.gz'], stdout=subprocess.PIPE)
    else:
        find = subprocess.Popen(['find', '-L', fastq_dir, '-name', '*.fastq.gz'], stdout=subprocess.PIPE)
    fastq_regex = re.compile("^(.*)_(S[0-9]+)_([L0-9\-]+)?_?(R[1|2])_001\.fastq\.gz$")
    for fastq_file in find.stdout:
        fastq_file = fastq_file.rstrip('\n').split()[0]
        fastq_file_basename = os.path.basename(fastq_file)
        if fastq_regex.match(fastq_file_basename):
            sample_name = fastq_regex.match(fastq_file_basename).group(1)
            sample_id = fastq_regex.match(fastq_file_basename).group(2)
            lane = fastq_regex.match(fastq_file_basename).group(3)
            read = fastq_regex.match(fastq_file_basename).group(4)
            
            if sample_name == "Undetermined":
                continue
            
            if not lane:
                lane = 'None'
                
            if lane not in fastq_files:
                fastq_files[lane] = {}
            
            sample = (sample_id, sample_name)
            if sample not in fastq_files[lane]:
                fastq_files[lane][sample] = {}
            
            if read in fastq_files[lane][sample]:
                raise Exception('already observed fastq file for lane %s, sample_id %s, sample_name %s, read %read in fastq_dir %s' % (lane, sample_id, sample_name, read, fastq_dir))
            
            fastq_files[lane][sample][read] = fastq_file
    
    return fastq_files

class starFlow(WorkflowRunner):
    """
    Pyflow workflow to call the STAR aligner.  We follow file formatting guidelines from the isaac workflow
    from bioxutils.
    """
    def __init__(self, star_path, star_index, fastq_dir, out_dir, max_job_cores=8):
        self.star_path = star_path
        self.star_index = star_index
        self.fastq_dir = fastq_dir
        self.out_dir = out_dir
        self.max_job_cores = max_job_cores

    def workflow(self):
        if not os.path.exists(self.out_dir):
            os.makedirs(self.out_dir)
        project_name = 'autogen_project'
        fastq_files = find_fastq_files(self.fastq_dir, rm_undetermined=True)
        
        mv_fastq_tasks = []
        fastq_cleanup_cmds = []
        sample_id_TO_sample_name = {}
        fastq_by_sample_1 = defaultdict(list)
        fastq_by_sample_2 = defaultdict(list)        
        for (lane, sample_data) in fastq_files.iteritems():
            for (sample, reads) in sample_data.iteritems():
                    self.flowLog(reads.keys())
                    fastq_by_sample_1[sample].append(reads['R1'])
                    try:
                        fastq_by_sample_2[sample].append(reads['R2'])
                    except KeyError:
                        no_r2 = True
        for sample in fastq_by_sample_1.keys():
            (sample_id, sample_name) = sample
            sample_path = os.path.join(self.out_dir, '_'.join(sample))
            if not os.path.exists(sample_path):
                os.makedirs(sample_path)
            # STAR alignment
            if no_r2:
                read_files_string = ','.join(fastq_by_sample_1[sample])
            else:
                read_files_string = ','.join(fastq_by_sample_1[sample])+" "+','.join(fastq_by_sample_2[sample])
            self.star_command = "{star_path} --genomeDir {star_index}  --readFilesIn {read_files_string} --runThreadN {max_job_cores}  --outFileNamePrefix {out_dir}  --outTmpDir {tmp_dir} --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within".format(
                star_path=self.star_path, star_index=self.star_index, 
                read_files_string=read_files_string, max_job_cores=self.max_job_cores,
                out_dir=sample_path+"/", tmp_dir = os.path.join(self.out_dir, 'tmp'+str(sample_id)))
            self.flowLog(self.star_command)
            star_task = self.addTask("star"+str(sample_id), self.star_command, nCores=self.max_job_cores, memMb=40000, dependencies=mv_fastq_tasks)
            # cleanup
            self.addTask("delete_tmp"+str(sample_id),  "rm -r {tmp_dir}".format(tmp_dir=os.path.join(self.out_dir, 'tmp'+str(sample_id))), dependencies=star_task)
            rename_task = self.addTask("rename_star"+str(sample_id),  "mv {path}/Aligned.sortedByCoord.out.bam {path}/{sample_id}_{sample_name}.raw.bam".format(path=sample_path, sample_id=sample_id, sample_name=sample_name), dependencies=star_task)
            # build bai
            self.addTask("bai"+str(sample_id), "samtools index {path}/{sample_id}_{sample_name}.raw.bam".format(path=sample_path, sample_id=sample_id, sample_name=sample_name), dependencies=rename_task)

class SingleStarFlow(WorkflowRunner):
    """
    Pyflow workflow to call the STAR aligner.  We follow file formatting guidelines from the isaac workflow
    from bioxutils.
    """
    def __init__(self, star_path, star_index, sample, fastqs, out_dir, max_job_cores=8, tmp_path='', command_args=''):
        self.star_path = star_path
        self.star_index = star_index
        self.sample = sample
        self.fastqs = fastqs
        self.out_dir = out_dir
        self.max_job_cores = max_job_cores
        if tmp_path !='':
            self.tmp_path = tmp_path
        else:
            self.tmp_path = os.path.join(self.out_dir, 'tmp'+self.sample)
        self.tmp_path+=''.join(random.choice(char_set) for _ in range(4))
        self.command_args = command_args

    def workflow(self):
        if not os.path.exists(self.out_dir):
            os.makedirs(self.out_dir)
        pre_out_prefix = os.path.join(self.out_dir, self.sample)
        fastq_files = self.fastqs
        fastq_cleanup_cmds = []
        fastq_by_sample_1 = [x for x in fastq_files if '_R1_' in x]
        fastq_by_sample_2 = [x for x in fastq_files if '_R2_' in x]
        if len(fastq_by_sample_2) > 0:
                read_files_string = ','.join(fastq_by_sample_1)+" "+','.join(fastq_by_sample_2)
        else:
                read_files_string = ','.join(fastq_by_sample_1)
        self.star_command = "{star_path} --genomeDir {star_index}  --readFilesIn {read_files_string} --runThreadN {max_job_cores}  --outFileNamePrefix {out_dir}  --outTmpDir {tmp_dir} --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within".format(
            star_path=self.star_path, star_index=self.star_index, 
            read_files_string=read_files_string, max_job_cores=self.max_job_cores,
            out_dir=pre_out_prefix, tmp_dir = self.tmp_path)
        if self.command_args != "":
            self.star_command+=" "+self.command_args
        self.flowLog(self.star_command)
        star_task = self.addTask("star", self.star_command, nCores=self.max_job_cores, memMb=100*1024, dependencies=None)
        #self.addTask("delete_tmp",  "rm -r {tmp_dir}".format(tmp_dir=os.path.join(self.out_dir, 'tmp'+self.sample)), dependencies=star_task)
        rename_task = self.addTask("rename_star"+str(self.sample),  "mv {prefix}Aligned.sortedByCoord.out.bam {path}/{sample}.raw.bam".format(prefix=pre_out_prefix, path=self.out_dir, sample=self.sample), dependencies=star_task)
        # build bai
        self.addTask("bai", "samtools index {path}/{sample}.raw.bam".format(path=self.out_dir, sample=self.sample), dependencies=rename_task)


    

if __name__ == '__main__':
    wflow = STARFlow()
    retval = wflow.run(mode='sge')
    sys.exit(retval)