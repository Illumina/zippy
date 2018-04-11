from __future__ import print_function
import re
import json
import copy
import sys
import subprocess
import imp
import os
import types
from functools import partial
from pyflow import WorkflowRunner
from .modular_runner import *
from .modular_helper import run_help
from .params import get_params, save_params_as_json
from .make_params import case_insensitive_list_match


def docker_task_reformat(params, current_stage_setup, command):
    '''
    Pulls docker information from the json file and wraps the command appropriately.
    '''
    stage_params = getattr(params, current_stage_setup)
    if hasattr(stage_params, 'use_docker'):
        docker_config = getattr(params.docker,stage_params.use_docker)
        mount_points = []
        for m in docker_config.mount_points:
            mount_points.append('{prefix} {mount_dir}'.format(
                    prefix=_get_docker_mount_prefix(docker_config), mount_dir=m))
            if getattr(docker_config, 'intercept_mounts', False):
                (from_path,to_path) = m.split(':')
                command = re.sub(r"(\s+)"+from_path, r'\1'+to_path, command)
                command = re.sub("^"+from_path, to_path, command)
        mount_points = " ".join(mount_points)
        if getattr(docker_config, 'pull', True):
            load_command = _get_docker_pull_command(docker_config)                
        else:
            load_command = _get_docker_load_command(docker_config)
        command = _get_docker_run_command(docker_config, mount_points, command)
        command = load_command+" ; "+command
    return command

def _get_docker_mount_prefix(docker_config):
    if docker_config.engine == 'singularity':
        return '-B'
    else:
        return '-v'

def _get_docker_load_command(docker_config):
    if docker_config.engine == 'singularity':
        raise NotImplementedError('must use "pull" mode with singularity.  Set "pull" to true in the config.')
    else:
        return "cat {} | docker load".format(docker_config.image_file_path)

def _get_docker_pull_command(docker_config):
    if docker_config.engine == 'singularity':
        return 'singularity pull docker://{}'.format(docker_config.image_file_path)
    else:
        return 'docker pull {}'.format(docker_config.image_file_path)

def _get_docker_run_command(docker_config, mount_points, command):
    if docker_config.engine == 'singularity':
        return 'singularity exec {mount_points} docker://{docker_image} {command}'.format(
            mount_points=mount_points, command=command, docker_image=docker_config.image)
    else:
        if hasattr(docker_config, 'uid'):
            uid = docker_config.uid
        else:
            uid = subprocess.check_output('id -u', shell=True).strip()
        return 'docker run --rm {mount_points} --user={uid} --entrypoint=/bin/bash {docker_image} -c "{command}"'.format(
            mount_points=mount_points, command=command, docker_image=docker_config.image, uid=uid)

def subworkflow_addTask(current_stage_setup, params, self, label, command=None, **kwargs):
    '''
    Used to monkey patch pyflow's subworkflows, so that we can inject docker wrappers into them.
    Very similar to the version used for the main workflow, but we have to hardcode current_stage_setup
    and params using partial functions, to fix the context for the subworkflow.
    '''
    if current_stage_setup is not None and command is not None:
        command = docker_task_reformat(params, current_stage_setup, command)
        self.flowLog(command)
    return WorkflowRunner.addTask(self, label, command=command, **kwargs)



class ModularMain(WorkflowRunner):
    """
    ZIPPY, the modular pipeline execution system.
    """
    def __init__(self, input_params, from_dict=False):
        self.modules_loaded = 0
        if from_dict: #initiate from api
            self.params = input_params
        else: #initiate from command line
            if not os.path.exists(input_params.defaults):
                print("Warning: Defaults file {} is not found.".format(input_params.defaults))
                self.params = get_params(input_params.workflow)
            else:
                self.params = get_params(input_params.workflow, input_params.defaults)
        self.stage_dict = {} # map from identifier to stage
        for x in getattr(self.params, 'imports', []):
            self.load_external_modules(x)

    def addWorkflowTask(self, label, workflowRunnerInstance, dependencies=None):
        '''
        We inject a docker-aware addTask function into subworkflows.
        '''
        workflowRunnerInstance.addTask = types.MethodType(partial(subworkflow_addTask, self.current_stage_setup, self.params), workflowRunnerInstance)
        return WorkflowRunner.addWorkflowTask(self, label, workflowRunnerInstance, dependencies=dependencies)


    def addTask(self, label, command=None, **kwargs):
        '''
        We wrap pyflow's addTask function to inject docker commands when applicable.
        self.current_stage_setup is set in run_stage, and is switched as each stage adds its tasks to the DAG.
        '''
        if self.current_stage_setup is not None and command is not None:
            command = docker_task_reformat(self.params, self.current_stage_setup, command)
            self.flowLog(command)
        return WorkflowRunner.addTask(self, label, command=command, **kwargs)

    def run_stage(self, stage, previous_stage):
        '''
        So this is some crazy stuff.  What this does is allow us to name classes as strings
        and then instantiate them as workflow tasks.  Each workflow task is self contained, and takes as
        input its previous stage, which knows where its outputs will be per-sample.  Thus, so long
        as the input/output formats line up, stages can be arbitrarily chained.
        @return the newly instantiated stage
        TODO: use __import__ to allow the importation of arbitrary dependencies.
        '''
        print(stage)
        class_name = case_insensitive_list_match('{}Runner'.format(stage.stage), globals().keys())
        print (class_name, stage, previous_stage)
        try:
            stage_out = globals()[class_name](stage.identifier, self.params, previous_stage)
        except KeyError:
            raise KeyError("Your workflow stage {} is not found.  If it is defined in a custom file, make sure that file is imported in your params file.  If it is built-in, make sure you are spelling it correctly!".format(stage.stage))
        self.current_stage_setup = stage.identifier
        stage_out.setup_workflow(self)
        self.stage_dict[stage.identifier] = stage_out
        self.current_stage_setup = None
        return stage_out

    def load_external_modules(self, module_path):
        modules_to_add = {}
        m = imp.load_source('user_module_import_{}'.format(self.modules_loaded), module_path)
        self.modules_loaded += 1
        for k in vars(m).keys():
            if 'Runner' in k:
                modules_to_add[k] = vars(m)[k]
        globals().update(modules_to_add)



    def workflow(self):
        '''
        For each stage, we add its workflow, passing along its explicit dependencies, if it has any.  If it does not,
        we assume it depends on the previous stage run.
        '''
        previous_stage = None
        for stage in self.params.stages:
            print("Adding stage {}".format(stage.identifier))
            if hasattr(stage, 'previous_stage'):
                if isinstance(stage.previous_stage, list):
                    previous_stages = [self.stage_dict[x] for x in stage.previous_stage]
                    previous_stage = self.run_stage(stage, previous_stages)
                elif isinstance(stage.previous_stage, str):
                    previous_stage = self.run_stage(stage, [self.stage_dict[stage.previous_stage]])
                elif isinstance(stage.previous_stage, unicode):
                    previous_stage = self.run_stage(stage, [self.stage_dict[str(stage.previous_stage)]])
                else:
                    raise TypeError('Previous stage for {} is neither list or string'.format(stage.identifier))
            else:
                previous_stage = self.run_stage(stage, [previous_stage])

    def run_zippy(self, mode='sge', mail_to=None):
        """Runs the workflow. Returns 0 upon successful completion, 1 otherwise"""
        # pyflow.WorkflowRunner's run function by default already returns 0/1 for success/fail
        return self.run(mode=mode, dataDirRoot=self.params.scratch_path, retryMax=0, mailTo=mail_to)

def build_zippy(params_dictionary):
    """
    Experimental! Call ZIPPY from python, specifying a params dictionary instead of going through makeparams.  Returns a zippy object,
    which can then be run using the call x.api_run_zippy()
    """
    return ModularMain(params_dictionary, from_dict=True)

def load_params(fname, defaults_fname=None):
    """
    Loads a valid zippy parameters file or template file from disk.  Represents the file as a native python object (c.f., the python json module)
    """
    if defaults_fname is not None:
        return get_params(fname)
    return get_params(fname, defaults_fname)

def save_params(params, fname):
    """
    Writes a python object (structured as json) to a file.  Used to write files which can then by run using 'python zippy.py your_file.json'
    """
    save_params_as_json(copy.deepcopy(params), fname)

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(usage='ZIPPY will run your workflow compiled by make_params and filled out by YOU!  Run "python zippy.py --help" for more.', formatter_class=argparse.ArgumentDefaultsHelpFormatter, add_help=False)
    parser.add_argument('workflow', help='The ZIPPY workflow json file.', nargs='?')
    parser.add_argument('--defaults', default='defaults.json', help='File with default parameter values (in json format)')
    parser.add_argument('--email', default=None, help='If specified, you will be emailed here upon job completion')
    parser.add_argument('--run_local', dest='run_local', action='store_true', help='Setting this flag will run things on the local node, rather than attempt to use SGE.')
    parser.add_argument('--help', dest='help', action='store_true', help='Our help system!')
    parser.set_defaults(help=False)
    parser.set_defaults(run_local=False)
    input_params = parser.parse_args()
    if input_params.help or input_params.workflow is None:
        run_help(input_params.workflow)
        if not input_params.help and input_params.workflow is not None:
            print(parser.get_help())
    else:
        wflow = ModularMain(input_params)
        if input_params.run_local:
            mode = 'local'
        else:
            mode = 'sge'
        retval = wflow.run(mode=mode, dataDirRoot=wflow.params.scratch_path, retryMax=0, mailTo=input_params.email)
        sys.exit(retval)
