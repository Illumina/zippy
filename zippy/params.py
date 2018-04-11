'''
params.py contains the code to parse a compiled ZIPPY .json file.  It performs several layers of functionality
compared to a pure json loader, such as setting up special keywords 'self' and 'optional', and filling in wildcards.
'''
from __future__ import print_function
import copy
import commentjson as json
import inspect
import sys
from .utils import sub_check_wilds
from .modular_runner import ModularRunner

class ObjectDict(object):
    '''
    It is required that all bound methods of ObjectDict be private (start with _)
    Note: need to figure out how to deal with nested function names.  Like an if not in dir(this object) or something?
    '''
    def __init__(self, **kwargs):
        self._self_flag = False
        for (k,v) in kwargs.iteritems():
            setattr(self, k, v)

    def _lookup_identifier_doublestack(self):
        try:
            stack = inspect.stack()
            parentframe = stack[2][0]
            print(parentframe.f_locals['self'])
            if isinstance(parentframe.f_locals['self'], ModularRunner):
                return parentframe.f_locals['self'].identifier
            else:
                return None
        except (AttributeError,KeyError):
            return None


    def _update_from_dict(self, d):
        if self._self_flag:
            identifier = self._lookup_identifier_doublestack()
            if identifier is not None:
                stage_params = object.__getattribute__(self, identifier)
                stage_params._update_from_dict(d)
            self._self_flag = False
        else:
            for (k,v) in d.iteritems():
                if not hasattr(self, k):
                    setattr(self, k, v)


    def __getattribute__(self, attr):
        '''
        When .self is encountered, a flag (self.self_mode) is set so that the next non-'self' getattribute will be in 'self mode'.
        Thus, we first check in the params.self namespace, and then (if not in self mode) look in the broader namespace.
        After every parameter lookup, it sets the self_mode flag to false
        Downside: lookups to the params object are not threadsafe.  This shouldn't be a big issue.
        '''
        if attr == 'self':
            self._self_flag = True
            return self
        elif attr == 'optional':
            return self
        elif attr[0] == '_': # no magic for private references
            return object.__getattribute__(self, attr)
        else:
            try:
                stack = inspect.stack()
                parentframe = stack[1][0]
                if isinstance(parentframe.f_locals['self'], ModularRunner):
                    identifier = parentframe.f_locals['self'].identifier
                else:
                    self._self_flag = False                
                    return object.__getattribute__(self, attr)
            except (AttributeError,KeyError): # This handles other failures where the calling frame of reference is not a stage
                self._self_flag = False                
                return object.__getattribute__(self, attr)
            stage_params = object.__getattribute__(self, identifier)
            try: #check in self.params.self
                value =  object.__getattribute__(stage_params, attr)
                self._self_flag = False
                return value
            except AttributeError:
                if not self._self_flag: # if we have not required self, we can check the top level
                    self._self_flag = False
                    return object.__getattribute__(self, attr)
                else:
                    self._self_flag = False
                    raise AttributeError('Stage level parameter {} not found'.format(attr))

    def __repr__(self):
        return str(self.__dict__)


def _json_object_hook(d, fname, proto):
    """
    Turns the dictionary into a simple object, emulating a Namespace object from argparse.
    """
    if not proto and 'stages' in d:
        for (i,stage) in enumerate(copy.copy(d['stages'])): 
            if stage['identifier'] in d:
                raise KeyError('Stage identifier {} is the same as a parameter name.  This is not allowed.'.format(stage['identifier']))
            print(stage)
            d[stage['identifier']] = ObjectDict(**stage)
            d['stages'][i] = d[stage['identifier']]
    if 'docker' in d:
        for docker in copy.copy(d['docker']).keys():
            d['docker'][docker] = ObjectDict(**d['docker'][docker])
        d['docker'] = ObjectDict(**d['docker'])
    d['params'] = fname
    return ObjectDict(**d)

def merge_defaults(defaults, params):
    """
    There's some subtlety to this.  Currently, we merge over:
    1.  All top-level keys
    2.  Every stage matched by it's .stage 
    """
    for k in defaults.keys():
        if k not in params:
            params[k] = defaults[k]
    for (i,stage) in enumerate(copy.copy(params['stages'])):
        if stage['stage'] not in defaults['stages']:
            continue
        for (k,v) in defaults['stages'][stage['stage']].iteritems():
            if k not in stage:
                params['stages'][i][k] = v
    return params

def resolve_wildcards(wildcard_map, params):
    new_params_map = {}
    for (k, v) in params.iteritems():
        k = sub_check_wilds(wildcard_map, k)
        if isinstance(v, dict):
            v = resolve_wildcards(wildcard_map, v)
        elif isinstance(v, list):
            v = resolve_wildcards_list(wildcard_map, v)
        elif isinstance(v, str) or isinstance(v, unicode):
            v = sub_check_wilds(wildcard_map, v)
        new_params_map[k] = v
    return new_params_map

def resolve_wildcards_list(wildcard_map, params_list):
    '''
    Helper function for resolve wildcards that iterates over a list, not a dict.
    '''
    new_members = []
    for x in params_list:
        if isinstance(x, dict):
            new_members.append(resolve_wildcards(wildcard_map, x))
        elif isinstance(x, list):
            new_members.append(resolve_wildcards_list(wildcard_map, x))
        else:
            new_members.append(sub_check_wilds(wildcard_map, x))
    return new_members


def old_get_params(fname, defaults_fname=None, proto=False):
    """
    Returns an object with fields defined for every value in the json file.  Hardcodes its own file location
    to the 'params' entry.
    proto: whether this is a prototype params file or not.  If it is, then we do not unroll stages (because it is flat)
    """
    with open(fname, 'rb') as f:
        params = json.load(f)
        if 'wildcards' in params:
            params = resolve_wildcards(params['wildcards'], params)
        if defaults_fname is not None:
            with open(defaults_fname, 'rb') as def_f:
                default_params = json.load(def_f)
                params = merge_defaults(default_params, params)
        params = _json_object_hook(params, fname, proto=proto)
        return params

def get_params(fname, proto=False):
    '''
    To build a params object:
    0.  Turn the json into a python object
    1.  Resolve wildcards
    2.  Build dockers
    3.  Build stages.
        3a.  Expand subworkflows
        3b.  Turn stage_dict into an ObjectDict
        3c.  Update stage_dict and docker_dict from stage environment
    4.  Build global
        4a.  Turn d into an ObjectDict
        4b.  Update global environment
    '''
    subworkflow_index = 0
    with open(fname, 'rb') as f:
        d = json.load(f)
    if proto:
        return ObjectDict(**d)
    if 'wildcards' in d:
        d = resolve_wildcards(d['wildcards'], d)
    if 'docker' in d:
        for docker in copy.copy(d['docker']).keys():
            d['docker'][docker] = ObjectDict(**d['docker'][docker])
        d['docker'] = ObjectDict(**d['docker'])
    if 'stages' in d:
        stages = copy.deepcopy(d['stages'])
        d['stages'] = []
        for stage in stages:
            if 'subworkflow' in stage:
                handle_subworkflow(stage, d, subworkflow_index)
                subworkflow_index += 1
            elif 'identifier' in stage:
                stage_dict = ObjectDict(**stage)
                handle_stage(stage_dict, d)
            else:
                raise ValueError('Stage must either define identifier (and be a stage) or subworkflow (and be a subworkflow)')
    d['params'] = fname
    d = ObjectDict(**d)
    if hasattr(d, 'environment'):
        env = get_params(d.environment)
        if hasattr(env, 'docker'):
            d.docker._update_from_dict(env.docker.__dict__)
        d._update_from_dict(env.__dict__)
    return d

def handle_stage(stage_dict, d):
    '''
    Adds a stage to the params workflow representation, and uses the environment
    file to fill in any missing parameters (if defined)
    '''
    if stage_dict.identifier in d:
        raise KeyError('Stage identifier {} is the same as a parameter name.  This is not allowed.'.format(stage['identifier']))
    # set up an alias for the stage at params.stage_identifier.
    d[stage_dict.identifier] = stage_dict
    d['stages'].append(stage_dict)
    if hasattr(stage_dict, 'environment'):
        env = get_params(stage_dict.environment)
        stage_dict._update_from_dict(getattr(env, stage_dict.identifier).__dict__)
        if hasattr(env, 'docker'):
            d['docker']._update_from_dict(env.docker.__dict__)

def handle_subworkflow(stage, d, subworkflow_index):
    '''
    Adds a subworkflow to the params representation by iteratively adding its stages.
    '''
    subworkflow = get_params(stage['subworkflow'])
    subworkflow_feedins = stage['previous_stage'] if 'previous_stage' in stage else {}
    subworkflow_identifier = stage['identifier'] if 'identifier' in stage else subworkflow_index
    for stage_dict in subworkflow.stages:
        # we add the subworkflow identifier to internal previous_stage references
        if hasattr(stage_dict, 'previous_stage') and stage_dict.previous_stage !='':
            if isinstance(stage_dict.previous_stage, list):
                stage_dict.previous_stage = [x+'{}'.format(subworkflow_identifier) for x in stage_dict.previous_stage]
            else:
                stage_dict.previous_stage+='{}'.format(subworkflow_identifier)
        # we add external previous_stages defined in the feedins
        if stage_dict.identifier in subworkflow_feedins:
            add_feedins_to_previous_stage(stage_dict, subworkflow_feedins[stage_dict.identifier])
        # by default, we append the subworkflow identifier to output paths.
        if not getattr(stage_dict, 'zippy_do_not_append_identifier_to_paths', False):
            if hasattr(stage_dict, 'output_dir'):
                stage_dict.output_dir+='{}'.format(subworkflow_identifier)
        # now we add the stage to the workflow
        stage_dict.identifier+='{}'.format(subworkflow_identifier)
        handle_stage(stage_dict, d)


def add_feedins_to_previous_stage(stage_dict, feedins):
    '''
    We add the new feedins into the dependencies of the curent stage.
    '''
    if feedins == '': #there's a feedin defined, but it's not an identifier.  We can ignore it.
        return ''
    if not isinstance(feedins, list):
        feedins = [feedins]
    if isinstance(stage_dict.previous_stage, list): # formerly, there were multiple previous stages
        stage_dict.previous_stage.extend(feedins)
    elif stage_dict.previous_stage != '': # formerly, there was one previous stage
        stage_dict.previous_stage = [stage_dict.previous_stage]+feedins
    else: # formerly, there was no prvious stage
        stage_dict.previous_stage = feedins

def save_params_as_json(params, fname):
    params = rejsonify(params)
    for stage in params['stages']: #we delete the duplicate copies of the stages, to match the input json style
        del params[stage['identifier']]
    with open(fname, 'w') as f:
        json.dump(params, f, indent=4)

def rejsonify(params):
    '''
    Used to convert ObjectDicts to serializable dicts
    '''
    if isinstance(params, ObjectDict):
        to_return = params.__dict__
        #strip out private members
        to_return = {k:v for (k,v) in to_return.iteritems() if k[0]!='_'}
        for (k,v) in copy.copy(to_return).iteritems():
            to_return[k] = rejsonify(v)
    elif isinstance(params, list):
        for (i,v) in enumerate(copy.copy(params)):
            params[i] = rejsonify(v)
        to_return = params
    else:
        to_return = params
    return to_return
