import argparse
import ast
import _ast
import imp
import inspect
try:
    import meta #meta is optional, and used for parameter detection if the source code is not found.
except ImportError:
    pass
import sys
import re
from collections import OrderedDict, namedtuple
import commentjson as json
from .modular_runner import *
from .params import get_params

param_tuple = namedtuple('param_tuple', ['param', 'is_optional'])
class ParamFinder(ast.NodeVisitor):
    '''
    AST graph walker designed to find all uses of the params namespace in the ast graph.
    There is some magic in handling the 'self' keyword, which is used to denote class-specific params.
    For example, params.self.output_dir in class RSEMRunner will be mapped to the tuple ('RSEMRunner', 'output_dir')
    '''
    def __init__(self):
        self.params_found = []

    def filter_private(self):
        print self.params_found
        self.params_found = [x for x in self.params_found if x.param[-1][0]!='_']

    def register_param(self, raw_tuple):
        '''
        Creates a list of params that keeps the order they were seen, but removes duplicates.
        '''
        if self.is_valid_param(raw_tuple):
            self.params_found.append(self.format_param(raw_tuple))

    def is_valid_param(self, x):
        '''
        Given a raw tuple formed from an AST path, we see if it's well formed.  A well formed
        tuple is currently defined as any param that is not wholly containing self and optional
        tags.  More sophisticated validation here is possible.
        '''
        if len(x) == 1 and x[0] in ['self', 'optional']:
            return False
        if len(x) == 2 and x[0] in ['self', 'optional'] and x[1] in ['self', 'optional']:
            return False
        return True

    def format_param(self, raw_tuple):
        '''
        Converts a raw ast trace from the ast graph into a param_tuple named tuple.  It:
        -strips out optional and puts it in its own field
        -converts self into the current class name. 
        '''
        try:
            raw_tuple.remove('optional')
            is_optional = True
        except ValueError: #optional is not in the param.
            is_optional = False
        try:
            while True:
                self_index = raw_tuple.index('self')
                raw_tuple[self_index] = self.current_class
        except ValueError:
            pass
        if len(raw_tuple) > 2:
            raise ValueError('Malformed parameter tuple: {}'.format(raw_tuple))
        return param_tuple(tuple(raw_tuple), is_optional)




    def uniqueify(self):
        seen = set()
        self.params_found = [x for x in self.params_found if x not in seen and not seen.add(x)]

    def old_visit_Attribute(self, node):
        '''
        TODO: check that it's read and not written
        '''
        if isinstance(node.value, _ast.Attribute) and node.value.attr == 'params' and node.attr != 'self':
            self.params_found.append(node.attr)
        #This next bit identifies lines two steps removed from params.  In that case, we know that the middle attribute is 'self' and we 
        #replace it with the relevant class name.
        elif isinstance(node.value, _ast.Attribute) and isinstance(node.value.value, _ast.Attribute) and node.value.value.attr == 'params':
            self.params_found.append((self.current_class, node.attr))
        self.generic_visit(node)

    def visit_Attribute(self, node):
        '''
        TODO: check that it's read and not written
        '''
        if isinstance(node, _ast.Attribute):
            current_node = node
            param_builder = []
            param_builder.append(node.attr)
            while isinstance(current_node.value, _ast.Attribute):
                if current_node.value.attr == 'params':
                    self.register_param(param_builder[::-1])
                    break
                param_builder.append(current_node.value.attr)
                current_node = current_node.value
        self.generic_visit(node)

def case_insensitive_list_match(query, l):
    """
    Returns the first case-insensitive matching item to the query
    """
    query = query.lower()
    l_low = [x.lower() for x in l]
    for (i,x) in enumerate(l_low):
        if query == x:
            return l[i]
    return None

def walk_all_methods(pf, class_object):
    '''
    '''
    for func in dir(class_object):
        try:
            try:
                #we get each method from source, deindent it, and then parse it
                source_lines = inspect.getsourcelines(getattr(class_object,func))[0]
                indent = len(source_lines[0]) - len(source_lines[0].lstrip())
                source_lines = [line[indent:] for line in source_lines]
                ast_tree = ast.parse(''.join(source_lines))
            except IOError:
                print('Module source code not found: ({}, {}).  Decompiling instead.'.format(class_object, func))
                try:
                    ast_tree = meta.decompile(getattr(class_object,func).__code__)
                except AssertionError:
                    print 'meta failed to decompile function {} in class {}.  Some parameters may be missing from the generated file.'.format(func, class_object.__name__)
                    continue
                except NameError:
                    print 'meta is not installed.  Parameters may be missing from the generated file.'
            except TypeError:
                continue
            pf.visit(ast_tree)
        except AttributeError:
            continue

class JsonBuilderMain():
    """
    Compiles a params proto file into a template that must be filled out for running ZIPPY.
    """
    def __init__(self, input_args):
        self.input_args = input_args
        if self.input_args.out:
            self.output_file = self.input_args.out
        else:
            arg_split = self.input_args.proto.split('.')
            self.output_file = '.'.join(arg_split[:-1]+['compiled']+ [arg_split[-1]])
        self.params = get_params(self.input_args.proto, proto=True)
        self.modules_loaded = 0
        for x in getattr(self.params, 'imports', []):
            self.load_external_modules(x)

    def load_external_modules(self, module_path):
        modules_to_add = {}
        m = imp.load_source('user_module_import_{}'.format(self.modules_loaded), module_path)
        self.modules_loaded += 1
        for k in vars(m).keys():
            if 'Runner' in k:
                modules_to_add[k] = vars(m)[k]
        globals().update(modules_to_add)

    def make_params_file(self):
        '''
        Given the params.stages list in the input params file, finds all the needed params to run the pipeline.
        Each stage will have an identifier, and may have a previous stage.  We also add a 'scratch_path' for the pyflow workspace.
        '''
        pf = ParamFinder()
        for stage in self.params.stages:
            class_name = case_insensitive_list_match('{}Runner'.format(stage), globals().keys())
            try:
                class_object = globals()[class_name]
            except KeyError:
                raise KeyError("One of your workflow stages, {}, is not found.  If it is defined in a custom file, make sure that file is imported in your params file.  If it is built-in, make sure you are spelling it correctly!".format(stage))
            pf.current_class = stage
            walk_all_methods(pf, class_object)
        pf.uniqueify()
        pf.filter_private()
        try:
            defaults = get_params(self.input_args.defaults, proto=True)
        except IOError:
            print 'Warning: No defaults file found.'
            defaults = None
        output_map = OrderedDict()
        if hasattr(self.params, 'imports'):
            output_map["imports"] = self.params.imports
        identifiers = set()
        optional_params = filter(lambda x: x.is_optional, pf.params_found)
        required_params = filter(lambda x: not x.is_optional, pf.params_found)
        output_map["stages"] = []
        #stage params
        for (i,stage) in enumerate(self.params.stages):
            identifier = get_identifier(stage, identifiers)
            stage_map = {"stage": stage, "identifier": identifier}
            if i > 0:
                print output_map
                stage_map['previous_stage'] = output_map["stages"][i-1]["identifier"]
            for param in required_params:
                param = param.param
                print param
                if param[0] == stage and len(param)>1:
                    if hasattr(defaults, 'stages') and stage in defaults.stages and param[1] in defaults.stages[stage]:
                        stage_map[param[1]] = defaults.stages[stage][param[1]]
                    else:
                        stage_map[param[1]] = ''
            for param in optional_params:
                if self.input_args.default_behavior == 'ignore':
                    continue
                param = param.param
                if param[0] == stage:
                    if self.input_args.default_behavior == 'include':
                        if stage in defaults.stages and param[1] in defaults.stages[stage]:
                            stage_map[param[1]] = defaults.stages[stage][param[1]]
                        else:
                            stage_map[param[1]] = ''
                    elif self.input_args.default_behavior == 'warn':
                        if not hasattr(defaults, 'stages') or stage not in defaults.stages or param[1] not in defaults.stages[stage]:
                            print "Warning: parameter {} is not included in stage {} defaults".format(param[1], stage)
            output_map["stages"].append(stage_map)
        #global params
        for param in required_params:
            if len(param.param) > 1:
                continue
            param = param.param[0]
            if hasattr(defaults, param):
                output_map[param] = getattr(defaults, param)
            else:
                output_map[param] = ''
        for param in optional_params:
            if len(param.param) > 1:
                continue
            if self.input_args.default_behavior == 'ignore':
                    continue
            param = param.param[0]
            if self.input_args.default_behavior == 'include':            
                if hasattr(defaults, param):
                    output_map[param] = getattr(defaults, param)
                else:
                    output_map[param] = ''
            elif self.input_args.default_behavior == 'warn':            
                if not hasattr(defaults, param):
                    print "Warning: global parameter {} is not included".format(param)
        output_map["scratch_path"] = ''
        with open(self.output_file, 'w') as f:
            json.dump(output_map, f, indent=4)

def get_identifier(stage, identifiers):
    '''
    Called when beginning the output for a new pipeline stage.  Returns a unique id for that pipeline stage.
    '''
    while stage in identifiers:
        instance = re.search('(.+_)(\d+)$', stage)
        if instance is not None: #we are an identifer of the form stagename.#
            stage = instance.group(1)+str(int(instance.group(2))+1)
        else: #we are an identifier of the form stagename
            stage = stage+'_1'
    identifiers.add(stage)
    return stage

def get_argparser():
    parser = argparse.ArgumentParser(usage='make_params is used to compile a proto-workflow into a form useable by ZIPPY', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('proto', help='The proto-workflow to compile.')
    parser.add_argument('out', help='The name of the compiled workflow to create.  If none, prepends compiled to the filetype as output (e.g., foo.compiled.json)', nargs='?')
    parser.add_argument('--defaults', default='defaults.json', help='File with default parameter values (in json format)')
    parser.add_argument('--default_behavior', default='warn', choices=['include', 'warn', 'ignore'], help='How to treat optional parameters.  <Include> them all in your file, <warn> you if they are not in the specified defaults file, or <ignore> them.')
    return parser
if __name__ == '__main__':    
    parser = get_argparser()
    input_params = parser.parse_args()
    jsb = JsonBuilderMain(input_params)
    jsb.make_params_file()
