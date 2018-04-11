from __future__ import print_function
from .modular_runner import *
from .make_params import case_insensitive_list_match, ParamFinder, walk_all_methods

def run_help(arg):
    print('Welcome to the ZIPPY help system.')
    if arg is None:
        print('ZIPPY is a powerful, easy-to-use NGS pipeline prototyping system, with batteries included.')
        print('This help system can be used to learn more about the various modules included with ZIPPY')
        print(' ')
        print('To see a list of all modules, use "python zippy.py --help all"')
        print('To get help for a specific module, use "python zippy.py --help <modulename>"')
        print(' ')
        print('ZIPPY command line arguments:')
        print(' ')
        print('zippy.py workflow [--defaults defaults_file]')
        print('workflow\tThe ZIPPY workflow json file.')
        print('--defaults\tDEFAULTS\tFile with default parameter values (in json format)\n\t\t(default: defaults.json)')
        print('--run_local\tIf set, we will run things on the local node, rather than attempt to use SGE.')
    elif arg == 'all':
        ignore_set = set(['ModularRunner', 'WorkflowRunner'])
        print('All built-in ZIPPY modules:')
        print(' ')
        print('\n'.join(sorted([x for x in globals().keys() if 'Runner' in x and x not in ignore_set])))
    else:
        if 'Runner' not in arg:
            arg = case_insensitive_list_match('{}Runner'.format(arg), globals().keys())
        module = globals()[arg]
        print(module)
        print(module.__doc__)
        print("\tParameters for {}:".format(arg))
        pf = ParamFinder()
        pf.current_class = arg
        walk_all_methods(pf, module)
        pf.uniqueify()
        for param in pf.params_found:
            if param.param[0] == arg:
                if param.is_optional:
                    print("\t\t{}\t(optional)".format(param.param[1]))
                else:
                    print("\t\t{}".format(param.param[1]))



if __name__ == '__main__':
    main()
