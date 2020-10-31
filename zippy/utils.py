import re
def sub_check_wilds(wildcard_map, s):
    '''
    Given a dictionary of substitutions to make, and a string <s> to find such substitutions, return
    the substituted string.
    '''
    for wildcard in wildcard_map.keys():
        s = s.replace('{'+wildcard+'}', wildcard_map[wildcard])
    return s

def sub_check_params(params, s):
    '''
    Given a params object of substitutions to make, and a string <s> to find such substitutions, return
    the substituted string.
    '''
    wildcards = re.findall('\{(\S+)\}', s)
    print(wildcards)
    for wildcard in wildcards:
        print(wildcard)
        if hasattr(params, wildcard):
            s = re.sub('{'+wildcard+'}', getattr(params, wildcard), s)
    return s
