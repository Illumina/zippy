def sub_check_wilds(wildcard_map, s):
    '''
    Given a dictionary of substitutions to make, and a string <s> to find such substitutions, return
    the substituted string.
    '''
    for wildcard in wildcard_map.keys():
        s = s.replace('{'+wildcard+'}', wildcard_map[wildcard])
    return s
