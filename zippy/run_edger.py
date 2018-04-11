'''
Runs edger from a set of RSEM inputs
'''
import os
import random
import subprocess
import csv

def parse_rsem(f_name):
    '''
    Returns a map from gene name to a list with one item, the expected count.
    '''
    rsem_map = {}
    with open(f_name, 'rb') as f:
        c = csv.reader(f, delimiter='\t')
        c.next()
        for line in c:
            gene = line[0]
            rsem_map[gene] = [float(line[4])]
    return rsem_map

def add_to_data_grid(data_grid, sample_map):
    if data_grid is None:
        data_grid = sample_map
    else:
        for k in sample_map.keys():
            try:
                data_grid[k].extend(sample_map[k])
            except KeyError:
                continue
    return data_grid

def write_data_grid(data_grid, key_ordering, temp_folder):
    #we write out every gene that has the 
    random_path_suffix = str(random.randint(0,100000))
    data_grid_path = os.path.join(temp_folder, 'datagrid'+random_path_suffix+'.csv')
    print 'DGP: {}'.format(data_grid_path)
    with open(data_grid_path, 'wb') as f:
        c = csv.writer(f)
        c.writerow(['Symbol']+key_ordering)
        for symbol in data_grid.keys():
            if len(data_grid[symbol]) == len(key_ordering):
                c.writerow([symbol]+data_grid[symbol])
            else:
                print len(data_grid[s])
    return data_grid_path

def run_edger(params):
    data_grid = None
    key_ordering = []
    for (group_index, samples) in enumerate([params.group1, params.group2]):
        for sample in samples:
            sample_map = parse_rsem(sample)
            data_grid = add_to_data_grid(data_grid, sample_map)
            key_ordering.append(os.path.basename(sample))
    data_grid_path = write_data_grid(data_grid, key_ordering, params.temp_folder)
    subprocess.call('{r_path} edger.r {data_grid_path} {group1_len} {group2_len} {out_file} {temp_path}'.format(
        data_grid_path=data_grid_path,
        group1_len=len(params.group1),
        group2_len=len(params.group2),
        out_file=params.out_file,
        r_path=params.r_path,
        temp_path=params.temp_folder), shell=True)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--group1', nargs='+')
    parser.add_argument('--group2', nargs='+')
    parser.add_argument('--temp_folder')
    parser.add_argument('--out_file')
    parser.add_argument('--r_path')
    params = parser.parse_args()
    run_edger(params)