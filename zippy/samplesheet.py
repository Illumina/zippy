'''
@author jtoung
'''
import re
from collections import defaultdict

valid_sections = set(('Header', 'Reads', 'Settings', 'Data'))

def check_valid_samplename(name):
    if "_R1_" in name or "_R2_" in name or name.endswith("_R1") or name.endswith("_R2"):
        raise RuntimeError('Sample name {} contains _R1_ or endswith _R1.  This interferes with how paired-end reads are formatted.  Please rename your samples'.format(name))

def check_valid_section(section):
    if section not in valid_sections:
        raise Exception('invalid section %s; valid sections are %s' % (section, " ".join(valid_sections)))

class SampleSheet(object):

    def __init__(self, sample_sheet, fix_dup_sample_names=False):
        self.unique_sample_name_map = {}
        if not sample_sheet:
            raise AttributeError('samplesheet is required for SampleSheet')
        self.sample_sheet = sample_sheet
        self.fix_dup_sample_names = fix_dup_sample_names
        
        self._parse_sample_sheet(self.sample_sheet)
        self._check_names()

    def _check_names(self):
        '''
        Ensures that names are unique (it either fails, or fixes them
        if fix_dup_sample_names is true.  Also ensures that names fit in
        our sample name guidelines (do not end with _R1/_R2 and do not
        contain _R1_/_R2_)
        '''
        name_to_ids = defaultdict(set)
        ids_to_names = defaultdict(set)
        for sample_id in self.get_sample_ids():
            sample_lines = self.get_samples_for_id(sample_id)
            check_valid_samplename(sample_id)
            for line in sample_lines:
                name = line.get("Sample_Name")
                name_to_ids[name].add(sample_id)
                ids_to_names[sample_id].add(name)
                check_valid_samplename(name)
        for (identifier, name_set) in ids_to_names.iteritems():
            if len(name_set) > 1:
                raise IOError('Sample ID {} has multiple sample names.  Sample name and sample ID must be unique'.format(identifier))        
        for (name, id_set) in name_to_ids.iteritems():
            if len(id_set) > 1:
                if self.fix_dup_sample_names:
                    for sample_id in id_set:
                        self.uniqueify_sample_names(sample_id)
                        if name == '':
                            self.overwrite_sample_name(sample_id)
                else:
                    raise IOError('Sample name {} has multiple sample IDs.  Sample name and sample ID must be unique'.format(name))
            else:
                self.unique_sample_name_map[list(id_set)[0]] = name

    def sample_id_to_sample_index(self, sample_id):
        sample_indices = 1
        samples_seen = set()
        for line in self._ss["Data"]:
            current_id = line.get("Sample_ID")
            if current_id == sample_id:
                return sample_indices
            else:
                if current_id not in samples_seen:
                    sample_indices += 1
                    samples_seen.add(current_id)
        raise ValueError('Sample id {} not found in samplesheet'.format(sample_id))


    def uniqueify_sample_names(self, sample_id):
        self.unique_sample_name_map[sample_id] = sample_id

    def overwrite_sample_name(self, sample_id):
        sample_lines = self.get_samples_for_id(sample_id)
        for line in sample_lines:
            line.set("Sample_Name", sample_id)

    def print_sample_sheet(self, output):
        """prints sample sheet to output"""
        fh_out = open(output, 'wb')
        for line in self._l:
            fh_out.write(",".join(line) + "\n")
        fh_out.close()
        
    def _parse_sample_sheet(self, sample_sheet):
        """parse the sample sheet into two objects: _l is a list of each line in the original file; _ss is a dictionary that lets us look up the index of each line conveniently"""
        """first, _ss points to sections Header, Reads, or Settings. the value of these keys are themselves dictionaries whose keys are the first column, and values are the line idx of the original file, which let's us get the corresponding line in _l"""
        """second, _ss points to sections Data, which contains a list of "Data" elements, each of which contain a link to the original sample sheet (self), as well as the line idx of the original file"""
        
        header_regex = re.compile("^\[(.*)\]$")
        section_name = None
        # _ss is a dictionary that stores keys called "Header", "Reads", or "Settings"
        # 
        self._ss = {}
        # _l is a list that of the original file
        self._l = []
        _data_header = None
        _data_i = 1
        with open(sample_sheet, 'rb') as f:
            for i, line in enumerate(f):
                line = line.rstrip('\n\r').split(',')
                self._l.append(line)

                if len(line) == 1 and line[0] == '': #skip blank lines
                    continue
                
                if header_regex.match(line[0]):
                    section_name = header_regex.match(line[0]).group(1)
                    continue

                if not section_name:
                    raise Exception('could not parse section name in sample sheet %s' % self.sample_sheet)
                
                if section_name not in self._ss:
                    if section_name in ["Data", "Reads"]:
                        self._ss[section_name] = []
                    else:
                        self._ss[section_name] = {}
                    
                if section_name in ["Header", "Settings"]:
                    (key, value) = (line[0], line[1])
                    if key in self._ss[section_name]:
                        raise Exception('key %s observed twice in section %s' % (key, section_name))
                    self._ss[section_name][key] = i
                    
                if section_name in ["Data"]:
                    if not _data_header:
                        _data_header = line
                    else:
                        self._ss[section_name].append(Data(_data_header, _data_i, self, i))
                        _data_i += 1

                if section_name in ["Reads"]:
                    self._ss[section_name].append(line)

    def get(self, section, attribute=None):
        
        check_valid_section(section)
        
        if section in ["Data",]:
            # returns a list of "Data" items
            return self._ss[section]
        else:
            if attribute in self._ss[section]:
                i = self.ss_[section][attribute]
                return self._l[i][1]
            elif attribute is None:
                return self._ss[section]
            else:
                return None
    
    def set(self, section, attribute, new_value):
        
        check_valid_section(section)
        
        if section in ["Data",]:
            raise Exception('set method not defined for Data')
        
        if not attribute in self._ss[section]:
            raise Exception('%s not found in section %s' % (attribute, section))
        
        i = self._ss[section][attribute]
        self._l[i][1] = new_value
    
    def calculate_number_of_lanes(self):
        lanes = set()
        for d in self._ss["Data"]:
            for lane in d.get("Lane").split("+"):
                lanes.add(lane)
        lanes = sorted([int(lane) for lane in lanes])
        return lanes

    def get_sample_ids(self):
        sample_ids = set()
        for line in self._ss["Data"]:
            sample_ids.add(line.get("Sample_ID"))
        return sample_ids
    
    def get_samples_for_id(self, sample_id):
        sample_lines = []
        for line in self._ss["Data"]:
            if line.get("Sample_ID") == sample_id:
                sample_lines.append(line)
        return sample_lines

    def sample_id_to_sample_name(self, sample_id):
        '''
        Will return the first sample name associated with this sample id
        '''
        sample_lines = self.get_samples_for_id(sample_id)
        return sample_lines[0].get("Sample_Name")

    def is_paired_end(self):
        count = 0
        for line in self._ss["Reads"]:
            if len(line) > 0 and line[0] != '':
                count += 1
        if count == 2:
            return True
        if count == 1:
            return False
        raise ValueError("Reads section of samplesheet is malformed")

class Data(object):
    
    def __init__(self,
                 data_header,
                 data_i, # the ith data object (order of this data object relative to other data objects in the spreadsheet)
                 sample_sheet, # the original sample sheet object
                 sample_sheet_i # the line number in the original sample sheet
                 ):
                
        for i in ('data_i', 'sample_sheet', 'sample_sheet_i'):
            setattr(self, i, locals()[i])
        
        self.data_header_idx = {}
        for (i, v) in enumerate(data_header):
            if v == '':
                continue
            if v in self.data_header_idx:
                raise Exception('observed Data header %s twice' % v)
            self.data_header_idx[v] = i
        
    def __repr__(self):
        return str(self.sample_sheet._l[self.sample_sheet_i])

    def set(self, attribute, new_value):
        # get the line from the original sample sheet
        line = self.sample_sheet._l[self.sample_sheet_i]
        if attribute not in self.data_header_idx:
            raise AttributeError('invalid attribute for Data %s' % attribute)
        line[self.data_header_idx[attribute]] = new_value
    
    def get(self, attribute):
        
        # row_idx is another synonym for the ith data object or 'data_i'
        if attribute == "row_idx":
            return self.data_i
        
        # get the line from the original sample sheet
        line = self.sample_sheet._l[self.sample_sheet_i]
        if attribute not in self.data_header_idx:
            raise AttributeError('invalid attribute for Data %s' % attribute)
        return line[self.data_header_idx[attribute]]

    def has(self, attribute):
        
        # row_idx is another synonym for the ith data object or 'data_i'
        if attribute == "row_idx":
            return True
        if attribute not in self.data_header_idx:
            return False
        return True