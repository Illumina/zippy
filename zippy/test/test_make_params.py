import os
import unittest
from ..make_params import get_argparser, JsonBuilderMain, get_identifier, case_insensitive_list_match
import commentjson as json

def compare_json_structures(a_path,b_path):
    with open(a_path, 'r') as a_handle:
        a = json.load(a_handle)
    with open(b_path, 'r') as b_handle:
        b = json.load(b_handle)
    return recursive_compare_keys(a,b) and recursive_compare_keys(b,a)

def recursive_compare_keys(a,b):
    for k in a.keys():
        if k in b:
            if isinstance(a[k], dict):
                if not recursive_compare_keys(a[k], b[k]):
                    return False
        else:
            print 'missing key {}'.format(k)
            return False
    return True

class ImportsTestCase(unittest.TestCase):
    '''
    Asserts that we pass through the imports
    '''
    def runTest(self):
        parser = get_argparser()
        params = parser.parse_args(['zippy/test/import_test.proto', 'test_output_12344.json', '--defaults', 'zippy/test/defaults.json'])
        jsb = JsonBuilderMain(params)
        jsb.make_params_file()
        self.assertTrue(compare_json_structures('test_output_12344.json','zippy/test/import_test.json'))
    
    def tearDown(self):
        try:
            pass
            os.remove('test_output_12344.json')
        except OSError:
            pass

class MacsTestCase(unittest.TestCase):
    '''
    Asserts that we produce a json file with the proper parameters exposed given the workflow
    '''
    def runTest(self):
        parser = get_argparser()
        params = parser.parse_args(['zippy/test/merge_and_macs.proto', 'test_output_12345.json', '--defaults', 'zippy/test/defaults.json'])
        jsb = JsonBuilderMain(params)
        jsb.make_params_file()
        self.assertTrue(compare_json_structures('test_output_12345.json','zippy/test/merge_and_macs.json'))
    
    def tearDown(self):
        try:
            pass
            os.remove('test_output_12345.json')
        except OSError:
            pass

class RSEMTestCase(unittest.TestCase):
    '''
    Asserts that we produce a json file with the proper parameters exposed given the workflow
    '''
    def runTest(self):
        parser = get_argparser()
        params = parser.parse_args(['zippy/test/rsem.proto', 'test_output_12346.json', '--defaults', 'zippy/test/defaults.json'])
        jsb = JsonBuilderMain(params)
        jsb.make_params_file()
        self.assertTrue(compare_json_structures('test_output_12346.json','zippy/test/rsem.json'))
    
    def tearDown(self):
        try:
            pass
            os.remove('test_output_12346.json')
        except OSError:
            pass

class GetIdentifiersTestCase(unittest.TestCase):
    def runTest(self):
        current_identifiers = set()
        get_identifier('test', current_identifiers)
        self.assertEqual(current_identifiers, set(['test']))
        get_identifier('test', current_identifiers)
        self.assertEqual(current_identifiers, set(['test', 'test_1']))
        get_identifier('test', current_identifiers)
        self.assertEqual(current_identifiers, set(['test', 'test_1', 'test_2']))
        print current_identifiers

class CaseInsensitiveListMatchTestCase(unittest.TestCase):
    def runTest(self):
        parser = get_argparser()
        params = parser.parse_args(['zippy/test/rsem.proto', 'test_output_12348.json'])
        jsb = JsonBuilderMain(params)
        l = ['hello', 'HELLO', 'HeLlO', 'hello1', 'ello', 'h', 'e', 'l' ,'l', 'o', 'heLLO']
        result = case_insensitive_list_match('hello', l)
        self.assertEqual(result, 'hello')
        l = ['HELLO', 'HeLlO', 'hello1', 'ello', 'h', 'e', 'l' ,'l', 'o', 'heLLO', 'hello']
        result = case_insensitive_list_match('hello', l)
        self.assertEqual(result, 'HELLO')
        l = ['ELLO', 'eLlO', 'helo1', 'ello', 'h', 'e', 'l' ,'l', 'o', 'heL', 'hllo']
        result = case_insensitive_list_match('hello', l)
        self.assertEqual(result, None)

class DefaultParameterTest(unittest.TestCase):
    def runTest(self):
        for behavior in ['include', 'warn', 'ignore']:
            parser = get_argparser()
            params = parser.parse_args(['zippy/test/defaults_test.proto', 'test_output_12347.json', '--defaults', 'zippy/test/test_defaults.json', '--default_behavior', behavior])
            jsb = JsonBuilderMain(params)
            jsb.make_params_file()
            self.assertTrue(compare_json_structures('test_output_12347.json','zippy/test/defaults_{}.json'.format(behavior)))

    def tearDown(self):
        try:
            pass
            os.remove('test_output_12347.json')
        except OSError:
            pass

if __name__ == '__main__':
    unittest.main()