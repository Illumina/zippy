import os
import unittest
import commentjson as json
from ..params import get_params
from ..modular_runner import ModularRunner

THIS_DIR = os.path.dirname(os.path.abspath(__file__))

class DummyRunner(ModularRunner):
    '''
    Does nothing.  Designed for doing asserts.
    '''

    def get_output(self, sample):
        pass

    def workflow(self, workflowRunner):
        pass

    def get_self_parameter(self):
        return self.params.self.args

    def get_global_parameter(self):
        return self.params.dummy_path

    def get_overridden_parameter(self):
        return self.params.overridden_arg

    def get_nonlocalized_local(self):
        return self.params.self.missing_arg

    def check_hasattr_self(self, arg):
        if hasattr(self.params.self, arg):
            return True
        return False

    def check_hasattr(self, arg):
        if hasattr(self.params, arg):
            return True
        return False

class ParamsTestCase(unittest.TestCase):
    '''
    Asserts that we properly handle self.params.self parameters
    '''
    def runTest(self):
        params = get_params(os.path.join(THIS_DIR, 'dummy_params_test.json'))
        test_stage = DummyRunner(params.dummy.identifier, params, None)
        self.assertTrue(test_stage.check_hasattr_self('args'))
        self.assertFalse(test_stage.check_hasattr_self('missing_args'))
        self.assertTrue(test_stage.check_hasattr('missing_args'))
        self.assertEqual(test_stage.get_self_parameter(), 'test_args')
        self.assertEqual(test_stage.get_global_parameter(), 'path/dummy')
        self.assertEqual(test_stage.get_overridden_parameter(), 'internal')
        self.assertRaises(AttributeError, test_stage.get_nonlocalized_local)

    def tearDown(self):
        pass


class SubworkflowTestCase(unittest.TestCase):
    '''
    Asserts that we properly handle self.params.self parameters
    '''
    def runTest(self):
        params = get_params(os.path.join(THIS_DIR, 'subworkflow_params_test.json'))
        test_stage = DummyRunner(params.super_dummy.identifier, params, None)
        test_stage2 = DummyRunner(params.dummy.identifier, params, test_stage)
        self.assertTrue(test_stage.check_hasattr_self('args'))
        self.assertFalse(test_stage.check_hasattr_self('missing_args'))
        self.assertTrue(test_stage.check_hasattr('missing_args'))
        self.assertEqual(test_stage.get_self_parameter(), 'test_args')
        self.assertEqual(test_stage.get_global_parameter(), 'path/dummy2')
        self.assertEqual(test_stage.get_overridden_parameter(), 'internal')
        self.assertTrue(test_stage2.check_hasattr_self('args'))
        self.assertFalse(test_stage2.check_hasattr_self('missing_args'))
        self.assertTrue(test_stage2.check_hasattr('missing_args'))
        self.assertEqual(test_stage2.get_self_parameter(), 'test_args')
        self.assertEqual(test_stage2.get_global_parameter(), 'path/dummy2')
        self.assertEqual(test_stage2.get_overridden_parameter(), 'internal')
        self.assertRaises(AttributeError, test_stage.get_nonlocalized_local)

    def tearDown(self):
        pass

class EnvironmentTestCase(unittest.TestCase):
    '''
    Asserts that we properly handle self.params.self parameters
    '''
    def runTest(self):
        params = get_params(os.path.join(THIS_DIR, 'environment_params_test.json'))
        test_stage = DummyRunner(params.super_dummy.identifier, params, None)
        self.assertTrue(test_stage.check_hasattr_self('args'))
        self.assertFalse(test_stage.check_hasattr_self('missing_args'))
        self.assertTrue(test_stage.check_hasattr('missing_args'))
        self.assertEqual(test_stage.get_self_parameter(), 'test_args')
        self.assertEqual(test_stage.get_global_parameter(), 'path/dummy2')
        self.assertEqual(test_stage.get_overridden_parameter(), 'internal')
        self.assertRaises(AttributeError, test_stage.get_nonlocalized_local)

    def tearDown(self):
        pass

if __name__ == '__main__':
    unittest.main()