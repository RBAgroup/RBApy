from rba.core import solver
import rba.cli.generate_rba_model
import rba.cli.solve_rba_model
import os
import rba
import shutil
import tempfile
import unittest
import unittest.mock


class CliTestCase(unittest.TestCase):
    def setUp(self):
        self.tmp_dirname = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmp_dirname)

    def test(self):
        # help
        with unittest.mock.patch('sys.argv', ['', '--help']):
            with self.assertRaises(SystemExit):
                rba.cli.generate_rba_model.main()

        with unittest.mock.patch('sys.argv', ['', '--help']):
            with self.assertRaises(SystemExit):
                rba.cli.solve_rba_model.main()

        # model generation
        parameters_file = os.path.join(os.path.dirname(__file__), '..', 'sample_input', 'params.in')
        with unittest.mock.patch('sys.argv', ['', parameters_file, '--model-dir', self.tmp_dirname]):
            rba.cli.generate_rba_model.main()

        self.assertTrue(os.path.join(self.tmp_dirname, 'model_file_index.in'))

        parameters_file = os.path.join(os.path.dirname(__file__), '..', 'sample_input', 'params.in')
        with unittest.mock.patch('sys.argv', ['', parameters_file, '--model-dir', self.tmp_dirname, '--verbose']):
            rba.cli.generate_rba_model.main()
