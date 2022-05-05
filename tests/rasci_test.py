from pyqchem.qchem_core import get_output_from_qchem, create_qchem_input, redefine_calculation_data_filename, QchemInput
from pyqchem.parsers.parser_rasci import parser_rasci
from pyqchem.parsers.parser_optimization import basic_optimization
from pyqchem.structure import Structure
from pyqchem.test import standardize_dictionary
import yaml
import unittest
import os, sys


redefine_calculation_data_filename('test_data.db')

#if 'USER' in os.environ and os.environ['USER'] == 'travis':
recalculate = False
remake_tests = False
dir_path = os.path.dirname(os.path.realpath(__file__))

cwd = os.getcwd()

class HydrogenTest(unittest.TestCase):

    def setUp(self):
        self.assertDictEqual.__self__.maxDiff = None

        # generate molecule
        self.molecule = Structure(coordinates=[[0.0, 0.0, 0.0],
                                               [0.0, 0.0, 1.5]],
                                  symbols=['H', 'H'],
                                  charge=0,
                                  multiplicity=1)

        print('work_dir:', cwd)


    def test_srdft(self):

        # create qchem input
        qc_input = create_qchem_input(self.molecule,
                                       jobtype='sp',
                                       exchange='hf',
                                       correlation='rasci',
                                       basis='6-31G(d,p)',
                                       ras_act=2,
                                       ras_elec=2,
                                       ras_spin_mult=0,
                                       ras_roots=2,
                                       ras_do_hole=True,
                                       ras_sts_tm=True,
                                       # rasci sr-dft
                                       ras_omega=400,
                                       ras_srdft_cor='srpw92',
                                       ras_srdft_exc='srpbe',
                                       ras_natorb=False,
                                       set_iter=30,
                                       ras_srdft_damp=0.5)


        from pyqchem.cache import SqlCache as CacheSystem

        cache = CacheSystem()
        cache.list_database()
        output = cache.retrieve_calculation_data(qc_input, 'fullout')
        print(output)

        # calculate and parse qchem output
        output = get_output_from_qchem(qc_input,
                                       processors=4,
                                       force_recalculation=recalculate,
                                       store_full_output=True)
        print(output)
        data = parser_rasci(output)

        filename = dir_path + '/' + self.__class__.__name__ + '_srdft.yaml'

        if remake_tests:
            with open(filename, 'w') as outfile:
                yaml.dump(data, outfile, default_flow_style=False, allow_unicode=True)

        data = standardize_dictionary(data, decimal=2)

        with open(filename, 'r') as stream:
            data_loaded = yaml.load(stream, Loader=yaml.Loader)

        print(data_loaded)
        data_loaded = standardize_dictionary(data_loaded, decimal=2)

        self.assertDictEqual(data, data_loaded)


    def _test_rasci(self):
        # create qchem input
        txt_input = create_qchem_input(self.molecule,
                                       jobtype='sp',
                                       exchange='hf',
                                       correlation='rasci',
                                       basis='sto-3g',
                                       ras_act=2,
                                       ras_elec=2,
                                       ras_spin_mult=0,
                                       ras_roots=2,
                                       ras_print=5,
                                       ras_do_hole=True,
                                       ras_sts_tm=True)

        # calculate and parse qchem output
        output = get_output_from_qchem(txt_input,
                                       processors=4,
                                       force_recalculation=recalculate,
                                       store_full_output=True)
        print(output)
        data = parser_rasci(output)

        filename = dir_path + '/' + self.__class__.__name__ + '_rasci.yaml'

        if remake_tests:
            with open(filename, 'w') as outfile:
                yaml.dump(data, outfile, default_flow_style=False, allow_unicode=True)

        data = standardize_dictionary(data, decimal=2)

        with open(filename, 'r') as stream:
            data_loaded = yaml.load(stream, Loader=yaml.Loader)

        print(data_loaded)
        data_loaded = standardize_dictionary(data_loaded, decimal=2)

        self.assertDictEqual(data, data_loaded)


class WaterTest(unittest.TestCase):

    def setUp(self):
        self.assertDictEqual.__self__.maxDiff = None

        # generate molecule
        molecule = Structure(coordinates=[[0.0, 1.0, 0.0],
                                          [0.0, 0.0, 1.0],
                                          [0.0, 0.0, -1.0]],
                             symbols=['O', 'H', 'H'],
                             charge=0,
                             multiplicity=1)

        # optimization
        txt_input = create_qchem_input(molecule,
                                       jobtype='opt',
                                       exchange='hf',
                                       basis='sto-3g',
                                       geom_opt_tol_gradient=300,
                                       geom_opt_tol_energy=100,
                                       geom_opt_coords=-1,
                                       geom_opt_tol_displacement=1200)

        parsed_data = get_output_from_qchem(txt_input,
                                            processors=4,
                                            parser=basic_optimization,
                                            force_recalculation=recalculate,
                                            store_full_output=True)

        self.molecule = parsed_data['optimized_molecule']

    def _test_rasci(self):
        # create qchem input
        txt_input = create_qchem_input(self.molecule,
                                       jobtype='sp',
                                       exchange='hf',
                                       correlation='rasci',
                                       basis='sto-3g',
                                       ras_act=4,
                                       ras_elec=4,
                                       ras_spin_mult=1,
                                       ras_roots=6,
                                       ras_print=5,
                                       ras_do_hole=True,
                                       ras_sts_tm=True)

        # calculate and parse qchem output
        output = get_output_from_qchem(txt_input,
                                       processors=4,
                                       force_recalculation=recalculate,
                                       store_full_output=True)
        print(output)
        data = parser_rasci(output)

        print(data)

        filename = dir_path + '/' + self.__class__.__name__ + '_rasci.yaml'

        if remake_tests:
            with open(filename, 'w') as outfile:
                yaml.dump(data, outfile, default_flow_style=False, allow_unicode=True)

        data = standardize_dictionary(data, decimal=2)

        with open(filename, 'r') as stream:
            data_loaded = yaml.load(stream, Loader=yaml.Loader)

        print(data_loaded)
        data_loaded = standardize_dictionary(data_loaded, decimal=2)

        self.assertDictEqual(data, data_loaded)


    def _test_srdft(self):

        # create qchem input
        txt_input = create_qchem_input(self.molecule,
                                       jobtype='sp',
                                       exchange='hf',
                                       correlation='rasci',
                                       basis='6-31G(d,p)',
                                       ras_act=4,
                                       ras_elec=4,
                                       ras_spin_mult=0,
                                       ras_roots=6,
                                       ras_do_hole=True,
                                       ras_sts_tm=True,
                                       # rasci sr-dft
                                       ras_omega=300,
                                       ras_srdft_cor='srpbe',
                                       ras_srdft_exc='srlsda',
                                       ras_natorb=False,
                                       set_iter=30,
                                       ras_srdft_damp=0.5)

        print(txt_input.get_txt())

        # calculate and parse qchem output
        output = get_output_from_qchem(txt_input,
                                       processors=4,
                                       force_recalculation=recalculate,
                                       store_full_output=True)

        print(output)
        data = parser_rasci(output)

        print(data)

        filename = dir_path + '/' + self.__class__.__name__ + '_srdft.yaml'

        if remake_tests:
            with open(filename, 'w') as outfile:
                yaml.dump(data, outfile, default_flow_style=False, allow_unicode=True)

        data = standardize_dictionary(data)

        with open(filename, 'r') as stream:
            data_loaded = yaml.load(stream, Loader=yaml.Loader)

        print(data_loaded)
        data_loaded = standardize_dictionary(data_loaded)

        self.assertDictEqual(data, data_loaded)
