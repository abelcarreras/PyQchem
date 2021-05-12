from pyqchem.qchem_core import get_output_from_qchem, create_qchem_input, redefine_calculation_data_filename
from pyqchem.parsers.parser_rasci_basic import basic_rasci
from pyqchem.structure import Structure
from pyqchem.test import standardize_dictionary

import yaml
import unittest

ctrl_print = False  # set to true to generate reference files
do_alpha_beta = False # explicitly defining alpha/beta number of electrons?


redefine_calculation_data_filename('test_data.db')


class Eth00(unittest.TestCase):

    def setUp(self):
        self.assertDictEqual.__self__.maxDiff = None

        # generate molecule
        self.molecule = Structure(coordinates=[
                                  [0.665000, 0.000000, 0.000000],
                                  [-0.665000, 0.000000, 0.000000],
                                  [-1.230407, 0.915473, 0.000000],
                                  [-1.230407, -0.915473, 0.000000],
                                  [1.230407, 0.915473, 0.000000],
                                  [1.230407, -0.915473, 0.000000]],
                  symbols=['C', 'C', 'H', 'H', 'H', 'H'],
                  charge=0,
                  multiplicity=1)

        self.rem = {'jobtype': 'sp',
                    'exchange': 'hf',
                    'basis': '6-31G(d,p)',
                    'thresh': 14,
                    'scf_convergence': 8,
                    'max_scf_cycles': 150,
                    'max_cis_cycles': 150,
                    'unrestricted': False,
                    # RASCI
                    'correlation': 'rasci',
                    'ras_act': 2,
                    'ras_elec': 2,
                    'ras_occ': 7,
                    'ras_spin_mult': 0,
                    'ras_roots': 7,
                    'ras_do_hole': True,
                    'ras_do_part': True,
                    'ras_sts_tm': True,
                    'ras_natorb': True,
                    # rasci sr-dft
                    # 'ras_srdft': False,
                    # 'ras_omega': 400,
                    # 'ras_srdft_cor': 'srpw92',
                    # 'ras_srdft_exc': 'srpbe',
                    # 'ras_srdft_damp': 0.5,
                    #  rasci print level
                    'ras_print': 2,
                    # Frozen
                    'n_frozen_core': 0,
                    'n_frozen_virt': 0,
                    'set_iter': 60}

    def test_eth_00_ras22(self):

        rasci = dict(self.rem)
        if do_alpha_beta:
            rasci.update({'ras_elec_alpha': 1,
                          'ras_elec_beta': 1})

        # create qchem input
        txt_input = create_qchem_input(self.molecule, **rasci)

        # calculate and parse qchem output
        output = get_output_from_qchem(txt_input, processors=4)
        print(output)
        data = basic_rasci(output)

        filename = self.__class__.__name__ + '_ras22.yaml'

        # create reference file
        if ctrl_print:
            with open(filename, 'w') as outfile:
                yaml.dump(data, outfile, default_flow_style=False, allow_unicode=True)

        with open(filename, 'r') as stream:
            data_loaded = yaml.safe_load(stream)

        data = standardize_dictionary(data)

        #print(data)

        print(data_loaded)
        data_loaded = standardize_dictionary(data_loaded)

        self.assertDictEqual(data, data_loaded)

    def test_eth_00_ras44(self):

        rasci = dict(self.rem)
        rasci.update({'ras_act': 4,
                      'ras_elec': 4,
                      'ras_occ': 6})
        if do_alpha_beta:
            rasci.update({'ras_elec_alpha': 2,
                          'ras_elec_beta': 2})


        # create qchem input
        txt_input = create_qchem_input(self.molecule, **rasci)

        # calculate and parse qchem output
        output = get_output_from_qchem(txt_input, processors=4)
        print(output)
        data = basic_rasci(output)

        print(data)

        filename = self.__class__.__name__ + '_ras44.yaml'

        # creat reference file
        if ctrl_print:
            with open(filename, 'w') as outfile:
                yaml.dump(data, outfile, default_flow_style=False, allow_unicode=True)

        with open(filename, 'r') as stream:
            data_loaded = yaml.safe_load(stream)

        data = standardize_dictionary(data)

        print(data_loaded)
        data_loaded = standardize_dictionary(data_loaded)

        self.assertDictEqual(data, data_loaded)

    def test_eth_00_ras66(self):

        rasci = dict(self.rem)
        rasci.update({'ras_act': 6,
                      'ras_elec': 6,
                      'ras_occ': 5})
        if do_alpha_beta:
            rasci.update({'ras_elec_alpha': 3,
                          'ras_elec_beta': 3})

        # create qchem input
        txt_input = create_qchem_input(self.molecule, **rasci)

        # calculate and parse qchem output
        output = get_output_from_qchem(txt_input, processors=4)
        print(output)
        data = basic_rasci(output)

        print(data)

        filename = self.__class__.__name__ + '_ras66.yaml'

        # creat reference file
        if ctrl_print:
            with open(filename, 'w') as outfile:
                yaml.dump(data, outfile, default_flow_style=False, allow_unicode=True)

        with open(filename, 'r') as stream:
            data_loaded = yaml.safe_load(stream)

        data = standardize_dictionary(data)

        print(data_loaded)
        data_loaded = standardize_dictionary(data_loaded)

        self.assertDictEqual(data, data_loaded)


class Eth90(unittest.TestCase):

    def setUp(self):
        self.assertDictEqual.__self__.maxDiff = None

        # generate molecule
        self.molecule = Structure(coordinates=[
            [0.665000, 0.000000, 0.000000],
            [-0.665000, 0.000000, 0.000000],
            [-1.230407, 0.915473, 0.000000],
            [-1.230407, -0.915473, 0.000000],
            [1.230407, 0.0, 0.915473],
            [1.230407, 0.0, -0.915473]],
            symbols=['C', 'C', 'H', 'H', 'H', 'H'],
            charge=0,
            multiplicity=3)

        self.rem = {'jobtype': 'sp',
                    'exchange': 'hf',
                    'basis': '6-31G(d)',
                    'thresh': 14,
                    'scf_convergence': 8,
                    'max_scf_cycles': 150,
                    'max_cis_cycles': 150,
                    'unrestricted': False,
                    # RASCI
                    'correlation': 'rasci',
                    'ras_act': 2,
                    'ras_elec': 2,
                    'ras_occ': 7,
                    'ras_spin_mult': 0,
                    'ras_roots': 7,
                    'ras_do_hole': True,
                    'ras_do_part': True,
                    'ras_sts_tm': True,
                    'ras_natorb': True,
                    # rasci sr-dft
                    # 'ras_srdft': False,
                    # 'ras_omega': 400,
                    # 'ras_srdft_cor': 'srpw92',
                    # 'ras_srdft_exc': 'srpbe',
                    # 'ras_srdft_damp': 0.5,
                    # rasci print level
                    'ras_print': 2,
                    # Frozen
                    'n_frozen_core': 0,
                    'n_frozen_virt': 0,
                    'set_iter': 60}

    def test_eth_90_ras22(self):

        rasci = dict(self.rem)
        if do_alpha_beta:
            rasci.update({'ras_elec_alpha': 1,
                          'ras_elec_beta': 1})

        # create qchem input
        txt_input = create_qchem_input(self.molecule, **rasci)

        # calculate and parse qchem output
        output = get_output_from_qchem(txt_input, processors=4)
        print(output)
        data = basic_rasci(output)

        print(data)

        filename = self.__class__.__name__ + '_ras22.yaml'

        # create reference file
        if ctrl_print:
            with open(filename, 'w') as outfile:
                yaml.dump(data, outfile, default_flow_style=False, allow_unicode=True)

        with open(filename, 'r') as stream:
            data_loaded = yaml.safe_load(stream)

        data = standardize_dictionary(data)

        print(data_loaded)
        data_loaded = standardize_dictionary(data_loaded)

        self.assertDictEqual(data, data_loaded)

    def test_eth_90_ras44(self):

        rasci = dict(self.rem)
        rasci.update({'ras_act': 4,
                      'ras_elec': 4,
                      'ras_occ': 6})
        if do_alpha_beta:
            rasci.update({'ras_elec_alpha': 2,
                          'ras_elec_beta': 2})

        # create qchem input
        txt_input = create_qchem_input(self.molecule, **rasci)

        # calculate and parse qchem output
        output = get_output_from_qchem(txt_input, processors=4)
        print(output)
        data = basic_rasci(output)

        print(data)

        filename = self.__class__.__name__ + '_ras44.yaml'

        # creat reference file
        if ctrl_print:
            with open(filename, 'w') as outfile:
                yaml.dump(data, outfile, default_flow_style=False, allow_unicode=True)

        with open(filename, 'r') as stream:
            data_loaded = yaml.safe_load(stream)

        data = standardize_dictionary(data, decimal=3)

        print(data_loaded)
        data_loaded = standardize_dictionary(data_loaded, decimal=3)

        self.assertDictEqual(data, data_loaded)

    def test_eth_90_ras66(self):

        rasci = dict(self.rem)
        rasci.update({'ras_act': 6,
                      'ras_elec': 6,
                      'ras_occ': 5})
        if do_alpha_beta:
            rasci.update({'ras_elec_alpha': 3,
                          'ras_elec_beta': 3})

        # create qchem input
        txt_input = create_qchem_input(self.molecule, **rasci)

        # calculate and parse qchem output
        output = get_output_from_qchem(txt_input, processors=4)
        print(output)
        data = basic_rasci(output)

        print(data)

        filename = self.__class__.__name__ + '_ras66.yaml'

        # create reference file
        if ctrl_print:
            with open(filename, 'w') as outfile:
                yaml.dump(data, outfile, default_flow_style=False, allow_unicode=True)

        with open(filename, 'r') as stream:
            data_loaded = yaml.safe_load(stream)

        data = standardize_dictionary(data, decimal=3)

        print(data_loaded)
        data_loaded = standardize_dictionary(data_loaded, decimal=3)

        self.assertDictEqual(data, data_loaded)


class EthDist(unittest.TestCase):

    def setUp(self):
        self.assertDictEqual.__self__.maxDiff = None

        # generate molecule
        self.molecule = Structure(coordinates=[
            [0.865000, 0.200000, 0.800000],
            [-0.665000, -0.300000, 0.500000],
            [-1.130407, 0.715473, 0.220000],
            [-1.430407, -0.915473, 0.300000],
            [1.230407, 0.815473, -0.4],
            [1.000407,-0.888888, 0.45]],
            symbols=['C', 'C', 'H', 'H', 'H', 'H'],
            charge=0,
            multiplicity=1)

        self.rem = {'jobtype': 'sp',
                    'exchange': 'hf',
                    'basis': '6-31G(d,p)',
                    'thresh': 14,
                    'scf_convergence': 8,
                    'max_scf_cycles': 150,
                    'max_cis_cycles': 150,
                    'unrestricted': False,
                    # RASCI
                    'correlation': 'rasci',
                    'ras_act': 2,
                    'ras_elec': 2,
                    'ras_occ': 7,
                    'ras_spin_mult': 0,
                    'ras_roots': 7,
                    'ras_do_hole': True,
                    'ras_do_part': True,
                    'ras_sts_tm': True,
                    'ras_natorb': True,
                    # rasci sr-dft
                    # 'ras_srdft': False,
                    # 'ras_omega': 400,
                    # 'ras_srdft_cor': 'srpw92',
                    # 'ras_srdft_exc': 'srpbe',
                    # 'ras_srdft_damp': 0.5,
                    # rasci print level
                    'ras_print': 2,
                    # Frozen
                    'n_frozen_core': 0,
                    'n_frozen_virt': 0,
                    'set_iter': 60}

    def test_eth_dist_ras22(self):

        rasci = dict(self.rem)
        if do_alpha_beta:
            rasci.update({'ras_elec_alpha': 1,
                          'ras_elec_beta': 1})

        # create qchem input
        txt_input = create_qchem_input(self.molecule, **rasci)

        # calculate and parse qchem output
        output = get_output_from_qchem(txt_input, processors=4)
        print(output)
        data = basic_rasci(output)

        print(data)

        filename = self.__class__.__name__ + '_ras22.yaml'

        # creat reference file
        if ctrl_print:
            with open(filename, 'w') as outfile:
                yaml.dump(data, outfile, default_flow_style=False, allow_unicode=True)

        with open(filename, 'r') as stream:
            data_loaded = yaml.safe_load(stream)

        data = standardize_dictionary(data)

        print(data_loaded)
        data_loaded = standardize_dictionary(data_loaded)

        self.assertDictEqual(data, data_loaded)

    def test_eth_dist_ras44(self):

        rasci = dict(self.rem)
        rasci.update({'ras_act': 4,
                      'ras_elec': 4,
                      'ras_occ': 6})
        if do_alpha_beta:
            rasci.update({'ras_elec_alpha': 2,
                          'ras_elec_beta': 2})

        # create qchem input
        txt_input = create_qchem_input(self.molecule, **rasci)

        # calculate and parse qchem output
        output = get_output_from_qchem(txt_input, processors=4)
        print(output)
        data = basic_rasci(output)

        print(data)

        filename = self.__class__.__name__ + '_ras44.yaml'

        # creat reference file
        if ctrl_print:
            with open(filename, 'w') as outfile:
                yaml.dump(data, outfile, default_flow_style=False, allow_unicode=True)

        with open(filename, 'r') as stream:
            data_loaded = yaml.safe_load(stream)

        data = standardize_dictionary(data)

        print(data_loaded)
        data_loaded = standardize_dictionary(data_loaded)

        self.assertDictEqual(data, data_loaded)

    def test_eth_dist_ras66(self):

        rasci = dict(self.rem)
        rasci.update({'ras_act': 6,
                      'ras_elec': 6,
                      'ras_occ': 5})
        if do_alpha_beta:
            rasci.update({'ras_elec_alpha': 3,
                          'ras_elec_beta': 3})

        # create qchem input
        txt_input = create_qchem_input(self.molecule, **rasci)

        # calculate and parse qchem output
        output = get_output_from_qchem(txt_input, processors=4)
        print(output)
        data = basic_rasci(output)

        print(data)

        filename = self.__class__.__name__ + '_ras66.yaml'

        # creat reference file
        if ctrl_print:
            with open(filename, 'w') as outfile:
                yaml.dump(data, outfile, default_flow_style=False, allow_unicode=True)

        with open(filename, 'r') as stream:
            data_loaded = yaml.safe_load(stream)

        data = standardize_dictionary(data)

        print(data_loaded)
        data_loaded = standardize_dictionary(data_loaded)

        self.assertDictEqual(data, data_loaded)
