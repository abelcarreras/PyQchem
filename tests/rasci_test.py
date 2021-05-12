from pyqchem.qchem_core import get_output_from_qchem, create_qchem_input, redefine_calculation_data_filename
from pyqchem.parsers.parser_rasci_basic import basic_rasci
from pyqchem.parsers.parser_optimization import basic_optimization
from pyqchem.structure import Structure
from pyqchem.test import standardize_dictionary
import yaml
import unittest
import os, sys


redefine_calculation_data_filename('test_data.db')

if 'USER' in os.environ and os.environ['USER'] == 'travis':
    recalculate = False
else:
    recalculate = True


class HydrogenTest(unittest.TestCase):

    def setUp(self):
        self.assertDictEqual.__self__.maxDiff = None

        # generate molecule
        self.molecule = Structure(coordinates=[[0.0, 0.0, 0.0],
                                               [0.0, 0.0, 1.5]],
                                  symbols=['H', 'H'],
                                  charge=0,
                                  multiplicity=1)

    def test_srdft(self):

        # create qchem input
        txt_input = create_qchem_input(self.molecule,
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

        # calculate and parse qchem output
        output = get_output_from_qchem(txt_input,
                                       processors=4,
                                       force_recalculation=recalculate,
                                       store_full_output=True)
        print(output)
        data = basic_rasci(output)

        filename = self.__class__.__name__ + '_srdft.yaml'

        #with open(filename, 'w') as outfile:
        #      yaml.dump(data, outfile, default_flow_style=False, allow_unicode=True)

        data = standardize_dictionary(data)

        with open(filename, 'r') as stream:
            data_loaded = yaml.safe_load(stream)

        print(data_loaded)
        data_loaded = standardize_dictionary(data_loaded)

        self.assertDictEqual(data, data_loaded)


    def test_rasci(self):
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
        data = basic_rasci(output)

        filename = self.__class__.__name__ + '_rasci.yaml'

        # with open(filename, 'w') as outfile:
        #      yaml.dump(data, outfile, default_flow_style=False, allow_unicode=True)

        data = standardize_dictionary(data)

        with open(filename, 'r') as stream:
            data_loaded = yaml.safe_load(stream)

        print(data_loaded)
        data_loaded = standardize_dictionary(data_loaded)

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

    def test_rasci(self):
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
        data = basic_rasci(output)

        print(data)

        filename = self.__class__.__name__ + '_rasci.yaml'

        # with open(filename, 'w') as outfile:
        #     yaml.dump(data, outfile, default_flow_style=False, allow_unicode=True)

        data = standardize_dictionary(data)

        with open(filename, 'r') as stream:
            data_loaded = yaml.safe_load(stream)

        print(data_loaded)
        data_loaded = standardize_dictionary(data_loaded)

        self.assertDictEqual(data, data_loaded)


    def test_srdft(self):

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
        data = basic_rasci(output)

        print(data)

        filename = self.__class__.__name__ + '_srdft.yaml'

        # with open(filename, 'w') as outfile:
        #     yaml.dump(data, outfile, default_flow_style=False, allow_unicode=True)

        data = standardize_dictionary(data)

        with open(filename, 'r') as stream:
            data_loaded = yaml.safe_load(stream)

        print(data_loaded)
        data_loaded = standardize_dictionary(data_loaded)

        self.assertDictEqual(data, data_loaded)
