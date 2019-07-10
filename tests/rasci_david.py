from pyqchem.qchem_core import get_output_from_qchem, create_qchem_input
from pyqchem.parsers.parser_rasci_basic import basic_rasci
from pyqchem.parsers.parser_optimization import basic_optimization
from pyqchem.structure import Structure
import yaml
import unittest

ctrl_print = True  # set to true to generate reference files

def trunc_dictionary_list(w, decimal=6):

    def round_float(num):
        if isinstance(num, float):
            return int(num * 10**decimal)
        else:
            return num

    def iterlist(d):
        lf = []
        for i, v in enumerate(d):
            if isinstance(v, dict):
                iterdict(v)
            elif isinstance(v, list):
                iterlist(v)
            else:
                d[i] = round_float(v)

    def iterdict(d):
        lf = []
        for k,v in d.items():
            if isinstance(v, dict):
                iterdict(v)
            elif isinstance(v, list):
                iterlist(v)
            else:
                d[k] = round_float(v)

    if isinstance(w, dict):
        iterdict(w)
    elif isinstance(w, list):
        iterlist(w)
    else:
        round_float(w)


def modify_dictionary(dic_data):
    for state in dic_data['excited states rasci']:
        for amp in state['amplitudes']:
            amp['amplitude'] = abs(amp['amplitude'])


    trunc_dictionary_list(dic_data)

    return dic_data


def assertDeepAlmostEqual(test_case, expected, actual, *args, **kwargs):
    """
    Assert that two complex structures have almost equal contents.

    Compares lists, dicts and tuples recursively. Checks numeric values
    using test_case's :py:meth:`unittest.TestCase.assertAlmostEqual` and
    checks all other values with :py:meth:`unittest.TestCase.assertEqual`.
    Accepts additional positional and keyword arguments and pass those
    intact to assertAlmostEqual() (that's how you specify comparison
    precision).

    :param test_case: TestCase object on which we can call all of the basic
    'assert' methods.
    :type test_case: :py:class:`unittest.TestCase` object
    """

    import numpy
    from numpy import long
    is_root = not '__trace' in kwargs
    trace = kwargs.pop('__trace', 'ROOT')
    try:
        if isinstance(expected, (int, float, long, complex)):
            test_case.assertAlmostEqual(expected, actual, *args, **kwargs)
        elif isinstance(expected, (list, tuple, numpy.ndarray)):
            test_case.assertEqual(len(expected), len(actual))
            for index in range(len(expected)):
                v1, v2 = expected[index], actual[index]
                assertDeepAlmostEqual(test_case, v1, v2,
                                      __trace=repr(index), *args, **kwargs)
        elif isinstance(expected, dict):
            test_case.assertEqual(set(expected), set(actual))
            for key in expected:
                assertDeepAlmostEqual(test_case, expected[key], actual[key],
                                      __trace=repr(key), *args, **kwargs)
        else:
            test_case.assertEqual(expected, actual)
    except AssertionError as exc:
        exc.__dict__.setdefault('traces', []).append(trace)
        if is_root:
            trace = ' -> '.join(reversed(exc.traces))
            exc = AssertionError("%s\nTRACE: %s" % (exc.message, trace))
        raise exc


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
                  atomic_elements=['C', 'C', 'H', 'H', 'H', 'H'],
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
                    'ras_elec_alpha': 1,
                    'ras_elec_beta': 1,
                    'ras_occ': 7,
                    'ras_spin_mult': 0,
                    'ras_roots': 7,
                    'ras_do_hole': True,
                    'ras_do_part': True,
                    'ras_sts_tm': True,
                    'ras_natorb': True,
                    # rasci sr-dft
                    'ras_srdft': False,
                    'ras_omega': 400,
                    'ras_srdft_cor': 'srpw92',
                    'ras_srdft_exc': 'srpbe',
                    'ras_srdft_damp': 0.5,
                    # rasci print level
                    'ras_print': 2,
                    # Frozen
                    'n_frozen_core': 0,
                    'n_frozen_virt': 0}


    def test_eth_00_ras22(self):

        rasci = dict(self.rem)

        # create qchem input
        txt_input = create_qchem_input(self.molecule, **rasci)

        # calculate and parse qchem output
        output = get_output_from_qchem(txt_input, processors=4)
        print(output)
        data = basic_rasci(output)
        data = modify_dictionary(data)

        print(data)

        filename = self.__class__.__name__ + '_ras22.yaml'

        # create reference file
        if ctrl_print:
            with open(filename, 'w') as outfile:
                yaml.dump(data, outfile, default_flow_style=False, allow_unicode=True)

        with open(filename, 'r') as stream:
            data_loaded = yaml.safe_load(stream)

        print(data_loaded)
        trunc_dictionary_list(data_loaded)

        self.assertDictEqual(data, data_loaded)


    def test_eth_00_ras44(self):

        rasci = dict(self.rem)
        rasci.update({'ras_act': 4,
                      'ras_elec': 4,
                      'ras_elec_alpha': 2,
                      'ras_elec_beta': 2,
                      'ras_occ': 6})

        # create qchem input
        txt_input = create_qchem_input(self.molecule, **rasci)

        # calculate and parse qchem output
        output = get_output_from_qchem(txt_input, processors=4)
        print(output)
        data = basic_rasci(output)
        data = modify_dictionary(data)

        print(data)

        filename = self.__class__.__name__ + '_ras44.yaml'

        # creat reference file
        if ctrl_print:
            with open(filename, 'w') as outfile:
                yaml.dump(data, outfile, default_flow_style=False, allow_unicode=True)

        with open(filename, 'r') as stream:
            data_loaded = yaml.safe_load(stream)

        print(data_loaded)
        trunc_dictionary_list(data_loaded)

        self.assertDictEqual(data, data_loaded)

    def test_eth_00_ras66(self):

        rasci = dict(self.rem)
        rasci.update({'ras_act': 6,
                      'ras_elec': 6,
                      'ras_elec_alpha': 3,
                      'ras_elec_beta': 3,
                      'ras_occ': 5})

        # create qchem input
        txt_input = create_qchem_input(self.molecule, **rasci)

        # calculate and parse qchem output
        output = get_output_from_qchem(txt_input, processors=4)
        print(output)
        data = basic_rasci(output)
        data = modify_dictionary(data)

        print(data)

        filename = self.__class__.__name__ + '_ras66.yaml'

        # creat reference file
        if ctrl_print:
            with open(filename, 'w') as outfile:
                yaml.dump(data, outfile, default_flow_style=False, allow_unicode=True)

        with open(filename, 'r') as stream:
            data_loaded = yaml.safe_load(stream)

        print(data_loaded)
        trunc_dictionary_list(data_loaded)

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
            atomic_elements=['C', 'C', 'H', 'H', 'H', 'H'],
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
                    'ras_elec_alpha': 1,
                    'ras_elec_beta': 1,
                    'ras_occ': 7,
                    'ras_spin_mult': 0,
                    'ras_roots': 7,
                    'ras_do_hole': True,
                    'ras_do_part': True,
                    'ras_sts_tm': True,
                    'ras_natorb': True,
                    # rasci sr-dft
                    'ras_srdft': False,
                    'ras_omega': 400,
                    'ras_srdft_cor': 'srpw92',
                    'ras_srdft_exc': 'srpbe',
                    'ras_srdft_damp': 0.5,
                    # rasci print level
                    'ras_print': 2,
                    # Frozen
                    'n_frozen_core': 0,
                    'n_frozen_virt': 0}


    def test_eth_90_ras22(self):

        rasci = dict(self.rem)

        # create qchem input
        txt_input = create_qchem_input(self.molecule, **rasci)

        # calculate and parse qchem output
        output = get_output_from_qchem(txt_input, processors=4)
        print(output)
        data = basic_rasci(output)
        data = modify_dictionary(data)

        print(data)

        filename = self.__class__.__name__ + '_ras22.yaml'

        # create reference file
        if ctrl_print:
            with open(filename, 'w') as outfile:
                yaml.dump(data, outfile, default_flow_style=False, allow_unicode=True)

        with open(filename, 'r') as stream:
            data_loaded = yaml.safe_load(stream)

        print(data_loaded)
        trunc_dictionary_list(data_loaded)

        self.assertDictEqual(data, data_loaded)

    def test_eth_90_ras44(self):

        rasci = dict(self.rem)
        rasci.update({'ras_act': 4,
                      'ras_elec': 4,
                      'ras_elec_alpha': 2,
                      'ras_elec_beta': 2,
                      'ras_occ': 6})

        # create qchem input
        txt_input = create_qchem_input(self.molecule, **rasci)

        # calculate and parse qchem output
        output = get_output_from_qchem(txt_input, processors=4)
        print(output)
        data = basic_rasci(output)
        data = modify_dictionary(data)

        print(data)

        filename = self.__class__.__name__ + '_ras44.yaml'

        # creat reference file
        if ctrl_print:
            with open(filename, 'w') as outfile:
                yaml.dump(data, outfile, default_flow_style=False, allow_unicode=True)

        with open(filename, 'r') as stream:
            data_loaded = yaml.safe_load(stream)

        print(data_loaded)
        trunc_dictionary_list(data_loaded)

        self.assertDictEqual(data, data_loaded)

    def test_eth_90_ras66(self):

        rasci = dict(self.rem)
        rasci.update({'ras_act': 6,
                      'ras_elec': 6,
                      'ras_elec_alpha': 3,
                      'ras_elec_beta': 3,
                      'ras_occ': 5})

        # create qchem input
        txt_input = create_qchem_input(self.molecule, **rasci)

        # calculate and parse qchem output
        output = get_output_from_qchem(txt_input, processors=4)
        print(output)
        data = basic_rasci(output)
        data = modify_dictionary(data)

        print(data)

        filename = self.__class__.__name__ + '_ras66.yaml'

        # create reference file
        if ctrl_print:
            with open(filename, 'w') as outfile:
                yaml.dump(data, outfile, default_flow_style=False, allow_unicode=True)

        with open(filename, 'r') as stream:
            data_loaded = yaml.safe_load(stream)

        print(data_loaded)
        trunc_dictionary_list(data_loaded)

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
            atomic_elements=['C', 'C', 'H', 'H', 'H', 'H'],
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
                    'ras_elec_alpha': 1,
                    'ras_elec_beta': 1,
                    'ras_occ': 7,
                    'ras_spin_mult': 0,
                    'ras_roots': 7,
                    'ras_do_hole': True,
                    'ras_do_part': True,
                    'ras_sts_tm': True,
                    'ras_natorb': True,
                    # rasci sr-dft
                    'ras_srdft': False,
                    'ras_omega': 400,
                    'ras_srdft_cor': 'srpw92',
                    'ras_srdft_exc': 'srpbe',
                    'ras_srdft_damp': 0.5,
                    # rasci print level
                    'ras_print': 2,
                    # Frozen
                    'n_frozen_core': 0,
                    'n_frozen_virt': 0}


    def test_eth_dist_ras22(self):

        rasci = dict(self.rem)

        # create qchem input
        txt_input = create_qchem_input(self.molecule, **rasci)

        # calculate and parse qchem output
        output = get_output_from_qchem(txt_input, processors=4)
        print(output)
        data = basic_rasci(output)
        data = modify_dictionary(data)

        print(data)

        filename = self.__class__.__name__ + '_ras22.yaml'

        # creat reference file
        if ctrl_print:
            with open(filename, 'w') as outfile:
                yaml.dump(data, outfile, default_flow_style=False, allow_unicode=True)

        with open(filename, 'r') as stream:
            data_loaded = yaml.safe_load(stream)

        print(data_loaded)
        trunc_dictionary_list(data_loaded)

        self.assertDictEqual(data, data_loaded)


    def test_eth_dist_ras44(self):

        rasci = dict(self.rem)
        rasci.update({'ras_act': 4,
                      'ras_elec': 4,
                      'ras_elec_alpha': 2,
                      'ras_elec_beta': 2,
                      'ras_occ': 6})

        # create qchem input
        txt_input = create_qchem_input(self.molecule, **rasci)

        # calculate and parse qchem output
        output = get_output_from_qchem(txt_input, processors=4)
        print(output)
        data = basic_rasci(output)
        data = modify_dictionary(data)

        print(data)

        filename = self.__class__.__name__ + '_ras44.yaml'

        # creat reference file
        if ctrl_print:
            with open(filename, 'w') as outfile:
                yaml.dump(data, outfile, default_flow_style=False, allow_unicode=True)

        with open(filename, 'r') as stream:
            data_loaded = yaml.safe_load(stream)

        print(data_loaded)
        trunc_dictionary_list(data_loaded)

        self.assertDictEqual(data, data_loaded)

    def test_eth_dist_ras66(self):

        rasci = dict(self.rem)
        rasci.update({'ras_act': 6,
                      'ras_elec': 6,
                      'ras_elec_alpha': 3,
                      'ras_elec_beta': 3,
                      'ras_occ': 5})

        # create qchem input
        txt_input = create_qchem_input(self.molecule, **rasci)

        # calculate and parse qchem output
        output = get_output_from_qchem(txt_input, processors=4)
        print(output)
        data = basic_rasci(output)
        data = modify_dictionary(data)

        print(data)

        filename = self.__class__.__name__ + '_ras66.yaml'

        # creat reference file
        if ctrl_print:
            with open(filename, 'w') as outfile:
                yaml.dump(data, outfile, default_flow_style=False, allow_unicode=True)

        with open(filename, 'r') as stream:
            data_loaded = yaml.safe_load(stream)

        print(data_loaded)
        trunc_dictionary_list(data_loaded)

        self.assertDictEqual(data, data_loaded)
