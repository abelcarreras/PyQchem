from pyqchem.qchem_core import get_output_from_qchem, create_qchem_input
from pyqchem.parsers.parser_rasci_basic import basic_rasci
from pyqchem.parsers.parser_optimization import basic_optimization
from pyqchem.structure import Structure
import yaml
import unittest


def trunc_dictionary_list(w, decimal=5):

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

class HydrogenTest(unittest.TestCase):

    def setUp(self):
        self.assertDictEqual.__self__.maxDiff = None

        # generate molecule
        self.molecule = Structure(coordinates=[[0.0, 0.0, 0.0],
                                               [0.0, 0.0, 1.5]],
                             atomic_elements=['H', 'H'],
                             charge=0,
                             multiplicity=1)

#        self.rasci = {'jobtype': 'sp',
#                      'exchange' : 'hf'}


    def test_srdft(self):

 #       rasci2 = dict(self.rasci)
 #       rasci2.update({'ras_act': 2,
 #                      'ras_hole': 4})

 #       txt_input = create_qchem_input(self.molecule, **rasci2)

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
                                       ras_srdft=True,
                                       ras_omega=400,
                                       ras_srdft_cor='srpw92',
                                       ras_srdft_exc='srpbe',
                                       ras_natorb=False,
                                       ras_print=0,
                                       set_iter=30,
                                       ras_srdft_damp=0.5)

        # calculate and parse qchem output
        output = get_output_from_qchem(txt_input, processors=4)
        print(output)
        data = basic_rasci(output)

        filename = self.__class__.__name__ + '_srdft.yaml'

        #with open(filename, 'w') as outfile:
        #      yaml.dump(data, outfile, default_flow_style=False, allow_unicode=True)

        data = modify_dictionary(data)

        with open(filename, 'r') as stream:
            data_loaded = yaml.safe_load(stream)

        print(data_loaded)
        data_loaded = modify_dictionary(data_loaded)

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
                                       ras_do_hole=True,
                                       ras_sts_tm=True)

        # calculate and parse qchem output
        output = get_output_from_qchem(txt_input, processors=4)
        print(output)

        data = basic_rasci(output)
        print(data)

        filename = self.__class__.__name__ + '_rasci.yaml'

        # with open(filename, 'w') as outfile:
        #      yaml.dump(data, outfile, default_flow_style=False, allow_unicode=True)

        data = modify_dictionary(data)

        with open(filename, 'r') as stream:
            data_loaded = yaml.safe_load(stream)

        print(data_loaded)
        data_loaded = modify_dictionary(data_loaded)

        self.assertDictEqual(data, data_loaded)


class WaterTest(unittest.TestCase):

    def setUp(self):
        self.assertDictEqual.__self__.maxDiff = None

        # generate molecule
        molecule = Structure(coordinates=[[0.0, 1.0, 0.0],
                                          [0.0, 0.0, 1.0],
                                          [0.0, 0.0, -1.0]],
                                  atomic_elements=['O', 'H', 'H'],
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
                                            parser=basic_optimization)

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
                                       ras_do_hole=True,
                                       ras_sts_tm=True)

        # calculate and parse qchem output
        output = get_output_from_qchem(txt_input, processors=4)
        print(output)
        data = basic_rasci(output)

        print(data)

        filename = self.__class__.__name__ + '_rasci.yaml'

        # with open(filename, 'w') as outfile:
        #     yaml.dump(data, outfile, default_flow_style=False, allow_unicode=True)

        data = modify_dictionary(data)

        with open(filename, 'r') as stream:
            data_loaded = yaml.safe_load(stream)

        print(data_loaded)
        data_loaded = modify_dictionary(data_loaded)

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
                                       ras_srdft=True,
                                       ras_omega=300,
                                       ras_srdft_cor='srpbe',
                                       ras_srdft_exc='srlsda',
                                       ras_natorb=False,
                                       ras_print=0,
                                       set_iter=30,
                                       ras_srdft_damp=0.5)

        # calculate and parse qchem output
        output = get_output_from_qchem(txt_input, processors=4)
        print(output)
        data = basic_rasci(output)

        print(data)

        filename = self.__class__.__name__ + '_srdft.yaml'

        # with open(filename, 'w') as outfile:
        #     yaml.dump(data, outfile, default_flow_style=False, allow_unicode=True)

        data = modify_dictionary(data)

        with open(filename, 'r') as stream:
            data_loaded = yaml.safe_load(stream)

        print(data_loaded)
        data_loaded = modify_dictionary(data_loaded)

        self.assertDictEqual(data, data_loaded)

