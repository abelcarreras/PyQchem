__author__ = 'Abel Carreras'

import os
from subprocess import Popen, PIPE
import numpy as np
import hashlib
import pickle
import warnings
from pyqchem.qc_input import QchemInput

__calculation_data_filename__ = 'calculation_data.pkl'
try:
    with open(__calculation_data_filename__, 'rb') as input:
        calculation_data = pickle.load(input)
        print('Loaded data from {}'.format(__calculation_data_filename__))
except FileNotFoundError:
    calculation_data = {}


# Layer of compatibility with old version
def create_qchem_input(*args, **kwargs):
    return QchemInput(*args, **kwargs)


def parse_output(get_output_function):
    """
    to be deprecated

    :param get_output_function:
    :return:
    """

    global calculation_data

    def func_wrapper(*args, **kwargs):
        parser = kwargs.pop('parser', None)
        parser_parameters = kwargs.pop('parser_parameters', {})
        store_output = kwargs.pop('store_output', None)

        force_recalculation = kwargs.pop('force_recalculation', False)

        if parser is not None:
            hash = get_input_hash(args[0]+'{}'.format(parser.__name__))
            if hash in calculation_data and not force_recalculation:
                print('already calculated. Skip')
                return calculation_data[hash]

        output, err = get_output_function(*args, **kwargs)

        if store_output is not None:
            with open('{}'.format(store_output), 'w') as f:
                f.write(output)

        if len(err) > 0:
            print(output[-800:])
            print(err)
            raise Exception('q-chem calculation finished with error')

        if parser is None:
            return output

        parsed_output = parser(output, **parser_parameters)

        calculation_data[hash] = parsed_output
        with open(__calculation_data_filename__, 'wb') as output:
            pickle.dump(calculation_data, output, pickle.HIGHEST_PROTOCOL)

        return parsed_output

    return func_wrapper


def get_output_from_qchem(input_qchem,
                          processors=1,
                          use_mpi=False,
                          scratch=None,
                          read_fchk=False,
                          parser=None,
                          parser_parameters={},
                          force_recalculation=False,
                          fchk_only=False):
    """
    Runs qchem and returns the output in the following format:
    1) If read_fchk is requested:
        [output, error, parsed_fchk]
    2) If read_fchk is not requested:
        [output, error]

    Note: if parser is set then output contains a dictionary with the parsed info
          else output contains the q-chem output in plain text

    error: contains the standard error data from the calculation (if all OK, then should contain nothing)
    read_fchk: contains a dictionary with the parsed info inside fchk file.

    :param input_qchem:
    :param processors:
    :param use_mpi:
    :param scratch:
    :param read_fchk:
    :param parser:
    :param parser_parameters:
    :param force_recalculation:
    :param fchk_only:
    :return output, error[, fchk_dict]:
    """

    # check gui > 2 if read_fchk
    if read_fchk:
        if input_qchem.gui is None or input_qchem.gui < 1:
            input_qchem.gui = 2

    if scratch is None:
        scratch = os.environ['QCSCRATCH']

    scratch_dir = '{}/qchem{}'.format(scratch, os.getpid())

    # check scf_guess if guess
    if input_qchem.mo_coefficients is not None:
        guess = input_qchem.mo_coefficients
        # set guess in place
        mo_coeffa = np.array(guess['alpha'], dtype=np.float)
        l = len(mo_coeffa)
        if 'beta' in guess:
            mo_coeffb = np.array(guess['beta'], dtype=np.float)
        else:
            mo_coeffb = mo_coeffa

        mo_ene = np.zeros(l)

        guess_file = np.vstack([mo_coeffa, mo_ene, mo_coeffb, mo_ene]).flatten()
        with open(scratch_dir + '53.0', 'w') as f:
            guess_file.tofile(f, sep='')

    input_txt = input_qchem.get_txt()

    # Check if already calculate
    hash_fchk = get_input_hash(input_txt + '__fchk__')

    if not force_recalculation:

        data_fchk = None
        if read_fchk and hash_fchk in calculation_data:
            # warnings.warn('already fchk calculated. Skip')
            data_fchk = calculation_data[hash_fchk]

        if parser is not None:
            hash = get_input_hash(input_txt + '{}'.format(parser.__name__))
            if hash in calculation_data:
                # hash = get_input_hash(input_txt + '{}'.format(parser.__name__))
                # warnings.warn('already calculated. Skip')
                data = calculation_data[hash]
                err = ''.encode()
                if data_fchk is None:
                    return data, err
                else:
                    return data, err, data_fchk
        else:
            if fchk_only and data_fchk is not None:
                return None, None, data_fchk

    temp_file_name = scratch_dir + '/qchem_temp_{}.inp'.format(os.getpid())

    try:
        os.mkdir(scratch_dir)
    except FileExistsError:
        pass

    qchem_input_file = open(temp_file_name, mode='w')
    qchem_input_file.write(input_txt)
    qchem_input_file.close()
    # print(temp_file_name)
    if use_mpi:
        flag = '-np'
    else:
        flag = '-nt'
        os.environ["QCTHREADS"] = "{}".format(processors)
        os.environ["OMP_NUM_THREADS"] = "{}".format(processors)
        os.environ["MKL_NUM_THREADS"] = "1"

    fchk_file = 'qchem_temp_{}.fchk'.format(os.getpid())
    os.environ["GUIFILE"] = fchk_file

    qc_dir = os.environ['QC']

    binary = "{}/exe/qcprog.exe".format(qc_dir)
    #command = binary + ' {} {} '.format(flag, processors) + ' {} '.format(temp_file_name)
    command = binary + ' {} '.format(temp_file_name) + ' {} '.format(scratch)

    # print(command)

    qchem_process = Popen(command, stdout=PIPE, stdin=PIPE, stderr=PIPE, shell=True)
    (output, err) = qchem_process.communicate(input=input_txt.encode())
    qchem_process.wait()
    os.remove(temp_file_name)

    output = output.decode()
    err = err.decode()

    if parser is not None:
        hash = get_input_hash(input_txt + '{}'.format(parser.__name__))
        output = parser(output, **parser_parameters)
        calculation_data[hash] = output
        with open(__calculation_data_filename__, 'wb') as f:
            pickle.dump(calculation_data, f, pickle.HIGHEST_PROTOCOL)

    if read_fchk:
        from pyqchem.parsers.parser_fchk import parser_fchk

        if not os.path.isfile(fchk_file):
            warnings.warn('fchk not found! Make sure the input generates it (gui 2)')
            return output, err, []

        with open('qchem_temp_{}.fchk'.format(os.getpid())) as f:
            fchk_txt = f.read()
        os.remove('qchem_temp_{}.fchk'.format(os.getpid()))

        data_fchk = parser_fchk(fchk_txt)
        calculation_data[hash_fchk] = data_fchk
        with open(__calculation_data_filename__, 'wb') as f:
            pickle.dump(calculation_data, f, pickle.HIGHEST_PROTOCOL)

        return output, err, data_fchk

    return output, err


def get_input_hash(data):
    return hashlib.md5(data.encode()).hexdigest()

