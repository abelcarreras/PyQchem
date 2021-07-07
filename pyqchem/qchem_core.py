from pyqchem.qc_input import QchemInput
from pyqchem.errors import ParserError, OutputError
from subprocess import Popen, PIPE
import os, shutil
import numpy as np
import hashlib
import pickle
import warnings
import sys

if sys.version_info[0] < 3 or sys.platform.startswith('win'):
    # For python 2.x or Windows use pickle based cache system
    from pyqchem.cache import SimpleCache as CacheSystem
else:
    # For python 3.x use SQL Database based cache system
    from pyqchem.cache import SqlCache as CacheSystem


# Backwards Compatibility
def redefine_calculation_data_filename(filename):
    cache = CacheSystem(filename=filename)


# Check if calculation finished ok
def finish_ok(output):
    return output[-1000:].find('Thank you very much for using Q-Chem') != -1


# Return qchem version
def get_version_output(output):

    class QChemVersion:
        def __init__(self, string):

            string_version = string.split()[1]
            self._major = string_version.split('.')[0]
            self._minor = string_version.split('.')[1]

            string_branch = string.split()[2]
            self._devel = True if '(devel)' in string_branch else False

        def __str__(self):
            dev = 'dev' if self.is_development else ''
            return '{}.{} {}'.format(self.major, self.minor, dev)

        def __eq__(self, other):

            """
            Here put the logic for more sophisticated comparison
            between versions

            :param other: string/QchemVersion
            :return:
            """

            if isinstance(other, QChemVersion):
                return True if self.__str__() == other.__str__() else False

            o_major = other.split('.')[0]
            o_minor = other.split('.')[1]

            if int(o_major) == self.major:

                # handle expresions like 2.3+
                if '+' in o_minor[-1]:
                    if self.minor >= int(o_minor[:-1]):
                        return True
                    else:
                        return False

                if int(o_minor) == self.minor:
                    return True

            return False

        @property
        def major(self):
            return int(self._major)

        @property
        def minor(self):
            return int(self._minor)

        @property
        def is_development(self):
            return self._devel

    index = output[:500].find('\n Q-Chem')
    string = output[index: index + 30]

    return QChemVersion(string)


def get_compatibility_list_from_parser(parser):
    docstring = parser.__doc__
    if docstring is None:
        return None

    lines = docstring.split('\n')
    for line in lines:
        if 'compatibility' in line.lower():
            try:
                return [version.strip() for version in line.split(':')[1].split(',')]
            except IndexError:
                continue

    return None


# Layer of compatibility with old version
def create_qchem_input(*args, **kwargs):
    return QchemInput(*args, **kwargs)


def parse_output(get_output_function):
    """
    to be deprecated

    :param get_output_function:
    :return: parsed output
    """

    cache = CacheSystem()

    def func_wrapper(*args, **kwargs):
        parser = kwargs.pop('parser', None)
        parser_parameters = kwargs.pop('parser_parameters', {})
        store_output = kwargs.pop('store_output', None)

        force_recalculation = kwargs.pop('force_recalculation', False)

        if parser is not None:
            hash_p = (args[0], parser.__name__)
            if hash_p in cache.calculation_data and not force_recalculation:
                print('already calculated. Skip')
                return cache.calculation_data[hash_p]

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

        cache.calculation_data[hash_p] = parsed_output
        with open(cache._calculation_data_filename, 'wb') as output:
            pickle.dump(cache.calculation_data, output, protocol=pickle.DEFAULT_PROTOCOL)

        return parsed_output

    return func_wrapper


def local_run(input_file_name, work_dir, fchk_file, use_mpi=False, processors=1):
    """
    Run Q-Chem locally

    :param input_file_name: Q-Chem input file in plain text format
    :param work_dir:  Scratch directory where calculation run
    :param fchk_file: filename of fchk
    :param use_mpi: use mpi instead of openmp

    :return: output, err: Q-Chem standard output and standard error
    """

    if not use_mpi:
        os.environ["QCTHREADS"] = "{}".format(processors)
        os.environ["OMP_NUM_THREADS"] = "{}".format(processors)
        os.environ["MKL_NUM_THREADS"] = "1"

    os.environ["GUIFILE"] = fchk_file
    qc_dir = os.getenv('QC')
    binary = "{}/exe/qcprog.exe".format(qc_dir)
    # command = binary + ' {} {} '.format(flag, processors) + ' {} '.format(temp_file_name)
    command = binary + ' {} '.format(os.path.join(work_dir, input_file_name)) + ' {} '.format(work_dir)

    qchem_process = Popen(command, stdout=PIPE, stdin=PIPE, stderr=PIPE, shell=True, cwd=work_dir)
    (output, err) = qchem_process.communicate()
    qchem_process.wait()
    output = output.decode()
    err = err.decode()

    return output, err


def remote_run(input_file_name, work_dir, fchk_file, remote_params, use_mpi=False, processors=1):
    """
    Run Q-Chem remotely

    :param input_file: Q-Chem input file in plain text format
    :param work_dir:  Scratch directory where calculation run
    :param fchk_file: filename of fchk
    :param remote_params: connection parameters for paramiko
    :param use_mpi: use mpi instead of openmp

    :return: output, err: Q-Chem standard output and standard error
    """
    import paramiko

    # get precommands
    commands = remote_params.pop('precommand', [])
    remote_scratch = remote_params.pop('remote_scratch', None)

    # Setup SSH connection
    ssh = paramiko.SSHClient()
    ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    ssh.connect(**remote_params)

    ssh.get_transport()
    sftp = ssh.open_sftp()
    print('connected to {}..'.format(remote_params['hostname']))

    # Define temp remote dir
    _, stdout, _ = ssh.exec_command('pwd', get_pty=True)

    if remote_scratch is None:
        remote_scratch = stdout.read().decode().strip('\n').strip('\r')

    remote_dir = '{}/temp_pyqchem_remote/'.format(remote_scratch)

    # Create temp directory in remote machine
    try:
        sftp.mkdir(remote_dir)
    except OSError:
        pass
    sftp.chdir(remote_dir)

    # Copy all files in local workdir to remote machine
    file_list = os.listdir(work_dir)
    for file in file_list:
        sftp.put(os.path.join(work_dir, file), '{}'.format(file))

    flag = '-np' if use_mpi else '-nt'

    # Define commands to run Q-Chem in remote machine
    commands += ['cd  {}'.format(remote_dir),  # go to remote work dir
                 'qchem {} {} {}'.format(flag, processors, input_file_name)]  # run qchem

    # Execute command in remote machine
    stdin, stdout, stderr = ssh.exec_command('bash -l -c "{}"'.format(';'.join(commands)), get_pty=True)

    # Reformat output/error files
    output = ''.join(stdout.readlines())
    error = ''.join(stderr.readlines())

    # get files and remove them from remote server
    for file in sftp.listdir():
        sftp.get(os.path.join(remote_dir, file), os.path.join(work_dir, file))
        sftp.remove(os.path.join(remote_dir, file))
    sftp.rmdir(remote_dir)

    sftp.close()
    ssh.close()

    # Rename fchk file to match expected name
    if input_file_name + '.fchk' in os.listdir(work_dir):
        os.rename(os.path.join(work_dir, input_file_name + '.fchk'), os.path.join(work_dir, fchk_file))

    return output, error



def generate_additional_files(input_qchem, work_dir):
    # Hessian
    if input_qchem.hessian is not None:
        # ndim = len(input_qchem.hessian)
        hessian_triu = np.array(input_qchem.hessian)
        with open(work_dir + '132.0', 'w') as f:
            hessian_triu.tofile(f, sep='')


def get_output_from_qchem(input_qchem,
                          processors=1,
                          use_mpi=False,
                          scratch=None,
                          read_fchk=False,
                          parser=None,
                          parser_parameters=None,
                          force_recalculation=False,
                          fchk_only=False,
                          store_full_output=False,
                          delete_scratch=True,
                          remote=None,
                          strict_policy=False):
    """
    Runs qchem and returns the output in the following format:

    1) If read_fchk is requested:
        [output, parsed_fchk]
    2) If read_fchk is not requested:
        [output]

    Note: if parser is set then output contains a dictionary with the parsed info
          else output contains the q-chem output in plain text

    read_fchk: contains a dictionary with the parsed info inside fchk file.

    :param input_qchem: QcInput object containing the Q-Chem input
    :param processors: number of threads/processors to use in the calculation
    :param use_mpi: If False use OpenMP (threads) else use MPI (processors)
    :param scratch: Full Q-Chem scratch directory path. If None read from $QCSCRATCH
    :param read_fchk: if True, returns the parsed FCHK file containing the electronic structure
    :param parser: function to use to parse the Q-Chem output
    :param parser_parameters: additional parameters that parser function may have
    :param force_recalculation: Force to recalculate even identical calculation has already performed
    :param fchk_only: If true, returns only the electronic structure data parsed from FCHK file
    :param remote: dictionary containing the data for remote calculation (beta)
    :param store_full_output: store full output in plain text in pkl file
    :param delete_scratch: delete all scratch files when calculation is finished

    :return: output [, fchk_dict]
    """
    from pyqchem.parsers.parser_fchk import parser_fchk
    cache = CacheSystem()

    # Always generate fchk
    if input_qchem.gui is None or input_qchem.gui < 1:
        input_qchem.gui = 2

    if scratch is None:
        scratch = os.getenv('QCSCRATCH')

    work_dir = '{}/pyqchem_{}/'.format(scratch, os.getpid())

    try:
        os.mkdir(work_dir)
    except OSError:
        pass

    # handle custom guess
    if input_qchem.mo_coefficients is not None:
        input_qchem.store_mo_file(work_dir)

    # set scf energy if skip_scfman (to not break)
    # TODO: now SCF energy is set to zero. This works for very little features.
    if input_qchem._skip_scfman:
        input_qchem.store_energy_file(work_dir)

    # check if parameters is None
    if parser_parameters is None:
        parser_parameters = {}

    # check if full output is stored
    output = cache.retrieve_calculation_data(input_qchem, 'fullout')

    #output, err = cache.calculation_data[(hash(input_qchem), 'fullout')] if (hash(input_qchem), 'fullout') in cache.calculation_data else [None, None]
    data_fchk = cache.retrieve_calculation_data(input_qchem, 'fchk')

    # check if repeated calculation
    if not force_recalculation and not store_full_output:  # store_full_output always force re-parsing
        if parser is not None:
            parsed_data = cache.retrieve_calculation_data(hash(input_qchem), parser.__name__)
            if parsed_data is not None:
                if read_fchk:
                    return parsed_data, data_fchk
                else:
                    return parsed_data
        else:
            if fchk_only and data_fchk is not None:
                return output, data_fchk

    fchk_filename = 'qchem_temp_{}.fchk'.format(os.getpid())
    temp_filename = 'qchem_temp_{}.inp'.format(os.getpid())

    # generate the input in TXT form
    input_txt = input_qchem.get_txt()
    qchem_input_file = open(os.path.join(work_dir, temp_filename), mode='w')
    qchem_input_file.write(input_txt)
    qchem_input_file.close()

    # generate extra files in calculation directory
    generate_additional_files(input_qchem, work_dir)

    # Q-Chem calculation
    if output is None or force_recalculation is True:
        if remote is None:
            output, err = local_run(temp_filename, work_dir, fchk_filename, use_mpi=use_mpi, processors=processors)
        else:
            output, err = remote_run(temp_filename, work_dir, fchk_filename, remote, use_mpi=use_mpi, processors=processors)

        if not finish_ok(output):
            raise OutputError(output, err)

        # parse fchk file
        if not os.path.isfile(os.path.join(work_dir, fchk_filename)):
            warnings.warn('fchk not found! something may be wrong in calculation')
        else:
            with open(os.path.join(work_dir, fchk_filename)) as f:
                fchk_txt = f.read()

            data_fchk = parser_fchk(fchk_txt)
            cache.store_calculation_data(input_qchem, 'fchk', data_fchk)

        if store_full_output:
            cache.store_calculation_data(input_qchem, 'fullout', output)

    if parser is not None:

        # Check parser compatibility
        version = get_version_output(output)
        compatibility_list = get_compatibility_list_from_parser(parser)
        if compatibility_list is not None:
            if version not in compatibility_list:
                warnings.warn('Parser "{}" may not be compatible with Q-Chem {}'.format(parser.__name__, version))

        # minimum functionality for parser error capture
        try:
            output = parser(output, **parser_parameters)
        except:
            raise ParserError(parser.__name__, 'Undefined error')

        cache.store_calculation_data(input_qchem, parser.__name__, output)

    if delete_scratch:
        shutil.rmtree(work_dir)

    if read_fchk:
        return output, data_fchk
    else:
        return output


def get_input_hash(data):
    return hashlib.md5(data.encode()).hexdigest()

