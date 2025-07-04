from pyqchem.qc_input import QchemInput
from pyqchem.errors import ParserError, OutputError
from pyqchem.utils import get_sdm
from subprocess import Popen, PIPE, STDOUT
from pathlib import Path
import os, shutil, sys
import numpy as np
import hashlib
import pickle
import warnings
import copy


if sys.version_info[0] < 3 or sys.platform in ["win32", "cygwin"] or os.getenv('PYQCHEM_CACHE') == '1':
    # For python 2.x or Windows use pickle based cache system
    from pyqchem.cache import SimpleCache as CacheSystem
    warnings.warn('Using SimpleCache')
else:
    # For python 3.x use SQL Database based cache system
    from pyqchem.cache import SqlCache as CacheSystem


# Backwards Compatibility
def redefine_calculation_data_filename(filename, compress=False):
    cache = CacheSystem(filename=filename, compress=compress)


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


def local_run(input_file_name, work_dir, fchk_file, use_mpi=False, processors=1, stream_output=False):
    """
    Run Q-Chem locally

    :param input_file_name: Q-Chem input file in plain text format
    :param work_dir:  Scratch directory where calculation run
    :param fchk_file: filename of fchk
    :param use_mpi: use mpi instead of openmp
    :param stream_output: print Q-Chem output live during execution
    :return: output, err: Q-Chem standard output and standard error
    """

    if not use_mpi:
        os.environ["QCTHREADS"] = "{}".format(processors)
        os.environ["OMP_NUM_THREADS"] = "{}".format(processors)
        os.environ["MKL_NUM_THREADS"] = "1"

    os.environ["GUIFILE"] = fchk_file
    qc_dir = os.getenv('QC')
    exe_dir = os.getenv('QC_EXE_DIR') if 'QC_EXE_DIR' in os.environ else 'exe'

    binary = Path(qc_dir).joinpath(exe_dir).joinpath('qcprog.exe')
    command = [binary, Path(work_dir).joinpath(input_file_name), Path(work_dir)]

    if stream_output:
        qchem_process = Popen(command, stdout=PIPE, stdin=PIPE, stderr=STDOUT, cwd=work_dir, text=True, bufsize=1)
        output = ''
        err = ''
        for line_out in qchem_process.stdout:
            print(line_out, end='')  # Print output as generated
            output += line_out
        qchem_process.wait()
    else:
        qchem_process = Popen(command, stdout=PIPE, stdin=PIPE, stderr=PIPE, cwd=work_dir)
        (output, err) = qchem_process.communicate()
        qchem_process.wait()
        output = output.decode(errors='ignore')
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
    remote_params = copy.deepcopy(remote_params)

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

    remote_dir = os.path.join(remote_scratch, 'temp_pyqchem_remote')

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
    """
    Generate additional files on scratch (work dir) for special calculations

    :param input_qchem: QChem input object
    :param work_dir: scratch directory
    """

    # handle custom guess
    if input_qchem.mo_coefficients is not None:
        input_qchem.store_mo_file(work_dir)

    if input_qchem.scf_density is not None:
        input_qchem.store_density_file(work_dir)

    # set scf energy if skip_scfman (to not break)
    # TODO: now SCF energy is set to zero. This may not work for all features.
    if input_qchem._skip_scfman:
        if input_qchem.mo_coefficients is None:
            raise Exception('Explicit MO guess has to be provided for scf_skip')
        input_qchem.store_energy_file(work_dir)

    # Write hessian
    if input_qchem.hessian is not None:
        input_qchem.store_hessian_file(work_dir)

    # write RAS-GUESS
    if input_qchem.ras_guess is not None:
        input_qchem.store_ras_guess_file(work_dir)


def retrieve_additional_files(input_qchem, data_fchk, work_dir, scratch_read_level=0):
    """
    retrieve data from files in scratch data (on development, currently for test only)

    :param input_qchem: QChem input object
    :param data_fchk: FCHK parsed dictionary
    :param work_dir: scratch directory
    :param scratch_read_level: defines what data to retrieve

    :return: dictionary with additional data
    """

    additional_data = {}

    natom = len(input_qchem.molecule.get_coordinates())
    file_list = os.listdir(work_dir)

    # OLD_DIMENSIONS
    if '819.0' in file_list:
        with open(work_dir + '819.0', 'r') as f:
            data = np.fromfile(f, dtype=np.int32)
            norb_alpha, norb_beta = data[0:2]
            norb = norb_alpha
            nbas = norb  # assumption
    else:
        norb = np.shape(data_fchk['coefficients']['alpha'])[0]
        nbas = np.shape(data_fchk['coefficients']['alpha'])[1]


    # MO_COEFS (Already in fchk) in internal order
    if '53.0' in file_list and 'coefficients' in data_fchk:
        with open(work_dir + '53.0', 'r') as f:
            data = np.fromfile(f, dtype=float)
            mo_alpha = data[:norb*nbas].reshape(-1, norb).tolist()
            mo_beta = data[norb*nbas: 2*norb_beta*nbas].reshape(-1, norb_beta).tolist()
            # additional_data['coefficients_internal'] = {'alpha': mo_alpha, 'beta': mo_beta}

            # obtain the order indices between fchk order and Q-Chem internal order of basis functions
            diff_square = get_sdm(data_fchk['coefficients']['alpha'], mo_alpha)

            # get non-repeating indices
            indices = []
            for row in diff_square.T:
                for i in np.argsort(row):
                    if i not in indices:
                        indices.append(int(i))
                        break

            # indices = np.argmin(diff_square, axis=0).tolist()

            # store q-chem index order for later use (e.g  guess)
            data_fchk['coefficients']['qchem_order'] = indices
    else:
        indices = list(range(nbas))

    # FOCK_MATRIX
    if '58.0' in file_list:
        with open(work_dir + '58.0', 'r') as f:
            data = np.fromfile(f, dtype=float)
            fock_alpha = data[:nbas*nbas].reshape(-1, nbas)
            fock_beta = data[nbas*nbas: 2*nbas*nbas].reshape(-1, nbas)

            # set basis functions in fchk order
            fock_alpha = fock_alpha[:, indices]
            fock_alpha = fock_alpha[indices, :]
            fock_beta = fock_beta[:, indices]
            fock_beta = fock_beta[indices, :]

            additional_data['fock_matrix'] = {'alpha': fock_alpha.tolist(), 'beta': fock_beta.tolist()}

    if scratch_read_level == -1:
        # FILE_ENERGY (Not really worth to read it)
        if '99.0' in file_list:
            with open(work_dir + '99.0', 'r') as f:
                data = np.fromfile(f, dtype=float)

        # FILE_DENSITY_MATRIX (Already in fchk)
        if '54.0' in file_list:
            with open(work_dir + '54.0', 'r') as f:
                data = np.fromfile(f, dtype=float)
                density_alpha = data[:nbas*nbas].reshape(-1, nbas)
                density_beta = data[nbas*nbas: 2*nbas*nbas].reshape(-1, nbas)
                # set basis functions in fchk order
                density_alpha = density_alpha[:, indices]
                density_alpha = density_alpha[indices, :]
                density_beta = density_beta[:, indices]
                density_beta = density_beta[indices, :]
                additional_data['scf_density_internal'] = {'alpha': density_alpha.tolist(), 'beta': density_beta.tolist()}

    # HESSIAN_MATRIX
    if '132.0' in file_list:
        with open(work_dir + '132.0', 'r') as f:
            data = np.fromfile(f, dtype=float)
            hessian = data.reshape(-1, natom*3)
            additional_data['hessian'] = hessian.tolist()

    # AO_INTS_DEBUG
    if '21.0' in file_list:
        with open(work_dir + '21.0', 'r') as f:
            data = np.fromfile(f, dtype=float)
            ao_integrals = data.reshape(-1, nbas, nbas, nbas)

            # set basis functions in fchk order
            ao_integrals = ao_integrals[:, :, :, indices]
            ao_integrals = ao_integrals[:, :, indices, :]
            ao_integrals = ao_integrals[:, indices, :, :]
            ao_integrals = ao_integrals[indices, :, :, :]

            additional_data['ao_integrals'] = ao_integrals.tolist()

    if scratch_read_level > 0:
        # FILE_RAS_AMP
        if '704.0' in file_list:
            with open(work_dir + '705.0', 'r') as f:
                ras_energies = np.fromfile(f, dtype=float)
                n_ras_roots = len(ras_energies)

            with open(work_dir + '704.0', 'r') as f:
                data = np.fromfile(f, dtype=float)
                ras_amplitudes = data.reshape(n_ras_roots, -1)
                additional_data['ras_amplitudes'] = ras_amplitudes.tolist()

    return additional_data


def get_output_from_qchem(input_qchem,
                          processors=1,
                          use_mpi=False,
                          scratch=None,
                          read_fchk=False,  # to be deprecated
                          return_electronic_structure=False,
                          parser=None,
                          parser_parameters=None,
                          force_recalculation=False,
                          fchk_only=False,
                          store_full_output=False,
                          delete_scratch=True,
                          stream_output=False,
                          remote=None,
                          scratch_read_level=0):

    """
    Runs qchem and returns the output in the following format:

    1) If return_electronic_structure is requested:
        [output, parsed_fchk]
    2) If return_electronic_structure is not requested:
        [output]

    Note: if parser is set then output contains a dictionary with the parsed info
          else output contains the q-chem output in plain text

    :param input_qchem: QcInput object containing the Q-Chem input
    :param processors: number of threads/processors to use in the calculation
    :param use_mpi: If False use OpenMP (threads) else use MPI (processors)
    :param scratch: Full Q-Chem scratch directory path. If None read from $QCSCRATCH
    :param return_electronic_structure: if True, returns the parsed FCHK file containing the electronic structure
    :param read_fchk: same as return_electronic_structure (to be deprecated)
    :param parser: function to use to parse the Q-Chem output
    :param parser_parameters: additional parameters that parser function may have
    :param force_recalculation: Force to recalculate even identical calculation has already performed
    :param fchk_only: If true, returns electronic structure data from cache ignoring output (to be deprecated)
    :param remote: dictionary containing the data for remote calculation (beta)
    :param store_full_output: store full output in plain text in pkl file
    :param delete_scratch: delete all scratch files when calculation is finished
    :param stream_output: print Q-Chem output live during execution

    :return: output [, electronic_structure]
    """
    from pyqchem.parsers.parser_fchk import parser_fchk
    cache = CacheSystem()

    # back-compatibility layer
    if read_fchk:
        warnings.warn("'read_fchk' will be deprecated, use return_electronic_structure instead",
                      DeprecationWarning, stacklevel=2)
        return_electronic_structure = read_fchk

    if fchk_only:
        warnings.warn("'fchk_only' option will be deprecated",
                      DeprecationWarning, stacklevel=2)

    # Always generate fchk
    if input_qchem.gui is None:
        input_qchem.gui = 2

    if scratch is None:
        scratch = os.getenv('QCSCRATCH')
        if scratch is None:
            warnings.warn('QCSCRATCH environment variable not defined, using workdir')
            scratch = '.'

    work_dir = '{}/pyqchem_{}/'.format(scratch, os.getpid())

    try:
        os.mkdir(work_dir)
    except OSError:
        pass

    # check if parameters is None
    if parser_parameters is None:
        parser_parameters = {}

    # check if full output is stored
    output = cache.retrieve_calculation_data(input_qchem, 'fullout')

    #output, err = cache.calculation_data[(hash(input_qchem), 'fullout')] if (hash(input_qchem), 'fullout') in cache.calculation_data else [None, None]
    elect_struct_data = cache.retrieve_calculation_data(input_qchem, 'fchk')

    # check if repeated calculation
    if not force_recalculation and not store_full_output:  # store_full_output always force re-parsing
        if parser is not None:
            parsed_data = cache.retrieve_calculation_data(hash(input_qchem), parser.__name__)
            if parsed_data is not None:
                if return_electronic_structure:
                    return parsed_data, elect_struct_data
                else:
                    return parsed_data
        else:
            if fchk_only and elect_struct_data is not None:
                return output, elect_struct_data

    # temp filenames generated in temp directory
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
            # print('local:')
            output, err = local_run(temp_filename, work_dir, fchk_filename,
                                    use_mpi=use_mpi, processors=processors, stream_output=stream_output)
        else:
            # print('Remote:')
            output, err = remote_run(temp_filename, work_dir, fchk_filename, remote, use_mpi=use_mpi, processors=processors)

        if not finish_ok(output):
            raise OutputError(output, err)

        # parse fchk file & and additional scratch dir files
        if not os.path.isfile(os.path.join(work_dir, fchk_filename)):
            warnings.warn('fchk not found! something may be wrong in calculation')
        else:
            with open(os.path.join(work_dir, fchk_filename)) as f:
                fchk_txt = f.read()

            elect_struct_data = parser_fchk(fchk_txt)
            elect_struct_data.update(retrieve_additional_files(input_qchem,
                                                               elect_struct_data,
                                                               work_dir,
                                                               scratch_read_level))
            cache.store_calculation_data(input_qchem, 'fchk', elect_struct_data)

        if store_full_output:
            cache.store_calculation_data(input_qchem, 'fullout', output)

    # plot directory
    plot_dir = work_dir + 'plots'
    if os.path.isdir(plot_dir):
        plots_dict = {}
        for file_name in os.listdir(plot_dir):
            plots_dict[file_name] = open(plot_dir + '/' + file_name, 'r').read()

        elect_struct_data.update({'plots': plots_dict})

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
            raise ParserError(parser.__name__, 'Undefined error', output)

        cache.store_calculation_data(input_qchem, parser.__name__, output)

    if delete_scratch:
        shutil.rmtree(work_dir)

    if return_electronic_structure:
        return output, elect_struct_data
    else:
        return output


def get_input_hash(data):
    return hashlib.md5(data.encode()).hexdigest()

