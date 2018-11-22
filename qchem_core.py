__author__ = 'Abel Carreras'

import os
from subprocess import Popen, PIPE
import numpy as np
import hashlib
import pickle

try:
    with open('calculation_data.pkl', 'rb') as input:
        calculation_data = pickle.load(input)
        print('Loaded data from calculation_data.pkl')
        # print(calculation_data)
except FileNotFoundError:
    calculation_data = {}


def create_qchem_input(molecule,
                       jobtype='sp',
                       exchange='HF',
                       correlation=None,
                       basis='6-31G',
                       thresh=14,
                       scf_convergence=8,
                       multiplicity=1,
                       ras_roots=1,
                       ras_do_hole=True,
                       ras_do_part=True,
                       ras_act=None,
                       ras_elec=None,
                       # ras_occ=None,
                       ras_spin_mult=1,
                       ras_omega=400,
                       ras_srdft=False,
                       ras_srdft_damp=0.5,
                       ras_srdft_exc=None,
                       ras_srdft_cor=None,
                       ras_srdft_spinpol=0,
                       ras_sts_tm=False,
                       ras_natorb=False,
                       ras_print=4,
                       # cis
                       cis_convergence=6,
                       cis_n_roots=None,
                       cis_singlets=False,
                       cis_triplets=False,
                       cis_ampl_anal=False,
                       loc_cis_ov_separate=False,
                       er_cis_numstate=0,
                       cis_diabath_decompose=False,
                       localized_diabatization=None,
                       RPA=False,
                       set_iter=30,
                       gui=0):

    if ras_elec is not None:
        ras_occ = np.sum(molecule.get_atomic_numbers()) - ras_elec - molecule.charge
    else:
        ras_occ = np.sum(molecule.get_atomic_numbers()) - molecule.charge
    print('ras_occ = {}'.format(ras_occ))

    input_file = ''

    # Molecule definition
    input_file += '$molecule\n'

    input_file += '{} {}\n'.format(molecule.charge, multiplicity)

    atomic_elements = molecule.get_atomic_elements()

    coordinates = molecule.get_coordinates()

    for index, element in enumerate(atomic_elements):
        input_file += (element + '\t' +
                       str(coordinates[index][0]) + '\t' +
                       str(coordinates[index][1]) + '\t' +
                       str(coordinates[index][2]) + '\n')

    input_file += '$end\n'

    # Rem variables
    input_file += '$rem\n'
    input_file += 'jobtype {}\n'.format(jobtype)
    input_file += 'exchange {}\n'.format(exchange)
    input_file += 'basis {}\n'.format(basis)
    input_file += 'thresh {}\n'.format(thresh)
    input_file += 'scf_convergence {}\n'.format(scf_convergence)
    input_file += 'gui {}\n'.format(gui)
    # input_file += 'purecart {}\n'.format(2)
    input_file += 'set_iter {}\n'.format(set_iter)
    input_file += 'RPA {}\n'.format(RPA)

    if correlation is not None:
        input_file += 'correlation {}\n'.format(correlation)

        # RasCI variables
        if correlation.upper() == 'RASCI':
            input_file += 'ras_roots {}\n'.format(ras_roots)
            input_file += 'ras_do_hole {}\n'.format(ras_do_hole)
            input_file += 'ras_do_part {}\n'.format(ras_do_part)
            input_file += 'ras_occ {}\n'.format(ras_occ)
            input_file += 'ras_omega {}\n'.format(ras_omega)
            input_file += 'ras_spin_mult {}\n'.format(ras_spin_mult)
            input_file += 'ras_print {}\n'.format(ras_print)
            input_file += 'ras_natorb {}\n'.format(ras_natorb)
            input_file += 'ras_sts_tm {}\n'.format(ras_sts_tm)
            # input_file += 'RAS_RESTR_TYPE {}\n'.format(True)
            if ras_act is not None:
                input_file += 'ras_act {}\n'.format(ras_act)
            else:
                print('test')
                raise Exception('{} not defined'.format('ras_act'))

            if ras_elec is not None:
                input_file += 'ras_elec {}\n'.format(ras_elec)
            else:
                raise Exception('{} not defined'.format('ras_elec'))

            if ras_act is not None:
                input_file += 'ras_act {}\n'.format(ras_act)
            else:
                raise Exception('{} not defined'.format('ras_act'))

            # Sr-DFT
            if ras_srdft:
                input_file += 'ras_srdft {}\n'.format('True')
                input_file += 'ras_srdft_damp {}\n'.format(ras_srdft_damp)
                input_file += 'ras_srdft_spinpol {}\n'.format(ras_srdft_spinpol)

                if ras_srdft_exc is not None:
                    input_file += 'ras_srdft_exc {}\n'.format(ras_srdft_exc)
                else:
                    raise Exception('{} not defined'.format('ras_srdft_exc'))

                if ras_srdft_cor is not None:
                    input_file += 'ras_srdft_cor {}\n'.format(ras_srdft_cor)
                else:
                    raise Exception('{} not defined'.format('ras_srdft_cor'))
    # CIS variables
    if cis_n_roots is not None:
        input_file += 'cis_convergence {}\n'.format(cis_convergence)
        input_file += 'cis_n_roots {}\n'.format(cis_n_roots)
        input_file += 'cis_singlets {}\n'.format(cis_singlets)
        input_file += 'cis_triplets {}\n'.format(cis_triplets)
        input_file += 'cis_ampl_anal {}\n'.format(cis_ampl_anal)
        input_file += 'loc_cis_ov_separate {}\n'.format(loc_cis_ov_separate)
        input_file += 'er_cis_numstate {}\n'.format(er_cis_numstate)
        input_file += 'cis_diabath_decompose {}\n'.format(cis_diabath_decompose)

    input_file += '$end\n'

    # localized diabatization
    if localized_diabatization is not None:
        input_file += '$localized_diabatization\nadiabatic states\n'
        input_file += ' '.join(np.array(localized_diabatization, dtype=str))
        input_file += '\n$end\n'

    return input_file + "\n"


def parse_output(get_output_function):

    global calculation_data

    def func_wrapper(*args, **kwargs):
        parser = kwargs.pop('parser', None)
        parser_parameters = kwargs.pop('parser_parameters', {})
        force_recalculation = kwargs.pop('force_recalculation', False)

        if parser is not None:
            hash = get_input_hash(args[0]+'{}'.format(parser.__name__))
            if hash in calculation_data and not force_recalculation:
                print('already calculated. Skip')
                return calculation_data[hash]

        output, err = get_output_function(*args, **kwargs)

        if len(err) > 0:
            print(output[-800:])
            raise Exception('q-chem calculation finished with error')

        if parser is None:
            return output

        parsed_output = parser(output, **parser_parameters)

        calculation_data[hash] = parsed_output
        with open('calculation_data.pkl', 'wb') as output:
            pickle.dump(calculation_data, output, pickle.HIGHEST_PROTOCOL)

        return parsed_output

    return func_wrapper


@parse_output
def get_output_from_qchem(input_data, processors=1, binary='qchem', use_mpi=False, scratch=None):

    if scratch is None:
        scratch = os.environ['QCSCRATCH']
    # print('scratch', scratch)

    temp_file_name = scratch + '/qchem_temp_{}'.format(os.getpid())
    # temp_file_name = tempfile.gettempdir() + '/qchem_temp_{}'.format(os.getpid())
    qchem_input_file = open(temp_file_name,mode='w')
    qchem_input_file.write(input_data)
    qchem_input_file.close()
    # print(temp_file_name)
    if use_mpi:
        flag = '-np'
    else:
        flag = '-nt'
    command = binary + ' {} {} '.format(flag, processors) + ' {} '.format(temp_file_name)
    # print(command)

    qchem_process = Popen(command, stdout=PIPE, stdin=PIPE, stderr=PIPE, shell=True)
    (output, err) = qchem_process.communicate(input=input_data.encode())
    qchem_process.wait()
    os.remove(temp_file_name)

    return output.decode(), err.decode()


def get_input_hash(data):
    return hashlib.md5(data.encode()).hexdigest()

