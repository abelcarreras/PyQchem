import numpy as np
from copy import deepcopy
from pyqchem.basis import basis_to_txt
import hashlib, json
import warnings
from pyqchem.errors import QchemInputWarning, QchemInputError


def normalize_values(value):
    """
    Set all string values (including keys and values of dictionaries) to lower case

    :param value: the values
    :return: normalized values
    """
    def normalize(value):
        if isinstance(value, str):
            return value.lower()
        return value
    if isinstance(value, dict):
        return {normalize(k): normalize(v) for k, v in value.items()}
    else:
        return normalize(value)


class QchemInput:
    """
    Handles the Q-Chem input info
    """

    def __init__(self, molecule,
                 jobtype='sp',
                 method=None,
                 exchange=None,
                 correlation=None,
                 unrestricted=None,
                 basis='6-31G',
                 thresh=14,
                 scf_convergence=8,
                 max_scf_cycles=50,
                 scf_algorithm='diis',
                 purecart=None,
                 # RASCI
                 ras_roots=1,
                 ras_do_hole=True,
                 ras_do_part=True,
                 ras_act=None,
                 ras_act_orb=None,
                 ras_elec=None,
                 ras_elec_alpha=None,
                 ras_elec_beta=None,
                 ras_occ=None,
                 ras_spin_mult=1,
                 ras_sts_tm=False,
                 ras_natorb=False,
                 ras_natorb_state=None,
                 ras_print=1,
                 ras_diabatization_scheme=None,
                 ras_diabatization_states=None,
                 # RASCI SrDFT
                 ras_omega=400,
                 ras_srdft=None,
                 ras_srdft_damp=0.5,
                 ras_srdft_exc=None,
                 ras_srdft_cor=None,
                 ras_srdft_spinpol=0,
                 # SOC
                 calc_soc=False,
                 state_analysis=False,
                 # EOM
                 ee_singlets=False,
                 ee_triplets=False,
                 cc_trans_prop=False,
                 cc_symmetry=True,
                 cc_e_conv=None,
                 cc_t_conv=None,
                 eom_davidson_conv=5,
                 # CIS
                 cis_convergence=6,
                 cis_n_roots=None,
                 cis_singlets=False,
                 cis_triplets=False,
                 cis_ampl_anal=False,
                 loc_cis_ov_separate=False,
                 er_cis_numstate=0,
                 cis_diabath_decompose=False,
                 max_cis_cycles=30,
                 localized_diabatization=None,
                 sts_multi_nroots=None,
                 cc_state_to_opt=None,
                 cis_state_deriv=None,
                 RPA=False,
                 set_iter=30,
                 gui=0,
                 # optimization
                 geom_opt_dmax=300,
                 geom_opt_update=-1,
                 geom_opt_linear_angle=165,
                 geom_opt_coords=-1,
                 geom_opt_tol_gradient=300,
                 geom_opt_tol_displacement=1200,
                 geom_opt_tol_energy=100,
                 geom_opt_max_cycles=50,
                 geom_opt_constrains=None,
                 # solvent
                 solvent_method=None,
                 solvent_params=None,
                 pcm_params=None,
                 # IRC
                 rpath_coords=0,
                 rpath_direction=1,
                 rpath_max_cycles=20,
                 rpath_max_stepsize=150,
                 rpath_tol_displacement=5000,
                 # symmetry
                 symmetry=True,
                 sym_ignore=False,
                 # other
                 n_frozen_core=None,
                 n_frozen_virt=0,
                 mom_start=False,
                 reorder_orbitals=None,
                 namd_nsurfaces=None,
                 scf_print=None,
                 scf_guess=None,
                 scf_energies=None,
                 scf_guess_mix=False,
                 hessian=None,
                 sym_tol=5,
                 mem_total=2000,
                 mem_static=64,
                 skip_scfman=False,
                 # special
                 extra_rem_keywords=None,
                 extra_sections=None,
                 ):

        # put to arguments self._* (will be written explicitly)
        self.purecart = purecart
        for name, value in vars().items():
            if name != 'self':
                # set keywords in lower case
                value = normalize_values(value)
                setattr(self, '_' + name, value)

        # set ras_occ
        if correlation is not None:
            if correlation.lower() == 'rasci':
                if exchange is None:
                    self._exchange = 'hf'
                if ras_occ is None:
                    if ras_elec is not None:
                        self._ras_occ = (np.sum(
                            molecule.get_atomic_numbers()) - ras_elec - molecule.charge) // 2
                    elif ras_elec_alpha is not None or ras_elec_beta is not None:
                        self._ras_occ = (np.sum(
                            molecule.get_atomic_numbers()) - ras_elec_alpha - ras_elec_beta - molecule.charge) // 2
                    else:
                        self._ras_occ = (np.sum(molecule.get_atomic_numbers()) - molecule.charge) // 2
                    self._ras_occ = int(self._ras_occ)
                    warnings.warn(QchemInputWarning('set ras_occ = {}'.format(self._ras_occ)))

        # Handle custom basis set
        if type(basis) is not str:
            self._basis = 'gen'
            self._custom_basis = basis

        # Handle cc_trans_prop
        if isinstance(self._cc_trans_prop, dict):
            self._trans_prop = self._cc_trans_prop
            self._cc_trans_prop = 1
        else:
            self._trans_prop = None

        # Handle reorder mo
        if self._reorder_orbitals is not None:
            if not 'beta' in self._reorder_orbitals:
                self._reorder_orbitals['beta'] = self._reorder_orbitals['alpha']

        # handle explicit guess (from MO coefficients)
        if scf_guess is not None and type(scf_guess) is not str:
            # print('set custom guess')
            self._scf_guess = 'read'
            self._mo_coefficients = scf_guess
        else:
            self._mo_coefficients = None

        if self._ras_srdft is not None:
            from warnings import warn
            warn('Warning! ras_srdft keyword is deprecated, this will be automatically '
                 'activated when using ras_srdft_exc and ras_srdft_cor')
        if ras_srdft_exc is not None or ras_srdft_cor is not None:
           self._ras_srdft = True
        else:
            self._ras_srdft = False

    def __hash__(self):

        # take all keywords defined in input
        keywords = dict(self.__dict__)

        # remove keywords that not affect the results (these keywords will be ignored in hash)
        for key in ['_mem_total', '_mem_static', '_gui', '_set_iter', '_max_scf_cycles',
                    '_geom_opt_max_cycles', '_max_cis_cycles']:
            keywords.pop(key, None)

        # Change molecule object by molecule coordinates (Structure class too complex for JSON)
        keywords['_molecule'] = hash(keywords['_molecule'])

        digest = hashlib.md5(json.dumps(keywords, sort_keys=True).encode()).hexdigest()
        return int(digest, 16)

    def get_txt(self):
        """
        get qchem input in plain text

        :return string: qchem input in plain text
        """

        input_file = ''

        ##############################
        # Molecule definition
        ##############################

        input_file += '$molecule\n'

        input_file += '{} {}\n'.format(self._molecule.charge, self._molecule.multiplicity)

        atomic_elements = self._molecule.get_symbols()

        coordinates = self._molecule.get_coordinates()

        for index, element in enumerate(atomic_elements):
            input_file += (element + '\t' + '{:20.10f} {:20.10f} {:20.10f}\n'.format(*coordinates[index]))

        input_file += '$end\n'

        #####################
        # Rem variables
        #####################

        input_file += '$rem\n'
        input_file += 'jobtype {}\n'.format(self._jobtype)

        if self._exchange is not None:
            input_file += 'exchange {}\n'.format(self._exchange)

        if self._method:
            input_file += 'method {}\n'.format(self._method)

        input_file += 'basis {}\n'.format(self._basis)
        input_file += 'thresh {}\n'.format(self._thresh)
        input_file += 'scf_convergence {}\n'.format(self._scf_convergence)
        input_file += 'scf_algorithm {}\n'.format(self._scf_algorithm)
        input_file += 'max_scf_cycles {}\n'.format(self._max_scf_cycles)
        input_file += 'gui {}\n'.format(self.gui)
        input_file += 'set_iter {}\n'.format(self._set_iter)
        input_file += 'RPA {}\n'.format(self._RPA)
        input_file += 'mem_total {}\n'.format(self._mem_total)
        input_file += 'mem_static {}\n'.format(self._mem_static)
        input_file += 'n_frozen_virtual {}\n'.format(self._n_frozen_virt)
        input_file += 'mom_start {}\n'.format(self._mom_start)
        input_file += 'skip_scfman {}\n'.format(self._skip_scfman)
        input_file += 'scf_guess_mix {}\n'.format(self._scf_guess_mix)

        if self._n_frozen_core is not None:
            input_file += 'n_frozen_core {}\n'.format(self._n_frozen_core)

        if self._solvent_method is not None:
            input_file += 'solvent_method {}\n'.format(self._solvent_method)

        if self._unrestricted is not None:
            input_file += 'unrestricted {}\n'.format(self._unrestricted)

        if self._purecart is not None:
            input_file += 'purecart {}\n'.format(self._purecart)

        if self._correlation is not None:
            input_file += 'correlation {}\n'.format(self._correlation)

            # RasCI variables
            if self._correlation.upper() == 'RASCI':

                if self._ras_natorb_state is not None:
                    input_file += 'ras_natorb_state {}\n'.format(self._ras_natorb_state)
                    self._ras_natorb = True

                input_file += 'ras_roots {}\n'.format(self._ras_roots)
                input_file += 'ras_do_hole {}\n'.format(self._ras_do_hole)
                input_file += 'ras_do_part {}\n'.format(self._ras_do_part)
                input_file += 'ras_occ {}\n'.format(self._ras_occ)
                input_file += 'ras_spin_mult {}\n'.format(self._ras_spin_mult)
                input_file += 'ras_print {}\n'.format(self._ras_print)
                input_file += 'ras_natorb {}\n'.format(self._ras_natorb)
                input_file += 'ras_sts_tm {}\n'.format(self._ras_sts_tm)
                # input_file += 'max_cis_cycles {}\n'.format(self._max_cis_cycles)
                # input_file += 'RAS_RESTR_TYPE {}\n'.format(True)

                if self._ras_elec is not None:
                    input_file += 'ras_elec {}\n'.format(self._ras_elec)
                else:
                    if self._ras_elec_alpha is None and self._ras_elec_beta is None:
                        raise QchemInputError('{} not defined'.format('ras_elec'))

                if self._ras_elec_alpha is not None:
                    input_file += 'ras_elec_alpha {}\n'.format(self._ras_elec_alpha)
                if self._ras_elec_beta is not None:
                    input_file += 'ras_elec_beta {}\n'.format(self._ras_elec_beta)

                if self._ras_act is not None:
                    input_file += 'ras_act {}\n'.format(self._ras_act)
                else:
                    raise QchemInputError('{} not defined'.format('ras_act'))

                if self._ras_act_orb is not None:
                    input_file += 'ras_act_orb [' + ','.join([str(num) for num in self._ras_act_orb]) + ']\n'

                # Sr-DFT
                if self._ras_srdft:
                    input_file += 'ras_srdft {}\n'.format('True')
                    input_file += 'ras_srdft_damp {}\n'.format(self._ras_srdft_damp)
                    input_file += 'ras_srdft_spinpol {}\n'.format(self._ras_srdft_spinpol)
                    input_file += 'ras_omega {}\n'.format(self._ras_omega)

                    if self._ras_srdft_exc is not None:
                        input_file += 'ras_srdft_exc {}\n'.format(self._ras_srdft_exc)

                    if self._ras_srdft_cor is not None:
                        input_file += 'ras_srdft_cor {}\n'.format(self._ras_srdft_cor)
                    else:
                        raise QchemInputError('{} not defined'.format('ras_srdft_cor'))

                # Diabatization
                diab_methods = {'ER': 1, 'Boys': 2, 'DQ': 3, 'Gamma': 4}
                if self._ras_diabatization_states is not None:
                    input_file += 'sts_multi_nroots {}\n'.format(len(self._ras_diabatization_states))
                    input_file += 'cis_diabath_decompose {}\n'.format(len(self._ras_diabatization_scheme))

                    input_file += 'ras_diab_seq_data ['
                    for seq in self._ras_diabatization_scheme:
                        input_file += '{} '.format([num for num in seq['states']] +
                                                   [diab_methods[seq['method']]] +
                                                   [seq['parameters'] if 'parameters' in seq else 0.0]).replace(' ', '')[1:-1]
                        input_file += ','
                    input_file = input_file[:-1] + ']\n'
                    input_file += 'ras_diab_seq_list ' + '{}\n'.format([len(seq['states']) for seq in self._ras_diabatization_scheme]).replace(' ', '')

                # Borrowed keywords
                input_file += 'cis_convergence {}\n'.format(self._cis_convergence)

        #if self._method.upper() in ['EOM-CCSD'] or self._correlation.upper() in ['CCSD']:
        if self._method is not None:

            # EOM
            if self._method.upper() in ['EOM-CCSD']:
                input_file += 'cc_trans_prop {}\n'.format(self._cc_trans_prop)
                input_file += 'cc_symmetry {}\n'.format(self._cc_symmetry)

                if self._cc_e_conv is not None:
                    input_file += 'cc_e_conv {}\n'.format(self._cc_e_conv)
                if self._cc_t_conv is not None:
                    input_file += 'cc_t_conv {}\n'.format(self._cc_t_conv)

                if self._ee_singlets is not False:
                    input_file += 'ee_singlets [' + ','.join([str(num) for num in self._ee_singlets]) + ']\n'

                if self._ee_triplets is not False:
                    input_file += 'ee_triplets [' + ','.join([str(num) for num in self._ee_triplets]) + ']\n'

                input_file += 'eom_davidson_conv {}\n'.format(self._eom_davidson_conv)

        # SOC
        if self._calc_soc is not False:
            input_file += 'calc_soc {}\n'.format(self._calc_soc)

        if self._state_analysis is not False:
            input_file += 'state_analysis {}\n'.format(self._state_analysis)

        # CIS variables
        if self._cis_n_roots is not None:
            input_file += 'cis_convergence {}\n'.format(self._cis_convergence)
            input_file += 'cis_n_roots {}\n'.format(self._cis_n_roots)
            input_file += 'cis_singlets {}\n'.format(self._cis_singlets)
            input_file += 'cis_triplets {}\n'.format(self._cis_triplets)
            input_file += 'cis_ampl_anal {}\n'.format(self._cis_ampl_anal)
            input_file += 'loc_cis_ov_separate {}\n'.format(self._loc_cis_ov_separate)
            input_file += 'er_cis_numstate {}\n'.format(self._er_cis_numstate)
            input_file += 'max_cis_cycles {}\n'.format(self._max_cis_cycles)
        # other
        if self._namd_nsurfaces is not None:
            input_file += 'namd_nsurfaces {}\n'.format(self._namd_nsurfaces)
        if self._sts_multi_nroots is not None:
            input_file += 'sts_multi_nroots {}\n'.format(self._sts_multi_nroots)
        if self._localized_diabatization is not None:
            input_file += 'cis_diabath_decompose {}\n'.format(self._cis_diabath_decompose)
        if self._cc_state_to_opt is not None:
            input_file += 'cc_state_to_opt [{},{}]\n'.format(self._cc_state_to_opt[0], self._cc_state_to_opt[1])
        if self._cis_state_deriv is not None:
            input_file += 'cis_state_deriv {}\n'.format(self._cis_state_deriv)

        if self._scf_print is not None:
            input_file += 'scf_print {}\n'.format(self._scf_print)

        if self._scf_guess is not None:
            input_file += 'scf_guess {}\n'.format(self._scf_guess)

        if self._symmetry is not None:
            input_file += 'symmetry {}\n'.format(self._symmetry)

        if self._sym_ignore is not None:
            input_file += 'sym_ignore {}\n'.format(self._sym_ignore)

        if self._sym_tol is not None:
            input_file += 'sym_tol {}\n'.format(self._sym_tol)

        # optimization
        if self._jobtype.lower() in ['opt', 'ts']:
            input_file += 'geom_opt_dmax {}\n'.format(self._geom_opt_dmax)
            input_file += 'geom_opt_update {}\n'.format(self._geom_opt_update)
            input_file += 'geom_opt_linear_angle {}\n'.format(self._geom_opt_linear_angle)
            input_file += 'geom_opt_coords {}\n'.format(self._geom_opt_coords)
            input_file += 'geom_opt_tol_gradient {}\n'.format(self._geom_opt_tol_gradient)
            input_file += 'geom_opt_tol_displacement {}\n'.format(self._geom_opt_tol_displacement)
            input_file += 'geom_opt_tol_energy {}\n'.format(self._geom_opt_tol_energy)
            input_file += 'geom_opt_max_cycles {}\n'.format(self._geom_opt_max_cycles)
        # IRC
        if self._jobtype.lower() in ['rpath']:
            input_file += 'rpath_coords {}\n'.format(self._rpath_coords)
            input_file += 'rpath_direction {}\n'.format(self._rpath_direction)
            input_file += 'rpath_max_cycles {}\n'.format(self._rpath_max_cycles)
            input_file += 'rpath_max_stepsize {}\n'.format(self._rpath_max_stepsize)
            input_file += 'rpath_tol_displacement {}\n'.format(self._rpath_tol_displacement)

        # Extra keywords
        if self._extra_rem_keywords is not None:
            for key, value in self._extra_rem_keywords.items():
                input_file += '{} {}\n'.format(key, value)

        input_file += '$end\n'

        ##########################
        # Additional sections
        ##########################

        # localized diabatization
        if self._localized_diabatization is not None:
            input_file += '$localized_diabatization\nadiabatic states\n'
            input_file += ' '.join(np.array(self._localized_diabatization, dtype=str))
            input_file += '\n$end\n'

        # Constrains section
        def modulate_angles(type, value):
            if type in ['tors', 'outp', 'linc', 'linp', 'bend']:
                value = np.mod(value + 180, 360) - 180
            if type in ['bend']:
                value = np.abs(value)
            return value

        if self._geom_opt_constrains is not None:
            input_file += '$opt\n'
            input_file += 'CONSTRAINT\n'
            for type, constrains in self._geom_opt_constrains.items():
                for constrain in constrains:
                    input_file += '{} {} {:15.6f}\n'.format(type,
                                                            ' '.join([str(num) for num in constrain['atoms']]),
                                                            modulate_angles(type, constrain['value']))
            input_file += 'ENDCONSTRAINT\n'
            input_file += '$end\n'

        # Diabatization section
        if self._ras_diabatization_states is not None:
            input_file += '$localized_diabatization\n'
            input_file += 'adiabatic states\n'
            input_file += ' '.join([str(num) for num in self._ras_diabatization_states])
            input_file += '\n$end\n'

        # custom basis section
        if self._basis == 'gen':
            input_file += '$basis\n'
            input_file += basis_to_txt(self._custom_basis)
            input_file += '$end\n'

        # reorder orbitals
        if self._reorder_orbitals is not None:
            input_file += '$reorder_mo\n'
            input_file += ' '.join([str(s) for s in self._reorder_orbitals['alpha']]) + '\n'
            input_file += ' '.join([str(s) for s in self._reorder_orbitals['beta']]) + '\n'
            input_file += '$end\n'

        # trans_prop section
        if self._trans_prop is not None:
            input_file += '$trans_prop\n'

            if 'state_list' in self._trans_prop:
                input_file += 'state_list\n'

                if 'ee_singlets' in self._trans_prop['state_list']:
                    for pair in self._trans_prop['state_list']['ee_singlets']:
                        input_file += 'ee_singlets {} {}\n'.format(*pair)

                if 'ee_triplets' in self._trans_prop['state_list']:
                    for pair in self._trans_prop['state_list']['ee_triplets']:
                        input_file += 'ee_triplets {} {}\n'.format(*pair)

                if 'ref' in self._trans_prop['state_list']:
                    input_file += 'ref {}\n'.format(self._trans_prop['state_list']['ref'])

                input_file += 'end_list\n'

            if 'state_pair_list' in self._trans_prop:
                for pair in self._trans_prop['state_pair_list']['ee_singlets']:
                    input_file += '{} {}\n'.format(*pair)
                input_file += 'end_pairs\n'

            if 'calc' in self._trans_prop:
                for calc in self._trans_prop['calc']:
                    input_file += 'CALC {}\n'.format(calc)

            if self._calc_soc is not False and self._calc_soc != 0:
                input_file += 'CALC soc\n'
            input_file += '$end\n'

        # solvent
        if self._solvent_params is not None:
            input_file += '$solvent\n'
            for prop, value in self._solvent_params.items():
                input_file += '{} {}\n'.format(prop, value)
            input_file += '$end\n'

        if self._pcm_params is not None:
            input_file += '$pcm\n'
            for prop, value in self._pcm_params.items():
                input_file += '{} {}\n'.format(prop, value)
            input_file += '$end\n'

        # extra sections
        if self._extra_sections is not None:
            if isinstance(self._extra_sections, list):
                for section in self._extra_sections:
                    input_file += section.get_txt()
            else:
                # Only one section
                input_file += self._extra_sections.get_txt()

        return input_file + "\n"

    def store_mo_file(self, path='.'):
        guess_coeff = self._mo_coefficients
        guess_energies = self._scf_energies

        # set guess in place
        mo_coeffa = np.array(guess_coeff['alpha'], dtype=np.float)

        if 'beta' in guess_coeff:
            mo_coeffb = np.array(guess_coeff['beta'], dtype=np.float)
        else:
            mo_coeffb = mo_coeffa

        if guess_energies is not None:
            mo_enea = np.array(guess_energies['alpha'], dtype=np.float)
            if 'beta' in guess_coeff:
                mo_coeffb = np.array(guess_coeff['beta'], dtype=np.float)
                mo_eneb = np.array(guess_energies['beta'], dtype=np.float)
            else:
                mo_eneb = mo_enea

        else:
            # here we set orbital energies to 0 (no problem if skip_scfman=False)
            # since they will be recalculated
            mo_enea = np.zeros(len(mo_coeffa))
            mo_eneb = np.zeros(len(mo_coeffb))

        guess_file = np.vstack([mo_coeffa, mo_coeffb, mo_enea, mo_eneb]).flatten()
        with open(path + '/53.0', 'w') as f:
            guess_file.tofile(f, sep='')

    def store_energy_file(self, path='.'):
        mo_enea = self._mo_coefficients['alpha']
        energy_file = np.zeros(len(mo_enea))
        warnings.warn('warining: FILE_ENERGY will be set to zeros, this may affect post HF methods')
        with open(path + '/99.0', 'w') as f:
            energy_file.tofile(f, sep='')

    # Access to properties (only a reduced set should be accessible/editable)
    @property
    def mo_coefficients(self):
        return self._mo_coefficients

    @property
    def mo_energies(self):
        return self._scf_energies

    @property
    def hessian(self):
        return self._hessian

    @property
    def gui(self):
        return self._gui

    @gui.setter
    def gui(self, value):
        value = int(value)
        if value < 0 or value > 10:
            raise ValueError('GUI value error')
        self._gui = value

    def get_copy(self):
        """
        Get a copy of the input

        :return:
        """
        return deepcopy(self)

    def update_input(self, dictionary):
        """
        Update the input from data in a dictionary
        Note: already existing parameters will be overwritten

        :param dictionary: parameters to add
        """
        # put to arguments self._* (will be written explicitly)
        for name, value in dictionary.items():
            setattr(self, '_' + name, value)


class CustomSection:
    def __init__(self, title, keywords):
        self._tile = title
        self._data = keywords

    def get_txt(self):
        txt_input = '${}\n'.format(self._tile)
        for prop, value in self._data.items():
            txt_input += ' {} {}\n'.format(prop, value)
        txt_input += '$end\n'
        return txt_input
