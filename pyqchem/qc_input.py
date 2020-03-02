import numpy as np
from copy import deepcopy
from pyqchem.basis import basis_to_txt
import hashlib, json
import warnings
from pyqchem.errors import QchemInputWarning


class QchemInput:
    """
    Handles the Q-Chem input info
    """

    def __init__(self, molecule,
                 jobtype='sp',
                 method='HF',
                 exchange=None,
                 correlation=None,
                 unrestricted=None,
                 basis='6-31G',
                 thresh=14,
                 scf_convergence=8,
                 max_scf_cycles=50,
                 purecart=None,
                 # RASCI
                 ras_roots=1,
                 ras_do_hole=True,
                 ras_do_part=True,
                 ras_act=None,
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
                 # cis
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
                 geom_opt_coords=-1,
                 geom_opt_tol_gradient=300,
                 geom_opt_tol_displacement=1200,
                 geom_opt_tol_energy=100,
                 geom_opt_max_cycles=50,
                 # other
                 n_frozen_core=None,
                 n_frozen_virt=None,
                 namd_nsurfaces=None,
                 scf_print=None,
                 scf_guess=None,
                 mem_total=2000,
                 mem_static=64
                 ):

        # put to arguments self._* (will be written explicitly)
        for name, value in vars().items():
            if name != 'self':
                # set keywords in lower case
                if type(value) is str:
                    value = value.lower()
                setattr(self, '_' + name, value)

        # set ras_occ
        if correlation is not None:
            if ras_occ is None and correlation.upper() == 'RASCI':
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

        # remove keywords that not affect the results
        for key in ['_mem_total', '_mem_static', '_gui', '_set_iter']:
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

        # Molecule definition
        input_file += '$molecule\n'

        input_file += '{} {}\n'.format(self._molecule.charge, self._molecule.multiplicity)

        atomic_elements = self._molecule.get_atomic_elements()

        coordinates = self._molecule.get_coordinates()

        for index, element in enumerate(atomic_elements):
            input_file += (element + '\t' + '{:20.10f} {:20.10f} {:20.10f}\n'.format(*coordinates[index]))

        input_file += '$end\n'

        # Rem variables
        input_file += '$rem\n'
        input_file += 'jobtype {}\n'.format(self._jobtype)

        if self._exchange is not None:
            input_file += 'exchange {}\n'.format(self._exchange)
        else:
            input_file += 'method {}\n'.format(self._method)

        input_file += 'basis {}\n'.format(self._basis)
        input_file += 'thresh {}\n'.format(self._thresh)
        input_file += 'scf_convergence {}\n'.format(self._scf_convergence)
        input_file += 'max_scf_cycles {}\n'.format(self._max_scf_cycles)
        input_file += 'gui {}\n'.format(self.gui)
        input_file += 'set_iter {}\n'.format(self._set_iter)
        input_file += 'RPA {}\n'.format(self._RPA)
        input_file += 'mem_total {}\n'.format(self._mem_total)
        input_file += 'mem_static {}\n'.format(self._mem_static)

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

                if self._ras_act is not None:
                    input_file += 'ras_act {}\n'.format(self._ras_act)
                else:
                    print('test')
                    raise Exception('{} not defined'.format('ras_act'))

                if self._ras_elec is not None:
                    input_file += 'ras_elec {}\n'.format(self._ras_elec)
                else:
                    if self._ras_elec_alpha is None and self._ras_elec_beta is None:
                        raise Exception('{} not defined'.format('ras_elec'))

                if self._ras_elec_alpha is not None:
                    input_file += 'ras_elec_alpha {}\n'.format(self._ras_elec_alpha)
                if self._ras_elec_beta is not None:
                    input_file += 'ras_elec_beta {}\n'.format(self._ras_elec_beta)

                if self._ras_act is not None:
                    input_file += 'ras_act {}\n'.format(self._ras_act)
                else:
                    raise Exception('{} not defined'.format('ras_act'))

                # Sr-DFT
                if self._ras_srdft:
                    input_file += 'ras_srdft {}\n'.format('True')
                    input_file += 'ras_srdft_damp {}\n'.format(self._ras_srdft_damp)
                    input_file += 'ras_srdft_spinpol {}\n'.format(self._ras_srdft_spinpol)
                    input_file += 'ras_omega {}\n'.format(self._ras_omega)

                    if self._ras_srdft_exc is not None:
                        input_file += 'ras_srdft_exc {}\n'.format(self._ras_srdft_exc)
                    else:
                        raise Exception('{} not defined'.format('ras_srdft_exc'))

                    if self._ras_srdft_cor is not None:
                        input_file += 'ras_srdft_cor {}\n'.format(self._ras_srdft_cor)
                    else:
                        raise Exception('{} not defined'.format('ras_srdft_cor'))

                # Diabatization
                diab_methods = {'ER': 1, 'Boys': 2, 'DQ': 3}
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
        # SOC
        if self._calc_soc is not False:
            input_file += 'calc_soc {}\n'.format(self._calc_soc)

        # CIS variables
        if self._cis_n_roots is not None:
            input_file += 'cis_convergence {}\n'.format(self._cis_convergence)
            input_file += 'cis_n_roots {}\n'.format(self._cis_n_roots)
            input_file += 'cis_singlets {}\n'.format(self._cis_singlets)
            input_file += 'cis_triplets {}\n'.format(self._cis_triplets)
            input_file += 'cis_ampl_anal {}\n'.format(self._cis_ampl_anal)
            input_file += 'loc_cis_ov_separate {}\n'.format(self._loc_cis_ov_separate)
            input_file += 'er_cis_numstate {}\n'.format(self._er_cis_numstate)
            input_file += 'cis_diabath_decompose {}\n'.format(self._cis_diabath_decompose)
            input_file += 'max_cis_cycles {}\n'.format(self._max_cis_cycles)
        # other
        if self._n_frozen_core is not None:
            input_file += 'n_frozen_core {}\n'.format(self._n_frozen_core)
        if self._n_frozen_virt is not None:
            input_file += 'n_frozen_virtual {}\n'.format(self._n_frozen_virt)
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

        input_file += '$end\n'

        # localized diabatization
        if self._localized_diabatization is not None:
            input_file += '$localized_diabatization\nadiabatic states\n'
            input_file += ' '.join(np.array(self._localized_diabatization, dtype=str))
            input_file += '\n$end\n'

        # optimization
        if self._jobtype.lower() == 'opt':
            input_file += 'geom_opt_coords {}\n'.format(self._geom_opt_coords)
            input_file += 'geom_opt_tol_gradient {}\n'.format(self._geom_opt_tol_gradient)
            input_file += 'geom_opt_tol_displacement {}\n'.format(self._geom_opt_tol_displacement)
            input_file += 'geom_opt_tol_energy {}\n'.format(self._geom_opt_tol_energy)
            input_file += 'geom_opt_max_cycles {}\n'.format(self._geom_opt_max_cycles)

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

        return input_file + "\n"

    # Access to properties (only a reduced set should be accessible/editable)
    @property
    def mo_coefficients(self):
        return self._mo_coefficients

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
