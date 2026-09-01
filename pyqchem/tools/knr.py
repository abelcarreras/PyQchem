import numpy as np
import math
import copy
from itertools import combinations_with_replacement
from pyqchem.tools.duschinsky import *
from pyqchem.units import CM1_TO_AU, AMU_TO_ELECTRONMASS, AU_TO_EV, CM_TO_S
from pyqchem.tools.morse_fit import fit_De_per_mode
from functools import cache
from scipy.special import gammaln, hyp2f1

class Knr:
    """
    Implementation of eq. 15 in Phys. Chem. Chem. Phys., 2026,28, 776-789
    """

    def __init__(self,
                 structure_initial, #angstrom
                 structure_final, # angstrom
                 modes_initial, #  read from output, not mass-weighted
                 modes_final, #  read from output, not mass-weighted
                 frequencies_initial, # cm-1
                 frequencies_final, # cm-1
                 derivative_coupling,
                 duschinsky,
                 excitation_energy,
                 ):

        self._modes_initial = NormalModes(structure_initial,
                                          modes_initial,
                                          frequencies_initial,
                                          is_mass_weighted=False)

        self._modes_final = NormalModes(structure_final,
                                        modes_final,
                                        frequencies_final,
                                        is_mass_weighted=False)

        self._modes_initial.trim_negative_frequency_modes()
        self._modes_final.trim_negative_frequency_modes()

        if len(self._modes_initial) != len(self._modes_final):
            n = np.min([len(self._modes_initial), len(self._modes_final)])
            self._modes_initial.trim_modes(list(range(n)))
            self._modes_final.trim_modes(list(range(n)))

        self.derivative_coupling = np.array(derivative_coupling).flatten() # dim 3N
        self._duschinsky = duschinsky
        self._excitation_energy = excitation_energy
        #self.modes_matrix_initial = np.array([np.array(m).flatten() for m in modes_initial])
        #self.modes_matrix_final = np.array([np.array(m).flatten() for m in modes_final])


    def compute_C(self, max_modes=-1):
        """
        Compute non-adiabatic coupling vector C_mu in eV using Atomic Units.
        C_mu = hbar * w_mu * D_mu
        """
        # 1. Obtenemos matrices y vectores base
        L = self._modes_initial.get_displacements().T  # dim (3N, 3N)
        mass_amu = np.array([[m] * 3 for m in self._modes_initial.get_atomic_masses()]).flatten()
        freqs_cm1 = np.array(self._modes_initial.get_frequencies())  # Frecuencias del estado inicial

        # Convertimos todo a a.u.
        mass_au = mass_amu * AMU_TO_ELECTRONMASS
        freqs_au = freqs_cm1 * CM1_TO_AU

        # 2. Proyectamos el acoplamiento a las coordenadas normales (mass-weighted)
        D_MW = np.dot(L, self.derivative_coupling / np.sqrt(mass_au))

        # 3. Convertimos en adimensional
        # D_adim = sqrt(hbar/w_mu) * D_MW. En a.u., hbar = 1.
        # C_mu = w_mu * D_adim = w_mu * D_MW / sqrt(w_mu) = sqrt(w_mu) * D_MW
        D_adim = D_MW / np.sqrt(freqs_au)

        # 4. Compute C_mu = hbar * w_mu * D_adim
        C_mu_au = freqs_au * D_adim

        num_modes = len(freqs_cm1)
        if 0 < max_modes < num_modes:
            C_mu_au = C_mu_au[:max_modes]

        # 5. Return C_mu in eV
        return C_mu_au * AU_TO_EV

    def compute_vibrational(self,
                            max_quanta=2,
                            max_modes=-1,
                            fcf_threshold=1e-8,
                            curvilinear=True,
                            HWHM=0.005):
        """
        Compute vibrational part of the rate

        max_quanta: maximum quanta inside every normal mode
        num_modes: number of modes involved in transition
        """

        hwhm = HWHM

        num_modes = len(self._modes_final)

        if 0 < max_modes < num_modes:
            num_modes = max_modes

        CM1_TO_EV = CM1_TO_AU * AU_TO_EV
        mode_energies = [freq * CM1_TO_EV for freq in self._modes_final.get_frequencies()]
        if curvilinear:
            s_i, s_f = self._duschinsky.get_huang_rys_curvilinear()
        else:
            s_i, s_f = self._duschinsky.get_huang_rys()

        # Slice para quedarnos con los primreos max_modes modos
        mode_energies = mode_energies[:num_modes]
        s_f = s_f[:num_modes]


        ## Decidir quanta maximo por modo
        max_q_per_mode = []
        for i in range(num_modes):
            S = s_f[i]
            allowed_q = 0
            for q in range(max_quanta +1):
                if _recursive_FCF(S,q) >= fcf_threshold:
                    allowed_q = q
            max_q_per_mode.append(allowed_q)

        # Save the remaining allowed quanta
        max_q_remaining = [sum(max_q_per_mode[i:]) for i in range(num_modes)]
        max_q_remaining.append(0)

        # Generador de configuraciones
        def _generar_configs(modo_idx, cuantos_restantes, config_actual):
            if modo_idx == num_modes:
                if cuantos_restantes == 0:
                    yield config_actual.copy()
                return

            if cuantos_restantes > max_q_remaining[modo_idx]:
                return

            limite = min(cuantos_restantes,max_q_per_mode[modo_idx])

            for q in range(limite +1):
                config_actual[modo_idx] = q
                yield from _generar_configs(modo_idx+1, cuantos_restantes-q,config_actual)
                config_actual[modo_idx] = 0

        # Loop on configurations
        vibration = 0.0
        gap_ev = self._excitation_energy * AU_TO_EV
        for total_quanta in range(max_quanta + 1): #
            vector_base = [0] * num_modes

            for configuracion in _generar_configs(0,total_quanta,vector_base):
                energy_conf = sum(quanta * mode_energies[i] for i, quanta in enumerate(configuracion)) # in eV

                fcf_total = 1.0
                for i, quanta in enumerate(configuracion):
                    fcf_total *= _recursive_FCF(s_f[i],quanta)

                x = gap_ev - energy_conf
                densidad = _lorentzian(x,hwhm)

                vibration += fcf_total * densidad

        return vibration

    def compute_FCWD(self,
                     max_quanta=10,
                     fcf_threshold=1e-8,
                     max_modes=-1,
                     curvilinear=True,
                     HWHM=0.005):
        """
        Compute FCWD using the convolution algorithm
        """

        hwhm = HWHM

        # Pasamos a cm-1 para reproducir el articulo
        num_modes = len(self._modes_final)
        if 0 < max_modes < num_modes:
            num_modes = max_modes

        mode_energies_cm1 = self._modes_final.get_frequencies()
        if curvilinear:
            s_i, s_f = self._duschinsky.get_huang_rys_curvilinear()
        else:
            s_i, s_f = self._duschinsky.get_huang_rys()

        # Paso (P): calculo de FCF y energias de configuracion
        modos_data = [] #diccionario

        for i in range(num_modes):
            S = s_f[i]
            freq = mode_energies_cm1[i]

            fcf_list = []
            energy_list = []

            for q in range(max_quanta + 1):
                fcf = _recursive_FCF(S,q)

                if q > S and fcf < fcf_threshold:
                    break

                fcf_list.append(fcf)
                energy_list.append(int(round(q * freq)))

            modos_data.append({
                'fcfs': np.array(fcf_list),
                'energies': np.array(energy_list)
            })

        # Paso A: generar grid de energia (debe cubrir el gap)
        gap_cm1 = self._excitation_energy / CM1_TO_AU
        max_energy = max(25000, int(gap_cm1 * 1.2) + 1000)  # cm-1
        step = 1           # cm-1

        energy_grid = np.arange(0,max_energy + step,step)
        num_puntos = len(energy_grid)

        # Paso B: Iniciar grid
        r_grid = np.zeros(num_puntos)
        r_grid[0] = 1

        # Paso C: convolucion
        for modo in modos_data:
            r_new = np.zeros(num_puntos)

            for fcf, e_shift in zip(modo['fcfs'],modo['energies']):
                if e_shift < num_puntos:
                    r_new[e_shift:] += r_grid[:num_puntos-e_shift] * fcf

            r_grid = r_new

        # Paso D: Average over gaussian function
        gamma = hwhm / AU_TO_EV / CM1_TO_AU

        lorentzian_array = 1 / np.pi * (gamma / ((energy_grid - gap_cm1)**2 + gamma**2))
        vibrational_part = np.sum(r_grid * lorentzian_array) * step / CM1_TO_AU / AU_TO_EV

        return vibrational_part

    def compute_FCWD_anharmonic(self, max_quanta=10, fcf_threshold=1e-8, D_e=30000, max_modes=-1,
                                curvilinear=True, HWHM=0.005, scan_dir='.'):
        """
        Compute FCWD using the convolution algorithms and anharmonic Morse FCF
        """
        hwhm = HWHM

        num_modes = len(self._modes_final)
        if 0 < max_modes < num_modes:
            num_modes = max_modes
        freqs_cm1 = self._modes_final.get_frequencies()
        freqs_au = np.array(freqs_cm1) * CM1_TO_AU

        if curvilinear:
            d_ini, d_fin = self._duschinsky.get_d_vector_curvilinear()
        else:
            d_ini, d_fin = self._duschinsky.get_d_vector()  # desplazamientos en a.u.
        sign_corr = _get_mode_sign_corrections(self._modes_final)
        d_fin = sign_corr * np.array(d_fin)  # signo físico por modo (sin abs: preservar signo)
        Delta = -np.sqrt(freqs_au) * d_fin  # Δ = Q_GS - Q_ES = -d_fin
        s_f = 0.5 * freqs_au * d_fin ** 2  # H-R factors (no cambia, depende de |d_fin|²)

        if D_e == 'fitting':
            D_e = fit_De_per_mode(scan_dir)

        # Paso (P): calculo de FCF y energias de configuracion
        modos_data = []

        for i in range(num_modes):
            freq = freqs_cm1[i]
            S = s_f[i]
            D = Delta[i]

            D_e_i = D_e.get(i, None) if isinstance(D_e, dict) else D_e

            fcf_list = []
            energy_list = []

            if D_e_i is None:
                # Modo sin scan: aprox. armonica
                for q in range(max_quanta + 1):
                    fcf = _recursive_FCF(S, q)
                    if q > S and fcf < fcf_threshold:
                        break
                    fcf_list.append(fcf)
                    energy_list.append(int(round(q * freq)))
            else:
                for q in range(max_quanta + 1):
                    if _morse_j(freq, q, D_e_i) <= 0:
                        #print(f"Warning: mode {i} stops at quanta {q} due to state not being bounded in its Morse Potential")
                        break
                    fcf = _anharmonic_FCF(freq, D, q, D_e_i)
                    if q > S and fcf < fcf_threshold:
                        break
                    fcf_list.append(fcf)
                    energy_list.append(int(round(_morse_energy(freq,q,D_e_i))))

            modos_data.append({
                'fcfs': np.array(fcf_list),
                'energies': np.array(energy_list)
            })

        # Paso A: generar grid de energia (debe cubrir el gap)
        gap_cm1 = self._excitation_energy / CM1_TO_AU
        max_energy = max(25000, int(gap_cm1 * 1.2) + 1000)  # cm-1
        step = 1           # cm-1

        energy_grid = np.arange(0,max_energy + step,step)
        num_puntos = len(energy_grid)

        # Paso B: Iniciar grid
        r_grid = np.zeros(num_puntos)
        r_grid[0] = 1

        # Paso C: convolucion
        for modo in modos_data:
            r_new = np.zeros(num_puntos)

            for fcf, e_shift in zip(modo['fcfs'],modo['energies']):
                if e_shift < num_puntos:
                    r_new[e_shift:] += r_grid[:num_puntos-e_shift] * fcf

            r_grid = r_new

        # Paso D: Average over gaussian function
        gamma = hwhm / AU_TO_EV / CM1_TO_AU

        lorentzian_array = 1 / np.pi * (gamma / ((energy_grid - gap_cm1) ** 2 + gamma ** 2))
        vibrational_part = np.sum(r_grid * lorentzian_array) * step / CM1_TO_AU / AU_TO_EV

        return vibrational_part

    def print_fcf(self,
                  max_quanta=10,
                  fcf_threshold=1e-3,
                  D_e=30000,
                  max_modes=-1,
                  curvilinear=True,
                  output_file=None):
        """
        Print harmonic and anharmonic FC factors above fcf_threshold for each mode.
        If output_file is given, writes to that file.
        """
        freqs_cm1 = self._modes_final.get_frequencies()
        if curvilinear:
            s_i, s_f = self._duschinsky.get_huang_rys_curvilinear()
            d_ini, d_fin = self._duschinsky.get_d_vector_curvilinear()
        else:
            s_i, s_f = self._duschinsky.get_huang_rys()
            d_ini, d_fin = self._duschinsky.get_d_vector()
        sign_corr = _get_mode_sign_corrections(self._modes_final)
        d_fin_signed = sign_corr * np.array(d_fin)
        freqs_au = np.array(freqs_cm1) * CM1_TO_AU
        Delta = -np.sqrt(freqs_au) * d_fin_signed

        num_modes = len(freqs_cm1)
        if 0 < max_modes < num_modes:
            num_modes = max_modes

        header = (f"\n{'Mode':>5}  {'freq(cm-1)':>10}  {'S':>7}  {'Delta':>8}  "
                  f"{'q':>3}  {'FCF_harm':>12}  {'FCF_anh':>12}\n" + "-" * 70)

        lines = [header]
        for i in range(num_modes):
            freq = freqs_cm1[i]
            S = s_f[i]
            D = Delta[i]

            harm_rows = []
            anh_rows = []
            for q in range(max_quanta + 1):
                fh = _recursive_FCF(S, q)
                if fh >= fcf_threshold:
                    harm_rows.append((q, fh))

                if _morse_j(freq, q, D_e) > 0:
                    fa = _anharmonic_FCF(freq, D, q, D_e)
                    if fa >= fcf_threshold:
                        anh_rows.append((q, fa))

            all_q = sorted(set(q for q, _ in harm_rows) | set(q for q, _ in anh_rows))
            if not all_q:
                continue

            harm_dict = dict(harm_rows)
            anh_dict = dict(anh_rows)
            first = True
            for q in all_q:
                fh = harm_dict.get(q, 0.0)
                fa = anh_dict.get(q, 0.0)
                if fh < fcf_threshold and fa < fcf_threshold:
                    continue
                if first:
                    lines.append(f"{i+1:>5}  {freq:>10.1f}  {S:>7.4f}  {D:>8.4f}  "
                                 f"{q:>3}  {fh:>12.6f}  {fa:>12.6f}")
                    first = False
                else:
                    lines.append(f"{'':>5}  {'':>10}  {'':>7}  {'':>8}  "
                                 f"{q:>3}  {fh:>12.6f}  {fa:>12.6f}")
            lines.append("")

        text = "\n".join(lines)
        if output_file is not None:
            with open(output_file, 'w') as f:
                f.write(text + "\n")
        else:
            print(text)

    def compute_knr(self,
                    max_quanta=2,
                    max_modes=-1,
                    fcf_threshold=1e-8,
                    anharmonic=False,
                    D_e=30000,
                    curvilinear=True,
                    HWHM=0.005,
                    scan_dir='.'):
        """
        Compute knr as the product of its two parts:
        """

        hbar = 6.582119569E-16 # eV * s
        # return (np.pi / hbar ) * np.sum(self.compute_C()**2) * self.compute_vibrational(max_quanta, max_modes,fcf_threshold)
        if anharmonic:
            fcwd = self.compute_FCWD_anharmonic(max_quanta, fcf_threshold, D_e, max_modes,curvilinear, HWHM, scan_dir)
        else:
            fcwd = self.compute_FCWD(max_quanta, fcf_threshold, max_modes,curvilinear, HWHM)
        return (np.pi / hbar ) * np.sum(self.compute_C(max_modes)**2) * fcwd


@cache
def _recursive_FCF(S,n):
    """
    Compute the FCF using the relation FCF(S,n+1) = S/(n+1) * F(S,n)
    """
    if n == 0:
        return np.exp(-S)
    else:
        return S/(n) * _recursive_FCF(S,n-1)

def _compute_FCF(S,n):
    """
    Compute F-C factor between 0 quanta S1 and n quanta S0, for a mode with Huang-Rhys S
    """
    return (S**n)/ math.factorial(n) * np.exp(-S)

def _lorentzian(E, sigma):
    """
    Devuelve el valor de una Lorentziana normalizada en el punto E con HWHM sigma
    """

    return (1.0 / np.pi) * (sigma / (E**2 + sigma**2))


def _gaussian(E, sigma):
    """
    Devuelve el valor de una función Gaussiana normalizada en el punto E.
    sigma: Desviación estándar de la distribución.
    """
    prefactor = 1.0 / (sigma * np.sqrt(2.0 * np.pi))
    exponencial = np.exp(-0.5 * (E / sigma) ** 2)

    return prefactor * exponencial

def _get_bonds(coords_angstrom, cutoff=1.7):
    """
    Detecta enlaces por distancia entre átomos.
    Devuelve lista de pares (i, j) con i < j.
    """
    N = len(coords_angstrom)
    bonds = []
    for i in range(N):
        for j in range(i + 1, N):
            dist = np.linalg.norm(coords_angstrom[i] - coords_angstrom[j])
            if dist < cutoff:
                bonds.append((i, j))
    return bonds


def _build_B_matrix(coords_angstrom, bonds):
    """
    Construye la Wilson B-matrix para longitudes de enlace.
    Dimensiones (N_bonds, 3*N_atoms).
    """
    N = len(coords_angstrom)
    B = np.zeros((len(bonds), 3 * N))
    for row, (i, j) in enumerate(bonds):
        vec = coords_angstrom[j] - coords_angstrom[i]
        e = vec / np.linalg.norm(vec)   # vector unitario i -> j
        B[row, 3*i : 3*i+3] = -e
        B[row, 3*j : 3*j+3] = +e
    return B


def _get_mode_sign_corrections(modes_final):
    """
    Determina la corrección de signo (+1 o -1) según el criterio de la SI Sec. S2.1
    """
    coords = np.array(modes_final.get_coordinates())   # (N_atoms, 3), Angstrom
    L = np.array(modes_final.get_displacements())      # (3N, N_modes)

    bonds = _get_bonds(coords)
    B = _build_B_matrix(coords, bonds)   # (N_bonds, 3N)

    # Cambio en cada longitud de enlace para cada modo: (N_bonds, N_modes)
    # B es (N_bonds, 3N), L es (3N, N_modes)
    delta_bonds = B @ L

    # Signo del modo = signo del cambio de enlace de mayor magnitud
    N_modes = L.shape[1]
    signs = np.ones(N_modes)
    for mu in range(N_modes):
        col = delta_bonds[:, mu]
        idx_max = np.argmax(np.abs(col))
        s = np.sign(col[idx_max])
        if s != 0:
            signs[mu] = s

    return signs

def _morse_j(freq_cm1, v, D_e =30000):
    # eq. 24 in SI
    return 4 * D_e / freq_cm1 - 2*v -1

def _morse_energy(freq_cm1, v, D_e =30000):
    # eq. 27 in SI
    
    # yield energy in cm-1
    xi = freq_cm1/(4 * D_e)
    return freq_cm1 * v - xi * freq_cm1 * (v + 0.5)**2 + 0.25 * xi * freq_cm1

def _anharmonic_FCF(freq_cm1, Delta, v, D_e =30000):
    # eq. 28 in SI
    # Delta = sqrt(omega_mu) * d_fin

    j0 = _morse_j(freq_cm1, 0, D_e)
    jv = _morse_j(freq_cm1, v, D_e)

    # Morse stifness a_mu
    a = np.sqrt(freq_cm1/(2*D_e))
    d = np.exp(-a*Delta) #eq.29

    # Prefactor de normalizacion 1/a * N(mu,0) * N(mu,nu)
    # N dado en eq.26
    log_norm = 0.5 * (np.log(j0)+np.log(jv)+ gammaln(v+1)-gammaln(j0+1)-gammaln(jv+v+1))

    # Argumentos de I_v
    a_p = (d+1)/2
    b_p = (jv +j0)/2 -1
    c_p = jv

    # Prefactor I_v
    log_I_pre = gammaln(1+b_p) + gammaln(c_p+v+1)-gammaln(v+1)-gammaln(c_p+1)

    # Polinomio hipergeometrico 2F1
    h = hyp2f1(1+b_p, -v, 1+c_p, 1/a_p)

    # log FCF (sin el factor hipergeometrico)
    log_F = (log_norm + (j0/2)*np.log(d)+log_I_pre+(-1-b_p)*np.log(a_p))

    # return FCF
    return (np.exp(log_F) * abs(h))**2

def _align_structure(source_output, target_output):
    """
    Rotate (Kabsch) the structure and normal modes of source_output onto the frame of target_output.
    """
    coor_s = np.array(source_output['structure'].get_coordinates())
    coor_t = np.array(target_output['structure'].get_coordinates())
    com_s, com_t = coor_s.mean(axis=0), coor_t.mean(axis=0)

    H = (coor_s - com_s).T @ (coor_t - com_t)
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T

    aligned = copy.deepcopy(source_output)
    aligned['structure'].set_coordinates(((R @ (coor_s - com_s).T).T + com_t).tolist())
    for mode in aligned['modes']:
        mode['displacement'] = (R @ np.array(mode['displacement']).T).T.tolist()
    return aligned


def get_knr(origin_frequency_output, target_frequency_output,
            origin_energy_output, target_energy_output,
            derivative_coupling_output):
    """
    build knr instance object from frequency parser dictionary and derivative coupling calculation

    :param origin_frequency_output: frequency parsed output of origin state (typically ground state)
    :param target_frequency_output: frequency parsed output of target state (typically excited state)
    :param derivative_coupling_output: calculation of the derivative coupling between states
    :return: knr object
    """

    # origin keeps its own frame (derivative_coupling lives there); target gets aligned onto it
    target_frequency_output = _align_structure(target_frequency_output, origin_frequency_output)

    return Knr(structure_initial=origin_frequency_output['structure'],
               structure_final=target_frequency_output['structure'],
               modes_initial=[mode['displacement'] for mode in origin_frequency_output['modes']],
               modes_final=[mode['displacement'] for mode in target_frequency_output['modes']],
               frequencies_initial=[mode['frequency'] for mode in origin_frequency_output['modes']],
               frequencies_final=[mode['frequency'] for mode in target_frequency_output['modes']],
               derivative_coupling=derivative_coupling_output['derivative_coupling'],
               duschinsky=get_duschinsky(origin_frequency_output, target_frequency_output),
               excitation_energy=origin_energy_output['energy']-target_energy_output['energy'],
               )