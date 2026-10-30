import re
import numpy as np
import glob
from scipy.optimize import curve_fit, minimize_scalar
from pyqchem.units import CM1_TO_AU

def load_morse_scan(filepath):
    header = {}
    Q, E = [], []
    with open(filepath) as f:
        first_line = f.readline()
        for key in ('mode_idx', 'freq_cm1', 'S', 'elong_sign'):
            m = re.search(rf'{key}=([-\d.eE+]+)', first_line)
            if m:
                header[key] = float(m.group(1))
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            q, e = line.split()
            Q.append(float(q))
            E.append(float(e))

    header['mode_idx'] = int(header['mode_idx'])
    return header, np.array(Q), np.array(E)

def morse_model(Q, D_e, Q0, V0, freq_cm1):
    a = np.sqrt(freq_cm1 / (2 * D_e))
    return V0 + D_e * (1 - np.exp(-a * (Q - Q0))) ** 2 * CM1_TO_AU

def fit_morse_De(Q, E_hartree, freq_cm1):
    Q = np.asarray(Q)
    E_hartree = np.asarray(E_hartree)

    i0 = np.argmin(np.abs(Q))
    V0 = E_hartree[i0]   # energia en Q=0, fija

    def model_log(Q, log_D_e):
        D_e = 10 ** log_D_e
        return morse_model(Q, D_e, 0.0, V0, freq_cm1)

    def ssr(log_D_e):
        residuals = E_hartree - model_log(Q, log_D_e)
        return np.sum(residuals ** 2)

    # barrido grueso en escala log para localizar la zona del minimo real
    # (el gradiente en D_e lineal se vuelve casi plano para D_e grande y
    # curve_fit para demasiado pronto si arranca directamente ahi)
    log_grid = np.linspace(2, 9, 71)  # D_e de 100 a 1e9 cm-1
    ssr_grid = [ssr(lg) for lg in log_grid]
    log_D_e0 = log_grid[np.argmin(ssr_grid)]

    popt, pcov = curve_fit(model_log, Q, E_hartree, p0=[log_D_e0])
    D_e_fit = 10 ** popt[0]

    residuals = E_hartree - model_log(Q, popt[0])
    ss_res = np.sum(residuals ** 2)
    ss_tot = np.sum((E_hartree - E_hartree.mean()) ** 2)
    r_squared = 1 - ss_res / ss_tot

    return {
        'D_e': D_e_fit,
        'V0': V0,
        'r_squared': r_squared,
    }

def fit_morse_scan_file(filepath):
    header, Q, E = load_morse_scan(filepath)
    fit = fit_morse_De(Q, E, header['freq_cm1'])
    return header, fit

def fit_De_per_mode(scan_dir):
    D_e_per_mode = {}
    for filepath in glob.glob(f'{scan_dir}/morse_scan_mode*.dat'):
        header, fit = fit_morse_scan_file(filepath)
        D_e_per_mode[header['mode_idx']] = fit['D_e']
    return D_e_per_mode

def _fit_De_fixed_V0(Q, E_hartree, V0, freq_cm1):
    """
    Least-squares D_e (cm-1) of a Morse potential with V0 and frequency fixed (log scale search)
    """
    ssr = lambda lg: np.sum((E_hartree - morse_model(Q, 10 ** lg, 0.0, V0, freq_cm1)) ** 2)
    grid = np.linspace(1, 20, 381)   # D_e de 10 a 1e20 cm-1
    g0 = grid[np.argmin([ssr(g) for g in grid])]
    if grid[0] < g0 < grid[-1]:
        g0 = minimize_scalar(ssr, bracket=(g0 - 0.05, g0, g0 + 0.05)).x
    return 10 ** g0


def fit_morse_De_oriented(Q, E_hartree, freq_cm1, Q_max=None, V0=None, De_cap=1e8, min_explained=0.5):
    """
    Fit the D_e of a Morse potential choosing the orientation of Q from the scan itself.
    """
    Q = np.asarray(Q, dtype=float)
    E_hartree = np.asarray(E_hartree, dtype=float)
    if Q_max is None:
        Q_max = min(-Q.min(), Q.max())
    sel = np.abs(Q) <= Q_max + 1e-8
    Q, E = Q[sel], E_hartree[sel]
    if V0 is None:
        V0 = E[np.argmin(np.abs(Q))]

    e_harm = V0 + 0.5 * freq_cm1 * Q ** 2 * CM1_TO_AU
    dE = (E - e_harm) / CM1_TO_AU
    c3, c4 = np.linalg.lstsq(np.vstack([Q ** 3, Q ** 4]).T, dE, rcond=None)[0]
    sign = 1 if c3 < 0 else -1

    D_e = _fit_De_fixed_V0(sign * Q, E, V0, freq_cm1)
    e_morse = morse_model(sign * Q, D_e, 0.0, V0, freq_cm1)
    ssr_harm = np.sum((E - e_harm) ** 2)
    ssr_morse = np.sum((E - e_morse) ** 2)
    explained = 1 - ssr_morse / ssr_harm if ssr_harm > 0 else 0.0

    if D_e >= De_cap:
        status = 'harmonic'
    elif explained < min_explained:
        status = 'non_morse'
    else:
        status = 'morse'

    return {'D_e': D_e if status == 'morse' else None,
            'D_e_fit': D_e,
            'sign': sign,
            'status': status,
            'c3': c3,
            'c4': c4,
            'explained': explained,
            'r_squared': 1 - ssr_morse / np.sum((E - E.mean()) ** 2),
            'V0': V0,
            'Q_max': Q_max}


def fit_morse_scan_file_oriented(filepath, **kwargs):
    header, Q, E = load_morse_scan(filepath)
    fit = fit_morse_De_oriented(Q, E, header['freq_cm1'], **kwargs)
    return header, fit


def fit_De_and_signs_per_mode(scan_dir, **kwargs):
    """
    Fit all morse_scan_mode*.dat files in scan_dir with fit_morse_De_oriented

    :return: D_e_per_mode {mode: D_e} (only modes with status 'morse'; the rest are harmonic) and
             signs_per_mode {mode: +1/-1} (orientation of Q, relative to the scan files)
    """
    D_e_per_mode, signs_per_mode = {}, {}
    for filepath in glob.glob('{}/morse_scan_mode*.dat'.format(scan_dir)):
        header, fit = fit_morse_scan_file_oriented(filepath, **kwargs)
        if fit['status'] == 'morse':
            D_e_per_mode[header['mode_idx']] = fit['D_e']
            signs_per_mode[header['mode_idx']] = fit['sign']
    return D_e_per_mode, signs_per_mode
