import re
import numpy as np
import glob
from scipy.optimize import curve_fit
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