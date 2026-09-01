import numpy as np
from collections import defaultdict
from pyqchem.units import ANGSTROM_TO_AU, AMU_TO_ELECTRONMASS

def get_bonds(coords, cutoff=1.7):
    """
    Builds a matrix that yields indexes
    """
    coords = np.array(coords)
    N = len(coords)
    bonds = []
    for i in range(N):
        for j in range(i + 1, N):
            if np.linalg.norm(coords[i] - coords[j]) < cutoff:
                bonds.append((i, j))
    return bonds

def measure_bond(coords, i, j):
    """
    Computes bond distance between two atoms
    """
    return np.linalg.norm(coords[j] - coords[i])

def measure_angle(coords, i, j, k):
    """
    Computes angle between three atoms in radians
    """
    v1 = coords[i] - coords[j]
    v2 = coords[k] - coords[j]
    cos_t = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
    return np.arccos(np.clip(cos_t, -1.0, 1.0)) # safe against numerical noise

def measure_dihedral(coords, i, j, k, l):
    """
    Computes dihedral angle between fourla atoms in radians
    """
    b0 = coords[i] - coords[j]
    b1 = coords[k] - coords[j]
    b2 = coords[l] - coords[k]
    b1n = b1 / np.linalg.norm(b1)
    v = b0 - np.dot(b0, b1n) * b1n
    w = b2 - np.dot(b2, b1n) * b1n
    x = np.dot(v, w) # |v||w|cos(phi)
    y = np.dot(np.cross(b1n, v), w) # |v||w|·sin(phi),
    return np.arctan2(y, x) # return phi between -pi and pi

def get_angles(coords, bonds, linear_cutoff_deg=179.0):
    """
    Build matrix with angles. It excludes angles with angle higher than the cutoff.
    """
    adj = defaultdict(set) # Adyacent atoms dictionary, adj[j] yields atoms connected to atom j
    for i, j in bonds:
        adj[i].add(j)
        adj[j].add(i)

    linear_cutoff = np.radians(linear_cutoff_deg)
    angles = []
    for j in adj:
        neighbors = sorted(adj[j])
        for a in range(len(neighbors)):
            for b in range(a + 1, len(neighbors)):
                i, k = neighbors[a], neighbors[b]
                theta = measure_angle(coords, i, j, k)
                if theta < linear_cutoff:
                    angles.append((i, j, k))
    return angles

def get_dihedrals(bonds, angles):
    """
    Builds matrix with dihedral angles
    """
    adj = defaultdict(set)
    for i, j in bonds:
        adj[i].add(j)
        adj[j].add(i)

    valid_angle_pairs = {(a, b, c) for a, b, c in angles} | {(c, b, a) for a, b, c in angles}

    dihedrals = []
    for j, k in bonds:
        for i in adj[j] - {k}:
            if (i, j, k) not in valid_angle_pairs:
                continue
            for l in adj[k] - {j, i}:
                if (j, k, l) not in valid_angle_pairs:
                    continue
                dihedrals.append((i, j, k, l))
    return dihedrals

def compute_internal_coordinates(coords, bonds, angles, dihedrals):
    """
    Builds vector with the values of the bonds, angles and dihedrals.
    This is vector z in J. Chem. Phys. 115, 9103–9109 (2001)
    """
    q_bonds = [measure_bond(coords, *b) for b in bonds]
    q_angles = [measure_angle(coords, *a) for a in angles]
    q_dihedrals = [measure_dihedral(coords, *d) for d in dihedrals]
    return np.concatenate([q_bonds, q_angles, q_dihedrals])

def _b_row_bond(coords, i, j, n_atoms):
    """
    Derivative of bond distance r_ij with respect to cartesian coordinates
    """
    row = np.zeros(3 * n_atoms)
    vec = coords[j] - coords[i]
    e = vec / np.linalg.norm(vec)
    row[3*i:3*i+3] = -e
    row[3*j:3*j+3] = e
    return row


def _b_row_angle(coords, i, j, k, n_atoms):
    """
    Derivative of angle ijk with respect to cartesian coordinates
    """
    row = np.zeros(3 * n_atoms)
    r_ji = coords[i] - coords[j]
    r_jk = coords[k] - coords[j]
    d_ji, d_jk = np.linalg.norm(r_ji), np.linalg.norm(r_jk)
    e_ji, e_jk = r_ji / d_ji, r_jk / d_jk
    cos_t = np.dot(e_ji, e_jk)
    sin_t = np.sqrt(max(1 - cos_t**2, 1e-8))  # evita division por 0 cerca de angulo lineal

    d_i = (cos_t * e_ji - e_jk) / (d_ji * sin_t)
    d_k = (cos_t * e_jk - e_ji) / (d_jk * sin_t)
    d_j = -(d_i + d_k)

    row[3*i:3*i+3] = d_i
    row[3*j:3*j+3] = d_j
    row[3*k:3*k+3] = d_k
    return row


def _b_row_dihedral(coords, i, j, k, l, n_atoms):
    """
    Derivative of dihedral angle ijkl with respect to cartesian coordinates
    """
    row = np.zeros(3 * n_atoms)
    b1 = coords[j] - coords[i]
    b2 = coords[k] - coords[j]
    b3 = coords[l] - coords[k]
    n1 = np.cross(b1, b2)
    n2 = np.cross(b2, b3)
    b2_norm = np.linalg.norm(b2)

    d_i = -(b2_norm / np.dot(n1, n1)) * n1
    d_l = (b2_norm / np.dot(n2, n2)) * n2
    d_j = -d_i - (np.dot(b1, b2) / np.dot(b2, b2)) * d_i + (np.dot(b3, b2) / np.dot(b2, b2)) * d_l
    d_k = -d_l + (np.dot(b1, b2) / np.dot(b2, b2)) * d_i - (np.dot(b3, b2) / np.dot(b2, b2)) * d_l

    row[3*i:3*i+3] = d_i
    row[3*j:3*j+3] = d_j
    row[3*k:3*k+3] = d_k
    row[3*l:3*l+3] = d_l
    return row

def compute_b_matrix(coords, bonds, angles, dihedrals):
    """
    Wilson B Matrix (n x 3N):
    row = derivative of internal coordinate with respect to cartesian coordinates
    """
    n_atoms = len(coords)
    rows = [_b_row_bond(coords, i, j, n_atoms) for i, j in bonds]
    rows += [_b_row_angle(coords, i, j, k, n_atoms) for i, j, k in angles]
    rows += [_b_row_dihedral(coords, i, j, k, l, n_atoms) for i, j, k, l in dihedrals]
    return np.array(rows)

def compute_g_matrix(B, masses):
    """
    Wilson G matrix: G = B * m^-1 * B^T   (Eq. 4, Reimers)
    """
    mass_vec = np.repeat(masses, 3)
    return B @ np.diag(1.0 / mass_vec) @ B.T

def symmetric_inverse_sqrt(M, tol=1e-8):
    """
    Inverse square root of a symmetric positive semi-definite matrix,
    M^(-1/2) = V * diag(1/sqrt(eigval)) * V^T
    """
    eigvals, eigvecs = np.linalg.eigh(M)
    eigvals = np.clip(eigvals, tol, None)
    return eigvecs @ np.diag(1.0 / np.sqrt(eigvals)) @ eigvecs.T

def get_nonredundant_projector(G, n_modes):
    """
    Matrix 'a' (Eq. 3, Reimers): the n_modes eigenvectors of G with nonzero eigenvalue
    """
    eigvals, eigvecs = np.linalg.eigh(G)
    order = np.argsort(eigvals)[::-1]
    return eigvecs[:, order][:, :n_modes]

def compute_nonredundant_b(B, masses, n_modes):
    """
    B'' (Eq. 8, Reimers): orthonormal projection of the redundant internal
    coordinates onto the vibrational subspace (n_v x 3N)
    """
    G = compute_g_matrix(B, masses)
    a = get_nonredundant_projector(G, n_modes)

    B_prime = a.T @ B                                          # Eq. 5
    mass_vec = np.repeat(masses, 3)
    G_prime = B_prime @ np.diag(1.0 / mass_vec) @ B_prime.T    # Eq. 7
    G_prime_inv_sqrt = symmetric_inverse_sqrt(G_prime)          # for Eq. 6/8

    return a, G_prime_inv_sqrt, B_prime

def curvilinear_displacement(coords_ref, L_ref, masses_ref, coords_other, bonds, angles, dihedrals):
    """
    Curvilinear-corrected displacement (Eq. 20, Reimers)
    """
    coords_ref_au = np.array(coords_ref) * ANGSTROM_TO_AU
    coords_other_au = np.array(coords_other) * ANGSTROM_TO_AU
    masses_au = np.array(masses_ref) * AMU_TO_ELECTRONMASS

    n_modes = L_ref.shape[1]
    B = compute_b_matrix(coords_ref_au, bonds, angles, dihedrals)
    a, G_prime_inv_sqrt, B_prime = compute_nonredundant_b(B, masses_au, n_modes)
    B_dprime = G_prime_inv_sqrt @ B_prime                            # B'' , Eq. 8

    mass_vec = np.repeat(masses_au, 3)
    c_dprime = B_dprime @ np.diag(1.0 / np.sqrt(mass_vec)) @ L_ref   # c'' , Eq. 11

    z_ref = compute_internal_coordinates(coords_ref_au, bonds, angles, dihedrals)
    z_other = compute_internal_coordinates(coords_other_au, bonds, angles, dihedrals)
    dz = z_other - z_ref

    n_dih = len(dihedrals)
    if n_dih > 0:
        dz[-n_dih:] = (dz[-n_dih:] + np.pi) % (2 * np.pi) - np.pi   # wrap dihedral jumps

    d = c_dprime.T @ G_prime_inv_sqrt @ a.T @ dz    # Eq. 20
    return d