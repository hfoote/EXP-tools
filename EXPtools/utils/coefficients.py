"""
Functionality to handle coefficients
"""

import numpy as np
from EXPtools.utils.indexing import list_states_range, find_max_order


def remove_terms(original_coefficients, n, l, m, verbose=False):
    """
    Remove coefficients 
    """ 
    assert len(l) == len(m) == len(n)
    copy_coeffcients = original_coefficients.deepcopy()
    coefs_matrix = original_coefficients.getAllCoefs()
    t_snaps = original_coefficients.Times()
    
    for t in range(len(t_snaps)):
        for i in range(len(l)):
            lm_idx = int(l[i]*(l[i]+1) / 2) + m[i]
            if verbose:
                print(n[i],l[i],m[i], lm_idx)
            try: coefs_matrix[lm_idx, n[i], t] = np.complex128(0)
            except IndexError: continue
        copy_coeffcients.setMatrix(mat=coefs_matrix[:,:, t], time=t_snaps[t])
    
    return copy_coeffcients
    
def reorder_nlm(coefficients, nmax, lmax):
    """
    return coefficients in order (n, l, l+1)
    """
    new_order = np.zeros((2, nmax+1, lmax+1, lmax+1))
    coefs_matrix = coefficients
    for n in range(nmax+1):
        print(n)
        for l in range(lmax+1):
            for m in range(l+1):
                lm_idx = int(l*(l+1) / 2) + m
                new_order[0][n][l][m] = coefs_matrix[n, lm_idx].real
                new_order[1][n][l][m] = coefs_matrix[n, lm_idx].imag
    return new_order

def truncate_expansion(original_coefficients, nmax=None, lmax=None):
    '''truncate_expansion removes terms higher than nmax, lmax

    Parameters
    ----------
    original_coefficients : pyEXP coefficients
        coefficient set to be truncated
    nmax : None or int, optional
        maximum n order to keep. If None, keeps all orders
    lmax : None or int, optional
        maximum l order to keep. If None, keeps all orders

    Returns
    -------
    pyEXP coefficients
        truncated coefficient set
    '''

    nmax_original, lmax_original = find_max_order(original_coefficients)

    if (nmax > nmax_original) or (nmax == None):
        print(f"Keeping all n up to original nmax = {nmax_original}")
        nmax = nmax_original

    if (lmax > nmax_original) or (lmax == None):
        print(f"Keeping all l up to original lmax = {lmax_original}")
        lmax = lmax_original

    n_remove, l_remove, m_remove = list_states_range(nmax, lmax, nmax_original, lmax_original)
    
    return remove_terms(original_coefficients, n_remove, l_remove, m_remove)

# Plot phase of coefficients! 

# 
