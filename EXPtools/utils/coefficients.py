"""
Functionality to handle coefficients
"""

import numpy as np
from EXPtools.utils.indexing import inverse_I, list_states_range


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

def truncate_expansion(original_coefficients, nmax, lmax):
    '''truncate_expansion removes terms higher than nmax, lmax

    Parameters
    ----------
    original_coefficients : pyEXP coefficients
        coefficient set to be truncated
    nmax : _type_
        maximum n order to keep
    lmax : _type_
        maximum l order to keep

    Returns
    -------
    pyEXP coefficients
        truncated coefficient set
    '''

    coefs_shape = original_coefficients.getAllCoefs().shape
    nmax_original = coefs_shape[1]
    lmax_original = int(np.sqrt(9./4. + 2*(coefs_shape[0] - 1)) - 3./2.)

    if (nmax > nmax_original) or (lmax > lmax_original):
        raise ValueError("New max orders must be lower than the size of the original expansion!")

    n_remove, l_remove, m_remove = list_states_range(nmax, lmax, nmax_original, lmax_original)
    
    return remove_terms(original_coefficients, n_remove, l_remove, m_remove)

# Plot phase of coefficients! 

# 
