
import numpy as np
import pandas as pd

def snphwe(obs_hets, obs_hom1, obs_hom2):
    # Convert cpp code from (Jeremy McRae) to python  
    # https://github.com/jeremymcrae/snphwe/blob/master/src/snp_hwe.cpp
    #/* (original comments)
    #// This code implements an exact SNP test of Hardy-Weinberg Equilibrium as
    #// described in Wigginton, JE, Cutler, DJ, and Abecasis, GR (2005) A Note on
    #// Exact Tests of Hardy-Weinberg Equilibrium. AJHG 76: 887-893
    #//
    #// Written by Jan Wigginton
    #*/
    
    obs_homr = min(obs_hom1, obs_hom2)
    obs_homc = max(obs_hom1, obs_hom2)
    
    rare = 2 * obs_homr + obs_hets
    genotypes = obs_hets + obs_homc + obs_homr
    
    probs = np.array([0.0 for i in range(rare +1)])
    
    mid = rare * (2 * genotypes - rare) // (2 * genotypes)
    
    if mid % 2 != rare%2:
        mid += 1
        
    probs[mid] = 1.0
    
    sum_p = 1 #probs[mid]
    curr_homr = (rare - mid) // 2
    curr_homc = genotypes - mid - curr_homr
    
    
    for curr_hets in range(mid, 1, -2):
        probs[curr_hets - 2] = probs[curr_hets] * curr_hets * (curr_hets - 1.0)/ (4.0 * (curr_homr + 1.0) * (curr_homc + 1.0))
        sum_p+= probs[curr_hets - 2]
        curr_homr += 1
        curr_homc += 1
    
    curr_homr = (rare - mid) // 2
    curr_homc = genotypes - mid - curr_homr
    
    for curr_hets in range(mid, rare-1, 2):
        probs[curr_hets + 2] = probs[curr_hets] * 4.0 * curr_homr * curr_homc/ ((curr_hets + 2.0) * (curr_hets + 1.0))
        sum_p += probs[curr_hets + 2]
        curr_homr -= 1
        curr_homc -= 1
    
    target = probs[obs_hets]
    p_hwe = 0.0
    
    for p in probs:
        if p <= target :
            p_hwe += p / sum_p  
    
    return min(p_hwe,1)