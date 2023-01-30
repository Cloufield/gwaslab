import pandas as pd
import numpy as np
from gwaslab.Log import Log

def adjust_text_position(positions, yspan, repel_force=0.01, max_iter=100,log=Log(),verbose=True):

    # check the number of variants to annotate 
    if repel_force>0:
        if 1/(repel_force*2 +0.01) < len(positions):
            if verbose: log.write(" -Too many variants to annotate; maybe it is better to reduce the number of variants")
    else:
        if len(positions)>30:
            if verbose: log.write(" -Too many variants to annotate; maybe it is better to reduce the number of variants")

    # calculate the steps
    step = int(yspan*repel_force)

    # start iteration
    for i in range(max_iter):
        # check overlap
        index = check_overlap(positions, step)

        # fix overlap if needed
        if index == len(positions)+1:
            # if no overlap
            return  np.floor(pd.to_numeric(positions, errors='coerce')).astype('Int64').copy()
        else:
            # if overlap exists
            move_position_from_center(positions, index, step)
    
    # when reaching maximum iteration, return anyway
    if verbose: log.write(" -Reaching maximum iteration: {}; Skipping...".format(max_iter))
    return np.floor(pd.to_numeric(positions, errors='coerce')).astype('Int64').copy()

def check_overlap(positions,step):
    # check overlap by walkthrough
    for i in range(1,len(positions)):
        if positions[i] - positions[i-1]< step:
            # if overlap: return the index
            return i
    #if no overlap , return maxmum length + 1
    return len(positions)+1

def move_position_from_center(positions, center_i, step):
    # left side
    #print("center_i",center_i)
    #print("before",positions)
    
    # if not pass left bound
    if positions[center_i-1] - step//2 > 0:
        positions[center_i-1] = positions[center_i-1] - step//2
    
    # if not the second element
    if positions[center_i-2] - step//2 > 0:
        for i in range(center_i-2, 0, -1):
            if abs(positions[i] - positions[i+1]) < step:
                 positions[i] = positions[i] - step//2
            else:
                break
    
    # right side            
    positions[center_i] = positions[center_i] + step//2
    for i in range(center_i+1, len(positions)):
        if abs(positions[i] - positions[i-1]) < step:
             positions[i] = positions[i] + step//2
        else:
            break
    #print("after",positions)
    return positions
    
    