import pandas as pd
import numpy as np
from gwaslab.info.g_Log import Log

def adjust_text_position(positions, yspan, repel_force=0.01, max_iter=100,amode="int",log=Log(),verbose=True, min_factor=None):
    # check the number of variants to annotate 
    #if repel_force>0:
    #    if 1/(repel_force*2 +0.01) < len(positions):
    #        log.write(" -Too many variants to annotate; maybe it is better to reduce the number of variants")
    #else:
    if len(positions)>30:
        log.write(" -Too many variants to annotate; maybe it is better to reduce the number of variants",verbose=verbose)

    # calculate the steps
    if amode=="int":
        step = int(yspan*repel_force) 
    elif amode=="log":
        if min_factor is None:
            min_factor = np.min(positions)
        #(1, max) -> (0, log(max)))
        positions = np.log2(positions/min_factor)
        step = max(positions)*repel_force 

    else:
        step = yspan*repel_force 

    # start iteration
    for i in range(max_iter):
        # check overlap
        index = check_overlap(positions, step)

        # fix overlap if needed
        if index == len(positions)+1:
            # if no overlap
            if amode=="int":
                return  np.floor(pd.to_numeric(positions, errors='coerce')).astype('Int64').copy()
            elif amode=="log":
                
                return  np.power(2, pd.to_numeric(positions, errors='coerce'))* min_factor
            else:
                return  pd.to_numeric(positions, errors='coerce')
        else:
            # if overlap exists
            if amode=="int":
                move_position_from_center_int(positions, index, step)
            else:
                move_position_from_center_float(positions, index, step)
    
    # when reaching maximum iteration, return anyway
    log.write(" -Reaching maximum iteration: {}; Skipping...".format(max_iter),verbose=verbose)
    if amode=="int":
        return np.floor(pd.to_numeric(positions, errors='coerce')).astype('Int64').copy()
    elif amode=="log":
        return np.exp(pd.to_numeric(positions, errors='coerce')) * min_factor
    else:
        return  pd.to_numeric(positions, errors='coerce')

def check_overlap(positions,step):
    # check overlap by walkthrough
    for i in range(1,len(positions)):
        if positions[i] - positions[i-1]< step:
            # if overlap (distance between two positions < step), return the index
            return i
    #if no overlap , return maxmum length + 1
    return len(positions)+1

def move_position_from_center_int(positions, center_i, step):
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
    
def move_position_from_center_float(positions, center_i, step):
    # left side
    #print("center_i",center_i)
    #print("before",positions)
    
    # if not pass left bound
    if positions[center_i-1] - step/2 > 0:
        positions[center_i-1] = positions[center_i-1] - step/2
    
    # if not the second element
    if positions[center_i-2] - step/2 > 0:
        for i in range(center_i-2, 0, -1):
            if abs(positions[i] - positions[i+1]) < step:
                 positions[i] = positions[i] - step/2
            else:
                break
    
    # right side            
    positions[center_i] = positions[center_i] + step/2
    for i in range(center_i+1, len(positions)):
        if abs(positions[i] - positions[i-1]) < step:
             positions[i] = positions[i] + step/2
        else:
            break
    #print("after",positions)
    return positions

