


repelforce = 0.02 
yspan=100
positions = [1 , 3, 10, 19, 26, 50,51,52, 62,71,99]



def text_reposition( positions, span, repel_force=0.01):
    # if interval > repleforce * span
    # need adjust 

def adjust_position(positions, yspan, repelforce, max_iter=100):
    step = yspan*repelforce
    print(step)
    for i in range(max_iter):
        index = check_overlap(positions, step)
        print(index)
        if index == len(positions)+1:
            return positions
        else:
            move_position_from_center(positions, index, step)
    return positions

def check_overlap(positions,step):
    for i in range(1,len(positions)):
        if positions[i] - positions[i-1]< step:
            return i
    return len(positions)+1

def move_position_from_center(positions, center_i, step):
    # left side
    print("center_i",center_i)
    print("before",positions)
    if center_i-1>0:
        positions[center_i-1] = positions[center_i-1] - step/2
    
    if center_i-2>0:
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
    print("after",positions)
    return positions
    
    