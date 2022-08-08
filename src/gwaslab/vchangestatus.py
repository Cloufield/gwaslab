import pandas as pd

def vchange_status(status,digit,before,after):
    dic={str(i):str(i) for i in range(10)}
    for i in range(len(before)):
        dic[before[i]]=after[i]
    #pattern= (digit-1) * r'\w' + before[i] + (7 - digit)* r'\w'     
    #to_change = status.str.match(pattern, case=False, flags=0, na=False)  
    #if sum(to_change)>0:
    if digit>1:
        status_pre = status.str[:digit-1]
    else:
        status_pre = ""
    status_end=status.str[digit:]
    status_new = status_pre+status.str[digit-1].map(dic)+status_end
    return status_new