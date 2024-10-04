import pandas as pd

CATEGORIES = {str(j+i) for j in [1300000,1800000,1900000,3800000,9700000,9800000,9900000] for i in range(0,100000)}

def vchange_status(status,digit,before,after):
    dic={}
    for i in range(len(before)):
        dic[before[i]]=after[i]
    if digit>1:
        return pd.Categorical(status.str[:digit-1]+status.str[digit-1].replace(dic)+status.str[digit:],categories=CATEGORIES)
    else:
        return pd.Categorical(status.str[digit-1].replace(dic)+status.str[digit:],categories=CATEGORIES)

def copy_status(from_status,to_status, digit):
    if digit>1:
        return to_status.str[:digit-1]+from_status.str[digit-1]+to_status.str[digit:]
    else:
        return to_status.str[digit-1]+to_status.str[digit:]


def change_status(status,digit,after):
    prefix= status // 10**(7-digit+1)
    #middle= (status // 10**(7-digit) ) % 10
    suffix= status % 10**(7-digit)
    status = prefix*10**(7-digit+1) + after*10**(7-digit) + suffix
    return status

def schange_status(status,digit,after):
    prefix= status.floordiv(10**(7-digit+1))
    suffix= status.mod(10**(7-digit))
    #prefix= pd.eval("status // 10**(7-digit+1)")
    #suffix= pd.eval("status % 10**(7-digit)")
    #status = prefix*10**(7-digit+1) + after*10**(7-digit) + suffix
    status = pd.eval("prefix*10**(7-digit+1) + after*10**(7-digit) + suffix")
    return status



def status_match(status,digit,to_match):   
    #middle = status.floordiv(10**(7-digit)).mod(10)
    middle = pd.eval('(status // 10**(7-digit) ) % 10')
    if len(to_match)==1:
        is_match = middle==to_match[0]
    else:
        is_match = middle.isin(set(to_match))
    return is_match
'''
def vchange_status(status,digit,before,after):
    for i in range(len(before)):
        matched = status_match(status,digit,[int(before[i])])
        if sum(matched)>0:
            status[matched] = schange_status(status[matched],digit,int(after[i]))
    return status
'''