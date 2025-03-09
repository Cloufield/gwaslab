import copy

def _list_func_args(func):
    return func.__code__.co_varnames

def _extract_kwargs(prefix:str, default:dict, kwargs:dict) -> dict:
    # prefix: keyword
    # default: default dict
    # kwargs: all local kwargs  + args + kwargs

    extracted = []
    extracted_single=dict()
    for key,value in kwargs.items():
        # kwargs or args
        if key=="kwargs" or key=="args":
            for key_nested,value_nested in kwargs[key].items():
                if prefix in key_nested and "arg" in key_nested:

                    if len(key_nested.split("_"))<3:
                        extracted.append(value_nested)
                    ##
                    ## prefix_arg_fontsize
                    else:
                        print(key_nested.split("_")[-1], value)
                        extracted_single[key_nested.split("_")[-1]] = value_nested
        else:
            # local kwargs
            if prefix in key and "arg" in key:
                extracted.append(value)
    if len(extracted_single.keys()) >0:
        extracted.append(extracted_single)
    merged_arg = _merge_and_sync_dic(extracted, default)
    return merged_arg

def _merge_and_sync_dic(list_of_dics:list, default:dict) -> dict:
    temp = copy.copy(default)
    for dic in list_of_dics:
        if isinstance(dic, dict):
            temp.update(dic)
    return temp

def _update_args(args=None, default_args=None):
    
    if default_args is None:
        default_args={}

    if args is None:
        # if None, return default dict
        return default_args
    else:
        # if not None, update default dict
        for key,value in args.items():
            default_args[key] = value
        return default_args

    

def _update_arg(arg=None, default_arg=None):
    if arg is None:
        # if None, return default
        return default_arg
    else:
        # if not None, return arg  
        return arg