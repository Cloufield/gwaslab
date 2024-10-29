def _merge_and_sync_dic(*dics):
    temp = {}
    for dic in dics:
        temp = temp| dic
    return temp