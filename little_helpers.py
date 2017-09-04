
class no_error_dictionary:
    def __init__(self,dic_in):
        self.dic = dic_in

    def __getitem__(self, key):
        try:
            return self.dic[key]
        except KeyError:
            return None


def flatten_dictionary(dictionary):
    new_dic = {}
    for key,value in dictionary.items():
        for key2,value in value.items():
            new_dic[key+'_'+key2] = value
    return new_dic


