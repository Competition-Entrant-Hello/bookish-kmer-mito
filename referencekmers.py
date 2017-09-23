import itertools

nucleotides = ['A', 'T', 'C', 'G']
nuleo = ['a', 't', 'c', 'g']



def list(k, upper=True):
    if upper == True:
        list = [''.join(i) for i in itertools.product(nucleotides, repeat = k)]
    else:
        list = [''.join(i) for i in itertools.product(nuleo, repeat = k)]
    return list


#print(list(5))
