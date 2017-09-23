import combine_second_strand

k_dict = {}


def list_kmers(string, k):
    kmers = []
    length = len(string)
    for i in range(0, length - k + 1):
        kmers.append(string[i:i+k])
    return kmers




def count_kmers(k_list):
    for kmer in k_list:
        try:
            k_dict[kmer] = k_dict[kmer] + 1
        except KeyError:
            k_dict[kmer] = 1
    return k_dict



def adjust2strand(k_dict):
    k_dict_adj = {}
    for kmer in k_dict:
        try:
            k_dict_adj[kmer] = k_dict[kmer] + k_dict[combine_second_strand.secondStrand(kmer)] #adds reverse complement to every value
        except KeyError:
            k_dict_adj[kmer] = k_dict[kmer]
    return k_dict_adj




if __name__ == "__main__":
    k = 5
	
    source = open('/Users/jtf/datasets/hgs37/chrM.fa', 'rU')
    chr_m = source.read()
    source.close()
    
    chr_m = chr_m.replace('\n', '').replace('>chrM', '')
    
    #print(combine_second_strand.secondStrand('GCAAT'))
    
    k_list_m = list_kmers(chr_m, k)
    
    k_dict = count_kmers(k_list_m)
    
    print(adjust2strand(k_dict))