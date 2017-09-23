#k-mer counts of mitochondrial DNA and non-transposon regions of chromosome 1

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


def rm_nontransposer(string): #removes all the lowercase and Ns
    regions = []
    region = ''
    for char in string:
        if char.isupper() and char != 'N':
            region = region + char
            #print(char)
        else:
            if len(region) > 0:
                regions.append(region)
                region = ''
    if len(region) > 0:
        regions.append(region)
    return regions




if __name__ == "__main__":
    k = 5
    source = open ('/Users/jtf/datasets/hgs37/chr1/chr1.fa', 'rU')
    chr1 = source.read()
    source.close()

    #chr1 = chr1.replace('\n', '').replace('>chr1', '')

    #chr1_regions = rm_nontransposer(chr1)

    #print(chr1_regions)
    #for r in chr1_regions:
    #    print(r)

    source = open('/Users/jtf/datasets/hgs37/chrM.fa', 'rU')
    chr_m = source.read()
    source.close()

    chr_m = chr_m.replace('\n', '').replace('>chrM', '')

    chrm_regions = rm_nontransposer(chr_m)

    for r in chrm_regions:
        print(r)
    #k_list_1 = list_kmers(chr1, k)
    #k_list_m = list_kmers(chr_m, k)
    #print count_kmers(k_list_1)
