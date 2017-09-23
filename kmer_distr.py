import combine_second_strand
import referencekmers

percentages = False





def nontransposer(string): #removes all the lowercase and Ns
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

def seulRepeat(string):
    regions = []
    region = ''
    for char in string:
        if char.islower() and char != 'n':
            region = region + char
        else:
            if len(region) > 0:
                regions.append(region)
                #print(region)
                region = ''
    if len(region) > 0:
        regions.append(region)
    return regions


def list_kmers(string, k):
    kmers = []
    length = len(string)
    for i in range(0, length - k + 1):
        kmers.append(string[i:i+k])
    return kmers



def count_kmers(k_list, k_ref): #makes dictionary with counts of each k-mer
    k_dict = {}
    for mer in k_ref:
        k_dict[mer] = 0
    for kmer in k_list:
        try:
            k_dict[kmer] = k_dict[kmer] + 1
        except KeyError:
            k_dict[kmer] = 1
    return k_dict





def adjust2strand(k_dict, total, percentages = True):
    k_dict_adj = {}
    for kmer in k_dict:
        rcomp = combine_second_strand.secondStrand(kmer)
        if kmer >= rcomp:
            greater = kmer
            lesser = rcomp
        else:
            greater = rcomp
            lesser = kmer
        try:
            if percentages == False:
                k_dict_adj[greater + '/' + lesser] = k_dict[kmer] + k_dict[rcomp]
            else:
                #print(k_dict[kmer] + k_dict[rcomp])
                #print(total)
                k_dict_adj[greater + '/' + lesser] = (k_dict[kmer] + k_dict[rcomp])/total * 100
            k_dict_adj.pop(lesser + '/' + greater, None) #remove duplicate key
        except KeyError: #reverse complement not included
            if percentages == False:
                k_dict_adj[greater + '/' + lesser] = k_dict[kmer]
            else:
                k_dict_adj[greater + '/' + lesser] = (k_dict[kmer])/total * 100
    return k_dict_adj


def kmer_distr(string, k, percentages = True):
     k_list = []
     k_ref = referencekmers.list(k)
     str_list = nontransposer(string)
     #str_list = [string]
     #str_list = seulRepeat(string)
     #print(str_list)
     length = 0.0
     for s in str_list:
          #print(s)
          length += len(s) - k + 1 #number of kmers
          k_list += list_kmers(s, k)
     dict_temp = count_kmers(k_list, k_ref)
     if length == 0:
         dictk = adjust2strand(dict_temp, 100, percentages)
     else:
         dictk = adjust2strand(dict_temp, length, percentages)
     return dictk

def temp_kmer_distr(string, k):
     k_list = []
     k_ref = referencekmers.list(k, False)
     #str_list = nontransposer(string)
     #str_list = [string]
     str_list = seulRepeat(string)
     #print(str_list)
     length = 0.0
     for s in str_list:
          #print(s)
          length += len(s) - k + 1 #number of kmers
          k_list += list_kmers(s, k)
     dict_temp = count_kmers(k_list, k_ref)
     if length == 0:
         dictk = adjust2strand(dict_temp, 100)
     else:
         dictk = adjust2strand(dict_temp, length)
     return dictk



if __name__ == "__main__":

   #  k_ref = referencekmers.list(5)
#
#     source = open('/Users/jtf/Documents/k-mer_proj/k=5/chr1_non-transposon_kmers_k=5.txt', 'rU')
#     chr = source.read()
#     source.close()
#
#     chr = chr.replace('[', '').replace(']', ',').replace('\n','').replace("'", '').replace(' ', '')
#
#     k_listA = chr.split(',')
#
#     #print(len(k_listA))
#
#     k_dictA = count_kmers(k_listA, k_ref)
#
#     k_dict_adjA = adjust2strand(k_dictA, 109921931.0)


#     source = open('/Users/jtf/Documents/k-mer_proj/k=5/chrm_kmers_k=5.txt', 'rU')
#     chrm = source.read()
#     source.close()
#
#     chrm = chrm.replace('[', '').replace(']', ',').replace('\n','').replace("'", '').replace(' ', '')
#
#     k_listB = chrm.split(',')
#
#     #print(len(k_listB))
#
#     k_dictB = count_kmers(k_listB, k_ref)
#
#     k_dict_adjB = adjust2strand(k_dictB, 16179.0)
#
#     print("kmer chr1 chrm")
#     for kmer in k_dict_adjA:
#         #try:
#             print(kmer + " " + str(k_dict_adjA[kmer]) + " " + str(k_dict_adjB[kmer]))
#         #except KeyError:
# #             print(kmer + " " + str(k_dict_adjA[kmer]) + " " + str(k_dict_adjB[kmer[6:] + "/" + kmer[:5]]))

    source = open('/Users/jtf/datasets/hg37/chr1.fa', 'rU')
    chr_1 = source.read()
    source.close()

    print(kmer_distr(chr_1, 5))
