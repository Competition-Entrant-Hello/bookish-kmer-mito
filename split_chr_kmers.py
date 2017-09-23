import kmer_distr
import splitter
import referencekmers

#split chr1
#for each ^: rm non_transposer
#list k_mers
#k_mer count
#adjust 2 strand
#format print

k = 5

ref_kmers = referencekmers.list(k)




source = open ('/Users/jtf/datasets/hgs37/chr1/chr1.fa', 'rU')
chr = source.read()
source.close()

total  = 109921931.0

splat = splitter.split(chr, 10000000) #pieces are each 10 Mbp

regions_list = []

for piece in splat:
    regions_list.append(kmer_distr.rm_nontransposer(piece))

#print(regions_list)

# k_list = []
# k_dict_list = []
k_dictadj_list = []


for piece in regions_list:#every 10 Mbp piece
    for region in piece:#every non-transposon region in that 10 Mbp piece
        list = kmer_distr.list_kmers(region, k)
        dict = kmer_distr.count_kmers(list, ref_kmers)
        k_dictadj_list.append(kmer_distr.adjust2strand(dict, total))
        
        

# print(k_list)
# 
# 
# 
# for list in k_list:
#     #print(list)
#     k_dict_list.append(kmer_distr.count_kmers(list, ref_kmers))#add kmer counts for each region to list of dictionaries
# 
# print(k_dict_list)
# 
# 
# 
# 
# for dict in k_dict_list:
#     k_dictadj_list.append(kmer_distr.adjust2strand(dict, total))



header = 'kmer '
n = 1
for dict in k_dictadj_list:
    header += (str(n) + " ")
    n += 1
header = header[:-1]

print(header)
for kmer in k_dictadj_list[0]:
    line = kmer + " "
    for dicti in k_dictadj_list:
        line += str(dicti[kmer])
    print(line)






    