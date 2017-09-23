import kmer_distr
import splitter
import ranwin


k = 5

source = open ('/Users/jtf/datasets/hg37/chr1.fa', 'rU')
#source = open('/Users/jtf/datasets/hgs37/chrM.fa', 'rU')
chr = source.read()
source.close()


chr2 = kmer_distr.nontransposer(chr)


# jump = 10000000 #pieces are each 10 Mbp
# #jump = 8500
# splat = splitter.split(chr, jump)


window_size = 20000
splat = ranwin.pick_random_windows(chr2, window_size)




dictL = []

for split in splat:
    #test = kmer_distr.temp_kmer_distr(split, k)
    test = kmer_distr.kmer_distr(split, k)
    dictL.append(test)
    #print (kmer_distr.kmer_distr(split, k, 10000000.0))
    #print(test)
    #print(sum(test.values()))
    # sum = 0
#     for key in test:
#         sum += test[key]
#     print(sum)
# print(dictL)

source = open('/Users/jtf/datasets/hg37/chrM.fa', 'rU')
chr = source.read()
#chr = chr.lower() #REMOVE WHEN APPROPRIATE
source.close()

jump = 8500
splat = splitter.split(chr, jump)
for split in splat:
    #test = kmer_distr.temp_kmer_distr(split, k)
    test = kmer_distr.kmer_distr(split,k)
    #print("SKABAMWF AWFAWFW        FWAFNWANWANF\n\n\n\n\n\n")
    #print(sum(test.values()))
    dictL.append(test)

# dictL2 = []
#
# for dicto in dictL:
#     dictm = {}
#     for kmer in sorted(dicto.items()):
#         dictm[kmer.upper()] = dicto[kmer]
#     dictL2.append(dictm)
#
# print(dictL2)

header = "kmer 10MbpJUMPS mthalves"
print(header)

for kmer in dictL[0]:
    out = kmer + " "
    for dict in dictL:
        out += str(dict[kmer]) + " "
    print(out[:-1])#removes last space from line
