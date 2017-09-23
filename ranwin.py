import random
import kmer_distr
import referencekmers
import splitter
import bpcontent

atwin = []
winstart = []
random.seed(42)


def add_dicts(x, y):
    for key in x:
        #print(key)
        x[key] = x[key] + y[key]
    return x

def pick_random_windows_s(string, size):
    windows = []
    for n in range(30):
        pos = random.randint(0, len(string) - size)
        windows.append(string[pos:(pos+size)])#add the sequence starting from pos
    return windows


def pick_random_kmers(regions, size, k, num = 30):
    windows = []

    for n in range(num):
        seq = ''
        kmers = {}
        for i in referencekmers.list(k):
            kmers[i] = 0
        kmers = kmer_distr.adjust2strand(kmers, 100)
        length = sum(len(s) for s in regions)
        pos = random.randint(0, length - size)
        winstart.append(pos)
        n = 0
        i = 0
        while n <= pos:
                n += len(regions[i])
                i += 1
        i -= 1
        n -= len(regions[i])
        index = pos - n
        length = 0.0
        at = 0
        cg = 0
        while length < size:
            if index >= len(regions[i]):
                #print ("index: " + str(index))
                #print("i: " + str(i))
                #print(seq)
                seq_kmer = kmer_distr.kmer_distr(seq, k, False)
                add_dicts(kmers, seq_kmer)
                at += bpcontent.findATcontent(seq)
                cg += bpcontent.findCGcontent(seq)
                seq = ''
                index = 0
                i += 1
            #print (i)
            #print (index)
            try:
                seq += regions[i][index]
                index += 1
                length += 1
                # print("\n\n\n\n")
                # print("LENGTH: " + str(length))
                # print("seq: " + seq)
                # print("index: " + str(index))
            except IndexError:#if empty region
                index += 1
        seq_kmer = kmer_distr.kmer_distr(seq, k, False)
        add_dicts(kmers, seq_kmer)
        at += bpcontent.findATcontent(seq)
        cg += bpcontent.findCGcontent(seq)
        at = at/length
        cg = cg/length
        atwin.append(at)
        # print (at)
        # print (cg)
        klength = sum(kmers.values()) * 1.0
        for kmer in kmers:
            expAT = bpcontent.findATcontent(kmer)
            expCG = bpcontent.findCGcontent(kmer)
            kmers[kmer] = (kmers[kmer] / klength * 100) / ((at ** expAT) * (cg ** expCG))
        windows.append(kmers)
    return windows


# test = ["ATGCGCTAGCTAGCTAGCT", "AGCTATATCGCTATATATCGC", "GGTACTAGCAGCGC", "GGCGCGATGCATG", "ATATATATGCGCGATCGACGGG", "", "ATGGGCTA"]
#
# print(pick_random_kmers(test, 10, 5, 3))

source = open ('/Users/jtf/datasets/hg38/Homo_sapiens.GRCh38.dna_rm.chromosome.1.fa', 'rU')
#source = open('/Users/jtf/datasets/hgs37/chrM.fa', 'rU')
chr = source.read()
source.close()
chr2 = kmer_distr.nontransposer(chr)

k_list = pick_random_kmers(chr2, 20000, 5)
# for d in k_list:
#     print (sum(d.values()))




source = open('/Users/jtf/datasets/hg38/Homo_sapiens.GRCh38.dna_rm.chromosome.MT.fa', 'rU')
chr = source.read()
#chr = chr.lower() #REMOVE WHEN APPROPRIATE
source.close()

jump = 8500
splat = splitter.split(chr, jump)
for split in splat:
    numN = 0
    for c in split:
        if c == 'N':
            numN += 1
    winstart.append(0)
    winstart.append(len(split))
    #test = kmer_distr.temp_kmer_distr(split, k)
    # atM = bpcontent.findATcontent(split)/(len(split) * 1.0)
    # cgM = bpcontent.findCGcontent(split)/(len(split) * 1.0)
    test = kmer_distr.kmer_distr(split,5, True)
    atM = bpcontent.findATcontent(split)/((len(split) - numN) * 1.0)
    cgM = bpcontent.findCGcontent(split)/((len(split) - numN) * 1.0)
    atwin.append(atM)
    print("\n\n\n\n\n\n\n\n\n\n\n")
    print(bpcontent.findATcontent(split))
    print(bpcontent.findCGcontent(split))
    print(len(split) - numN)
    print(atM + cgM)
    #print("SKABAMWF AWFAWFW        FWAFNWANWANF\n\n\n\n\n\n")
    #print(sum(test.values()))
    for kmer in test:
        expATM = bpcontent.findATcontent(kmer)
        expCGM = bpcontent.findCGcontent(kmer)
        test[kmer] = test[kmer] / ((atM ** expATM) * (cgM ** expCGM))
    k_list.append(test)

del winstart[-1]
del winstart[-1]
#print(k_list)



print("kmer w1 w2 w3 w4 w5 w6 w7 w8 w9 w10 w11 w12 w13 w14 w15 w16 w17 w18 w19 w20 w21 w22 w23 w24 w25 w26 w27 w28 w29 w30 m1 m2")
for kmer in k_list[0]:
    out = kmer + " "
    for dict in k_list:
        out += str(dict[kmer]) + " "
    print(out[:-1])#removes last space from line

last = "atcontent "
for window in atwin:
    last += str(window) + " "
print(last[:-1])

rLast = "winstart "
for pos in winstart:
    rLast += str(pos) + " "
print(rLast[:-1])
