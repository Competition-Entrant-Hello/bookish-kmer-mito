import k_genes
import kmer_distr
import bpcontent
import re

LINELENGTH=60
GOLDEN_PATH_DIR='/Users/jtf/datasets/hg38/unmasked/'


def findDist(dictA, dictB):
    total = 0
    for key in dictA:
        term = (dictA[key] - dictB[key]) ** 2
        total += term
    return total ** 0.5


def getMito(k):
    source = open(GOLDEN_PATH_DIR + "Homo_sapiens.GRCh38.dna.chromosome.MT.fa", 'rU')
    m = source.read()
    return kmer_distr.kmer_distr(m, k)

def getWhole(chromosome):
    source = open(GOLDEN_PATH_DIR+'Homo_sapiens.GRCh38.dna.chromosome.'+str(chromosome)+'.fa','rU')
    return source.read().replace("\n", "")



def moving_win(chromosome, chrlength, winlength, k):
    #chro = getWhole(chromosome)
    winstart = 0
    winend = winlength
    out = "winS winE dist AT\n"
    mdict = getMito(k)
    seq = 'temp'
    while chrlength - winstart > winlength:
        seq = k_genes.get_sequence(winstart, winend, chromosome)
        N_count = seq.count('N')
        #print(N_count)
        kmers = kmer_distr.kmer_distr(seq, k)
        #print(kmers)
        dist = findDist(mdict, kmers)
        toAdd = str(winstart) + " " + str(winend) + " " + str(dist) + " "

        #only mito hard code because ugh
        # if winstart == 640000 or winstart == 239456000 or winstart == 241904000:
        #print(len(seq) - N_count)
        leng = len(seq)
        if leng - N_count != 0:
            toAdd += str(bpcontent.findATcontent(seq)/((leng - N_count) * 1.0)) + '\n'
        # else:
        # toAdd += "NA\n"

        print(toAdd)
        if N_count < 100:
            out += toAdd
        winstart = winend
        winend = winstart + winlength
    return out


def marKer(k, chrocomp):
    mit = getMito(k)
    LMito = []
    for key in mit:
        if mit[key] > 0.01:
            LMito.append(key)
    print(LMito)
    start = 0
    end = 3000
    LPos = []
    leng = len(getWhole(chrocomp))
    while end < leng:
        count = 0
        seq = k_genes.get_sequence(start, end, chrocomp)
        for mer in LMito:
            if seq.find(mer) != -1:
                count += 1
        print(count)
        print(start)
        if count > 5:
            LPos.append(start)
        start += 1500
        end += 1500
    return LPos


if __name__ == "__main__":
    # chromosome = 1
    # chrlength = len(getWhole(1))
    # chromosome = 'MT'
    # chrlength = len(getWhole('MT'))
    # output = moving_win(chromosome, chrlength, 3000, 9)
    # print(output)
    #print(marKer(7, 1))
    print(kmer_distr.kmer_distr(getWhole('MT'), 13))
