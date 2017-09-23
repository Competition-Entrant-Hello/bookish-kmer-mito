import kmer_distr

#euclidean distance between first half of mito reference and second half of chondria reference

GOLDEN_PATH_DIR='/Users/jtf/datasets/hg38/unmasked/'

def findDist(dictA, dictB):
    total = 0
    for key in dictA:
        term = (dictA[key] - dictB[key]) ** 2
        total += term
    return total ** 0.5

def getDist(k):
    source = open(GOLDEN_PATH_DIR+'Homo_sapiens.GRCh38.dna.chromosome.MT.fa','rU')
    mito = source.read()
    first = mito[:len(mito)/2]
    second = mito[len(mito)/2:]
    f = kmer_distr.kmer_distr(first, k)
    s = kmer_distr.kmer_distr(second, k)
    return findDist(f, s)



print(getDist(5))
