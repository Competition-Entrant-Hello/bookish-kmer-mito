import referencekmers
import k_genes
import combine_second_strand

LINELENGTH=60
GOLDEN_PATH_DIR='/Users/jtf/datasets/hg38/unmasked/'

def create_kmers(string, k, perc = True):
    merdict = {}
    #k_ref = referencekmers.list(k)
    length = len(string)
    # for kmer in k_ref:
    #     merdict[kmer] = 0
    for i in range(0, length - k + 1):
        if 'N' not in string[i:i+k]:
            try:
                merdict[string[i:i+k]] += 1
            except KeyError:
                merdict[string[i:i+k]] = 1
    if perc:
        total = 0.0
        for key in merdict:
            total += merdict[key]
        for key in merdict:
            merdict[key] = merdict[key] / total * 100
    return merdict

def getMito(k):
    source = open(GOLDEN_PATH_DIR + "Homo_sapiens.GRCh38.dna.chromosome.MT.fa", 'rU')
    m = k_genes.get_sequence(1,16569,'MT')
    return create_kmers(m, k, True)

def getRMito(k):
    source = open(GOLDEN_PATH_DIR + "Homo_sapiens.GRCh38.dna.chromosome.MT.fa", 'rU')
    m = combine_second_strand.secondStrand(k_genes.get_sequence(1, 16569, 'MT'))
    #print(m)
    return create_kmers(m, k, True)

def findDist(dictA, dictB):
    total = 0
    for key in dictA:#insert try/except
        try:
            term = (dictA[key] - dictB[key]) ** 2
            total += term
        except KeyError:
            term = (dictA[key] - 0) ** 2
            total += term
    diff = set(dictB.keys()) - set(dictA.keys())
    for key in diff:
        term = (dictB[key] - 0) ** 2
        total += term
    return total ** 0.5


def comp_mit(chromosome, chrlength, winlength, k):
    winstart = 0
    winend = winlength
    out = "winS winE dist strand\n"
    mf = getMito(k)
    mr = getRMito(k)
    while chrlength - winstart > winlength:
        shortstep = False
        seq = k_genes.get_sequence(winstart, winend, chromosome)
        kmers = create_kmers(seq, k)
        #print(kmers)
        # if len(kmers) < 10:
        #     dist = -1
        #else:
        distF = findDist(kmers, mf)
        distR = findDist(kmers,mr)
        if distF <= distR:
            if distF < 3:
                shortstep = True
            toAdd = str(winstart) + " " + str(winend) + " " + str(distF) + " +\n"
        else:
            if distR < 3:
                shortstep = True
            toAdd = str(winstart) + " " + str(winend) + " " + str(distR) + ' -\n'
        #print(toAdd)
        if seq.count('N') < 100:
            out += toAdd
        if shortstep == True:
            winstart += winlength / 10
            winend += winlength / 10
        else:
            winstart += winlength / 2
            winend += winlength / 2
    return out

#def whole_genome():



if __name__ == "__main__":
    #chromo = k_genes.get_sequence(1, 16569, 'MT')
    #print(chromo)
    #print(create_kmers(chromo, 51, False))
    k = 5
    chroms = [1,2,4,5,6,7,8,9,10,11,13,15,16,17,20,'X']
    #chrlengths = [248956422, 242193529, 190214555, ]
    #dists = comp_mit(1, 248956422, 3000, k)
    #dists = comp_mit(2, 242193529, 3000, k)
    #dists = comp_mit(3, 198295559, 3000, k)
    #dists = comp_mit(4, 190214555, 3000, k)
    #dists = comp_mit(5, 181538259, 3000, k)
    #dists = comp_mit(6, 170805979, 3000, k)
    #dists = comp_mit(7, 159345973, 3000, k)
    #dists = comp_mit(8, 145138636, 3000, k)
    dists = comp_mit(9, 138394717, 3000, k)
    #dists = comp_mit(10, 133797422, 3000, k)
    #dists = comp_mit(11, 135086622, 3000, k)
    #dists = comp_mit(12, 133275309, 3000, k)
    #dists = comp_mit(13, 114364328, 3000, k)
    #dists = comp_mit(14, 107043718, 3000, k)
    #dists = comp_mit(15, 101991189, 3000, k)
    #dists = comp_mit(16, 90338345, 3000, k)
    #dists = comp_mit(17, 83257441, 3000, k)
    #dists = comp_mit(18, 80373285, 3000, k)
    #dists = comp_mit(19, 58617616, 3000, k)
    #dists = comp_mit(20, 64444167, 3000, k)
    #dists = comp_mit(21, 46709983, 3000, k)
    #dists = comp_mit('MT', 16569, 3000, k)
    #dists = comp_mit('X', 156040895, 3000, k)
    print(dists)
