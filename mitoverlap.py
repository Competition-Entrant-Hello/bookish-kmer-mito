import kmer_distr


# 1 kb overlapping windows of mtdna. columns are kmers of each window

LINELENGTH=60
GOLDEN_PATH_DIR='/Users/jtf/datasets/hg38/unmasked/'


def mv_winO(chromosome, winlength, k):
    source = open(GOLDEN_PATH_DIR+'Homo_sapiens.GRCh38.dna.chromosome.'+str(chromosome)+'.fa','rU')
    chro = source.read()
    winstart = 0
    winend = winlength
    out = "" #make header
    while winend < len(chro):
        window = chro[winstart:winend]
        kmers = kmer_distr.kmer_distr(window, k)
        out += str(winstart) + " " + str(winend) + " "
        header = "win_S win_E "
        total = 0
        for kmer in kmers:
            #print(kmer)
            total += kmers[kmer]
            out += str(kmers[kmer]) + " "
            header += kmer + " "
        #print(total)
        header += "\n"
        out += "\n"
        #print(out)
        winstart += 500
        winend += 500
    return header + out


def mitovars(file, winlength):
    source = open(file)
    chro = open(GOLDEN_PATH_DIR + 'Homo_sapiens.GRCh38.dna.chromosome.MT.fa', 'rU').read()
    line = source.readline()
    assert line.find('#') == 0
    tokens = line.split()
    anno = source.read()
    annoL = anno.split("\n")
    posL = []
    for line in annoL:
        temp = line.split()
        if len(temp) != 0:
            posL.append(temp[2])

    winstart = 0
    winend = winlength

    L = []
    while winend < len(chro):
        varwin = 0
        for pos in posL:
            if int(pos) > winstart and int(pos) < winend:
                varwin += 1
        L.append(varwin)
        winstart += 500
        winend += 500
    return L




if __name__ == "__main__":
    toPrint = mv_winO("MT", 1000, 9)
    print(toPrint)
    # toL = mitovars("/Users/jtf/Documents/k-mer_proj/mitovar.txt", 1000)
    # for pos in toL:
    #     print(pos)
