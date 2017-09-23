import sys

PATH_DIR = '/Users/jtf/Documents/k-mer_proj/chr'


def getTopDist(chrom, toList = False):
    source = open(PATH_DIR + chrom + '/dist' + chrom + '_k5_sort.txt')
    dists = source.read()
    numlines = dists.count('\n')
    lines = dists.split('\n')
    top = int(0.00015 * numlines) + 1
    out = ""
    L = []
    for i in range(0, top):
        out += lines[i] + '\n'
        L.append(lines[i])
    if toList:
        return L
    return out

def getNumts(chrom):
    source = open('/Users/jtf/Documents/k-mer_proj/numtpos.txt')
    numts = source.read()
    lines = numts.split('\n')
    posL = []
    truePosL = []
    for line in lines:
        temp = line.split()
        posL.append(temp[6])
    for pos in posL:
        sep = pos.find('-')
        colon = pos.find(':')
        theChro = pos[:colon]
        start = pos[colon + 1:sep]
        end = pos[sep + 1:]
        try:
            # print(theChro)
            # print(chrom)
            if (theChro == chrom):
                truePosL.append(theChro + ':'+ str((int(start) + int(end)) / 2))
        except ValueError:
            pass
    return truePosL

def getNew(chrom):
    peaks = getTopDist(chrom, True)
    del peaks[0]
    theNewbies = []
    for peak in peaks:
        peako = peak.split()
        avg = (int(peako[0]) + int(peako[1])) / 2
        new = True
        for known in getNumts('chr' + chrom):
            pos = int(known[known.find(':') + 1:])
            if avg > pos - 10000 and avg < pos + 10000:
                # print('avg:' + str(avg))
                # print('pos:' + str(pos))
                new = False
        if new == True:
            theNewbies.append('chr' + chrom + ':' + peako[0] + '-' + peako[1] + " " + peako[2])
    return theNewbies


def getNewAll():
    chroms = ['01','02','04','05','06','07','08','09','10','11','13','15','16','17','20','X']
    allNewbies = []
    for chro in chroms:
        newbies = getNew(chro)
        allNewbies += newbies
    return allNewbies




if __name__ == "__main__":
    #print(getTopDist(sys.argv[1]))
    #print(getNumts('chr01'))
    #print(getNew('01'))
    for thing in getNewAll():
        print(thing)
    #print(getNewAll())
