import kmer_distr
import referencekmers

#!/usr/local/bin/python
LINELENGTH=60
GOLDEN_PATH_DIR='/Users/jtf/datasets/hg38/unmasked/'

def get_sequence(genomic_start,genomic_end,chromosome,repeat_masked=0):
  if repeat_masked==0:
    gp_file=open(GOLDEN_PATH_DIR+'Homo_sapiens.GRCh38.dna.chromosome.'+str(chromosome)+'.fa','rU')
  else:
    gp_file=open(GOLDEN_PATH_DIR+'Homo_sapiens.GRCh38.dna.chromsome.'+str(chromosome)+'.masked.fa','rU')
  line=gp_file.readline()
  assert line.find(">")==0
  tokens=line.split()
  chromosome_id=tokens[0][1:]
  line_file_pos=(genomic_start/LINELENGTH)*(LINELENGTH+1)
  line_pos=genomic_start%LINELENGTH
  if line_pos==0:
    file_pos=line_file_pos+line_pos-2
  else:
    file_pos=line_file_pos+line_pos-1
  search_length=genomic_end-genomic_start+1
  gp_file.seek(file_pos,1)
  seq=''
  while len(seq)<search_length:
    s=gp_file.readline().strip()
    seq=seq+s
  gp_file.close()
  seq=seq[:search_length]
  return seq


def add_dicts(x, y):
    for key in x:
        #print(key)
        x[key] = x[key] + y[key]
    return x



def kmer_specific(chromosome, k, size, coding=True):
    gene_file = open("/Users/jtf/datasets/hg38/genTable.txt")
    line = gene_file.readline()
    assert line.find("#") == 0
    tokens = line.split()
    listseq = []
    kmers = {}
    #print ("BAM")
    #print(referencekmers.list(k))
    for i in referencekmers.list(k):
        kmers[i] = 0
    kmers = kmer_distr.adjust2strand(kmers, 100)
    # print(tokens)
    #txst i = 3; txend i = 4
    while tokens[1] != "chr" + chromosome:#read through file until hit desired chromosome
        line = gene_file.readline()
        tokens = line.split()
        #print(tokens[1])
    start = 0
    while tokens[1] == "chr" + chromosome:
        if coding:
            start = int(tokens[5])
            end = int(tokens[6])
        else:
            end = int(tokens[5])
        #print(tokens[1])
        if chromosome == 'M':
            seq = get_sequence(start, end, 'MT')
            # print(start)
            # print(end)
            listseq.append(seq)
        else:
            seq = get_sequence(start, end, chromosome)
            listseq.append(seq)
        if not coding:
            start = int(tokens[6])
        line = gene_file.readline()
        tokens = line.split()
    # print(listseq)
    if coding == False:
        if chromosome == 'M':
            seq = get_sequence(start, size, 'MT')
            listseq.append(seq)
        else:
            seq = get_sequence(start, size, chromosome)
            listseq.append(seq)
    length = 0.0
    for seq in listseq:
        kdict = kmer_distr.kmer_distr(seq, k)
        length += sum(kdict.values())
        # print("SEQ:\n\n" + seq + "\n\n\n")
        # print("LENACT: " + str(len(seq)))
        # print("seq len: " + str(len(seq) - k + 1))
        # print("seq len: " + str(sum(kdict.values())))
        # # for key in kdict:
        # #     print(key + str(kdict[key]))
        # print("cur len: " + str(length))
        #length += len(seq) - k + 1
        #sums = 0
        # for i in kdict:
        #     sums += kdict[i]
        # print(sums)
        #print('BOM')
        #print(kmers)
        #print('BAM')
        #print(kdict)
        add_dicts(kmers, kdict)
        # print(kmers)
    #print("LENGTH: " + str(length))
    for kmer in kmers:
        kmers[kmer] = kmers[kmer] / length * 100
    return kmers


def mito_k_toText(dictC, dictN):
    out = "kmer coding non-coding\n"
    for kmer in dictC:
        out += (kmer + " " +  str(dictC[kmer]) + " " + str(dictN[kmer]) + "\n")
    return out

if __name__ == "__main__":
    #print(kmer_distr.kmer_distr(get_sequence(1000000, 2000000, 10), 5))
    k_thing = kmer_specific('M', 5, 16569, True)
    sums = 0
    for i in k_thing:
        sums += k_thing[i]
    #print(k_thing)
    print(sums)
    dictC = kmer_specific('M', 5, 16569, True)
    dictN = kmer_specific('M', 5, 16569, False)
    print(mito_k_toText(dictC, dictN))
