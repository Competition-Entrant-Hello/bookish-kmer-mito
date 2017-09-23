import kmer_distr



LINELENGTH=60
GOLDEN_PATH_DIR='/Users/jtf/datasets/hg38/unmasked/'

def get_sequence(genomic_start,genomic_end,chromosome,repeat_masked=0):
  if repeat_masked==0:
    gp_file=open(GOLDEN_PATH_DIR+'Homo_sapiens.GRCh38.dna.chromosome.'+str(chromosome)+'.fa','rU')
  else:
    gp_file=open(GOLDEN_PATH_DIR+'Homo_sapiens.GRCh38.dna.chromsome.'+str(chromosome)+'.masked.fa','rU')
  line=gp_file.readline()
  print('\n\n\n')
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

def getMito(k):
    source = open(GOLDEN_PATH_DIR + "Homo_sapiens.GRCh38.dna.chromosome.MT.fa", 'rU')
    m = source.read()
    return kmer_distr.kmer_distr(m, k)

def findDist(dictA, dictB):
    total = 0
    for key in dictA:
        term = (dictA[key] - dictB[key]) ** 2
        total += term
    return total ** 0.5




if __name__ == "__main__":
    print(get_sequence(9085000, 9098000, 1))
    # gen = get_sequence(106802629, 106806090, 1)
    # print(gen)
    # mit = get_sequence(5855, 9328, 'MT')
    # print(mit)
    # print("\n\n\n")
    # print(findDist(kmer_distr.kmer_distr(gen, 5), kmer_distr.kmer_distr(mit, 5)))
    # genfake = get_sequence(107001106, 107004602, 1)
    # print(findDist(kmer_distr.kmer_distr(genfake, 5), kmer_distr.kmer_distr(mit, 5)))
    # print(findDist(kmer_distr.kmer_distr(gen, 5), (getMito(5))))
    #
    # print("\n\n\n")
    # test = get_sequence(106803000, 106806000, 1)
    # print(findDist(kmer_distr.kmer_distr(test, 5), (getMito(5))))
