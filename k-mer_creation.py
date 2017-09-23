#k-mer counts of mitochondrial DNA and non-transposer regions of chromosome 1

k_dict = {}

def list_kmers(string, k):#writes down every k-mer of length k in string
    kmers = []
    length = len(string)
    for i in range(0, length - k + 1):
        kmers.append(string[i:i+k])
    return kmers



    



if __name__ == "__main__":
    k = 5
    
    source = open('/Users/jtf/Documents/k-mer_proj/chr1_non-transposer.txt', 'rU')
    chr1_nonT = source.read()
    source.close()
    
    chr1 = chr1_nonT.split('\n')
    for n in chr1:
    	t_list = list_kmers(n, k)
        if t_list != []:
            print (t_list)