
def secondStrand(string):
    other = ''
    for c in string:
        if c == 'A':
            other += 'T'
        elif c == 'a':
            other += 't'
        elif c == 'T':
            other += 'A'
        elif c == 't':
            other += 'a'
        elif c == 'G':
            other += 'C'
        elif c == 'g':
            other += 'c'
        elif c == 'C':
            other += 'G'
        elif c == 'c':
            other += 'g'
    return other[::-1]





if __name__ == "__main__":
    print(secondStrand('GCAAT'))
