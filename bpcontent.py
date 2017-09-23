

def findATcontent(sequence):
    AT = 0
    for c in sequence:
        if c == 'A' or c == 'T' or c == 'a' or c == 't':
            AT += 1
    return AT


def findCGcontent(sequence):
    CG = 0
    for c in sequence:
        if c == 'C' or c == 'G' or c == 'c' or c == 'g':
            CG += 1
    return CG



if __name__ == "__main__":
    print(findATcontent("ATCGTA"))
    print(findCGcontent("ATCGTA"))
