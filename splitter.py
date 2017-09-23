

def split(string, jump):
    pieces = []
    while len(string) >= jump:
        #print("test")
        #print(string[:jump])
        pieces.append(string[:jump])
        #print(pieces)
        string = string[jump:]
    if len(string) > 0:
        pieces.append(string)
    return pieces

#print(split("potato", 2))
