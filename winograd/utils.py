"""
Methods to convert strings to numbers and vice versa.
"""

alphabet1 = [" ","A","B","C","D","E",
"F","G","H","I","J","K","L",
"M","N","O","P","Q","R","S",
"T","U","V","W", "X","Y","Z",
":", ".", "'"]

alphabet_ascii = [chr(i) for i in range(256)]

ALPHABET = alphabet_ascii



def S2N(s, wl):
    """String to list of numbers of worm length w"""
    worms = _subdivide_list(_zero_padding(_string2list(s), wl), wl)
    return _worms2ints(worms)
    

def _string2list(string):
    """From string to corresponding letters in the alphabet"""
    return [ALPHABET.index(letter) for letter in string]


def _zero_padding(v, wl):
    """Makes the input list a list of length multiple of wl, padding with zeros """
    r = len(v) % wl    
    return v + [0] * ((wl - r) % wl)


def _subdivide_list(lis, wl):
    """From a list of number to a list of lists of len wl. """
    return [lis[j*wl:(j+1)*wl] for j in range(int(len(lis)/wl))]


def _worms2ints(list_of_worms):
    """From a list of worms to the list of corresponding integers"""
    return [_worm2int(w) for w in list_of_worms]


def _worm2int(worm):
    return sum([w*len(ALPHABET)**(len(worm)- i - 1) for i, w in enumerate(worm)])


def N2S(v, lw):
    """ From a list of integers to the corresponding word"""
    list_nums = _trim(_lifter(_ints2worms(v, lw)))
    return _list2string(list_nums)
    

def _ints2worms(list_of_ints, wl):
    """From list of integers to list of worms"""
    return [_int2worm(n, wl) for n in list_of_ints]


def _int2worm(num, wl):
    """From integer to corresponding worm"""
    worm = [0] * wl
    worm[-1] = num % len(ALPHABET)
    for i in range(wl-2, -1, -1):
        num = (num - worm[i + 1]) / len(ALPHABET)
        worm[i] = int(num % len(ALPHABET))
    return worm
        

def _lifter(li):
    """From list of lists to single list with all their elements"""
    return (
        _lifter(li[0]) + (_lifter(li[1:]) if len(li) > 1 else []) if type(li) is list else [li]
    )


def _trim(v):
    """Eliminates the zeros at the end of the list """
    while v[-1] == 0:
        v = v[:-1]
    return v


def _list2string(list_of_numbers):
    """From a list of numbers to corresponding letters """
    list_of_letters = [ALPHABET[number] for number in list_of_numbers]
    result = ""
    for s in list_of_letters:
        result += s
    return result


def list2binary(list_of_integers):
    max_len = len(ALPHABET) + 1  # TODO not so simple!
    return  ["{0:b}".format(k).zfill(max_len) for k in list_of_integers] 


def list2decimals(list_of_binaries):
    return [int(k, 2) for k in list_of_binaries]
