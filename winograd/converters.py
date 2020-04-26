"""
Class to convert strings to numbers and vice versa.
"""


class Convert:
    def __init__(self, alphabet=None, worm_length=None):
        self.A = alphabet
        self.wl = worm_length

        if self.A is None:
            self.A = [chr(i) for i in range(256)]
        if self.wl is None:
            self.wl = 5

        self._largest_int = len(self.A) ** self.wl - 1
        self.len_binaries = "{:b}".format(self._largest_int)

    def S2N(self, input_string):
        """String to list of numbers of worm length w"""
        worms = self._subdivide_list(
            self._zero_padding(self._string2list(input_string))
        )
        return self._worms2ints(worms)

    def _string2list(self, input_string):
        """From string to corresponding letters in the alphabet"""
        return [self.A.index(letter) for letter in input_string]

    def _zero_padding(self, v):
        """Makes the input list a list of length multiple of wl, padding with zeros """
        r = len(v) % self.wl
        return v + [0] * ((self.wl - r) % self.wl)

    def _subdivide_list(self, input_list):
        """From a list of number to a list of lists of len wl. """
        return [
            input_list[j * self.wl : (j + 1) * self.wl]
            for j in range(int(len(input_list) / self.wl))
        ]

    def _worms2ints(self, list_of_worms):
        """From a list of worms to the list of corresponding integers"""
        return [self._worm2int(w) for w in list_of_worms]

    def _worm2int(self, worm):
        return sum([w * len(self.A) ** (len(worm) - i - 1) for i, w in enumerate(worm)])

    def N2S(self, v):
        """ From a list of integers to the corresponding word"""
        list_nums = self._trim(self._lifter(self._ints2worms(v)))
        return self._list2string(list_nums)

    def _ints2worms(self, list_of_ints):
        """From list of integers to list of worms"""
        return [self._int2worm(n) for n in list_of_ints]

    def _int2worm(self, num):
        """From integer to corresponding worm"""
        worm = [0] * self.wl
        worm[-1] = num % len(self.A)
        for i in range(self.wl - 2, -1, -1):
            num = (num - worm[i + 1]) / len(self.A)
            worm[i] = int(num % len(self.A))
        return worm

    def _lifter(self, input_list):
        """From list of lists to single list with all their elements"""
        return (
            self._lifter(input_list[0])
            + (self._lifter(input_list[1:]) if len(input_list) > 1 else [])
            if type(input_list) is list
            else [input_list]
        )

    def _trim(self, input_list):
        """Eliminates the zeros at the end of the list """
        while input_list[-1] == 0:
            input_list = input_list[:-1]
        return input_list

    def _list2string(self, list_of_numbers):
        """From a list of numbers to corresponding letters """
        list_of_letters = [self.A[number] for number in list_of_numbers]
        result = ""
        for s in list_of_letters:
            result += s
        return result

    def listTen2ListBin(self, list_of_integers_base_ten):
        return [
            "{0:b}".format(k).zfill(self.len_binaries)
            for k in list_of_integers_base_ten
        ]

    def listBin2ListDec(self, list_of_base_two):
        return [int(k, 2) for k in list_of_base_two]
