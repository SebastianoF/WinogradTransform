from numpy import testing

from winograd.converters import Convert

ALPHABET = [chr(i) for i in range(256)]
WORM = 5


# --- Inverse consistency tests --- #


def test_string2ints_ints2string():
    c = Convert(ALPHABET, WORM)
    input_string = "Sting & Sons,"
    ints = c.string2ints("Sting & Sons,")
    output_string = c.ints2string(ints)
    testing.assert_equal(input_string, output_string)


def test_ints2string_string2ints():
    c = Convert(ALPHABET, WORM)
    input_ints = [1234, 4321, 1234, 4321, 2356]
    string = c.ints2string(input_ints)
    output_ints = c.string2ints(string)
    testing.assert_equal(input_ints, output_ints)


def test_string2intsb_intsb2string():
    c = Convert(ALPHABET, WORM)
    input_string = "This string will become a worms of binaries"
    worms = c.string2intsb(input_string)
    output_string = c.intsb2string(worms)
    testing.assert_equal(input_string, output_string)


def test_intsb2string_string2intsb():
    c = Convert(ALPHABET, WORM)
    input_binary_worms = [
        "0101010001101000011010010111001100100000",
        "0111001101110100011100100110100101101110",
        "0110011100100000011101110110100101101100",
        "0110110000100000011000100110010101100011",
        "0110111101101101011001010010000001100001",
        "0010000001110111011011110111001001101101",
        "0111001100100000011011110110011000100000",
        "0110001001101001011011100110000101110010",
        "0110100101100101011100110000000000000000",
    ]
    string = c.intsb2string(input_binary_worms)
    output_intsb = c.string2intsb(string)
    testing.assert_equal(input_binary_worms, output_intsb)


def test_string2list_list2string():
    c = Convert(ALPHABET, WORM)
    input_string = "Input string"
    output_list = c._string2list(input_string)
    output_string = c._list2string(output_list)
    testing.assert_equal(input_string, output_string)


def test_list2string_string2list():
    c = Convert(ALPHABET, WORM)
    input_list = [73, 110, 112, 117]
    output_string = c._list2string(input_list)
    output_list = c._string2list(output_string)
    testing.assert_equal(input_list, output_list)


def test_worm2int_int2worm():
    c = Convert(ALPHABET, WORM)
    input_worm = [32, 12, 3, 0, 0]
    output_int = c._worm2int(input_worm)
    output_worm = c._int2worm(output_int)
    testing.assert_equal(input_worm, output_worm)


def test_int2worm_worm2int():
    c = Convert(ALPHABET, WORM)
    input_int = 3234
    output_worm = c._int2worm(input_int)
    output_int = c._worm2int(output_worm)
    testing.assert_equal(input_int, output_int)


# --- Numerical values tests --- #


def test_strings_to_worms_sentence():
    c = Convert(ALPHABET, WORM)
    worms = c.string2ints("My funny valentine!")
    expected_worms = [504224577869, 507350969966, 500068346977, 560295529]
    testing.assert_equal(worms, expected_worms)


def test_single_letter_to_worm():
    c = Convert(ALPHABET, WORM)
    worms = c.string2ints("A")
    expected_worms = [65]
    testing.assert_equal(worms, expected_worms)


def test_single_letter_string2list():
    c = Convert(ALPHABET, WORM)
    output_list = c._string2list("A")
    expected_list = [65]
    testing.assert_equal(output_list, expected_list)


def test_single_letter_string2list_and_padding():
    c = Convert(ALPHABET, WORM)
    output_list = c._zero_padding(c._string2list("A"))
    expected_list = [65, 0, 0, 0, 0]
    print(output_list)
    testing.assert_equal(output_list, expected_list)


def test_single_letter_string2list_and_padding_and_subdivide():
    c = Convert(ALPHABET, WORM)
    output_list = c._subdivide_list(c._zero_padding(c._string2list("A")))
    expected_list = [[65, 0, 0, 0, 0]]
    testing.assert_equal(output_list, expected_list)


def test_convert_single_worm_to_num():
    c = Convert(ALPHABET, WORM)
    intput_worms = [[65, 0, 0, 0, 0]]
    output_ints = c._worms2ints(intput_worms)
    testing.assert_equal(output_ints, [65])

    intput_worms = [[0, 65, 0, 0, 0]]
    output_ints = c._worms2ints(intput_worms)
    testing.assert_equal(output_ints, [65 * len(c.A)])

    intput_worms = [[0, 0, 65, 0, 0]]
    output_ints = c._worms2ints(intput_worms)
    testing.assert_equal(output_ints, [65 * len(c.A) ** 2])


# --- extreme cases tests --- #


def test_strings_to_worm_empty_string():
    c = Convert(ALPHABET, WORM)
    worms = c.string2ints("")
    testing.assert_equal(worms, [])
