# HD QKD Schemes and Randomness Analysis
# Dunbar Birnie IV
# Rutgers ECE

import random
import math

# time is partitioned in frames consisting of n time units. We take n to be a power
# of two.
FRAME_COUNT = 4
TIME_UNITS_PER_FRAME = 8 # n

# Each frame is divided into n/k bins, each consisting of
# k ≤ n time units, and we are free to choose k. Note that k also needs to be a power of two in order
# for n to be divisible by k.
TIME_UNITS_PER_BIN = 2 # k
TOTAL_BIN_COUNT = TIME_UNITS_PER_FRAME*FRAME_COUNT/TIME_UNITS_PER_BIN # n/k

def entropy(string):
        # get probability of chars in string
        prob = [ float(string.count(c)) / len(string) for c in dict.fromkeys(list(string)) ]

        # calculate the entropy
        entropy = - sum([ p * math.log(p) / math.log(2.0) for p in prob ])

        return entropy

# Freq Analysis
def freq_anal(rawKey):
    zeroes = 0
    ones = 1
    for b in rawKey:
        if b == '0':
            zeroes+=1
        else:
            ones+=1
    return ones/zeroes
# better spdc prob p of being in bin
def better_spdc(p):
    r = random.random()
    if r < p:
        return 1
    else:
        return 0

# Returns time unit in which photon has arrived. 
def spdc_per_bin(time_bins_range):
    return random.randint(0,time_bins_range-1)

# Returns time unit in which photon has arrived. 
def spdc_per_unit(time_units_range):
    unit = random.randint(0,time_units_range-1)
    bin_num = int(unit / TIME_UNITS_PER_BIN)
    return bin_num

# Bins are labeled by log(n/k) bit strings
def getBitString(binNum):
    return ('{0:0'+str(int(math.log(TOTAL_BIN_COUNT,2)))+'b}').format(binNum)

# gets the partial kkey from the time window
def simple_binning(time_window):
    for i in range(0, FRAME_COUNT): # 0 1 2 3
        partial_frame_key = ''
        for j in range(0, TIME_UNITS_PER_FRAME/TIME_UNITS_PER_BIN): # 0 1 2 3
            detected = 0
            if time_window[i*FRAME_COUNT+j] or time_window[i*FRAME_COUNT + j + 1]: # only works for bins size 2 TODO
                detected+=1
                partial_key = getBitString(i*FRAME_COUNT+j)
                partial_frame_key += partial_key

    return 'a'


        
# experiment 1 w old spdc
def exp_1():
    print("Testing 4 Frames 2 units per bin 8 units per frame, w spdc per unit")
    key_length = 128 # bits
    generated_bits = 0
    spdc_source = ''
    raw_key = ''
    count = 0
    while generated_bits < key_length:
        count+=1
        a = spdc_per_bin(TOTAL_BIN_COUNT)
        spdc_source+=str(a)
        raw_key_partial = getBitString(a)
        print("Iter: " + str(count) + " a: " + str(a) + " Raw Partial Key: " + raw_key_partial)
        raw_key += raw_key_partial
        generated_bits += len(raw_key_partial)

    print("\nGenerated Key: " + raw_key)
    print("\nratio of ones:zeores in raw key: " + str(freq_anal(raw_key)))
    print("Shannon Entropy of raw key: " + str(entropy(raw_key)))

    print("SPDC output: " + spdc_source)
    print("Shannon Entropy of spdc sim souorce: " + str(entropy(spdc_source)))
    print("")

# experiment 2 w old spdc
def exp_2():
    print("Testing 4 Frames 2 units per bin 8 units per frame, w spdc per bin")
    key_length = 128 # bits
    generated_bits = 0
    spdc_source = ''
    raw_key = ''
    count = 0
    while generated_bits < key_length:
        count+=1
        a = spdc_per_unit(TIME_UNITS_PER_FRAME*FRAME_COUNT)
        spdc_source+=str(a)
        raw_key_partial = getBitString(a)
        print("Iter: " + str(count) + " a: " + str(a) + " Raw Partial Key: " + raw_key_partial)
        raw_key += raw_key_partial
        generated_bits += len(raw_key_partial)

    print("\nGenerated Key: " + raw_key)
    print("\nratio of ones:zeores in raw key: " + str(freq_anal(raw_key)))
    print("Shannon Entropy of raw key: " + str(entropy(raw_key)))

    print("SPDC output: " + spdc_source)
    print("Shannon Entropy of spdc sim souorce: " + str(entropy(spdc_source)))
    print("")


#trying out the new method for spdc simulation
def exp_3_new_method():
    print("Testing 4 Frames 2 units per bin 8 units per frame, w spdc per bin")
    p = 0.15
    key_length = 128 # bits
    generated_bits = 0
    spdc_source = ''
    raw_key = ''
    count = 0
    while generated_bits < key_length:
        count+=1
        time_window = []
        for i in range(0,TIME_UNITS_PER_FRAME*FRAME_COUNT):
            a = better_spdc(p)
            time_window = time_window + a
        
        raw_key_partial = simple_binning(time_window)

        print("Iter: " + str(count) + " a: " + str(a) + " Raw Partial Key: " + raw_key_partial)
        raw_key += raw_key_partial
        generated_bits += len(raw_key_partial)

    print("\nGenerated Key: " + raw_key)
    print("\nratio of ones:zeores in raw key: " + str(freq_anal(raw_key)))
    print("Shannon Entropy of raw key: " + str(entropy(raw_key)))

    print("SPDC output: " + spdc_source)
    print("Shannon Entropy of spdc sim souorce: " + str(entropy(spdc_source)))
    print("")


# just to test the output
def test_bit_string_gen():
    for i in range(4):
        a = spdc_per_bin(TOTAL_BIN_COUNT)
        print(str(a) + " : " + getBitString(a))


def main():
    exp_3_new_method()

    

    
if __name__ == "__main__":
    main()