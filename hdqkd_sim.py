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
BINS_PER_FRAME = TIME_UNITS_PER_FRAME / TIME_UNITS_PER_BIN
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
def p_spdc(p):
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
    return ('{0:0'+str(int(math.log(int(BINS_PER_FRAME),2)))+'b}').format(binNum)

# returns the indexes of the occupied bins in this frame
def get_occupied_bins(frame):
    bin_indexes = []
    for i,u in enumerate(frame):
        if u == 1:
            # unit occupied
            bin_indexes.append(int(i/TIME_UNITS_PER_BIN))
    
    return list(set(bin_indexes)) # removes duplicates. important as it allows multiple photons to arrive in the same bin

# returns index of bin that is left out
def get_unoccupied_bin(occupied_bins):
    for i in range(0, int(BINS_PER_FRAME)):
        if i not in occupied_bins:
            return i

# generate time window... poll SPDC
def get_photon_data(p):
    time_window = []
    for i in range(0, int(TIME_UNITS_PER_FRAME * FRAME_COUNT)):
        spdc_result = p_spdc(p)
        time_window.append(spdc_result)
    return time_window


# gets the partial kkey from the time window
def simple_binning(time_window):
    raw_key = ''
    for i in range(0, FRAME_COUNT):
        # calc frame  by frame
        frame = time_window[i*TIME_UNITS_PER_FRAME:(i+1)*TIME_UNITS_PER_FRAME]
        
        # find the bins in the frame which are occupied
        occupied_bins = get_occupied_bins(frame)
        
        # figure out if we keep this frame
        num_occupied = len(occupied_bins)
        if num_occupied != 1 and num_occupied != BINS_PER_FRAME-1:
            # drop frame.. include analytics later TODO
            print("Frame "+str(i)+" Dropped: " + print_time_window(frame))
        else:
            # use frame
            if num_occupied == 1:
                raw_key += getBitString(occupied_bins[0])
            else:
                raw_key += getBitString(get_unoccupied_bin(occupied_bins))
    return raw_key
        
# prints out the time window data
def print_time_window(time_window):
    s = ''
    for i in time_window:
        s+=str(i)
    return s

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
def simple_binning_experiment():
    print("Testing 4 Frames 2 units per bin 8 units per frame, w spdc per unit")
    p = 0.15
    key_length = 128 # bits
    generated_bits = 0
    spdc_source = ''
    raw_key = ''
    count = 0
    while generated_bits < key_length:
        count+=1
        time_window = get_photon_data(p)
        
        raw_key_partial = simple_binning(time_window)

        print("Iter: " + str(count) + " Time Window: \n" + print_time_window(time_window) + "\nRaw Partial Key: " + raw_key_partial)
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
    simple_binning_experiment()
    

    
if __name__ == "__main__":
    main()
