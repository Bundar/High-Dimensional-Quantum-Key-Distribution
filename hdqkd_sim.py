# HD QKD Schemes and Randomness Analysis
# Dunbar Birnie IV
# Rutgers ECE

import numpy as np
import random
import math
import matplotlib.pyplot as plt

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
            # drop frame..
            framedropped = 1
            # print("Frame "+str(i)+" Dropped: " + print_time_window(frame))
        else:
            # use frame
            if num_occupied == 1:
                raw_key += getBitString(occupied_bins[0])
            else:
                raw_key += getBitString(get_unoccupied_bin(occupied_bins))
    return raw_key # TODO what if all frames are dropped
        
# prints out the time window data
def print_time_window(time_window):
    s = ''
    for i in time_window:
        s+=str(i)
    return s

# calculates the max # of bits per unit that can be extracted from mthe timestamps 
def calc_h(p):
    return -p*math.log(p,2) - (1-p)*math.log(1-p,2)

# function to calculate photon utilization for simple binning
def calc_photon_utilization(p, n, k):
    FRAME_COUNT = 1
    TIME_UNITS_PER_FRAME = n # n 
    TIME_UNITS_PER_BIN = k # k

    BINS_PER_FRAME = TIME_UNITS_PER_FRAME / TIME_UNITS_PER_BIN
    TOTAL_BIN_COUNT = TIME_UNITS_PER_FRAME*FRAME_COUNT/TIME_UNITS_PER_BIN # n/k

    key_length = 32
    generated_bits = 0
    raw_key = ''
    count = 0

    total_time_units = 0

    while generated_bits < key_length:
        count+=1
        time_window = get_photon_data(p)
        total_time_units += len(time_window)
        raw_key_partial = simple_binning(time_window)
        raw_key += raw_key_partial
        generated_bits += len(raw_key_partial)

    raw_key_rate =  key_length / total_time_units
    optimal_entropy = calc_h(p)
    photon_utilization = raw_key_rate / optimal_entropy

    return photon_utilization

#######################################################################

# function to calculate photon utilization for simple binning. n units per frame, k units per bin
def calc_photon_utilization_single_frame(p, n, k):
    key_length = 1024
    generated_bits = 0
    raw_key = ''
    total_time_units = 0

    while generated_bits < key_length:
        # generate SPDC photon arrivals in frame
        frame = get_photon_data_single_frame(p, n)

        # increment time passed to calculate raw key rate later
        total_time_units += len(frame)
        
        # calculate raw key bits
        raw_key_partial = simple_binning_single_frame(frame, n, k)
        
        # construct total key from partial raw keys
        raw_key += raw_key_partial
        
        # increment stopping condition
        generated_bits += len(raw_key_partial)

    raw_key_rate =  generated_bits / total_time_units
    optimal_entropy = calc_h(p)
    photon_utilization = raw_key_rate / optimal_entropy

    return photon_utilization

# calculates the max # of bits per unit that can be extracted from mthe timestamps 
def calc_h(p):
    return -p*math.log(p,2) - (1-p)*math.log(1-p,2)

# gets the partial key from the time window
def simple_binning_single_frame(frame, n, k):
    key = ''

    # find the bins in the frame which are occupied
    occupied_bins = get_occupied_bins_in_frame(frame, k)
        
    # figure out if we keep this frame
    num_occupied = len(occupied_bins)
    if num_occupied != 1 and num_occupied != BINS_PER_FRAME-1:
        return key # drop frame
    
    unique_bin_index = -1

    # use frame
    if num_occupied == 1:
        unique_bin_index = occupied_bins[0]
    else:
        unique_bin_index = get_unoccupied_bin(occupied_bins, n, k)
    
    return getBitString(unique_bin_index, n, k)

# returns the indexes of the occupied bins in this frame
def get_occupied_bins_in_frame(frame, k):
    bin_indexes = []
    for i,u in enumerate(frame):
        if u == 1:
            # unit occupied
            bin_indexes.append(int(i/k))
    
    return list(set(bin_indexes)) # removes duplicates. important as it allows multiple photons to arrive in the same bin

# returns index of bin that is left out
def get_unoccupied_bin(occupied_bins, n, k):
    for i in range(0, int(n/k)):
        if i not in occupied_bins:
            return i

# Bins are labeled by log(n/k) bit strings
def getBitString(binNum, n, k):
    return ('{0:0'+str(int(math.log(int(n/k),2)))+'b}').format(binNum)

# generate time window... poll SPDC
def get_photon_data_single_frame(p, n):
    return [p_spdc(p) for i in range(0,n)]

# Generate specific graph for case n = 8, k = 1 ... 
def graph_a():
    p_range = np.arange(0, 0.94, 0.01)
    plt.title("N = 8, K = 1")
    plt.xlabel("Probability(p)")
    plt.ylabel("Photon Utilization")
    photon_utilization_per_p = [0]

    for p in p_range:
        print(p)
        if p == 0:
            continue
        else:
            pu = calc_photon_utilization_single_frame(p, 8, 1)
            photon_utilization_per_p.append(pu)
    
    plt.plot(p_range, photon_utilization_per_p, label = "k = 1")
    plt.show()

####








#trying out the new method for spdc simulation
def simple_binning_experiment():
    print("Testing 4 Frames 2 units per bin 8 units per frame, w spdc per unit")
    p = 0.15
    key_length = 128 # bits
    generated_bits = 0
    raw_key = ''
    count = 0

    total_time_units = 0 # for analysis of raw key rate.. r key rate = key_length / total_time_units

    while generated_bits < key_length:
        count+=1
        time_window = get_photon_data(p)
        total_time_units += len(time_window)
        raw_key_partial = simple_binning(time_window)

        print("Iter: " + str(count) + " Time Window: \n" + print_time_window(time_window) + "\nRaw Partial Key: " + raw_key_partial)
        raw_key += raw_key_partial
        generated_bits += len(raw_key_partial)
    
    print("\nGenerated Key: " + raw_key)
    print("\nratio of ones:zeores in raw key: " + str(freq_anal(raw_key)))
    print("Shannon Entropy of raw key: " + str(entropy(raw_key)))
    
    print("\nn = "+ str(TIME_UNITS_PER_FRAME) + ", k = " + str(TIME_UNITS_PER_BIN))
    print("Probability p: " + str(p))
    raw_key_rate =  key_length / total_time_units
    print("Raw Key Rate: " + str(raw_key_rate))
    optimal_entropy = calc_h(p)
    print("h(p) = " +str(optimal_entropy))
    print("Photon Utilization: " +  str(raw_key_rate / optimal_entropy))






# generate graphs fig 2
def fig_2_graphs():
    # genreate graphs from the paper figure 2
    p_step = np.arange(0.01,1,0.01)#[x*0.1 for x in range(1,10, 1)]
    for n in [8, 16, 64]:
        # data x,y 
        plt.title("N = " + str(n))
        plt.xlabel("Probability(p)")
        plt.ylabel("Photon Utilization")
        
        data_k1 = []
        data_k2 = []
        data_k4 = []
        #y = [x*0.05 for x in range(0,12)]

        for k in [1,2,4]:
            photon_u_per_k = []
            for p in p_step:
                photon_u_per_k.append(calc_photon_utilization(p, n, k))
            if k == 1:
                data_k1 = photon_u_per_k
            if k == 2:
                data_k2 = photon_u_per_k
            else:
                data_k4 = photon_u_per_k
        
        # add data sets to graph
        plt.plot(p_step, data_k1, label = "k = 1")
        plt.plot(p_step, data_k2, label = "k = 2")
        plt.plot(p_step, data_k4, label = "k = 4")
        # plot
        plt.show()


def main():
    #simple_binning_experiment()
    #fig_2_graphs()
    graph_a()
    
if __name__ == "__main__":
    main()
