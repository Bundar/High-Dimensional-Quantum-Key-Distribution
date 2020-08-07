# HD QKD Schemes and Randomness Analysis
# Dunbar Birnie IV
# Rutgers ECE

import numpy as np
import random
import math
import matplotlib.pyplot as plt

# better spdc prob p of being in bin
def p_spdc(p):
    r = random.random()
    if r < p:
        return 1
    else:
        return 0

# function to calculate photon utilization for simple binning. n units per frame, k units per bin
def calc_photon_utilization_single_frame(p, n, k):
    key_length = 128
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
    if num_occupied != 1 and num_occupied != (n/k)-1:
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

# Generate general graph for case n, k=1,2,4... 
def graph(n):
    p_range = np.arange(0, 0.94, 0.05)
    plt.title("N = "+str(n))
    plt.xlabel("Probability(p)")
    plt.ylabel("Photon Utilization")
    k1 = [0]
    k2 = [0]
    k4 = [0]

    for p in p_range:
        print(p)
        if p == 0:
            continue
        else:
            pu1 = calc_photon_utilization_single_frame(p, n, 1)
            pu2 = calc_photon_utilization_single_frame(p, n, 2)
            pu4 = calc_photon_utilization_single_frame(p, n, 4)

            k1.append(pu1)
            k2.append(pu2)
            k4.append(pu4)
    
    plt.plot(p_range, k1, label = "k = 1")
    plt.plot(p_range, k2, label = "k = 2")
    plt.plot(p_range, k4, label = "k = 4")

    plt.show()

def main():
    for n in [64]:
        graph(n)
    
if __name__ == "__main__":
    main()