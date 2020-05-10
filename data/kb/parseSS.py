#!/usr/bin/env python3

import sys

def main():
    Plate_ID = 0
    index2 = 1
    Twist_RNA_copies = 2
    ATCC_RNA_copies = 3
    ATCC_virus_copies = 4
    spike_copies = 5
    lysate = 6
    nCoV_amplicon = 7
    nCoV_primer_nM = 8
    RPP30_primer_nM = 9
    RPP30_inner_primer_nM = 10
    bc_set = 11
    RT_temp = 12
    PCR_cycles = 13
    Sample_Well = 14
    index = 15
    Sample_ID = 16


    plate_idx = 0
    index2_idx = 1
    index_idx = -2
    well_idx = -3
    rand_idx1 = -10
    rand_idx2 = -11    

    for line in sys.stdin:
        if "Plate" not in line:
            continue
        line = line.strip()
        data = line.split(",")
        if data[0] == "Plate_ID":
            continue
        
        sys.stdout.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(data[index] + data[index2], data[index], data[index2], data[Plate_ID], data[Sample_Well], data[nCoV_amplicon], data[lysate], data[Twist_RNA_copies], data[ATCC_RNA_copies], data[ATCC_virus_copies]))

    return 1

if __name__=="__main__":
    main()
