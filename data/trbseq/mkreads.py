#!/usr/bin/env python3

import sys
import random


bases = "ACGT"

def mutate(orig_string, mutation_rate=0.005):
    result = []
    mutations = []
    for base in orig_string:
        if random.random() < mutation_rate:
            new_base = bases[bases.index(base) - random.randint(1, 3)] # negatives are OK
            result.append(new_base)
            mutations.append((base, new_base))
        else:
            result.append(base)
    return "".join(result)

def make_read(seq):
    h = "NB552046:38:HLHM2BGXF:1:11101:23150:1029 1:N:0:0"
    l = len(seq)
    qual = l*"I"
    r = "@{}\n{}\n+\n{}\n".format(h, mutate(seq), qual)
    #r = "@{}\n{}\n+\n{}\n".format(h, seq, qual)
    return r

gene_ref = {"N1": "CACCCCGCATTACGTTTGGTGGACCCTCAGATTCAACTGGCAGTAACCAGAATGGAGAACGC",
           "N2": "CGCAAATTGCACAATTTGCCCCCAGCGCTTCAGCGTTCTTCGGAATGTCGCGCATTGGCATGGAAGTCA",
           "RPP30": "GGTTCTGACCTGAAGGCTCTGCGCGGACTTGTGGAGACAGCCGCTC"}

genes = ['N1', 'N1_spikein', 'RPP30']
constant = "GGTTAATTCTGATGAGCGTACC"

def main():
    bcs = []
    i = 0
    with open("trb_bcs.txt", "r") as f, open("whitelist.txt", "w") as bar:
        for lidx, line in enumerate(f):
            if lidx%100==0:
                i+=1
                bar.write(line)
                line = line.strip()
                bcs.append(line)
                if i==96:
                    break
    with open("R1.fastq", "w") as r1, open("R2.fastq","w") as r2:
        for lidx, l in enumerate(sys.stdin):
            l = l.strip()
            data = l.split(",")
            data = list(map(int, data))
            print(data)
            for gidx, cnt in enumerate(data):
                cnt = cnt
                g = genes[gidx]
                if g == "N1":
                    for i in range(cnt):
                        seq = gene_ref["N1"][:50]
                        b = bcs[lidx] + constant
                        # print(seq)
                        # print(b)
                        r2.write(make_read(seq))
                        r1.write(make_read(b))
                elif g == "N1_spikein":
                    for i in range(cnt):
                        seq = gene_ref["N2"][:50]
                        b = bcs[lidx] + constant
                        r2.write(make_read(seq))
                        r1.write(make_read(b))
                elif g=="RPP30":
                    for i in range(cnt):
                        seq = gene_ref[g][:50]
                        b = bcs[lidx] + constant
                        r2.write(make_read(seq))
                        r1.write(make_read(b))



if __name__=="__main__":
    main()
