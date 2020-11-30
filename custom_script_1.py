#!/usr/bin/env python3

from typing import List
import re
from collections import Counter
import gzip

import argparse as arp

import multiprocessing as mp

import os.path as osp




MAPPER = dict(A = "T",
              T = "A",
              G = "C",
              C = "G",
              N = "N",
              )

def reverse_complement(seq : str) -> str:
    new_seq = list()
    for s in seq[::-1]:
        new_seq.append(MAPPER.get(s,s))
    new_seq = ''.join(new_seq)
    return new_seq

def split_head(header : List[List[str]],
               i1_len : int = 8,
               i2_len : int = 28,
               )->List[List[str]]:

    idx_pattern = "(?<=\:)[CGTAN]{{{}}}\+[CGTAN]{{{}}}"\
        .format(i1_len,i2_len)
    break_pattern = "\:(?=[CGTAN]{{{}}}\+[CGTAN]{{{}}})"\
            .format(i1_len,i2_len)

    break_point = re.search(break_pattern,header).span()[1]
    stem = header[0:break_point]

    i1_tag = re.search(idx_pattern,header).group()
    n_plus = Counter(i1_tag)["+"]
    assert n_plus == 1, \
        "more + signs than expected"
    i1,tag = i1_tag.split("+")

    return (stem,i1,tag)


def run(PTH : str,
         OPTH : str,
         start_umi : int,
         i1_len : int,
         i2_len : int,
         )->None:

    reads = [[],[]]
    barcodes = []
    toggle = -1
    keep = True

    with gzip.open(PTH,"r+") as f:
        for line in f:
            dline = line.decode("utf-8").rstrip("\n")
            if dline[0] == "@":
                toggle += 1
                rx = toggle % 2
                stem,i1,tag = split_head(dline,i1_len,i2_len)
                rtag = reverse_complement(tag)
                quid = stem + i1
                squid = quid.split(" ")
                p1 = squid[0]
                p3 = "_".join(squid[1::])
                bc = rtag[0:start_umi]
                umi = rtag[start_umi::]
                quid = p1 +  "_BARCODE_" + umi  + " " + p3
                reads[rx].append([quid])
                if rx == 0:
                    barcodes.append(bc)
            else:
                reads[rx][-1].append(dline)

    with open(OPTH,"w+") as f:
        for r1,r2,bc in zip(*reads,barcodes):
            chunk = [r1[0],
                    bc + "_" + r1[1] + "_" + r2[1],
                    "+",
                    len(bc) * "A" + "_" + r1[3] + "_" + r2[3],
                    "\n"
            ]
            f.write('\n'.join(chunk))



if __name__ == "__main__":

    prs = arp.ArgumentParser()

    aa = prs.add_argument

    aa("-i","--input",
       nargs = "+",
       required = True,
       )

    aa("-od","--out_dir",
       required = True,
       )


    aa("-su","--start_umi",
       type = int,
       default = 16,
       )


    aa("-i1","--i1_length",
       type = int,
       default = 8,
       )

    aa("-i2","--i2_length",
       type = int,
       default = 28,
       )

    aa("-t","--tag",
       default = "mod_1"
       )


    args = prs.parse_args()


    out_pths = [osp.join(args.out_dir,
                         args.tag + "_" + osp.basename(x).rstrip(".gz")) for\
                x in args.input]

    nprocs = mp.cpu_count()
    print(f"Number of CPU cores: {nprocs}")

    pool = mp.Pool(processes=nprocs)

    pool.starmap(run,[(pth,
                        opth,
                        args.start_umi,
                        args.i1_length,
                        args.i2_length) for pth,opth\
                       in zip(args.input,out_pths)])
