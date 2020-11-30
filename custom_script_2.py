#!/usr/bin/env python3

import gzip
import os.path as osp
import argparse as arp
import multiprocessing as mp
import re

from typing import List



def run(PTH : str,
         OPTH : List,
         barcode_length : int,
         )->None:

    chunks_1 = []
    chunks_2 = []
    chunk_1 = []
    chunk_2 = []

    if PTH.endswith("gz"):
        open_fun = lambda x: gzip.open(x,"r+")
        line_fun = lambda x: x.decode("utf-8").rstrip("\n")
    else:
        open_fun = lambda x: open(x,"r+")
        line_fun = lambda x: x.rstrip("\n")

    if not OPTH.endswith("gz"):
        OPTH = OPTH + ".gz"


    with open_fun(PTH) as f:
        for line in f:
            dline = line_fun(line)
            if dline[0] == "@":

                bc = re.search("(?<=B0:Z:)[ACTG]{{{}}}".format(barcode_length),
                               dline).group(0)

                head = dline.replace("BARCODE",bc)
                chunk_1 = [head]
                chunk_2 = [head]
            else:
                if len(chunk_1) in [1,3]:
                    _,p_1,p_2 = dline.split("_")
                    chunk_1.append(p_1)
                    chunk_2.append(p_2)
                    if len(chunk_1) == 4:
                        chunks_1.append("\n".join(chunk_1))
                        chunks_2.append("\n".join(chunk_2))

                else:
                    chunk_1.append(dline)
                    chunk_2.append(dline)

    out_pth_1 = OPTH.replace("R@X","R1")
    out_pth_2 = OPTH.replace("R@X","R2")

    for out_pth,chunk in zip([out_pth_1,out_pth_2],
                            [chunks_1,chunks_2]):

        with gzip.open(out_pth,"wt") as f:
            f.write("\n".join(chunk) + "\n")



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


    aa("-bcl","--barcode_length",
       type = int,
       default = 16,
       )

    aa("-t","--tag",
       default = "mod_2",
       )


    args = prs.parse_args()


    out_pths = [osp.join(args.out_dir,
                         args.tag + "_R@X_" + osp.basename(x)) for\
                x in args.input]

    nprocs = mp.cpu_count()
    print(f"Number of CPU cores: {nprocs}")

    pool = mp.Pool(processes=nprocs)

    pool.starmap(run,[(pth,opth,args.barcode_length) for pth,opth\
                       in zip(args.input,out_pths)])
