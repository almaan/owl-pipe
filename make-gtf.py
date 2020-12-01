#!/usr/bin/env python3

import pandas as pd
import csv
import os.path as osp

import argparse as arp


def read_from_star_index(reference_fn : str,
                         )->pd.DataFrame:

    reference = pd.read_csv(reference_fn,
                            sep = "\t",
                            header = None,
                            index_col = 0,
                            )

    reference.columns = ["size"]

    return reference


def make_gtf(reference : pd.DataFrame,
             bin_size : int = 200000,
             feature : str = "chromosome",
             source : str = "ENSEMBL",
             include_end : bool = True,
             )->pd.DataFrame:

    gtf = dict(name = [],
            source = [],
            feature = [],
            start = [],
            end = [],
            score = [],
            strand = [],
            frame = [],
            attribute =[],
            )

    for chrom in reference.index.values:

        chrom_size = reference.loc[chrom,"size"]
        if chrom_size >= bin_size:
            start = list(range(1,reference.loc[chrom,"size"],bin_size))
            end = [x -1 for x in start[1::]]
            start = start[0:-1]
            if include_end and end[-1] < chrom_size:
                start.append(end[-1])
                end.append(chrom_size)
        elif include_end:
            start = [1]
            end = [chrom_size]

        else:
            continue

        n_bins = len(start)

        gtf["name"] += [chrom] * n_bins
        gtf["start"] += start
        gtf["end"] += end
        gtf["feature"] += [feature] * n_bins

        gtf["source"] += [source] * n_bins
        gtf["score"] += ["."] * n_bins
        gtf["strand"] += ["+"] * n_bins
        gtf["frame"] += ["."] * n_bins
        gtf["attribute"] += ['gene_id "{}_{}";'.format(chrom,n+1) for n in range(n_bins)]

    gtf = pd.DataFrame(gtf)

    return gtf


if __name__ == "__main__":

    prs = arp.ArgumentParser()

    aa = prs.add_argument


    aa("-i",
       "--input",
       required = True,
       )

    aa("-o","--output",
       required = True,
       )

    aa("-bs","--bin_size",
       default = 200000,
       type = int,
       required = False,
       )

    aa("-ie","--include_end",
       default = False,
       action = "store_true",
       )

    aa("-f","--feature",
       default = "chromosome",
       type = str,
       )

    aa("-s","--source",
       default = "ENSEMBL",
       required = False,
       type = str,
       )

    args = prs.parse_args()

    if osp.isdir(args.output):
        base = ".".join(osp.basename(args.input).split("."))[0:-1]
        new_fn = "bs_{}_{}.gtf".format(args.bin_size,base)
        out_pth = osp.join(args.output,new_fn)

    else:
        out_pth = args.output



    reference = read_from_star_index(args.input)

    gtf = make_gtf(reference,
                   bin_size = args.bin_size,
                   feature = args.feature,
                   source = args.source,
                   include_end = args.include_end,
                   )

    gtf.to_csv(out_pth,
               sep = "\t",
               header = None,
               index = None,
               quoting=csv.QUOTE_NONE,
               )

