#!/usr/bin/env python3

import subprocess as sp
from subprocess import Popen, PIPE

import multiprocessing as mp

import shutil
import sys

import custom_script_1 as cs1
import custom_script_2 as cs2

import argparse as arp

import yaml

import os.path as osp
import os

# demuxbyname
# reformat 1
# modify 1
# taggd
# modify 2
# STAR 
# nFeatureCounts
# umi_tools



pth = "mod_2_fj.fq.gz"

def demuxbyname(infile : str,
                outdir : str,
                read : str,
                index_names_file : str,
                index_length : int = 8,
                )->None:

    sp.call([DEMUXBYNAME_PROGRAM,
             "in={}".format(infile),
             "out={}".format(osp.join(outdir,"%","R{}.fq.gz".format(read))),
             "substringmode",
             "length={}".format(index_length),
             "names={}".format(index_names_file),
             ])

def reformat(read_1 : str,
             read_2 : str,
             out_pth : str,
             )->None:

    sp.call([REFORMAT_PROGRAM,
             "in={}".format(read_1),
             "in2={}".format(read_2),
             "out={}".format(out_pth)])


def taggd(reads_infile : str ,
          barcode_infile,
          prefix : str = "demux",
          barcode_start : int = 0,
          )->None:

    if isinstance(barcode_start,int):
        barcode_start = str(barcode_start)
    sp.call(["taggd_demultiplex.py",
             barcode_infile,
             reads_infile,
             prefix,
             "--start-position",barcode_start,
             ])


def star_align(read_1 : str,
               read_2 : str,
               out_dir : str,
               genome_dir : str,
               nthreads : int = None,
               filter_multimap_max : int = 1,
               )->None:
    # last settings from
    # https://github.com/alexdobin/STAR/issues/169
    # should revise this

    if nthreads is None:
        nthreads = str(mp.cpu_count())
    elif isinstance(nthreads,int):
        nthreads = str(nthreads)

    if isinstance(filter_multimap_max,int):
        filter_multimap_max = str(filter_multimap_max)



    sp.call(["STAR",
             "--runThreadN",
             nthreads,
             "--genomeDir",
             genome_dir,
             "--readFilesIn",
             read_1,
             read_2,
             "--readFilesCommand",
             "zcat",
             "--outFilterMultimapNmax",
             filter_multimap_max,
             "--outFileNamePrefix",
             out_dir,
             "--outSAMtype",
             "BAM",
             "SortedByCoordinate",
             "--outFilterScoreMinOverLread",
             "0",
             "--outFilterMatchNminOverLread",
             "0",
             "--outFilterMatchNmin",
             "0",
        ])

def annotate_and_count(infile : str,
                       annotation_file : str,
                       )->None:

    sp.call([FEATURECOUNTS_PROGRAM,
             "-a",
             annotation_file,
             "-o",
             osp.join(osp.dirname(infile),"gene_assigned"),
             infile,
             "-R",
             "BAM",
             "-T",
             "4",
             "-p",
             ])

    sp.call(["samtools",
             "sort",
             infile + ".featureCounts.bam",
             "-o",
             osp.join(osp.dirname(infile),"assigned_sorted.bam"),
             ])

    sp.call(["samtools",
             "index",
             osp.join(osp.dirname(infile),"assigned_sorted.bam"),
    ])


    sp.call(["umi_tools",
             "count",
             "--per-gene",
             "--gene-tag=XT",
             "--assigned-status-tag=XS",
             "--per-cell",
             "-I",
             osp.join(osp.dirname(infile),"assigned_sorted.bam"),
             "-S",
             osp.join(osp.dirname(infile),"counts.tsv.gz"),
             "--wide-format-cell-counts",
             ])


def parse_args():

    prs = arp.ArgumentParser()
    aa = prs.add_argument

    aa("-r1","--read_1",
       required = True,
       )
    aa("-r2","--read_2",
       required = True,
       )
    aa("-o","--out_dir",
       required = True,
       )
    aa("-in","--index_names",
       required = True,
       )

    aa("-stu","--start_umi",
       default = 16,
       help = "start_umi"
       )

    aa("-i1l","--i1_length",
       default = 8,
       help = "index 1 length"
       )

    aa("-i2l","--i2_length",
       default = 28,
       help = "index 2 length"
       )

    aa("-bc","--barcode_file",
       required = True,
       help = "barcode file"
       )

    aa("-tgp","--taggd_prefix",
       default = "demux",
       help = "prefix on taggd files"
       )

    aa("-bcl","--barcode_length",
       default = 16,
       help = "length of barcode",
       )

    aa("-ow","--overwrite",
       default = False,
       action = "store_true",
       )

    aa("-fmm","--filter_multimap_max",
       default = 1,
       )

    aa("-gd","--genome_directory",
       required = True,
       )

    aa("-an","--annotation_file",
       required = True,
       )



    args = prs.parse_args()

    return args


def run_sample():

    args = parse_args()

    with open(args.index_names,"r+") as f:
        index_names = f.readlines()
    index_names = [x.rstrip("\n") for x in index_names]


    if osp.exists(args.out_dir):
        if args.overwrite:
            shutil.rmtree(args.out_dir)
            os.mkdir(args.out_dir)
        else:
            print("[ERROR] : OUTPUT directory exists. Abort.")
            sys.exit(-1)
    else:
        os.mkdir(args.out_dir)

    sub_dirs = {"DEMUX":True,
                "TMP":True,
                "COUNT":True}

    for sd,make_idx_dirs in sub_dirs.items():
        os.mkdir(osp.join(args.out_dir,sd))
        for idx in index_names:
            if make_idx_dirs:
                os.mkdir(osp.join(args.out_dir,sd,idx))

    # demultiplex read 1 w.r.t. sample
    demuxbyname(args.read_1,
                osp.join(args.out_dir,"DEMUX"),
                "1",
                args.index_names,
                args.i1_length,
                )
    # demultiplex read 1 w.r.t. sample
    demuxbyname(args.read_2,
                osp.join(args.out_dir,"DEMUX"),
                "2",
                args.index_names,
                args.i1_length,
                )

    nprocs = mp.cpu_count()
    pool = mp.Pool(processes=nprocs)

    # join R1 and R2 with reformat
    reformat_in_reads = [(osp.join(args.out_dir,"DEMUX",idx,"R1.fq.gz"),
                          osp.join(args.out_dir,"DEMUX",idx,"R2.fq.gz")) for\
                         idx in index_names ]


    reformat_out_reads = [osp.join(args.out_dir,"TMP",idx,"R1R2.fq.gz") for\
                          idx in index_names ]

    reformat_args = [[rs[0],rs[1],od] for\
                     rs,od in zip(reformat_in_reads,
                                  reformat_out_reads)]


    pool.starmap(reformat,reformat_args)


    # move around barcode and prep for taggd

    cs1_outfiles = [osp.join(args.out_dir,"TMP",idx,"mod1_R1R2.fastq") for\
                    idx in index_names]

    cs1_args = [[infile,
                 outfile,
                 args.start_umi,
                 args.i1_length,
                 args.i2_length,
                 ] for infile,outfile in\
                zip(reformat_out_reads,cs1_outfiles)]

    pool.starmap(cs1.run,cs1_args)

    # execute taggd 

    taggd_infiles = cs1_outfiles

    taggd_args = [[infile,
                   args.barcode_file,
                   osp.join(osp.dirname(infile),args.taggd_prefix),
                   "0",
                  ] for infile in taggd_infiles]

    pool.starmap(taggd,taggd_args)


    # split in two files for STAR

    cs2_infiles = [osp.join(osp.dirname(infile),
                            args.taggd_prefix + "_matched.fastq") for\
    infile in taggd_infiles]

    cs2_outfiles = [osp.join(args.out_dir,"TMP",idx,"mod2_R@X.fastq.gz") \
                    for idx in index_names]

    cs2_args = [(infile,outfile,args.barcode_length) for\
                infile,outfile in zip(cs2_infiles,cs2_outfiles)]

    pool.starmap(cs2.run,cs2_args)

    # run star

    star_r1 = [infile.replace("R@X","R1") for infile in cs2_outfiles]
    star_r2 = [infile.replace("R@X","R2") for infile in cs2_outfiles]

    star_outdirs = [osp.join(osp.dirname(infile),"") for infile in cs2_outfiles]

    star_args = [(r1,
                  r2,
                  od,
                  args.genome_directory,
                  nprocs,
                  args.filter_multimap_max,
                  ) for r1,r2,od in\
                 zip(star_r1,
                     star_r2,
                     star_outdirs,
                     )
                 ]

    pool.starmap(star_align,
                 star_args)


    # annotate and count with FeatureCounts

    ann_infiles = [osp.join(args.out_dir,
                            "TMP",
                            idx,
                            "Aligned.sortedByCoord.out.bam") for\
                   idx in index_names
                   ]

    ann_args = [(infile,
                 args.annotation_file,
                 ) for infile in ann_infiles]

    pool.starmap(annotate_and_count,
                 ann_args)


if __name__ == "__main__":
    global DEMUXBYNAME_PROGRAM
    global REFORMAT_PROGRAM
    global FEATURECOUNTS_PROGRAM

    with open("configs/config.yaml","r+") as f:
        config = yaml.load(f,Loader=yaml.FullLoader)

    DEMUXBYNAME_PROGRAM = config["DEMUXBYNAME"]
    REFORMAT_PROGRAM = config["REFORMAT"]
    FEATURECOUNTS_PROGRAM = config["FEATURECOUNTS"]

    run_sample()
