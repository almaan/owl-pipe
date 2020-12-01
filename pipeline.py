#!/usr/bin/env python3

import subprocess as sp
from subprocess import Popen, PIPE
import tempfile

import multiprocessing as mp

import shutil
import sys

import custom_script_1 as cs1
import custom_script_2 as cs2
from parser import parse_args

import argparse as arp

import yaml

import os.path as osp
import os

import logging


def make_anndata(cnt_pth : str,
                 bc_pth : str,
                 index : str,
                 ):
    import anndata as ad
    import pandas as pd
    from numpy import array

    adata = pd.read_csv(cnt_pth,
                        sep = "\t",
                        header = 0,
                        index_col = 0).T

    if adata.shape[0] * adata.shape[1] < 1:
        return ad.AnnData(array([[]]),
                          obs = pd.DataFrame([[]]),
                          var = pd.DataFrame([]),
                          uns = dict(),
                          )
    else:
        obs = pd.read_csv(bc_pth,
                        sep = "\t",
                        header = None,
                        index_col = 0)

        inter = adata.index.intersection(obs.index)
        print(inter)

        adata = adata.loc[inter,:]
        obs = obs.loc[inter,:]
        obs.columns = pd.Index(["x","y"])
        obs.index = pd.Index(obs.index.astype(str))

        var = pd.DataFrame(adata.columns.values,
                        index = adata.columns.values,
                        columns = ["gene"],
        )

        uns = dict(index = index)

        adata = ad.AnnData(adata,
                        obs = obs,
                        var = var,
                        uns = uns,
        )
        return adata


def proc_call(cmd : list,
              verbose : bool = False,
              )->None:

    if verbose:
        stdout = sb.STDOUT
        stderr = sb.STDERR
    else:
        stdout = tempfile.TemporaryFile()
        stderr = stdout

    sp.call(cmd,
            stdout = stdout,
            stderr = stderr,
            )

    if not verbose:
        stdout.close()


def demuxbyname(infile : str,
                outdir : str,
                read : str,
                index_names_file : str,
                index_length : int = 8,
                verbose : bool = False,
                )->None:


    proc_call([DEMUXBYNAME_PROGRAM,
               "in={}".format(infile),
               "out={}".format(osp.join(outdir,"%","R{}.fq.gz".format(read))),
               "substringmode",
               "length={}".format(index_length),
               "names={}".format(index_names_file),
               ],
              verbose = verbose,
              )

def reformat(read_1 : str,
             read_2 : str,
             out_pth : str,
             verbose : bool = False,
             )->None:

    proc_call([REFORMAT_PROGRAM,
               "in={}".format(read_1),
               "in2={}".format(read_2),
               "out={}".format(out_pth)],
              verbose = verbose,
              )


def taggd(reads_infile : str ,
          barcode_infile,
          prefix : str = "demux",
          barcode_start : int = 0,
          verbose : bool = False,
          )->None:

    if isinstance(barcode_start,int):
        barcode_start = str(barcode_start)
    proc_call(["taggd_demultiplex.py",
             barcode_infile,
             reads_infile,
             prefix,
             "--start-position",barcode_start,
             ],verbose = verbose)


def star_align(read_1 : str,
               read_2 : str,
               out_dir : str,
               genome_dir : str,
               nthreads : int = None,
               filter_multimap_max : int = 1,
               verbose : bool = False,
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


    proc_call(["STAR",
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
             ],
              verbose = verbose,
             )

def annotate_and_count(infile : str,
                       out_dir : str,
                       annotation_file : str,
                       nthreads : int,
                       verbose : bool = False,
                       )->None:

    proc_call([FEATURECOUNTS_PROGRAM,
             "-a",
             annotation_file,
             "-o",
             osp.join(out_dir,"gene_assigned"),
             infile,
             "-R",
             "BAM",
             "--Rpath",
             out_dir,
             "-T",
             str(nthreads),
             "-p",
             "-t",
             "chromosome",
             ],
            verbose = verbose,
            )

    proc_call(["samtools",
             "sort",
             osp.join(out_dir,osp.basename(infile) + ".featureCounts.bam"),
             "-o",
             osp.join(out_dir,"annot_sorted.bam"),
             ],
              verbose = verbose,
              )

    proc_call(["samtools",
             "index",
             osp.join(out_dir,"annot_sorted.bam"),
             ],
              verbose = verbose,
            )

    proc_call(["umi_tools",
             "count",
             "--per-gene",
             "--gene-tag=XT",
             "--assigned-status-tag=XS",
             "--per-cell",
             "-I",
             osp.join(out_dir,"annot_sorted.bam"),
             "-S",
             osp.join(out_dir,"counts.tsv.gz"),
             "--wide-format-cell-counts",
             ],
              verbose = verbose,
            )

def run_sample():

    args = parse_args()

    logging.basicConfig(format = '[%(levelname)s:%(asctime)s] : %(message)s',
                        datefmt='%m:%d:%Y:%I:%M:%S',
                        level=logging.DEBUG,
                        )

    # ------- Settings for pipeline entrance point
    start_pos = {"TOP":0,
                 "MOD1":1,
                 "TAGGD":2,
                 "MOD2":3,
                 "MAP":4,
                 "ANN":5,
                 }

    args.start_from = args.start_from.upper()

    assert args.start_from in start_pos.keys(),\
        "Invalid starting point choose one from :\n {}"\
        .format("\n".join(list(start_pos.keys())))

    if args.start_from != "TOP":
        assert args.old_out_dir is not None,\
            "need to provide a directory to start from"

    if start_pos[args.start_from] != 0:
        logging.info("Entering Pipeline at : {}".format(args.start_from))
        logging.info("Will read data from : {}".format(args.old_out_dir))
    else:
        logging.info("Start pipeline from scratch")

    with open(args.index_names,"r+") as f:
        index_names = f.readlines()
    index_names = [x.rstrip("\n") for x in index_names]

    logging.info("Using indices :\n {}".format("\n ".join(index_names)))

    # ------- For Parallel Processing
    nprocs = mp.cpu_count()
    pool = mp.Pool(processes=nprocs)

    logging.info("Using {} Processes".format(nprocs))

    # ------- Prepare Directories for output
    if osp.exists(args.out_dir):
        if args.overwrite:
            shutil.rmtree(args.out_dir)
            logging.info("Will overwrite already existing output folder {}".\
                         format(args.out_dir))
            os.mkdir(args.out_dir)
        else:
            logging.error("OUTPUT directory {} already exists. Abort."\
                          .format(args.out_dir))
            logging.info('To overwrite old directory,'\
            'include the "--overwrite" flag')
            sys.exit(-1)
    else:
        os.mkdir(args.out_dir)

    sub_dirs = {"DEMUX":True,
                "TMP":True,
                "COUNT":False}

    for sd,make_idx_dirs in sub_dirs.items():
        os.mkdir(osp.join(args.out_dir,sd))
        for idx in index_names:
            if make_idx_dirs:
                os.mkdir(osp.join(args.out_dir,sd,idx))


    # ------- Demuultiplex by sample
    if start_pos[args.start_from] == start_pos["TOP"]:

        logging.info("Initiated Demultiplexing by name")

        # demultiplex read 1 w.r.t. sample
        demuxbyname(args.read_1,
                    osp.join(args.out_dir,"DEMUX"),
                    "1",
                    args.index_names,
                    args.i1_length,
                    verbose = args.verbose_backend,
                    )
        # demultiplex read 1 w.r.t. sample
        demuxbyname(args.read_2,
                    osp.join(args.out_dir,"DEMUX"),
                    "2",
                    args.index_names,
                    args.i1_length,
                    verbose = args.verbose_backend,
                    )

        reformat_in_reads = [(osp.join(args.out_dir,"DEMUX",idx,"R1.fq.gz"),
                              osp.join(args.out_dir,"DEMUX",idx,"R2.fq.gz")) for\
                             idx in index_names ]

        reformat_out_reads = [osp.join(args.out_dir,"TMP",idx,"R1R2.fq.gz") for\
                              idx in index_names ]

        reformat_args = [[rs[0],rs[1],od] for\
                         rs,od in zip(reformat_in_reads,
                                      reformat_out_reads)]


        pool.starmap(reformat,reformat_args)

        logging.info("Completed Demultiplexing by name")

    # ------- Move SB and UMI elements, prepare for Taggd
    if start_pos[args.start_from] <= start_pos["MOD1"]:
        logging.info("Started first modification")
        if start_pos[args.start_from] == start_pos["MOD1"]:
            reformat_out_reads = [osp.join(args.old_out_dir,"TMP",idx,"R1R2.fq.gz") for\
                                  idx in index_names ]


        cs1_outfiles = [osp.join(args.out_dir,"TMP",idx,"mod1_R1R2.fastq") for\
                        idx in index_names]

        cs1_args = [[infile,
                    outfile,
                    args.start_umi,
                    args.i1_length,
                    args.i2_length,
                    args.verbose_backend,
                    ] for infile,outfile in\
                    zip(reformat_out_reads,
                        cs1_outfiles)]

        pool.starmap(cs1.run,cs1_args)
        logging.info("Completed first modification")
        if args.clean_tmp:
            for fn in reformat_out_reads:
                os.remove(fn)

    # ------- Match SB using Taggd
    if start_pos[args.start_from] <= start_pos["TAGGD"]:
        logging.info("Initiated matching of Spatial Barcodes")
        if start_pos[args.start_from] == start_pos["TAGGD"]:
            taggd_infiles = [osp.join(args.old_out_dir,
                                     "TMP",
                                      idx,
                                      "mod1_R1R2.fastq") for\
                             idx in index_names]
        else:
            taggd_infiles = cs1_outfiles

        taggd_args = [[infile,
                       args.barcode_file,
                       osp.join(osp.dirname(infile),args.taggd_prefix),
                       "0",
                       args.verbose_backend,
                       ] for infile in taggd_infiles]

        pool.starmap(taggd,taggd_args)

        logging.info("Completed matching of Spatial Barcodes")
        if args.clean_tmp:
            for fn in taggd_infiles:
                os.remove(fn)

    # ------- Split Read 1 and Read 2 into separate files, remove observed SB
    if start_pos[args.start_from] <= start_pos["MOD2"]:
        logging.info("Initiated second modification")
        if start_pos[args.start_from] == start_pos["MOD2"]:
            taggd_infiles = [osp.join(args.old_out_dir,
                                     "TMP",
                                      idx,
                                      "mod1_R1R2.fastq") for\
                             idx in index_names]

        cs2_infiles = [osp.join(osp.dirname(infile),
                           args.taggd_prefix + "_matched.fastq") for\
                           infile in taggd_infiles]

        cs2_outfiles = [osp.join(args.out_dir,"TMP",idx,"mod2_R@X.fastq.gz") \
                        for idx in index_names]


        cs2_args = [(infile,outfile,args.barcode_length) for\
                    infile,outfile in zip(cs2_infiles,cs2_outfiles)]

        pool.starmap(cs2.run,cs2_args)
        logging.info("Completed second modification")
        if args.clean_tmp:
            for fn in cs2_infiles:
                os.remove(fn)

    # ------- Map to reference using STAR
    if start_pos[args.start_from] <= start_pos["MAP"]:
        logging.info("Initiated Mapping with STAR")
        logging.info("Using Genome Directory : {}"\
        .format(args.genome_directory))

        if start_pos[args.start_from] == start_pos["MAP"]:

            cs2_outfiles = [osp.join(args.old_out_dir,
                                     "TMP",
                                     idx,
                                     "mod2_R@X.fastq.gz") \
                            for idx in index_names]

        star_r1 = [infile.replace("R@X","R1") for infile in cs2_outfiles]

        star_r2 = [infile.replace("R@X","R2") for infile in cs2_outfiles]

        star_outdirs = [osp.join(osp.dirname(infile),"") for infile in cs2_outfiles]

        star_args = [(r1,
                      r2,
                      od,
                      args.genome_directory,
                      nprocs,
                      args.filter_multimap_max,
                      args.verbose_backend,
                      ) for r1,r2,od in\
                     zip(star_r1,
                         star_r2,
                         star_outdirs,
                        )]

        pool.starmap(star_align,
                    star_args)

        logging.info("Completed Mapping with STAR")
        if args.clean_tmp:
            for fn1,fn2 in zip(star_r1,
                               star_r2,
                               ):
                os.remove(fn1)
                os.remove(fn2)



    # ------- Annotate and count using FeatureCounts and umi_tools
    if start_pos[args.start_from] <= start_pos["ANN"]:
        logging.info("Initiated Annotation and Counting with featureCounts and umi_tools")
        if start_pos[args.start_from] == start_pos["ANN"]:
            ann_infiles_dir = args.old_out_dir
        else:
            ann_infiles_dir = args.out_dir

        ann_infiles = [osp.join(ann_infiles_dir,
                                "TMP",
                                idx,
                                "Aligned.sortedByCoord.out.bam") for\
                    idx in index_names
                    ]

        ann_outdirs = [osp.join(args.out_dir,"TMP",
                                index,
                                ) for index in index_names]

        ann_args = [(infile,
                     od,
                    args.annotation_file,
                     nprocs,
                     args.verbose_backend,
                    ) for infile,od in zip(ann_infiles,ann_outdirs)]

        pool.starmap(annotate_and_count,
                     ann_args)
        logging.info("Completed Annotation and Counting with featureCounts and umi_tools")
        if args.clean_tmp:
            for fn in ann_infiles:
                os.remove(fn)

    # ------- Move counts to COUNTS folder
    for index,dir_name in zip(index_names,ann_outdirs):
        source = osp.join(dir_name,"counts.tsv.gz")
        dest_dir = osp.join(args.out_dir,"COUNT")
        if args.h5ad_output:
            adata = make_anndata(source,
                                 args.barcode_file,
                                 index,
                                 )

            dest = osp.join(dest_dir,
                            "{}_counts.h5ad".format(index))
            adata.write_h5ad(dest)
            if osp.exists(dest): os.remove(source)
        else:
            dest = osp.join(dest_dir,
                            "{}_counts.tsv.gz".format(index))
            shutil.move(source,dest)

    if args.clean_tmp:
        shutil.rmtree(osp.join(args.out_dir,
                               "TMP",
                               ))

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
