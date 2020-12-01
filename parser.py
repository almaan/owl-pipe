#!/usr/bin/env python3

import argparse as arp

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

    aa("-sf","--start_from",
       default = "top",
       required = False,
       )

    aa("-old","--old_out_dir",
       default = None,
       required = False,
       )

    aa("-h5ad","--h5ad_output",
       default = False,
       action = "store_true",
       required = False,
       )

    aa("-cl","--clean_tmp",
       default = False,
       action = "store_true",
       )

    aa("-vb","--verbose_backend",
       default = False,
       action = "store_true",
       )





    args = prs.parse_args()

    return args

