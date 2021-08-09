#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__author__ = "Alex Stewart"

from datetime import datetime
from sys import platform
from . import init_config


class ProgDefaults:
    # configuration values
    CONFIG = init_config.validate_config_load(init_config.get_config())
    USER_DB_NAME = CONFIG["db_name"]
    OUTDIR_NAME = CONFIG["outdir_name"]

    # common filenames and default values #
    PROG_NAME = "FACET"
    PROG_VERSION = "v0.1.7"

    # parser aliases
    DB_ALIAS = ['database', 'db']
    FREE_ALIAS = ['db_free', 'database_free', 'dbf']
    OUTFILE_ALIAS = ['outfile', 'o']
    MASKER_ALIAS = ['masker', 'm']
    VC_ALIAS = ['vc', 'vcf', 'v']
    ACCEPTABLE_SUBPARSERS = DB_ALIAS + FREE_ALIAS + OUTFILE_ALIAS + MASKER_ALIAS + VC_ALIAS

    # blast evalues
    BLASTN_EVAL = "1e-15"
    BLASTN_SHORT_EVAL = "1e-5"

    # max_target seqs
    BLASTN_MAX_TARGET_SEQS = 20000

    # buffer for determining whether alignments 'overlap'
    ALN_CLEANING_BUFFER = 25

    # buffer for determining if an alignment spans an entire contig
    TIG_SPAN_BUFFER = 10

    # preset outfmt strings #
    BLASTN_DEFAULT_OUTFMT_SIX = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"
    BASE_OUTFMT_STR = "6 sseqid sstart send qseqid qstart qend sstrand pident"
    VERBOSE_OUTFMT_STR = "6 sseqid sstart send qseqid qstart qend sstrand pident btop"
    MASKER_OUTFMT = "6 sseqid sstart send"

    # outfmt type caster #
    BLASTN_TYPES = dict(qseqid=str, qgi=str, qacc=str, qaccver=str, qlen=int, sseqid=str, sallseqid=str, sgi=str,
                        sallgi=str, sacc=str, saccver=str, sallacc=str, slen=int, qstart=int, qend=int,
                        sstart=int, send=int, qseq=str, sseq=str, evalue=float, bitscore=float, score=int,
                        length=int, pident=float, nident=int, mismatch=int, positive=int, gapopen=int, gaps=int,
                        ppos=float, frames=str, qframe=str, sframe=str, btop=str, staxids=str, sscinames=str,
                        scomnames=str, sblastnames=str, sskingdoms=str, stitle=str, salltitles=str, sstrand=str,
                        qcovs=int, qcovhsp=int)

    # FACET_OUTFMT #
    BASE_INDX = dict(sseqid=0, sstart=1, send=2, qseqid=3, qstart=4, qend=5, sstrand=6, pident=7, btop=8)

    # OUTFMT_SIX #
    OUT_SIX_INDX = dict(qseqid=0, sseqid=1, pident=2, length=3, mismatch=4, gapopen=5, qstart=6, qend=7, sstart=8,
                        send=9, evalue=10, bitscore=11)

    # GFF colors #
    BOTH_ENDS_COLOR = "#56B4E9"
    ONLY_BEGINNING_COLOR = "#E69F00"
    ONLY_END_COLOR = "#F0E442"
    NEITHER_END_COLOR = "#999999"

    # TIC file coverage change needed to add to the list#
    TIC_DELTA = 5

    # Number of alignments needed in list before it flushes out #
    FLUSH_LEN = 100000

    # String of current time to make temp folders that don't recur #
    TIME_STR = datetime.now().strftime('%Y%m%d%H%M%S')

    # Temporary directory to store files in
    PROG_TEMP_DIR = "%s_temp_%s" % (PROG_NAME, TIME_STR)

    # column values in vcf files
    VCF_FMT_DICT = dict(chrom=0, pos=1, id=2, ref=3, alt=4, qual=5, filter=6, info=7, format=8)

    # vcf tempfile names
    VCF_HEADER_TEMPFILE = "%s/vcf_header" % PROG_TEMP_DIR
    VCF_VARIANT_TEMPFILE = "%s/variants.vcf" % PROG_TEMP_DIR
    VCF_TEMPDIR = "%s/VCF_TEMPDIR" % PROG_TEMP_DIR
    VCF_TRUE_SNP_FILE = "UniqueRegion_SNPs.vcf"
    VCF_FALSE_SNP_FILE = "RepetitiveRegion_SNPs.vcf"
    VCF_SUMM_FILE = "summary.txt"

    # RIP MUTATIONS
    RIP_MUT = dict(G="A", C="T", T="C", A="G")

    # rplot parameters
    RPLOT_WINDOW_STEP = 2
    RPLOT_WINDOW_SIZE = 20
    RPLOT_SEARCH_TERM = "AD"
    RPLOT_REF_AGREE_CTOFF = 0

    RPLOT_FILTERED_VCF = "%s/filtered_variamts.vcf" % PROG_TEMP_DIR
    RPLOT_HEADERS = ["sseqid", "start_pos_of_window", "RIP:non", "VCF_coverage", "BLAST_coverage"]
    RPLOT_TEMP_FILE = "%s/RPLOT_TEMP_FILE" % PROG_TEMP_DIR
    RPLOT_FINAL_OUT = "rplot_%swin_%sstep" % (RPLOT_WINDOW_SIZE, RPLOT_WINDOW_STEP)

    # Multiprocessing
    if platform.startswith('linux'):
        from os import sched_getaffinity
        AVAIL_CORES = len(sched_getaffinity(0))
        if AVAIL_CORES - 2 > 0:
            AVAIL_CORES -= 2
    else:  # allows multiprocessing on macOS and others
        from os import cpu_count
        AVAIL_CORES = cpu_count()
        if AVAIL_CORES - 2 > 0:
            AVAIL_CORES -= 2

    # Chunking large files
    BIG_CHUNK_LEN = 250000
    SMALL_CHUNK_LEN = 50000
