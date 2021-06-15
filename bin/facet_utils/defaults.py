#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__author__ = "Alex Stewart"


class ProgDefaults:
    # common filenames and default values #
    PROG_NAME = "FACET"
    PROG_VERSION = "v0.1.3"

    ACCEPTABLE_SUBPARSERS = ['database', 'db', 'db_free', 'database_free', 'dbf', 'outfile', 'masker']

    BLASTN_EVAL = "1e-15"
    BLASTN_SHORT_EVAL = "1e-5"
    ALN_CLEANING_BUFFER = 25

    TIG_SPAN_BUFFER = 10

    # outfmt strings #
    BLASTN_DEFAULT_OUTFMT_SIX = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"
    BASE_OUTFMT_STR = "6 sseqid sstart send qseqid qstart qend sstrand pident"
    DERIP_OUTFMT_STR = "6 sseqid sstart send qseqid qstart qend sstrand pident btop"
    VERBOSE_OUTFMT_STR = "6 sseqid sstart send qseqid qstart qend sstrand pident btop qseq"

    # BLAST indices #
    # BASE_OUTFMT
    SSEQID_INDEX = 0
    SSTART_INDEX = 1
    SEND_INDEX = 2
    QSEQID_INDEX = 3
    QSTART_INDEX = 4
    QEND_INDEX = 5
    SSTRAND_INDEX = 6
    PIDENT_INDEX = 7
    # VERBOSE_OUTFMT
    BTOP_INDEX = 8

    # OUTFMT_SIX
    FMTSIX_QSEQID = 0
    FMTSIX_SSEQID = 1
    FMTSIX_PIDENT = 2
    FMTSIX_LENGTH = 3
    FMTSIX_MISMATCH = 4
    FMTSIX_GAPOPEN = 5
    FMTSIX_QSTART = 6
    FMTSIX_QEND = 7
    FMTSIX_SSTART = 8
    FMTSIX_SEND = 9
    FMTSIX_EVALUE = 10
    FMTSIX_BITSCORE = 11

    # GFF colors #
    BOTH_ENDS_COLOR = "#56B4E9"
    ONLY_BEGINNING_COLOR = "#E69F00"
    ONLY_END_COLOR = "#F0E442"
    NEITHER_END_COLOR = "#999999"

    # String Masking Values #
    #      DEPRECATED       #
    MASKER_CEILING_CHAR = '~'
    MASKER_CEILING_INT = ord(MASKER_CEILING_CHAR)

    MASKER_MID_CHAR = 'O'
    MASKER_MID_INT = ord(MASKER_MID_CHAR)

    MASKER_FLOOR_CHAR = '!'
    MASKER_FLOOR_INT = ord(MASKER_FLOOR_CHAR)

    # TIC file coverage change #
    TIC_DELTA = 5
