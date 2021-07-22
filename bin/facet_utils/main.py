#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__author__ = "Alex Stewart"

import os
from pathlib import Path
from shutil import rmtree

from . import blast_utils
from . import masker_utils
from .defaults import ProgDefaults as dv
from . import file_export
from . import variant_utils


class ProgParser:
    def __init__(self):
        pass

    def database_run(self, args):
        if args.blaste is None or args.buffer is None or args.writefasta is None or \
                args.num_threads is None:
            print("FATAL: Option specified with no value")
            exit()

        user_db_name = dv.USER_DB_NAME
        outdir_name = dv.OUTDIR_NAME

        # set up names for subject and query #
        contigs = args.subject
        query = args.query
        contigsFilename = str(Path(Path(contigs).resolve().name).with_suffix(''))
        queryFilename = str(Path(Path(query).resolve().name).with_suffix(''))
        facet_prog_dir = str(Path(__file__).parent.absolute().resolve())

        # creates a BLAST database #
        if args.verbose:
            print("\nCreating a BLAST database....")
        dbfilepath = blast_utils.make_blast_db(user_db_name, contigs, facet_prog_dir, args)

        # creates a folder for data exports #
        datafilepath = outdir_name + "/" + contigsFilename

        try:
            os.makedirs(datafilepath)
        except FileExistsError:
            if args.verbose:
                print("\nDirectory \'" + datafilepath + '\' already exists; using it')

        # sets up a filepath for file export #
        filepath = contigsFilename + "_" + queryFilename + "_" + args.task

        if args.notigcov:
            filepath += "_ntc"
        if args.noclean:
            filepath += "_noclean"
        elif args.nocat:
            filepath += "_nocat"
        else:
            filepath += "_clean"

        # checks if the program is going to overwrite any files. If it will, exits the program #
        if not file_export.check_file_overwrite(datafilepath, filepath, args) and not args.force:
            print("FATAL: cowardly refusing to overwrite user files, use \'--force\' to override")
            exit()

        # directory containing temporary output files #
        aln_dir = blast_utils.blast_driver(dbfilepath, contigs, query, args)

        # data export #
        file_export.data_export_driver_new(contigs, aln_dir, datafilepath, filepath, query, args)

    def db_free_run(self, args):
        if args.blaste is None or args.buffer is None or args.writefasta is None:
            print("FATAL: Option specified with no value")
            exit()

        # initialize configuration for database and output directory names
        outdir_name = dv.OUTDIR_NAME

        # set up names for subject and query #
        contigs = args.subject
        query = args.query
        contigsFilename = str(Path(Path(contigs).resolve().name).with_suffix(''))
        queryFilename = str(Path(Path(query).resolve().name).with_suffix(''))

        if args.verbose:
            print("\nRunning BLASTn w/o a database....")
        dbfilepath = contigs

        # creates a folder for data exports #
        datafilepath = outdir_name + "/" + contigsFilename

        try:
            os.makedirs(datafilepath)
        except FileExistsError:
            if args.verbose:
                print("\nDirectory \'" + datafilepath + '\' already exists; using it')

        # sets up a filepath for file export #
        filepath = contigsFilename + "_" + queryFilename + "_" + args.task
        if args.notigcov:
            filepath += "_ntc"
        if args.noclean:
            filepath += "_noclean"
        elif args.nocat:
            filepath += "_nocat"
        else:
            filepath += "_clean"

        # checks if the program is going to overwrite any files. If it will, exits the program #
        if not file_export.check_file_overwrite(datafilepath, filepath, args) and not args.force:
            print("FATAL: cowardly refusing to overwrite user files, use \'--force\' to override")
            exit()

        # directory containing temporary output files #
        aln_dir = blast_utils.blast_driver(dbfilepath, contigs, query, args)

        # data export #
        file_export.data_export_driver_new(contigs, aln_dir, datafilepath, filepath, query, args)

    def outfile_run(self, args):
        if args.buffer is None or args.blaste is None or args.writefasta is None:
            # print(args.buffer, args.eval_filter, args.writefasta)
            print("FATAL: Option specified with no value")
            exit()

        # initialize configuration for database and output directory names
        outdir_name = dv.OUTDIR_NAME

        outfile = Path(args.blast_outfile).resolve()
        outfile_name = str(Path(Path(args.blast_outfile).resolve().name).with_suffix(''))

        contigs = args.subject

        query = args.query

        # creates a folder for data exports #
        datafilepath = outdir_name + "/" + outfile_name
        try:
            os.makedirs(datafilepath)
        except FileExistsError:
            if args.verbose:
                print("\nDirectory \'" + datafilepath + '\' already exists; using it')

        # creates base filename for data exports #
        filepath = outfile_name
        if args.noclean:
            filepath += "_noclean"
        elif args.nocat:
            filepath += "_nocat"
        else:
            filepath += "_clean"

        # checks if the program is going to overwrite any files. If it will, exits the program #
        if not file_export.check_file_overwrite(datafilepath, filepath, args) and not args.force:
            print("FATAL: cowardly refusing to overwrite user files, use \'--force\' to override")
            exit()

        outfile = blast_utils.infile_driver(outfile, args)

        # data export #
        file_export.data_export_driver_new(contigs, outfile, datafilepath, filepath, query, args)

    def masker_run(self, args):
        if args.cov_depth is None or args.mask_char is None or args.blaste is None or args.num_threads is None:
            print("FATAL: Option specified with no value")
            exit()
        if Path("%s_%smasked.fasta" % (Path(args.genome).stem, dv.PROG_NAME.lower())).exists() and not args.force:
            print("FATAL: cowardly refusing to overwrite user files, use \'--force\' to override")
            exit()
        elif Path("%s_%smasked.fasta" % (Path(args.genome).stem, dv.PROG_NAME.lower())).exists() and args.force:
            os.remove(Path("%s_%smasked.fasta" % (Path(args.genome).stem, dv.PROG_NAME.lower())))

        user_db_name = dv.USER_DB_NAME

        facet_prog_dir = str(Path(__file__).parent.absolute().resolve())

        # creates a BLAST database #
        if args.verbose:
            print("\nCreating a BLAST database....")
        dbfilepath = blast_utils.make_blast_db(user_db_name, Path(args.genome).resolve(), facet_prog_dir, args)

        masker_utils.masker_driver(blast_utils.blast_driver(dbfilepath, args.genome, args.genome, args), args)

    def vc_run(self, args):
        if args.genome == "" and args.outfile == "":
            print("FATAL: please provide a genome file or a self-blast outfile using \'--genome <genome.fasta>\' or "
                  "\'--outfile <file.out>\'")
            exit()

        # Checks file overwrites
        if args.verbose:
            print("Checking to see if %s will overwrite files...." % dv.PROG_NAME)
        variant_utils.check_vc_overwrite(args)

        # initialize configuration for database and output directory names
        user_db_name = dv.USER_DB_NAME
        facet_prog_dir = str(Path(__file__).parent.absolute().resolve())

        try:
            os.makedirs(dv.OUTDIR_NAME)
        except FileExistsError:
            if args.verbose:
                print("\nDirectory \'" + dv.OUTDIR_NAME + '\' already exists; using it')

        # the user has provided an outfile and has also provided a genome used to generate the outfile
        if args.outfile != "":
            blast_utils.infile_driver(Path(args.outfile).resolve(), args)

        # the user has not provided an outfile, but has provided a genome that can be used to perform a self-blast
        else:
            # creates a BLAST database #
            if args.verbose:
                print("\nCreating a BLAST database....")
            dbfilepath = blast_utils.make_blast_db(user_db_name, args.genome, facet_prog_dir, args)

            # runs blast on created database
            blast_utils.blast_driver(dbfilepath, args.genome, args.genome, args)

        variant_utils.vc_driver(args)

    def interpret_parser(self, args):
        try:
            os.mkdir(dv.PROG_TEMP_DIR)
        except FileExistsError:
            print("FATAL: temp directory already exists! How did this even happen!?")
            exit()
        except Exception:
            print("FATAL: Something went wrong with creating the temp directory.")
            print("Make sure you have permission to create directories in the current working directory!")
            exit()

        if args.verbose:
            print('.:| ' + dv.PROG_NAME + ' ' + dv.PROG_VERSION + ' |:.''')

        # creates outfmt dictionary and sets args.outfmt to the dict
        blast_utils.outfmt_parser(args)

        if args.subparser_id in dv.DB_ALIAS:
            if args.verbose:
                print("Running %s by creating a BLAST database....\n" % dv.PROG_NAME)
            self.database_run(args)
        elif args.subparser_id in dv.FREE_ALIAS:
            if args.verbose:
                print("Running %s without creating a BLAST database....\n" % dv.PROG_NAME)
            self.db_free_run(args)
        elif args.subparser_id in dv.OUTFILE_ALIAS:
            if args.verbose:
                print("Running %s on a BLAST output file....\n" % dv.PROG_NAME)
            self.outfile_run(args)
        elif args.subparser_id in dv.MASKER_ALIAS:
            if args.verbose:
                print("Using %s to mask repetitive regions in the provided genome....\n" % dv.PROG_NAME)
            self.masker_run(args)
        elif args.subparser_id in dv.VC_ALIAS:
            if args.verbose:
                print("Using %s to find SNP calls in repetitive regions of the provided genome....\n" % dv.PROG_NAME)
            self.vc_run(args)
        try:
            if args.verbose:
                print("\nRemoving temp directory....")
                rmtree(dv.PROG_TEMP_DIR)
                print("Successfully removed temp directory!")
        except Exception:
            print("Temp directory \'%s\' was not successfully removed." % dv.PROG_TEMP_DIR)