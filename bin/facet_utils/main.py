#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__author__ = "Alex Stewart"

import os
from pathlib import Path

from . import blast_utils
from . import masker_utils
from .defaults import ProgDefaults as dv
from . import file_export
from . import init_config


class ProgParser:
    def __init__(self):
        pass

    def database_run(self, args):
        if args.blaste is None or args.shorte is None or args.buffer is None or args.writefasta is None or \
                args.num_threads is None:
            print("FATAL: Option specified with no value")
            exit()

        # initialize configuration for database and output directory names
        config = init_config.validate_config_load(init_config.get_config())
        user_db_name = config["db_name"]
        outdir_name = config["outdir_name"]

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
        if args.shortonly:
            filepath = contigsFilename + "_" + queryFilename + "_shortonly"
        elif args.runshort:
            filepath = contigsFilename + "_" + queryFilename + "_short"
        else:
            filepath = contigsFilename + "_" + queryFilename
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

        # runs normal blast with the newly created database, cleans each list, and stores cleaned list in nested dict
        if not args.shortonly:
            if args.verbose:
                print("\nRunning BLASTn....")
            blastdata = blast_utils.blastn_driver(dbfilepath, contigs, query, args)
        else:
            blastdata = {}

        # runs short blast with the newly created database and cleans each repeat
        if args.runshort or args.shortonly:
            if args.verbose:
                print("\nRunning BLASTn-short....")
            shortblastdata = blast_utils.short_driver(dbfilepath, contigs, query, blastdata, args)
        else:
            shortblastdata = {}

        # combines shortblast hits with blast hits
        combinedhits = blast_utils.cat_megablast_short(blastdata, shortblastdata, args)

        # data export #
        file_export.data_export_driver(contigs, combinedhits, datafilepath, filepath, query, args)

    def db_free_run(self, args):
        if args.blaste is None or args.shorte is None or args.buffer is None or args.writefasta is None:
            print("FATAL: Option specified with no value")
            exit()

        # initialize configuration for database and output directory names
        config = init_config.validate_config_load(init_config.get_config())
        outdir_name = config["outdir_name"]

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
        if args.shortonly:
            filepath = contigsFilename + "_" + queryFilename + "_shortonly"
        elif args.runshort:
            filepath = contigsFilename + "_" + queryFilename + "_short"
        else:
            filepath = contigsFilename + "_" + queryFilename
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

        # runs normal blast with the newly created database, cleans each list, and stores cleaned list in nested dict
        if not args.shortonly:
            if args.verbose:
                print("\nRunning BLASTn....")
            blastdata = blast_utils.blastn_driver(dbfilepath, contigs, query, args)
        else:
            blastdata = {}

        # runs short blast with the newly created database and cleans each repeat
        if args.runshort or args.shortonly:
            if args.verbose:
                print("\nRunning BLASTn-short....")
            shortblastdata = blast_utils.short_driver(dbfilepath, contigs, query, blastdata, args)
        else:
            shortblastdata = {}

        # combines shortblast hits with blast hits
        combinedhits = blast_utils.cat_megablast_short(blastdata, shortblastdata, args)

        # data export #
        file_export.data_export_driver(contigs, combinedhits, datafilepath, filepath, query, args)

    def outfile_run(self, args):
        if args.buffer is None or args.eval_filter is None or args.writefasta is None:
            # print(args.buffer, args.eval_filter, args.writefasta)
            print("FATAL: Option specified with no value")
            exit()

        # initialize configuration for database and output directory names
        config = init_config.validate_config_load(init_config.get_config())
        outdir_name = config["outdir_name"]

        outfile = list(filter(None, [i.rstrip() for i in open(args.blast_outfile, 'r').readlines()]))
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
        file_export.data_export_driver(contigs, outfile, datafilepath, filepath, query, args)

    def masker_run(self, args):
        if args.cov_depth is None or args.mask_char is None or args.blaste is None or args.num_threads is None:
            print("FATAL: Option specified with no value")
            exit()
        if Path("%s_%smasked.fasta" % (Path(args.genome).stem, dv.PROG_NAME.lower())).exists() and not args.force:
            print("FATAL: cowardly refusing to overwrite user files, use \'--force\' to override")
            exit()
        # masker_utils.masker_str_driver(args) # DO NOT USE STRINGS IN PYTHON, LISTS ARE MUCH MUCH MUCH MORE EFFICIENT
        masker_utils.masker_list_driver(args)

    def interpret_parser(self, args):
        if args.verbose:
            print('.:| ' + dv.PROG_NAME + ' ' + dv.PROG_VERSION + ' |:.''')

        if args.subparser_id == "database" or args.subparser_id == "db":
            if args.verbose:
                print("Running %s by creating a BLAST database....\n" % dv.PROG_NAME)
            self.database_run(args)
        elif args.subparser_id == 'db_free' or args.subparser_id == 'database_free' or args.subparser_id == 'dbf':
            if args.verbose:
                print("Running %s without creating a BLAST database....\n" % dv.PROG_NAME)
            self.db_free_run(args)
        elif args.subparser_id == 'outfile':
            if args.verbose:
                print("Running %s on a BLAST output file....\n" % dv.PROG_NAME)
            self.outfile_run(args)
        elif args.subparser_id == 'masker':
            if args.verbose:
                print("Using %s to mask repetitive regions in the provided genome....\n" % dv.PROG_NAME)
            self.masker_run(args)
