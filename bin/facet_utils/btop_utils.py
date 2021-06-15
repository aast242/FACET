#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__author__ = "Alex Stewart"

import re

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


# using the qseq id's and sequences, makes a dictionary that will be used to analyze changes to seq in hits #
def make_rip_dict(query_info: list):
    ripDict = {}
    for query in query_info:
        tempPos = {}
        for pos in range(0, len(query[1])):
            tempPos[str(pos)] = {"A": 0, "G": 0, "T": 0, "C": 0, "-": 0, "match_ref": [0], "ins": [], "N": 0}
        tempPos["el_seq"] = query[1]
        ripDict[query[0]] = tempPos
    return ripDict


# takes a btop string and returns an easily readable list #
def parse_btop(btop: str):
    btop_list = list(filter(None, re.split(r'(\d+)', btop)))  # splits btop when numbers are encountered
    j = 0
    for i in range(0, len(btop_list)):
        if btop_list[j].isdigit():  # if the element is a digit
            btop_list[j] = int(btop_list[j])  # casts element to int
            j += 1
        else:  # if the element isn't a digit, must be a string
            temp = re.findall('..', btop_list[j])  # splits string list into 2's (GAGAGA becomes GA GA GA)
            if len(temp) != 1:  # if there is more than 1 pair of letters
                btop_list[j:j] = temp  # inserts the split list into the main btop list
                j += (len(temp))  # adds the length of the list to the counter
                btop_list.pop(j)  # removes the original element
            else:  # if the element is only one pair, no need to edit it; just keep going
                j += 1
    return btop_list


# adds a btop list to the rip dictionary #
def add_btop_to_dict(btop_list: list, ripDict: dict, qstart: int, qend: int, el_id: str):
    qloc = qstart - 1  # defines where the alignment starts
    stop_add = qend - 1  # defines where the alignment breaks down
    insString = ""
    while qloc <= stop_add:  # while we are in the alignment
        for ele in range(0, len(btop_list)):  # iterate through the btop
            if isinstance(btop_list[ele], int):  # this means the alignment perfectly matches the query
                for count in range(0, btop_list[ele]):  # does this loop as many times as the number in ele
                    if len(ripDict[el_id][str(qloc)]["match_ref"]) == 1:
                        ripDict[el_id][str(qloc)]["match_ref"].append(ripDict[el_id]["el_seq"][qloc])
                    ripDict[el_id][str(qloc)]["match_ref"][0] += 1  # base matches ref, increments that
                    qloc += 1  # increments the query location

            else:
                if btop_list[ele][0] != "-":  # if there isn't a gap in the query
                    # print("\n"+str(qloc))
                    # print(stop_add)
                    if len(ripDict[el_id][str(qloc)]["match_ref"]) == 1:
                        ripDict[el_id][str(qloc)]["match_ref"].append(ripDict[el_id]["el_seq"][qloc])
                    try:
                        ripDict[el_id][str(qloc)][btop_list[ele][1]] += 1  # adds one to the subject base
                    except KeyError:
                        """print("\nelement ID: %s" % el_id)
                        print("query location: %s" % qloc)
                        print("btop_list[ele][1]: %s" % btop_list[ele][1])
                        print("btop_list[ele]: %s" % btop_list[ele])
                        print("btop_list: %s" % btop_list)"""
                        pass
                    qloc += 1  # increments the query location
                else:  # there is a gap in the query
                    pass
                    # if you want to include insertions, uncomment block below
                    if isinstance(btop_list[ele + 1], str) and btop_list[ele + 1][0] != "-":
                        insString += btop_list[ele][1]
                        # print(qloc)
                        ripDict[el_id][str(qloc)]["ins"].append(insString)
                        # print(insString)
                        insString = ""
                    else:
                        insString += btop_list[ele][1]


# converts a btop sequence to a cigar sequence for storage in sam/bam files #
def btop_to_cigar(btop_list: list):
    cigar = ""
    btop_list.append([])  # adds an empty list to the end to allow for easy checks against ele+1
    count = 0
    for ele in range(0, len(btop_list) - 1):
        if isinstance(btop_list[ele], int):  # this covers matches
            cigar += str(btop_list[ele]) + "M"
        elif isinstance(btop_list[ele], str):  # this covers mismatches, deletions, and insertions
            if "-" not in btop_list[ele]:  # mismatches  (X)
                if isinstance(btop_list[ele + 1], str) and "-" not in str(btop_list[ele + 1]):
                    count += 1
                else:
                    count += 1
                    cigar += str(count) + "X"
                    count = 0
                pass
            elif btop_list[ele][0] == "-":  # deletions  (D)
                if isinstance(btop_list[ele + 1], str) and btop_list[ele + 1][0] == "-":
                    count += 1
                else:
                    count += 1
                    cigar += str(count) + "D"
                    count = 0
            elif btop_list[ele][1] == "-":  # insertions  (I)
                if isinstance(btop_list[ele + 1], str) and btop_list[ele + 1][1] == "-":
                    count += 1
                else:
                    count += 1
                    cigar += str(count) + "I"
                    count = 0
    del btop_list[-1]
    return cigar


# Runs a script to reverse RIP of ancient fungal transposable elements (C->T and G->A mutations)     #
# This parses a BTOP from BLAST alignments and uses it to reverse the point mutations and add indels #
def derip_query(query, flatcombinedhits, datafilepath, queryFilename, args):
    if args.verbose:
        print("\n\nBeginning de-RIP process....")
    derip_fasta = []
    cons_fasta = []
    query_info = [[records.id, str(records.seq)] for records in list(SeqIO.parse(query, "fasta"))]
    # print(query_info)
    ripDict = make_rip_dict(query_info)
    for hits in flatcombinedhits:  # adds each btop from each nr blast hit to dict
        btop_seq = parse_btop(hits[9])
        add_btop_to_dict(btop_seq, ripDict, hits[5], hits[6], hits[4])
    for repeats in ripDict.keys():  # goes through each repeat in the query file
        if args.verbose:
            print("\n\n" + repeats + "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
        el_len = len(ripDict[repeats]["el_seq"])  # finds element length
        derip_seq = ""  # blank derip seq
        cons_seq = ""  # blank consensus
        atog = 0  # counts A's converted back to G's
        ttoc = 0  # counts T's converted back to C's
        removedBases = 0  # counts number of bases removed
        insertions = 0  # counts number of bases inserted
        for bp in range(0, el_len):  # for each base pair in the element
            temp = list(ripDict[repeats][str(bp)].values())  # converts dict results to list
            max_cov = max(temp[0], temp[1], temp[2], temp[3], temp[4])
            if max_cov > temp[5][0]:
                if max_cov == temp[0]:
                    cons_seq += "A"
                elif max_cov == temp[1]:
                    cons_seq += "G"
                elif max_cov == temp[2]:
                    cons_seq += "T"
                elif max_cov == temp[3]:
                    cons_seq += "C"
            else:
                cons_seq += temp[5][1]
            # print(temp)
            coverage = temp[0] + temp[1] + temp[2] + temp[3] + temp[4] + temp[5][0]  # number of times base was aligned
            if 0 < coverage < 3:  # if coverage is very low (between 1 and 2), just add the consensus base
                derip_seq += temp[5][1]
            elif coverage == 0:
                derip_seq += ripDict[repeats]["el_seq"][bp]
            else:
                # doesn't remove base unless more than specified fraction of alignments are gaps
                """if float(temp[4]) / coverage >= args.deripdel:
                    removedBases += 1
                    pass"""  # change what's below to elif when you change back
                if temp[5][1] == "A":  # if the consensus is A
                    if float(temp[1] / coverage) > args.deripswitch:  # check to see if G/coverage is greater than opt
                        derip_seq += "G"  # if it is, add a G to the derip seq
                        atog += 1
                    else:
                        derip_seq += "A"  # if it isn't, add an A to the derip seq
                elif temp[5][1] == "T":  # if the consensus is T
                    if float(temp[3] / coverage) > args.deripswitch:  # check to see if T/coverage is greater than opt
                        derip_seq += "C"  # if it is, add a C to the derip seq
                        ttoc += 1
                    else:
                        derip_seq += "T"  # if it isn't, add a T to the derip seq
                else:
                    derip_seq += temp[5][1]  # if coverage is high and consensus isn't A or T, add consensus base
            """if len(temp[6]) > 0:
                maxinsertion = sorted(Counter(temp[6]).items(), key=lambda x: x[1], reverse=True)[0]
                if maxinsertion[1] >= args.deripins:  # adds insertion if there are >= opt occurrences of said insertion
                    derip_seq += maxinsertion[0]
                    cons_seq += maxinsertion[0]
                    insertions += len(maxinsertion[0])"""
        if args.verbose:
            print(str(atog) + " bases changed from A to G\n" + str(ttoc) + " bases changed from T to C\n" +
                  str(removedBases) + " bases removed\n" + str(insertions) + " bases inserted")
            print("INSERTIONS AND DELETIONS HAVE BEEN DISABLED")
        derip_fasta.append(SeqRecord(Seq(derip_seq),
                                     id=repeats + "_derip_" + str(int(args.deripswitch * 100)), description="",
                                     annotations={"molecule_type": "DNA"}))
        cons_fasta.append(SeqRecord(Seq(cons_seq),
                                    id=repeats + "_new_consensus", description="",
                                    annotations={"molecule_type": "DNA"}))
    if args.verbose:
        print("\nWriting de-ripped query file....\n")
    # export deripped sequences
    SeqIO.write(derip_fasta, datafilepath + "/" + queryFilename +
                "_derip_" + str(int(args.deripswitch * 100)) + "perc.fasta", "fasta")
    """SeqIO.write(cons_fasta, datafilepath + "/" + queryFilename +
                "_new_consensus.fasta", "fasta")"""


# Takes a btop and qseq and generates the subject sequence used to generate the BTOP
def sseq_from_btop(btop: str, qseq: str):
    mod_qseq = qseq.replace("-", "")  # maybe do this without modifying the qseq string?
    sseq = ""
    qseq_index = 0
    for i in parse_btop(btop):
        if type(i) is int:
            sseq += mod_qseq[qseq_index:qseq_index+i]
            qseq_index += i
        else:
            if i[0] == "-":  # gap in the query, base in subject; ADD SUBJ BASE, DON'T INCREMENT
                sseq += i[1]
            elif i[1] == "-":  # gap in subject, base in query; DONT ADD BASE, DO INCREMENT
                pass
                qseq_index += 1
            else:  # mismatch, take subject base
                sseq += i[1]
                qseq_index += 1
    return sseq


# Takes a btop and sseq and generates the query sequence used to generate the BTOP
def qseq_from_btop(btop: str, sseq: str):
    mod_sseq = sseq.replace("-", "")  # maybe do this without modifying the sseq string?
    qseq = ""
    sseq_index = 0
    for i in parse_btop(btop):
        if type(i) is int:
            qseq += mod_sseq[sseq_index:sseq_index+i]
            sseq_index += i
        else:
            if i[0] == "-":  # gap in the query, base in subject; DON'T ADD BASE, DO INCREMENT
                pass
                sseq_index += 1
            elif i[1] == "-":  # gap in subject, base in query; ADD BASE, DON'T INCREMENT
                qseq += i[0]
            else:  # mismatch, take query base
                qseq += i[0]
                sseq_index += 1
    return qseq


# takes a btop string and provides the reverse compliment of the btop #
def reverse_compliment(btop: str):
    btop = list(reversed(parse_btop(btop)))
    for element in range(0, len(btop)):
        # if the element is not an integer (matches), compliment the mismatch/gap string
        if not isinstance(btop[element], int):
            btop[element] = str(Seq(btop[element]).complement())
    return "".join(str(x) for x in btop)
