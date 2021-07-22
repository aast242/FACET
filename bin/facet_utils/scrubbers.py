from pprint import pprint
import argparse

parser = argparse.ArgumentParser(description='Program: BLASTrec (BLAST redundancy cleaner)\n'
                                             'Version: 1.5\n',
                                 formatter_class=argparse.RawDescriptionHelpFormatter,
                                 usage='%(prog)s contigs query [options]')
parser.add_argument("--buffer", default=15, type=int,
                    help='length in bp of a buffer to help with finding correct hits [default: 15]', nargs='?')
args = parser.parse_args()

cleaned = [['Chr6', 32030, 39301, 'lazarus_original', 1, 7284, 'minus', 93.249],
           ['Chr6', 280505, 280788, 'lazarus_original', 7003, 7288, 'plus', 90.21],
           ['Chr6', 1210488, 1217770, 'lazarus_original', 1, 7286, 'plus', 91.72],
           ['Chr6', 1228376, 1232815, 'lazarus_original', 1, 4453, 'plus', 86.437],
           ['Chr6', 1240333, 1240953, 'lazarus_original', 4500, 5115, 'plus', 84.244],
           ['Chr6', 2125100, 2127585, 'lazarus_original', 4780, 7282, 'minus', 89.505],
           ['Chr6', 2623322, 2623607, 'lazarus_original', 6997, 7282, 'plus', 87.762],
           ['Chr6', 2623343, 2627257, 'lazarus_original', 9, 3951, 'plus', 89.291],
           ['Chr6', 2627789, 2631108, 'lazarus_original', 3945, 7284, 'plus', 86.936],
           ['Chr6', 2633894, 2639547, 'lazarus_original', 1, 5640, 'plus', 86.569],
           ['Chr6', 2639539, 2641050, 'lazarus_original', 5776, 7288, 'plus', 87.657],
           ['Chr6', 3281903, 3282188, 'lazarus_original', 7003, 7288, 'minus', 90.559],
           ['Chr6', 3405330, 3412614, 'lazarus_original', 1, 7288, 'plus', 92.35],
           ['Chr6', 5787342, 5792572, 'lazarus_original', 93, 5339, 'minus', 75.878],
           ['Chr6', 5792569, 5792847, 'lazarus_original', 1, 279, 'minus', 89.964]]

dirty = [['Chr6', 32028, 32304, 'lazarus_original', 1, 277, 'minus', 88.809],
         ['Chr6', 32030, 39301, 'lazarus_original', 1, 7284, 'minus', 93.249],
         ['Chr6', 39023, 39313, 'lazarus_original', 6999, 7288, 'minus', 91.409],
         ['Chr6', 280505, 280788, 'lazarus_original', 7003, 7288, 'plus', 90.21],
         ['Chr6', 280512, 280798, 'lazarus_original', 1, 287, 'plus', 91.349],
         ['Chr6', 1210482, 1210766, 'lazarus_original', 7004, 7288, 'plus', 88.07],
         ['Chr6', 1210488, 1217770, 'lazarus_original', 1, 7286, 'plus', 91.72],
         ['Chr6', 1217494, 1217782, 'lazarus_original', 1, 287, 'plus', 85.813],
         ['Chr6', 1228372, 1228648, 'lazarus_original', 7006, 7282, 'plus', 83.032],
         ['Chr6', 1228376, 1232815, 'lazarus_original', 1, 4453, 'plus', 86.437],
         ['Chr6', 1240333, 1240953, 'lazarus_original', 4500, 5115, 'plus', 84.244],
         ['Chr6', 2125100, 2127585, 'lazarus_original', 4780, 7282, 'minus', 89.505],
         ['Chr6', 2623322, 2623607, 'lazarus_original', 6997, 7282, 'plus', 87.762],
         ['Chr6', 2623343, 2627257, 'lazarus_original', 9, 3951, 'plus', 89.291],
         ['Chr6', 2627789, 2631108, 'lazarus_original', 3945, 7284, 'plus', 86.936],
         ['Chr6', 2630834, 2631110, 'lazarus_original', 1, 277, 'plus', 85.199],
         ['Chr6', 2633883, 2634172, 'lazarus_original', 6999, 7288, 'plus', 89.655],
         ['Chr6', 2633894, 2639547, 'lazarus_original', 1, 5640, 'plus', 86.569],
         ['Chr6', 2639539, 2641050, 'lazarus_original', 5776, 7288, 'plus', 87.657],
         ['Chr6', 2640772, 2641053, 'lazarus_original', 1, 282, 'plus', 91.489],
         ['Chr6', 3281903, 3282188, 'lazarus_original', 7003, 7288, 'minus', 90.559],
         ['Chr6', 3281903, 3282181, 'lazarus_original', 1, 279, 'minus', 90.681],
         ['Chr6', 3405324, 3405608, 'lazarus_original', 7004, 7288, 'plus', 88.772],
         ['Chr6', 3405330, 3412614, 'lazarus_original', 1, 7288, 'plus', 92.35],
         ['Chr6', 3412336, 3412616, 'lazarus_original', 1, 281, 'plus', 92.171],
         ['Chr6', 5787342, 5792572, 'lazarus_original', 93, 5339, 'minus', 75.878],
         ['Chr6', 5792569, 5792847, 'lazarus_original', 1, 279, 'minus', 89.964],
         ['Chr6', 5792569, 5792853, 'lazarus_original', 7004, 7288, 'minus', 86.667]]


# takes a list of lists and retains elements that have no complete overlap with other elements  #
# Input looks like: [[9,20], [2,3],[5,7],[2,8],[1,9],[1,10],[8,15]], returns [[1,10],[9,20]]    #
# index1 and index2 allow you to have lists that have more than two elements but you still want #
# to compare them as if they did. Input ranges must be [min, max]                               #
def clean_list(cluttered, index1, index2, pident_index, args):
    # sorts cluttered list by first part of each sublist
    cluttered = sorted(cluttered, key=lambda x: x[index1])
    # reduces function calls #
    clutteredremove = cluttered.remove
    exitcon = True
    while exitcon:  # cleans the list until no ranges contained within other ranges remain
        exitcon = False  # sets exit condition to false each time the loop starts
        # for loop checks to see if either list is containing the other #
        for i in range(1, len(cluttered)):
            if cluttered[i][index1] >= cluttered[i - 1][index1] - args.buffer and \
                    cluttered[i][index2] <= cluttered[i - 1][index2] + args.buffer:
                if pident_index >= 0:  # checks percent identity to retain matches that match the best
                    # checks to see if hits are approx. the same length #
                    if -args.buffer <= (cluttered[i][index2] - cluttered[i][index1]) - \
                            (cluttered[i - 1][index2] - cluttered[i - 1][index1]) <= args.buffer:
                        # removes the hit with lower pident #
                        if cluttered[i][pident_index] - cluttered[i - 1][pident_index] > 0:
                            clutteredremove(cluttered[i - 1])
                        else:
                            clutteredremove(cluttered[i])
                    else:
                        clutteredremove(cluttered[i])
                    exitcon = True
                    break
                else:  # if we don't check by pident, remove the shorter element
                    if (cluttered[i][index2] - cluttered[i][index1]) > \
                            (cluttered[i - 1][index2] - cluttered[i - 1][index1]):
                        clutteredremove(cluttered[i - 1])
                    else:
                        clutteredremove(cluttered[i])
                    exitcon = True
                    break
            elif cluttered[i - 1][index1] >= cluttered[i][index1] - args.buffer and \
                    cluttered[i - 1][index2] <= cluttered[i][index2] + args.buffer:
                if pident_index >= 0:  # checks percent identity to retain matches that match the best
                    # checks to see if hits are approx. the same length #
                    if -args.buffer <= (cluttered[i][index2] - cluttered[i][index1]) - \
                            (cluttered[i - 1][index2] - cluttered[i - 1][index1]) <= args.buffer:
                        # removes the hit with lower pident #
                        if cluttered[i][pident_index] - cluttered[i - 1][pident_index] > 0:
                            clutteredremove(cluttered[i - 1])
                        else:
                            clutteredremove(cluttered[i])
                    else:
                        clutteredremove(cluttered[i - 1])
                    exitcon = True
                    break
                else:  # if we don't check by pident, remove the shorter element
                    if (cluttered[i][index2] - cluttered[i][index1]) > \
                            (cluttered[i - 1][index2] - cluttered[i - 1][index1]):
                        clutteredremove(cluttered[i - 1])
                    else:
                        clutteredremove(cluttered[i])
                    exitcon = True
                    break

    exitcon = True
    while exitcon:  # cleans the list until no ranges contained within surrounding ranges remain
        exitcon = False  # sets exit condition to false each time the loop starts
        for i in range(1, len(cluttered) - 1):
            # if the range in question is overlapping both the one before and after it
            if max(cluttered[i][index1], cluttered[i - 1][index1]) <= min(cluttered[i][index2],
                                                                          cluttered[i - 1][index2]):
                if max(cluttered[i][index1], cluttered[i + 1][index1]) <= min(cluttered[i][index2],
                                                                              cluttered[i + 1][index2]):
                    # if the range in question is contained by the range before and after
                    if cluttered[i - 1][index2] >= cluttered[i + 1][index1]:
                        exitcon = True
                        clutteredremove(cluttered[i])
                        break
    return cluttered


def seth_scrubber(cluttered, sindex, eindex, pidentindex, args):
    starts, ends, pidents, active, limbo, saved = [], [], [], [], [], []
    # creates start, stop, and pident lists [[id, start], [1, 32304], ...]
    for alignment in range(0, len(cluttered)):
        starts.append([alignment, cluttered[alignment][sindex]])
        ends.append([alignment, cluttered[alignment][eindex]])
        pidents.append([alignment, cluttered[alignment][pidentindex]])
    # sorts ends by stop site (starts are already sorted by BLASTn)
    ends = sorted(ends, key=lambda x: x[1])
    ###### TAKE CARE OF START/END ON SAME VALUE AFTER THIS ######

    scount, ecount, activestartindex = 1, 0, 0
    # adds first start to the active list
    active.append(starts[0])
    pprint(starts)
    print()
    pprint(ends)
    while len(ends) > 0:
        # if a stop is the next significant entry
        if scount == len(starts) or starts[scount][1] > ends[ecount][1]:
            # checks to see where the significant end ID occurs in the active list
            activestartindex = next(i for i, v in enumerate(active) if v[0] == ends[ecount][0])

            # creates a new active list where each of the new elements except the active one are artificially
            # extended by the buffer value, sorts it, and returns where the active ID shows up in the list
            bufferasi = next(i for i, v in enumerate(
                sorted([[pair[0], pair[1] - args.buffer] if pair[0] != ends[ecount][0] else
                        [pair[0], pair[1]] for pair in active], key=lambda x: x[1])) if v[0] == ends[ecount][0])
            # removes hit if the indices differ
            if activestartindex != bufferasi:
                print("%s %s" % (activestartindex, bufferasi))
                limbo.append([ends[ecount][0], active[activestartindex][1], ends[ecount][1]])
                del ends[ecount]
                del active[activestartindex]
                continue
            """
            limbocount = 0
            while limbocount < len(limbo):
                # length check #
                if abs(active[activestartindex][1]-limbo[limbocount][1]) + abs(ends[ecount][1]-limbo[limbocount][2])\
                        <= 2*args.buffer:
                    pass
            """
            # if it's at the beginning, add the index to the saved list, send the range to limbo, and remove the range
            if activestartindex == 0:
                saved.append(ends[ecount][0])
                limbo.append([ends[ecount][0], active[activestartindex][1], ends[ecount][1]])
            del ends[ecount]
            del active[activestartindex]

            # clean up the limbo list
            limbocount = 0
            while limbocount < len(limbo) and len(ends) > 0:
                if ends[ecount][1] > limbo[limbocount][2] + args.buffer:
                    del limbo[limbocount]
                else:
                    limbocount += 1
        # if a start is the next significant entry
        else:
            # adds start to the active list, increments the start count
            active.append(starts[scount])
            scount += 1
    print(saved)


# takes a list of lists and retains elements that have no complete overlap with other elements  #
# Input looks like: [[9,20], [2,3],[5,7],[2,8],[1,9],[1,10],[8,15]], returns [[1,10],[9,20]]    #
# index1 and index2 allow you to have lists that have more than two elements but you still want #
# to compare them as if they did. Input ranges must be [min, max]                               #
# TODO: add argument for checking by pident, always have pident_index, and if hits are exactly the same length,
# TODO: take the one with higher pident
def clean_list_backup(cluttered, index1, index2, pident_index, args):
    # sorts cluttered list by first part of each sublist
    cluttered = sorted(cluttered, key=lambda x: x[index1])
    # reduces function calls #
    clutteredremove = cluttered.remove
    exitcon = True
    while exitcon:  # cleans the list until no ranges contained within other ranges remain
        exitcon = False  # sets exit condition to false each time the loop starts
        # for loop checks to see if either list is containing the other #
        for i in range(1, len(cluttered)):
            if cluttered[i][index1] >= cluttered[i - 1][index1] - args.buffer and \
                    cluttered[i][index2] <= cluttered[i - 1][index2] + args.buffer:
                if pident_index >= 0:  # checks percent identity to retain matches that match the best
                    # checks to see if hits are approx. the same length #
                    if -args.buffer <= (cluttered[i][index2] - cluttered[i][index1]) - \
                            (cluttered[i - 1][index2] - cluttered[i - 1][index1]) <= args.buffer:
                        # removes the hit with lower pident #
                        if cluttered[i][pident_index] - cluttered[i - 1][pident_index] > 0:
                            clutteredremove(cluttered[i - 1])
                        # if the pidents are the same, keep longer hit
                        elif cluttered[i][pident_index] - cluttered[i - 1][pident_index] == 0:
                            if (cluttered[i][index2] - cluttered[i][index1]) > \
                                    (cluttered[i - 1][index2] - cluttered[i - 1][index1]):
                                clutteredremove(cluttered[i - 1])
                            else:
                                clutteredremove(cluttered[i])
                        else:
                            clutteredremove(cluttered[i])
                    else:
                        clutteredremove(cluttered[i])
                    exitcon = True
                    break
                else:  # if we don't check by pident, remove the shorter element
                    if (cluttered[i][index2] - cluttered[i][index1]) > \
                            (cluttered[i - 1][index2] - cluttered[i - 1][index1]):
                        clutteredremove(cluttered[i - 1])
                    else:
                        clutteredremove(cluttered[i])
                    exitcon = True
                    break
            elif cluttered[i - 1][index1] >= cluttered[i][index1] - args.buffer and \
                    cluttered[i - 1][index2] <= cluttered[i][index2] + args.buffer:
                if pident_index >= 0:  # checks percent identity to retain matches that match the best
                    # checks to see if hits are approx. the same length #
                    if -args.buffer <= (cluttered[i][index2] - cluttered[i][index1]) - \
                            (cluttered[i - 1][index2] - cluttered[i - 1][index1]) <= args.buffer:
                        # removes the hit with lower pident #
                        if cluttered[i][pident_index] - cluttered[i - 1][pident_index] > 0:
                            clutteredremove(cluttered[i - 1])
                        # if pidents are the same, keep the longer hit
                        elif cluttered[i][pident_index] - cluttered[i - 1][pident_index] == 0:
                            if (cluttered[i][index2] - cluttered[i][index1]) > \
                                    (cluttered[i - 1][index2] - cluttered[i - 1][index1]):
                                clutteredremove(cluttered[i - 1])
                            else:
                                clutteredremove(cluttered[i])
                        else:
                            clutteredremove(cluttered[i])
                    else:
                        clutteredremove(cluttered[i - 1])
                    exitcon = True
                    break
                else:  # if we don't check by pident, remove the shorter element
                    if (cluttered[i][index2] - cluttered[i][index1]) > \
                            (cluttered[i - 1][index2] - cluttered[i - 1][index1]):
                        clutteredremove(cluttered[i - 1])
                    else:
                        clutteredremove(cluttered[i])
                    exitcon = True
                    break

    exitcon = True
    while exitcon:  # cleans the list until no ranges contained within surrounding ranges remain
        exitcon = False  # sets exit condition to false each time the loop starts
        for i in range(1, len(cluttered) - 1):
            # if the range in question is overlapping both the one before and after it
            if max(cluttered[i][index1], cluttered[i - 1][index1]) <= min(cluttered[i][index2],
                                                                          cluttered[i - 1][index2]):
                if max(cluttered[i][index1], cluttered[i + 1][index1]) <= min(cluttered[i][index2],
                                                                              cluttered[i + 1][index2]):
                    # if the range in question is contained by the range before and after
                    if cluttered[i - 1][index2] >= cluttered[i + 1][index1]:
                        exitcon = True
                        clutteredremove(cluttered[i])
                        break
    return cluttered

seth_scrubber(dirty, 1, 2, 7, args)
