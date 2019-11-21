#!/usr/bin/env python3

from distutils.spawn import find_executable
import subprocess
import statistics
import uuid
import sys
import re
import os

"""

OPEN FOR BUSINESS


Mean stdev.
set_avg = statistics.mean(fvl)
set_stdev = statistics.stdev(fvl)


"""

################################################################################

def is_tool(name):
    """Check whether tool "name" is in PATH."""
    return find_executable(name) is not None


################################################################################

def bed_check_six_col_format(bed_file):
    """
    Check whether given .bed file has 6 columns.

    >>> test_bed = "test_data/test1.bed"
    >>> bed_check_six_col_format(test_bed)
    True
    >>> test_bed = "test_data/empty_file"
    >>> bed_check_six_col_format(test_bed)
    False

    """

    six_col_format = False
    with open(bed_file) as f:
        for line in f:
            cols = line.strip().split("\t")
            if len(cols) == 6:
                six_col_format = True
            break
    f.closed
    return six_col_format


################################################################################

def bed_check_unique_ids(bed_file):
    """
    Check whether .bed file (6 column format with IDs in column 4) 
    has unique column 4 IDs.
    
    >>> test_bed = "test_data/test1.bed"
    >>> bed_check_unique_ids(test_bed)
    True
    >>> test_bed = "test_data/test2.bed"
    >>> bed_check_unique_ids(test_bed)
    False

    """

    check_cmd = "cut -f 4 " + bed_file + " | sort | uniq -d"
    output = subprocess.getoutput(check_cmd)
    if output:
        return False
    else:
        return True


################################################################################

def bed_get_length_stats(bed_file):
    """
    Get length statistics for a .bed file.
    Return dictionary with following keys:
    min, max, mean, median, mode, stdev
    
    >>> test_bed = "test_data/test1.bed"
    >>> stats_dic = bed_get_length_stats(test_bed)
    >>> stats_dic["size"]
    7
    >>> stats_dic["min"]
    500
    >>> stats_dic["max"]
    1500
    >>> stats_dic["median"]
    1000
    >>> stats_dic["mode"]
    1000
    >>> stats_dic["mean"]
    1000
    >>> stats_dic["stdev"]
    288.67513459481285

    """
    # Stats dictionary.
    stats_dic = {}
    # Site lengths list.
    sll = []
    # Read in .bed file.
    with open(bed_file) as f:
        for line in f:
            cols = line.strip().split("\t")
            site_s = int(cols[1])
            site_e = int(cols[2])
            site_l = site_e - site_s
            sll.append(site_l)
    f.closed
    # Calcualte stats.
    stats_dic["size"] = len(sll)
    stats_dic["min"] = min(sll)
    stats_dic["max"] = max(sll)
    stats_dic["mean"] = statistics.mean(sll)
    stats_dic["median"] = statistics.median(sll)
    stats_dic["mode"] = statistics.mode(sll)
    stats_dic["stdev"] = statistics.stdev(sll)
    return stats_dic


################################################################################

def bed_filter_by_col5_score(in_bed_file, out_bed_file,
                             score_threshold=1,
                             generate_unique_ids=False,
                             center_sites=False,
                             ext_lr=False,
                             reverse_filter=False,
                             id_prefix="CLIP"):

    """
    Filter .bed file by column 5 score (assuming the higher the better).
    Output to new .bed file. Optionally do reverse filtering, and generate 
    new unique column 4 IDs.

    >>> in_bed_file = "test_data/test1.bed"
    >>> out_bed_file = "test_data/out.bed"
    >>> bed_filter_by_col5_score(in_bed_file, out_bed_file)
    >>> count_file_rows(out_bed_file)
    5
    >>> bed_filter_by_col5_score(in_bed_file, out_bed_file, reverse_filter=True)
    >>> count_file_rows(out_bed_file)
    4

    """

    # Output .bed file.
    OUTBED = open(out_bed_file,"w")

    # New site IDs.
    site_id_pref = id_prefix
    c_sites = 0

    # Open input .bed file.
    with open(in_bed_file) as f:
        for line in f:
            cols = line.strip().split("\t")
            seq_id = cols[0]
            site_s = int(cols[1])
            site_e = int(cols[2])
            site_id = cols[3]
            site_sc = float(cols[4])
            site_pol = cols[5]
            # Sanity checking .bed file.
            if site_s >= site_e:
                print ("ERROR: invalid region coordinates in .bed file \"%s\" (start >= end: %i >= %i)" % (in_bed_file, site_s, site_e))
                sys.exit()
            if site_s < 0 or site_e < 1:
                print ("ERROR: invalid region coordinates in .bed file \"%s\" (start < 0 or end < 1)" % (in_bed_file))
                sys.exit()
            # Filter by score.
            if reverse_filter:
                if site_sc > score_threshold:
                    continue
            else:
                if site_sc < score_threshold:
                    continue
            # Update start + end positions.
            new_s = site_s
            new_e = site_e
            # Center sites (get center position).
            if center_sites:
                new_e = get_center_position(site_s, site_e)
                new_s = new_e - 1
            # If site extension is specified.
            if ext_lr:
                new_s = new_s - ext_lr
                new_e = new_e + ext_lr
            # New site ID.
            if generate_unique_ids:
                c_sites += 1
                site_id = "%s_%i" % (site_id_pref, c_sites)
            # Output to new file.
            OUTBED.write("%s\t%i\t%i\t%s\t%f\t%s\n" % (seq_id,new_s,new_e,site_id,site_sc,site_pol))
    f.closed


################################################################################

def get_center_position(start, end):
    """
    Get center position (1-based), given a (genomic) start (0-based) and 
    end coordinate (1-based).

    >>> get_center_position(10, 11)
    11
    >>> get_center_position(1000,2000)
    1501
    >>> get_center_position(11, 20)
    17

    """
    # If region has length of 1, return end position.
    center_pos = end
    # Otherwise calculate new center position.
    if not end - start == 1:
        center_pos = round( ( (end - start) / 2 ) + start ) + 1
    return center_pos


################################################################################

def count_file_rows(in_file):
    """
    Count number of file rows for given input file.
    
    >>> test_file = "test_data/test1.bed"
    >>> count_file_rows(test_file)
    7
    >>> test_file = "test_data/empty_file"
    >>> count_file_rows(test_file)
    0

    """
    check_cmd = "cat " + in_file + " | wc -l"
    output = subprocess.getoutput(check_cmd)
    row_count = int(output.strip())
    return row_count


################################################################################

def read_ids_into_dic(ids_file,
                      ids_dic=False):
    """
    Read in IDs list file, where each line stores one ID.

    >>> test_ids_file = "test_data/test.ids"
    >>> ids_dic = read_ids_into_dic(test_ids_file)
    >>> print(ids_dic)
    {'CLIP1': 1, 'CLIP2': 1, 'CLIP3': 1}

    """
    if not ids_dic:
        ids_dic = {}
    # Read in file content.
    with open(ids_file) as f:
        for line in f:
            row_id = line.strip()
            ids_dic[row_id] = 1
    f.closed
    return ids_dic


################################################################################

def read_fasta_into_dic(fasta_file,
                        seqs_dic=False,
                        ids_dic=False,
                        skip_n_seqs=True):
    """
    Read in FASTA sequences, store in dictionary and return dictionary.
    
    >>> test_fasta = "test_data/test.fa"
    >>> d = read_fasta_into_dic(test_fasta)
    >>> print(d)
    {'seq1': 'acguACGUacgu', 'seq2': 'ugcaUGCAugcaACGUacgu'}
    >>> test_fasta = "test_data/test2.fa"
    >>> d = read_fasta_into_dic(test_fasta)
    >>> print(d)
    {}

    """
    if not seqs_dic:
        seqs_dic = {}
    seq_id = ""
    seq = ""
    # Go through FASTA file, extract sequences.
    with open(fasta_file) as f:
        for line in f:
            if re.search(">.+", line):
                m = re.search(">(.+)", line)
                seq_id = m.group(1)
                # If there is a ".", take only first part of header.
                # This assumes ENSEMBL header format ">ENST00000631435.1 cdna ..."
                if re.search(".+\..+", seq_id):
                    m = re.search("(.+?)\..+", seq_id)
                    seq_id = m.group(1)
                if seq_id in seqs_dic:
                    print ("ERROR: non-unique FASTA header \"%s\" in \"%s\"" % (seq_id, fasta_file))
                    sys.exit()
                else:
                    if ids_dic:
                        if seq_id in ids_dic:
                            seqs_dic[seq_id] = ""
                    else:
                        seqs_dic[seq_id] = ""
            elif re.search("[ACGTUN]+", line, re.I):
                m = re.search("([ACGTUN]+)", line, re.I)
                if seq_id in seqs_dic:
                    # Convert to RNA, concatenate sequence.
                    seqs_dic[seq_id] += m.group(1).replace("T","U").replace("t","u")
                # If sequences with N nucleotides should be skipped.
                if skip_n_seqs:
                    if "n" in m.group(1) or "N" in m.group(1):
                        print ("WARNING: sequence with seq_id \"%s\" in file \"%s\" contains N nucleotides. Discarding sequence ... " % (seq_id, fasta_file))
                        del seqs_dic[seq_id]
                        continue
    f.closed
    return seqs_dic


################################################################################



"""
    # Generate .tmp file.
    random_id = uuid.uuid1()
    tmp_bed = str(random_id) + ".tmp.bed"
    random_id = uuid.uuid1()
    tmp_tab = str(random_id) + ".tmp.tab"
"""

