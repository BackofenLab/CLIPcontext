#!/usr/bin/env python3

from distutils.spawn import find_executable
import subprocess
import statistics
import random
import gzip
import uuid
import sys
import re
import os

"""

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
~~~~~~ OPEN FOR BUSINESS ~~~~~~
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


AuthoR: uhlm [at] informatik.uni-freiburg.de


~~~~~~~~~~~~~
Run doctests
~~~~~~~~~~~~~

python3 -m doctest cliplib.py



"""

################################################################################

def is_tool(name):
    """Check whether tool "name" is in PATH."""
    return find_executable(name) is not None


################################################################################

def dir_get_files(file_dir,
                  file_ending=False):
    """
    Return list of files from given directory file_dir.
    E.g. file_ending="bed" to filter for .bed files.

    >>> test_dir = "test_data"
    >>> dir_get_files(test_dir, file_ending="profile")
    ['test2.profile', 'test.profile']

    """

    from os import listdir
    from os.path import isfile, join
    dir_files = [f for f in listdir(file_dir) if isfile(join(file_dir, f))]
    assert dir_files, "ERROR: given directory \"%s\" contains no files" %(file_dir)
    # If filter for file ending true.
    if file_ending:
        new_files = []
        for df in dir_files:
            if re.search(".+\.%s" %(file_ending), df):
                new_files.append(df)
        assert new_files, "ERROR: no files left after filtering by file ending \"%s\"" %(file_ending)
        return new_files
    else:
        return dir_files


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

def bed_get_region_id_scores(in_bed_file):
    """
    Read in .bed file, and store scores for each region in dictionary
    (unique column 4 ID and column 5 score have to be present).
    Return dictionary with mappings region ID -> region score

    >>> test_bed = "test_data/test3.bed"
    >>> bed_get_region_id_scores(test_bed)
    {'CLIP2': 2.57, 'CLIP1': 1.58, 'CLIP3': 3.11}

    """
    id2sc_dic = {}
    # Open input .bed file.
    with open(in_bed_file) as f:
        for line in f:
            cols = line.strip().split("\t")
            site_id = cols[3]
            site_sc = float(cols[4])
            id2sc_dic[site_id] = site_sc
    f.closed
    return id2sc_dic


################################################################################

def bed_map_region_id_to_seq_id(in_bed_file):
    """
    Read in .bed file, and store for each region ID (column 4) the sequence
    ID (column 1)
    Return dictionary with mappings region ID -> sequence ID

    >>> test_bed = "test_data/test3.bed"
    >>> bed_map_region_id_to_seq_id(test_bed)
    {'CLIP2': 'chr1', 'CLIP1': 'chr2', 'CLIP3': 'chr1'}

    """
    regid2seqid_dic = {}
    # Open input .bed file.
    with open(in_bed_file) as f:
        for line in f:
            cols = line.strip().split("\t")
            seq_id = cols[0]
            site_id = cols[3]
            regid2seqid_dic[site_id] = seq_id
    f.closed
    return regid2seqid_dic


################################################################################

def bed_process_bed_file(in_bed_file, out_bed_file,
                         score_thr=None,
                         min_len=False,
                         max_len=False,
                         generate_unique_ids=False,
                         center_sites=False,
                         ext_lr=False,
                         seq_len_dic=False,
                         siteids2keep_dic=False,
                         seqids2keep_dic=False,
                         siteseqids2keep_dic=False,
                         zero_scores=False,
                         int_whole_nr=True,
                         rev_filter=False,
                         id_prefix="CLIP"):

    """
    Process .bed file in various ways:
    - Filter by region length (min_len, max_len) or region score (column 5)
      (score_thr). By default no score or length filtering is applied.
    - Option to reverse-filter scores (the lower score the better)
    - Center regions (center_sites=True)
    - Extend sites up- downstream (ext_lr=value)
    - Generate new region IDs (column 4, generate_unique_ids=True),
      optionally providing an id_prefix
    - Filter by given dictionary of region IDs (ids2keep_dic)
    - Print "0" scores to column 5 (zero_scores)

    Output processed .bed file (in_bed_file) to new bed file (out_bed_file)

    >>> in_bed = "test_data/test1.bed"
    >>> out_bed = "test_data/out.bed"
    >>> bed_process_bed_file(in_bed, out_bed, score_thr=1)
    >>> count_file_rows(out_bed)
    5
    >>> bed_process_bed_file(in_bed, out_bed, rev_filter=True, score_thr=1)
    >>> count_file_rows(out_bed)
    4
    >>> in_bed = "test_data/test3.bed"
    >>> out_bed = "test_data/out.bed"
    >>> out_bed_exp = "test_data/test3_out.centered_zero_sc.bed"
    >>> bed_process_bed_file(in_bed, out_bed, zero_scores=True, center_sites=True)
    >>> diff_two_files_identical(out_bed, out_bed_exp)
    True

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
            site_l = site_e - site_s
            # Sanity checking .bed file.
            assert site_s < site_e, "ERROR: invalid region coordinates in .bed file \"%s\" (start >= end: %i >= %i)" % (in_bed_file, site_s, site_e)
            assert site_s >= 0 and site_e >= 1, "ERROR: invalid region coordinates in .bed file \"%s\" (start < 0 or end < 1)" % (in_bed_file)
            # Filter by IDs to keep dictionary.
            if siteids2keep_dic:
                if not site_id in siteids2keep_dic:
                    continue
            if seqids2keep_dic:
                if not seq_id in seqids2keep_dic:
                    continue
            if siteseqids2keep_dic:
                if site_id in siteseqids2keep_dic:
                    if not seq_id == siteseqids2keep_dic[site_id]:
                        continue
            # Filter by score.
            if score_thr is not None:
                if rev_filter:
                    if site_sc > score_thr:
                        continue
                else:
                    if site_sc < score_thr:
                        continue
            # Check whether score is whole number.
            if int_whole_nr:
                if not site_sc % 1:
                    site_sc = int(site_sc)
            # Filter by minimum site length.
            if min_len:
                if min_len > site_l:
                    continue
            # Filter by maximum site length.
            if max_len:
                if max_len < site_l:
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
                if new_s < 0:
                    new_s = 0
            # New site ID.
            if generate_unique_ids:
                c_sites += 1
                site_id = "%s_%i" % (site_id_pref, c_sites)
            new_sc = str(site_sc)
            if zero_scores:
                new_sc = "0"
            if seq_len_dic:
                assert seq_id in seq_len_dic, "sequence ID \"%s\" missing in given sequence lengths dictionary" %(seq_id)
                if new_e > seq_len_dic[seq_id]:
                    new_e = seq_len_dic[seq_id]
            # Output to new file.
            OUTBED.write("%s\t%i\t%i\t%s\t%s\t%s\n" % (seq_id,new_s,new_e,site_id,new_sc,site_pol))
    f.closed
    OUTBED.close()


################################################################################

def bed_read_rows_into_dic(in_bed,
                           two_ids_dic=False):
    """
    Read in .bed file rows into dictionary.
    Mapping is region ID -> bed row.
    two_ids_dic  : dictionary with site ID -> sequence ID, used for filtering sites.
                   Thus, row has to have both site and sequence ID to be kept.

    >>> test_bed = "test_data/test4.bed"
    >>> bed_read_rows_into_dic(test_bed)
    {'CLIP1': 'chr1\\t10\\t20\\tCLIP1\\t0\\t+', 'CLIP2': 'chr1\\t30\\t40\\tCLIP2\\t0\\t+'}

    """
    id2row_dic = {}
    with open(in_bed) as f:
        for line in f:
            row = line.strip()
            cols = line.strip().split("\t")
            seq_id = cols[0]
            site_id = cols[3]
            if two_ids_dic:
                if site_id in two_ids_dic:
                    if not seq_id == two_ids_dic[site_id]:
                        continue
                else:
                    continue
            id2row_dic[site_id] = row
    f.closed
    assert id2row_dic, "No IDs read into dictionary (input file \"%s\" empty or malformatted?)" % (in_bed)
    return id2row_dic


################################################################################

def get_seq_lengths_from_seqs_dic(seqs_dic):
    """
    Given a dictionary of sequences, return dictionary of sequence lengths.
    Mapping is sequence ID -> sequence length.
    """
    seq_len_dic = {}
    assert seqs_dic, "sequence dictionary seems to be empty"
    for seq_id in seqs_dic:
        seq_l = len(seqs_dic[seq_id])
        seq_len_dic[seq_id] = seq_l
    return seq_len_dic


################################################################################

def bed_get_region_lengths(bed_file):
    """
    Read in .bed file, store and return region lengths in dictionary.
    key   :  region ID (.bed col4)
    value :  region length (.bed col3-col2)

    >>> test_file = "test_data/test4.bed"
    >>> bed_get_region_lengths(test_file)
    {'CLIP1': 10, 'CLIP2': 10}

    """
    id2len_dic = {}
    with open(bed_file) as f:
        for line in f:
            row = line.strip()
            cols = line.strip().split("\t")
            site_s = int(cols[1])
            site_e = int(cols[2])
            site_id = cols[3]
            site_l = site_e - site_s
            assert site_id not in id2len_dic, "column 4 IDs not unique in given .bed file \"%s\"" %(bed_file)
            id2len_dic[site_id] = site_l
    f.closed
    assert id2len_dic, "No IDs read into dictionary (input file \"%s\" empty or malformatted?)" % (in_bed)
    return id2len_dic


################################################################################

def extract_transcript_sequences(bed_dic, seq_dic,
                                 ext_lr=False,
                                 full_hits_only=False):
    """
    Given a dictionary with bed regions (region ID -> BED row) and a
    sequence dictionary (Sequence ID -> sequence), extract the BED region
    sequences and return in new dictionary (region ID -> region sequence).
    Optionally, extend regions by ext_lr nt (up- and downstream).
    In case full extension is not possible, use maximum extension possible.
    Set full_hits_only=True to only recover full hits.

    >>> seq_dic = {"T1" : "AAAACCCCGGGGTTTT", "T2" : "ATATACACAGAGCGCGCTCTGTGT"}
    >>> bed_dic = {"S1" : "T1\\t4\\t8\\tS1\\t0\\t+", "S2" : "T2\\t6\\t8\\tS2\\t0\\t+"}
    >>> extract_transcript_sequences(bed_dic, seq_dic, ext_lr=2)
    {'S1': 'AACCCCGG', 'S2': 'ACACAG'}
    >>> extract_transcript_sequences(bed_dic, seq_dic, ext_lr=5, full_hits_only=True)
    {'S2': 'TATACACAGAGC'}

    """
    id2seq_dic = {}
    # Process .bed regions.
    for reg_id in bed_dic:
        cols = bed_dic[reg_id].split("\t")
        seq_id = cols[0]
        reg_s = int(cols[1])
        reg_e = int(cols[2])
        assert seq_id in seq_dic, "sequence ID \"%s\" not found in given sequence dictionary" %(seq_id)
        seq = seq_dic[seq_id]
        # Update region borders.
        new_s = reg_s
        new_e = reg_e
        exp_l = new_e - new_s
        # Adjust if given start or end is out of bounds.
        if new_s < 0:
            new_s = 0
        if new_e > len(seq):
            new_e = len(seq)
        # If region should be extended up- and downstream by ext_lr.
        if ext_lr:
            new_s = new_s - ext_lr
            new_e = reg_e + ext_lr
            exp_l = new_e - new_s
            # If start or end is out of bounds after extension.
            if new_s < 0:
                new_s = 0
            if new_e > len(seq):
                new_e = len(seq)
        reg_seq = seq[new_s:new_e]
        reg_l = len(reg_seq)
        if full_hits_only:
            if not reg_l == exp_l:
                continue
        id2seq_dic[reg_id] = reg_seq
    return id2seq_dic


################################################################################

def bed_extract_sequences_from_2bit(in_bed, out_fa, in_2bit):
    """
    Extract sequences from genome (provide genome .2bit file).
    twoBitToFa executable needs to be in PATH. Store extracted
    sequences in out_fa.

    """
    # Check for twoBitToFa.
    assert is_tool("twoBitToFa"), "twoBitToFa not in PATH"

    # Run twoBitToFa and check.
    check_cmd = "twoBitToFa -noMask -bed=" + in_bed + " " + in_2bit + " " + out_fa
    output = subprocess.getoutput(check_cmd)
    error = False
    if output:
        error = True
    assert error == False, "twoBitToFa is complaining:\n%s\n%s" %(check_cmd, output)


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

def bed_sort_file(in_bed, out_bed,
                  custom_params_str=False):
    """
    Sort input .bed file and output sorted .bed.
    Use command line sort, with default parameters: "-k1,1 -k2,2n"

    """
    params_str = "-k1,1 -k2,2n"
    if custom_params_str:
        params_str = custom_params_str
    check_cmd = "sort " + params_str + " " + in_bed + " > " + out_bed
    output = subprocess.getoutput(check_cmd)
    error = False
    if output:
        error = True
    assert error == False, "sort is complaining:\n%s\n%s" %(check_cmd, output)


################################################################################

def bed_merge_file(in_bed, out_bed,
                   custom_params_str=False):
    """
    Use mergeBed from bedtools to merge overlapping .bed entries, storing
    the region IDs to later pick one region for each set of overlapping
    regions.

    >>> in_bed = "test_data/test.sorted.bed"
    >>> out_bed = "test_data/test.sorted.merged.tmp.bed"
    >>> out_exp_bed = "test_data/test.sorted.merged.exp.bed"
    >>> bed_merge_file(in_bed, out_bed)
    >>> diff_two_files_identical(out_bed, out_exp_bed)
    True

    """
    # Check for bedtools.
    assert is_tool("bedtools"), "bedtools not in PATH"
    # Parameter string.
    params_str = '-s -c 4 -o distinct -delim ";"'
    if custom_params_str:
        params_str = custom_params_str
    check_cmd = "mergeBed -i " + in_bed + " " + params_str + " > " + out_bed
    output = subprocess.getoutput(check_cmd)
    error = False
    if output:
        error = True
    assert error == False, "mergeBed is complaining:\n%s\n%s" %(check_cmd, output)


################################################################################

def graphprot_predictions_read_in_ids(predictions_file):
    """
    Given a GraphProt .predictions file, read in column 1 IDs in order
    appearing in file, and store in list.

    >>> test_file = "test_data/test.predictions"
    >>> graphprot_predictions_read_in_ids(test_file)
    ['SERBP1_K562_rep01_2175', 'SERBP1_K562_rep01_544', 'SERBP1_K562_rep01_316']

    """
    ids_list = []
    with open(predictions_file) as f:
        for line in f:
            cols = line.strip().split("\t")
            seq_id = cols[0]
            ids_list.append(seq_id)
    f.close()
    return ids_list


################################################################################

def fasta_read_in_ids(fasta_file):
    """
    Given a .fa file, read in header IDs in order appearing in file,
    and store in list.

    >>> test_file = "test_data/test3.fa"
    >>> fasta_read_in_ids(test_file)
    ['SERBP1_K562_rep01_544', 'SERBP1_K562_rep02_709', 'SERBP1_K562_rep01_316']

    """
    ids_list = []
    with open(fasta_file) as f:
        for line in f:
            if re.search(">.+", line):
                m = re.search(">(.+)", line)
                seq_id = m.group(1)
                ids_list.append(seq_id)
    f.close()
    return ids_list


################################################################################

def graphprot_profile_extract_peak_regions(in_file, out_file,
                                           max_merge_dist=0,
                                           sc_thr=0):
    """
    Extract peak regions from GraphProt .profile file.
    Store the peak regions (defined as regions with scores >= sc_thr)
    as to out_file in 6-column .bed.

    TODO:
    Add option for genomic coordinates input (+ - polarity support).
    Output genomic regions instead of sequence regions.

    >>> in_file = "test_data/test4.avg_profile"
    >>> out_file = "test_data/test4_out.peaks.bed"
    >>> exp_file = "test_data/test4_out_exp.peaks.bed"
    >>> exp2_file = "test_data/test4_out_exp2.peaks.bed"
    >>> empty_file = "test_data/empty_file"
    >>> graphprot_profile_extract_peak_regions(in_file, out_file)
    >>> diff_two_files_identical(out_file, exp_file)
    True
    >>> graphprot_profile_extract_peak_regions(in_file, out_file, sc_thr=10)
    >>> diff_two_files_identical(out_file, empty_file)
    True
    >>> graphprot_profile_extract_peak_regions(in_file, out_file, max_merge_dist=2)
    >>> diff_two_files_identical(out_file, exp2_file)
    True

    """
    # Check files.
    assert in_file != out_file, "input file == output file (\"%s\" == \"%s\")" %(in_file, out_file)
    OUTPEAKS = open(out_file, "w")
    # Old site ID.
    old_id = ""
    # Current site ID.
    cur_id = ""
    # Scores list.
    scores_list = []
    site_starts_dic = {}
    with open(in_file) as f:
        for line in f:
            cols = line.strip().split("\t")
            cur_id = cols[0]
            pos = int(cols[1]) # 0-based.
            score = float(cols[2])
            # Store first position of site.
            if cur_id not in site_starts_dic:
                # If first position != zero, we assume positions are 1-based.
                if pos != 0:
                    # Make index 0-based.
                    site_starts_dic[cur_id] = pos - 1
                else:
                    site_starts_dic[cur_id] = pos
            # Case: new site (new column 1 ID).
            if cur_id != old_id:
                # Process old id scores.
                if scores_list:
                    # Extract peaks from region.
                    peak_list = list_extract_peaks(scores_list,
                                                   max_merge_dist=max_merge_dist,
                                                   coords="bed",
                                                   sc_thr=sc_thr)
                    start_pos = site_starts_dic[old_id]
                    # Print out peaks in .bed format.
                    for l in peak_list:
                        peak_s = start_pos + l[0]
                        peak_e = start_pos + l[1]
                        site_id = "%s,%i" %(old_id, l[2])
                        OUTPEAKS.write("%s\t%i\t%i\t%s\t%f\t+\n" %(old_id, peak_s, peak_e, site_id, l[3]))
                    # Reset list.
                    scores_list = []
                old_id = cur_id
                scores_list.append(score)
            else:
                # Add to scores_list.
                scores_list.append(score)
    f.close()
    # Process last block.
    if scores_list:
        # Extract peaks from region.
        peak_list = list_extract_peaks(scores_list,
                                       max_merge_dist=max_merge_dist,
                                       coords="bed",
                                       sc_thr=sc_thr)
        start_pos = site_starts_dic[old_id]
        # Print out peaks in .bed format.
        for l in peak_list:
            peak_s = start_pos + l[0]
            peak_e = start_pos + l[1]
            site_id = "%s,%i" %(old_id, l[2]) # best score also 1-based.
            OUTPEAKS.write("%s\t%i\t%i\t%s\t%f\t+\n" %(old_id, peak_s, peak_e, site_id, l[3]))
    OUTPEAKS.close()


################################################################################

def bed_peaks_to_genomic_peaks(peak_file, genomic_peak_file, genomic_sites_bed):
    """
    Given a .bed file of sequence peak regions (possible coordinates from
    0 to length of s), convert peak coordinates to genomic coordinates.
    Do this by taking genomic regions of sequences as input.

    >>> test_in = "test_data/test.peaks.bed"
    >>> test_exp = "test_data/test_exp.peaks.bed"
    >>> test_out = "test_data/test_out.peaks.bed"
    >>> gen_in = "test_data/test.peaks_genomic.bed"
    >>> bed_peaks_to_genomic_peaks(test_in, test_out, gen_in)
    >>> diff_two_files_identical(test_out, test_exp)
    True

    """
    # Read in genomic region info.
    id2row_dic = {}

    with open(genomic_sites_bed) as f:
        for line in f:
            row = line.strip()
            cols = line.strip().split("\t")
            site_id = cols[3]
            assert site_id not in id2row_dic, "column 4 IDs not unique in given .bed file \"%s\"" %(genomic_sites_bed)
            id2row_dic[site_id] = row
    f.close()

    # Read in peaks file and convert coordinates.
    OUTPEAKS = open(genomic_peak_file, "w")
    with open(peak_file) as f:
        for line in f:
            cols = line.strip().split("\t")
            site_id = cols[0]
            site_s = int(cols[1])
            site_e = int(cols[2])
            site_id2 = cols[3]
            site_sc = float(cols[4])
            assert re.search(".+,.+", site_id2), "regular expression failed for ID \"%s\"" %(site_id2)
            m = re.search(".+,(\d+)", site_id2)
            sc_pos = int(m.group(1)) # 1-based.
            assert site_id in id2row_dic, "site ID \"%s\" not found in genomic sites dictionary" %(site_id)
            row = id2row_dic[site_id]
            rowl = row.split("\t")
            gen_chr = rowl[0]
            gen_s = int(rowl[1])
            gen_e = int(rowl[2])
            gen_pol = rowl[5]
            new_s = site_s + gen_s
            new_e = site_e + gen_s
            new_sc_pos = sc_pos + gen_s
            if gen_pol == "-":
                new_s = gen_e - site_e
                new_e = gen_e - site_s
                new_sc_pos = gen_e - sc_pos + 1 # keep 1-based.
            new_row = "%s\t%i\t%i\t%s,%i\t%f\t%s" %(gen_chr, new_s, new_e, site_id, new_sc_pos, site_sc, gen_pol)
            OUTPEAKS.write("%s\n" %(new_row))
    OUTPEAKS.close()


################################################################################

def graphprot_profile_calculate_avg_profile(in_file, out_file,
                                            ap_extlr=5,
                                            seq_ids_list=False,
                                            method=1):
    """
    Given a GraphProt .profile file, calculate average profiles and output
    average profile file.
    Average profile means that the position-wise scores will get smoothed
    out by calculating for each position a new score, taking a sequence
    window -ap_extlr to +ap_extlr relative to the position
    and calculate the mean score over this window. The mean score then
    becomes the new average profile score at this position.
    Two different implementations of the task are given:
    method=1 (new python implementation, slower + more memory but easy to read)
    method=2 (old perl implementation, faster and less memory but more code)

    >>> in_file = "test_data/test2.profile"
    >>> out_file1 = "test_data/test2_1.avg_profile"
    >>> out_file2 = "test_data/test2_2.avg_profile"
    >>> out_file4 = "test_data/test2_3.avg_profile"
    >>> graphprot_profile_calculate_avg_profile(in_file, out_file1, ap_extlr=2, method=1)
    >>> graphprot_profile_calculate_avg_profile(in_file, out_file2, ap_extlr=2, method=2)
    >>> diff_two_files_identical(out_file1, out_file2)
    True
    >>> test_list = ["s1", "s2", "s3", "s4"]
    >>> out_file3_exp = "test_data/test3_added_ids_exp.avg_profile"
    >>> out_file3 = "test_data/test3_added_ids_out.avg_profile"
    >>> graphprot_profile_calculate_avg_profile(in_file, out_file3, ap_extlr=2, method=1, seq_ids_list=test_list)
    >>> diff_two_files_identical(out_file3_exp, out_file3)
    True

    """
    # Check input files.
    assert in_file != out_file, "input file == output file (\"%s\" == \"%s\")" %(in_file, out_file)
    if method == 1:
        # Dictionary of lists, with list of scores (value) for each site (key).
        lists_dic = {}
        site_starts_dic = {}
        with open(in_file) as f:
            for line in f:
                cols = line.strip().split("\t")
                site_id = int(cols[0])
                pos = int(cols[1]) # 0-based.
                score = float(cols[2])
                # Store first position of site.
                if site_id not in site_starts_dic:
                    site_starts_dic[site_id] = pos
                if site_id in lists_dic:
                    lists_dic[site_id].append(score)
                else:
                    lists_dic[site_id] = []
                    lists_dic[site_id].append(score)
        f.close()
        # Check number of IDs (# FASTA sequence IDs has to be same as # site IDs).
        if seq_ids_list:
            c_seq_ids = len(seq_ids_list)
            c_site_ids = len(site_starts_dic)
            assert c_seq_ids == c_site_ids, "# sequence IDs != # site IDs (%i != %i)" %(c_seq_ids, c_site_ids)
        OUTPROF = open(out_file, "w")
        # For each site, calculate average profile scores list.
        max_list = []
        for site_id in lists_dic:
            # Convert profile score list to average profile scores list.
            aps_list = list_moving_window_average_values(lists_dic[site_id],
                                                         win_extlr=ap_extlr)
            start_pos = site_starts_dic[site_id]
            # Get original FASTA sequence ID.
            if seq_ids_list:
                site_id = seq_ids_list[site_id]
            for i, sc in enumerate(aps_list):
                pos = i + start_pos + 1 # make 1-based.
                OUTPROF.write("%s\t%i\t%f\n" %(site_id, pos, sc))
        OUTPROF.close()
    elif method == 2:
        OUTPROF = open(out_file, "w")
        # Old site ID.
        old_id = ""
        # Current site ID.
        cur_id = ""
        # Scores list.
        scores_list = []
        site_starts_dic = {}
        with open(in_file) as f:
            for line in f:
                cols = line.strip().split("\t")
                cur_id = int(cols[0])
                pos = int(cols[1]) # 0-based.
                score = float(cols[2])
                # Store first position of site.
                if cur_id not in site_starts_dic:
                    site_starts_dic[cur_id] = pos
                # Case: new site (new column 1 ID).
                if cur_id != old_id:
                    # Process old id scores.
                    if scores_list:
                        aps_list = list_moving_window_average_values(scores_list,
                                                                     win_extlr=ap_extlr)
                        start_pos = site_starts_dic[old_id]
                        seq_id = old_id
                        # Get original FASTA sequence ID.
                        if seq_ids_list:
                            seq_id = seq_ids_list[old_id]
                        for i, sc in enumerate(aps_list):
                            pos = i + start_pos + 1 # make 1-based.
                            OUTPROF.write("%s\t%i\t%f\n" %(seq_id, pos, sc))
                        # Reset list.
                        scores_list = []
                    old_id = cur_id
                    scores_list.append(score)
                else:
                    # Add to scores_list.
                    scores_list.append(score)
        f.close()
        # Process last block.
        if scores_list:
            aps_list = list_moving_window_average_values(scores_list,
                                                         win_extlr=ap_extlr)
            start_pos = site_starts_dic[old_id]
            seq_id = old_id
            # Get original FASTA sequence ID.
            if seq_ids_list:
                seq_id = seq_ids_list[old_id]
            for i, sc in enumerate(aps_list):
                pos = i + start_pos + 1 # make 1-based.
                OUTPROF.write("%s\t%i\t%f\n" %(seq_id, pos, sc))
        OUTPROF.close()


################################################################################

def graphprot_profile_get_top_scores_median(profile_file,
                                            profile_type="profile",
                                            avg_profile_extlr=5):

    """
    Given a GraphProt .profile file, extract for each site (identified by
    column 1 ID) the top (= highest) score. Then return the median of these
    top scores.

    profile_type can be either "profile" or "avg_profile".
    "avg_profile means that the position-wise scores will first get smoothed
    out by calculating for each position a new score through taking a
    sequence window -avg_profile_extlr to +avg_profile_extlr of the position
    and calculate the mean score over this window and assign it to the position.
    After that, the maximum score of each site is chosen, and the median over
    all maximum scores is returned.
    "profile" leaves the position-wise scores as they are, directly extracting
    the maximum for each site and then reporting the median.

    >>> test_file = "test_data/test.profile"
    >>> graphprot_profile_get_top_scores_median(test_file)
    3.2

    """
    # Dictionary of lists, with list of scores (value) for each site (key).
    lists_dic = {}
    with open(profile_file) as f:
        for line in f:
            cols = line.strip().split("\t")
            seq_id = cols[0]
            score = float(cols[2])
            if seq_id in lists_dic:
                lists_dic[seq_id].append(score)
            else:
                lists_dic[seq_id] = []
                lists_dic[seq_id].append(score)
    f.close()
    # For each site, extract maximum and store in new list.
    max_list = []
    for seq_id in lists_dic:
        if profile_type == "profile":
            max_sc = max(lists_dic[seq_id])
            max_list.append(max_sc)
        elif profile_type == "avg_profile":
            # Convert profile score list to average profile scores list.
            aps_list = list_moving_window_average_values(lists_dic[seq_id],
                                                         win_extlr=avg_profile_extlr)
            max_sc = max(aps_list)
            max_list.append(max_sc)
        else:
            assert 0, "invalid profile_type argument given: \"%s\"" %(profile_type)
    # Return the median.
    return statistics.median(max_list)


################################################################################

def list_extract_peaks(in_list,
                       max_merge_dist=0,
                       coords="list",
                       sc_thr=0):
    """
    Extract peak regions from list.
    Peak region is defined as region >= score threshold.

    coords=bed  :  peak start 0-based, peak end 1-based.
    coords=list :  peak start 0-based, peak end 0-based.

    >>> test_list = [-1, 0, 2, 4.5, 1, -1, 5, 6.5]
    >>> list_extract_peaks(test_list)
    [[1, 4, 3, 4.5], [6, 7, 7, 6.5]]
    >>> list_extract_peaks(test_list, sc_thr=2)
    [[2, 3, 3, 4.5], [6, 7, 7, 6.5]]
    >>> list_extract_peaks(test_list, sc_thr=2, coords="bed")
    [[2, 4, 4, 4.5], [6, 8, 8, 6.5]]
    >>> list_extract_peaks(test_list, sc_thr=10)
    []
    >>> test_list = [2, -1, 3, -1, 4, -1, -1, 6, 9]
    >>> list_extract_peaks(test_list, max_merge_dist=2)
    [[0, 4, 4, 4], [7, 8, 8, 9]]
    >>> list_extract_peaks(test_list, max_merge_dist=3)
    [[0, 8, 8, 9]]

    """
    # Check.
    assert len(in_list), "Given list is empty"
    # Peak regions list.
    peak_list = []
    # Help me.
    inside = False
    pr_s = 0
    pr_e = 0
    pr_top_pos = 0
    pr_top_sc = -100000
    for i, sc in enumerate(in_list):
        # Part of peak region?
        if sc >= sc_thr:
            # At peak start.
            if not inside:
                pr_s = i
                pr_e = i
                inside = True
            else:
                # Inside peak region.
                pr_e = i
            # Store top position.
            if sc > pr_top_sc:
                pr_top_sc = sc
                pr_top_pos = i
        else:
            # Before was peak region?
            if inside:
                # Store peak region.
                #peak_infos = "%i,%i,%i,%f" %(pr_s, pr_e, pr_top_pos, pr_top_sc)
                peak_infos = [pr_s, pr_e, pr_top_pos, pr_top_sc]
                peak_list.append(peak_infos)
                inside = False
                pr_top_pos = 0
                pr_top_sc = -100000
    # If peak at the end, also report.
    if inside:
        # Store peak region.
        peak_infos = [pr_s, pr_e, pr_top_pos, pr_top_sc]
        peak_list.append(peak_infos)
    # Merge peaks.
    if max_merge_dist and len(peak_list) > 1:
        iterate = True
        while iterate:
            merged_peak_list = []
            added_peaks_dic = {}
            peaks_merged = False
            for i, l in enumerate(peak_list):
                if i in added_peaks_dic:
                    continue
                j = i + 1
                # Last element.
                if j == len(peak_list):
                    if i not in added_peaks_dic:
                        merged_peak_list.append(peak_list[i])
                    break
                # Compare two elements.
                new_peak = []
                if (peak_list[j][0] - peak_list[i][1]) <= max_merge_dist:
                    peaks_merged = True
                    new_top_pos = peak_list[i][2]
                    new_top_sc = peak_list[i][3]
                    if peak_list[i][3] < peak_list[j][3]:
                        new_top_pos = peak_list[j][2]
                        new_top_sc = peak_list[j][3]
                    new_peak = [peak_list[i][0], peak_list[j][1], new_top_pos, new_top_sc]
                # If two peaks were merged.
                if new_peak:
                    merged_peak_list.append(new_peak)
                    added_peaks_dic[i] = 1
                    added_peaks_dic[j] = 1
                else:
                    merged_peak_list.append(peak_list[i])
                    added_peaks_dic[i] = 1
            if not peaks_merged:
                iterate = False
            peak_list = merged_peak_list
            peaks_merged = False
    # If peak coordinates should be in .bed format, make peak ends 1-based.
    if coords == "bed":
        for i in range(len(peak_list)):
            peak_list[i][1] += 1
            peak_list[i][2] += 1 # 1-base best score position too.
    return peak_list


################################################################################

def list_moving_window_average_values(in_list,
                                      win_extlr=5,
                                      method=1):
    """
    Take a list of numeric values, and calculate for each position a new value,
    by taking the mean value of the window of positions -win_extlr and
    +win_extlr. If full extension is not possible (at list ends), it just
    takes what it gets.
    Two implementations of the task are given, chose by method=1 or method=2.

    >>> test_list = [2, 3, 5, 8, 4, 3, 7, 1]
    >>> list_moving_window_average_values(test_list, win_extlr=2, method=1)
    [3.3333333333333335, 4.5, 4.4, 4.6, 5.4, 4.6, 3.75, 3.6666666666666665]
    >>> list_moving_window_average_values(test_list, win_extlr=2, method=2)
    [3.3333333333333335, 4.5, 4.4, 4.6, 5.4, 4.6, 3.75, 3.6666666666666665]

    """
    l_list = len(in_list)
    assert l_list, "Given list is empty"
    new_list = [0] * l_list
    if win_extlr == 0:
        return l_list
    if method == 1:
        for i in range(l_list):
            s = i - win_extlr
            e = i + win_extlr + 1
            if s < 0:
                s = 0
            if e > l_list:
                e = l_list
            # Extract portion and assign value to new list.
            new_list[i] = statistics.mean(in_list[s:e])
    elif method == 2:
        for i in range(l_list):
            s = i - win_extlr
            e = i + win_extlr + 1
            if s < 0:
                s = 0
            if e > l_list:
                e = l_list
            l = e-s
            sc_sum = 0
            for j in range(l):
                sc_sum += in_list[s+j]
            new_list[i] = sc_sum / l
    else:
        assert 0, "invalid method ID given (%i)" %(method)
    return new_list


################################################################################

def seqs_dic_count_uc_nts(seqs_dic):
    """
    Count number of uppercase nucleotides in sequences stored in sequence
    dictionary.

    >>> seqs_dic = {'seq1': "acgtACGTacgt", 'seq2': 'acgtACacgt'}
    >>> seqs_dic_count_uc_nts(seqs_dic)
    6
    >>> seqs_dic = {'seq1': "acgtacgt", 'seq2': 'acgtacgt'}
    >>> seqs_dic_count_uc_nts(seqs_dic)
    0

    """
    assert seqs_dic, "Given sequence dictionary empty"
    c_uc = 0
    for seq_id in seqs_dic:
        c_uc += len(re.findall(r'[A-Z]', seqs_dic[seq_id]))
    return c_uc


################################################################################

def seqs_dic_count_lc_nts(seqs_dic):
    """
    Count number of lowercase nucleotides in sequences stored in sequence
    dictionary.

    >>> seqs_dic = {'seq1': "gtACGTac", 'seq2': 'cgtACacg'}
    >>> seqs_dic_count_lc_nts(seqs_dic)
    10
    >>> seqs_dic = {'seq1': "ACGT", 'seq2': 'ACGTAC'}
    >>> seqs_dic_count_lc_nts(seqs_dic)
    0

    """
    assert seqs_dic, "Given sequence dictionary empty"
    c_uc = 0
    for seq_id in seqs_dic:
        c_uc += len(re.findall(r'[a-z]', seqs_dic[seq_id]))
    return c_uc


################################################################################

def graphprot_predictions_get_median(predictions_file):
    """
    Given a GraphProt .predictions file, read in site scores and return
    the median value.

    >>> test_file = "test_data/test.predictions"
    >>> graphprot_predictions_get_median(test_file)
    0.571673

    """
    # Site scores list.
    sc_list = []
    with open(predictions_file) as f:
        for line in f:
            cols = line.strip().split("\t")
            score = float(cols[2])
            sc_list.append(score)
    f.close()
    # Return the median.
    return statistics.median(sc_list)


################################################################################

def graphprot_filter_predictions_file(in_file, out_file,
                                      sc_thr=0):
    """
    Filter GraphProt .predictions file by given score thr_sc.
    """
    assert in_file != out_file, "input file == output file (\"%s\" == \"%s\")" %(in_file, out_file)
    OUTPRED = open(out_file, "w")
    with open(in_file) as f:
        for line in f:
            row = line.strip()
            cols = line.strip().split("\t")
            score = float(cols[2])
            if score < sc_thr:
                continue
            OUTPRED.write("%s\n" %(row))
    f.close()
    OUTPRED.close()


################################################################################

def graphprot_get_param_dic(params_file):
    """
    Read in GraphProt .params file and store in dictionary.
    key = parameter
    value = parameter value

    >>> params_file = "test_data/test.params"
    >>> graphprot_get_param_dic(params_file)
    {'epochs': '40', 'lambda': '0.001', 'R': '4', 'D': '6', 'bitsize': '14', 'model_type': 'sequence', 'pos_train_ws_pred_median': '1.033690', 'pos_train_profile_median': '8.680340', 'pos_train_avg_profile_median_1': '4.027981', 'pos_train_avg_profile_median_2': '3.027981'}

    """
    param_dic = {}
    with open(params_file) as f:
        for line in f:
            cols = line.strip().split(" ")
            param = cols[0]
            setting = cols[1]
            if re.search(".+:", param):
                m = re.search("(.+):", line)
                par = m.group(1)
                param_dic[par] = setting
    f.close()
    return param_dic


################################################################################

def graphprot_get_param_string(params_file):
    """
    Get parameter string from GraphProt .params file.

    >>> test_params = "test_data/test.params"
    >>> graphprot_get_param_string(test_params)
    '-epochs 40 -lambda 0.001 -R 4 -D 6 -bitsize 14 -onlyseq '

    """
    param_string = ""
    with open(params_file) as f:
        for line in f:
            cols = line.strip().split(" ")
            param = cols[0]
            setting = cols[1]
            if re.search(".+:", param):
                m = re.search("(.+):", line)
                par = m.group(1)
                if re.search("pos_train.+", line):
                    continue
                if par == "model_type":
                    if setting == "sequence":
                        param_string += "-onlyseq "
                else:
                    param_string += "-%s %s " %(par, setting)
            else:
                assert 0, "pattern matching failed for string \"%s\"" %(param)
    return param_string


################################################################################

def echo_add_to_file(echo_string, out_file):
    """
    Add a string to file, using echo command.

    """
    check_cmd = 'echo "%s" >> %s' % (echo_string, out_file)
    output = subprocess.getoutput(check_cmd)
    error = False
    if output:
        error = True
    assert error == False, "echo is complaining:\n%s\n%s" %(check_cmd, output)


################################################################################

def bed_merge_file_select_top_ids(in_merged_bed, id2sc_dic,
                                  rev_filter=False):
    """
    Given a merged .bed file (using mergeBed or bed_merge_file() ),
    select the top scoring region IDs from the merged file, where in case
    of overlaps the best scoring region ID is picked.
    rev_filter=True leads to lower scores regarded as better, e.g. in case
    of p-values.

    >>> test_merged = "test_data/test.sorted.merged.bed"
    >>> id2sc_dic = {'r1': 1, 'r2': 2, 'r4': 4, 'r7' : 3}
    >>> bed_merge_file_select_top_ids(test_merged, id2sc_dic)
    {'r4': 4, 'r7': 3}
    >>> bed_merge_file_select_top_ids(test_merged, id2sc_dic, rev_filter=True)
    {'r1': 1, 'r7': 3}

    """
    ids2keep_dic = {}
    assert id2sc_dic, "given ID to score dictionary seems to be empty"
    with open(in_merged_bed) as f:
        for line in f:
            cols = line.strip().split("\t")
            ids = cols[3].split(";")
            best_id = "-"
            best_sc = 0
            if rev_filter:
                best_sc = 1000000
            for site_id in ids:
                assert site_id in id2sc_dic, "site ID \"%s\" not found in given site ID to score dictionary" % (site_id)
                site_sc = id2sc_dic[site_id]
                if rev_filter:
                    if site_sc < best_sc:
                        best_sc = site_sc
                        best_id = site_id
                else:
                    if site_sc > best_sc:
                        best_sc = site_sc
                        best_id = site_id
            ids2keep_dic[best_id] = best_sc
    f.closed
    assert ids2keep_dic, "No IDs read into dictionary (input file \"%s\" empty or malformatted?)" % (in_merged_bed)
    return ids2keep_dic


################################################################################

def intersect_bed_files(a_file, b_file, params, out_file,
                        sorted_out=False):
    """
    Intersect two .bed files, using intersectBed.

    """

    check_cmd = "intersectBed -a " + a_file + " -b " + b_file + " " + params + " > " + out_file
    if sorted_out:
        check_cmd = "intersectBed -a " + a_file + " -b " + b_file + " " + params + " | " + "sort -k1,1 -k2,2n > " + out_file
    output = subprocess.getoutput(check_cmd)
    error = False
    if output:
        error = True
    assert error == False, "intersectBed has problems with your input:\n%s\n%s" %(check_cmd, output)


################################################################################

def read_ids_into_dic(ids_file,
                      ids_dic=False):
    """
    Read in IDs list file, where each line stores one ID.

    >>> test_ids_file = "test_data/test.ids"
    >>> read_ids_into_dic(test_ids_file)
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

def count_fasta_headers(fasta_file):
    """
    Count number of FASTA headers in fasta_file using grep.

    >>> test_file = "test_data/test.fa"
    >>> count_fasta_headers(test_file)
    2
    >>> test_file = "test_data/empty_file"
    >>> count_fasta_headers(test_file)
    0

    """
    check_cmd = 'grep -c ">" ' + fasta_file
    output = subprocess.getoutput(check_cmd)
    row_count = int(output.strip())
    return row_count


################################################################################

def random_order_dic_keys_into_list(in_dic):
    """
    Read in dictionary keys, and return random order list of IDs.

    """
    id_list = []
    for key in in_dic:
        id_list.append(key)
    random.shuffle(id_list)
    return id_list


################################################################################

def split_fasta_into_test_train_files(in_fasta, test_out_fa, train_out_fa,
                                      test_size=500):
    """
    Split in_fasta .fa file into two files (e.g. test, train).

    """
    # Check files.
    assert in_fasta != test_out_fa, "input file == output file (\"%s\" == \"%s\")" %(in_fasta, test_out_fa)
    assert in_fasta != train_out_fa, "input file == output file (\"%s\" == \"%s\")" %(in_fasta, train_out_fa)
    # Read in in_fasta.
    seqs_dic = read_fasta_into_dic(in_fasta)
    # Shuffle IDs.
    rand_ids_list = random_order_dic_keys_into_list(seqs_dic)
    c_out = 0
    TESTOUT = open(test_out_fa, "w")
    TRAINOUT = open(train_out_fa, "w")
    for seq_id in rand_ids_list:
        seq = seqs_dic[seq_id]
        if (c_out >= test_size):
            TRAINOUT.write(">%s\n%s\n" % (seq_id, seq))
        else:
            TESTOUT.write(">%s\n%s\n" % (seq_id, seq))
        c_out += 1
    TESTOUT.close()
    TRAINOUT.close()


################################################################################

def read_fasta_into_dic(fasta_file,
                        seqs_dic=False,
                        ids_dic=False,
                        read_dna=False,
                        reject_lc=False,
                        convert_to_uc=False,
                        skip_n_seqs=True):
    """
    Read in FASTA sequences, convert to RNA, store in dictionary
    and return dictionary.

    >>> test_fasta = "test_data/test.fa"
    >>> read_fasta_into_dic(test_fasta)
    {'seq1': 'acguACGUacgu', 'seq2': 'ugcaUGCAugcaACGUacgu'}
    >>> test_fasta = "test_data/test2.fa"
    >>> read_fasta_into_dic(test_fasta)
    {}
    >>> test_fasta = "test_data/test.ensembl.fa"
    >>> read_fasta_into_dic(test_fasta, read_dna=True)
    {'ENST00000415118': 'GAAATAGT', 'ENST00000448914': 'ACTGGGGGATACGAAAA'}

    """
    if not seqs_dic:
        seqs_dic = {}
    seq_id = ""
    seq = ""

    # Go through FASTA file, extract sequences.
    if re.search(".+\.gz$", fasta_file):
        f = gzip.open(fasta_file, 'rt')
    else:
        f = open(fasta_file, "r")
    for line in f:
        if re.search(">.+", line):
            m = re.search(">(.+)", line)
            seq_id = m.group(1)
            # If there is a ".", take only first part of header.
            # This assumes ENSEMBL header format ">ENST00000631435.1 cdna ..."
            if re.search(".+\..+", seq_id):
                m = re.search("(.+?)\..+", seq_id)
                seq_id = m.group(1)
            assert seq_id not in seqs_dic, "non-unique FASTA header \"%s\" in \"%s\"" % (seq_id, fasta_file)
            if ids_dic:
                if seq_id in ids_dic:
                    seqs_dic[seq_id] = ""
            else:
                seqs_dic[seq_id] = ""
        elif re.search("[ACGTUN]+", line, re.I):
            if seq_id in seqs_dic:
                m = re.search("([ACGTUN]+)", line, re.I)
                seq = m.group(1)
                if reject_lc:
                    assert not re.search("[a-z]", seq), "lowercase characters detected in sequence \"%i\" (reject_lc=True)" %(seq_id)
                if convert_to_uc:
                    seq = seq.upper()
                # If sequences with N nucleotides should be skipped.
                if skip_n_seqs:
                    if "n" in m.group(1) or "N" in m.group(1):
                        print ("WARNING: \"%s\" contains N nucleotides. Discarding sequence ... " % (seq_id))
                        del seqs_dic[seq_id]
                        continue
                # Convert to RNA, concatenate sequence.
                if read_dna:
                    seqs_dic[seq_id] += m.group(1).replace("U","T").replace("u","t")
                else:
                    seqs_dic[seq_id] += m.group(1).replace("T","U").replace("t","u")
    f.close()
    return seqs_dic


################################################################################

def generate_random_fn(file_ending):
    """
    Generate a random file name for temporary files.

    """
    random_id = uuid.uuid1()
    random_fn = str(random_id) + ".tmp." . file_ending
    return random_fn


################################################################################

def make_file_copy(in_file, out_file):
    """
    Make a file copy by copying in_file to out_file.

    """
    check_cmd = "cat " + in_file + " > " + out_file
    assert in_file != out_file, "cat does not like to cat file into same file (%s)" %(check_cmd)
    output = subprocess.getoutput(check_cmd)
    error = False
    if output:
        error = True
    assert error == False, "cat did not like your input (in_file: %s, out_file: %s):\n%s" %(in_file, out_file, output)


################################################################################

def diff_two_files_identical(file1, file2):
    """
    Check whether two files are identical. Return true if diff reports no
    differences.

    >>> file1 = "test_data/file1"
    >>> file2 = "test_data/file2"
    >>> diff_two_files_identical(file1, file2)
    True
    >>> file1 = "test_data/test1.bed"
    >>> diff_two_files_identical(file1, file2)
    False

    """
    same = True
    check_cmd = "diff " + file1 + " " + file2
    output = subprocess.getoutput(check_cmd)
    if output:
        same = False
    return same


################################################################################

def bed_get_incomplete_overlaps(site_bed, region_bed, out_bed):
    """
    Overlap sites .bed with regions .bed, and return .bed file containing
    only incomplete (not full length matching) sites .bed.

    """

    # intersect with -f 1
    intersect_params = "-s -f 1 -v"
    intersect_bed_files(site_bed, region_bed, intersect_params, out_bed)


################################################################################

def bed_read_ids_into_dic(in_bed):
    """
    Read in .bed file IDs (column 4) into dictionary.

    >>> test_bed = "test_data/test4.bed"
    >>> bed_read_ids_into_dic(test_bed)
    {'CLIP1': 1, 'CLIP2': 1}

    """
    ids_dic = {}
    with open(in_bed) as f:
        for line in f:
            cols = line.strip().split("\t")
            site_id = cols[3]
            ids_dic[site_id] = 1
    f.closed
    assert ids_dic, "No IDs read into dictionary (input file \"%s\" empty or malformatted?)" % (in_bed)
    return ids_dic


################################################################################

def bed_count_region_overlaps(site_bed, region_bed, overlap_out,
                              only_incomplete=False,
                              intersectBed_f=False):
    """
    Overlap sites with regions, and return number of sites overlapping for
    each region. Also return the list of site IDs for each region.
    only_incomplete : If True, count only incompletely overlapping sites.

    >>> site_bed = "test_data/overlap_sites_in.bed"
    >>> region_bed = "test_data/overlap_regions_in.bed"
    >>> overlap_out = "test_data/intersectbed.overlap.tmp.bed"
    >>> bed_count_region_overlaps(site_bed, region_bed, overlap_out)
    ({'i1': 2, 'i2': 2}, {'i1': ['CLIP1', 'CLIP2'], 'i2': ['CLIP3', 'CLIP4']})

    """
    # Intersect.
    intersect_params = "-s -wb" # Any overlaps okay.
    if intersectBed_f:
        intersect_params = "-s -wb -f %f" %(intersectBed_f)
    intersect_bed_files(site_bed, region_bed, intersect_params, overlap_out)
    # If only_incomplete, first read in site lengths.
    site_len_dic = {}
    if only_incomplete:
        site_len_dic = bed_get_region_lengths(site_bed)
    c_all = 0
    reg_olc_dic = {}
    reg_site_ids_dic = {}
    # Read in overlaps.
    with open(overlap_out) as f:
        for line in f:
            cols = line.strip().split("\t")
            site_match_s = int(cols[1])
            site_match_e = int(cols[2])
            site_id = cols[3]
            reg_id = cols[9]
            if only_incomplete:
                site_match_l = site_match_e - site_match_s
                if site_match_l == site_len_dic[site_id]:
                    continue
            c_all += 1
            if reg_id in reg_olc_dic:
                reg_olc_dic[reg_id] += 1
            else:
                reg_olc_dic[reg_id] = 1
            if reg_id in reg_site_ids_dic:
                reg_site_ids_dic[reg_id].append(site_id)
            else:
                reg_site_ids_dic[reg_id] = [site_id]
    f.close()
    assert c_all, "ERROR: no overlapping sites when overlapping \"%s\" with \"%s\"" %(site_bed, region_bed)
    return reg_olc_dic, reg_site_ids_dic


################################################################################

def gtf_extract_exon_bed(in_gtf, out_bed,
                         out_intron_bed=False,
                         tr_ids_dic=False):
    """
    Given a .gtf file with exon features, extract exon regions and store in
    .bed file. Optionally, a dictionary of transcript IDs can be provided,
    meaning that only exon regions from the given transcripts will be extracted.
    If out_intron_bed is set, an intronic regions .bed file will also be
    extracted, based on the exonic regions .bed information.

    Output .bed will look like this (note column 4 ID format with transcript
    ID followed by _e+exon_number):
    chr1	1000	2000	ENST001_e1	0	+
    chr1	3000	4000	ENST001_e2	0	+
    chr1	8000	9000	ENST002_e1	0	-
    chr1	6000	7000	ENST002_e2	0	-
    ...

    NOTE that function has been tested with .gtf files from Ensembl. .gtf files
    from different sources sometimes have a slightly different format, which
    could lead to incompatibilities / errors. See test files for format that
    works.

    Some tested Ensembl GTF files:
    Homo_sapiens.GRCh38.97.gtf.gz
    Mus_musculus.GRCm38.81.gtf.gz
    Mus_musculus.GRCm38.79.gtf.gz

    >>> in_gtf = "test_data/map_test_in.gtf"
    >>> exp_out_bed = "test_data/gtf_exon_out_exp.bed"
    >>> exp_out_intron_bed = "test_data/gtf_intron_out_exp.bed"
    >>> out_bed = "test_data/gtf_exon_out.bed"
    >>> out_intron_bed = "test_data/gtf_intron_out.bed"
    >>> gtf_extract_exon_bed(in_gtf, out_bed, out_intron_bed=out_intron_bed)
    >>> diff_two_files_identical(out_bed, exp_out_bed)
    True
    >>> diff_two_files_identical(out_intron_bed, exp_out_intron_bed)
    True

    """

    # Output genomic exon regions.
    OUTBED = open(out_bed, "w")

    # Read in exon features from GTF file.
    c_gtf_ex_feat = 0
    # Start end coordinates of exons.
    exon_e_dic = {}
    exon_s_dic = {}
    # Transcript stats.
    tr2pol_dic = {}
    tr2chr_dic = {}
    # dic for sanity checking exon number order.
    tr2exon_nr_dic = {}

    # Open GTF either as .gz or as text file.
    if re.search(".+\.gz$", in_gtf):
        f = gzip.open(in_gtf, 'rt')
    else:
        f = open(in_gtf, "r")
    for line in f:
        # Skip header.
        if re.search("^#", line):
            continue
        cols = line.strip().split("\t")
        chr_id = cols[0]
        feature = cols[2]
        feat_s = int(cols[3])
        feat_e = int(cols[4])
        feat_pol = cols[6]
        infos = cols[8]
        if not feature == "exon":
            continue

        # Restrict to standard chromosomes.
        if re.search("^chr", chr_id):
            if not re.search("^chr[\dMXY]", chr_id):
                continue
        else:
            # Convert to "chr" IDs.
            if not re.search("^[\dMXY]", chr_id):
                continue
            if chr_id == "MT":
                chr_id == "M"
            chr_id = "chr" + chr_id

        # Make start coordinate 0-base (BED standard).
        feat_s = feat_s - 1

        # Extract transcript ID.
        m = re.search('transcript_id "(.+?)"', infos)
        assert m, "transcript_id entry missing in GTF file \"%s\", line \"%s\"" %(in_gtf, line)
        transcript_id = m.group(1)
        # Extract exon number.
        m = re.search('exon_number "(\d+?)"', infos)
        assert m, "exon_number entry missing in GTF file \"%s\", line \"%s\"" %(in_gtf, line)
        exon_nr = int(m.group(1))

        # Check if transcript ID is in transcript dic.
        if tr_ids_dic:
            if not transcript_id in tr_ids_dic:
                continue

        # Store transcript stats.
        tr2pol_dic[transcript_id] = feat_pol
        tr2chr_dic[transcript_id] = chr_id

        # Check whether exon numbers are incrementing for each transcript ID.
        if not transcript_id in tr2exon_nr_dic:
            tr2exon_nr_dic[transcript_id] = exon_nr
        else:
            assert tr2exon_nr_dic[transcript_id] < exon_nr, "transcript ID \"%s\" without increasing exon number order in GTF file \"%s\"" %(transcript_id, in_gtf)
            tr2exon_nr_dic[transcript_id] = exon_nr

        # Count exon entry.
        c_gtf_ex_feat += 1

        # Construct exon ID.
        exon_id = transcript_id + "_e" + str(exon_nr)
        # Store infos.
        exon_s_dic[exon_id] = feat_s
        exon_e_dic[exon_id] = feat_e

        # Output genomic exon region.
        OUTBED.write("%s\t%i\t%i\t%s\t0\t%s\n" % (chr_id,feat_s,feat_e,exon_id,feat_pol))

    OUTBED.close()
    f.close()

    # Check for read-in features.
    assert c_gtf_ex_feat, "no exon features read in from \"%s\"" %(in_gtf)

    # Output intron .bed.
    if out_intron_bed:
        tr2intron_nr_dic = {}
        OUTBED = open(out_intron_bed, "w")
        for tr_id in tr2pol_dic:
            tr_pol = tr2pol_dic[tr_id]
            chr_id = tr2chr_dic[tr_id]
            tr_c = tr2exon_nr_dic[tr_id]
            intron_c = 0
            tr2intron_nr_dic[tr_id] = 0
            # 1-exon transcripts, no introns.
            if tr_c == 1:
                continue
            ex_list = []
            for i in range(tr_c):
                ex_nr = i + 1
                ex_id = tr_id + "_e" + str(ex_nr)
                ex_list.append(ex_id)
            for i in range(len(ex_list)):
                ex1i = i
                ex2i = i + 1
                # At last exon, no more introns to add.
                if ex2i == len(ex_list):
                    break
                ex1id = ex_list[ex1i]
                ex2id = ex_list[ex2i]
                ex1s = exon_s_dic[ex1id]
                ex2s = exon_s_dic[ex2id]
                ex1e = exon_e_dic[ex1id]
                ex2e = exon_e_dic[ex2id]
                # Plus case.
                intron_s = ex1e
                intron_e = ex2s
                if tr_pol == "-":
                    intron_s = ex2e
                    intron_e = ex1s
                intron_id = tr_id + "_i" + str(ex2i)
                intron_c += 1
                OUTBED.write("%s\t%i\t%i\t%s\t0\t%s\n" % (chr_id,intron_s,intron_e,intron_id,tr_pol))
            tr2intron_nr_dic[tr_id] = intron_c
        OUTBED.close()
        # Sanity check exon + intron numbers.
        for tr_id in tr2exon_nr_dic:
            exon_nr = tr2exon_nr_dic[tr_id]
            intron_nr = tr2intron_nr_dic[tr_id]
            assert (exon_nr-1) == intron_nr, "ERROR: intron number != exon number - 1 for \"%s\" (%i != %i - 1)" %(tr_id, intron_nr, exon_nr)


################################################################################

def convert_genome_positions_to_transcriptome(in_bed, out_folder,
                                              in_gtf, tr_ids_dic,
                                              intersectBed_f=1,
                                              int_whole_nr=True,
                                              ignore_ids_dic=False):
    """
    Converts a BED file with genomic coordinates into a BED file with
    transcriptome coordinates. A GTF file with exon features needs to be
    supplied. A dictionary of transcript IDs defines to which transcripts
    the genomic regions will be mapped to. Note that input BED file column
    4 is used as region ID and should be unique.

    Output files are:
    genomic_exon_coordinates.bed
        Genomic exon coordinates extracted from .gtf file
    transcript_exon_coordinates.bed
        Transcript exon coordinates calculated from genomic ones
    hit_exon_overlap.bed
        Overlap between genomic input and exon regions
    transcript_matches_complete.bed
        Input regions that fully (completely) map to transcript regions
    transcript_matches_incomplete.bed
        Input regions that partly (incompletely) map to transcript regions
    transcript_matches_complete_unique.bed
        Unique + complete (full-length) matches
    transcript_matches_all_unique.bed
        All unique (complete+incomplete) matches
    hit_transcript_exons.bed
        Genomic coordinates of exons of transcripts with mapped regions
    hit_transcript_stats.out
        Various statistics for transcripts with hits
        e.g. gene ID, gene name, gene biotype, # unique complete hits,
        # unique all hits, # complete hits, # all hits

    NOTE that function has been tested with .gtf files from Ensembl. .gtf files
    from different sources sometimes have a slightly different format, which
    could lead to incompatibilities / errors. See test files for format that
    works.

    Some tested Ensembl GTF files:
    Homo_sapiens.GRCh38.97.gtf.gz
    Mus_musculus.GRCm38.81.gtf.gz
    Mus_musculus.GRCm38.79.gtf.gz

    Requirements:
    bedTools (tested with version 2.25.0)
    GTF file needs to have exons sorted (minus + plus strand exons, see test.gtf
    below as an example). Sorting should be the default (at least for tested
    Ensembl GTF files)

    >>> tr_ids_dic = {"ENST001" : 1, "ENST002" : 1}
    >>> in_bed = "test_data/map_test_in.bed"
    >>> in_gtf = "test_data/map_test_in.gtf"
    >>> comp_uniq_exp = "test_data/map_test_out_all_unique.bed"
    >>> comp_uniq_out = "test_data/map_out/transcript_hits_all_unique.bed"
    >>> tr_stats_exp = "test_data/map_test_out_transcript_stats.out"
    >>> tr_stats_out = "test_data/map_out/hit_transcript_stats.out"
    >>> out_folder = "test_data/map_out"
    >>> convert_genome_positions_to_transcriptome(in_bed, out_folder, in_gtf, tr_ids_dic, intersectBed_f=0.5)
    >>> diff_two_files_identical(comp_uniq_exp, comp_uniq_out)
    True
    >>> diff_two_files_identical(tr_stats_exp, tr_stats_out)
    True

    """
    # Check for bedtools.
    assert is_tool("bedtools"), "bedtools not in PATH"

    # Results output folder.
    if not os.path.exists(out_folder):
        os.makedirs(out_folder)
    # Output files.
    genome_exon_bed = out_folder + "/" + "genomic_exon_coordinates.bed"
    transcript_exon_bed = out_folder + "/" + "transcript_exon_coordinates.bed"
    overlap_out = out_folder + "/" + "hit_exon_overlap.bed"
    complete_transcript_hits_bed = out_folder + "/" + "transcript_hits_complete.bed"
    incomplete_transcript_hits_bed = out_folder + "/" + "transcript_hits_incomplete.bed"
    uniq_complete_out = out_folder + "/" + "transcript_hits_complete_unique.bed"
    uniq_all_out = out_folder + "/" + "transcript_hits_all_unique.bed"
    hit_tr_exons_bed = out_folder + "/" + "hit_transcript_exons.bed"
    hit_tr_stats_out = out_folder + "/" + "hit_transcript_stats.out"

    # Check for unique .bed IDs.
    assert bed_check_unique_ids(in_bed), "in_bed \"%s\" column 4 IDs not unique" % (in_bed)

    # Remove IDs to ignore from transcript IDs dictionary.
    if ignore_ids_dic:
        for seq_id in ignore_ids_dic:
            del tr_ids_dic[seq_id]

    # Output genomic exon regions.
    OUTBED = open(genome_exon_bed, "w")

    # Read in exon features from GTF file.
    tr2gene_id_dic = {}
    tr2gene_name_dic = {}
    tr2gene_biotype_dic = {}
    c_gtf_ex_feat = 0
    # dic for sanity checking exon number order.
    tr2exon_nr_dic = {}
    # dic of lists, storing exon lengths and IDs.
    tr_exon_len_dic = {}
    tr_exon_id_dic = {}
    exon_id_tr_dic = {}

    # Open GTF either as .gz or as text file.
    if re.search(".+\.gz$", in_gtf):
        f = gzip.open(in_gtf, 'rt')
    else:
        f = open(in_gtf, "r")
    for line in f:
        # Skip header.
        if re.search("^#", line):
            continue
        cols = line.strip().split("\t")
        chr_id = cols[0]
        feature = cols[2]
        feat_s = int(cols[3])
        feat_e = int(cols[4])
        feat_pol = cols[6]
        infos = cols[8]
        if not feature == "exon":
            continue

        # Restrict to standard chromosomes.
        if re.search("^chr", chr_id):
            if not re.search("^chr[\dMXY]", chr_id):
                continue
        else:
            # Convert to "chr" IDs.
            if not re.search("^[\dMXY]", chr_id):
                continue
            if chr_id == "MT":
                chr_id == "M"
            chr_id = "chr" + chr_id

        # Make start coordinate 0-base (BED standard).
        feat_s = feat_s - 1

        # Extract transcript ID and from infos.
        m = re.search('gene_id "(.+?)"', infos)
        assert m, "gene_id entry missing in GTF file \"%s\", line \"%s\"" %(in_gtf, line)
        gene_id = m.group(1)
        # Extract transcript ID.
        m = re.search('transcript_id "(.+?)"', infos)
        assert m, "transcript_id entry missing in GTF file \"%s\", line \"%s\"" %(in_gtf, line)
        transcript_id = m.group(1)
        # Extract exon number.
        m = re.search('exon_number "(\d+?)"', infos)
        assert m, "exon_number entry missing in GTF file \"%s\", line \"%s\"" %(in_gtf, line)
        exon_nr = int(m.group(1))
        # Extract gene name.
        m = re.search('gene_name "(.+?)"', infos)
        assert m, "gene_name entry missing in GTF file \"%s\", line \"%s\"" %(in_gtf, line)
        gene_name = m.group(1)
        # Extract gene biotype.
        m = re.search('gene_biotype "(.+?)"', infos)
        assert m, "gene_biotype entry missing in GTF file \"%s\", line \"%s\"" %(in_gtf, line)
        gene_biotype = m.group(1)

        # Check if transcript ID is in transcript dic.
        if not transcript_id in tr_ids_dic:
            continue

        # Check whether exon numbers are incrementing for each transcript ID.
        if not transcript_id in tr2exon_nr_dic:
            tr2exon_nr_dic[transcript_id] = exon_nr
        else:
            assert tr2exon_nr_dic[transcript_id] < exon_nr, "transcript ID \"%s\" without increasing exon number order in GTF file \"%s\"" %(transcript_id, in_gtf)
            tr2exon_nr_dic[transcript_id] = exon_nr

        # Make exon count 3-digit.
        #add = ""
        #if exon_nr < 10:
        #    add = "00"
        #if exon_nr >= 10 and exon_nr < 100:
        #    add = "0"

        # Count exon entry.
        c_gtf_ex_feat += 1

        # Construct exon ID.
        exon_id = transcript_id + "_e" + str(exon_nr)

        # Store more infos.
        tr2gene_name_dic[transcript_id] = gene_name
        tr2gene_biotype_dic[transcript_id] = gene_biotype
        tr2gene_id_dic[transcript_id] = gene_id
        exon_id_tr_dic[exon_id] = transcript_id

        # Store exon lengths in dictionary of lists.
        feat_l = feat_e - feat_s
        if not transcript_id in tr_exon_len_dic:
            tr_exon_len_dic[transcript_id] = [feat_l]
        else:
            tr_exon_len_dic[transcript_id].append(feat_l)
        # Store exon IDs in dictionary of lists.
        if not transcript_id in tr_exon_id_dic:
            tr_exon_id_dic[transcript_id] = [exon_id]
        else:
            tr_exon_id_dic[transcript_id].append(exon_id)

        # Output genomic exon region.
        OUTBED.write("%s\t%i\t%i\t%s\t0\t%s\n" % (chr_id,feat_s,feat_e,exon_id,feat_pol))

    OUTBED.close()
    f.close()

    # Check for read-in features.
    assert c_gtf_ex_feat, "ERROR: no exon features read in from \"%s\"" %(in_gtf)

    # Output transcript exon regions.
    OUTBED = open(transcript_exon_bed, "w")
    tr_exon_starts_dic = {}

    # Calculate transcript exon coordinates from in-order exon lengths.
    for tr_id in tr_exon_len_dic:
        start = 0
        for exon_i, exon_l in enumerate(tr_exon_len_dic[tr_id]):
            exon_id = tr_exon_id_dic[tr_id][exon_i]
            new_end = start + exon_l
            # Store exon transcript start positions (0-based).
            tr_exon_starts_dic[exon_id] = start
            OUTBED.write("%s\t%i\t%i\t%s\t0\t+\n" % (tr_id,start,new_end,exon_id))
            # Set for next exon.
            start = new_end
    OUTBED.close()

    # Get input .bed region lengths and scores.
    id2site_sc_dic = {}
    id2site_len_dic = {}
    with open(in_bed) as f:
        for line in f:
            cols = line.strip().split("\t")
            site_s = int(cols[1])
            site_e = int(cols[2])
            site_id = cols[3]
            site_sc = float(cols[4])
            site_pol = cols[5]
            assert site_pol == "+" or site_pol == "-", "invalid strand (in_bed: %s, site_pol: %s)" %(in_bed, site_pol)
            # Check whether score is whole number.
            if int_whole_nr:
                if not site_sc % 1:
                    site_sc = int(site_sc)
            # Store score and length of each genomic input site.
            id2site_sc_dic[site_id] = site_sc
            id2site_len_dic[site_id] = site_e - site_s
    f.close()

    # Number input sites.
    c_in_bed_sites = len(id2site_len_dic)

    # Calculate overlap between genome exon .bed and input .bed.
    intersect_params = "-s -wb -f %f" %(intersectBed_f)
    intersect_bed_files(in_bed, genome_exon_bed, intersect_params, overlap_out)

    # Calculate hit region transcript positions.

    # Store complete and incomplete hits in separate .bed files.
    OUTINC = open(incomplete_transcript_hits_bed, "w")
    OUTCOM = open(complete_transcript_hits_bed, "w")

    c_complete = 0
    c_incomplete = 0
    c_all = 0
    # Count site hits dic.
    c_site_hits_dic = {}
    # ID to site stats dic of lists.
    id2stats_dic = {}
    # ID to hit length dic.
    id2hit_len_dic = {}
    # Transcripts with hits dic.
    match_tr_dic = {}
    # Transcript ID to unique complete hits dic.
    tr2uniq_com_hits_dic = {}
    # Transcript ID to unique all hits dic.
    tr2uniq_all_hits_dic = {}
    # Transcript ID to complete hits dic.
    tr2com_hits_dic = {}
    # Transcript ID to all hits dic.
    tr2all_hits_dic = {}
    # Site ID to transcript ID dic.
    site2tr_id_dic = {}

    with open(overlap_out) as f:
        for line in f:
            cols = line.strip().split("\t")
            s_gen_hit = int(cols[1])
            e_gen_hit = int(cols[2])
            site_id = cols[3]
            s_gen_exon = int(cols[7])
            e_gen_exon = int(cols[8])
            exon_id = cols[9]
            exon_pol = cols[11]
            c_all += 1
            # Count how many transcriptome matches site has.
            if site_id in c_site_hits_dic:
                c_site_hits_dic[site_id] += 1
            else:
                c_site_hits_dic[site_id] = 1
            # Exon transcript start position (0-based).
            s_tr_exon = tr_exon_starts_dic[exon_id]
            # Hit length.
            l_gen_hit = e_gen_hit - s_gen_hit
            # Site length.
            l_site = id2site_len_dic[site_id]
            # Hit region transcript positions (plus strand).
            hit_tr_s_pos = s_gen_hit - s_gen_exon + s_tr_exon
            hit_tr_e_pos = hit_tr_s_pos + l_gen_hit
            # If exon on reverse (minus) strand.
            if exon_pol == "-":
                hit_tr_s_pos = e_gen_exon - e_gen_hit + s_tr_exon
                hit_tr_e_pos = hit_tr_s_pos + l_gen_hit
            # Site score.
            site_sc = id2site_sc_dic[site_id]
            # Transcript ID.
            tr_id = exon_id_tr_dic[exon_id]
            # Store transcript ID for each site ID.
            # In case site ID has several transcript hits, this will be overwritten,
            # but we just use this dic for unique hits so no problem.
            site2tr_id_dic[site_id] = tr_id
            # Store transcript ID has a match.
            match_tr_dic[tr_id] = 1
            # If site score round number, make integer.
            if not site_sc % 1:
                site_sc = int(site_sc)
            site_sc = str(site_sc)
            # Store hit stats list (bed row) for each site.
            bed_row = "%s\t%i\t%i\t%s\t%s\t+" %(tr_id, hit_tr_s_pos, hit_tr_e_pos, site_id, site_sc)
            if not site_id in id2stats_dic:
                id2stats_dic[site_id] = [bed_row]
            else:
                id2stats_dic[site_id].append(bed_row)
            id2hit_len_dic[site_id] = l_gen_hit
            if l_gen_hit == l_site:
                # Output complete hits.
                OUTCOM.write("%s\n" % (bed_row))
                # Count complete hits per transcript.
                if tr_id in tr2com_hits_dic:
                    tr2com_hits_dic[tr_id] += 1
                else:
                    tr2com_hits_dic[tr_id] = 1
            else:
                # Output incomplete hits.
                OUTINC.write("%s\n" % (bed_row))
            # Count all hits per transcript.
            if tr_id in tr2all_hits_dic:
                tr2all_hits_dic[tr_id] += 1
            else:
                tr2all_hits_dic[tr_id] = 1

    OUTCOM.close()
    OUTINC.close()
    f.close()

    # Output unique hits (two files, one for complete hits, other for all).
    OUTUNIALL = open(uniq_all_out, "w")
    OUTUNICOM = open(uniq_complete_out, "w")

    for site_id in c_site_hits_dic:
        c_hits = c_site_hits_dic[site_id]
        if c_hits != 1:
            continue
        l_hit = id2hit_len_dic[site_id]
        l_site = id2site_len_dic[site_id]
        tr_id = site2tr_id_dic[site_id]
        bed_row = id2stats_dic[site_id][0]
        if l_hit == l_site:
            # Store unique + complete hit.
            OUTUNICOM.write("%s\n" % (bed_row))
            if tr_id in tr2uniq_com_hits_dic:
                tr2uniq_com_hits_dic[tr_id] += 1
            else:
                tr2uniq_com_hits_dic[tr_id] = 1
        # Store unique hit (complete or incomplete).
        OUTUNIALL.write("%s\n" % (bed_row))
        if tr_id in tr2uniq_all_hits_dic:
            tr2uniq_all_hits_dic[tr_id] += 1
        else:
            tr2uniq_all_hits_dic[tr_id] = 1

    OUTUNICOM.close()
    OUTUNIALL.close()

    # For all transcripts with mapped regions, store exons.bed + stats.out.
    OUTEXBED = open(hit_tr_exons_bed, "w")
    OUTSTATS = open(hit_tr_stats_out, "w")
    # Statistics out file header.
    OUTSTATS.write("tr_id\tchr\tgen_s\tgen_e\tpol\tgene_id\tgene_name\tgene_biotype\ttr_len\tcomp_hits\tall_hits\tuniq_comp_hits\tuniq_all_hits\n")
    # transcript stats.
    tr2len_dic = {}
    tr2gen_s_dic = {}
    tr2gen_e_dic = {}
    tr2gen_chr_dic = {}
    tr2gen_pol_dic = {}

    with open(genome_exon_bed) as f:
        for line in f:
            cols = line.strip().split("\t")
            chr_id = cols[0]
            ex_s = int(cols[1])
            ex_e = int(cols[2])
            ex_id = cols[3]
            ex_pol = cols[5]
            ex_l = ex_e - ex_s
            # Print out exons of transcripts with hits.
            tr_id = exon_id_tr_dic[ex_id]
            # Store transcripts lengths.
            if tr_id in tr2len_dic:
                tr2len_dic[tr_id] += ex_l
            else:
                tr2len_dic[tr_id] = ex_l
            # Store more transcript stats.
            if tr_id in tr2gen_s_dic:
                if ex_s < tr2gen_s_dic[tr_id]:
                    tr2gen_s_dic[tr_id] = ex_s
            else:
                tr2gen_s_dic[tr_id] = ex_s
            if tr_id in tr2gen_e_dic:
                if ex_e > tr2gen_e_dic[tr_id]:
                    tr2gen_e_dic[tr_id] = ex_e
            else:
                tr2gen_e_dic[tr_id] = ex_e
            tr2gen_chr_dic[tr_id] = chr_id
            tr2gen_pol_dic[tr_id] = ex_pol
            if tr_id in match_tr_dic:
                bed_row = "%s\t%i\t%i\t%s\t0\t%s" %(chr_id, ex_s, ex_e, ex_id, ex_pol)
                OUTEXBED.write("%s\n" % (bed_row))
    OUTEXBED.close()

    # Transcript hit statistics.
    for tr_id in match_tr_dic:
        gene_id = tr2gene_id_dic[tr_id]
        gene_biotype = tr2gene_biotype_dic[tr_id]
        gene_name = tr2gene_name_dic[tr_id]
        tr_l = tr2len_dic[tr_id]
        tr_chr = tr2gen_chr_dic[tr_id]
        tr_pol = tr2gen_pol_dic[tr_id]
        tr_gen_s = tr2gen_s_dic[tr_id]
        tr_gen_e = tr2gen_e_dic[tr_id]
        c_com_hits = 0
        c_all_hits = 0
        c_uniq_com_hits = 0
        c_uniq_all_hits = 0
        if tr_id in tr2com_hits_dic:
            c_com_hits = tr2com_hits_dic[tr_id]
        if tr_id in tr2all_hits_dic:
            c_all_hits = tr2all_hits_dic[tr_id]
        if tr_id in tr2uniq_com_hits_dic:
            c_uniq_com_hits = tr2uniq_com_hits_dic[tr_id]
        if tr_id in tr2uniq_all_hits_dic:
            c_uniq_all_hits = tr2uniq_all_hits_dic[tr_id]
        stats_row = "%s\t%s\t%i\t%i\t%s\t%s\t%s\t%s\t%i\t%i\t%i\t%i\t%i" %(tr_id, tr_chr, tr_gen_s, tr_gen_e, tr_pol, gene_id, gene_name, gene_biotype, tr_l, c_com_hits, c_all_hits, c_uniq_com_hits, c_uniq_all_hits)
        OUTSTATS.write("%s\n" % (stats_row))
    OUTSTATS.close()


################################################################################

def bed_merge_identical_regions(in_bed, out_bed,
                                out_uniq = False):
    """
    Merge identical regions inside in_bed .bed file and output to out_bed.
    Store IDs of merged regions in column 4.

    out_uniq:
    Only output unique regions

    >>> in_bed = "test_data/test5.bed"
    >>> out_bed = "test_data/test5.tmp.bed"
    >>> exp_bed = "test_data/test5.exp.bed"
    >>> bed_merge_identical_regions(in_bed, out_bed)
    >>> diff_two_files_identical(out_bed, exp_bed)
    True

    """
    assert os.path.isfile(in_bed), "ERROR: Cannot open in_bed \"%s\"" % (in_bed)
    # Read in IDs.
    reg2ids_dic = {}
    reg2c_dic = {}
    with open(in_bed) as f:
        for line in f:
            row = line.strip()
            cols = line.strip().split("\t")
            chr_id = cols[0]
            site_s = cols[1]
            site_e = cols[2]
            site_id = cols[3]
            site_pol = cols[5]
            reg_id = "%s,%s,%s,%s" %(chr_id, site_s, site_e, site_pol)
            if reg_id in reg2ids_dic:
                reg2ids_dic[reg_id] += ",%s" %(site_id)
                reg2c_dic[reg_id] += 1
            else:
                reg2ids_dic[reg_id] = site_id
                reg2c_dic[reg_id] = 1
    f.closed
    assert reg2ids_dic, "ERROR: reg2ids_dic empty. in_bed empty?"
    # Write merged regions to out_bed.
    OUTBED = open(out_bed,"w")
    c_out = 0
    for reg in reg2ids_dic:
        reg_c = reg2c_dic[reg]
        reg_ids = reg2ids_dic[reg]
        if out_uniq:
            if reg_c != 1:
                continue
        cols = reg.strip().split(",")
        c_out += 1
        OUTBED.write("%s\t%s\t%s\t%s\t0\t%s\n" % (cols[0],cols[1],cols[2],reg_ids,cols[3]))
    OUTBED.close()
    assert c_out, "ERROR: nothing was output into out_bed"


################################################################################

def bed_convert_transcript_to_genomic_sites(in_bed, in_gtf, out_bed,
                                            site2hitc_dic=None,
                                            out_folder=False):
    """
    Dependencies:
    bedtools (tested with 2.26.0)
    gzip

    Convert in_bed .bed file with transcript sites into genomic coordinates
    sites file. in_bed column 1 transcript IDs have to be present in
    in_gtf GTF file, from which genomic exon coordinates of the transcript
    will get extracted.

    site2hitc_dic:
    A site2hitc_dic can be given, where site ID to hit count will be stored
    for usage outside the function.

    Output:
    By default output to out_bed file, using id_p1, id_p2 IDs.
    If out_folder=True, use out_bed name as folder name.
    In this case, output these files to folder:
    exon_regions_genome.bed
    exon_regions_transcript.bed
    complete_hits.bed
    split_hits.bed
    all_hits.bed

    >>> test_gtf = "test_data/test_tr2gen.gtf"
    >>> test_in_bed = "test_data/test_tr2gen.bed"
    >>> test_out_exp_bed = "test_data/test_tr2gen.exp.bed"
    >>> test_out_tmp_bed = "test_data/test_tr2gen.tmp.bed"
    >>> bed_convert_transcript_to_genomic_sites(test_in_bed, test_gtf, test_out_tmp_bed)
    >>> diff_two_files_identical(test_out_exp_bed, test_out_tmp_bed)
    True
    >>> test_out = "test_data/tr2gen_tmp_out"
    >>> test_out_tmp_bed = "test_data/tr2gen_tmp_out/all_hits.bed"
    >>> bed_convert_transcript_to_genomic_sites(test_in_bed, test_gtf, test_out, out_folder=True)
    >>> diff_two_files_identical(test_out_exp_bed, test_out_tmp_bed)
    True

    """

    # Generate .tmp files.
    random_id = uuid.uuid1()
    tmp_bed = str(random_id) + ".tmp.bed"
    random_id = uuid.uuid1()
    tmp_out = str(random_id) + ".tmp.out"

    # Output files if output_folder=True.
    if out_folder:
        if not os.path.exists(out_bed):
            os.makedirs(out_bed)
    out_exon_regions_genome_bed = out_bed + "/" + "exon_regions_genome.bed"
    out_exon_regions_transcript_bed = out_bed + "/" + "exon_regions_transcript.bed"
    out_unique_hits_bed = out_bed + "/" + "unique_hits.bed"
    out_split_hits_bed = out_bed + "/" + "split_hits.bed"
    out_all_hits_bed = out_bed + "/" + "all_hits.bed"

    # Transcript IDs dic.
    tr_ids_dic = bed_get_chromosome_ids(in_bed)

    # Extract transcript exon regions from GTF and store as BED.
    gtf_extract_exon_bed(in_gtf, tmp_bed, tr_ids_dic=tr_ids_dic)
    if out_folder:
        make_file_copy(tmp_bed, out_exon_regions_transcript_bed)

    # Get exon region lengths.
    exid2len_dic = bed_get_region_lengths(tmp_bed)

    # Get exon numbers for each transcript.
    tr_exc_dic = bed_get_transcript_exon_numbers(tmp_bed)

    # Read in exon region stats.
    id2chr_dic = {}
    id2s_dic = {}
    id2e_dic = {}
    id2pol_dic = {}
    exid2trid_dic = {}
    with open(tmp_bed) as f:
        for line in f:
            row = line.strip()
            cols = line.strip().split("\t")
            chr_id = cols[0]
            site_s = int(cols[1])
            site_e = int(cols[2])
            site_id = cols[3]
            site_pol = cols[5]
            id2chr_dic[site_id] = chr_id
            id2s_dic[site_id] = site_s
            id2e_dic[site_id] = site_e
            id2pol_dic[site_id] = site_pol
            if re.search(".+_e\d", site_id):
                m = re.search("(.+)_e\d", site_id)
                tr_id = m.group(1)
                exid2trid_dic[site_id] = tr_id
            else:
                assert False, "ERROR: site ID \"%s\" missing added _e exon number" %(site_id)
    f.close()

    # Output exon regions with transcript coordinates.
    OUTBED = open(tmp_bed, "w")
    for tr_id in tr_exc_dic:
        ex_c = tr_exc_dic[tr_id]
        new_s = 0
        for i in range(ex_c):
            i += 1
            ex_id = tr_id + "_e" + str(i)
            gen_s = id2s_dic[ex_id]
            gen_e = id2e_dic[ex_id]
            ex_len = gen_e - gen_s
            tr_s = new_s
            tr_e = new_s + ex_len
            OUTBED.write("%s\t%i\t%i\t%s\t0\t+\n" % (tr_id,tr_s,tr_e,ex_id))
            new_s = tr_e
    OUTBED.close()

    if out_folder:
        make_file_copy(tmp_bed, out_exon_regions_genome_bed)

    # Overlap in_bed with tmp_bed.
    params = "-wb"
    intersect_bed_files(in_bed, tmp_bed, params, tmp_out,
                        sorted_out=True)

    # Read in transcript site overlaps with transcript exon regions.
    site2c_dic = {}
    # Dictionaries for later outputting unique + split hits separately.
    siteid2pol_dic = {}
    siteid2sc_dic = {}
    partid2chrse_dic = {}
    with open(tmp_out) as f:
        for line in f:
            row = line.strip()
            cols = line.strip().split("\t")
            tr_id = cols[0]
            part_s = int(cols[1])
            part_e = int(cols[2])
            site_id = cols[3]
            site_sc = cols[4]
            ex_s = int(cols[7])
            ex_e = int(cols[8])
            ex_id = cols[9]
            ex_pol = id2pol_dic[ex_id]
            siteid2pol_dic[site_id] = ex_pol
            siteid2sc_dic[site_id] = site_sc
            if site_id in site2c_dic:
                site2c_dic[site_id] += 1
            else:
                site2c_dic[site_id] = 1
            # Hit part number.
            hit_c = site2c_dic[site_id]
            # Calculate genomic hit coordinates.
            # Plus strand case.
            gen_s = id2s_dic[ex_id] + part_s - ex_s
            gen_e = id2s_dic[ex_id] + part_e - ex_s
            # Minus strand case.
            if ex_pol == "-":
                gen_s = id2e_dic[ex_id] - part_e + ex_s
                gen_e = id2e_dic[ex_id] - part_s + ex_s
            # part ID.
            part_id = site_id + "_p" + str(hit_c)
            # Store chrse for each part ID.
            chrse = "%s\t%i\t%i" %(id2chr_dic[ex_id],gen_s,gen_e)
            partid2chrse_dic[part_id] = "%s\t%i\t%i" %(id2chr_dic[ex_id],gen_s,gen_e)

    # Produce seperate output files for unique + split hits.
    all_hits_bed = out_bed
    if out_folder:
        all_hits_bed = out_all_hits_bed
    ALLBED = open(all_hits_bed, "w")
    if out_folder:
        UNIBED = open(out_unique_hits_bed, "w")
        SPLBED = open(out_split_hits_bed, "w")
        for site_id in site2c_dic:
            hit_c = site2c_dic[site_id]
            if site2hitc_dic is not None:
                site2hitc_dic[site_id] = hit_c
            site_pol = siteid2pol_dic[site_id]
            site_sc = siteid2sc_dic[site_id]
            # For unique hit use site ID, for split hits use part IDs.
            if hit_c == 1:
                # Unique hits.
                part_id = site_id + "_p1"
                UNIBED.write("%s\t%s\t%s\t%s\n" %(partid2chrse_dic[part_id],site_id,site_sc,site_pol))
            else:
                # Split hits.
                for i in range(hit_c):
                    i += 1
                    part_id = site_id + "_p" + str(i)
                    SPLBED.write("%s\t%s\t%s\t%s\n" %(partid2chrse_dic[part_id],part_id,site_sc,site_pol))
    # Output all hits.
    for site_id in site2c_dic:
        hit_c = site2c_dic[site_id]
        if site2hitc_dic is not None:
            site2hitc_dic[site_id] = hit_c
        site_pol = siteid2pol_dic[site_id]
        site_sc = siteid2sc_dic[site_id]
        # For unique hit use site ID, for split hits use part IDs.
        if hit_c == 1:
            # Unique hits.
            part_id = site_id + "_p1"
            ALLBED.write("%s\t%s\t%s\t%s\n" %(partid2chrse_dic[part_id],site_id,site_sc,site_pol))
        else:
            # Split hits.
            for i in range(hit_c):
                i += 1
                part_id = site_id + "_p" + str(i)
                ALLBED.write("%s\t%s\t%s\t%s\n" %(partid2chrse_dic[part_id],part_id,site_sc,site_pol))

    # Delete tmp files.
    if os.path.exists(tmp_bed):
        os.remove(tmp_bed)
    if os.path.exists(tmp_out):
        os.remove(tmp_out)


################################################################################

def bed_get_chromosome_ids(bed_file):
    """
    Read in .bed file, return chromosome IDs (column 1 IDs).
    Return dic with chromosome ID -> count mapping.

    >>> test_file = "test_data/test6.bed"
    >>> bed_get_chromosome_ids(test_file)
    {'chr1': 2, 'chr2': 2, 'chr3': 1}

    """
    ids_dic = {}
    with open(bed_file) as f:
        for line in f:
            row = line.strip()
            cols = line.strip().split("\t")
            chr_id = cols[0]
            if chr_id in ids_dic:
                ids_dic[chr_id] += 1
            else:
                ids_dic[chr_id] = 1
    f.closed
    assert ids_dic, "No chromosome IDs read into dictionary (input file \"%s\" empty or malformatted?)" % (bed_file)
    return ids_dic


################################################################################

def bed_get_transcript_exon_numbers(in_bed):
    """
    Get number of exons for each transcript from in_bed BED file with
    transcript exon regions, with ID format:
    transcriptid_e1 (exon 1), transcriptid_e1 (exon 2)
    This is the output format from gtf_extract_exon_bed(), so both can
    be used in combination.

    >>> in_bed = "test_data/test6.bed"
    >>> bed_get_transcript_exon_numbers(in_bed)
    {'ENST1': 2, 'ENST2': 2, 'ENST3': 1}

    """
    tr_exc_dic = {}
    # Open input .bed file.
    with open(in_bed) as f:
        for line in f:
            cols = line.strip().split("\t")
            site_id = cols[3]
            if re.search(".+_e\d", site_id):
                m = re.search("(.+)_e\d", site_id)
                tr_id = m.group(1)
                if tr_id not in tr_exc_dic:
                    tr_exc_dic[tr_id] = 1
                else:
                    tr_exc_dic[tr_id] += 1
            else:
                assert False, "ERROR: site ID \"%s\" missing added _e exon number" %(site_id)
    f.close()
    assert tr_exc_dic, "ERROR: nothing was read in (\"%s\" empty or malformatted?)" %(in_bed)
    return tr_exc_dic


################################################################################

def get_transcript_sequences_from_gtf(in_gtf, in_2bit,
                                      tr_ids_dic=False):
    """
    Get spliced transcript sequences based on in_gtf annotations. For
    transcripts with > 1 exon, concatenate the exon sequences to build
    the transcript sequence. If one exon is missing / not extracted or
    if extracted lengths don't fit, the transcript sequence will be
    skipped / not output.
    Return dictionary with transcript_id -> sequence mapping.

    tr_ids_dic:
    Defines transcript IDs for which sequence should be extracted.

    """
    # Generate .tmp files.
    random_id = uuid.uuid1()
    tmp_bed = str(random_id) + ".tmp.bed"
    random_id = uuid.uuid1()
    tmp_fa = str(random_id) + ".tmp.fa"

    # Transcript sequences dic.
    tr_seqs_dic = {}

    # Extract transcript exon regions from GTF and store as BED.
    gtf_extract_exon_bed(in_gtf, tmp_bed, tr_ids_dic=tr_ids_dic)

    # Extract exon region sequences from .2bit.
    bed_extract_sequences_from_2bit(tmp_bed, tmp_fa, in_2bit)

    # Get transcript lengths from tmp_bed for comparison.
    tr_len_dic = bed_get_transcript_lengths_from_exon_regions(tmp_bed)
    # Get exon numbers for each transcript.
    tr_exc_dic = bed_get_transcript_exon_numbers(tmp_bed)

    # Read in sequences.
    exon_seqs_dic = read_fasta_into_dic(tmp_fa)

    # Concatenate exon region sequences.
    for tr_id in tr_exc_dic:
        ex_c = tr_exc_dic[tr_id]
        for i in range(ex_c):
            i += 1
            ex_id = tr_id + "_e" + str(i)
            if ex_id in exon_seqs_dic:
                ex_seq = exon_seqs_dic[ex_id]
                if tr_id not in tr_seqs_dic:
                    tr_seqs_dic[tr_id] = ex_seq
                else:
                    tr_seqs_dic[tr_id] += ex_seq
            else:
                print("WARNING: no sequence extracted for exon ID \"%s\". Skipping \"%s\" .. " %(ex_id, tr_id))
                if tr_id in tr_seqs_dic:
                    del tr_seqs_dic[tr_id]
                break
    # Checks.
    assert tr_seqs_dic, "ERROR: tr_seqs_dic empty (no FASTA sequences extracted?)"
    for tr_id in tr_seqs_dic:
        tr_len = len(tr_seqs_dic[tr_id])
        exp_len = tr_len_dic[tr_id]
        assert tr_len == exp_len, "ERROR: BED transcript length != FASTA transcript length for \"%s\"" %(tr_id)

    # Delete tmp files.
    if os.path.exists(tmp_bed):
        os.remove(tmp_bed)
    if os.path.exists(tmp_fa):
        os.remove(tmp_fa)

    # Return transcript sequences dic constructed from exon sequences.
    return tr_seqs_dic


################################################################################

def bed_get_transcript_lengths_from_exon_regions(in_bed):
    """
    Get spliced transcript lengths from in_bed BED file with transcript
    exon regions, with ID format:
    transcriptid_e1 (exon 1), transcriptid_e1 (exon 2)
    This is the output format from gtf_extract_exon_bed(), so both can
    be used in combination.

    >>> in_bed = "test_data/test6.bed"
    >>> bed_get_transcript_lengths_from_exon_regions(in_bed)
    {'ENST1': 4000, 'ENST2': 1500, 'ENST3': 2500}

    """
    tr_len_dic = {}
    # Open input .bed file.
    with open(in_bed) as f:
        for line in f:
            cols = line.strip().split("\t")
            site_s = int(cols[1])
            site_e = int(cols[2])
            site_id = cols[3]
            site_len = site_e - site_s
            if re.search(".+_e\d", site_id):
                m = re.search("(.+)_e\d", site_id)
                tr_id = m.group(1)
                if tr_id not in tr_len_dic:
                    tr_len_dic[tr_id] = site_len
                else:
                    tr_len_dic[tr_id] += site_len
            else:
                assert False, "ERROR: site ID \"%s\" missing added _e exon number" %(site_id)
    f.close()
    assert tr_len_dic, "ERROR: nothing was read in (\"%s\" empty or malformatted?)" %(in_bed)
    return tr_len_dic


################################################################################

def get_chromosome_lengths_from_2bit(in_2bit, out_lengths):
    """
    Get chromosome lengths from in_2bit .2bit file. Write lengths
    to out_lengths, with format:
    chr1	248956422
    chr10	133797422
    chr11	135086622
    ...
    Also return a dictionary with key=chr_id and value=chr_length.

    """

    # Check for twoBitInfo.
    assert is_tool("twoBitInfo"), "twoBitInfo not in PATH"

    # Run twoBitInfo and check.
    check_cmd = "twoBitInfo " + in_2bit + " " + out_lengths
    output = subprocess.getoutput(check_cmd)
    error = False
    if output:
        error = True
    assert error == False, "twoBitInfo is complaining:\n%s\n%s" %(check_cmd, output)

    # Read in lengths into dictionary.
    chr_len_dic = {}
    with open(out_lengths) as f:
        for line in f:
            row = line.strip()
            cols = line.strip().split("\t")
            chr_id = cols[0]
            chr_l = int(cols[1])
            assert chr_id not in chr_len_dic, "ERROR: non-unique chromosome ID \"%s\" encountered in \"%s\"" %(chr_id, out_lengths)
            chr_len_dic[chr_id] = chr_l
    f.closed
    assert chr_len_dic, "ERROR: chr_len_dic empty (\"%s\" empty?)" %(out_lengths)

    # Return chromosome lengths dic.
    return chr_len_dic


################################################################################

def move_rename_file(in_file, out_file):
    """
    Move / rename in_file to out_file.

    """
    check_cmd = "mv " + in_file + " " + out_file
    assert in_file != out_file, "mv does not like to mv file into same file (%s)" %(check_cmd)
    output = subprocess.getoutput(check_cmd)
    error = False
    if output:
        error = True
    assert error == False, "mv did not like your input (in_file: %s, out_file: %s):\n%s" %(in_file, out_file, output)


################################################################################

def merge_files(files_list, out_file):
    """
    Merge list of files into one output file.

    """
    assert files_list, "ERROR: given files_list is empty"
    for f in files_list:
        assert os.path.isfile(f), "ERROR: list file \"%s\" not found" % (f)
        assert f != out_file, "cat does not like to cat file into same file (%s)" %(check_cmd)
        check_cmd = "cat " + f + " >> " + out_file
        output = subprocess.getoutput(check_cmd)
        error = False
        if output:
            error = True
        assert error == False, "cat did not like your input (in_file: %s, out_file: %s):\n%s" %(f, out_file, output)


################################################################################

def fasta_output_dic(fasta_dic, fasta_out,
                     split=False,
                     split_size=60):
    """
    Output FASTA sequences dictionary (sequence_id -> sequence) to fasta_out.

    split        Split FASTA sequence for output to file
    split_size   Split size

    >>> fasta_dic = {'seq1': 'ACGTACGTACGTAC', 'seq2': 'ACGT'}
    >>> split_size = 4
    >>> fasta_exp = "test_data/test5.exp.fa"
    >>> fasta_out = "test_data/test5.tmp.fa"
    >>> fasta_output_dic(fasta_dic, fasta_out, split=True, split_size=split_size)
    >>> diff_two_files_identical(fasta_exp, fasta_out)
    True

    """
    # Check.
    assert fasta_dic, "ERROR: given fasta_dic empty"
    # Write sequences to FASTA file.
    OUTFA = open(fasta_out,"w")
    for seq_id in fasta_dic:
        seq = fasta_dic[seq_id]
        if split:
            OUTFA.write(">%s\n" %(seq_id))
            for i in range(0, len(seq), split_size):
                OUTFA.write("%s\n" %((seq[i:i+split_size])))
        else:
            OUTFA.write(">%s\n%s\n" %(seq_id, seq))
    OUTFA.close()


################################################################################

def gtf_get_transcript_lengths(in_gtf,
                               tr2exc_dic=None):
    """
    Get transcript lengths (= length of their exons, not unspliced length!)
    from GTF file.

    >>> in_gtf = "test_data/map_test_in.gtf"
    >>> gtf_get_transcript_lengths(in_gtf)
    {'ENST001': 2000, 'ENST002': 2000}

    """
    # Transcript ID to exonic length dictionary.
    tr2len_dic = {}
    # Open GTF either as .gz or as text file.
    if re.search(".+\.gz$", in_gtf):
        f = gzip.open(in_gtf, 'rt')
    else:
        f = open(in_gtf, "r")
    for line in f:
        # Skip header.
        if re.search("^#", line):
            continue
        cols = line.strip().split("\t")
        feature = cols[2]
        feat_s = int(cols[3])
        feat_e = int(cols[4])
        infos = cols[8]
        if not feature == "exon":
            continue
        # Extract transcript ID.
        m = re.search('transcript_id "(.+?)"', infos)
        assert m, "transcript_id entry missing in GTF file \"%s\", line \"%s\"" %(in_gtf, line)
        tr_id = m.group(1)
        # Sum up length.
        ex_len = feat_e - feat_s + 1
        if not tr_id in tr2len_dic:
            tr2len_dic[tr_id] = ex_len
        else:
            tr2len_dic[tr_id] += ex_len
        if tr2exc_dic is not None:
            if not tr_id in tr2exc_dic:
                tr2exc_dic[tr_id] = 1
            else:
                tr2exc_dic[tr_id] += 1
    f.close()
    assert tr2len_dic, "No IDs read into dictionary (input file \"%s\" empty or malformatted?)" % (in_gtf)
    return tr2len_dic


################################################################################

def gtf_extract_most_prominent_transcripts(in_gtf, out_file,
                                           strict=False,
                                           min_len=False,
                                           add_infos=False):
    """
    Extract most prominent transcripts list from in_gtf.

    in_gtf:
    Genomic annotations (hg38) GTF file (.gtf or .gtf.gz)
    NOTE: tested with Ensembl GTF files, expects transcript
    support level (TSL) information.
    out_file:
    File to output transcript IDs (optionally with add_infos)
    min_len:
    Accept only transcripts with length >= --min-len
    strict:
    Accept only transcripts with transcript support level (TSL) 1-5
    add_infos:
    Add additional information columns (gene ID, TSL, length) to out_list
    output file.

    """

    # Comparison dictionary.
    id2sc = {}
    for i in range(5):
        pos = i + 1
        pos_str = "%i" %(pos)
        id2sc[pos_str] = pos
    id2sc["NA"] = 6

    if strict:
        print("Strict transcript selection enabled ... ")
    if add_infos:
        print("Additional transcript infos in output file enabled ... ")

    # Read in transcript length (exonic regions).
    print("Read in transcript lengths (exonic lengths) from GTF ... ")
    tr2exc_dic = {}
    tr2len_dic = gtf_get_transcript_lengths(in_gtf, tr2exc_dic=tr2exc_dic)
    assert tr2len_dic, "ERROR: no transcript lengths read in from --gtf (invalid file format?)"
    print("# transcripts read in:  %i" %(len(tr2len_dic)))

    # Store most prominent transcript.
    g2tr_id = {}
    g2tr_tsl = {}
    g2tr_len = {}
    g2tr_bt = {}
    g2gn = {}
    g2gbt = {}

    print("Extract most prominent transcripts ... ")

    # Open GTF either as .gz or as text file.
    if re.search(".+\.gz$", in_gtf):
        f = gzip.open(in_gtf, 'rt')
    else:
        f = open(in_gtf, "r")
    for line in f:
        # Skip header.
        if re.search("^#", line):
            continue
        cols = line.strip().split("\t")
        chr_id = cols[0]
        feature = cols[2]
        feat_s = int(cols[3])
        feat_e = int(cols[4])
        feat_pol = cols[6]
        infos = cols[8]
        if not feature == "transcript":
            continue

        # Extract gene ID.
        m = re.search('gene_id "(.+?)"', infos)
        assert m, "gene_id entry missing in GTF file \"%s\", line \"%s\"" %(in_gtf, line)
        gene_id = m.group(1)
        # Extract transcript ID.
        m = re.search('transcript_id "(.+?)"', infos)
        assert m, "transcript_id entry missing in GTF file \"%s\", line \"%s\"" %(in_gtf, line)
        tr_id = m.group(1)
        # Extract gene name.
        m = re.search('gene_name "(.+?)"', infos)
        assert m, "gene_name entry missing in GTF file \"%s\", line \"%s\"" %(in_gtf, line)
        gene_name = m.group(1)
        # Extract gene biotype.
        m = re.search('gene_biotype "(.+?)"', infos)
        assert m, "gene_biotype entry missing in GTF file \"%s\", line \"%s\"" %(in_gtf, line)
        gene_biotype = m.group(1)
        # Extract transcript biotype.
        m = re.search('transcript_biotype "(.+?)"', infos)
        assert m, "transcript_biotype entry missing in GTF file \"%s\", line \"%s\"" %(in_gtf, line)
        tr_biotype = m.group(1)

        # Transcript length.
        tr_len = tr2len_dic[tr_id]
        # Gene name.
        g2gn[gene_id] = gene_name
        # Gene biotype.
        g2gbt[gene_id] = gene_biotype

        # Look for basic tag.
        m = re.search('tag "basic"', infos)
        if not m:
            continue
        # Get transcript support level (TSL).
        m = re.search('transcript_support_level "(.+?)"', infos)
        tsl_id = "NA"
        if m:
            tsl_id = m.group(1)
            if re.search("assigned to previous", tsl_id):
                m = re.search("(.+?) \(", tsl_id)
                tsl_id = m.group(1)

        # More filtering.
        if strict:
            if tsl_id == "NA":
                continue
        if min_len:
            if tr_len < min_len:
                continue

        # Update most prominent transcript.
        if not gene_id in g2tr_id:
            g2tr_id[gene_id] = tr_id
            g2tr_len[gene_id] = tr_len
            g2tr_tsl[gene_id] = tsl_id
            g2tr_bt[gene_id] = tr_biotype
        else:
            if id2sc[tsl_id] < id2sc[g2tr_tsl[gene_id]]:
                g2tr_id[gene_id] = tr_id
                g2tr_len[gene_id] = tr_len
                g2tr_tsl[gene_id] = tsl_id
                g2tr_bt[gene_id] = tr_biotype
            elif id2sc[tsl_id] == id2sc[g2tr_tsl[gene_id]]:
                if tr_len > g2tr_len[gene_id]:
                    g2tr_id[gene_id] = tr_id
                    g2tr_len[gene_id] = tr_len
                    g2tr_tsl[gene_id] = tsl_id
                    g2tr_bt[gene_id] = tr_biotype
    f.close()

    assert g2tr_id, "No IDs read into dictionary (input file \"%s\" empty or malformatted?)" % (in_gtf)
    c_prom_tr = len(g2tr_id)
    print("Number of selected transcripts: %i" %(c_prom_tr))

    # Output transcript IDs list.
    OUT = open(out_file, "w")
    if add_infos:
        OUT.write("gene_id\tgene_name\tgene_biotype\ttr_id\ttr_biotype\ttr_len\ttr_exc\ttsl\n")
    for gene_id in g2tr_id:
        tr_id = g2tr_id[gene_id]
        tr_len = g2tr_len[gene_id]
        tsl_id = g2tr_tsl[gene_id]
        tr_bt = g2tr_bt[gene_id]
        tr_exc = tr2exc_dic[tr_id]
        gene_name = g2gn[gene_id]
        gene_bt = g2gbt[gene_id]
        if add_infos:
            OUT.write("%s\t%s\t%s\t%s\t%s\t%i\t%i\t%s\n" % (gene_id,gene_name,gene_bt,tr_id,tr_bt,tr_len,tr_exc,tsl_id))
        else:
            OUT.write("%s\n" % (tr_id))
    OUT.close()

    if add_infos:
        print("%i transcript IDs + additional infos written to:\n%s" %(c_prom_tr, out_file))
    else:
        print("%i transcript IDs written to:\n%s" %(c_prom_tr, out_file))


################################################################################
