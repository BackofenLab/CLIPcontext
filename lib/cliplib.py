#!/usr/bin/env python3

from distutils.spawn import find_executable
import subprocess
import statistics
import gzip
import uuid
import sys
import re
import os

"""

=====================
  OPEN FOR BUSINESS
=====================

Run doctests:

python3 -m doctest cliplib.py

AssertionError only, no ERROR

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

def bed_get_region_id_scores(in_bed_file):
    """
    Read in .bed file, and store scores for each region in dictionary 
    (unique column 4 ID and column 5 score have to be present).
    Return dictionary with mappings region ID -> region score

    >>> test_bed = "test_data/test3.bed"
    >>> id2sc_dic = bed_get_region_id_scores(test_bed)
    >>> id2sc_dic
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

def bed_filter_by_col5_score(in_bed_file, out_bed_file,
                             score_threshold=1,
                             disable_filter=False,
                             generate_unique_ids=False,
                             center_sites=False,
                             ext_lr=False,
                             reverse_filter=False,
                             id_prefix="CLIP"):

    """
    Filter .bed file by column 5 score (assuming the higher the better).
    Output to new .bed file. Optionally do reverse filtering, and generate 
    new unique column 4 IDs. You can also use the function to just extend 
    the sites, by disable_filter=True + set ext_lr to extension you want.
    Moreover, centering of sites before extension is possible by 
    center_sites=True.

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
            if not disable_filter:
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

def intersect_bed_files(a_file, b_file, params, out_file):
    """
    Intersect two .bed files, using intersectBed.

    """

    check_cmd = "intersectBed -a " + a_file + " -b " + b_file + " " + params + " > " + out_file
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


def convert_genome_positions_to_transcriptome(in_bed, out_folder,
                                              in_gtf, tr_ids_dic,
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
    >>> convert_genome_positions_to_transcriptome(in_bed, out_folder, in_gtf, tr_ids_dic)
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
        add = ""
        if exon_nr < 10:
            add = "00"
        if exon_nr >= 10 and exon_nr < 100:
            add = "0"

        # Count exon entry.
        c_gtf_ex_feat += 1
        
        # Construct exon ID.
        exon_id = transcript_id + "_e" + add + str(exon_nr)

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
    assert c_gtf_ex_feat, "no exon features read in from \"%s\"" %(in_gtf)

    # Output genomic exon regions.
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
            id2site_sc_dic[site_id] = float(site_sc)
            id2site_len_dic[site_id] = site_e - site_s
    f.close()

    # Number input sites.
    c_in_bed_sites = len(id2site_len_dic)

    # Calculate overlap between genome exon .bed and input .bed.
    intersect_params = "-s -wb"
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
            site2tr_id_dic[site_id] = tr_id
            # Store match for transcript ID yes.
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

    # Output unique hits (two files, one only complete hits, other all).
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
    OUTSTATS.write("tr_id\tgene_id\tgene_name\tgene_biotype\tcomp_hits\tall_hits\tuniq_comp_hits\tuniq_all_hits\n")

    with open(genome_exon_bed) as f:
        for line in f:
            cols = line.strip().split("\t")
            chr_id = cols[0]
            ex_s = int(cols[1])
            ex_e = int(cols[2])
            ex_id = cols[3]
            ex_pol = cols[5]
            # Print out exons of transcripts with hits.
            tr_id = exon_id_tr_dic[ex_id]
            if tr_id in match_tr_dic:
                bed_row = "%s\t%i\t%i\t%s\t0\t%s" %(chr_id, ex_s, ex_e, ex_id, ex_pol)
                OUTEXBED.write("%s\n" % (bed_row))
            # Transcript hit statistics.
            gene_id = tr2gene_id_dic[tr_id]
            gene_biotype = tr2gene_biotype_dic[tr_id]
            gene_name = tr2gene_name_dic[tr_id]
            c_com_hits = tr2com_hits_dic[tr_id]
            c_all_hits = tr2all_hits_dic[tr_id]
            c_uniq_com_hits = tr2uniq_com_hits_dic[tr_id]
            c_uniq_all_hits = tr2uniq_all_hits_dic[tr_id]
            stats_row = "%s\t%s\t%s\t%s\t%i\t%i\t%i\t%i" %(tr_id, gene_id, gene_name, gene_biotype, c_com_hits, c_all_hits, c_uniq_com_hits, c_uniq_all_hits)
            OUTSTATS.write("%s\n" % (stats_row))
    OUTEXBED.close()
    OUTSTATS.close()



"""


ENSG001
ENST001
ENSG002
ENST002

1       havana  exon    1001   2000   .       +       .       gene_id "ENSG001"; gene_version "1"; transcript_id "ENST001"; transcript_version "1"; exon_number "1"; gene_name "ABC1"; gene_source "havana"; gene_biotype "transcribed_unprocessed_pseudogene"; transcript_name "ABC1-101"; transcript_source "havana"; transcript_biotype "lncRNA"; exon_id "ENSE001"; exon_version "1"; tag "basic"; transcript_support_level "1";
1       havana  exon    3001   4000   .       +       .       gene_id "ENSG001"; gene_version "1"; transcript_id "ENST001"; transcript_version "1"; exon_number "2"; gene_name "ABC1"; gene_source "havana"; gene_biotype "transcribed_unprocessed_pseudogene"; transcript_name "ABC1-101"; transcript_source "havana"; transcript_biotype "lncRNA"; exon_id "ENSE002"; exon_version "1"; tag "basic"; transcript_support_level "1";

1       havana  exon    8001   9000   .       -       .       gene_id "ENSG002"; gene_version "1"; transcript_id "ENST002"; transcript_version "1"; exon_number "1"; gene_name "DEF1"; gene_source "havana"; gene_biotype "transcribed_unprocessed_pseudogene"; transcript_name "DEF1-101"; transcript_source "havana"; transcript_biotype "lncRNA"; exon_id "ENSE003"; exon_version "1"; tag "basic"; transcript_support_level "1";
1       havana  exon    6001   7000   .       -       .       gene_id "ENSG002"; gene_version "1"; transcript_id "ENST002"; transcript_version "1"; exon_number "2"; gene_name "DEF1"; gene_source "havana"; gene_biotype "transcribed_unprocessed_pseudogene"; transcript_name "DEF1-101"; transcript_source "havana"; transcript_biotype "lncRNA"; exon_id "ENSE004"; exon_version "1"; tag "basic"; transcript_support_level "1";





1       havana  exon    12613   12721   .       +       .       gene_id "ENSG00000223972"; gene_version "5"; transcript_id "ENST00000456328"; transcript_version "2"; exon_number "2"; gene_name "DDX11L1"; gene_source "havana"; gene_biotype "transcribed_unprocessed_pseudogene"; transcript_name "DDX11L1-202"; transcript_source "havana"; transcript_biotype "lncRNA"; exon_id "ENSE00003582793"; exon_version "1"; tag "basic"; transcript_support_level "1";
1       havana  exon    13221   14409   .       +       .       gene_id "ENSG00000223972"; gene_version "5"; transcript_id "ENST00000456328"; transcript_version "2"; exon_number "3"; gene_name "DDX11L1"; gene_source "havana"; gene_biotype "transcribed_unprocessed_pseudogene"; transcript_name "DDX11L1-202"; transcript_source "havana"; transcript_biotype "lncRNA"; exon_id "ENSE00002312635"; exon_version "1"; tag "basic"; transcript_support_level "1";





/home/uhlm/Data/ensembl_data/Homo_sapiens.GRCh38.97.gtf.gz
1	ensembl_havana	exon	28887124	28887210	.	+	.	gene_id "ENSG00000159023"; gene_version "21"; transcript_id "ENST00000373800"; transcript_version "7"; exon_number "1"; gene_name "EPB41"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "EPB41-207"; transcript_source "ensembl_havana"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS331"; exon_id "ENSE00001900393"; exon_version "1"; tag "basic"; transcript_support_level "1";
1	ensembl_havana	exon	28987448	28987905	.	+	.	gene_id "ENSG00000159023"; gene_version "21"; transcript_id "ENST00000373800"; transcript_version "7"; exon_number "2"; gene_name "EPB41"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "EPB41-207"; transcript_source "ensembl_havana"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS331"; exon_id "ENSE00001429864"; exon_version "1"; tag "basic"; transcript_support_level "1";


1       havana  exon    11869   12227   .       +       .       gene_id "ENSG00000223972"; gene_version "5"; transcript_id "ENST00000456328"; transcript_version "2"; exon_number "1"; gene_name "DDX11L1"; gene_source "havana"; gene_biotype "transcribed_unprocessed_pseudogene"; transcript_name "DDX11L1-202"; transcript_source "havana"; transcript_biotype "lncRNA"; exon_id "ENSE00002234944"; exon_version "1"; tag "basic"; transcript_support_level "1";
1       havana  exon    12613   12721   .       +       .       gene_id "ENSG00000223972"; gene_version "5"; transcript_id "ENST00000456328"; transcript_version "2"; exon_number "2"; gene_name "DDX11L1"; gene_source "havana"; gene_biotype "transcribed_unprocessed_pseudogene"; transcript_name "DDX11L1-202"; transcript_source "havana"; transcript_biotype "lncRNA"; exon_id "ENSE00003582793"; exon_version "1"; tag "basic"; transcript_support_level "1";
1       havana  exon    13221   14409   .       +       .       gene_id "ENSG00000223972"; gene_version "5"; transcript_id "ENST00000456328"; transcript_version "2"; exon_number "3"; gene_name "DDX11L1"; gene_source "havana"; gene_biotype "transcribed_unprocessed_pseudogene"; transcript_name "DDX11L1-202"; transcript_source "havana"; transcript_biotype "lncRNA"; exon_id "ENSE00002312635"; exon_version "1"; tag "basic"; transcript_support_level "1";


"""



"""
GTF more labels to exploit:

CDS (several per transcript possible of course since on different exons often)
start_codon
stop_codon
five_prime_utr
three_prime_utr

"""




################################################################################




