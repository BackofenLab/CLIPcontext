#!/usr/bin/env python3

from distutils.spawn import find_executable
import subprocess
import argparse
import sys
import os

"""

TOOL DEPENDENCIES
=================

Tool:     bedtools intersect (aka intersectBed)
Version:  v2.26.0


TEST CALLS
==========

python bed_get_regions_near_exon_borders.py -in ../data/SERBP1_K562_rep1_sites_chr1_hg38.bed -ex ../data/transcript_exons_chr1_hg38.bed -out test_out.bed --thr 2

"""

################################################################################

def setup_argument_parser():
    """Setup argparse parser."""
    help_description = """
    Given genomic regions (CLIP peak regions + exons), report portion of peaks 
    within specified distance to exon borders.
    """
    # Define argument parser.
    p = argparse.ArgumentParser(add_help=False,
                                prog="bed_get_regions_near_exon_borders.py",
                                #usage="%(prog)s -i [input.bed] ... ",
                                description=help_description,
                                formatter_class=argparse.MetavarTypeHelpFormatter)

    # Argument groups.
    p_man = p.add_argument_group("REQUIRED ARGUMENTS")
    p_opt = p.add_argument_group("OPTIONAL ARGUMENTS")

    # Required arguments.
    p_opt.add_argument("-h", "--help",
                   action="help",
                   help="Print help message")
    p_man.add_argument("-in",
                   dest="in_bed",
                   type=str,
                   required = True,
                   help = "Genomic CLIP peak regions input .bed file (6-column format)")
    p_man.add_argument("-out",
                   dest="out_bed",
                   type=str,
                   required = True,
                   help = "Genomic CLIP peak regions near exon borders output .bed file")
    p_man.add_argument("-ex",
                   dest="in_exon_bed",
                   type=str,
                   required = True,
                   help = "Genomic exon regions .bed file for calculating distances to exon borders")
    # Optional arguments.
    p_opt.add_argument("--max-dist",
                   dest = "max_dist",
                   type = int,
                   default = 50,
                   help = "Maximum distance of CLIP region end to nearest exon end for CLIP region to still be output (default: 50)")
    p_opt.add_argument("--min-len",
                   dest = "min_len",
                   type = int,
                   default = 1,
                   help = "Minimum input site length for filtering -i bed_file (default: 1)")
    p_opt.add_argument("--max-len",
                   dest = "max_len",
                   type = int,
                   default = 200,
                   help = "Maximum input site length for filtering -i bed_file (default: 200)")
    p_opt.add_argument("--thr",
                   dest = "score_thr",
                   type = float,
                   default = None,
                   help = "Filter out -i .bed regions < --thr column 5 score (default: no filtering)")
    p_opt.add_argument("--rev-filter",
                   dest = "rev_filter",
                   default = False,
                   action = "store_true",
                   help = "Reverse filtering (keep values <= --thr and prefer sites with smaller values) (default: false)")
    return p


################################################################################

def is_tool(name):
    """Check whether tool "name" is in PATH."""
    return find_executable(name) is not None


################################################################################

def bed_check_unique_ids(bed_file):
    """
    Check whether .bed file (6 column format with IDs in column 4) 
    has unique column 4 IDs.

    """
    check_cmd = "cut -f 4 " + bed_file + " | sort | uniq -d"
    output = subprocess.getoutput(check_cmd)
    if output:
        return False
    else:
        return True


################################################################################

def bed_check_six_col_format(bed_file):
    """
    Check whether given .bed file has 6 columns.

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

def count_file_rows(in_file):
    """
    Count number of file rows for given input file.

    """
    check_cmd = "cat " + in_file + " | wc -l"
    output = subprocess.getoutput(check_cmd)
    row_count = int(output.strip())
    return row_count


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

if __name__ == '__main__':

    # Setup argparse.
    parser = setup_argument_parser()
    # Read in command line arguments.
    args = parser.parse_args()

    # Do some checking.
    if not "linux" in sys.platform:
        print("ERROR: please use Linux")
        sys.exit()
    # Check tool availability.
    if not is_tool("bedtools"):
        print("ERROR: bedtools not in PATH")
        sys.exit()
    # Check file inputs.
    if not os.path.exists(args.in_bed):
        print("ERROR: input .bed file \"%s\" not found" %(args.in_bed))
    if not os.path.exists(args.in_exon_bed):
        print("ERROR: input exon .bed file \"%s\" not found" %(args.in_exon_bed))
    # Check .bed for 6-column format.
    if not bed_check_six_col_format(args.in_bed):
        print("ERROR: input .bed file \"%s\" appears to be not in 6-column .bed format" %(args.in_bed))
        sys.exit()
    if not bed_check_six_col_format(args.in_exon_bed):
        print("ERROR: input exon .bed file \"%s\" appears to be not in 6-column .bed format" %(args.in_exon_bed))
        sys.exit()
    # Check for unique column 4 IDs.
    if not bed_check_unique_ids(args.in_bed):
        print("ERROR: input .bed file \"%s\" column 4 IDs not unique" %(args.in_bed))
        sys.exit()
    # Check .bed for content.
    c_in = count_file_rows(args.in_bed)
    if not c_in:
        print("ERROR: input .bed file \"%s\" is empty" %(args.in_bed))
        sys.exit()

    # First get regions inside exons (overlapping >= 90 % with them).
    tmp_bed1 = "sites_overlapping_with_exons.tmp.bed"
    params = "-s -u -wa -f 0.90"
    intersect_bed_files(args.in_bed, args.in_exon_bed, params, tmp_bed1)

    # Filter and extend overlapping sites.
    tmp_bed2 = "extended_sites.tmp.bed"
    # Output .bed file.
    TMPOUT = open(tmp_bed2,"w")
    c_ol = 0
    id2len_dic = {}
    id2stats_dic = {}
    with open(tmp_bed1) as f:
        for line in f:
            cols = line.strip().split("\t")
            seq_id = cols[0]
            site_s = int(cols[1])
            site_e = int(cols[2])
            site_id = cols[3]
            site_sc = float(cols[4])
            site_pol = cols[5]
            site_l = site_e - site_s
            # Filter by site score.
            if args.score_thr is not None:
                if args.rev_filter:
                    if site_sc > args.score_thr:
                        continue
                else:
                    if site_sc < args.score_thr:
                        continue
            # Filter by site length.
            if site_l > args.max_len:
                continue
            if site_l < args.min_len:
                continue
            new_s = site_s - args.max_dist - 1
            new_e = site_e + args.max_dist + 1
            new_l = new_e - new_s
            id2len_dic[site_id] = new_l
            c_ol += 1
            # Store original region.
            id2stats_dic[site_id] = "%s\t%i\t%i\t%s\t%f\t%s" %(seq_id,site_s,site_e,site_id,site_sc,site_pol)
            # Output extended region.
            TMPOUT.write("%s\t%i\t%i\t%s\t%f\t%s\n" % (seq_id,new_s,new_e,site_id,site_sc,site_pol))
    f.close()
    TMPOUT.close()

    # Overlap sites with exons, get bases overlapping.
    tmp_bed3 = "base_overlaps.tmp.bed"
    intersect_bed_files(tmp_bed2, args.in_exon_bed, "-s", tmp_bed3)

    # Output .bed.
    OUT = open(args.out_bed,"w")
    seen_dic = {}
    # Number of sites close to exon ends.
    c_close = 0

    # Get sites within border range.
    with open(tmp_bed3) as f:
        for line in f:
            cols = line.strip().split("\t")
            site_s = int(cols[1])
            site_e = int(cols[2])
            site_id = cols[3]
            if site_id in seen_dic:
                continue
            site_l = site_e - site_s
            full_l = id2len_dic[site_id]
            bed_row = id2stats_dic[site_id]
            if not full_l == site_l:
                c_close += 1
                OUT.write("%s\n" %(bed_row))
            seen_dic[site_id] = 1

    clean_up = True
    if clean_up:
        # Remove tmp files.
        if os.path.exists(tmp_bed1):
            os.remove(tmp_bed1)
        if os.path.exists(tmp_bed2):
            os.remove(tmp_bed2)
        if os.path.exists(tmp_bed3):
            os.remove(tmp_bed3)


    # Report results.
    print("Region stats (post-filtering)")
    print("=============================")
    print("Number of -in regions:                        %i" %(c_in))
    print("Number of -in regions overlapping with -ex:   %i" %(c_ol))
    print("Number of -in regions close to exon ends:     %i" %(c_close))




