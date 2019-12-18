#!/usr/bin/env python3

import subprocess
import argparse
import gzip
import sys
import os
import re

"""

TEST CALLS
==========

python gtf_extract_exon_regions.py --gtf ../data/Homo_sapiens.GRCh38.98.gtf.gz --out test_out.bed --tr test_list.out

"""

################################################################################

def setup_argument_parser():
    """Setup argparse parser."""
    help_description = """
    Extract exon regions in BED format for a given list of transcript IDs.
    """
    # Define argument parser.
    p = argparse.ArgumentParser(add_help=False,
                                prog="gtf_extract_exon_regions.py",
                                description=help_description,
                                formatter_class=argparse.MetavarTypeHelpFormatter)

    # Arguments.
    p.add_argument("-h", "--help",
                   action="help",
                   help="Print help message")
    p.add_argument("--gtf",
                   dest="in_gtf",
                   type=str,
                   required = True,
                   help = "Input GTF file with genomic annotations")
    p.add_argument("--tr",
                   dest="in_tr_list",
                   type=str,
                   required = True,
                   help = "Transcript sequence IDs list file")
    p.add_argument("--out",
                   dest="out_bed",
                   type=str,
                   required = True,
                   help = "Output BED file (6-column format) with exon regions")
    return p


################################################################################

def read_ids_into_dic(ids_file,
                      ids_dic=False):
    """
    Read in IDs list file, where each line stores one ID.

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

if __name__ == '__main__':

    # Setup argparse.
    parser = setup_argument_parser()
    # Read in command line arguments.
    args = parser.parse_args()
    # Do some checking.
    if not "linux" in sys.platform:
        print("ERROR: please use Linux")
        sys.exit()
    # Check file inputs.
    if not os.path.exists(args.in_gtf):
        print("ERROR: --gtf input file \"%s\" not found" %(args.in_gtf))
    if not os.path.exists(args.in_tr_list):
        print("ERROR: --tr input file \"%s\" not found" %(args.in_tr_list))

    # Read in transcript ID list.
    print("Read in transcript IDs ... ")
    tr_ids_dic = read_ids_into_dic(args.in_tr_list)
    tr_ids_c = len(tr_ids_dic)
    if tr_ids_c == 0:
        print("ERROR: no transcript IDs read in from \"%s\"" %(args.in_tr_list))

    print("# transcript IDs read in:  %i" %(tr_ids_c))

    OUTBED = open(args.out_bed, "w")
    c_out = 0

    print("Extract exon regions ... ")

    # Open GTF either as .gz or as text file.
    if re.search(".+\.gz$", args.in_gtf):
        f = gzip.open(args.in_gtf, 'rt')
    else: 
        f = open(args.in_gtf, "r")
    for line in f:
        # Skip header.
        if re.search("^#", line):
            continue
        cols = line.strip().split("\t")
        chr_id = cols[0]
        feature = cols[2]
        feat_s = int(cols[3]) - 1 # make 0-based.
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

        # Extract transcript ID.
        m = re.search('transcript_id "(.+?)"', infos)
        assert m, "transcript_id entry missing in GTF file \"%s\", line \"%s\"" %(in_gtf, line)
        tr_id = m.group(1)
        # Extract exon ID.
        m = re.search('exon_id "(.+?)"', infos)
        assert m, "exon_id entry missing in GTF file \"%s\", line \"%s\"" %(in_gtf, line)
        exon_id = m.group(1)
        # Extract exon number.
        m = re.search('exon_number "(.+?)"', infos)
        assert m, "exon_number entry missing in GTF file \"%s\", line \"%s\"" %(in_gtf, line)
        exon_nr = m.group(1)

        # Only exon regions of given transcript IDs.
        if not tr_id in tr_ids_dic:
            continue

        # Output exon region.
        bed_id = "%s,%s,e%s" %(tr_id,exon_id,exon_nr)
        c_out += 1
        OUTBED.write("%s\t%i\t%i\t%s\t0\t%s\n" % (chr_id,feat_s,feat_e,bed_id,feat_pol))
    f.close()
    OUTBED.close()

    print("# exon regions output:     %i" %(c_out))



