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

python gtf_extract_most_prominent_transcripts.py --gtf ../data/Homo_sapiens.GRCh38.98.gtf.gz --out test_list.out --strict

python gtf_extract_most_prominent_transcripts.py --gtf ../data/Homo_sapiens.GRCh38.98.gtf.gz --out test.out --add-infos --strict

"""

################################################################################

def setup_argument_parser():
    """Setup argparse parser."""
    help_description = """
    Extract the most prominent transcript for each gene from the given GTF file.
    For all transcripts with tag "basic", filter based on transcript support 
    level (TSL, highest priority) and transcript length (longer transcripts 
    preferred). Output transcript IDs to --out.
    """
    # Define argument parser.
    p = argparse.ArgumentParser(add_help=False,
                                prog="gtf_extract_most_prominent_transcripts.py",
                                description=help_description,
                                formatter_class=argparse.MetavarTypeHelpFormatter)

    # Argument groups.
    p_man = p.add_argument_group("REQUIRED ARGUMENTS")
    p_opt = p.add_argument_group("OPTIONAL ARGUMENTS")
    
    # Arguments.
    p_opt.add_argument("-h", "--help",
                   action="help",
                   help="Print help message")
    p_man.add_argument("--gtf",
                   dest="in_gtf",
                   type=str,
                   required = True,
                   help = "Input GTF file with genomic annotations")
    p_man.add_argument("--out",
                   dest="out_file",
                   type=str,
                   required = True,
                   help = "Output transcript IDs list file")
    p_opt.add_argument("--min-len",
                   dest="min_len",
                   type = int,
                   default = False,
                   help = "Accept only transcripts with length >= --min-len (default: false)")
    p_opt.add_argument("--strict",
                   dest="strict",
                   default = False,
                   action = "store_true",
                   help = "Accept only transcripts with TSL 1-5 (default: false)")
    p_opt.add_argument("--add-infos",
                   dest="add_infos",
                   default = False,
                   action = "store_true",
                   help = "Add additional information columns (gene ID, TSL, length) to output file (default: false)")
    return p


################################################################################

def gtf_get_transcript_lengths(in_gtf):
    """
    Get transcript lengths (= length of their exons, not unspliced length!) 
    from GTF file.
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
    f.close()
    assert tr2len_dic, "No IDs read into dictionary (input file \"%s\" empty or malformatted?)" % (in_gtf)
    return tr2len_dic


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
        print("ERROR: input .gtf file \"%s\" not found" %(args.in_gtf))

    # Comparison dictionary.
    id2sc = {}
    for i in range(5):
        pos = i + 1
        pos_str = "%i" %(pos)
        id2sc[pos_str] = pos
    id2sc["NA"] = 6

    # Read in transcript length (exonic regions).
    print("Read in transcript lengths (exonic lengths) ... ")
    tr2len_dic = gtf_get_transcript_lengths(args.in_gtf)

    # Store most prominent transcript.
    g2tr_id = {}
    g2tr_tsl = {}
    g2tr_len = {}
    g2tr_bt = {}
    g2gn = {}
    g2gbt = {}

    print("Extract most prominent transcripts ... ")

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
        if args.strict:
            if tsl_id == "NA":
                continue
        if args.min_len:
            if tr_len < args.min_len:
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

    assert g2tr_id, "No IDs read into dictionary (input file \"%s\" empty or malformatted?)" % (args.in_gtf)
    c_prom_tr = len(g2tr_id)
    print("Number of selected transcripts: %i" %(c_prom_tr))

    # Output transcript IDs list.
    OUT = open(args.out_file, "w")
    if args.add_infos:
        OUT.write("gene_id\tgene_name\tgene_biotype\ttr_id\ttr_biotype\ttr_len\ttsl\n")
    for gene_id in g2tr_id:
        tr_id = g2tr_id[gene_id]
        tr_len = g2tr_len[gene_id]
        tsl_id = g2tr_tsl[gene_id]
        tr_bt = g2tr_bt[gene_id]
        gene_name = g2gn[gene_id]
        gene_bt = g2gbt[gene_id]
        if args.add_infos:
            OUT.write("%s\t%s\t%s\t%s\t%s\t%i\t%s\n" % (gene_id,gene_name,gene_bt,tr_id,tr_bt,tr_len,tsl_id))
        else:
            OUT.write("%s\n" % (tr_id))
    OUT.close()

"""

1       havana  transcript      11869   14409   .       +       .       gene_id "ENSG00000223972"; gene_version "5"; transcript_id "ENST00000456328"; transcript_version "2"; gene_name "DDX11L1"; gene_source "havana"; gene_biotype "transcribed_unprocessed_pseudogene"; transcript_name "DDX11L1-202"; transcript_source "havana"; transcript_biotype "lncRNA"; tag "basic"; transcript_support_level "1";


1       havana  gene    11869   14409   .       +       .       gene_id "ENSG00000223972"; gene_version "5"; gene_name "DDX11L1"; gene_source "havana"; gene_biotype "transcribed_unprocessed_pseudogene";
1       havana  transcript      11869   14409   .       +       .       gene_id "ENSG00000223972"; gene_version "5"; transcript_id "ENST00000456328"; transcript_version "2"; gene_name "DDX11L1"; gene_source "havana"; gene_biotype "transcribed_unprocessed_pseudogene"; transcript_name "DDX11L1-202"; transcript_source "havana"; transcript_biotype "lncRNA"; tag "basic"; transcript_support_level "1";

"""

