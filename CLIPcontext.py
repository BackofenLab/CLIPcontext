#!/usr/bin/env python3


from lib import cliplib
import argparse
import gzip
import sys
import os


"""

TOOL DEPENDENCIES
=================

python3
perl
bedtools (tested with version v2.26.0)
twoBitToFa


Take input .bed files, map to transcriptome.
By default 90% overlap is required for mapping.


python CLIPcontext.py -i 

-g /home/uhlm/Data/genome_2bit/hg38.2bit
-a /home/uhlm/Data/ensembl_data/Homo_sapiens.GRCh38.97.gtf.gz


"""

################################################################################

def setup_argument_parser():
    """Setup argparse parser."""
    help_description = """
    CLIPcontext maps genomic regions to the transcriptome and outputs 
    transcript and genomic region sequences.
    """
    # Define argument parser.
    p = argparse.ArgumentParser(add_help=False,
                                prog="CLIPcontext.py",
                                #usage="%(prog)s -i [input.bed] ... ",
                                description=help_description,
                                formatter_class=argparse.ArgumentDefaultsHelpFormatter)
                                # ArgumentDefaultsHelpFormatter
                                # RawTextHelpFormatter
    # Required arguments.
    p.add_argument("-h", "--help",
                   action="help",
                   help="Print help message")
    p.add_argument("-i",
                   dest="in_bed",
                   type=str,
                   required = True,
                   help = "Geomic regions (hg38) BED file (6 column format)")
    p.add_argument("-o",
                   dest="out_folder",
                   type=str,
                   required = True,
                   help = "Output results folder")
    p.add_argument("-f",
                   dest="in_fasta",
                   type=str,
                   required = True,
                   help = "Transcript sequences FASTA file")
    p.add_argument("-t",
                   dest="in_tr_list",
                   type=str,
                   required = True,
                   help = "Transcript sequence IDs list file")
    p.add_argument("-a",
                   dest="in_gtf",
                   type=str,
                   required = True,
                   help = "Genomic annotation (hg38) GTF file (.gtf.gz also supported)")
    p.add_argument("-g",
                   dest="in_genome_2bit",
                   type=str,
                   required = True,
                   help = "Genome sequence (hg38) .2bit file")
    # Optional arguments.
    p.add_argument("--thr",
                   dest = "score_thr",
                   type = float,
                   default = False,
                   help = "Site score threshold for filtering -i bed_file")
    p.add_argument("--seq-ext",
                   dest="us_ds_ext",
                   type = int,
                   default = 30,
                   help = "Up- and downstream extension of center positions on transcripts")
    p.add_argument("--rev-filter",
                   dest = "rev_filter",
                   default = False,
                   action = "store_true",
                   help = "Reverse filtering (keep values <= threshold)")
    p.add_argument("--merge-ext",
                   dest="merge_ext",
                   type = int,
                   default = 20,
                   help = "Extend regions mapped to transcripts by --merge-ext before running mergeBed to merge overlapping regions")
    p.add_argument("--min-overlap",
                   dest = "min_overlap",
                   type = float,
                   default = 0.9,
                   help = "Minimum exon overlap required for site to be added to transcript context set (= intersectBed -f parameter)")
    p.add_argument("--gen-uniq-ids",
                   dest = "gen_uniq_ids",
                   default = False,
                   action = "store_true",
                   help = "Generate unique column 4 IDs for -i .bed file entries")
    return p


################################################################################

if __name__ == '__main__':

    # Setup argparse.
    parser = setup_argument_parser()
    # Read in command line arguments.
    args = parser.parse_args()

    # Check for Linux.
    assert ('linux' in sys.platform), "Sorry but this tool runs on Linux only"
    # Check tool availability.
    if not is_tool("bedtools"):
        print("ERROR: bedtools not in PATH")
        sys.exit()
    if not is_tool("twoBitToFa"):
        print("ERROR: twoBitToFa not in PATH")
        sys.exit()
    # Check .bed for content.
    c_in_sites = cliplib.count_file_rows(args.in_bed)
    if not c_in_sites:
        print("ERROR: input .bed file \"%s\" is empty" %(args.in_bed))
        sys.exit()
    # Check .bed for 6 column format.
    if not cliplib.bed_check_six_col_format(bed_file):
        print("ERROR: input .bed file \"%s\" appears to be not in 6-column .bed format" %(args.in_bed))
        sys.exit()
    # Check for unique column 4 IDs.
    if not args.gen_uniq_ids:
        if not cliplib.bed_check_unique_ids(args.in_bed):
            print("ERROR: input .bed file \"%s\" column 4 IDs not unique (change or use --gen-uniq-ids option)" %(args.in_bed))
            sys.exit()

    # Generate results output folder.
    if not os.path.exists(args.out_folder):
        os.makedirs(args.out_folder)

    # Read in transcript ID list.
    tr_ids_dic = cliplib.read_ids_into_dic(args.in_tr_list)
    if tr_ids_c == 0:
        print("ERROR: no transcript IDs read in from \"%s\"" %(args.in_tr_list))
        sys.exit()
    print("# transcript IDs read in:                  %i" %(tr_ids_c))
    print("# input .bed sites:                        %i" %(c_in_sites))

    # Read in transcript sequences.
    tr_seqs_dic = cliplib.read_fasta_into_dic(args.in_fasta, ids_dic=tr_ids_dic)

    # Check for IDs present in list but not in FASTA file.
    lost_ids_dic = {}
    for tr_id in tr_ids_dic:
        if not tr_id in tr_seqs_dic:
            lost_ids_dic[tr_id] = 1
            print("WARNING: transcript ID \"%s\" not found in -f .fa file" %(tr_id))
    if lost_ids_dic:
        print("WARNING: sites on missing transcripts will be skipped!")

    # Filter sites by threshold.
    tmp_bed1 = generate_random_fn("bed")
    if args.score_thr:
        cliplib.bed_filter_by_col5_score(args.in_bed, tmp_bed1,
                                         score_threshold=args.score_thr,
                                         generate_unique_ids=args.gen_uniq_ids,
                                         reverse_filter=args.rev_filter)
    else:
        cliplib.make_file_copy(args.in_bed, tmp_bed1)

    # Number of remaining sites.
    c_filt_sites = cliplib.count_file_rows(tmp_bed1)
    print("# input .bed sites after --thr filtering:  %i" %(c_filt_sites))
    if not c_filt_sites:
        print("ERROR: no remaining BED regions after --thr filtering"
        sys.exit()




print "Read in -bed regions:                         $c_in_bed\n";
print "Remaining -bed regions after -thr filtering:  $c_filt_in_bed\n";

"""
import gzip

with gzip.open('input.gz','rt') as f:
    for line in f:
        print('got line', line)
"""

    
    # Get region IDs, lengths, strand polarities.
    id2pol_dic = {}
    id2len_dic = {}
    with open(args.in_bed) as f:
        for line in f:
            cols = line.strip().split("\t")
            start = int(cols[1])
            end = int(cols[2])
            region_id = cols[3]
            pol = cols[5]
            region_l = end - start
            id2pol_dic[region_id] = pol
            id2len_dic[region_id] = region_l
    f.closed
    

