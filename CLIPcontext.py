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
    if not "linux" in sys.platform:
        print("ERROR: please use Linux")
        sys.exit()
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
    # Remove lost IDs from transcript IDs dic.
    for seq_id in lost_ids_dic:
        del tr_ids_dic[seq_id]

    # Filter sites by threshold.
    tmp_bed1 = args.out_folder + "/" + "filtered_input_sites.tmp.bed"
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

    # # First map full length .bed to transcriptome.
    full_length_out = args.out_folder + "/" + "transcript_map_full_length_out"
    cliplib.convert_genome_positions_to_transcriptome(tmp_bed1, full_length_out, 
                                                      args.in_gtf, tr_ids_dic)
    # All unique matches (complete + incomplete).
    uniq_hits_full_length_map = full_length_out + "/" + "transcript_hits_all_unique.bed"

    # Count unique hits.
    c_uniq_fl_hits = cliplib.count_file_rows(uniq_hits_full_length_map)
    print("# unique transcript hits:                  %i" %(c_uniq_fl_hits))
    if not c_uniq_fl_hits:
        print("ERROR: no unique transcript hits for given genomic .bed and transcripts"
        sys.exit()

    # Get scores for unique hits.
    id2sc_dic = cliplib.bed_get_region_id_scores(unique_hits_full_length_map)

    # Prolong unique hits.
    tmp_bed2 = args.out_folder + "/" + "prolonged_unique_transcript_hits.tmp.bed";
    tmp_bed3 = args.out_folder + "/" + "prolonged_unique_transcript_hits.sorted.tmp.bed";
    tmp_bed4 = args.out_folder + "/" + "prolonged_unique_transcript_hits.merged.tmp.bed";

    # Prolong unique hits by set merge extension.
    cliplib.bed_filter_by_col5_score(unique_hits_full_length_map, tmp_bed2,
                                     disable_filter=True,
                                     ext_lr=args.merge_ext,
                                     center_sites=False)
    

"""


# Sort BED file.
qx/sort -k1,1 -k2,2n $tmp_bed2 > $tmp_bed3/;

# mergeBed the file.
qx/mergeBed -i $tmp_bed3 -s -c 4 -o distinct -delim ";" >  $tmp_bed4/;

# Store selected regions IDs.
my %ids2keep;

open(IN, $tmp_bed4) or die "Cannot open $tmp_bed4: $!";

while (<IN>) {
    chomp;
    my ($ids) = (split /\t/)[4];
    $ids .= ";";
    # Select best ID from cluster.
    my $best_id = "-";
    my $best_sc = 0;
    while ($ids =~ /(.+?);/g) {
        my $id = $1;
        if ($id2sc{$id} > $best_sc) {
            $best_sc = $id2sc{$id};
            $best_id = $id;
        }
    }
    $ids2keep{$best_id} = $best_sc;
}
close IN;

my $c_keep = keys %ids2keep;

print "Unique transcript matches to keep:            $c_keep\n";

unless ($c_keep) {
    die "ERROR: no best IDs selected";
}


"""





