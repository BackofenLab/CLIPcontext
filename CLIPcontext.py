#!/usr/bin/env python3


from lib import cliplib
import argparse
import shutil
import gzip
import sys
import os


"""

TOOL DEPENDENCIES
=================

python3 (tested with version 3.7.3)

bedtools (tested with version v2.26.0)
https://github.com/arq5x/bedtools2/releases

twoBitToFa
http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/twoBitToFa


REQUIRED DATASETS
=================

Currently only works with / tested on human datasets (hg38)

.2bit:
https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.2bit
Go here to get latest datasets (transcript FASTA, GTF):
http://www.ensembl.org/info/data/ftp/index.html


EXAMPLE CALLS
=============

cd data/
wget https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.2bit
wget ftp://ftp.ensembl.org/pub/release-98/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
wget ftp://ftp.ensembl.org/pub/release-98/fasta/homo_sapiens/ncrna/Homo_sapiens.GRCh38.ncrna.fa.gz
cat Homo_sapiens.GRCh38.cdna.all.fa.gz Homo_sapiens.GRCh38.ncrna.fa.gz > Homo_sapiens.GRCh38.cdna.ncrna.fa.gz
wget ftp://ftp.ensembl.org/pub/release-98/gtf/homo_sapiens/Homo_sapiens.GRCh38.98.gtf.gz
cd ..

python CLIPcontext.py --in data/SERBP1_K562_rep1_sites_chr1_hg38.bed --out test_out --tr data/GRCh38.p12.prominent_isoforms_chr1.out --fa data/Homo_sapiens.GRCh38.cdna.ncrna.fa.gz --gtf data/Homo_sapiens.GRCh38.98.gtf.gz --gen data/hg38.2bit






Add mapping count and list of transcript IDs per site ID

Supply a script that selects most prominent transcripts from .gtf file.
Option to extract most prominent isoforms from GTF file.
generate transcript IDs file script,
taking GTF and output list file (as input to CLIPcontext)



"""

################################################################################

def setup_argument_parser():
    """Setup argparse parser."""
    help_description = """
    CLIPcontext takes genomic RBP binding regions identified by CLIP-seq, 
    maps them to the transcriptome, and retrieves the region sequences 
    with both genomic and transcript sequence context.
    """
    # Define argument parser.
    p = argparse.ArgumentParser(add_help=False,
                                prog="CLIPcontext.py",
                                description=help_description,
                                formatter_class=argparse.MetavarTypeHelpFormatter)

    # Argument groups.
    p_man = p.add_argument_group("REQUIRED ARGUMENTS")
    p_opt = p.add_argument_group("OPTIONAL ARGUMENTS")

    # Required arguments.
    p_opt.add_argument("-h", "--help",
                   action="help",
                   help="Print help message")
    p_man.add_argument("--in",
                   dest="in_bed",
                   type=str,
                   required = True,
                   help = "Genomic regions (hg38) BED file (6-column format)")
    p_man.add_argument("--out",
                   dest="out_folder",
                   type=str,
                   required = True,
                   help = "Output results folder")
    p_man.add_argument("--fa",
                   dest="in_fasta",
                   type=str,
                   required = True,
                   help = "Transcript sequences FASTA file (.fa or .fa.gz)")
    p_man.add_argument("--tr",
                   dest="in_tr_list",
                   type=str,
                   required = True,
                   help = "Transcript sequence IDs list file")
    p_man.add_argument("--gtf",
                   dest="in_gtf",
                   type=str,
                   required = True,
                   help = "Genomic annotation (hg38) GTF file (.gtf or .gtf.gz)")
    p_man.add_argument("--gen",
                   dest="in_genome_2bit",
                   type=str,
                   required = True,
                   help = "Genome sequence (hg38) .2bit file")
    # Optional arguments.
    p_opt.add_argument("--thr",
                   dest = "score_thr",
                   type = float,
                   default = None,
                   help = "Site score threshold for filtering -i bed_file (default: None)")
    p_opt.add_argument("--min-len",
                   dest = "min_site_len",
                   type = int,
                   default = False,
                   help = "Minimum input site length for filtering -i bed_file (default: False)")
    p_opt.add_argument("--rev-filter",
                   dest = "rev_filter",
                   default = False,
                   action = "store_true",
                   help = "Reverse filtering (keep values <= threshold and prefer sites with smaller values) (default: false)")
    p_opt.add_argument("--max-len",
                   dest = "max_site_len",
                   type = int,
                   default = False,
                   help = "Maximum input site length for filtering -i bed_file (default: False)")
    p_opt.add_argument("--min-exon-ol",
                   dest = "min_exon_ovlp",
                   type = float,
                   default = 0.9,
                   help = "Minimum exon overlap of a site to be reported as transcript hit (intersectBed -f parameter) (default: 0.9)")
    p_opt.add_argument("--merge-ext",
                   dest="merge_ext",
                   type = int,
                   default = 10,
                   help = "Extend regions mapped to transcripts by --merge-ext before running mergeBed to merge overlapping regions (default: 10)")
    p_opt.add_argument("--merge-all",
                   dest = "merge_all",
                   default = False,
                   action = "store_true",
                   help = "Merge all overlapping transcript sites extended by --merge-ext (default: only merge sites overlapping at exon borders) (default: False)")
    p_opt.add_argument("--seq-ext",
                   dest="us_ds_ext",
                   type = int,
                   default = 30,
                   help = "Up- and downstream extension of centered sites for context sequence extraction (default: 30)")
    p_opt.add_argument("--gen-uniq-ids",
                   dest = "gen_uniq_ids",
                   default = False,
                   action = "store_true",
                   help = "Generate unique column 4 IDs for -i .bed file entries (default: False)")
    return p


################################################################################

if __name__ == '__main__':

    # Setup argparse.
    parser = setup_argument_parser()
    # Read in command line arguments.
    args = parser.parse_args()
    # Some banner.
    print("   _______   _______                __          __   ")
    print("  / ___/ /  /  _/ _ \_______  ___  / /______ __/ /_  ")
    print(" / /__/ /___/ // ___/ __/ _ \/ _ \/ __/ -_) \ / __/  ")
    print(" \___/____/___/_/   \__/\___/_//_/\__/\__/_\_\\__/   ")
    print("                                                     ")
    # Check for Linux.
    if not "linux" in sys.platform:
        print("ERROR: please use Linux")
        sys.exit()
    # Check tool availability.
    if not cliplib.is_tool("bedtools"):
        print("ERROR: bedtools not in PATH")
        sys.exit()
    if not cliplib.is_tool("twoBitToFa"):
        print("ERROR: twoBitToFa not in PATH")
        sys.exit()
    # Check file inputs.
    if not os.path.exists(args.in_bed):
        print("ERROR: input .bed file \"%s\" not found" %(args.in_bed))
    if not os.path.exists(args.in_gtf):
        print("ERROR: input .gtf file \"%s\" not found" %(args.in_gtf))
    if not os.path.exists(args.in_fasta):
        print("ERROR: input .fa file \"%s\" not found" %(args.in_fasta))
    if not os.path.exists(args.in_tr_list):
        print("ERROR: input transcript ID list file \"%s\" not found" %(args.in_tr_list))
    if not os.path.exists(args.in_genome_2bit):
        print("ERROR: input .2bit file \"%s\" not found" %(args.in_genome_2bit))
    # Check .bed for content.
    c_in_sites = cliplib.count_file_rows(args.in_bed)
    if not c_in_sites:
        print("ERROR: input .bed file \"%s\" is empty" %(args.in_bed))
        sys.exit()
    # Check .bed for 6-column format.
    if not cliplib.bed_check_six_col_format(args.in_bed):
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

    # Parameter log file.
    pars_file = args.out_folder + "/parameters.log"
    PARS = open(pars_file, "w")
    for arg in vars(args):
        PARS.write("%s\t%s\n" %(arg, str(getattr(args, arg))))
    PARS.close()

    # Output files.
    filt_in_bed = args.out_folder + "/" + "genomic_sites.bed"
    gen_cp_bed = args.out_folder + "/" + "genomic_sites.cp.bed"
    gen_cp_ext_bed = args.out_folder + "/" + "genomic_sites.cp.ext.bed"
    gen_cp_ext_fa = args.out_folder + "/" + "genomic_sites.cp.ext.fa"
    tr_uniq_cp_bed = args.out_folder + "/" + "transcript_sites.unique_hits.cp.bed"
    tr_uniq_cp_ext_bed = args.out_folder + "/" + "transcript_sites.unique_hits.cp.ext.bed"
    tr_uniq_cp_ext_fa = args.out_folder + "/" + "transcript_sites.unique_hits.cp.ext.fa"
    # Full length files.
    tr_fl_comp_bed = args.out_folder + "/" + "transcript_hits_complete.bed"
    tr_fl_incomp_bed = args.out_folder + "/" + "transcript_hits_incomplete.bed"
    tr_fl_uniq_comp_bed = args.out_folder + "/" + "transcript_hits_complete_unique.bed"
    tr_fl_uniq_all_bed = args.out_folder + "/" + "transcript_hits_all_unique.bed"
    # Full lengths hit transcript exons + stats.
    hit_tr_exons_bed = args.out_folder + "/" + "hit_transcript_exons.bed"
    hit_tr_stats_out = args.out_folder + "/" + "hit_transcript_stats.out"

    # Read in transcript ID list.
    tr_ids_dic = cliplib.read_ids_into_dic(args.in_tr_list)
    tr_ids_c = len(tr_ids_dic)
    if tr_ids_c == 0:
        print("ERROR: no transcript IDs read in from \"%s\"" %(args.in_tr_list))
        sys.exit()
    print("# transcript IDs read in:                  %i" %(tr_ids_c))
    print("# input .bed sites:                        %i" %(c_in_sites))

    # Read in transcript sequences.
    print("Reading in -f .fa sequences ...")
    tr_seqs_dic = cliplib.read_fasta_into_dic(args.in_fasta, ids_dic=tr_ids_dic)
    tr_seqs_c = len(tr_seqs_dic)
    print("# transcript sequences read in:            %i" %(tr_seqs_c))

    # Check for IDs present in list but not in FASTA file.
    lost_ids_dic = {}
    for tr_id in tr_ids_dic:
        if not tr_id in tr_seqs_dic:
            lost_ids_dic[tr_id] = 1
            print("WARNING: transcript ID \"%s\" missing in -f .fa file" %(tr_id))
    if lost_ids_dic:
        print("WARNING: sites that map on missing transcripts will be skipped!")
    # Remove lost IDs from transcript IDs dic.
    for seq_id in lost_ids_dic:
        del tr_ids_dic[seq_id]

    # Filter sites by threshold.
    if args.score_thr is not None:
        cliplib.bed_process_bed_file(args.in_bed, filt_in_bed,
                                     score_thr=args.score_thr,
                                     min_len=args.min_len,
                                     max_len=args.max_len,
                                     generate_unique_ids=args.gen_uniq_ids,
                                     rev_filter=args.rev_filter)
    else:
        cliplib.make_file_copy(args.in_bed, filt_in_bed)

    # Number of remaining sites.
    c_filt_sites = cliplib.count_file_rows(filt_in_bed)
    print("# input .bed sites after --thr filtering:  %i" %(c_filt_sites))
    if not c_filt_sites:
        print("ERROR: no remaining BED regions after --thr filtering")
        sys.exit()

    """
    First map full length filtered input .bed to transcriptome.
    Do this to get regions at exon borders, which will later be merged.
    """
    print("Mapping full-length genomic input sites to transcriptome ... ")
    full_length_out = args.out_folder + "/" + "transcript_map_full_length_out"
    cliplib.convert_genome_positions_to_transcriptome(filt_in_bed, full_length_out, 
                                                      args.in_gtf, tr_ids_dic,
                                                      intersectBed_f=args.min_exon_ovlp)
    # All unique transcript hits (complete + incomplete hits).
    fl_all_uniq_tr_hits_bed = full_length_out + "/" + "transcript_hits_all_unique.bed"
    # Complete and incomplete hits, both unique and multi hits.
    fl_comp_transcript_hits_bed = full_length_out + "/" + "transcript_hits_complete.bed"
    fl_incomp_transcript_hits_bed = full_length_out + "/" + "transcript_hits_incomplete.bed"
    # Unique hits, only complete matches.
    fl_uniq_comp_hits_bed = full_length_out + "/" + "transcript_hits_complete_unique.bed"

    # Genomic exon .bed of transcripts provided by -t.
    genome_exon_bed = full_length_out + "/" + "genomic_exon_coordinates.bed"
    # Exons .bed for all transcripts with full length hits.
    fl_hit_tr_exons_bed = full_length_out + "/" + "hit_transcript_exons.bed"
    # Stats for transcripts with full length hits.
    fl_hit_tr_stats_out = full_length_out + "/" + "hit_transcript_stats.out"
    # Copy to output files.
    cliplib.make_file_copy(fl_hit_tr_exons_bed, hit_tr_exons_bed)
    cliplib.make_file_copy(fl_hit_tr_stats_out, hit_tr_stats_out)
    cliplib.make_file_copy(fl_comp_transcript_hits_bed, tr_fl_comp_bed)
    cliplib.make_file_copy(fl_incomp_transcript_hits_bed, tr_fl_incomp_bed)
    cliplib.make_file_copy(fl_all_uniq_tr_hits_bed, tr_fl_uniq_all_bed)
    cliplib.make_file_copy(fl_uniq_comp_hits_bed, tr_fl_uniq_comp_bed)

    # Count unique transcript hits.
    c_uniq_fl_hits = cliplib.count_file_rows(fl_all_uniq_tr_hits_bed)
    print("# unique transcript hits:                  %i" %(c_uniq_fl_hits))
    if not c_uniq_fl_hits:
        print("ERROR: no unique transcript hits for given genomic .bed and transcripts")
        sys.exit()

    # Get IDs for unique transcript hits (site ID -> site score dic).
    uniq_tr_hit_dic = cliplib.bed_get_region_id_scores(fl_all_uniq_tr_hits_bed)
    # Map site IDs to sequence IDs.
    uniq_tr_seqid2siteid_dic = cliplib.bed_map_region_id_to_seq_id(fl_all_uniq_tr_hits_bed)

    # Get IDs for genomic sites that completely overlap with exons (after merge_ext).
    tmp_bed1 = args.out_folder + "/" + "genomic_sites.merge_ext.tmp.bed"
    cliplib.bed_process_bed_file(filt_in_bed, tmp_bed1,
                                 ids2keep_dic=uniq_tr_hit_dic,
                                 ext_lr=args.merge_ext)
    tmp_bed2 = args.out_folder + "/" + "genomic_sites.merge_ext.comp_exon_ovlp.tmp.bed"
    cliplib.intersect_bed_files(tmp_bed1, genome_exon_bed, "-s -f 1 -u", tmp_bed2)
    comp_ex_hit_dic = cliplib.bed_get_region_id_scores(tmp_bed2)

    # Get unique transcript sites that only partially overlap with exons (after merge_ext).
    uniq_part_tr_hit_dic = {}
    uniq_comp_tr_hit_dic = {}
    for site_id in uniq_tr_hit_dic:
        if site_id in comp_ex_hit_dic: # if completely overlapping.
            uniq_comp_tr_hit_dic[site_id] = uniq_tr_hit_dic[site_id]
        else:
            uniq_part_tr_hit_dic[site_id] = uniq_tr_hit_dic[site_id]
    fl_uniq_comp_c = len(uniq_comp_tr_hit_dic)
    fl_uniq_incomp_c = len(uniq_part_tr_hit_dic)
    print("# unique complete transcript hits:         %i" %(fl_uniq_comp_c))
    print("# unique incomplete transcript hits:       %i" %(fl_uniq_incomp_c))

    """
    Merge overlapping sites, keep only highest-scoring site for each set of 
    overlapping regions.
    By default, merging takes place only for sites near exon borders 
    (more precisely, sites that do not fully overlap with exons after 
    extending them by --merge-ext). However, if --merge-all is set, 
    merging will be done for all overlapping sites.
    """
    # Temp files for merging operations.
    tmp_bed3 = args.out_folder + "/" + "prolonged_unique_transcript_hits.tmp.bed";
    tmp_bed4 = args.out_folder + "/" + "prolonged_unique_transcript_hits.sorted.tmp.bed";
    tmp_bed5 = args.out_folder + "/" + "prolonged_unique_transcript_hits.merged.tmp.bed";

    # IDs to keep after merging.
    ids2keep_dic = {}
    # If all overlapping sets of sites should be merged.
    if args.merge_all:
        print("Merging overlapping sites ... ")
        # Prolong unique hits on transcripts by set merge extension.
        cliplib.bed_process_bed_file(fl_all_uniq_tr_hits_bed, tmp_bed3,
                                     zero_scores=True,
                                     ext_lr=args.merge_ext)
        # Sort .bed file.
        cliplib.bed_sort_file(tmp_bed3, tmp_bed4)
        # Merge .bed file (mergeBed overlapping regions).
        cliplib.bed_merge_file(tmp_bed4, tmp_bed5)
        # Select region IDs (in case of overlaps choose highest-scoring region).
        ids2keep_dic = cliplib.bed_merge_file_select_top_ids(tmp_bed5, uniq_tr_hit_dic)
    else:
        print("Merging overlapping sites at exon borders ... ")
        if uniq_part_tr_hit_dic:
            # Do merging only for sites near exon borders (unique hits partially overlapping with exons).
            cliplib.bed_process_bed_file(fl_all_uniq_tr_hits_bed, tmp_bed3,
                                         ids2keep_dic=uniq_part_tr_hit_dic,
                                         ext_lr=args.merge_ext)
            # Sort .bed file.
            cliplib.bed_sort_file(tmp_bed3, tmp_bed4)
            # Merge .bed file (mergeBed overlapping regions).
            cliplib.bed_merge_file(tmp_bed4, tmp_bed5)
            # Select region IDs (in case of overlaps choose highest-scoring region).
            ids2keep_dic = cliplib.bed_merge_file_select_top_ids(tmp_bed5, uniq_part_tr_hit_dic)
        # Add complete exon hits to ids2keep dictionary.
        ids2keep_dic.update(uniq_comp_tr_hit_dic)
    # Number of transcript sites after merging.
    c_ids2keep = len(ids2keep_dic)
    print("# IDs after merging overlapping sites      %i" %(c_ids2keep))

    """
    Now map center position genomic sites to transcriptomes.
    Use only merged sites that uniquely mapped in the first step,
    with their IDs stored in ids2keep_dic. After center-position mapping,
    extension by --seq-ext will be done for both transcript sites and 
    genomic sites to get both transcript and genomic context sites.
    """

    print("Mapping center positions of genomic input sites to transcriptome ... ")

    # Convert filtered input .bed to genomic center positions .bed.
    cliplib.bed_process_bed_file(filt_in_bed, gen_cp_bed,
                                 ids2keep_dic=ids2keep_dic,
                                 center_sites=True)
    # Make a second .bed file without "0" scores, for twoBitToFa compatibility.
    tmp_bed6 = args.out_folder + "/" + "genomic_sites.cp.no_sc.tmp.bed"
    cliplib.bed_process_bed_file(filt_in_bed, tmp_bed6,
                                 ids2keep_dic=ids2keep_dic,
                                 zero_scores=True,
                                 center_sites=True)

    # Now map center position regions to transcriptome.
    center_pos_out = args.out_folder + "/" + "transcript_center_pos_out"
    cliplib.convert_genome_positions_to_transcriptome(gen_cp_bed, center_pos_out, 
                                                      args.in_gtf, tr_ids_dic)
    # All unique hits (complete + incomplete).
    cp_all_uniq_tr_hits_bed = center_pos_out + "/" + "transcript_hits_all_unique.bed"
    cp_comp_tr_hits_bed = center_pos_out + "/" + "transcript_hits_complete.bed"

    """
    Get site ID / sequence ID combinations to keep.
    We need these since mapping center positions to transcriptome can result 
    in sites which mapped uniquely using their full lengths before, but with 
    center positions could map to more than one exon / transcript. This 
    is due to the intersectBed_f overlap parameter set to 0.9, which does not 
    report hits if they do not overlap 90%+ with exons, while for center 
    position mapping this is always satisfied.
    """
    # Get site ID / sequence ID combinations to keep.
    ids2keep2_dic = {}
    for site_id in uniq_tr_seqid2siteid_dic:
        seq_id = uniq_tr_seqid2siteid_dic[site_id]
        if site_id in ids2keep_dic:
            ids2keep2_dic[site_id] = seq_id
    # Filter center position complete hits by site ID -> sequence ID combinations.
    cliplib.bed_process_bed_file(cp_comp_tr_hits_bed, tr_uniq_cp_bed,
                                 ids2keep_dic=ids2keep2_dic,
                                 seq_id_filter=True)

    # Transcript sequence lengths dictionary.
    tr_seq_len_dic = cliplib.get_seq_lengths_from_seqs_dic(tr_seqs_dic)

    # Extend center position transcript sites.
    cliplib.bed_process_bed_file(tr_uniq_cp_bed, tr_uniq_cp_ext_bed,
                                 seq_len_dic=tr_seq_len_dic,
                                 ext_lr=args.us_ds_ext)

    # Get .bed rows.
    id2row_dic = cliplib.bed_read_rows_into_dic(tr_uniq_cp_ext_bed)
    id2cprow_dic = cliplib.bed_read_rows_into_dic(tr_uniq_cp_bed)

    # Count unique hits.
    c_uniq_cp_hits = len(id2row_dic)
    print("# unique transcript hits (center pos)      %i" %(c_uniq_cp_hits))
    if not c_uniq_cp_hits:
        print("ERROR: no unique transcript hits for given genomic .bed (center positions) and transcripts")
        sys.exit()

    # Get transcript sequences for sites (site ID -> transcript sequence).
    id2seq_dic = cliplib.extract_transcript_sequences(id2row_dic, tr_seqs_dic,
                                                      full_hits_only=False)
    # Count # of full-length transcript sequences.
    fl_l = 1 + args.us_ds_ext*2
    c_fl = 0
    c_all = len(id2seq_dic)
    for site_id in id2seq_dic:
        seq_l = len(id2seq_dic[site_id])
        if seq_l == fl_l:
            c_fl += 1
    print("# extracted transcript region sequences    %i" %(c_all))
    print("# full-length transcript region sequences  %i" %(c_fl))

    # Output sequences to FASTA file.
    cliplib.output_seqs_dic_to_fasta(id2seq_dic, tr_uniq_cp_ext_fa)

    """
    Next, do a sanity check and extract genomic context sequences.
    1) Sanity check:
       Compare extracted transcript center position nucleotides with
       with genomic center position nucleotides.
    2) Extract genomic sequences for given input regions,
       using up- and downstream extension used for transcript sequences.
    """

    print("Checking transcript against genomic center position nucleotides ... ")

    # Get center position nucleotides from transcript regions.
    id2tcp_dic = cliplib.extract_transcript_sequences(id2cprow_dic, tr_seqs_dic,
                                                      ext_lr=False,
                                                      full_hits_only=False)
    # Get center position nucleotides from genome.
    tmp_fa1 = args.out_folder + "/" + "genomic_sites.cp.tmp.fa";
    cliplib.bed_extract_sequences_from_2bit(tmp_bed6, tmp_fa1, args.in_genome_2bit)
    id2gcp_dic = cliplib.read_fasta_into_dic(tmp_fa1)
    # Compare center position nucleotides.
    for site_id in id2tcp_dic:
        tcp = id2tcp_dic[site_id]
        assert site_id in id2gcp_dic, "site ID \"%s\" not in genomic center positions dictionary" %(site_id)
        gcp = id2gcp_dic[site_id]
        l_tcp = len(tcp)
        l_gcp = len(gcp)
        assert l_tcp == 1, "Transcript center position sequence length != 1 for site ID \"%s\" (%i != 1)" %(site_id, l_tcp)
        assert l_gcp == 1, "Genomic center position sequence length != 1 for site ID \"%s\" (%i != 1)" %(site_id, l_gcp)
        if not tcp == gcp:
            print("WARNING: differing center position nucleotides for \"%s\" (%s != %s)" % (site_id, tcp, gcp))

    # Center genomic sites.
    cliplib.bed_process_bed_file(filt_in_bed, gen_cp_bed,
                                 center_sites=True)
    # Extend genomic center position sites.
    cliplib.bed_process_bed_file(filt_in_bed, gen_cp_ext_bed,
                                 ext_lr=args.us_ds_ext,
                                 center_sites=True)
    # Extend genomic center position sites with zero scores for twoBitToFa.
    tmp_bed7 = args.out_folder + "/" + "genomic_sites.cp.ext.no_sc.tmp.bed"
    cliplib.bed_process_bed_file(filt_in_bed, tmp_bed7,
                                 ext_lr=args.us_ds_ext,
                                 zero_scores=True,
                                 center_sites=True)
    # Get extended genomic regions.
    cliplib.bed_extract_sequences_from_2bit(tmp_bed7, gen_cp_ext_fa, args.in_genome_2bit)

    # Clean up.
    clean_up = True
    if clean_up:
        print("Removing temporary files ... ")
        # Remove tmp files.
        if os.path.exists(tmp_bed1):
            os.remove(tmp_bed1)
        if os.path.exists(tmp_bed2):
            os.remove(tmp_bed2)
        if os.path.exists(tmp_bed3):
            os.remove(tmp_bed3)
        if os.path.exists(tmp_bed4):
            os.remove(tmp_bed4)
        if os.path.exists(tmp_bed5):
            os.remove(tmp_bed5)
        if os.path.exists(tmp_bed6):
            os.remove(tmp_bed6)
        if os.path.exists(tmp_bed7):
            os.remove(tmp_bed7)
        if os.path.exists(tmp_fa1):
            os.remove(tmp_fa1)
        # Remove mapping folders.
        if os.path.exists(full_length_out):
            shutil.rmtree(full_length_out, ignore_errors=True)
        if os.path.exists(center_pos_out):
            shutil.rmtree(center_pos_out, ignore_errors=True)

    # Output files report.
    print("")
    print("GENOMIC FILES")
    print("=============")
    print("Filtered genomic input sites .bed:\n%s" %(filt_in_bed))
    print("Filtered genomic input sites center positions .bed:\n%s" %(gen_cp_bed))
    print("Filtered genomic input sites extended .bed:\n%s" %(gen_cp_ext_bed))
    print("Filtered genomic input sites extended .fa:\n%s" %(gen_cp_ext_fa))
    print("Genomic exon regions .bed for all transcripts with hits:\n%s\n" %(hit_tr_exons_bed))

    print("TRANSCRIPT SITES (FROM MAPPING OF FULL-LENGTH GENOMIC SITES)")
    print("============================================================")
    print("Complete matches on transcripts .bed:\n%s" %(tr_fl_comp_bed))
    print("Incomplete matches on transcripts .bed:\n%s" %(tr_fl_incomp_bed))
    print("Unique complete matches on transcripts .bed:\n%s" %(tr_fl_uniq_comp_bed))
    print("All unique (complete + incomplete) matches on transcripts .bed:\n%s" %(tr_fl_uniq_comp_bed))
    print("Mapping statistics for all transcripts with hits:\n%s\n" %(hit_tr_stats_out))

    print("MERGED + EXTENDED TRANSCRIPT SITES")
    print("==================================")
    print("Up- and downstream context sequence extension: %i" %(args.us_ds_ext))
    print("Unique matches on transcripts center positions .bed:\n%s" %(tr_uniq_cp_bed))
    print("Unique matches on transcripts center positions extended .bed:\n%s" %(tr_uniq_cp_ext_bed))
    print("Unique matches on transcripts center positions extended .fa:\n%s\n" %(tr_uniq_cp_ext_fa))



