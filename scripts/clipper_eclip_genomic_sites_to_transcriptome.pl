#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;


=head1 NAME

=head1 SYNOPSIS

This script maps genomic sites identified by the CLIPper peak caller 
(input has to be a 6 column .bed file containing the peak regions with 
their log2 fold change values in column 5) to transcriptome. The script 
takes care of several things, e.g. merging peaks called at 
adjacent exon ends, and also extracts extended sequences from transcript 
context (using additional scripts stored in subfolder scripts/).



care of several 
things (e.g. adjacent peaks called on two neighbor exons). It also 
extracts extended sequences from transcriptomic context (using 
array of existing scripts, so scripts folder needed). A two-step 
approach is used to remove double peaks called by CLIPper, first 
mapping full binding sites 


Arguments:

    -help|h          display help page
    -bed             Genomic eCLIP peak regions input .bed file
                     with log2 fold change values in column 5
    -data-id         Data ID for naming output files
    -out             Results output folder
    -tr-list         Transcript IDs list file, one ID per line
                     Genomic regions will be mapped on these 
                     transcripts
    -gtf             GTF file with exon information for mapping 
                     genomic sites to transcripts specified 
                     in -tr-list file
    -fa              FASTA file with transcript sequences for 
                     -tr-list IDs
    -thr             Filter out -in regions < -thr column 5 score
                     For CLIPper this sets the minimum log2 fold 
                     change a site has to have to be kept
                     default: 1
    -merge-ext       Merge extension, extend regions mapped to 
                     transcripts by -merge-ext before running 
                     mergeBed to merge overlapping regions
                     default: 20
    -vp-ext          Viewpoint extension of center position on 
                     transcripts
                     default: 30 (i.e. 61 nt viewpoint)
    -con-ext         Context extension (lower case characters)
                     default: 150 (i.e. with 61 nt vp = 361 total)
    -genome-2bit     Genome .2bit file to compare center positions 
                     extracted from transcripts with positions 
                     from genome as sanity check
                     If set this also requires twoBitToFa in path
                     default: false

=head1 DESCRIPTION

Tools needed:

Tool:    bedtools intersect (aka intersectBed)
Tool:    bedtools merge (aka mergeBed)
Version: v2.26.0
twoBitToFa (if -genome-2bit is set)


From scripts repo:
bed_extract_sequences_from_fasta.pl
convert_genome_positions_to_transcriptome.pl

clipper_eclip_genomic_sites_to_transcriptome.pl 
-fa /home/uhlm/Data/ensembl_data/GRCh38.97_cdna_ncrna.new_header.fa
-bed bed_minfold1_rmncc_uniqcol4_hg38_merged_rep12_near_exons_maxd50_out/FXR1_K562.bed
-data-id FXR1_K562_near_exbo_maxd50
-out FXR1_K562_near_exbo_maxd50_out
-gtf /home/uhlm/Data/ensembl_data/Homo_sapiens.GRCh38.97.gtf.gz
-tr-list GRCh38.p12.prominent_isoforms.tp.exons.enst_list.out
-genome-2bit /home/uhlm/Data/genome_2bit/hg38.2bit

# Normal call.
perl clipper_eclip_genomic_sites_to_transcriptome.pl -fa /home/uhlm/Data/ensembl_data/GRCh38.97_cdna_ncrna.new_header.fa -bed bed_minfold1_rmncc_uniqcol4_hg38_merged_rep12_near_exons_maxd50_out/FXR1_K562.bed -data-id FXR1_K562_near_exbo_maxd50 -out FXR1_K562_near_exbo_maxd50_out -gtf /home/uhlm/Data/ensembl_data/Homo_sapiens.GRCh38.97.gtf.gz -tr-list GRCh38.p12.prominent_isoforms.tp.exons.enst_list.out

# Plus sanity check extracted sequences.
perl clipper_eclip_genomic_sites_to_transcriptome.pl -fa /home/uhlm/Data/ensembl_data/GRCh38.97_cdna_ncrna.new_header.fa -bed bed_minfold1_rmncc_uniqcol4_hg38_merged_rep12_near_exons_maxd50_out/FXR1_K562.bed -data-id FXR1_K562_near_exbo_maxd50 -out FXR1_K562_near_exbo_maxd50_out -gtf /home/uhlm/Data/ensembl_data/Homo_sapiens.GRCh38.97.gtf.gz -tr-list GRCh38.p12.prominent_isoforms.tp.exons.enst_list.out -genome-2bit /home/uhlm/Data/genome_2bit/hg38.2bit


=head1 EXAMPLES

=cut

my ($i_help, $i_bed, $i_gtf, $i_tr_list, $i_out, $i_thr, $i_data_id, $i_vp_ext, $i_con_ext, $i_merge_ext, $i_fa, $i_genome_2bit);

GetOptions ( "help|h" => \$i_help,
             "bed:s" => \$i_bed,
             "out:s" => \$i_out,
             "gtf:s" => \$i_gtf,
             "fa:s" => \$i_fa,
             "genome-2bit:s" => \$i_genome_2bit,
             "merge-ext:i" => \$i_merge_ext,
             "vp-ext:i" => \$i_vp_ext,
             "con-ext:i" => \$i_con_ext,
             "data-id:s" => \$i_data_id,
             "tr-list:s" => \$i_tr_list,
             "thr:f" => \$i_thr ) or pod2usage(1);

# Check input.
pod2usage(1) if $i_help;
($i_bed and $i_out and $i_gtf and $i_tr_list and $i_data_id and $i_fa) or pod2usage("ERROR: -bed, -out, -gtf, -tr-list, -data-id, -fa required.");

# Check for programs in PATH.
my $ch1 = qx/which twoBitToFa/;
unless ($ch1) {
    die "ERROR: twoBitToFa not in PATH\n";
}
my $ch2 = qx/which intersectBed/;
unless ($ch2) {
    die "ERROR: intersectBed not in PATH\n";
}
my $ch3 = qx/which mergeBed/;
unless ($ch3) {
    die "ERROR: intersectBed not in PATH\n";
}


# Input.
my $con_ext = 150;
my $vp_ext = 30;
my $merge_ext = 20;
my $thr = 1;
if ($i_con_ext) {
    $con_ext = $i_con_ext;
}
if ($i_vp_ext) {
    $vp_ext = $i_vp_ext;
}
if ($i_merge_ext) {
    $merge_ext = $i_merge_ext;
}
if ($i_thr) {
    $thr = $i_thr;
}

unless (-d $i_out) {
    qx/mkdir $i_out/;
}

# Output files.
my $out_map1 = $i_out . "/" . "full_length_transcriptome_map_out";
my $out_map2 = $i_out . "/" . "centerpos_transcriptome_map_out";
my $out_tr_seqs = $i_out . "/" . "transcript_sites_fasta_out";
my $out_tr_cp_seqs = $i_out . "/" . "transcript_sites_centerpos_fasta_out";
my $ignore_list_file = $i_out . "/" . "ignore_ids.out";

# Create empty ignore file.
qx/touch $ignore_list_file/;

# Check IDs.
my $c_ids = qx/cut -f 4 $i_bed | sort | uniq -d/;
if ($c_ids) {
    die "ERROR: non-unique col4 IDs in -in $i_bed file";
}

# Check whether given transcript ids are present in .fa file.
my %tr_ids;
open(IN, $i_tr_list) or die "Cannot open $i_tr_list: $!";
while(<IN>) {
    if ($_ =~ /(.+)\n/) {
        $tr_ids{$1}++;
    }
}
close IN;
my %fa_ids;
open(IN, $i_fa) or die "Cannot open $i_fa: $!";
while(<IN>) {
    if ($_ =~ />(.+)\n/) {
        $fa_ids{$1}++;
    }
}
close IN;
my %lost_tr_ids;
foreach (keys %tr_ids) {
    my $tr_id = $_;
    unless(exists $fa_ids{$tr_id}) {
        print "WARNING: transcript id \"$tr_id\" not found in given -fa file\n";
        $lost_tr_ids{$tr_id}++;
    }
}
if (keys %lost_tr_ids) {
    print "WARNING: sites on missing transcripts will be skipped!\n";
    open(OUT, '>', $ignore_list_file) or die "Cannot open $ignore_list_file: $!";
    foreach (keys %lost_tr_ids) {
        my $id = $_;
        print OUT "$id\n";
    }
    close OUT;
}

my $tmp_bed1 = $i_out . "/" . "filtered_sites.tmp.bed";
open(OUT, '>', $tmp_bed1) or die "Cannot open $tmp_bed1: $!";
open(IN, $i_bed) or die "Cannot open $i_bed: $!";
my $c_in_bed = 0;
my $c_filt_in_bed = 0;

while(<IN>) {
    my ($sc) = (split /\t/)[4];
    $c_in_bed++;
    next if ($sc < $thr);
    $c_filt_in_bed++;
    print OUT $_;
}
close IN;
close OUT;

print "Read in -bed regions:                         $c_in_bed\n";
print "Remaining -bed regions after -thr filtering:  $c_filt_in_bed\n";
unless ($c_filt_in_bed) {
    die "ERROR: no remaining sites after filtering";
}
if ($c_filt_in_bed < 4000) {
    print "WARNING: less than 4000 sites remaining\n";
}

print "Mapping full-length regions to transcriptome \n";

# First map full length -bed to transcriptome.
qx/convert_genome_positions_to_transcriptome.pl -gtf $i_gtf -bed $tmp_bed1 -out $out_map1 -data-id $i_data_id -transcript-list $i_tr_list -ignore-list $ignore_list_file -merge-uniq/;

# Unique transcript region results.
my $uniq_matches_map1 = $out_map1 . "/" . $i_data_id . "_transcript_unique_matches.bed";

# Prolong unique matches.
my $tmp_bed2 = $i_out . "/" . "prolonged_unique_transcript_matches.tmp.bed";
my $tmp_bed3 = $i_out . "/" . "prolonged_unique_transcript_matches.sorted.tmp.bed";
my $tmp_bed4 = $i_out . "/" . "prolonged_unique_transcript_matches.merged.tmp.bed";
open(OUT, '>', $tmp_bed2) or die "Cannot open $tmp_bed2: $!";
open(IN, $uniq_matches_map1) or die "Cannot open $uniq_matches_map1: $!";
my $c_uniq_m1 = 0;
my %id2sc;

while(<IN>) {
    chomp;
    my ($chr, $s, $e, $id, $sc, $pol) = (split /\t/)[0,1,2,3,4,5];
    my $new_s = $s - $merge_ext;
    my $new_e = $e + $merge_ext;
    if ($new_s < 0) {
        $new_s = 0;
    }
    $id2sc{$id} = $sc;
    $c_uniq_m1++;
    print OUT "$chr\t$new_s\t$new_e\t$id\t$sc\t$pol\n";
}
close IN;
close OUT;

print "Unique transcript matches full-length input:  $c_uniq_m1\n";

unless ($c_uniq_m1) {
    die "ERROR: no unique matches on full-length inputs";
}

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

# Now map full length -bed to transcriptome.
my $tmp_bed5 = $i_out . "/" . "genomic_center_positions.tmp.bed";
my $tmp_bed6 = $i_out . "/" . "genomic_center_positions.no_sc.tmp.bed";

open(OUT, '>', $tmp_bed5) or die "Cannot open $tmp_bed5: $!";
open(OUT2, '>', $tmp_bed6) or die "Cannot open $tmp_bed6: $!";
open(IN, $i_bed) or die "Cannot open $i_bed: $!";

while(<IN>) {
    chomp;
    my ($chr, $s, $e, $id, $sc, $pol) = (split /\t/)[0,1,2,3,4,5];
    next unless (exists $ids2keep{$id});
    my $new_start = $s;
    my $new_end = $e;
    unless ($e - $s == 1) {
        $new_end = sprintf("%.0f", ((($e - $s) / 2) + $s)) + 1;
        $new_start = $new_end - 1;
    }
    print OUT "$chr\t$new_start\t$new_end\t$id\t$sc\t$pol\n";
    print OUT2 "$chr\t$new_start\t$new_end\t$id\t0\t$pol\n";
}
close IN;
close OUT;
close OUT2;

print "Mapping center positions regions to transcriptome ... \n";

# First center position regions -bed to transcriptome.
qx/convert_genome_positions_to_transcriptome.pl -gtf $i_gtf -bed $tmp_bed5 -out $out_map2 -data-id $i_data_id -transcript-list $i_tr_list -ignore-list $ignore_list_file -merge-uniq/;

# Unique transcript region results.
my $uniq_matches_map2 = $out_map2 . "/" . $i_data_id . "_transcript_unique_matches.bed";

my $c_uniq_m2 = qx/cat $uniq_matches_map2 | wc -l/;
chomp $c_uniq_m2;

print "Unique transcript matches center positions:   $c_uniq_m2\n";
unless ($c_keep == $c_uniq_m2) {
    die "ERROR: c_keep != c_uniq_m2 ($c_keep != $c_uniq_m2)";
}

print "Extracting transcript region sequences ... \n";

# Now extract vp, context sequences from transcript center pos.
qx/bed_extract_sequences_from_fasta.pl -i $uniq_matches_map2 -o $out_tr_seqs -pos-ext $vp_ext -con-ext $con_ext -fasta $i_fa/;

# Output files.
my $tr_sites_fa = $out_tr_seqs . "/" . "extracted_sites.fa";
my $tr_sites_bed = $out_tr_seqs . "/" . "extracted_sites.bed";
my $tr_compl_sites_fa = $out_tr_seqs . "/" . "complete_sites.fa";
my $tr_compl_sites_bed = $out_tr_seqs . "/" . "complete_sites.bed";

unless (-e $tr_sites_fa) {
    die "ERROR: $tr_sites_fa not found";
}
unless (-e $tr_sites_bed) {
    die "ERROR: $tr_sites_bed not found";
}
unless (-e $tr_compl_sites_fa) {
    die "ERROR: $tr_compl_sites_fa not found";
}
unless (-e $tr_compl_sites_bed) {
    die "ERROR: $tr_compl_sites_bed not found";
}

my $c_compl_sites = qx/cat $tr_compl_sites_bed | wc -l/;
chomp $c_compl_sites;
my $c_all_sites = qx/cat $tr_sites_bed | wc -l/;
chomp $c_all_sites;

print "# all extracted transcript region sequences:  $c_all_sites\n";
print "# complete transcript region sequences:       $c_compl_sites\n";


# .2bit comparison.
if ($i_genome_2bit) {

    print "Extracting center position nucleotides ... \n";

    # Extract center position transcript nucleotides.
    qx/bed_extract_sequences_from_fasta.pl -i $uniq_matches_map2 -o $out_tr_cp_seqs -no-pos-ext -no-context -fasta $i_fa/;

    # Center position transcript FASTA.
    my $tr_cp_sites_fa = $out_tr_cp_seqs . "/" . "extracted_sites.fa";
    # Center position genomic FASTA.
    my $gen_cp_sites_fa = $i_out . "/" . "genomic_cp_sites.tmp.fa";

    # Extract center position genome nucleotides.
    qx/twoBitToFa -noMask -bed=$tmp_bed6 $i_genome_2bit $gen_cp_sites_fa/;


    # Read in extracted center position sequences.
    my $id;
    my %tr_cp_seqs;
    my %gen_cp_seqs;

    open(IN, $tr_cp_sites_fa) or die "Cannot open $tr_cp_sites_fa: $!";
    while(<IN>) {
        if ($_ =~ />(.+)\n/) {
            $id = $1;
        } elsif ($_ =~ /([ACGTUN]+)\n/i) {
            $tr_cp_seqs{$id} .= $1;
        }
    }
    close IN;

    open(IN, $gen_cp_sites_fa) or die "Cannot open $gen_cp_sites_fa: $!";
    while(<IN>) {
        if ($_ =~ />(.+)\n/) {
            $id = $1;
        } elsif ($_ =~ /([ACGTUN]+)\n/i) {
            $gen_cp_seqs{$id} .= $1;
        }
    }
    close IN;

    my $c_tr_cp_seqs = keys %tr_cp_seqs;
    my $c_gen_cp_seqs = keys %gen_cp_seqs;

    unless($c_tr_cp_seqs) {
        die "ERROR: no sequences extracted from $uniq_matches_map2";
    }
    unless($c_gen_cp_seqs) {
        die "ERROR: no sequences extracted from $tmp_bed5";
    }

    print "# center position genomic sequences:          $c_gen_cp_seqs\n";
    print "# center position transcript sequences:       $c_tr_cp_seqs\n";

    unless ($c_gen_cp_seqs == $c_tr_cp_seqs) {
        die "ERROR: number of center position sequences differs ($c_gen_cp_seqs != $c_tr_cp_seqs)";
    }

    # Compare the two sequences.
    foreach (keys %gen_cp_seqs) {
        my $id = $_;
        my $g_seq = $gen_cp_seqs{$id};
        my $l_g_seq = length($g_seq);
        unless (exists $tr_cp_seqs{$id}) {
            die "ERROR: missing gen_cp_seqs sequence ID \"$id\" in tr_cp_seqs";
        }
        my $tr_seq = $tr_cp_seqs{$id};
        my $l_tr_seq = length($tr_seq);
        unless ($l_g_seq == 1) {
            die "ERROR: l_g_seq != 1 ($l_g_seq)";
        }
        unless ($l_tr_seq == 1) {
            die "ERROR: l_tr_seq != 1 ($l_tr_seq, $tr_seq, $id)";
        }
        unless ($tr_seq eq $g_seq) {
            print "WARNING: differing center position nucleotides for \"$id\" ($tr_seq != $g_seq)\n";
        }
    }
}

# All exons of transcripts with mapped sits.
my $mapped_tr_exons_bed = $out_map2 . "/" . $i_data_id . "_mapped_transcript_exons.bed";
# Bundle files.
my $out_gen_cp_bed = $i_out . "/" . $i_data_id . "_hg38.centerpos.bed";
my $out_tr_ext_bed = $i_out . "/" . $i_data_id . "_transcript_sites.bed";
my $out_tr_ext_fa = $i_out . "/" . $i_data_id . "_transcript_sites.fa";
my $out_tr_exons_bed = $i_out . "/" . $i_data_id . "_mapped_transcript_exons.bed";
qx/cat $tmp_bed5 > $out_gen_cp_bed/;
qx/cat $tr_sites_bed > $out_tr_ext_bed/;
qx/cat $tr_sites_fa > $out_tr_ext_fa/;
qx/cat $mapped_tr_exons_bed > $out_tr_exons_bed/;

print "Output center position genomic sites .bed:\n$out_gen_cp_bed\n";
print "Output mapped transcript exons genomic .bed:\n$out_tr_exons_bed\n";
print "Output extended transcript sites .bed:\n$out_tr_ext_bed\n";
print "Output extended transcript sites .fa:\n$out_tr_ext_fa\n";




