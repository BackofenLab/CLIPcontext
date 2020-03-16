# CLIPcontext
CLIPcontext is a tool suite that offers several modes to map RBP binding regions to the transcriptome or genome. The following modes are currently available:

### G2T

In **G2T** mode, CLIPcontext takes genomic RBP binding regions or sites identified by CLIP-seq and maps them to the transcriptome. This way, region sequences are retrieved with both genomic and transcript sequence context. Depending on the location of the binding regions and the set context length, this leads to the extraction of two different sequence contexts:

<img src="doc/gen_tr_context.png" alt="Site with genomic and transcript context"
	title="Site with genomic and transcript context" width="700" />

(A) illustrates the usual way to extract CLIP-seq binding region sequences, after mapping of CLIP-seq reads to the genome and peak calling to identify the binding regions. The context sequence is obtained directly from the genome. In contrast, (B) shows the region mapped to the underlying transcript, from which CLIPcontext then takes the sequence to extract the (possibly more authentic) transcript context.

In **G2T** mode, CLIPcontext essentially takes care of the following tasks:

- Mapping of genomic RBP binding regions to underlying transcripts
- Optionally merge adjacent peak regions at exon borders (keep site with highest score)
- Extract both genomic and transcript context sets (BED regions + FASTA sequences) for comparative analysis

### T2G

In **T2G** mode, CLIPcontext maps transcript binding sites back to the genome. Again, context sequence sets (BED regions + FASTA sequences) are extracted for both genomic and transcript context. Both full and split matches (if sites overlap exon borders) are output.

### lst

In **LST** mode, CLIPcontext extracts the most prominent transcript for each gene from a given GTF file, producing a list of transcript IDs. The output transcript IDs list file can then be used as input (--tr) for the other modes. Most prominent here is defined as the transcript that is part of the GENCODE basic dataset + having the highest transcript support level (TSL). If there are two or more transcripts of one gene having the same TSL, the longest transcript will be selected. Genes or transcripts that do not meet the filtering criteria will not be output. Note that Ensembl GTF files were used for testing, and that due to possible differences in formatting, other GTF files might not work. 


### INT

In **INT** mode, CLIPcontext maps input sites to intron regions, returning only sites that overlap with intron regions. This can be useful to e.g. observe properties of intron-binding sites in a CLIP set of interest.

### EXB

In **EXB** mode, CLIPcontext extracts binding regions near exon borders from a set of input BED regions. It can be used to e.g. create an input dataset for CLIPcontext, focussing only on regions near exon borders for downstream analysis.

### EIR

In **EIR** mode, CLIPcontext creates exon + intron regions BED files for a given list of transcript IDs. Exon + intron regions are extracted for each input transcript ID.


## Installation

To install CLIPcontext, simply clone the repository and use the Python script within the folder:

```
git clone https://github.com/michauhl/CLIPcontext.git
cd CLIPcontext
python clipcontext -h
```

### Dependencies
Dependencies for CLIPcontext are as follows:

- python3 (tested with version 3.7.3)
- [bedtools](https://github.com/arq5x/bedtools2/releases)  (tested with version 2.26.0)
- [twoBitToFa](http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/twoBitToFa) executable in PATH
- [twoBitInfo](http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/twoBitInfo) executable in PATH

### Dataset requirements
CLIPcontext was implemented to work with and tested on human datasets retrieved from [Ensembl](http://www.ensembl.org/index.html). If you want to use CLIPcontext for different organisms or with datasets from different ressources, feel free to open up an issue.

Depending on the set mode, the following datasets need to be obtained for running CLIPcontext:
- A BED file (6-column format) with (genomic) RBP binding regions (e.g. eCLIP CLIPper peak regions obtained from [ENCODE](https://www.encodeproject.org/))
- A GTF file with genomic annotations from Ensembl (see [download page](http://www.ensembl.org/info/data/ftp/index.html))
- A list of transcript IDs defining the transcriptome to map to (from sequencing data or generate with **LST** mode)
- A genome .2bit file for extracting genomic sequences (for hg38 assembly click [here](https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.2bit))

### Test run

A small BED file of genomic RBP binding regions as well as a list of transcript IDs is already provided in the data/ subfolder. First we download and store all the remaining necessary datasets in the same folder:
```
cd data/
wget https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.2bit
wget ftp://ftp.ensembl.org/pub/release-98/gtf/homo_sapiens/Homo_sapiens.GRCh38.98.gtf.gz
cd ..
```
Now we can run CLIPcontext in **G2T** mode on the datasets in the data/ subfolder:
```
python clipcontext g2t --in data/SERBP1_K562_rep1_sites_chr1_hg38.bed --out test_out --tr data/GRCh38.p12.prominent_isoforms_chr1.out --gtf data/Homo_sapiens.GRCh38.98.gtf.gz --gen data/hg38.2bit
```

CLIPcontext will output an overview of the files produced at the end of the run (all stored in set --out folder). For more details see the documentation section below.

## Documentation

### Command line arguments

CLIPcontext command line arguments are grouped into required (mandatory) and optional arguments:

```
python CLIPcontext.py -h
usage: CLIPcontext.py [-h] --in str --out str --fa str --tr str --gtf str
                      --gen str [--thr float] [--min-len int] [--max-len int]
                      [--min-exon-ol float] [--seq-ext int] [--rev-filter]
                      [--merge-ext int] [--merge-all] [--gen-uniq-ids]

CLIPcontext takes genomic RBP binding regions identified by CLIP-seq, maps
them to the transcriptome, and retrieves the region sequences with both
genomic and transcript sequence context.

REQUIRED ARGUMENTS:
  --in str             Genomic regions (hg38) BED file (6-column format)
  --out str            Output results folder
  --fa str             Transcript sequences FASTA file (.fa or .fa.gz)
  --tr str             Transcript sequence IDs list file
  --gtf str            Genomic annotation (hg38) GTF file (.gtf or .gtf.gz)
  --gen str            Genome sequence (hg38) .2bit file

OPTIONAL ARGUMENTS:
  -h, --help           Print help message
  --thr float          Site score threshold for filtering -i bed_file
                       (default: None)
  --rev-filter         Reverse filtering (keep values <= threshold and prefer
                       sites with smaller values) (default: false)
  --min-len int        Minimum input site length for filtering -i bed_file
                       (default: False)
  --max-len int        Maximum input site length for filtering -i bed_file
                       (default: False)
  --min-exon-ol float  Minimum exon overlap of a site to be reported as
                       transcript hit (intersectBed -f parameter) (default:
                       0.9)
  --seq-ext int        Up- and downstream extension of centered sites for
                       context sequence extraction (default: 30)
  --merge-ext int      Extend regions mapped to transcripts by --merge-ext
                       before running mergeBed to merge overlapping regions
                       (default: 10)
  --merge-all          Merge all overlapping transcript sites extended by
                       --merge-ext (default: only merge sites overlapping at
                       exon borders) (default: False)
  --gen-uniq-ids       Generate unique column 4 IDs for -i .bed file entries
                       (default: False)
```
The required arguments comprise the necessary input files (--in genomic RBP binding regions, --fa transcript sequences, --gtf genomic annotations, --gen .2bit genomic sequences, --tr transcript IDs list) and the output folder (--out) to store all files produced during the run.
Optional arguments are provided for filtering the input .bed sites (--thr, --rev-filter, --min-len, --max-len), controlling the exon overlap amount (--min-exon-ol), the context sequence extension (--seq-ext), and the merging of nearby sites (--merge-ext, --merge-all). In case the input .bed file does not have unique column 4 IDs, new ones can be generated by --gen-uniq-ids.

### Algorithm description

CLIPcontext first reads in the transcript IDs (--tr) and extracts the corresponding sequences from the --fa FASTA file. Transcript IDs without corresponding sequence will be ignored, thus also genomic regions that map to these transcripts. 
After optionally filtering the input sites (--in) by score or length, the full-length sites are mapped to the given set of transcripts. This is mainly done to identify regions at exon borders. CLIPcontext performs mapping to the transcriptome only for genomic sites that show considerable overlap with an exonic region (>= 90%, --min-exon-ol 0.9). Identified border regions are merged if they overlap after applying --merge-ext (i.e., for each overlapping set of sites select the site with the highest score). Merging is done since for RBPs that bind to exonic regions, peaks called at adjacent exon ends often originate from the same binding event. In case all overlapping sites should be merged (not only sites at exon borders), use --merge-all. Only sites that map to a single exonic region (uniquely mapped sites) are kept.
Merged and uniquely mapped sites are then mapped again to the transcriptome, this time only taking their center positions for mapping. After center-position mapping, extension by --seq-ext is performed for both transcript sites and 
genomic sites to get both transcript and genomic context sequences. Transcript regions without full extension indicate their location near transcript ends.
A number of sanity checks are performed throughout the script. For example, extracted transcript center position nucleotides are compared with genomic center position nucleotides to check for correct sequence extraction.

At the end of the run, CLIPcontext prints an overview with short discriptions for each output file (--out test_out):

```
....

VARIOUS OUTPUT FILES
====================
Filtered genomic input sites .bed:
test_out/genomic_sites.bed
Genomic exon regions .bed for all transcripts with hits:
test_out/hit_transcript_exons.bed

TRANSCRIPT SITES (FROM MAPPING OF FULL-LENGTH GENOMIC SITES)
============================================================
Complete matches (unique + non-unique) on transcripts .bed:
test_out/transcript_hits_complete.bed
Incomplete matches (unique + non-unique) on transcripts .bed:
test_out/transcript_hits_incomplete.bed
Complete matches (unique only) on transcripts .bed:
test_out/transcript_hits_complete_unique.bed
All matches (complete + incomplete, unique only) on transcripts .bed:
test_out/transcript_hits_all_unique.bed
Mapping statistics for all transcripts with hits:
test_out/hit_transcript_stats.out

MERGED + EXTENDED TRANSCRIPT SITES (UNIQUE MATCHES ONLY)
========================================================
Unique transcript matches center positions .bed:
test_out/transcript_sites.unique_hits.cp.bed
Unique transcript matches center positions extended .bed:
test_out/transcript_sites.unique_hits.cp.ext.bed
Unique transcript matches center positions extended .fa:
test_out/transcript_sites.unique_hits.cp.ext.fa

GENOMIC SITES (CORRESPONDING TO MERGED + EXTENDED TRANSCRIPT SITES)
===================================================================
Genomic sites center positions .bed:
test_out/genomic_sites.cp.bed
Genomic sites center positions extended .bed:
test_out/genomic_sites.cp.ext.bed
Genomic sites center positions extended .fa:
test_out/genomic_sites.cp.ext.fa

```
Notice the naming conventions of the output files (cp : center-positioned site, ext : sites extended by --seq-ext). Additional mapping statistics and information for each transcript are stored in the file hit_transcript_stats.out (see below for format).


### Dataset formats

#### BED file (--in)

The --in BED file should be in 6-column format, with unique region (site) IDs in column 4:
```
head -5 SERBP1_K562_rep1_sites_chr1_hg38.bed
chr1	944348	944389	SERBP1_K562_rep01_1323	3.837688284099	-
chr1	944756	944801	SERBP1_K562_rep01_1242	4.25272578337784	-
chr1	951145	951177	SERBP1_K562_rep01_1327	3.837688284099	-
chr1	1020139	1020235	SERBP1_K562_rep01_1014	5.51576018921164	+
chr1	1040657	1040700	SERBP1_K562_rep01_1049	4.25272578337784	+
```
Note that column 5 is expected to store the region score (in this case the log2 fold change output by the CLIPper peak caller). By default, a higher score is preferred when filtering the sites by --thr. In case this should be reversed (e.g. if p-values are given), use --rev-filter in combination with --thr.

#### Transcript IDs file (--tr)

The transcript IDs list (--tr) is simply a file with one transcript ID per row:
```
head -5 GRCh38.p12.prominent_isoforms_chr1.out
ENST00000456328
ENST00000417324
ENST00000635159
ENST00000445118
ENST00000446136
```
The transcript IDs defined in this file are used to define the transcript set onto which the genomic regions (--in) are mapped to.

#### Mapping statistics output file

The mapping statistics output file stores the different hit counts and additional information for each transcript:
```
head -5 test_out/hit_transcript_stats.out 
tr_id	chr	gen_s	gen_e	pol	gene_id	gene_name	gene_biotype	tr_len	comp_hits	all_hits	uniq_comp_hits	uniq_all_hits
ENST00000361427	chr1	26472439	26476642	+	ENSG00000198830HMGN2	protein_coding	1940	1	1	1	1
ENST00000327300	chr1	32013867	32043877	+	ENSG00000121774KHDRBS1	protein_coding	2713	6	6	6	6
ENST00000368300	chr1	156114710	156140081	+	ENSG00000160789LMNA	protein_coding	3178	11	12	11	12
ENST00000379370	chr1	1020119	1056116	+	ENSG00000188157	AGRN	protein_coding	7326	11	12	11	12
```

Additional information includes the gene ID, gene name, gene biotype, transcript length, and the genomic region coordinates of the transcript. The different hit counts are: # complete (full-length matching) hits, # all hits (complete and incomplete), # unique (matching to one exon/transcript only) + complete hits, and # all unique hits.

### Data preprocessing scripts

Additional data preprocessing scripts are available in the scripts/ subfolder:

#### bed_get_regions_near_exon_borders.py
This script extracts binding regions near exon borders from a set of input BED regions. It can be used to create an input dataset for CLIPcontext, focussing only on regions near exon borders for further analysis.

#### gtf_extract_exon_regions.py
This script creates an exon regions BED file for a given list of transcript IDs. For each input transcript ID the exon regions get extracted. You will need this BED file e.g. as an input file to **bed_get_regions_near_exon_borders.py**.

#### gtf_extract_most_prominent_transcripts.py
This script extracts the most prominent transcript for each gene from a GTF file, producing a list of transcript IDs. Most prominent here is defined as the transcript that is part of the GENCODE basic dataset + having the highest transcript support level (TSL). If there are two or more transcripts of one gene having the same TSL, the longest transcript will be selected. Genes or transcripts that do not meet the filtering criteria will not be output. See the script help page for more details / extraction options. The most prominent transcript IDs output list can be used as input for CLIPcontext (--tr).

