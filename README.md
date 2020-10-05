# CLIPcontext

CLIPcontext is a tool suite that offers several modes to map RBP binding regions to the transcriptome or genome. The following modes are currently available:

### G2T

In **G2T** mode, CLIPcontext takes genomic RBP binding regions or sites identified by CLIP-seq and maps them to the transcriptome. This way, region sequences are retrieved with both genomic and transcript sequence context. Depending on the location of the binding regions and the set context length, this leads to the extraction of two different sequence contexts:

<img src="docs/img/gen_tr_context.png" alt="Site with genomic and transcript context"
	title="Site with genomic and transcript context" width="700" />

(A) illustrates the usual way to extract CLIP-seq binding region sequences, after mapping of CLIP-seq reads to the genome and peak calling to identify the binding regions. The context sequence is obtained directly from the genome. In contrast, (B) shows the region mapped to the underlying transcript, from which CLIPcontext then takes the sequence to extract the (possibly more authentic) transcript context.

In **G2T** mode, CLIPcontext essentially takes care of the following tasks:

- Mapping of genomic RBP binding regions to underlying transcripts
- Optionally merge adjacent peak regions at exon borders (keep site with highest score)
- Extract both genomic and transcript context sets (BED regions + FASTA sequences) for comparative analysis

### T2G

In **T2G** mode, CLIPcontext maps transcript binding sites back to the genome. Again, context sequence sets (BED regions + FASTA sequences) are extracted for both genomic and transcript context. Both full and split matches (if sites overlap exon borders) are output.

### LST

In **LST** mode, CLIPcontext extracts the most prominent transcript for each gene from a given GTF file, producing a list of transcript IDs. The output transcript IDs list file can then be used as input (--tr) for the other modes. Most prominent here is defined as the transcript that is part of the GENCODE basic dataset + having the highest transcript support level (TSL). If there are two or more transcripts of one gene having the same TSL, the longest transcript will be selected. Genes or transcripts that do not meet the filtering criteria will not be output. Note that Ensembl GTF files were used for testing, and that due to possible differences in formatting, other GTF files might not work. 


### INT

In **INT** mode, CLIPcontext maps input sites to intron regions, returning only sites that overlap with intron regions. This can be useful to e.g. observe properties of intron-binding sites in a CLIP set of interest.

### EXB

In **EXB** mode, CLIPcontext extracts binding regions near exon borders from a set of input BED regions. It can be used to e.g. create an input dataset for CLIPcontext, focussing only on regions near exon borders for downstream analysis.

### EIR

In **EIR** mode, CLIPcontext creates exon + intron regions BED files for a given list of transcript IDs. Exon + intron regions are extracted for each input transcript ID.

### MTF

In **MTF** mode, CLIPcontext searches for a motif or regular expression inside the extracted transcript and genomic context sequences (output folder of clipcontext g2t or clipcontext t2g as --in input), and reports counts and frequencies in the two sets. In addition, CLIPcontext mtf can be used to search for a motif in the transcriptome (or any specified set of transcripts), with the hits also being mapped back to the genome (full and split hits).


## Installation

To install CLIPcontext, simply clone the repository and use the Python script within the folder:

```
git clone https://github.com/BackofenLab/CLIPcontext.git
cd CLIPcontext
python -m pip install . --ignore-installed --no-deps -vv
```

CLIPcontext can also be installed via [conda](https://anaconda.org/bioconda/clipcontext). This is the most convenient way to install CLIPcontext, since conda takes care of all the dependencies. Note however that the conda version might not always be the latest release.


### Dependencies
Dependencies for CLIPcontext are as follows:

- python3 (tested with version 3.7.3)
- python libraries: seaborn>=0.10.0, matplotlib>=3.1.3, markdown>=3.2.1, pandas>=1.0.3
- [bedtools](https://github.com/arq5x/bedtools2/releases)  (tested with version 2.29.0) executables in PATH
- [twoBitToFa](http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/twoBitToFa) executable in PATH
- [twoBitInfo](http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/twoBitInfo) executable in PATH


### Dataset requirements
CLIPcontext was implemented to work with and tested on human datasets retrieved from [Ensembl](http://www.ensembl.org/index.html). If you want to use CLIPcontext for different organisms or with datasets from different ressources, feel free to open up an issue.

Depending on the set mode, the following datasets need to be obtained for running CLIPcontext:
- A BED file (6-column format) with (genomic) RBP binding regions (e.g. eCLIP CLIPper peak regions obtained from [ENCODE](https://www.encodeproject.org/))
- A GTF file with genomic annotations from Ensembl (see [download page](http://www.ensembl.org/info/data/ftp/index.html))
- A list of transcript IDs defining the transcriptome to map to (from sequencing data or generate with **LST** mode)
- A genome .2bit file for extracting genomic and transcript sequences (for hg38 assembly click [here](https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.2bit))

### Test run

A small BED file of genomic RBP binding regions as well as a list of transcript IDs is already provided in the test/ subfolder. First we download and store all the remaining necessary datasets in the same folder:
```
cd test/
wget https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.2bit
wget ftp://ftp.ensembl.org/pub/release-98/gtf/homo_sapiens/Homo_sapiens.GRCh38.98.gtf.gz
```
Now we can run CLIPcontext in **G2T** mode on the datasets in the test/ subfolder:
```
clipcontext g2t --in SERBP1_K562_rep1_sites_chr1_hg38.bed --out test_out --tr GRCh38.p12.prominent_isoforms_chr1.out --gtf Homo_sapiens.GRCh38.98.gtf.gz --gen hg38.2bit
```

CLIPcontext will output an overview of the files produced at the end of the run (all stored in set --out folder). For more details see the documentation section below.

## Documentation


An overview of the modes offered by CLIPcontext can be obtained by:

```
clipcontext -h
usage: clipcontext [-h] [-v] {g2t,t2g,lst,int,exb,eir,mtf} ...

Tool suite for mapping RBP binding regions to transcriptome or genome.

positional arguments:
  {g2t,t2g,lst,int,exb,eir,mtf}
                        Program modes
    g2t                 Map genomic sites to transcript sites
    t2g                 Map transcript sites to genomic sites
    lst                 Get list of most prominent transcripts
    int                 Get sites overlapping with introns
    exb                 Get sites near exon borders
    eir                 Get exon and intron regions
    mtf                 Search for motif in g2t or t2g sets

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit

```

### G2T

The following command line arguments are available in **G2T** mode:

```
clipcontext g2t -h
usage: clipcontext g2t [-h] --in str --out str --tr str --gtf str --gen str
                       [--thr float] [--rev-filter] [--min-len int]
                       [--max-len int] [--min-exon-ol float]
                       [--merge-mode {1,2,3}] [--merge-ext int] [--add-out]
                       [--seq-ext int] [--all-gen-out] [--gen-uniq-ids]

optional arguments:
  -h, --help            show this help message and exit
  --thr float           Site score threshold for filtering --in BED file
                        (default: None)
  --rev-filter          Reverse filtering (keep values <= threshold and prefer
                        sites with smaller values) (default: False)
  --min-len int         Minimum input site length for filtering --in BED file
                        (default: False)
  --max-len int         Maximum input site length for filtering --in BED file
                        (default: False)
  --min-exon-ol float   Minimum exon overlap of a site to be reported as
                        transcript hit (intersectBed -f parameter) (default:
                        0.9)
  --merge-mode {1,2,3}  Defines how to merge overlapping transcript sites
                        (overlap controlled by --merge-ext). (1) only merge
                        sites overlapping at exon borders, (2) merge all
                        overlapping sites, (3) do NOT merge overlapping sites
                        (default: 1)
  --merge-ext int       Extend regions mapped to transcripts by --merge-ext
                        before running mergeBed to merge overlapping regions
                        (default: 10)
  --add-out             Output centered and extended sites and sequences for
                        all transcript matches (unique + non-unique) (default:
                        False)
  --seq-ext int         Up- and downstream extension of centered sites for
                        context sequence extraction (default: 30)
  --all-gen-out         Output all centered and extended genomic regions,
                        instead of only the ones with unique transcript
                        matches (default: False)
  --gen-uniq-ids        Generate unique column 4 IDs for --in BED file entries
                        (default: False)
  --report              Output an .html report with statistics and plots
                        comparing transcript and genomic sequences (default:
                        False)

required arguments:
  --in str              Genomic regions (hg38) BED file (6-column format)
  --out str             Output results folder
  --tr str              Transcript sequence IDs list file to define
                        transcripts to map on
  --gtf str             Genomic annotations (hg38) GTF file (.gtf or .gtf.gz)
  --gen str             Genomic sequences (hg38) .2bit file

```

The required arguments comprise the necessary input files (--in genomic RBP binding regions, --gtf genomic annotations, --gen .2bit genomic sequences, --tr transcript IDs list) and the output folder (--out) to store all files produced during the run.
Optional arguments are provided for filtering the input .bed sites (--thr, --rev-filter, --min-len, --max-len), controlling the exon overlap amount (--min-exon-ol), the context sequence extension (--seq-ext), and the merging of nearby sites (--merge-ext, --merge-mode). In case the input .bed file does not have unique column 4 IDs, unique ones can be generated by --gen-uniq-ids. 

#### G2T algorithm description

In **G2T** mode, CLIPcontext first reads in the transcript IDs (--tr) and extracts the corresponding sequences using --gtf and --gen.
After optionally filtering the input sites (--in) by score or length, the full-length sites are mapped to the given set of transcripts. This is mainly done to identify regions at exon borders. CLIPcontext performs mapping to the transcriptome only for genomic sites that show considerable overlap with an exonic region (by default >= 90%, --min-exon-ol 0.9). Identified border regions are merged by default (--merge-mode 1) if they overlap after applying --merge-ext (i.e., for each overlapping set of sites select the site with the highest score). Merging is done since for RBPs that bind to exonic regions, peaks called at adjacent exon ends often originate from the same binding event. In case all overlapping sites should be merged (not only sites at exon borders), use --merge-mode 2. If merging should be skipped altogether, use --merge-mode 3. Only sites that map to a single exonic region (uniquely mapped sites) are kept for subsequent sequence extraction. Non-unique matches are output only in BED format (see unique + non-unique matches on transcripts BED below).
Merged and uniquely mapped sites are then mapped again to the transcriptome, this time only taking their center positions for mapping. After center-position mapping, extension by --seq-ext is performed for both transcript sites and genomic sites to get both their transcript and genomic context sequences. Transcript regions without full extension indicate their location near transcript ends.

At the end of the run, an overview with short discriptions for each output file is printed out (--out test_out, --report enabled):

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
test_out/transcript_sites.cp.bed
Unique transcript matches center positions extended .bed:
test_out/transcript_sites.cp.ext.bed
Unique transcript matches center positions extended .fa:
test_out/transcript_sites.cp.ext.fa

GENOMIC SITES (CORRESPONDING TO MERGED + EXTENDED TRANSCRIPT SITES)
===================================================================
Genomic sites center positions .bed:
test_out/genomic_sites.cp.bed
Genomic sites center positions extended .bed:
test_out/genomic_sites.cp.ext.bed
Genomic sites center positions extended .fa:
test_out/genomic_sites.cp.ext.fa

G2T CONTEXT SET COMPARISON HTML REPORT
======================================
G2T context set comparison report .html:
test_out/report.g2t.html

```
Notice the naming conventions of the output files (cp : center-positioned site, ext : sites extended by --seq-ext). 


#### Output statistics & HTML report

If --report is set, an **HTML report** is output as well, providing comparative statistics for the generated transcript context and genomic context sequences. These currently include: length statistics, sequence complexity distribution, 2-5 mer distributions, target gene biotype statistics, and site overlap statistics. Note that the HTML report only takes into account the centered and extended sequences (.cp.ext.fa). Additional mapping statistics and information for each transcript are stored in the file hit_transcript_stats.out (see below for format).


compare sites containing genomic context with sites containing transcript context

### T2G

The following command line arguments are available in **T2G** mode:

```
clipcontext t2g -h
usage: clipcontext t2g [-h] --in str --out str --gtf str --gen str
                       [--thr float] [--rev-filter] [--min-len int]
                       [--max-len int] [--seq-ext int] [--gen-uniq-ids]
                       [--report]

optional arguments:
  -h, --help      show this help message and exit
  --thr float     Site score threshold for filtering --in BED file (default:
                  None)
  --rev-filter    Reverse filtering (keep values <= threshold and prefer sites
                  with smaller values) (default: False)
  --min-len int   Minimum input site length for filtering --in BED file
                  (default: False)
  --max-len int   Maximum input site length for filtering --in BED file
                  (default: False)
  --seq-ext int   Up- and downstream extension of centered sites for context
                  sequence extraction (default: 30)
  --gen-uniq-ids  Generate unique column 4 IDs for --in BED file entries
                  (default: False)
  --report        Output an .html report with statistics and plots comparing
                  transcript and genomic sequences (default: False)

required arguments:
  --in str        Transcript regions BED file (6-column format) (transcript
                  IDs need to be in --gtf)
  --out str       Output results folder
  --gtf str       Genomic annotations (hg38) GTF file (.gtf or .gtf.gz)
  --gen str       Genomic sequences (hg38) .2bit file

```

A test run (still inside the test/ subfolder) that also produces an HTML report can be evoked by:

```
clipcontext t2g --in test_tr2gen.bed --gtf test_tr2gen.gtf --out test_out_t2g --gen hg38.2bit --report

```

Transcript sites that span exon borders lead to split mappings, where the two (or more in case of long sites or short exons) parts map to different genomic locations. Here the following site IDs get assigned (original site ID: siteid): siteid_p1, siteid_p2. Again, an overview with short discriptions for each output file is printed out (--out test_out_t2g) at the end of the run:


```
....

TRANSCRIPT FILES
================
Filtered transcript sites .bed:
test_out_t2g/transcript_sites.bed
Filtered transcript sites center positions .bed:
test_out_t2g/transcript_sites.cp.bed
Filtered transcript sites extended .bed:
test_out_t2g/transcript_sites.cp.ext.bed
Filtered transcript sites extended .fa:
test_out_t2g/transcript_sites.cp.ext.fa
Transcript exon regions .bed:
test_out_t2g/exon_regions_transcript.bed

GENOMIC FILES (FROM MAPPING FULL-LENGTH TRANSCRIPT SITES)
=========================================================
All genomic matches .bed:
test_out_t2g/genomic_hits.all.bed
Full-length genomic matches .bed:
test_out_t2g/genomic_hits.unique.bed
Split genomic matches .bed:
test_out_t2g/genomic_hits.split.bed
Transcript genomic exon regions .bed:
test_out_t2g/exon_regions_genome.bed

GENOMIC FILES (FROM MAPPING CENTERED TRANSCRIPT SITES)
======================================================
All genomic matches .bed:
test_out_t2g/genomic_hits.all.bed
Genomic matches center positions .bed:
test_out_t2g/genomic_hits.cp.bed
Genomic matches extended .bed:
test_out_t2g/genomic_sites.cp.ext.bed
Genomic matches extended .fa:
test_out_t2g/genomic_sites.cp.ext.fa

T2G CONTEXT SET COMPARISON HTML REPORT
======================================
T2G context set comparison report .html:
test_out_t2g/report.t2g.html

```

The **HTML report** provides the same comparative statistics for the extracted genomic and transcript sequences (.cp.ext.fa) as the report of clipcontext g2t (see details above).


### Additional modes (LST, INT, EXB, EIR, MTF)

Executing the additional modes should be self-explanatory. Here are a few example runs (again inside the test/ subfolder) for the individual modes:

```
clipcontext lst --gtf Homo_sapiens.GRCh38.98.gtf.gz --out prominent_transcripts_gtf.out --strict --add-infos

```

This command extracts a list of prominent transcript IDs from the given --gtf file. Setting --strict leads to a more strict selection, accepting only transcripts with TSL of 1-5. Setting --add-infos adds additional information to the output file, which then becomes a tabular file storing for each transcript ID: gene ID, gene name, gene biotype, transcript ID, transcript biotype, transcript length, number of transcript exons, transcript support level. See clipcontext lst -h for available arguments.

```
clipcontext int --in g2t_test_in.bed --tr g2t_test_in.tr_list --gtf g2t_test_in.gtf --out sites_on_introns.bed

```

This command gets the input sites overlapping with intron regions. Intron regions to consider are defined by a list of transcript IDs (--tr) and a GTF file (--gtf). See clipcontext int -h for available arguments.

```
clipcontext exb --gtf g2t_test_in.gtf --tr g2t_test_in.tr_list --in g2t_test_in.bed --out sites_near_exon_borders.bed

```

This command returns input sites near exon borders (distance to borders is controlled by --max-dist). Sites are returned in BED format (set output file name with --out). See clipcontext exb -h for available arguments.

```
clipcontext eir --tr g2t_test_in.tr_list --gtf g2t_test_in.gtf --exon-out extracted_exon_regions.bed --intron-out extracted_intron_regions.bed

```

This commands extracts exon + intron regions in BED format for a given set of transcripts (--tr).

```
clipcontext mtf --in test_out --motif '[AC]GCGC'

```
This command searches for the motif or regular expression '[AC]GCGC' inside the extracted transcript and genomic context sequences (output of clipcontext g2t), and reports counts and frequencies in the two sets.

```
clipcontext mtf --in pum2_mtf_test.tr_list --motif 'UGUA[ACGU]AUA' --gen hg38.2bit --gtf Homo_sapiens.GRCh38.98.gtf.gz --out pum2_tr_hits_out

```
This command searches for the PUM2 motif termed Pumilio Response Element inside a specified set of transcripts (pum2_mtf_test.tr_list), and also maps the transcript hits back to the genome (reporting both full and split hits).


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
The transcript IDs defined in this file are used to define the transcript set to be used in the different modes.

#### Mapping statistics output file

The mapping statistics output file (**G2T** mode only) stores the different hit counts and additional information for each transcript:
```
head -5 test_out/hit_transcript_stats.out 
tr_id	chr	gen_s	gen_e	pol	gene_id	gene_name	gene_biotype	tr_len	comp_hits	all_hits	uniq_comp_hits	uniq_all_hits
ENST00000361427	chr1	26472439	26476642	+	ENSG00000198830HMGN2	protein_coding	1940	1	1	1	1
ENST00000327300	chr1	32013867	32043877	+	ENSG00000121774KHDRBS1	protein_coding	2713	6	6	6	6
ENST00000368300	chr1	156114710	156140081	+	ENSG00000160789LMNA	protein_coding	3178	11	12	11	12
ENST00000379370	chr1	1020119	1056116	+	ENSG00000188157	AGRN	protein_coding	7326	11	12	11	12
```

Additional information includes the gene ID, gene name, gene biotype, transcript length, and the genomic region coordinates of the transcript. The different hit counts are: # complete (full-length matching) hits, # all hits (complete and incomplete), # unique (matching to one exon/transcript only) + complete hits, and # all unique hits.


