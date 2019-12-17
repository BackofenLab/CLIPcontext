# CLIPcontext
CLIPcontext takes genomic RBP binding regions or sites identified by CLIP-seq, maps them to the transcriptome, 
and retrieves the region sequences with both genomic and transcript sequence context. Depending on the location of the binding regions, this leads to the extraction of two different sequence contexts:

<img src="doc/gen_tr_context.png" alt="Site with genomic and transcript context"
	title="Site with genomic and transcript context" width="700" />

(A) illustrates the usual way to extract CLIP-seq binding region sequences, after mapping of CLIP-seq reads to the genome and peak calling to identify the binding regions. The context sequence is obtained directly from the genome. In contrast, (B) shows the region mapped to the underlying transcript, from which CLIPcontext then takes the sequence to extract the (possibly more authentic) transcript context.

CLIPcontext essentially takes care of the following tasks:

- Mapping of genomic RBP binding regions to underlying transcripts
- Merge genomic regions adjacent to each other at exon borders (keep site with highest score)
- Extract both genomic and transcript context sequences for comparative analysis


## Installation

To install CLIPcontext, simply clone the repository and use the Python script within the folder:

```
git clone https://github.com/michauhl/CLIPcontext.git
cd CLIPcontext
python CLIPcontext.py -h
```

### Dependencies
Dependencies for CLIPcontext are as follows:

- python3 (tested with version 3.7.3)
- [bedtools](https://github.com/arq5x/bedtools2/releases)  (tested with version 2.26.0)
- [twoBitToFa](http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/twoBitToFa) executable in PATH

### Dataset requirements
CLIPcontext was implemented to work with and tested on human datasets retrieved from [Ensembl](http://www.ensembl.org/index.html). If you want to use CLIPcontext for different organisms or with datasets from different ressources, feel free to open up an issue.

The following datasets need to be obtained for running CLIPcontext:
- A BED file (6-column format) with genomic RBP binding regions (e.g. obtained from [ENCODE](https://www.encodeproject.org/))
- A GTF file with genomic annotations from Ensembl (see [download page](http://www.ensembl.org/info/data/ftp/index.html))
- A transcript sequences FASTA file from Ensembl (see [download page](http://www.ensembl.org/info/data/ftp/index.html))
- A list of transcript IDs which defining the transcriptome to map to
- A genome .2bit file for extracting genomic sequences (for hg38 assembly click [here](https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.2bit))

### Test run
A small BED file of genomic RBP binding regions as well as a list of transcript IDs is already provided in the data/ subfolder. First we download and store all the remaining necessary datasets in the same folder:
```
cd data/
wget https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.2bit
wget ftp://ftp.ensembl.org/pub/release-98/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
wget ftp://ftp.ensembl.org/pub/release-98/fasta/homo_sapiens/ncrna/Homo_sapiens.GRCh38.ncrna.fa.gz
cat Homo_sapiens.GRCh38.cdna.all.fa.gz Homo_sapiens.GRCh38.ncrna.fa.gz > Homo_sapiens.GRCh38.cdna.ncrna.fa.gz
wget ftp://ftp.ensembl.org/pub/release-98/gtf/homo_sapiens/Homo_sapiens.GRCh38.98.gtf.gz
cd ..
```
Now we can run CLIPcontext on the datasets stored in the data/ subfolder:
```
python CLIPcontext.py --in data/SERBP1_K562_rep1_sites_chr1_hg38.bed --out test_out --tr data/GRCh38.p12.prominent_isoforms_chr1.out --fa data/Homo_sapiens.GRCh38.cdna.ncrna.fa.gz --gtf data/Homo_sapiens.GRCh38.98.gtf.gz --gen data/hg38.2bit
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
The required arguments comprise the necessary input files (--in genomic RBP binding regions, --fa transcript sequences, --gtf genomic annotations, --gen .2bit genomic sequences, --tr transcript IDs list) and the output folder to store all files produced during the run.
Optional arguments are provided for filtering the input .bed sites (--thr, --rev-filter, --min-len, --max-len), controlling the exon overlap amount (--min-exon-ol), the context sequence extension (--seq-ext), and the merging of nearby sites (--merge-ext, --merge-all). In case the input .bed file does not have unique column 4 IDs, new ones can be generated by --gen-uniq-ids.

### Algorithm description

CLIPcontext first reads in the transcript IDs (--tr), and extracts the corresponding sequences from the --fa FASTA file. Transcript IDs without corresponding sequence will be ignored, thus also genomic regions that map to these transcripts. 
After the filtering the input sites (--in) by score or length, the full-length sites are mapped to the given set of transcripts. This is done to get regions at exon borders. Border regions are then merged 

By default transcript regions near exon borders 
    are merged if they overlap after applying --merge-ext 
    (i.e. for each overlapping set of sites select the site with the 
    highest score). In case all overlapping sites should be merged, 
    use --merge-all.

Filter
Read in transcript IDs, get sequences 




Merge overlapping sites, keep only highest-scoring site for each set of 
overlapping regions.
By default, merging takes place only for sites near exon borders 
(more precisely, sites that do not fully overlap with exons after 
extending them by --merge-ext). However, if --merge-all is set, 
merging will be done for all overlapping sites.


Now map center position genomic sites to transcriptomes.
Use only merged sites that uniquely mapped in the first step,
with their IDs stored in ids2keep_dic. After center-position mapping,
extension by --seq-ext will be done for both transcript sites and 
genomic sites to get both transcript and genomic context sites.


Get site ID / sequence ID combinations to keep.
We need these since mapping center positions to transcriptome can result 
in sites which mapped uniquely using their full lengths before, but with 
center positions could map to more than one exon / transcript. This 
is due to the intersectBed_f overlap parameter set to 0.9, which does not 
report hits if they do not overlap 90%+ with exons, while for center 
position mapping this is always satisfied.


Next, do a sanity check and extract genomic context sequences.
1) Sanity check:
Compare extracted transcript center position nucleotides with
with genomic center position nucleotides.
2) Extract genomic sequences for given input regions,
using up- and downstream extension used for transcript sequences.


CLIPcontext maps genomic regions (--in) to a defined transcriptome (--tr) and outputs transcript and genomic region sequences. CLIPcontext performs mapping to the transcriptome for genomic sites that show considerable overlap with an exonic region (>= 90%, --min-exon-ol 0.9). 




CLIPcontext performs mapping to the transcriptome for genomic sites that show considerable overlap with an exonic region (>= 90%, --min-exon-ol 0.9). 



Consider a region overlapping with exons




The idea of CLIPcontext is to get both the genomic and transcript sequence context for a given set of RBP binding regions (--in). It is based on the assumption that binding regions that map to exons (minimum overlap set by --min-exon-ol)





    CLIPcontext maps genomic regions to a defined transcriptome and outputs 
    transcript and genomic region sequences. 
    
    
    Note that only regions
    uniquely mapped to transcripts (i.e. in total one transcript hit) 
    are reported. Moreover, the exon overlap of a genomic region has 
    to be >= --min-exon-ovlp for the region to be reported. 
    Both transcript and genomic region sequences are extracted by 
    centering the input regions and extending them by --seq-ext. 
    Transcript regions without full extension indicate their location 
    near exon ends. By default transcript regions near exon borders 
    are merged if they overlap after applying --merge-ext 
    (i.e. for each overlapping set of sites select the site with the 
    highest score). In case all overlapping sites should be merged, 
    use --merge-all.



Here's how it works:



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
