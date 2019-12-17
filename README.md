# CLIPcontext
CLIPcontext takes genomic RBP binding regions identified by CLIP-seq, maps them to the transcriptome, 
and retrieves the region sequences with both genomic and transcript sequence context.



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

Here's how it works:

<img src="doc/gen_tr_context.png" alt="Site with genomic and transcript context"
	title="Site with genomic and transcript context" width="700" />


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
The transcript IDs defined in this file are used for mapping the input genomic regions to transcript regions.
