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
First download and store the datasets necessary for conducting a test run in the data/ subfolder:
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
CLIPcontext will output an overview of the files produced at the end of the run. For more details see the documentation section below.

## Documentation

Here's how it works:

<img src="doc/gen_tr_context.png" alt="Site with genomic and transcript context"
	title="Site with genomic and transcript context" width="700" />
