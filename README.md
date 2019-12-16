# CLIPcontext
CLIPcontext takes genomic RBP binding regions identified by CLIP-seq, maps them to the transcriptome, 
and retrieves the region sequences with both genomic and transcript sequence context.



## Installation

To install CLIPcontext, simply clone the repository and use the tool within the folder:

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



## Documentation

Here's how it works:

<img src="doc/figures/gen_tr_context.png" alt="Site with genomic and transcript context"
	title="Site with genomic and transcript context" width="700" />
