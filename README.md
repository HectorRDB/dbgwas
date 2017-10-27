# About DBGWAS
DBGWAS is a tool for quick and efficient bacterial GWAS. It uses a compacted De Bruijn Graph (cDBG) structure to represent the variability within all bacterial genome assemblies given as input. Then cDBG nodes are tested for association with a phenotype of interest and the resulting associated nodes are then re-mapped on the cDBG. The output of DBGWAS consists of regions of the cDBG around statistically significant nodes with several informations related to the phenotypes, offering a representation helping in the interpretation. The output can be viewed with any modern web browser, and thus easily shared.

**IMPORTANT: DBGWAS ONLY WORKS ON LINUX FOR THE MOMENT.**

# DBGWAS in a nutshell

For a quick example on how DBGWAS works, we provide the output
obtained from 282 bacterial genomes along with their drug sensitivity
phenotype.

Genomes are taken from van Belkum A et al., "Phylogenetic Distribution
of CRISPR-Cas Systems in Antibiotic-Resistant Pseudomonas
aeruginosa.", MBio, 2015 Nov 24;6(6):e01796-15

Phenotypes are taken from Jaillard M at al. ,"Correlation between
phenotypic antibiotic susceptibility and the resistome in Pseudomonas
aeruginosa". IJAA , 50 (2017) 210-218.

The graph nodes were annotated with a resistance database composed by four sources (ResFinders, MEGAres, ArgAnnot and CARD) and the UniProt database of bacterial proteins (http://www.uniprot.org/).

The output can be found here: http://leoisl.gitlab.io/DBGWAS_support/full_dataset_visualization/

For reproducibility reasons, in the following you have easily the input data, annotation databases and command-line to generate the above output:

1. Strains' sequences along with DBGWAS' required files (input): https://www.dropbox.com/s/0g1llvdbfv1jys6/pseudomonas_aeruginosa_full_dataset.zip?dl=1
2. Resistance database: https://www.dropbox.com/s/mt3g4oh0bt5jwmr/Resistance_DB_for_DBGWAS.fasta?dl=1
3. UniProt database: https://www.dropbox.com/s/9y1p0yw918ips6k/uniprot_sprot_bacteria_for_DBGWAS.fasta?dl=1
4. Extract everything in the bin folder and execute DBGWAS as:
```
./DBGWAS -strains pseudomonas_aeruginosa_full_dataset/strains -newick pseudomonas_aeruginosa_full_dataset/strains.newick -nc_db Resistance_DB_for_DBGWAS.fasta -pt_db uniprot_sprot_bacteria_for_DBGWAS.fasta
```

# Downloading, installing and running
## Downloading the precompiled binaries
This is the easiest way to run the tool since it is already precompiled for Linux AMD64 machines.
Download the latest binary here (v0.3.7): https://www.dropbox.com/s/cl3ai2bybqyor2g/DBGWAS-0.3.7-Linux.tar.gz?dl=1

Pre-compiled versions history:

DBGWAS v0.3.6: https://www.dropbox.com/s/opwncrbhb0uexzo/DBGWAS-0.3.6-Linux.tar.gz?dl=1

DBGWAS v0.3.5: https://www.dropbox.com/s/b2jgqqicymxb93x/DBGWAS-0.3.5-Linux.tar.gz?dl=1

DBGWAS v0.3.4: https://www.dropbox.com/s/0wvxgv6s2tapbk0/DBGWAS-0.3.4-Linux.tar.gz?dl=1

DBGWAS v0.3.3: https://www.dropbox.com/s/538o8p8vyalcsia/DBGWAS-0.3.3-Linux.tar.gz?dl=1

DBGWAS v0.3.2: https://www.dropbox.com/s/vd4kjl1obrkfmey/DBGWAS-0.3.2-Linux.tar.gz?dl=1




## Compiling
If you still want to compile, clone the repository and execute inside the repository directory:
```
mkdir build && cd build && cmake .. && make && make package
```

## Dependencies installation
DBGWAS uses several thirdparty libraries, but most of them were already packed and were statically linked during compilation, so almost no dependencies are needed. However, you still need to:

1. Install bugwas (an R package). Execute these commands:

```
R
install.packages("ape")
install.packages("phangorn")
install.packages("https://raw.githubusercontent.com/sgearle/bugwas/master/build/bugwas_1.0.tar.gz", repos=NULL, type="source")
```

## Running on a sample example

Now that everything is installed, let's try running the tool in a
sample example comprising 50 bacterial genomes (subset of the 282
described in the "DBGWAS in a nutshell" section):

1. Go to the binary folder:
```
cd bin/
```
2. Execute the program, using demo files:
```
./DBGWAS -strains ../sample_example/strains -newick ../sample_example/strains.newick
```
3. The main output, which are subgraphs that can be visualised with any modern web browser, can be found in ```bin/output/visualisations/index.html```
4. For help and understanding the parameters:
```
./DBGWAS -h
```
5. See also the directory ```sample_example``` to understand better this example;
Check at least the file ```sample_example/strains``` to know how to build the input to the program.

## Lineage vs locus effect

DBGWAS is based on bugwas, which is described in the following paper:

Sarah Earle et al., "Identifying lineage effects when controlling for
population structure improves power in bacterial association studies",
Nature Microbiology 2016.

Bugwas offers association tests at the "locus" level (SNPs, kmers, or
in our case unitigs) or at the lineage level, i.e. groups of
SNPs/kmers/unitigs represented by a principal component of the full
design matrix. The latter option needs more memory as it requires a
PCA step. By default, DBGWAS skips lineage effect estimation to save
memory. If you would like to activate lineage effect estimation,
change do.lineage to TRUE in bin/DBGWAS.R.

## Memory and CPU requirements

The analysis on the full dataset with 280 P. aeruginosa strains
described in section "DBGWAS in a nutshell" and in our paper (ref at
the end of the README) produces 54,397,312 kmers which are compressed
into 2,356,052 unitigs corresponding to 1,141,877 unique patterns.

The analysis runs in about 40 minutes on a PC with 8 Intel(R) Xeon(R)
CPU E5-1620 0 @ 3.60GHz cores, although most of the operations only
use a single core at the moment.

With do.lineage=TRUE, memory usage peaks at ~15Gb, vs ~9Gb with
do.lineage=FALSE.

Scaling on new datasets will essentially depend on the number of
strains and the number of unique presence/absence patterns obtained
across unitigs, i.e., on the genetic variability in the dataset.

# Thirdparties
DBGWAS makes use of several thirdparty libraries:
1. GATB (https://github.com/GATB/gatb-core)
2. Boost C++ Libraries (http://www.boost.org/)
3. Bugwas (https://github.com/sgearle/bugwas)
4. GEMMA (https://github.com/genetics-statistics/GEMMA)
5. Blast (https://blast.ncbi.nlm.nih.gov/Blast.cgi)
6. Cytoscape.js (http://js.cytoscape.org/) and these extensions:
   1. cytoscape.js-cxtmenu (https://github.com/cytoscape/cytoscape.js-cxtmenu)
   2. cytoscape-ngraph.forcelayout (https://github.com/Nickolasmv/cytoscape-ngraph.forcelayout)
   3. cytoscape.js-panzoom (https://github.com/cytoscape/cytoscape.js-panzoom)
7. PhantomJS (http://phantomjs.org/)
8. PStreams (http://pstreams.sourceforge.net/)
9. Alasql (https://github.com/agershun/alasql)
10. Handsontable (https://github.com/handsontable/handsontable)
11. Bootstrap (http://getbootstrap.com/javascript/)
12. jQuery (https://jquery.com/)
13. jQuery QueryBuilder (http://querybuilder.js.org/)
14. jQuery Layout (http://layout.jquery-dev.com/)
15. jQuery BlockUI Plugin (https://github.com/malsup/blockui/)
16. Fastclick (https://github.com/ftlabs/fastclick)


# How to cite
Magali Jaillard, Maud Tournoud, Leandro Lima, Vincent Lacroix, Jean-Baptiste Veyrieras and Laurent Jacob, "Representing Genetic Determinants in Bacterial GWAS with Compacted De Bruijn Graphs", 2017,  Cold Spring Harbor Labs Journals, doi:10.1101/113563.

# License
Copyright (C) <2017>  <bioMerieux, Universite Claude Bernard Lyon 1,
Centre National de la Recherche Scientifique>

1. This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as published
by the Free Software Foundation version 3 of the  License and under the
terms of article 2 below.

2. This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE. See below the GNU Affero General Public License for more details.
You should have received a copy of the GNU Affero General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

3. Communication to the public by any means, in particular in the form of
a scientific paper, a poster, a slideshow, an internet page, or a patent,
of a result obtained directly or indirectly by running this program must
cite the following paper :

Magali Jaillard, Maud Tournoud, Leandro Lima, Vincent Lacroix, Jean-Baptiste Veyrieras and Laurent Jacob, "Representing Genetic Determinants in Bacterial GWAS with Compacted De Bruijn Graphs", 2017, Cold Spring Harbor Labs Journals, doi:10.1101/113563.(url: http://www.biorxiv.org/content/early/2017/03/03/113563).

# Affiliations:
bioMérieux (http://www.biomerieux.fr/)
LBBE (https://lbbe.univ-lyon1.fr/?lang=en)
Erable (https://team.inria.fr/erable/en/)