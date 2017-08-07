# About DBGWAS
DBGWAS is a tool for quick and efficient bacterial GWAS. It uses a compacted De Bruijn Graph (cDBG) structure to represent the variability within all bacterial genome assemblies given as input. Then cDBG nodes are tested for association with a phenotype of interest and the resulting associated nodes are then re-mapped on the cDBG. The output of DBGWAS consists of regions of the cDBG around statistically significant nodes with several informations related to the phenotypes, offering a representation helping in the interpretation. The output can be viewed with any modern web browser, and thus easily shared.

**Important: DBGWAS only works on Linux for the moment.**

# DBGWAS in a nutshell
For a quick example on how DBGWAS works, consider the following input example:
288 bacterial genomes testing a phenotype if each strain is resistant or sensitive to a drug (TODO - rephrase this, give link to the input).

The output can be found here: TODO

# Downloading, installing and running
## Downloading the latest release
This is the easiest way to run the tool since it is already pre-compiled for Linux AMD64 machines.
Download the latest release here: TODO

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
Now that everything is installed, let's try running the tool in a sample example comprising 50 bacterial genomes:
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

# Thirdparties
DBGWAS makes use of several thirdparties libraries:
1. GATB (https://github.com/GATB/gatb-core)
2. Boost C++ Libraries (http://www.boost.org/)
3. Bugwas (https://github.com/sgearle/bugwas)
4. GEMMA (https://github.com/genetics-statistics/GEMMA)
5. Cytoscape.js (http://js.cytoscape.org/) and these extensions:
   5.1: cytoscape.js-cxtmenu (https://github.com/cytoscape/cytoscape.js-cxtmenu)
   5.2: cytoscape-ngraph.forcelayout (https://github.com/Nickolasmv/cytoscape-ngraph.forcelayout)
   5.3: cytoscape.js-panzoom (https://github.com/cytoscape/cytoscape.js-panzoom)
6. PhantomJS (http://phantomjs.org/)
7. Alasql (https://github.com/agershun/alasql)
8. Handsontable (https://github.com/handsontable/handsontable)
9. Bootstrap (http://getbootstrap.com/javascript/)
10. JQuery (https://jquery.com/)
11. jQuery QueryBuilder (http://querybuilder.js.org/)
12. Fastclick (https://github.com/ftlabs/fastclick)


# How to cite
Magali Jaillard, Maud Tournoud, Leandro Lima, Vincent Lacroix, Jean-Baptiste Veyrieras and Laurent Jacob, "Representing Genetic Determinants in Bacterial GWAS with Compacted De Bruijn Graphs", 2017,  Cold Spring Harbor Labs Journals, doi:10.1101/113563.

# License
Copyright (C) <2017>  <bioMerieux, Universite Claude Bernard Lyon 1, 
Centre National de la Recherche Scientifique> 

1. This program is free software: you can redistribute it and/or modify 
it under the terms of the GNU Affero General Public License as published 
by the Free Software Foundation version 3 of the  License and under the 
terms of article 2 below.

This program is distributed in the hope that it will be useful, but 
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY 
or FITNESS FOR A PARTICULAR PURPOSE. See below the GNU Affero General  
Public License for more details.

You should have received a copy of the GNU Affero General Public License 
along with this program.  If not, see <http://www.gnu.org/licenses/>.

3. Communication to the public by any means, in particular in the form of 
a scientific paper, a poster, a slideshow, an internet page, or a patent, 
of a result obtained directly or indirectly by running this program must 
cite the following paper :   

Magali Jaillard, Maud Tournoud, Leandro Lima, Vincent Lacroix, Jean-Baptiste Veyrieras and Laurent Jacob, "Representing Genetic Determinants in Bacterial GWAS with Compacted De Bruijn Graphs", 2017, Cold Spring Harbor Labs Journals, doi:10.1101/113563.(url: http://www.biorxiv.org/content/early/2017/03/03/113563).
