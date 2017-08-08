## Copyright (C) <2017>  <bioMerieux, Universite Claude Bernard Lyon 1, 
## Centre National de la Recherche Scientifique> 

## 1. This program is free software: you can redistribute it and/or modify 
## it under the terms of the GNU Affero General Public License as published 
## by the Free Software Foundation version 3 of the  License and under the 
## terms of article 2 below.
## 2. This program is distributed in the hope that it will be useful, but 
## WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY 
## or FITNESS FOR A PARTICULAR PURPOSE. See below the GNU Affero General  
## Public License for more details.
## You should have received a copy of the GNU Affero General Public License 
## along with this program.  If not, see <http://www.gnu.org/licenses/>.
## 3. Communication to the public by any means, in particular in the form of 
## a scientific paper, a poster, a slideshow, an internet page, or a patent, 
## of a result obtained directly or indirectly by running this program must 
## cite the following paper :   
##  Magali Jaillard, Maud Tournoud, Leandro Lima, Vincent Lacroix, 
##  Jean-Baptiste Veyrieras and Laurent Jacob, "Representing Genetic 
##  Determinants in Bacterial GWAS with Compacted De Bruijn Graphs", 2017,
##  Cold Spring Harbor Labs Journals, doi:10.1101/113563.
##  (url: http://www.biorxiv.org/content/early/2017/03/03/113563)
## -------------------------------------------------------------------------

## Authors (alphabetically): Jacob L., Jaillard M., Lima L.

## Writing pheno file (strainsRS option) in [data.dir]/[drug]_pheno.txt
## and tree file in [data.dir]/bmx.newick

library(ape)

data.dir <- '../../../data/bmx'
strainIds <- file.path(data.dir, 'strainsIds') # Sorted list of used samples ids
strains <- file.path(data.dir, 'strains') # Sorted list of used samples path
full.tree.file <- file.path(data.dir, 'VTK_simuPheno-subtree.newick') # Input tree
tree.file <- file.path(data.dir, 'bmx.newick') # Output tree (genbank names)

drug <- 'Amikacin'
## drug <- 'Levofloxacin'

##---------------------------------------------------------------
## Read strain ids, just to select and sort the phenotype vector
##---------------------------------------------------------------

sample.id <- unlist(read.table(strainIds, as.is=TRUE))

##------------------------------
## Read phenotypes
##------------------------------

pheno <- read.csv(file=file.path(data.dir, 'mmc6.csv'))
bmx.id <- pheno$Isolate.id
rownames(pheno) <- bmx.id
## Supp table 3 of
## http://mbio.asm.org/content/6/6/e01796-15.full#sec-16. Necessary to
## join the pheno id (bioMerieux id) to the genotype id (genbank ids)
junction <- read.csv(file=file.path(data.dir, 'mbo006152563st1.csv'), colClasses=c(rep('character', 2), rep('NULL', 156)))
rownames(junction) <- junction[, 'name']
bmxID.from.SampleID <- junction[sample.id, 'bm_name']
pheno <- pheno[bmxID.from.SampleID, ]

annotated.sample <- !is.na(pheno[, drug])
if(!all(pheno[annotated.sample, drug] %in% c('S', 'I', 'R'))){
    stop('Phenotypes should be either NA, S, I or R')
}

##-----------------------------------------
## Build phenotype vector and strains file
##-----------------------------------------

pheno.vec <- rep(NA, nrow(pheno))
pheno.vec[annotated.sample] <- as.numeric(pheno[annotated.sample, drug] %in% c('I', 'R')) # Encode cases (I or R) as 1, controls (S) as 0
## pheno.mat <- cbind(sample.id, pheno.vec)
## colnames(pheno.mat) <- c("ID", "pheno")
## pheno.file <- file.path(data.dir, sprintf('%s_pheno.txt', drug))
## write.table(pheno.mat, file=pheno.file, row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t')
sample.path <- unlist(read.table(strains, as.is=TRUE))
strains.mat <- cbind(sample.id, pheno.vec, sample.path)
colnames(strains.mat) <- c("ID", "pheno", "Path")
strains.file <- file.path(data.dir, sprintf('%s_strains.txt', drug))
write.table(strains.mat, file=strains.file, row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t')



##-----------------------------------------------------------------
## Read tree, re-write it using names consistent genbank ids rather
## thank bioMerieux ids
##-----------------------------------------------------------------

bmx.tree <- read.tree(full.tree.file)
rownames(junction) <- junction[, 'bm_name']
bmx.tree$tip.label <- junction[bmx.tree$tip.label, 'name']
write.tree(bmx.tree, file=tree.file)

## ## Restrict tree to annotated strains
## restr.tree <- drop.tip(bmx.tree, setdiff(bmx.tree$tip.label, sample.id[annotated.sample]))
## write.tree(restr.tree, file=tree.file)

