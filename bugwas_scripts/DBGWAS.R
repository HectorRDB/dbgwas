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


library(bugwas)

source('cdbg_lin_loc.R')
source('cdbg_all_plots.R')

cleanMem <- function(n=50) { for (i in 1:n) gc() }

## step1.output <- './output'
## pheno.file <- 'output/bugwas_input.id_phenotype'
## tree.file <- '../sample_example/strains.newick'
## ## tree.file <- '../../../data/bmx/bmx.newick'
## prefix <- 'output/bugwas_out_'
## gem.path <- './gemma.0.93b'
## ## bh.thr <- 0.02
## maf.filter <- 0.01

args = commandArgs(trailingOnly=TRUE)
if (length(args) < 6) {
  stop("Six arguments must be supplied.", call.=FALSE)
}

step1.output <- args[1]
pheno.file <- args[2]
tree.file <- args[3]
prefix <- args[4]
gem.path <- args[5]
maf.filter <- as.numeric(args[6])


## If this is TRUE, estimate lineage effect and generate related
## plots. Otherwise skip tree management (the input file can be
## empty), svd, pca and plots. Could be an input of the script.

do.lineage <- FALSE
    
output.dir <- '.' # must be '.' as gemma automatically writes in ./output.

##--------------------------------
## Read gen file output in step 1
##--------------------------------
gen.file <- file.path(step1.output, 'bugwas_input.unique_rows.binary')
message(sprintf('[DBGWAS] Reading unitigs from %s', gen.file))
gen <- read.table(file=gen.file, header=TRUE,
                  row.names=1, check.names=FALSE, colClasses='integer')
cleanMem()


##------------------------------
## Read phenotypes
##------------------------------
message(sprintf('[DBGWAS] Reading phenotypes from %s', pheno.file))
pheno.mat <- read.table(file=pheno.file, header=TRUE)
if((ncol(gen) != nrow(pheno.mat)) || any(colnames(gen) != pheno.mat['ID'])){
    stop('Mismatch between genotype and phenotype ids (should be in the same order)')
}

##-------------------------------------------------------------------
## Restrict genotype, phenotype and tree to annotated strains (gemma
## deals with missing phenotypes but bugwas does not, first bug
## encountered in bugwas:::ridge_regression)
##-------------------------------------------------------------------

## annotated.sample <- !is.na(pheno.mat['pheno'])
## message(sprintf('[DBGWAS] Restricting genotype, phenotype and tree to %d/%d annotated strains.', sum(annotated.sample), length(annotated.sample)))
## pheno.mat <- pheno.mat[annotated.sample, ] # Restrict phenotype
restr.pheno.file <- paste0(prefix, "_restricted_pheno.txt")
write.table(pheno.mat, file=restr.pheno.file, row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t')
## gen <- gen[, annotated.sample] # Restrict genotype
## ## Need to re-switch minor allele to 1
## rmg.mask <- (rowMeans(gen) > 0.5)
## cleanMem()
## gen[rmg.mask, ] <- 1L - gen[rmg.mask, ]
XX.ID <- colnames(gen)

if(do.lineage){
    ## Restrict tree
    ## library(ape)
    ## bmx.tree <- read.tree(tree.file)
    ## restr.tree <- drop.tip(bmx.tree, setdiff(bmx.tree$tip.label, XX.ID))
    ## restr.tree.file <- paste0(prefix, "_restricted_tree.txt")
    ## write.tree(restr.tree, file=restr.tree.file)
    restr.tree.file <- tree.file
}else{
    restr.tree.file <- NULL
}

##-----------------
## Filter patterns
##-----------------

## Takes up memory (~3x design size)

## Remove mono-allelic (constant) patterns (even if maf.filter=0)
polyallelic.mask <- (apply(gen, 1, FUN=function(v) length(unique(v))) > 1)
cleanMem()
## Filter on MAF
maf.mask <- (rowMeans(gen) >= maf.filter)
pattern.mask <- polyallelic.mask & maf.mask
cleanMem()
message(sprintf('[DBGWAS] Restricting genotype %d/%d patterns with MAF >= %g.', sum(pattern.mask), length(pattern.mask), maf.filter))
gen <- gen[pattern.mask, ]

##---------------------------------------------------------------------
## Build the objects required by gemma (to build relmatrix) and bugwas
## (for the GWAS)
##---------------------------------------------------------------------
snps.by.pattern <- file.path(step1.output, 'gemma_input.pattern_to_nb_of_unitigs.binary')
snp.to.pattern <- file.path(step1.output, 'gemma_input.unitig_to_pattern.binary')

XX.all <- list()
XX.all$XX <- gen
## gc(verbose=FALSE)

bippat <- read.table(file=snps.by.pattern)
XX.all$bippat <- bippat[, -1]
names(XX.all$bippat) <- as.character(bippat[, 1])
XX.all$bippat <- XX.all$bippat[pattern.mask]
rm(gen, polyallelic.mask, pattern.mask)
cleanMem()

pattern <- read.table(file=snp.to.pattern, colClasses=rep('character', 2))
XX.all$pattern <- pattern[, -1]
names(XX.all$pattern) <- pattern[, 1]
## Remove unitigs corresponding to filtered out patterns
XX.all$pattern <- XX.all$pattern[XX.all$pattern %in% names(XX.all$bippat)]
## Because run_lmm_bi calls gemma with hard coded pattern names
## 1:nrow(XX), XX.all$pattern needs to index rows of XX.all$XX, not
## arbitrary pattern IDs. We will need to convert back the result to
## know which patterns were selected using names(XX.all$bippat).
XX.all$ps <- names(XX.all$pattern)
XX.all$pattern <- match(XX.all$pattern, names(XX.all$bippat))
names(XX.all$pattern) <- XX.all$ps

sample_ID <- XX.ID
npcs <- length(sample_ID)
y <- pheno.mat[, 'pheno']

SNPdata <- list(XX.all=XX.all, sample_ID=sample_ID, npcs=npcs, y=y, XX.ID=XX.ID)
rm(XX.all)
cleanMem()

## Gemma apparently filters out patterns with maf <
## 2*maf.filter. Deactivate filtering since we did it upstream anyway.
maf.filter <- 0

##-----------------
## Build relmatrix
##-----------------
message('[DBGWAS] Building kinship matrix')

## gemma automatically writes the matrix in a ./output directory, no
## matter what dir and relmatrix are. For consistency, dir should be
## . (to which get_kinship automatically adds /output).
gemma.pheno.file <- bugwas:::write_pheno(pheno = y, prefix = prefix)
relmatrix <- bugwas:::get_kinship(XX = SNPdata$XX.all$XX,
                                  pattern = SNPdata$XX.all$pattern,
                                  prefix = prefix,
                                  path = gem.path,
                                  dir = '.',
                                  maf=maf.filter, 
                                  pheno.file = gemma.pheno.file)

##---------------------------------------------------------------------
## Call a modified version of lin_loc which tests for both lineage and
## locus effects and generates a bunch of plots.
##---------------------------------------------------------------------
message('[DBGWAS] Performing association tests')

data <- cdbg_lin_loc(SNPdata=SNPdata,
                    pheno=restr.pheno.file,
                    phylo=restr.tree.file,
                    prefix=prefix,
                    gem.path=gem.path,
                    relmatrix=relmatrix,
                    output.dir=output.dir,
                    maf=maf.filter,
                    creatingAllPlots=do.lineage)

## svd.XX and pca can be recovered from data$biallelic and
## used as arguments for future calls of cdbg_lin_loc (eg with
## different phenotypes):
## svd.XX <- data$biallelic$svd.XX
## pca <- list(pca=data$biallelics$pca)

## NOTE:
## Warning in cor(mtp, pca$pca$x[, 1:npcs]) : the standard deviation
## is zero. Not surprising according to S. Earle (arises from XX being
## rank deficient and since we are keeping all PCs some of them
## correspond to 0 variance).

##-----------------------------------------------
## Extract and write list of patterns
##-----------------------------------------------

## Could also use data$lmm, but would need to de-duplicate to apply
## multiple testing correction procedure.

lmm.result <- read.table(file=file.path(step1.output, sprintf('%s_biallelic_lmmout_patterns.assoc.txt', prefix)), as.is=TRUE, header=TRUE)
bh.pv <- p.adjust(lmm.result$p_lrt, method='BH')

## bh.pv.mask <- (bh.pv < bh.thr)
## save.image('DBGWAS.RData')
## write.table(file=sprintf('%s_significant_patterns_BH%g.txt', prefix, bh.thr), cbind(names(SNPdata$XX.all$bippat)[lmm.result$ps[bh.pv.mask]],bh.pv[lmm.result$ps[bh.pv.mask]]), quote=FALSE, row.names=FALSE, col.names=FALSE)

write.table(file=sprintf('%s_DBGWAS_patterns.txt', prefix),
            cbind(names(SNPdata$XX.all$bippat)[lmm.result$ps],
                  bh.pv[lmm.result$ps], # q-values
                  lmm.result$beta[lmm.result$ps], # effect in the linear model
                  (lmm.result$beta/lmm.result$se)[lmm.result$ps]), # Wald statistic
            quote=FALSE, row.names=FALSE, col.names=FALSE)

