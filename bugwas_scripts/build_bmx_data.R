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

