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

## Modified (2017) from bugwas (copyright, authors and license below)
## file BUGWAS_modular.R. Appending "cdbg" to the modified functions.

## # Authors: Earle, S. G., Wu, C.-H. and Wilson, D. J.
## #
## # Copyright (C) 2015 University of Oxford
## #
## # This file is part of the bugwas R package.
## # See the NOTICE file distributed with this work for additional
## # information regarding copyright ownership and licensing.
## #
## # This is free software; you can redistribute it and/or modify
## # it under the terms of the GNU Lesser General Public License as
## # published by the Free Software Foundation; either version 2
## # of the License, or (at your option) any later version.
## #
## #  This is distributed in the hope that it will be useful,
## #  but WITHOUT ANY WARRANTY; without even the implied warranty of
## #  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## #  GNU Lesser General Public License for more details.
## #
## # You should have received a copy of the GNU Lesser General Public
## # License along with this software package; if not, write to the
## # Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
## # Boston, MA  02110-1301  USA

############################################################
## BUGWAS				   						
############################################################

##################################################################################################################
#### Library dependencies
##################################################################################################################

## library(ape)
## library(phangorn)

## source("/dipro/mmm/ana/Saur/AbxGWAS/bugwas/BUGWAS_functions.R")
## source("/home/wilson/R/ridge regression.R")
## source("/dipro/mmm/ana/Saur/AbxGWAS/bugwas/all_plots_new.R")
## source("/dipro/mmm/ana/Saur/AbxGWAS/bugwas/jessie_plots_code.R")


##' Lineage and locus tests for bacterial GWAS
##'
##' This function tests for locus effects using GEMMA and lineage effects using a bayesian wald test for haploid data
##' @param gen A file name specified by either a variable of mode character, or a double-quoted string, containing imputed haploid SNP data. Rows are SNPs, and columns are samples, with the first column being SNP positions. Column headers must contain 'ps' for the SNP positions with the others being the sample names. This must contain biallelic SNPs, but can also contain tri and tetra-allelic SNPs. Required argument.
##' @param pheno A file name specified by either a variable of mode character, or a double-quoted string, containing a column of sample names with header 'ID' and a column of the binary phenotype (coded by 0s and 1s) with column header 'pheno'. Required argument.
##' @param phylo A file name specified by either a variable of mode character, or a double-quoted string, containing a phylogeny of the samples, with the same names matching with arguments gen and pheno. Required argument.
##' @param prefix Output file prefix. Required argument.
##' @param gem.path A file path specified by either a variable of mode character, or a double-quoted string. gem.path is the file path to the software GEMMA (version >= ?). Required argument.
##' @param pcs A file name specified by either a variable of mode character, or a double-quoted string, containing the principle components of the data. Column names should be 'PC1' to 'PCn' and row names should be the sample names.
##' @param lmm.bi A file name specified by either a variable of mode character, or a double-quoted string, containing GEMMA results (ending '.assoc.txt') for the biallelic SNPs in argument 'gen'.
##' @param logreg.bi A file name specified by either a variable of mode character, or a double-quoted string, containing logistic regression -log10(p) for the biallelic SNPs with column names 'ps' for SNP positions/IDs and 'negLog10' for -log10(p).
##' @param cutOffCor Correlation cut-off for assigning and colouring variants by Principal Components (Default = 0, variants are coloured by the PC they are most correlated with).
##' @param run.lmm Whether to run GEMMA (Default = TRUE).
##' @param maf Minor allele frequency for GEMMA (Default = 0, all varaints are tested).
##' @param relmatrix A file name specified by either a variable of mode character, or a double-quoted string of a file containing the GEMMA relatedness matrix of the samples created from biallelic SNPs. The individual ordering must be in the same order as the column names in argument 'gen'.
##' @param lognull The log likelihood under the null from GEMMA.
##' @param lambda Lambda from GEMMA.
##' @param output.dir Output file directory.
##' @param creatingAllPlots Whether to create all bugwas plots. Default = TRUE.
##' @param allBranchAndPCCor Whether or not to retreive correlation matrix between branches and PCs. Default = FALSE.
##' @param svd.XX The output of svd(XX), where XX is a rescaled version of SNPdata$XX.all$XX. Optional, only used if creatingAllPlots is TRUE. Default = NULL (computed if creatingAllPlots is TRUE).
##' @param pca The output of bugwas:::do_pca(pcs = pcs, XX = XX, XX.ID = SNPdata$XX.ID), where XX is a rescaled version of SNPdata$XX.all$XX. Optional, only used if creatingAllPlots is TRUE. Default = NULL (computed if creatingAllPlots is TRUE).
##' @keywords bacteria GWAS locus lineage wald GEMMA
##' @export
##' @examples
##' lin_loc()
##' #### An example of running lin_loc with the minimum required inputs
##' #### Assuming gemma is installed in the present working directory
##' gen <- system.file("extdata", "gen.txt", package = "bugwas")
##' phylo <- system.file("extdata", "tree.txt", package = "bugwas")
##' prefix <- "test_bugwas"
##' gem.path <- "./gemma"
##' data <- lin_loc(gen = gen, pheno = pheno, phylo = phylo, prefix = prefix, gem.path = gem.path)


cdbg_lin_loc <- function(SNPdata = NULL,
                         phylo = NULL,                         
                         cov.file = NULL,
                         prefix = NULL,
                         gem.path = NULL,
                         pcs = NULL,
                         lmm.bi = NULL,
                         logreg.bi = NULL,
                         cutOffCor = 0,
                         run.lmm = TRUE,
                         maf = 0,
                         relmatrix = NULL,
                         lognull = NULL,
                         lambda = NULL,
                         output.dir = getwd(),
                         creatingAllPlots = TRUE,
                         allBranchAndPCCor = FALSE,
                         svd.XX = NULL,
                         pca = NULL){
    
    phylo = bugwas:::extractInputArgument(arg = phylo, canBeNULL=TRUE, checkExist = TRUE)
    prefix = bugwas:::extractInputArgument(arg = prefix)
    gem.path = bugwas:::extractInputArgument(arg = gem.path, checkExist = TRUE)
    pcs = bugwas:::extractInputArgument(arg = pcs, canBeNULL = TRUE)
    lmm.bi = bugwas:::extractInputArgument(arg = lmm.bi, canBeNULL = TRUE, checkExist = TRUE)
    logreg.bi = bugwas:::extractInputArgument(arg = logreg.bi, canBeNULL = TRUE, checkExist = TRUE)
    cutOffCor = bugwas:::extractInputArgument(arg = cutOffCor, default = 0)
    run.lmm = bugwas:::extractInputArgument(arg = run.lmm, default = TRUE)
    maf = bugwas:::extractInputArgument(arg = maf, default = 0)
    relmatrix = bugwas:::extractInputArgument(arg = relmatrix, canBeNULL = TRUE, checkExist = TRUE)
    lognull = bugwas:::extractInputArgument(arg = lognull, canBeNULL = TRUE)
    lambda = bugwas:::extractInputArgument(arg = lambda, canBeNULL = TRUE)
    output.dir = bugwas:::extractInputArgument(arg = output.dir, default = getwd())
    creatingAllPlots = bugwas:::extractInputArgument(arg = creatingAllPlots, default = TRUE)
    allBranchAndPCCor = bugwas:::extractInputArgument(arg = allBranchAndPCCor, default = FALSE)

    XX.all <- SNPdata$XX.all
    sample_ID <- SNPdata$sample_ID
    y <- SNPdata$y
    npcs <- min(SNPdata$npcs, sum(!is.na(y)))
    XX.ID <- SNPdata$XX.ID
    rm(SNPdata)
    cleanMem()
    
    bugwas:::get_log_file(XX.all = XX.all, prefix = prefix)
    
    ## Read in kinship matrix, if no kinship matrix to read in, compute kinship matrix
    if(run.lmm){
        pheno.file <- bugwas:::write_pheno(pheno = y, prefix = prefix)
        if(is.null(relmatrix)){
            message("Calculating kinship matrix.")
            relmatrix <- bugwas:::get_kinship(XX = XX.all$XX,
                                              pattern = XX.all$pattern,
                                              prefix = prefix,
                                              path = gem.path,
                                              dir = output.dir,
                                              maf = maf,
                                              pheno.file = pheno.file)
            message("Kinship matrix calculated successfully.")
        }
    }
    

    if(creatingAllPlots){## XX <- bugwas:::rescale_variants(var = XX.all$XX[, !is.na(y)], varpat = XX.all$bippat)
        
        ## We choose to retain strains with no phenotype for PCA (could better highlight population structure)
        
        XX <- bugwas:::rescale_variants(var = XX.all$XX, varpat = XX.all$bippat)
        cleanMem()
        message("Rescaled variants.")

    }else{
        XX <- NULL
    }

    ## Since XX is column-centered, SVD info is contained in PCA    
    ## if(creatingAllPlots && is.null(svd.XX)){
    ##     svd.XX <- svd(XX)
    ##     ## Since XX is column-centered, PCA is the same as SVD
    ##     ## pca <- list()
    ##     ## pca$sdev <- svd.XX$d / sqrt(max(1, nrow(XX) - 1))
    ##     ## pca$rotation <- svd.XX$v
    ##     ## pca$x <- svd.XX$u %*% diag(svd.XX$d)
    ##     ## pca <- list(pca=pca)
    ##     cleanMem()
    ##     message("Singular value decomposition complete.")
    ## }

    if(creatingAllPlots && is.null(pca)){## PCA on the bips
        ## Since XX is column-centered, PCA is the same as SVD
        pca <- svd(XX)
        names(pca)[names(pca) == 'v'] <- 'rotation'
        names(pca)[names(pca) == 'u'] <- 'x'        
        pca$x <- pca$x %*% diag(pca$d)
        pca <- list(pca=pca)
        ## pca <- bugwas:::do_pca(pcs = pcs, XX = XX, XX.ID = XX.ID)        
        ## Singular values of XX will be useful later, differ from
        ## sdev by a factor sqrt(nrow(XX)-1)        
        ## pca$pca$d <- pca$pca$sdev * sqrt(max(1, nrow(XX) - 1))        
        cleanMem()
        message("Principal component analysis complete.")
    }    
    if(!is.null(pca)){
        ## Samples with missing phenotype will not be represented in the tree.
        pca$pca$x <- pca$pca$x[!is.na(y), 1:sum(!is.na(y))]
        pca$pca$d <- pca$pca$d[1:sum(!is.na(y))]
        pca$pca$rotation <- pca$pca$rotation[, 1:sum(!is.na(y))]
    }
    
    biallelic <- cdbg_get_biallelic(logreg.bi = logreg.bi,
                                    XX.all = XX.all,
                                    XX = XX[!is.na(y), ],
                                    lmm.bi = lmm.bi,
                                    lognull = lognull,
                                    lambda = lambda,
                                    relmatrix = relmatrix,
                                    pheno.file = pheno.file,
                                    cov.file = cov.file,
                                    maf = maf,
                                    gem.path = gem.path,
                                    output.dir = output.dir,
                                    prefix = prefix,
                                    run.lmm = run.lmm,
                                    XX.ID = XX.ID,
                                    pca = pca$pca,
                                    npcs = npcs)
    message("Biallelic data processed successfully.")
    
    if(creatingAllPlots){
        ## Get list of all tree info
        treeInfo <- cdbg_get_tree(phylo = phylo,
                                  prefix = prefix,
                                  XX.ID = XX.ID[!is.na(y)],
                                  pca = pca$pca,
                                  npcs = npcs,
                                  allBranchAndPCCor = allBranchAndPCCor)
        message("Tree data processed successfully.")
        
        ## Ridge regression
        XX.all$XX <- XX.all$XX[, !is.na(y)]
        wald <- cdbg_wald_test(y = y[!is.na(y)],
                               XX = XX[!is.na(y), ],
                               lambda = biallelic$lambda,
                               XX.all = XX.all,
                               prefix = prefix,
                               npcs = npcs,
                               pca = pca$pca)
        
        
        ## rm(list=c("XX", "svd.XX"))
        rm(XX)
        cleanMem()
    }else{
        treeInfo <- NULL
        wald <- NULL
    }
    
    ## Put everything into lists
    biallelic <- list("pattern" = XX.all$pattern,
                      "cor.XX" = biallelic$cor.XX,
                      "npcs" = npcs,
                      "pheno" = y,
                      "logreg" = biallelic$logreg.bi,
                      "lmm" = biallelic$lmm.bi,
                      "ps" = XX.all$ps,
                      "pred" = wald$pred,
                      "pca" = pca$pca,
                      "pc_order" = wald$pc_order,
                      "p.pca.bwt" = wald$p.pca.bwt,
                      "bippat" = XX.all$bippat,
                      "id" = XX.ID, ## New  - Jessie need to add
                      "svd.XX" = svd.XX)
    
    config <- list("prefix" = prefix,
                   "signif_cutoff" = wald$signif_cutoff,
                   "cutoffCor" = cutOffCor)
    
    if(creatingAllPlots){
        ## Plot will only involve phenotyped strains (tree was
        ## restricted to phenotyped strains too).
        biallelic.plot <- biallelic
        keep.mask <- !is.na(biallelic.plot$pheno)
        biallelic.plot$pheno <- biallelic.plot$pheno[keep.mask]
        biallelic.plot$id <- biallelic.plot$id[keep.mask]
        cdbg_all_plots(biallelic = biallelic.plot, triallelic = NULL,
                       genVars = NULL, treeInfo = treeInfo, config = config)
        rm(biallelic.plot)
        cleanMem()
    }
    
    return(list("biallelic" = biallelic, "config" = config, "treeInfo" = treeInfo))
}


## all_plots(biallelic = biallelic, triallelic = triallelic, config = config, treeInfo = treeInfo, genVars = genVars)


################################################################################################
##
## Adapted from BUGWAS_functions.R
##
## Changelog:
##   - calls gemma with option "-lmm 4" rather than "-lmm 2" in order
##    to get beta (weights in the linear model) and standard errors
##    rather than just lrt p-values.##
##   - optional cov.file argument, providing the name of a file which
##    contains covariates to be added to the linear model.
##
## Run GEMMA software on binary data.
## @XX: Binary variant patterns
## @relmatrix: Path to a relatedness matrix built from biallelic SNP data
## @pattern: For each variant, which pattern (row of XX) is it
## @ps: For each variant, what position is it in a reference genome (if mapped data)
##		Otherwise, an ID number
## @pheno.file: File containing phenotype to be tested (in same order as columns of XX)
## @maf: Minor allele frequency to test using GEMMA (Default: 0, all variants are tested)
## @prefix: Output file prefix
## @path: Path to where GEMMA is installed
## @dir: Working directory
## @process.results: Whether to process the p-value output of LMM or just return lambda and
## 					 the log likelihood under the null
## 					 Defaults to TRUE, which processes LMM p-value results
##
## Outputs:
## Likelihood ratio test results, the log likelihood under the null, lambda
################################################################################################

cdbg_run_lmm_bi <- function(XX = NULL,
                            relmatrix = NULL,
                            pattern = NULL,
                            ps = NULL,
                            pheno.file = NULL,
                            cov.file = NULL,
                            maf = NULL,
                            prefix = NULL,
                            path = NULL,
                            dir = NULL,
                            process.results = TRUE){
    
                                        #message("cdbg_run_lmm_bi")				   	
    
                                        # Output file names
    gen.output.file <- paste0(prefix, "_gemma_genfile.txt")
    snp.output.file <- paste0(prefix, "_gemma_snpfile.txt")

    n.pat <- nrow(XX)
    
    if(is.null(dim(XX))){
        XX <- matrix(XX,nrow=1)
    }
    
    ## Check if the variants are binary. Takes up a lot of memory, so we comment it for now.
    ## num.alleles <- apply(XX, 1, function(data) length(unique(data)))    
    ## if(length(which(num.alleles>2))!=0){
    ##     stop("\nError: function cdbg_run_lmm_bi requires binary variants\n")
    ## }

    
    ## gen.file <- cbind(paste0("pattern",1:n.pat),rep(1,n.pat),rep(0,n.pat),XX)
    ## write.table(gen.file, file = gen.output.file, row=F, col=F, sep="\t", quote=F)

    XX <- cbind(rep(1,n.pat),rep(0,n.pat), XX)
    write.table(XX, file = gen.output.file, row.names=paste0("pattern",1:n.pat), col.names=F, sep="\t", quote=F)
    
    snp.file <- cbind(paste0("pattern",1:n.pat), 1:n.pat, rep(24,n.pat))
    write.table(snp.file, file = snp.output.file, row=F, col=F, sep="\t", quote=F)

    ## system(paste0(path, " -g ", gen.output.file, " -p ", pheno.file, " -a ", snp.output.file,
    ##               " -k ", relmatrix," -lmm 4 -o ", prefix, "_lmmout_patterns"," -maf ", maf))

    if(!is.null(cov.file)){
        system(paste0(path, " -g ", gen.output.file, " -p ", pheno.file, " -a ", snp.output.file,
                      " -k ", relmatrix," -c ", cov.file," -lmm 4 -o ", prefix, "_lmmout_patterns"," -maf ", maf))
    }else{
        system(paste0(path, " -g ", gen.output.file, " -p ", pheno.file, " -a ", snp.output.file,
                      " -k ", relmatrix," -lmm 4 -o ", prefix, "_lmmout_patterns"," -maf ", maf))
    }
    
                                        ##message("LMM calculations completed successfully.")
    lmm.log <- paste0(dir, "/output/", prefix, "_lmmout_patterns.log.txt")
    
                                        ##message(paste(c("GEMMA log file:", lmm.log), collapse=""))
    
    lambda.lognull <- bugwas:::extract_lambda_lognull(lmm.log)
    
                                        ##message("Lambda estimations extracted successfully.")
                                        ##message(paste(c("Extracted lambda: ", lambda.lognull), collapse=""))
    
    if(process.results == FALSE){
        
        return(lambda.lognull)
        
    } else {
        
        assocFile = paste0(dir, "/output/", prefix, "_lmmout_patterns.assoc.txt")
                                        ##message(paste(c("assocFile:", assocFile), collapse=" "))
        lmm <- read.table(assocFile, header=T, sep="\t", as.is=T)

                                        #message(paste(c("Header:", names(lmm)), collapse=" "))
        
        LH1 <- lmm$logl_H1
        
        
        D <- sapply(LH1, bugwas:::get_deviance, lognull=lambda.lognull["lognull"], USE.NAMES=FALSE)
        pvals <- pchisq(as.numeric(D), 1, low=F)
        negLog10 <- -log10(pvals)	
        
                                        ##message("flag1")
        
        ## !!! lmm$ps is not the same thing as ps. It is indexing
        ## !!! unique patterns in the output file of gemma so its
        ## !!! length is at most the number of unique
        ## !!! patterns. The following operation therefore simply
        ## !!! restrict the "pattern" variable to the subset of
        ## !!! unique patterns which were tested by gemma.  ##
        ## !!! Important: since these patterns are arbitrarilly
        ## !!! named 1:nrow(XX), patterns should contain row
        ## !!! numbers in XX, not pattern IDs.
        m <- match(pattern, lmm$ps) 
        pattern <- pattern[!is.na(m)]
        ## This restricts ps (which has one value per variant) to
        ## the set of variants corresponding to tested patterns.
        ps <- ps[!is.na(m)]
        
                                        ##message("flag2")
        
        lmm <- cbind(lmm, negLog10)
        lmm <- lmm[match(pattern,lmm$ps), ]
        lmm$ps <- ps
        
                                        ##message("flag3")
        
        write.table(lmm, file = paste0(prefix, "_lmmout_allSNPs.txt"), sep="\t",
                    row=F, col=T, quote=F)
        
        return(list("lmm" = lmm, "lognull" = lambda.lognull["lognull"],
                    "lambda" = lambda.lognull["lambda"]))
    }
    
}

################################################################################################
## Adapted from BUGWAS_functions.R
## Changelog:
##   - calls cdbg_run_lmm_bi rather than bugwas:::run_lmm_bi
##   - calls cdbg_get_correlations rather than bugwas:::get_correlations
## Functions for biallelic data
################################################################################################

cdbg_get_biallelic <- function(logreg.bi = NULL,
                               XX.all = NULL,
                               XX = NULL,
                               lmm.bi = NULL,
                               lognull = NULL,
                               lambda = NULL,
                               relmatrix = NULL,
                               pheno.file = NULL,
                               cov.file = NULL,
                               maf = NULL,
                               gem.path = NULL,
                               output.dir = NULL,
                               prefix = NULL,
                               run.lmm = NULL,
                               XX.ID = NULL,
                               pca = NULL,
                               npcs = NULL){

    if(!is.null(logreg.bi)){
        logreg.bi <- read.table(logreg.bi, header=T, sep="\t", as.is=T)
        if(nrow(logreg.bi)!=length(XX.all$pattern)){
            cat(paste0("logreg.bi = ", nrow(logreg.bi), " variants"),"\n")
            cat(paste0("gen = ", length(XX.all$pattern), " biallelic variants"),"\n")
            stop("\nError: number of variants in logreg.bi does not match number of biallelic SNPs in gen")
        }
        if(any(is.na(match(logreg.bi$ps, XX.all$ps)))){
            stop("\nError: variant positions/IDs do not match between logreg.bi and biallelic variants in gen\n")
        }
        if(any(logreg.bi$ps != XX.all$ps)){
            m <- match(XX.all$ps, logreg.bi$ps)
            logreg.bi <- logreg.bi[m,]
        }
    }
    
    
    if(!is.null(lmm.bi)){
        lmm.bi <- read.table(lmm.bi, header=T, sep="\t", as.is=T)
        if(nrow(lmm.bi)!=length(XX.all$pattern)){
            cat(paste0("lmm.bi = ", nrow(lmm.bi), " variants"),"\n")
            cat(paste0("gen = ", length(XX.all$pattern), " biallelic variants"),"\n")
            stop("\nError: number of variants in lmm.bi does not match number of biallelic SNPs in gen")
        }
        if(any(is.na(match(lmm.bi$ps, XX.all$ps)))){
            stop("\nError: variant positions/IDs do not match between lmm.bi and biallelic variants in gen\n")
        }
        if(any(lmm.bi$ps != XX.all$ps)){
            m <- match(XX.all$ps, lmm.bi$ps)
            lmm.bi <- lmm.bi[m,]
        }
        
        if(is.null(lognull) | is.null(lambda)){
            lambda.lognull <- cdbg_run_lmm_bi(XX = XX.all$XX[1,], relmatrix = relmatrix,
                                              pheno.file = pheno.file, cov.file = cov.file, maf = maf,
                                              prefix = paste0(prefix, "_getlognull"), path = gem.path,
                                              dir = output.dir, process.results = FALSE)
            if(is.null(lambda)){
                lambda <- as.numeric(lambda.lognull$lambda)
                cat(paste0("## Lambda = ", lambda), file = paste0(prefix,"_logfile.txt"),
                    sep="\n", append = TRUE)
            }
            if(is.null(lognull)){
                lognull <- as.numeric(lambda.lognull$lognull)
                cat(paste0("## Log likelihood under the null = ", lognull),
                    file = paste0(prefix,"_logfile.txt"), sep="\n", append = TRUE)
            }
        }
        
    } else if(is.null(lmm.bi) & run.lmm){	
        lmm.bi <- cdbg_run_lmm_bi(XX = XX.all$XX, relmatrix = relmatrix, pattern = XX.all$pattern,
                                  ps = XX.all$ps, pheno.file = pheno.file, cov.file = cov.file, maf = maf,
                                  prefix = paste0(prefix, "_biallelic"), path = gem.path, dir = output.dir)
        
        lognull <- as.numeric(lmm.bi$lognull)
        cat(paste0("## Log likelihood under the null = ", lognull),
            file = paste0(prefix,"_logfile.txt"), sep="\n", append = TRUE)
        lambda <- as.numeric(lmm.bi$lambda)
        cat(paste0("## Lambda = ", lambda), file = paste0(prefix,"_logfile.txt"),
            sep="\n", append = TRUE) 
        lmm.bi <- lmm.bi$lmm
    }

    if(!is.null(pca)){
        message('Computing correlations between patterns and PCs')
        cor.XX <- cdbg_get_correlations(XX = XX, pca = pca$x, npcs = npcs, id = XX.ID)
        cleanMem()
    }else{
        cor.XX <- NULL
    }

    
    return(list("logreg.bi" = logreg.bi, "lmm.bi" = lmm.bi, "lognull" = lognull,
                "lambda" = lambda, "cor.XX" = cor.XX))

}

################################################################################################
## Get correlations between variant patterns and principal components.
## @XX: variant matrix scaled by how many variants each pattern represents
## @pca: principal component analysis
## @npcs: number of principal components (= number of individuals)
##
## Outputs:
## which.pc: for each variant pattern, which PC it is most correlated to
## max.cor.pc: for each variant pattern, the maximum correlation value
## Adapted from BUGWAS_functions.R
##
## Changelog:
## Using more efficient (memory+cpu) technique for which.pc and max.cor.pc computation.
################################################################################################


cdbg_get_correlations <- function (XX = NULL, 
                                   pca = NULL,
                                   npcs = NULL,
                                   id = NULL,
                                   all.cor  = FALSE){

    cor.XX.pca <- cor(XX,pca[, 1:npcs])
    cleanMem()
    cor.XX.pca[is.na(cor.XX.pca)] = 0
    which.pc <- max.col(abs(cor.XX.pca), ties.method='first')
    max.cor.pc <- abs(cor.XX.pca[cbind(1:nrow(cor.XX.pca), which.pc)])
    names(which.pc) <- names(max.cor.pc) <- rownames(cor.XX.pca)
    if(all.cor){
        return(list("which.pc" = which.pc, "max.cor.pc" = max.cor.pc, "all.cor.pc" = cor.XX.pca))
    }else{
        return(list("which.pc" = which.pc, "max.cor.pc" = max.cor.pc))
    }
    
}

################################################################################################
## Changelog:
##
## - Use pca instead of svd.XX (original function used both but they
##  contain the same information, taking up space and requiring
##  unnecessary computation).
## - Call cdbg_ridge_regression instead of bugwas:::ridge_regression.
## - Call cdbg_get_wald_input instead of bugwas:::get_wald_input.
##
## Do Wald test
################################################################################################

cdbg_wald_test <- function(y = NULL,
                           XX = NULL,
                           lambda = NULL,
                           XX.all = NULL,
                           prefix = NULL,
                           npcs = NULL,
                           pca = NULL){
    

    fit.lmm <- cdbg_ridge_regression(y, XX, pca=pca,
                                     lambda_init=as.numeric(lambda)/sum(XX.all$bippat),
                                     maximize=FALSE, skip.var=TRUE)
    
    ## Fit the grand null model
    
    fit.0 <- lm(y~1)
    
    ## LRT for the LMM null vs grand null
    LRTnullVgrand <- -log10(pchisq(2*(fit.lmm$ML - as.numeric(logLik(fit.0))), 1, low=F)/2)
    cat(paste0("## LRT for the LMM null vs grand null = ", LRTnullVgrand),
    	file = paste0(prefix, "_logfile.txt"), sep="\n", append = TRUE)
    
    ## Heritability
    fit.lmm.ypred <- XX %*% fit.lmm$Ebeta
    cat(paste0("## Heritability (R^2) = ", cor(fit.lmm.ypred,y, use='complete.obs')^2),
    	file=paste0(prefix, "_logfile.txt"), sep="\n", append=TRUE)
    
    ## Get full posterior covariance matrix for Bayesian Wald Test
    ## Need the full posterior covariance matrix for the Bayesian Wald test,
    ## to get the posterior uncertainty for each point estimate
    
    wald_input <- cdbg_get_wald_input(fit.lmm = fit.lmm, pca = pca,
                                 y = y)
    
    ## Bayesian Wald Test
    pca.bwt <- wald_input$Ebeta^2/diag(wald_input$Vbeta)
    
    p.pca.bwt <- -log10(exp(1))*pchisq(pca.bwt, 1, low=F, log=T)
    cat(paste0("## Bayesian Wald Test for PCs range = ", paste(range(p.pca.bwt), collapse=" ")),
    	file=paste0(prefix, "_logfile.txt"), sep="\n", append=TRUE)
    write.table(p.pca.bwt, file = paste0(prefix, "_Bayesian_Wald_Test_negLog10.txt"),
    	        sep="\t", row=T, col = F, quote=F)
    
    ## Get order of PCs by Wald test results
    signif_cutoff <- -log10(0.05/npcs)
    pc_order <- bugwas:::get_pc_order(p.pca.bwt = p.pca.bwt, signif_cutoff = signif_cutoff)
    ## Predict phenotype using effect sizes
    ## effect <- t(t(XX) * as.vector(fit.lmm$Ebeta))
    ## pred <- rowSums(effect)
    
    pred <- wald_input$XE
    
    return(list("pc_order" = pc_order, "p.pca.bwt" = p.pca.bwt, "pred" = pred,
                "signif_cutoff" = signif_cutoff))
    
}

################################################################################################
##
## Changelog:
## - Use pca instead of pca and svd (redundant since scaled XX was column centered).
## - Use algebraic simplifications (caution: assuming npcs = ncol(pca$rotation))
##
## Get Bayesian Wald Test inputs.
## @fit.lmm: ridge regression results
## @pca: principal component analysis
##
## Outputs:
## pca.Ebeta
## pca.Vbeta
################################################################################################

cdbg_get_wald_input <- function(fit.lmm = NULL,
                                pca = NULL,
                                y = NULL){

    ## Get full posterior covariance matrix for Bayesian Wald Test
    ## Need the full posterior covariance matrix for the Bayesian Wald test,
    ## to get the posterior uncertainty for each point estimate
    lambda = fit.lmm$lambda_MLE
    ## Cstar = diag(lambda * pca$d^2 / (lambda * pca$d^2 + 1))
    ## For the null model, the posterior mean and variance (may be slow!)
    ## astar = t(y) %*% y - t(y) %*% XX %*% fit.lmm$Ebeta

    ## Deal with NA by restricting post.mean computation to phenotyped
    ## samples
    missing.y <- is.na(y)
    pca.Ebeta <- diag(1/(1/(fit.lmm$prefered.lambda) + pca$d^2)) %*% t(pca$x[!missing.y, ]) %*% y[!missing.y]
    XE <- pca$x %*% pca.Ebeta
    astar <- crossprod(y[!missing.y]) - t(y[!missing.y]) %*% XE[!missing.y, ]
    ## csqxty <-  diag(1/sqrt(1/(fit.lmm$prefered.lambda) + pca$d^2)) %*% t(pca$x) %*% y
    ## astar = crossprod(y) - crossprod(csqxty)
    dstar = sum(!missing.y) #length(y)
    tau = as.numeric((dstar-2)/astar)
    
    ## rotation = t(pca$rotation[,1:npcs])
    ## Based on the PCA rotations of raw genetic diversity
    ## pca.Ebeta = rotation %*% fit.lmm$Ebeta
    ## pca.Ebeta <- diag(1/sqrt(1/(fit.lmm$prefered.lambda) + pca$d^2)) %*% csqxty
    
    ## rtr = tcrossprod(rotation, rotation); # Should be n (sample size) by n
    ## rv = rotation %*% svd.XX$v; # Should be n by n
    ## rv = crossprod(pca$rotation); # Should be n by n
    ## rtr <- rv <- diag(nrow(rotation)) # Since rotation is an orthogonal matrix
    ## pca.Vbeta = rv %*% Cstar; # Should be n by n
    ## pca.Vbeta = tcrossprod(pca.Vbeta,rv); # Should be n by n
    ## pca.Vbeta = lambda/tau * (rtr - pca.Vbeta); # Should be n by n
    pca.Vbeta <- lambda/tau * diag(1 - (lambda * pca$d^2 / (lambda * pca$d^2 + 1)))
    
    return(list("Ebeta" = pca.Ebeta, "Vbeta" = pca.Vbeta, "XE"=XE))    
}

################################################################################################
## Changelog:
##
## - Call cdbg_get_correlation instead of bugwas:::get_correlation.
##
################################################################################################

cdbg_get_tree <- function (phylo = NULL, prefix = NULL, XX.ID = NULL, pca = NULL, 
                           npcs = NULL, allBranchAndPCCor = FALSE) 
{
    tree <- ape::read.tree(phylo)
    if (any(is.na(match(XX.ID, tree$tip.label)))) {
        stop("\nError: Phylogeny sample names do not match the sample names of the SNP data\n")
    }
    tree <- phangorn::midpoint(tree)
    ape::write.tree(tree, paste0(prefix, "_midpointrooted_tree.txt"))
    tree <- ape::read.tree(paste0(prefix, "_midpointrooted_tree.txt"))
    treepat <- bugwas:::tree2patterns(tree = tree, tiporder = XX.ID)
    mtp <- treepat$pat
    message("Retrieve all correlations between branches and PCs: ", 
        allBranchAndPCCor)
    cor.tree <- cdbg_get_correlations(XX = mtp, pca = pca$x, npcs = npcs, 
        id = XX.ID, all.cor = allBranchAndPCCor)
    if (allBranchAndPCCor) {
        branchPCCorFileName = paste(prefix, "allBranchAndPCCor.txt", 
            sep = "_")
        write.table(signif(cor.tree$all.cor.pc, digits = 3), 
            branchPCCorFileName, row.names = T, col.names = T, 
            quote = F, sep = "\t")
        ape::write.tree(treepat$labelled.tree, paste0(prefix, 
            "_node_labelled_tree.txt"))
    }
    return(list(tree = tree, pattern = treepat, cor.tree = cor.tree))
}
