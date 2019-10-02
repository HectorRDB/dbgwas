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
## files all_plots.R, createAllplots.R and
## snpManhattanPlot.R. Appending "cdbg" to the modified functions.

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

## #' Generates all plots.
## #' 
## #' This function generates all the plots
## #' @param biallelic A list called 'biallelic' created from the lin_loc function
## #' @param triallelic A list called 'triallelic' created from the lin_loc function
## #' @param genVars A list called 'genVars' created from the lin_loc function
## #' @param treeInfo A list called 'treeInfo' created from the lin_loc function
## #' @param config A list called 'config' created from the lin_loc function
## #' @keywords plot
## #' @export
## #' @examples
## #' data <- lin_loc(gen = gen, pheno = pheno, phylo = phylo, 
## #'  prefix = prefix, gem.path = gem.path)
## #' all_plots(biallelic = data$biallelic, triallelic = data$triallelic, 
## #' 	genVars = data$genVars, treeInfo = data$treeInfo, config = data$config)

cdbg_all_plots = function(biallelic = NULL,
                     triallelic = NULL,
                     genVars = NULL,
                     treeInfo = NULL,
                     config = NULL){
  
  
  cdbg_createAllPlots(prefix = config$prefix,
                 genVars = genVars,
                 cutoffCor = config$cutoffCor,
                 npcs = biallelic$npcs,
                 phenotype = biallelic$pheno,
                 pca = biallelic$pca,
                 fit.lm = biallelic$logreg,
                 fit.lmm = biallelic$lmm,
                 fit.lm.tritetra = triallelic$logreg,
                 fit.lmm.tritetra = triallelic$lmm,
                 gemma = biallelic$lmm,
                 gemma.tritetra = triallelic$lmm,
                 ipat = biallelic$pattern,
                 ipat.snps = triallelic$pattern,
                 pred2 = biallelic$pred,
                 cor.XX = biallelic$cor.XX,
                 cor.tritetra = triallelic$cor.tritetra,
                 which.pc =  biallelic$cor.XX$which.pc,
                 max.cor.pc = biallelic$cor.XX$max.cor.pc,
                 which.pc.tritetra = triallelic$cor.tritetra$which.pc,
                 max.cor.pc.tritetra = triallelic$cor.tritetra$max.cor.pc,
                 which.mtp.pc = treeInfo$cor.tree$which.pc,
                 max.mtp.cor.pc = treeInfo$cor.tree$max.cor.pc,
                 XX.comid = biallelic$id,
                 pc_order = biallelic$pc_order,
                 o = biallelic$pc_order$pc_order,
                 pc.lim = biallelic$pc_order$pc.lim,
                 p.pca.bwt = biallelic$p.pca.bwt,
                 signifCutOff = config$signif_cutoff,
                 bippat = biallelic$bippat,
                 snppat = triallelic$snppat,
                 tree = treeInfo$tree,
                 treepat = treeInfo$pattern);
  
}


cdbg_createAllPlots = function(prefix = NULL, 
genVars = NULL, cutoffCor = NULL,  npcs = NULL, phenotype = NULL,
                          pca = NULL, fit.lm = NULL, fit.lmm = NULL, fit.lm.tritetra = NULL,
                          fit.lmm.tritetra = NULL, gemma = NULL, gemma.tritetra = NULL, 
                          ipat = NULL, ipat.snps = NULL, pred2 = NULL, cor.XX = NULL, cor.tritetra = NULL,
                          which.pc = NULL, max.cor.pc = NULL, 
                          which.pc.tritetra = NULL, max.cor.pc.tritetra = NULL,
                          which.mtp.pc = NULL, max.mtp.cor.pc = NULL,
                          XX.comid = NULL, pc_order = NULL, o = NULL, pc.lim = NULL, p.pca.bwt = NULL,
                          signifCutOff = NULL, bippat = NULL, snppat = NULL,
                          tree = NULL, treepat = NULL){
              
  
  # Get a random sample of colours to use for all plots, equal to the number of significant PCs
  colourPalette = bugwas:::getColourPalette(p.pca.bwt = p.pca.bwt, signifCutOff = signifCutOff, pc.lim = pc.lim)
  
  sampleCount = length(phenotype)
  m = match(o[pc.lim], which.mtp.pc)
    
    ##Bayesian Wald test for genome-wide PCs
    p.genomewidepc = bugwas:::.testGenomeWidePCs(prefix = prefix,
                                                 pc.lim = pc.lim,
                                                 pca = pca,
                                                 bippat = bippat,
                                                 ipat = ipat,
                                                 o = o)
  message("Bayesian Wald test for genome-wide PCs has been completed successfully.")
    
  #The barplot for the Bayesian wald test for genome-wide PCs
  bugwas:::.BayesianWaldTestPCsBarplot(prefix = prefix,
                              p.pca.bwt = p.pca.bwt,
                              colourPalette = colourPalette,
                              o = o,
                              m = m,
                              p.genomewidepc = p.genomewidepc,
                              pc.lim = pc.lim)
   message("The barplot for the Bayesian wald test for genome-wide PCs has been completed successfully.")
       
    ## snpColours = cdbg_getSNPColours(sampleCount = sampleCount,
    ##                                colourPalette = colourPalette,
    ##                                colouredPCs =  o[1:20],
    ##                                which.pc = which.pc,
    ##                                ipat = ipat,
    ##                                max.cor.pc = max.cor.pc,
    ##                                cutoffCor = cutoffCor,
    ##                                which.pc.tritetra = which.pc.tritetra,
    ##                                ipat.snps = ipat.snps,
    ##                                max.cor.pc.tritetra = max.cor.pc.tritetra,
    ##                                cor.XX = cor.XX)
    
  ## bipCount = 0
  ## ttpCount = 0
  ## if(!is.null(fit.lmm)){
  ## 	bipCount = nrow(fit.lmm)
  ## }
  ## if(!is.null(fit.lmm.tritetra)){
  ## 	ttpCount = nrow(fit.lmm.tritetra)
  ## }
  
  ## snpType = c(rep(1,bipCount),rep(2,ttpCount))

    ## Not sure it makes sense when using unitigs rather than SNPs: no
    ## linear chromosome structure. Skip for now.

    ## #The Manhattan plot for SNP GWAS using logistic regerssion
    ## logregPvalues = NULL
    ## logregPos = NULL
    ## if(!is.null(fit.lm)){
        
    	##logregPos = c(fit.lm$ps,fit.lm.tritetra$ps)
        ##logregPvalues = c(fit.lm$negLog10,fit.lm.tritetra$negLog10)
        
        ## bugwas:::.manhattanPlot(prefix = paste0(prefix,"_ManhattanRawPvalues"),
        ##                         snpPos = logregPos,
        ##                         pValues = logregPvalues,
        ##                         snpType = snpType,
        ##                         main = "Logistic Regression SNPs Manhattan Plot",
        ##                         col = c(snpColours$bip, snpColours$ttp))
        
        ## message("The Manhattan plot for SNP GWAS using logistic regerssion has been completed successfully.")
    
    ##}
    
    
    ## #The Manhattan plot for SNP GWAS using LMM
    ##lmmPos = c(fit.lmm$ps,fit.lmm.tritetra$ps)
    ##lmmPvalues = c(-log10(fit.lmm$p_lrt),-log10(as.numeric(fit.lmm.tritetra$pvals)))
    
    ## bugwas:::.manhattanPlot(prefix = paste0(prefix,"_ManhattanLMMPvalues"),
    ##               snpPos = lmmPos,
    ##               pValues = lmmPvalues,
    ##               snpType = snpType,
    ##               main = "LMM SNPs Manhattan Plot",
  ##               col = c(snpColours$bip, snpColours$ttp))
    
    ## message("The Manhattan plot for SNP GWAS using LMM has been completed successfully.")
  

    ## Skip this too for now since we are never doing a logistic
    ## regression. Putting it back requires uncommenting the
    ## computation of lmmPvalues and logregPvalues above, as well as
    ## snpColours computation.
    
  #The plot of logistic regression P-values vs. LMM P-values for SNP GWAS
  ## if(!is.null(logregPvalues)){
  ##   bugwas:::.logregVsLMM(prefix = prefix,
  ##                        col = c(snpColours$bip, snpColours$ttp),
  ##                        logregPvalues = logregPvalues,
  ##                        lmmPvalues =  lmmPvalues,
  ##                        snpType = snpType,
  ##                        pcOrder = o,
  ##                        pc.lim = pc.lim,
  ##                        colourPalette = colourPalette)
                         
  ##    message("The plot of logistic regression P-values vs. LMM P-values for SNP GWAS has been completed successfully.")                    
  ## }
  
  
  
  # Plot the individuals by their top two significant additive PCs
    .plotIndividualBy2PCs(pc1 = o[1], pc2 = o[2],
                         pc1.scores = pca$x[,o[1]], pc2.scores = pca$x[,o[2]],
                         prefix= prefix, phenotype)
                       
    message("The reduced space plot of the sample on the top two significant additive PCs been completed successfully.")                    
                       
  #The plot with true and predicted phenotype on the tree
    .trueAndPredPhenoOnTreePlot(prefix = prefix, tree = tree, which.mtp.pc = unlist(which.mtp.pc), #Check with SGE
                                max.mtp.cor.pc = max.mtp.cor.pc, cutoffCor = cutoffCor, treepat = treepat,
                                pcOrder = o, p.genomewidepc = p.genomewidepc, phenotype = phenotype, 
                                XX.comid = XX.comid, colourPalette = colourPalette,
                                pc.lim = pc.lim, pred2 = pred2)
                             
    message("The plot with true and predicted phenotype on the tree has been completed successfully.")                         
  
    ## The plots of PCs loadings    
    bugwas:::.pcLoadingsPlot(prefix = prefix, pca = pca, pc.lim = pc.lim, 
                             pcOrder = o, ipat = ipat, bippat = bippat, 
                             bipPos = fit.lmm$ps)

  message("The plots of PCs loadings have been completed successfully.")
  
  ## To run for all SNPs (biallelic and tri and tetra allelic)
  new.pat <- ipat.snps
  new.pat <- sapply(new.pat, function(x, pat){ x + length(pat)}, pat = bippat)
  new.pat <- c(ipat, unlist(new.pat))
  
  #The Manhattan plot organised by PCs for SNP GWAS
    .plot_pc_manhattan(o = o, 
                       which.pc = c(cor.XX$which.pc, cor.tritetra$which.pc), 
                       pattern = new.pat, 
                       p.pca.bwt = p.pca.bwt, 
                       pc.lim = pc_order$pc.lim, 
                       negLog10 = c(fit.lmm$negLog10, fit.lmm.tritetra$negLog10), 
                       pat.weight = c(bippat, snppat), 
                       prefix=paste0(prefix,"_SNPs"),
                       colourPalette = colourPalette,
                       npcs = npcs)
    message("The Manhattan plot organised by PCs for SNP GWAS has been completed successfully.")                
  
  
  #The plots for general genetic variants
  if(!is.null(genVars)){
    bugwas:::.genVarPlots(genVars = genVars, o = o, p.pca.bwt = p.pca.bwt, pc.lim = pc.lim, 
                colourPalette = colourPalette, prefix = prefix, npcs = npcs, 
                sampleCount = sampleCount,cutoffCor = cutoffCor)
     message("The plots for general genetic variants have been completed successfully.")            
    
   
  }
  
}

## Get the colouring of the SNPs for manhattanp plot
cdbg_getSNPColours = function(sampleCount = NULL,
                         colourPalette = NULL,#colourPalette
                         colouredPCs = NULL, # o[1:20]
                         which.pc = NULL,
                         ipat = NULL,
                         max.cor.pc = NULL,
                         cutoffCor = NULL,
                         which.pc.tritetra = NULL,
                         ipat.snps = NULL,
                         max.cor.pc.tritetra = NULL,
                         cor.XX = NULL){
 
    ##COL = rep("grey50", sampleCount)  
    ##COL[colouredPCs] = colourPalette
    ##COL = COL[which.pc][ipat]
    COL = rep("grey50", length(ipat))
    col.keep <- (colourPalette != "grey50")
    colourPalette <- colourPalette[col.keep]
    colouredPCs <- colouredPCs[col.keep]
    names(colourPalette) <- colouredPCs  
    col.mask <- (ipat %in% names(which.pc[which.pc %in% colouredPCs]))
    COL[col.mask] <- colourPalette[as.character(which.pc[ipat[col.mask]])]
    if(cutoffCor > 0){
        COL[max.cor.pc[ipat] < cutoffCor] = "grey50"
    }
  
  return(list(bip = COL ,ttp = NULL))
  
}

###################################################################
## Plot manhattan plot ordered on the x-axis by PCs
## @o: PCs in order of significance
## @which.pc: For each pattern, which PC is it most correlated to
## @pattern: For each variant, what pattern it is
## @p.pca.bwt: Bayesian Wald test results for each PC
## @pc.lim: Which PCs are significant by the BWT
## @lmm: LMM results matrix
## @pat.weight: How many variants does each pattern represent
## @prefix: output prefix
## Outputs:
## Manhattan plot with the x-axis ordered by PCs
##
## Changelog:
## Deal with infinite negLog10 which sometimes arise when pvalues
## are 0.
###################################################################
.plot_pc_manhattan <- function(o = NULL,
                               which.pc = NULL,
                               pattern = NULL,
                               p.pca.bwt = NULL,
                               pc.lim = NULL,
                               negLog10 = NULL,
                               pat.weight = NULL,
                               prefix = NULL,
                               colourPalette = NULL,
                               npcs = NULL){
    ## Find which PCs have a variant that is most correlated to it
    m <- match(o,unique(which.pc))
    o.pats <- o[which(is.na(m)==FALSE)]
    ## Subset the Bayesian Wald Test -log10(p) to those PCs which have kmers most correlated to them
    p.pca.bwt.pats <- p.pca.bwt[o.pats]
    
    
    num.variants <- sapply(o.pats, function(x, pat.weight=NULL, which.pc = NULL)
        sum(pat.weight[which(which.pc==x)]),pat.weight = pat.weight,
        which.pc = which.pc, USE.NAMES=FALSE)
    
    cols.pcs <- rep_len(c("#5a5a5a","#c6c6c6"), length.out = length(o.pats))
    if(!is.null(pc.lim)){
        cols.pcs[1:length(pc.lim)] <- colourPalette[1:length(pc.lim)]
    }
    m <- match(which.pc,o.pats)
    cols <- cols.pcs[m[pattern]]
    ##cols[which(max.cor.pc[which(which.pc==o.pats[i])]<0.3)]="grey50"
    
    
    pos.gap <- sapply(num.variants, function(x) 10000/x, USE.NAMES=FALSE)
    pos.gap[num.variants > 10000] <- 1
    if(!is.null(pc.lim)){
        pos.gap[is.na(match(o.pats, o[pc.lim]))] <- pos.gap[is.na(match(o.pats, o[pc.lim]))]/ (npcs/10)
    } else {
        pos.gap[21:length(pos.gap)] <- pos.gap[21:length(pos.gap)] / (npcs/10)
    }
    
    pos <- rep(0,length(negLog10))
    max.pos <- 0
    plot.lines <- matrix(rep(0, length(o.pats)*2), ncol=2)
    
    which.pc <- which.pc[pattern]
    
    for(i in 1:length(o.pats)){
        
        s = seq(from = (max.pos+pos.gap[i]), by = pos.gap[i], length.out = num.variants[i])
        if(num.variants[i]>1){
            pos[which(which.pc==o.pats[i])] <- sample(s, num.variants[i], replace=FALSE)
        } else {
            pos[which(which.pc==o.pats[i])] <- s
        }
        max.pos <- max(s)
        plot.lines[i,] <- c((max.pos-pos.gap[i]*num.variants[i]+((num.variants[i]*0.15)*pos.gap[i])), (max.pos-((num.variants[i]*0.15)*pos.gap[i])))
        
    }
    
    neg.log.pv <- c(p.pca.bwt.pats, negLog10)    
    ## Replace infinite values by maximal finite value. Keep track of
    ## their indices to mark them by a special symbol in the plot.    
    ymax = max(neg.log.pv[!is.infinite(neg.log.pv)])
    pca.inf <- which(is.infinite(p.pca.bwt.pats))
    negLog10.inf <- which(is.infinite(negLog10))
    p.pca.bwt.pats[pca.inf] <- negLog10[negLog10.inf] <- ymax
    
    bugwas:::.pl2(paste0(prefix,"_PC_manhattan"),{
        plot(x = pos, y = negLog10, col = cols, pch = c(1,18)[1 + (1:length(pos) %in% negLog10.inf)],
             xlab = "", ylab = "LMM -log10(p)",
             ylim = c(0, ymax), xaxt = "n")
        
        for(i in 1:length(o.pats)){
            lines(x = c(plot.lines[i, ]), y = c(p.pca.bwt.pats[i], p.pca.bwt.pats[i]),
                  type = "l", col = cols.pcs[i], lty = (1 + (i %in% pca.inf)), lwd=2)
            if(i<=20){
                text(x = c(plot.lines[i, 1]+((plot.lines[i, 2]-plot.lines[i, 1])/2)),
                     y = c(p.pca.bwt.pats[i]+(max(negLog10)/90)), labels = o.pats[i],
                     col = cols.pcs[i], font = 2, cex = 0.9)
            }
        }
    })
    

}

###################################################################
#' Plot of the sample on the first two principal components
#' 
#' This function generates a plot of the sample on the first two
#' principal components.
#'
#' Changelog: deal with continuous phenotypes
###################################################################

.plotIndividualBy2PCs <- function (pc1 = NULL, pc2 = NULL, pc1.scores = NULL, pc2.scores = NULL, 
    prefix = NULL, phenotype = NULL) 
{

    if(length(unique(phenotype[!is.na(phenotype)])) == 2){
        pheno.colors <- c("blue", "red")[1 + phenotype]
    }else{
        ## Use n equally spaced breaks to assign each value to n-1 equal sized bins 
        ii <- cut(phenotype, breaks = seq(min(phenotype), max(phenotype), len = 100), 
                  include.lowest = TRUE)
        ## Use bin indices, ii, to select color from vector of n-1 equally spaced colors
        pheno.colors <- colorRampPalette(c("#aaaaff", "navy"))(99)[ii]
    }
    bugwas:::.pl(paste0(prefix, "_indiv_first2signifPCs"), {
        plot(pc1.scores + rnorm(length(pc1.scores), 0, sd(pc1.scores)/10), 
            pc2.scores + rnorm(length(pc2.scores), 0, sd(pc2.scores)/10), 
            col = pheno.colors, lwd = 0.5, 
            cex = 1.5, xlab = paste0("PC ", pc1), ylab = paste0("PC ", 
                pc2))
    })
}

###################################################################
#' 
#' This function plots the true and predicted phenotypes on a tree.
#' 
#' Changelog: deal with continuous phenotypes
###################################################################

.trueAndPredPhenoOnTreePlot <- function (prefix = NULL, tree = NULL, which.mtp.pc = NULL, max.mtp.cor.pc = NULL, 
    cutoffCor = NULL, treepat = NULL, pcOrder = NULL, p.genomewidepc = NULL, 
    phenotype = NULL, XX.comid = NULL, colourPalette = NULL, 
    pc.lim = NULL, pred2 = NULL) 
{
    m = match(pcOrder[pc.lim], which.mtp.pc)
    branch.col.pc = bugwas:::.getBranchColors(treepat = treepat, tree = tree, 
        phenotype = phenotype, pcOrder = pcOrder, colourPalette = colourPalette, 
        which.mtp.pc = which.mtp.pc, max.mtp.cor.pc = max.mtp.cor.pc, 
        cutoffCor = cutoffCor)
    tree.eq = bugwas:::.get.tree.eq(tree)
    pred3 = (pred2 - min(phenotype - mean(phenotype, na.rm=TRUE), na.rm=TRUE))/diff(range(phenotype - 
        mean(phenotype, na.rm=TRUE), na.rm=TRUE))
    edgeLabels = rep("", length(tree$edge))
    edgeLabelsCOL = rep("white", length(tree$edge))
    add.signif = NULL
    if (length(pc.lim) > 0) {
        sigInd = bugwas:::.getSigIndicators(pc.lim = pc.lim, pcOrder = pcOrder, 
            p.genomewidepc = p.genomewidepc, tree.eq = tree.eq, 
            branch.col.pc = branch.col.pc, colourPalette = colourPalette)
        edgeLabels = sigInd$edgeLabels
        edgeLabelsCOL = sigInd$edgeLabelsCOL
        add.signif = sigInd$add.signif
    }
    predCOLS = colorRamp(c("grey", "black"))
    predCOLS = predCOLS(bugwas:::.st(pred3))
    predCOLS = rgb(predCOLS, maxColorValue = 256)
    .trueAndPredPhenoOnTree(prefix = prefix, phenotype = phenotype, 
        tree.eq = tree.eq, branch.col.pc = branch.col.pc, tree = tree, 
        XX.comid = XX.comid, predCOLS = predCOLS, edgeLabels = edgeLabels, 
        edgeLabelsCOL = edgeLabelsCOL, add.signif = add.signif, 
        pc.lim = pc.lim, colourPalette = colourPalette, matchedPCs = m)
}

.trueAndPredPhenoOnTree <- function (prefix = NULL, phenotype = NULL, tree.eq = NULL, branch.col.pc = NULL, 
    tree = NULL, XX.comid = NULL, predCOLS = NULL, edgeLabels = NULL, 
    edgeLabelsCOL = NULL, add.signif = NULL, pc.lim = NULL, colourPalette = NULL, 
    matchedPCs = NULL) 
{
    bugwas:::.pl(paste0(prefix, "_tree_branchescolouredbyPC"), {
        n = length(phenotype)
        par(oma = c(1, 1, 1, 5))
        if(length(unique(phenotype[!is.na(phenotype)])) == 2){
            pheno.colors <- c("gray", "black")[1 + phenotype - min(phenotype)][match(tree$tip.label, XX.comid)]
        }else{
            ## Use n equally spaced breaks to assign each value to n-1 equal sized bins 
            ii <- cut(phenotype, breaks = seq(min(phenotype), max(phenotype), len = 100), 
                  include.lowest = TRUE)
            ## Use bin indices, ii, to select color from vector of n-1 equally spaced colors
            pheno.colors <- colorRampPalette(c("grey", "black"))(99)[ii]
        }
        ape::plot.phylo(tree.eq, type = "fan", show.tip.label = TRUE, 
                        edge.col = branch.col.pc, tip.col = pheno.colors, lwd = 3, adj = 0.5, lab4ut = "axial")
        ## cdbg_plot_phylo(tree.eq, type = "fan", show.tip.label = TRUE, 
        ##                 edge.col = branch.col.pc, tip.col = pheno.colors, lwd = 3, adj = 0.5, lab4ut = "axial")
        xx = get("last_plot.phylo", envir = ape::.PlotPhyloEnv)$xx[1:n][match(XX.comid, 
            tree$tip.label)]
        yy = get("last_plot.phylo", envir = ape::.PlotPhyloEnv)$yy[1:n][match(XX.comid, 
            tree$tip.label)]
        points(1.08 * xx, 1.08 * yy, col = predCOLS, xpd = TRUE, 
            lwd = 0.5, cex = 1.5)
        ape::edgelabels(edgeLabels, frame = "none", bg = "none", 
            col = edgeLabelsCOL)
        par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 
            0, 0, 0), new = TRUE)
        plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
        if (length(pc.lim) > 0) {
            colourPalette[which(is.na(matchedPCs) == TRUE)] = "grey50"
            legend("right", c(add.signif, "Other"), fill = c(colourPalette[pc.lim], 
                "grey50"), bty = "n", cex = 1.8, xpd = TRUE, 
                inset = c(0, 0))
        }
    })
}
