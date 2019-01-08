gc(reset = TRUE)

library(mlr)

# ---------------- #
# input arguments
# ---------------- #
args = commandArgs(trailingOnly=TRUE)
if (length(args) != 4) {
  stop("Number of arguments must be 4.", call.=FALSE)
}

pathToRScripts  <- args[1]
node_desc 		<- args[2]
subgraph_desc 	<- args[3]
output_file 	<- args[4]

# ---------------- #
# load models and scaling parameters
# ---------------- #
load(paste0(pathToRScripts, '/prediction.Rdata'))
# load mu, sigma and the models: binary_modl and multiv_modl

# ---------------- #
# functions
# ---------------- #
distribution <- function(covar, prefix) {
  dist_names 	<- c("Median", "Qu.0.05", "Qu.0.95", "sd")
  vect 			<- c(median(covar),quantile(covar,0.05, na.rm=T),quantile(covar,0.95, na.rm=T),sd(covar))
  names(vect)   <- paste0(prefix, dist_names)
  return(vect) 
}

set2zero <- function(tab){

	tab[which(tab$X.sig_node_number == 1),which(colnames(tab) == "X.sign_node_effect_sd" | colnames(tab) == "X.sign_node_degree_sd") ] <- 0
	tab[which(tab$X.length_neg_sum == 0),which(colnames(tab) == "X.length_neg_mean") ] <- 0
	tab[which(tab$X.length_pos_sum == 0),which(colnames(tab) == "X.length_pos_mean") ] <- 0

	return(tab);
}

standardize <- function(z,mu,sigma) {
	rv <- sweep(z, 2, mu,"-")  #subtracting median expression
	rv <- sweep(rv, 2, sigma, "/")  # dividing by median absolute deviation
	return(rv)
}

read_desc_files <- function(node_f, subgraph_f){

	if(!file.exists(subgraph_f)) message("Error file does not exists.")  	# !!! do something
    subgraphs <- read.table(subgraph_f, header=T) 
    if (nrow(subgraphs) ==0) message("Error file empty.")					# !!! do something
	 	# colnames(subgraphs)
	 	# 1. subgraph_id     2. node_number     3. sig_node_number 4. sig_node_ratio  5.pos_effect_ratio        
	 	# 6. nb_Pheno0    7. nb_Pheno1       

	subgraph_id <- subgraphs$subgraph_id

	if(!file.exists(node_f)) message("Error file does not exists.")  	# !!! do something
    nodes <- read.table(node_f, header=T) 
    if (nrow(nodes) ==0) message("Error file empty.")					# !! do something
		# colnames(nodes)
		# 1. subgraph_id     2. node_id 3. node_sign       4. node_effect     5. node_length     
		# 6. node_degree     7. node_in_Pheno0  8. node_in_Pheno1

    allele_0 	<- nodes$node_in_Pheno0 / ( subgraphs$nb_Pheno0[1] + subgraphs$nb_Pheno1[1] )
    allele_1 	<- nodes$node_in_Pheno1 / ( subgraphs$nb_Pheno0[1] + subgraphs$nb_Pheno1[1] )
    allele_fq 	<- allele_0 + allele_1

	node_list <- lapply(unique(nodes$subgraph_id), function(x){
	    # x is the subgraph id
	    subgraph_level <- c()

	    # 1. distribution of node effect, within signif 
	    covar <- nodes$node_effect[which(nodes$subgraph_id == x)]
	    subgraph_level  <- c(subgraph_level, distribution(covar, "node_effect_"))

	    covar <- nodes$node_effect[which(nodes$subgraph_id == x & nodes$node_sign == 1)]
	    subgraph_level  <- c(subgraph_level, distribution(covar, "sign_node_effect_"))

	    # 2. sum of unitig length (signif, pos effect)
	    k <- 31
	    length_pos_sum  <- sum(nodes$node_length[which(nodes$subgraph_id == x
	            & nodes$node_sign == 1 & nodes$node_effect > 0)])
	    length_neg_sum  <- sum(nodes$node_length[which(nodes$subgraph_id == x
	            & nodes$node_sign == 1 & nodes$node_effect < 0)])
	    length_pos_mean  <- mean(nodes$node_length[which(nodes$subgraph_id == x
	            & nodes$node_sign == 1 & nodes$node_effect > 0)] - k )
	    length_neg_mean  <- mean(nodes$node_length[which(nodes$subgraph_id == x
	            & nodes$node_sign == 1 & nodes$node_effect < 0)] - k )

	    sum_length      <- length_neg_sum + length_pos_sum
	    if(sum_length > 0){
	            diff_length <- abs(length_neg_sum - length_pos_sum) / max(length_pos_sum, length_neg_sum, na.rm = T)
	    }else{
	            diff_length <- 0
	    }
	    vectL           <- c(length_pos_mean, length_neg_mean, length_pos_sum, length_neg_sum, diff_length, sum_length)
	    names(vectL)    <- c("length_pos_mean", "length_neg_mean", "length_pos_sum", "length_neg_sum", "diff_length", "sum_length")
	    subgraph_level  <- c(subgraph_level, vectL)

	    # 3. Node degree within sign. 
	    covar <- nodes$node_degree[which(nodes$subgraph_id == x)]
	    subgraph_level  <- c(subgraph_level, distribution(covar, "node_degree_"))

	    covar <- nodes$node_degree[which(nodes$subgraph_id == x & nodes$node_sign == 1)]
	    subgraph_level  <- c(subgraph_level, distribution(covar, "sign_node_degree_"))

	    # 4. Allele frequency 
	    allele_fq_mean  <- mean( allele_fq[which(nodes$subgraph_id == x )])
	    allele_fq_sd    <- sd( allele_fq[which(nodes$subgraph_id == x )])
	    diff_alleles    <- abs(mean( allele_0[which(nodes$subgraph_id == x )]) -
	            mean( allele_1[which(nodes$subgraph_id == x )]))
	    vectA                   <- c(allele_fq_mean, allele_fq_sd, diff_alleles)
	    names(vectA)    <- c("allele_fq_mean", "allele_fq_sd", "diff_alleles")
	    subgraph_level  <- c(subgraph_level, vectA)

	    return(subgraph_level)
	})

	df <- data.frame(matrix(unlist(node_list), nrow=length(node_list), byrow=T))
    colnames(df) <- names(node_list[[1]])
    rownames(df) <- subgraph_id
    X_tab <- cbind(subgraphs[,2:5], df)
    colnames(X_tab) <- paste0("X.",colnames(X_tab))
	
	return(X_tab)
}

# ---------------- #
# main
# ---------------- #
X_mat <- read_desc_files(node_desc, subgraph_desc)
X_mat <- set2zero(X_mat)
x <- standardize(X_mat, mu,sigma)
bin <- predict(binary_modl,newdata = x)
# mul <- predict(multiv_modl,newdata = x) # no need right now

tab <- data.frame(subrgaph_id=rownames(bin$data), bin$data[,1:2])
write.table(tab, file=output_file, quote=F, row.names=F)

# ---------------- #
# leave
# ---------------- #
rm(list = ls(all = TRUE))


