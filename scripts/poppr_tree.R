#!/usr/bin/env Rscript
library(vcfR)
library(poppr)
library(ape)
library(optparse)
 
option_list = list(
  make_option(c("-i", "--vcf"), type="character", default=NULL, 
              help="vcf file", metavar="character"),
  make_option(c("-t", "--tree"), type="character", default=NULL, 
	      help="output tree name", metavar="character"),
  make_option(c("-m","--method"), type="character",default="bionj", 
		help="tree building method ('nj' or 'upgma') [default= %default]",metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
if ( is.null(opt$vcf) || is.null(opt$tree) ) {
	print_help(opt_parser)
	stop("No vcf input file or tree file via -i and -t options")
}

bd.VCF <- read.vcfR(opt$vcf)
gl.bd <- vcfR2genlight(bd.VCF)
ploidy(gl.bd) <- 2

tree <- aboot(gl.bd, tree = opt$method, distance = bitwise.dist, sample = 100, missing="loci",
              showtree = F, cutoff = 50, quiet = F)
write.tree(tree,opt$tree)
