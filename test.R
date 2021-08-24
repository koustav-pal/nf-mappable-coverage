# test script for testing mappable coverage of junctions
setwd("/Users/koustavpal/Repositories/nf-mappable-coverage")
source("funcrions.R")
annotation_file <- "Homo_sapiens.GRCh38.104.chromosome.19.gff3"
# format <- "gff3"
require("rtracklayer")

gene_ranges <- produce_splice_junctions(annotation_file = annotation_file, 
	format = "gff3")

annotation_file <- "Homo_sapiens.GRCh38.104.chr19.gtf"
gene_ranges <- produce_splice_junctions(annotation_file = annotation_file, 
	format = "gtf")
