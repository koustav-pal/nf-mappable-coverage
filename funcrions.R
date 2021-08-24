.read_annotation <- function(annotation_file, ...){
    arg_list <- list(...)
    if(arg_list$format == "gff3"){
        annotation_ranges <- import.gff3(annotation_file)
    }else if(arg_list$format == "gtf"){
        annotation_ranges <- import.gff2(annotation_file)
    }
    return(annotation_ranges)
}

.filter_annotation <- function(a_ranges, ...){
    arg_list <- list(...)
    genes_ranges <- a_ranges[a_ranges$type == "gene"]
    all_chromosomes <- seqlevels(genes_ranges)
    if(!is.na(arg_list$read_chromosomes)){
        all_chromosomes <- all_chromosomes[all_chromosomes %in% arg_list$read_chromosomes]
    }
    genes_ranges <- genes_ranges[seqnames(genes_ranges) %in% all_chromosomes]
    
    if(arg_list$format == "gff3"){
        transcript_ranges <- a_ranges[as.vector(a_ranges$Parent) %in% unique(genes_ranges$ID)]
        exon_ranges <- a_ranges[a_ranges$type == "exon" & as.vector(a_ranges$Parent) %in% unique(transcript_ranges$ID)]
    }else if(arg_list$format == "gtf"){
        exon_ranges <- a_ranges[a_ranges$type == "exon" & a_ranges$gene_id %in% genes_ranges$gene_id]
    }
    # filter ranges for genes greater than length
    if(!is.na(arg_list$min_length)){
        genes_ranges <- genes_ranges[width(genes_ranges) > arg_list$min_length]
    }
    if(arg_list$non_overlapping){
        ol_object <- findOverlaps(genes_ranges, genes_ranges, ignore.strand = T)
        ol_object <- ol_object[queryHits(ol_object) != subjectHits(ol_object)]
        keep_indices <- unique(c(queryHits(ol_object), subjectHits(ol_object)))
        genes_ranges <- genes_ranges[!(seq_along(genes_ranges) %in% keep_indices)]
    }
    if(!is.na(arg_list$min_tss_tes_gap)){
        tss_ranges <- genes_ranges
        end(tss_ranges[strand(tss_ranges) != "-"]) <- start(tss_ranges[strand(tss_ranges) != "-"])
        start(tss_ranges[strand(tss_ranges) == "-"]) <- end(tss_ranges[strand(tss_ranges) == "-"])
        start(tss_ranges) <- start(tss_ranges) - arg_list$min_tss_tes_gap
        end(tss_ranges) <- end(tss_ranges) + arg_list$min_tss_tes_gap
        tes_ranges <- genes_ranges
        start(tss_ranges[strand(tss_ranges) != "-"]) <- end(tss_ranges[strand(tss_ranges) != "-"])
        end(tss_ranges[strand(tss_ranges) == "-"]) <- start(tss_ranges[strand(tss_ranges) == "-"])
        start(tes_ranges) <- start(tes_ranges) - arg_list$min_tss_tes_gap
        end(tes_ranges) <- end(tes_ranges) + arg_list$min_tss_tes_gap
        ol_object <- findOverlaps(tss_ranges, tes_ranges)
        ol_object <- ol_object[queryHits(ol_object) != subjectHits(ol_object)]
        keep_indices <- unique(c(queryHits(ol_object), subjectHits(ol_object)))
        genes_ranges <- genes_ranges[!(seq_along(genes_ranges) %in% keep_indices)]
    }
    if(arg_list$format == "gff3"){
        filtered_transcript_ranges <- transcript_ranges[as.vector(transcript_ranges$Parent) %in% genes_ranges$ID]
        filtered_exon_ranges <- exon_ranges[as.vector(exon_ranges$Parent) %in% unique(filtered_transcript_ranges$ID)]
        filtered_exon_ranges$gene_id <- genes_ranges$gene_id[match(
            as.vector(filtered_transcript_ranges$Parent[match(as.vector(filtered_exon_ranges$Parent), 
            filtered_transcript_ranges$ID)]), genes_ranges$ID)]
    }else if (arg_list$format == "gtf"){
        filtered_exon_ranges <- exon_ranges[exon_ranges$gene_id %in% genes_ranges$gene_id]
    }
    filtered_ranges <- list(genes = genes_ranges, exons = filtered_exon_ranges)
    return(filtered_ranges)
}

.create_exon_intron_pair <- function(exon_ranges, gene_ranges) {
    cur_strand <- as.vector(strand(gene_ranges))
    merged_exons <- reduce(exon_ranges, ignore.strand = FALSE)
    names(merged_exons) <- paste(unique(gene_ranges$gene_id), seq_along(merged_exons), sep = "_")
    merged_exons$type <- "exon"
    start(gene_ranges) <- min(start(merged_exons))
    end(gene_ranges) <- max(end(merged_exons))
    if(length(merged_exons) == 1){
        merged_exons$intron_pos <- 0
        merged_exons$gene_id <- unique(gene_ranges$gene_id)
        return(merged_exons)
    }
    intron_ranges <- setdiff(gene_ranges, exon_ranges)
    names(intron_ranges) <- paste(unique(gene_ranges$gene_id), 
        seq_along(merged_exons)[-length(merged_exons)], 
        seq_along(merged_exons)[-1], sep = "_")
    intron_ranges$type <- "intron"
    all_ranges <- c(merged_exons, intron_ranges)
    all_ranges <- all_ranges[order(start(all_ranges))]
    junction_pos <- which(all_ranges$type == "intron")
    shifter <- ifelse(cur_strand == "-", 1, -1)
    upstream_exon <- junction_pos + shifter
    index <- seq_along(junction_pos)
    all_ranges$intron_pos <- NA
    all_ranges$intron_pos[junction_pos] <- index
    all_ranges$intron_pos[upstream_exon] <- index
    all_ranges$intron_pos[all_ranges$type == "exon"] <- NA
    all_ranges$gene_id <- unique(gene_ranges$gene_id)
    return(all_ranges)
}

.create_splice_junctions <- function(exon_ranges, gene_ranges) {
    cur_strand <- as.vector(strand(gene_ranges))
    merged_exons <- reduce(exon_ranges, ignore.strand = FALSE)
    merged_exons$type <- "exon"
    start(gene_ranges) <- min(start(merged_exons))
    end(gene_ranges) <- max(end(merged_exons))
    if(length(merged_exons) == 1){
        merged_exons$intron_pos <- 0
        merged_exons$gene_id <- unique(gene_ranges$gene_id)
        return(merged_exons)
    }
    intron_ranges <- setdiff(gene_ranges, exon_ranges)
    intron_ranges$type <- "intron"
    all_ranges <- c(merged_exons, intron_ranges)
    all_ranges <- all_ranges[order(start(all_ranges))]
    junction_pos <- which(all_ranges$type == "intron")
    if (cur_strand == "-") {
        junction_pos <- rev(junction_pos)
    }
    splice_junction_ranges_list <- lapply(seq_along(junction_pos),
        function(x_of_x) {
        x <- junction_pos[x_of_x]
        current_intron <- all_ranges[x]
        exon_up <- all_ranges[x - 1]
        start(exon_up) <- end(exon_up) - ceiling(width(exon_up) / 2) + 1
        exon_down <- all_ranges[x + 1]
        end(exon_down) <- (end(exon_down) - ceiling(width(exon_down) / 2))
        temp_ranges <- c(exon_up, current_intron, exon_down)
        temp_ranges$intron_pos <- x_of_x
        return(temp_ranges)
    })
    splice_junction_ranges <- do.call(c, unlist(splice_junction_ranges_list,
        use.names = FALSE))
    splice_junction_ranges <- splice_junction_ranges[
        order(start(splice_junction_ranges))]
    splice_junction_ranges$gene_id <- unique(gene_ranges$gene_id)
    return(splice_junction_ranges)
}

produce_splice_junctions <- function(annotation_file, read_chromosomes = NA, 
    format = c("gff3", "gtf"), min_length = 1000, non_overlapping = T,
    min_tss_tes_gap = NA){
    format <- match.arg(format)
    # message(annotation_file, " ", format)
    annotation_ranges <- .read_annotation(annotation_file, format = format)
    filtered_ranges_list <- .filter_annotation(annotation_ranges, 
        format = format, 
        read_chromosomes = read_chromosomes, 
        min_length = min_length,
        min_tss_tes_gap = min_tss_tes_gap,
        non_overlapping = non_overlapping)
    exon_ranges <- filtered_ranges_list[["exons"]]
    gene_ranges <- filtered_ranges_list[["genes"]]
    exon_ranges_split <- split(exon_ranges, exon_ranges$gene_id)
    exon_intron_pair_ranges_list <- lapply(exon_ranges_split,
        function(an_exon_range) {
        gene_id <- unique(an_exon_range$gene_id)
        a_gene_ranges <- gene_ranges[gene_ranges$gene_id == gene_id]
        temp_range <- .create_exon_intron_pair(an_exon_range, a_gene_ranges)
    })
    exon_intron_pair_ranges <- do.call(c, unlist(exon_intron_pair_ranges_list, 
        use.names = FALSE))
    return(exon_intron_pair_ranges)
}

.extend_regions_beyond_ranges <- function(a_ranges, ...){
    arg_list <- list(...)
    preceded_by <- precede(a_ranges, ignore.strand = T)
    preceded_by_seq <- seq_along(preceded_by)
    which_not_na <- which(!is.na(preceded_by)) 
    a_ranges$upstream_distance <- NA
    a_ranges$downstream_distance <- NA
    a_ranges$upstream_distance[preceded_by[which_not_na]] <- floor((start(a_ranges[preceded_by[which_not_na]]) - end(a_ranges[preceded_by_seq[which_not_na]]))/2)
    a_ranges$downstream_distance[preceded_by_seq[which_not_na]] <- floor((start(a_ranges[preceded_by[which_not_na]]) - end(a_ranges[preceded_by_seq[which_not_na]]))/2)
    a_ranges <- a_ranges[!is.na(a_ranges$upstream_distance) & !is.na(a_ranges$downstream_distance)]
    if(!is.na(min_region_width)){
        a_ranges <- a_ranges[a_ranges$upstream_distance > arg_list$min_region_width & 
            a_ranges$downstream_distance > arg_list$min_region_width]
    }
    if(!is.na(arg_list$max_region_width)){
        a_ranges$upstream_distance[a_ranges$upstream_distance > arg_list$max_region_width] <- arg_list$max_region_width
        a_ranges$downstream_distance[a_ranges$downstream_distance > arg_list$max_region_width] <- arg_list$max_region_width
    }
    upstream_region_start <- start(a_ranges) - a_ranges$upstream_distance + 1
    upstream_region_end <- start(a_ranges) - 1
    downstream_region_start <- end(a_ranges) + 1
    downstream_region_end <- end(a_ranges) + a_ranges$downstream_distance
    upstream_ranges <- GRanges(seqnames = as.vector(seqnames(a_ranges)), IRanges(upstream_region_start, upstream_region_end))
    upstream_ranges$gene_id <- paste(a_ranges$gene_id, "1", sep = "_")
    downstream_ranges <- GRanges(seqnames = as.vector(seqnames(a_ranges)), IRanges(downstream_region_start, downstream_region_end))
    downstream_ranges$gene_id <- paste(a_ranges$gene_id, "3", sep = "_")
    a_ranges$gene_id <- paste(a_ranges$gene_id, "2", sep = "_")
    all_ranges <- c(upstream_ranges, a_ranges, downstream_ranges)
    all_ranges <- all_ranges[order(all_ranges$gene_id)]
    names(all_ranges) <- all_ranges$gene_id
    return(all_ranges)
}

produce_regions_sorrounding_genes <- function(annotation_file, read_chromosomes = NA, 
    format = c("gff3", "gtf"), min_length = 1000, non_overlapping = T,
    min_tss_tes_gap = NA, min_region_width = 1000, max_region_width = 5000){
    format <- match.arg(format)
    annotation_ranges <- .read_annotation(annotation_file = annotation_file, format = format)
    filtered_ranges_list <- .filter_annotation(annotation_ranges, 
        format = format,
        read_chromosomes = read_chromosomes,
        min_length = min_length, 
        min_tss_tes_gap = min_tss_tes_gap,
        non_overlapping = T)
    gene_ranges <- filtered_ranges_list[["genes"]]
    gene_ids <- gene_ranges$gene_id
    elementMetadata(gene_ranges) <- NULL
    elementMetadata(gene_ranges)[["gene_id"]] <- gene_ids
    all_ranges <- .extend_regions_beyond_ranges(a_ranges = gene_ranges, 
        min_region_width = min_region_width,
        max_region_width = max_region_width)
    return(all_ranges)
}

