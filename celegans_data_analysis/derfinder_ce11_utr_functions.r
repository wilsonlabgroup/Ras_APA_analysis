# load libraries
# suppressMessages(library(edgeR))
# library(derfinder)
# library(derfinderPlot)
# library(regionReport)
# library(GenomicRanges)
# library(GenomicFeatures)
# library(Biostrings)
# library(bumphunter)
# library(gplots)
# library(RColorBrewer)

used_libraries <- c("ggpubr", "edgeR", "RUVSeq", "derfinder", "derfinderPlot", "regionReport", "GenomicRanges", "GenomicFeatures", "Biostrings", "bumphunter", "gplots","ggplot2","ggrepel", "gridExtra", "RColorBrewer", "tidyverse", "UpSetR","reshape2", "pcaMethods", "pheatmap","DESeq2", "openxlsx")

lapply(used_libraries, require, character.only = T, quietly = T)
#ggcolors <- c("N2"="#F8766D", "cfim1"="#00BFC4", "daf2-cfim1"="#619CFF" )
#library(TxDb.)
suppressMessages(library(BSgenome.Celegans.UCSC.ce11))

# set factor levels of conditions
con_order <- c("N2", "daf2", "daf16", "daf2-daf16", "cfim1", "daf2-cfim1")

cols <- colorRampPalette(brewer.pal(8, "RdYlBu"))(8)
condition_cols <- cols[c(1:4, 8,6)]
names(condition_cols) <- con_order

#comparison_labels <- c("notSig", "others", "longest_N2", "shorter_mut", "longest_N2:shorter_mut", "longest_mut:shorter_N2")

comparison_labels <- rev(c("distal_N2:proximal_mut","distal_mut:proximal_N2","proximal_mut","distal_N2", "others", "notSig"))


comparison_cols <- c("grey60", cols[4], cols[2:1], cols[6], cols[8])
names(comparison_cols) <- comparison_labels


hmcol<- colorRampPalette(brewer.pal(9, 'GnBu'))(100)
mono_hue_cols <- colorRampPalette(brewer.pal(9, 'Blues'))(100)

get_count_table <- function(counts_file){
  # load count table from featureCount (one run), also process colnames to get a sample condition table
  count_table <- read.table(counts_file, sep="\t", header=T, as.is = T, stringsAsFactors = F, check.names = FALSE)
  row.names(count_table) <- count_table[,1]
  gene_info <- count_table[,1:6]
  count_table <- count_table[,7:ncol(count_table)]
  # obtain sample information
  sampleFiles <- names(count_table)
  samplenames <- unname(sapply(basename(sampleFiles), function(x) {
    temp <- unlist(strsplit(x, "\\."))
    return(temp[grepl("^WL", temp)])}))
  samplenames <- unname(sapply(samplenames, function(x) paste(unlist(strsplit(basename(x), "_"))[c(1,7)], collapse="_")))
  #samplenames <- unname(sapply(samplenames, fixname))
  
  # rename columns of the count table
  names(count_table) <- samplenames
  
  # build a sample condition table from sample names 
  sampleCondition <- transform_df_samplename(samplenames, field = 2, list = FALSE)
  row.names(sampleCondition) <- samplenames
  #sampleCondition$rep <- sapply(row.names(sampleCondition), function(x) substring(x, nchar(x)))
  
  return(list(count_table, sampleCondition))
}

condition_col_list <- list(condition = condition_cols)

plotHM <- function(counts, method="pearson", cols=colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100), anno_df, show_col_names=T, show_row_names=T,new_col_names=colnames(counts), clust_method="complete"){
  mat<- as.matrix(cor(counts, method=method))
  
  # hc <- hclust(distsRL)
  #anno_df <- anno[, c("age","sex")]
  
  rownames(anno_df) <- rownames(mat) <- colnames(mat) <- new_col_names
  pheatmap(mat, 
           clustering_distance_cols = "correlation",
           clustering_distance_rows = "correlation",
           clustering_method = clust_method,
           annotation_col = anno_df,
           #annotation_colors = anno_cols,
           show_rownames = show_row_names,
           show_colnames = show_col_names,
           color = cols,
           annotation_colors = condition_col_list,
           border_color = NA)
}


# plotHM <- function(counts, method="pearson"){
#   distsRL <- as.dist(1-cor(counts, method=method))
#   #distsRL <- dist(t(counts))
#   mat<- as.matrix(cor(counts, method=method))
#   disfun <- function(x) {as.dist(1-x)}
#   #mat <- as.matrix(distsRL)
#   rownames(mat) <- colnames(mat) <- colnames(counts)
#   hc <- hclust(distsRL, method="average")
#   
#   #row_order <- hc$labels[hc$order]
#   # condition_df <- plename(colnames(mat), field = 2, list = F)
#   # sex_colors <- sapply(condition_df$sex, function(x) color_list[[as.character(x)]])
#   # age_colors <- sapply(condition_df$age, function(x) color_list[[as.character(x)]])
#   # batch_colors <- sapply(condition_df$batch, function(x) color_list[[as.character(x)]])
#   # colors <- cbind(sex_colors, age_colors, batch_colors)
#   
#   heatmap.2(mat, Rowv=as.dendrogram(hc),
#             Colv = as.dendrogram(hc),
#             distfun = disfun,
#             symm=TRUE, trace='none',
#             col = hmcol, margin=c(10, 13))#,
#   #           ColSideColors = colors,
#   #           ColSideColorsSize = 2)
#   # used_colors <- c(unique(as.character(condition_df$sex)),unique(as.character(condition_df$batch)),unique(as.character(condition_df$age))) 
#   # 
#   # used_color_list <- color_list[used_colors]
#   # 
#   # legend("topright",legend=names(used_color_list),
#   #        fill=unlist(used_color_list), border=TRUE, bty="n", y.intersp = 0.7, cex=0.7)
# }

get_samplesfiles <- function(sample_dir, strand="str1"){
  samplefiles <- list.files(path=sample_dir,
                            pattern = glob2rx(paste0("*", strand, "*bw")),
                            full.names = T)
 # files <- rawFiles(datadir = NULL,
  #                  sampledirs = samplefiles,
  #                  fileterm = NULL)
  samplenames <- sapply(basename(samplefiles), function(x) paste(unlist(strsplit(x, "_"))[c(1,7)], collapse = "_"))
  names(samplefiles) <- samplenames
  return(samplefiles)
}



# transform_df_samplename <- function(samplenames_all, field=7, list=TRUE){
#   transformed_list <- lapply(samplenames_all, function(samplenames){
#     print(samplenames)
#     name_elements <- lapply(as.character(samplenames), function(x) unlist(strsplit(x, "_")))
#     samplenames <- sapply(name_elements, "[[", field)
#     
#     samplefields <- strsplit(samplenames, "-")
#     
#     #condition <- sapply(samplefields, function(x) paste(x[-length(x)], collapse="-"))
#     #rep <- sapply(samplefields, function(x) x[length(x)])
#     ID <- sapply(name_elements, "[[", 1)
#     condition <- sapply(name_elements, "[[", 2)
#     rep <- sapply(name_elements, "[[", 3)
#     if(list){
#       return(list(condition, rep, ID))
#     } else {
#       return(data.frame(condition=condition, rep=rep, id=ID))
#     }
#   })
#   outdf <- do.call("rbind", transformed_list)
# #  outdf$batch <- batches[as.character(outdf$id)]
#   return(outdf)
# }

transform_df_samplename <- function(samplenames, field=7, list=TRUE){
    print(samplenames)
    name_elements <- lapply(as.character(samplenames), function(x) unlist(strsplit(x, "_")))
    samplenames <- sapply(name_elements, "[[", field)
    
    #condition <- sapply(samplefields, function(x) paste(x[-length(x)], collapse="-"))
    #rep <- sapply(samplefields, function(x) x[length(x)])
    ID <- sapply(name_elements, "[[", 1)
    condition <- sapply(name_elements, "[[", 2)
    rep <- sapply(name_elements, "[[", 3)
    if(list){
      return(list(condition, rep, ID))
    } else {
      return(data.frame(condition=condition, rep=rep, id=ID))
    }
}




get_txid <- function(anno, type="transcript_id"){
  anno_list <- unlist(strsplit(anno, "; "))
  output <- gsub(paste0(type, " "), "", anno_list[grepl(type, anno_list)])
  if(identical(output, character(0))){
    return(NA)
  } else{
    print(output)
    return(output)
  }
}

grange_add <- function(grange, add){
  meta1 <- as.data.frame(elementMetadata(grange))
  meta2 <- as.data.frame(add)
  meta <- cbind(meta1, meta2)
  elementMetadata(grange) <- meta
  return(grange)
}

get_pca_var <- function(res){
  all_var <- R2cum(res)
  var1 <- round(all_var[1], 3)
  var2 <- round(all_var[2] - all_var[1], 3)
  var3 <- round(all_var[3] - all_var[2], 3)
  return(c(var1, var2))
}

get_pca_plotdata <- function(cpm_table, top=500, npc=15){
  sds <- apply(cpm_table, 1, sd)
  
   variable_genes <- names(sort(-abs(sds)))[1:top]
 # variable_genes <- names(sds)
  variable_counts <- cpm_table[row.names(cpm_table) %in% variable_genes, ]
  
  center_data <- prep(t(variable_counts), scale="none", center=TRUE)
  resPCA <- pca(center_data, method="svd", center=FALSE, nPcs=npc)
  
  #plot_data <- scores(resPCA)[,1:2]
  plot_data <- scores(resPCA)
  plot_data <- cbind(plot_data, transform_df_samplename(row.names(plot_data), field = 2, list = FALSE))
  #plot_data$age <- factor(plot_data$age, levels=c("d12","d22","d27","d32","d37"))
  return(list(plot_data, resPCA))
}


get_pca_plot <- function(cpm_table, npc=5, top=500){
  PCA <- get_pca_plotdata(cpm_table,npc=npc,top=top)
  var <- get_pca_var(PCA[[2]])
  plotdata <- PCA[[1]]
  plotdata$condition <- factor(plotdata$condition, levels=con_order)
  p <- ggplot(plotdata, aes(x=PC1, y=PC2)) +
    geom_point(aes(color=condition), size=6, alpha=0.7) +
    theme_bw() +
    geom_text_repel(aes(label=id)) +
    # scale_color_manual(values = c("red", "blue"))+
    xlab(paste("PC1 (", var[1], ")")) +
    ylab(paste("PC2 (", var[2], ")")) +
    scale_color_manual(values=condition_cols)
  return(list(p, PCA[[2]]))
}

reformat_anno <- function(map_to_ends, txends, type="ends", nonUTR_regions=nonUTR_regions){
  anno_ends <- nonUTR_regions[from(map_to_ends)]
  ends_meta <- txends[to(map_to_ends)]
  
  anno_ends <- cbind(as.data.frame(anno_ends, row.names=NULL, stringsAsFactors=F), as.data.frame(ends_meta, row.names = NULL, stringsAsFactors=F)[, c("tx_id","tx_name")]) 
  #anno_end_txnames <- as.character(sapply(anno_ends$tx_name, function(x) unlist(strsplit(x, ";"))[1]))
  anno_ends$gene_id <- txgenemap[as.character(anno_ends$tx_id)]
  anno_ends$gene_name <- genename[anno_ends$gene_id]
  anno_ends$type <- rep(type, nrow(anno_ends))
  anno_ends$distance <- mcols(map_to_ends)$distance
  # anno_ends2 <- anno_ends %>%
  #   group_by(seqnames, start, end, strand) %>%
  #   summarise(tx_name=paste(sort(tx_name), collapse=";"), gene_id = gene_id[1]) %>%
  #   mutate(gene_name = genename[gene_id], type=type) %>%
  #   as.data.frame()
  anno_ends <- makeGRangesFromDataFrame(anno_ends, keep.extra.columns = T)
  return(anno_ends)
}

grange_add <- function(grange, add){
  meta1 <- as.data.frame(elementMetadata(grange))
  meta2 <- as.data.frame(add)
  meta <- cbind(meta1, meta2)
  elementMetadata(grange) <- meta
  return(grange)
}



grange_extend <- function(x, upstream=0, downstream=0){
  if (any(strand(x) == "*"))
    warning("'*' ranges were treated as '+'")
  on_plus <- strand(x) == "+" | strand(x) == "*"
  new_start <- start(x) - ifelse(on_plus, upstream, downstream)
  new_end <- end(x) + ifelse(on_plus, downstream, upstream)
  ranges(x) <- IRanges(new_start, new_end)
  trim(x)
}

extend_reduce <- function(anno_additional, left=40, right=20){
  #function to reduce the extanded regions and match the gene id back
  anno_extended <- grange_extend(anno_additional, left, right)
  anno_ends_extended_reduced <- reduce(anno_extended)
  
  temp <- findOverlaps(anno_ends_extended_reduced, anno_extended, select="first")
  mcols(anno_ends_extended_reduced)$gene_id <- mcols(anno_extended[temp])$gene_id
  return(anno_ends_extended_reduced)
}

grange_to_gtf <- function(refsequtr, type="RefSeq", outname){
  outgtf <- data.frame(chr=seqnames(refsequtr), 
                       type=type,
                       feature="exon",
                       start=start(refsequtr),
                       end=end(refsequtr),
                       noname1=".",
                       strand=strand(refsequtr),
                       noname2=".",
                       anno=paste('gene_id "', refsequtr$gene_id, '";', sep=""))
  
  write.table(unique(outgtf), outname, quote=F, row.names=F, col.names = F, sep="\t")
}

reformat_anno <- function(map_to_ends, txends, type="ends", nonUTR_regions=nonUTR_regions){
  anno_ends <- nonUTR_regions[from(map_to_ends)]
  ends_meta <- txends[to(map_to_ends)]
  
  anno_ends <- cbind(as.data.frame(anno_ends, row.names=NULL, stringsAsFactors=F), as.data.frame(ends_meta, row.names = NULL, stringsAsFactors=F)[, c("tx_id","tx_name")]) 
  #anno_end_txnames <- as.character(sapply(anno_ends$tx_name, function(x) unlist(strsplit(x, ";"))[1]))
  anno_ends$gene_id <- txgenemap[as.character(anno_ends$tx_id)]
  anno_ends$gene_name <- genename[anno_ends$gene_id]
  anno_ends$type <- rep(type, nrow(anno_ends))
  anno_ends$distance <- mcols(map_to_ends)$distance
  # anno_ends2 <- anno_ends %>%
  #   group_by(seqnames, start, end, strand) %>%
  #   summarise(tx_name=paste(sort(tx_name), collapse=";"), gene_id = gene_id[1]) %>%
  #   mutate(gene_name = genename[gene_id], type=type) %>%
  #   as.data.frame()
  anno_ends <- makeGRangesFromDataFrame(anno_ends, keep.extra.columns = T)
  return(anno_ends)
}

grange_to_gtf <- function(refsequtr, type="RefSeq", outname){
  outgtf <- data.frame(chr=seqnames(refsequtr), 
                       type=type,
                       feature="exon",
                       start=start(refsequtr),
                       end=end(refsequtr),
                       noname1=".",
                       strand=strand(refsequtr),
                       noname2=".",
                       anno=paste('gene_name "', refsequtr$gene_name, '";', sep=""))
  
  write.table(unique(outgtf), outname, quote=F, row.names=F, col.names = F, sep="\t")
}



# explore sequence composiiton of intronic signals
downstream_seq_sum_plot <- function(tx_regions, title="", start=-20, end=250){
  # tx_regions should have "region_name" column which is the index from all regions
  tx_downstream_seq <- downstream_seq[as.numeric(tx_regions$region_name)]
  letter_freq <- consensusMatrix(tx_downstream_seq)[1:4,]
  letter_freq <- as.data.frame(melt(letter_freq), stringAsFactors=F)
  names(letter_freq) <- c("Nucleotide", "pos","counts")
  letter_freq$pos <- as.numeric(letter_freq$pos) + start
  
  # p <- ggplot(letter_freq) +
  #   geom_bar(aes(x=pos, y=counts, fill=letter), stat="identity", position="stack") +
  #   scale_x_continuous(breaks = seq(-20, 250, by = 20), expand = c(0,0)) +
  #   scale_y_continuous(expand = c(0,0)) +
  #   theme_linedraw() +
  #   xlab("relative position from ends of regions")
  
  p <- ggplot(letter_freq) +
    geom_line(aes(x=pos, y=counts, color=Nucleotide, group=Nucleotide)) +
    scale_x_continuous(breaks = seq(start, end, by = floor(0.1*(end-start))), expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    theme_bw() +
    xlab("relative position from ends of regions") +
    ggtitle(title) +
    theme(legend.position = c(0.9, 0.8),
          plot.title = element_text(hjust=0.5, margin=margin(t=10, b=-20)))
  #theme(legend.position = c(max(letter_freq$pos)-10, max(letter_freq$counts)-10))
  
  return(p)
}

# a function to get number of reads mapped to a region from bigwig signal.
get_signal <- function(region,full_signal, return_df = F){
  chr <- as.character(seqnames(region))
  chr_signal <- full_signal[[chr]]
  signal_df <- as.data.frame(chr_signal[start(region):end(region),])
  signal <- colSums(signal_df)
  if(return_df){
    return(signal_df)
  }else {
    return (signal)
  }
}

reformat_regions <- function(tx_df){
  index <<- index +1
  
  gene <- tx_df$gene_id[1]
  print(gene)
  tx_df <- as.data.frame(tx_df)
  strand = as.character(tx_df[1, "strand"])
  
  # reorder, putting the longest UTR at first row
  tx_df = tx_df[with(tx_df, order(start, end)),]
  if(strand == "+"){
    tx_df = tx_df[rev(1:nrow(tx_df)),]
  }
  row.names(tx_df) <- seq(1, nrow(tx_df))
  count_df <- tx_df[grepl("WL", names(tx_df))]
  
  # sum of all samples for each PA site
  pasums <- rowSums(count_df)
  
  # if the internal polyA site is the longest, or second longest followed by a real polyA site, and also has significant proportion of reads mapped, make it a valid site
  if(tx_df[1,"valid"]=="no" & pasums[1]/sum(pasums) > 0.5){
    tx_df[1, "valid"]="yes_internal"
  }
  if(tx_df[2,"valid"]=="no" & tx_df[1, "valid"] == "yes" & pasums[2]/sum(pasums) > 0.5){
    tx_df[2, "valid"]="yes_internal"
  }
  
  # if(nrow(subset(tx_df, valid=="no"))==nrow(tx_df)){
  #   return("only_internalAs")
  # }
  
  # check if the internal polyA site is the longest and a major site  
  
  # add up all exonic signals and return a count table
  # include exonic internal As; only excluding non-exonic internal As. 
  tx_df_for_gene_counting <- subset(tx_df, !(type != "exon" & valid == "no"))
  count_df_for_gene_counting <- tx_df_for_gene_counting[grepl("WL", names(tx_df_for_gene_counting))]
  count_table_gene[index, ] <<- colSums(count_df_for_gene_counting)
  
  # if the majority of reads come from internal polyA events, discard the gene
  internalA_perc <- sum(pasums[which(tx_df$valid == "no")])/sum(pasums)
  if(internalA_perc > 0.75){
    return("Majority polyAs")
  }
  
  # After excluding internal polyA events that are not the longest signal, if there is less than 2 region left, discard
  tx_df_fil <- subset(tx_df, valid != "no")
  
  if(nrow(tx_df_fil) < 2){
    return("single PA after internalAs")
  }
  
  
  # merge regions that are within 25bp of each other 
  gr <- makeGRangesFromDataFrame(tx_df_fil, keep.extra.columns = T)
  gr_reduce <- GenomicRanges::reduce(gr, min.gapwidth=25, with.revmap=T)
  
  revmap <- gr_reduce$revmap
  count_list <- extractList(mcols(gr)[grepl("WL", names(mcols(gr)))], revmap)
  count_list <- lapply(count_list, function(x) colSums(as.data.frame(x)))
  count_df_temp <- do.call("rbind", count_list)
  
  
  # after filtering, only leaving with regions that contain more than 5% of all reads and have more than 10 reads in at least 3 samples
  
  #count_df2 <- tx_df_fil[grepl("WL", names(tx_df_fil))]
  pasums2 <-  rowSums(count_df_temp)
  regions_perc <- pasums2/sum(pasums2) 
  include1 <- which(regions_perc > 0.05)
  
  include2 <- which(rowSums(count_df_temp > 10) > 3)
  include <- intersect(include1, include2)
  if(length(include) == 0){
    return("Too few reads")
  }
  # calculate the percentage of reads that ended up being included 
  include_perc <- sum(pasums2[include])/sum(pasums2)
  
  #count_df_final <- count_df_temp[include,]
  tx_df_final <- as.data.frame(gr_reduce[include,])
  
  if(nrow(tx_df_final) < 2){
    return("single PA after merging and filtering for low counts")
  }
  
  tx_df_final$gene_id <- rep(gene, nrow(tx_df_final))
  tx_df_final <- cbind(tx_df_final, count_df_temp[include,])
  
  
  
  # sort the data frame again
  tx_df_final = tx_df_final[with(tx_df_final, order(start, end)),]
  
  # 20200903 update: order dataframe differently
  if(strand == "-"){
    tx_df_final = tx_df_final[rev(1:nrow(tx_df_final)),]
  }
  row.names(tx_df_final) <- 1:nrow(tx_df_final)
  
  # tx_df_gr <- makeGRangesFromDataFrame(tx_df_final)
  count_df_final <- tx_df_final[,grepl("WL", colnames(tx_df_final))]
  
  
  prop_df <- t(apply(count_df_final, 1, function(x) x/colSums(count_df_final)))
  # get the proportion of the longest and shortest UTRs
  
  # 20200903 update: get the proportion of the shortest
  #prop_df <- prop_df %>% map_df(rev)
  utrindex_short <- as.data.frame(prop_df)[1,]
  utrindex_long <- as.data.frame(prop_df)[nrow(prop_df),]
  row.names(count_df_final) <- paste0(unique(tx_df_final$gene), ":E00", row.names(tx_df_final))
  outgr <- makeGRangesFromDataFrame(tx_df_final, keep.extra.columns = T)
  names(outgr) <- row.names(count_df_final)
  # ifvalidPA <- nrow(subset(tx_df_final, type=="validPA"))
  return(list(utrindex_long, count_df_final, outgr, tx_df_fil, include_perc, internalA_perc, utrindex_short))
}


get_counts_multi <- function(tx_df){
  index <<- index +1
  
  gene <- tx_df$gene_id[1]
  print(gene)
  tx_df <- as.data.frame(tx_df)
  strand = as.character(tx_df[1, "strand"])
  
  # reorder, putting the longest UTR at first row
  tx_df = tx_df[with(tx_df, order(start, end)),]
  if(strand == "+"){
    tx_df = tx_df[rev(1:nrow(tx_df)),]
  }
  row.names(tx_df) <- seq(1, nrow(tx_df))
  count_df <- tx_df[grepl("WL", names(tx_df))]
  
  # sum of all samples for each PA site
  pasums <- rowSums(count_df)
  
  # if the internal polyA site is the longest, or second longest followed by a real polyA site, and also has significant proportion of reads mapped, make it a valid site
  if(tx_df[1,"valid"]=="no" & pasums[1]/sum(pasums) > 0.5){
    tx_df[1, "valid"]="yes_internal"
  }
  if(tx_df[2,"valid"]=="no" & tx_df[1, "valid"] == "yes" & pasums[2]/sum(pasums) > 0.5){
    tx_df[2, "valid"]="yes_internal"
  }
  
  # if(nrow(subset(tx_df, valid=="no"))==nrow(tx_df)){
  #   return("only_internalAs")
  # }
  
  # check if the internal polyA site is the longest and a major site  
  
  # add up all exonic signals and return a count table
  # include exonic internal As; only excluding non-exonic internal As. 
  tx_df_for_gene_counting <- subset(tx_df, !(type != "exon" & valid == "no"))
  count_df_for_gene_counting <- tx_df_for_gene_counting[grepl("WL", names(tx_df_for_gene_counting))]
  count_table_gene[index, ] <<- colSums(count_df_for_gene_counting)
}


wilcox_test <- function(data1, data2){
   alt <- ifelse(median(data1, na.rm=T) > median(data2, na.rm=T), "greater", "less")
   test <- wilcox.test(data1, data2, alternative = alt)
   return(test$p.value)
}
temp_function <- function(used_data, con){
  if(con!="notSig"){
    data1 <- subset(used_data, APAresult_simple == "notSig")$meanlogCPM
    data2 <- subset(used_data, APAresult_simple == con)$meanlogCPM
    return(wilcox_test(data1, data2))
  }
}


my_plot <- function (object, region, annotation = NULL, assay = NULL, fragment.path = NULL, 
                     group.by = NULL, window = 100, downsample = 0.1, height.tracks = 2, 
                     extend.upstream = 0, extend.downstream = 0, cells = NULL, 
                     idents = NULL, sep = c("-", "-")) 
{
  # cells <- cells %||% colnames(x = object)
  if (!is.null(x = idents)) {
    ident.cells <- WhichCells(object = object, idents = idents)
    cells <- intersect(x = cells, y = ident.cells)
  }
  if (!is(object = region, class2 = "GRanges")) {
    region <- StringToGRanges(regions = region, sep = sep)
  }
  region <- suppressWarnings(expr = Extend(x = region, upstream = extend.upstream, 
                                           downstream = extend.downstream))
  reads <- GetReadsInRegion(object = object, assay = assay, 
                            region = region, cells = cells, group.by = group.by, 
                            fragment.path = fragment.path, verbose = FALSE)
  cells.per.group <- CellsPerGroup(object = object, group.by = group.by)
  reads.per.group <- AverageCounts(object = object, group.by = group.by, 
                                   verbose = FALSE)
  coverages <- suppressWarnings(CalculateCoverages(reads = reads, 
                                                   cells.per.group = cells.per.group, reads.per.group = reads.per.group, 
                                                   window = window, verbose = FALSE))
  if (downsample > 1) {
    warning("Requested downsampling <0%, retaining all positions")
    downsample <- 1
  }
  chromosome <- as.character(x = seqnames(x = region))
  start.pos <- start(x = region)
  end.pos <- end(x = region)
  stepsize <- 1/downsample
  total_range <- end.pos - start.pos
  steps <- ceiling(x = (total_range/stepsize))
  retain_positions <- seq(from = start.pos, to = end.pos, by = stepsize)
  downsampled_coverage <- coverages[coverages$position %in% 
                                      retain_positions, ]
  ymax <- signif(x = max(downsampled_coverage$coverage, na.rm = TRUE), 
                 digits = 2)
  ymin <- 0
  downsampled_coverage <- downsampled_coverage[!is.na(x = downsampled_coverage$coverage), 
                                               ]
  p <- ggplot(data = downsampled_coverage, mapping = aes(x = position, 
                                                         y = coverage, fill = group)) + geom_area(stat = "identity") + 
    geom_hline(yintercept = 0, size = 0.1) + facet_wrap(facets = ~group, 
                                                        strip.position = "right", ncol = 1) + xlab(label = paste0(chromosome, 
                                                                                                                  " position (bp)")) + ylab(label = paste0("Normalized accessibility \n(range ", 
                                                                                                                                                           as.character(x = ymin), " - ", as.character(x = ymax), 
                                                                                                                                                           ")")) + ylim(c(ymin, ymax)) + theme_classic() + theme(axis.text.y = element_blank(), 
                                                                                                                                                                                                                 legend.position = "none", strip.text.y = element_text(angle = 0))
  if (!is.null(x = annotation)) {
    gr <- GRanges(seqnames = gsub(pattern = "chr", replacement = "", 
                                  x = chromosome), IRanges(start = start.pos, end = end.pos))
    filters <- AnnotationFilterList(GRangesFilter(value = gr), 
                                    GeneBiotypeFilter(value = "protein_coding"))
    if (suppressMessages(expr = nrow(x = ensembldb::select( annotation, 
                                                            filters)) > 0)) {
      genes <- suppressMessages(expr = autoplot(object = annotation, 
                                                filters, names.expr = "gene_name"))
      gene.plot <- genes@ggplot + xlim(start.pos, end.pos) + 
        xlab(label = paste0(chromosome, " position (bp)")) + 
        theme_classic()
      p <- p + theme(axis.title.x = element_blank(), axis.text.x = element_blank(), 
                     axis.line.x.bottom = element_blank(), axis.ticks.x.bottom = element_blank())
      p <- suppressWarnings(plot_grid(p, gene.plot, ncol = 1, 
                                      axis = "btlr", rel_heights = c(height.tracks, 
                                                                     1), align = "v"))
    }
  }
  return(p)
}


get_sem <- function(x){
  return(round(sd(x)/sqrt(length(x)),2))
}

plot_utrindex_bar <- function(utrindex_used = proportion_df, plot_gene, plot_cons = c("N2", "cfim1", "daf2-cfim1")){
  #plot_geneid <- names(genename[which(genename == plot_gene)])
  
  plot_samples <- as.character(subset(sampleCondition, condition %in% plot_cons)$sample)
  
  utrindex_counts <- utrindex_used %>%
    filter(genename == plot_gene) %>%
    gather(sample, utrindex, -genename) %>%
    filter(sample %in% plot_samples) %>%
    left_join(sampleCondition, by = "sample") %>% 
    mutate(utrindex = 100*utrindex)
  
  
  
  utrindex_counts_mean <- utrindex_counts %>%
    group_by(genename, condition) %>%
    summarise(mean = mean(utrindex), sem = get_sem(utrindex))
  
  print("to get plot")
  
  p_utrindex_bar <- ggplot() +
    geom_bar(data = utrindex_counts_mean, aes(x = genename, y = mean, fill = condition), stat = "identity", position = position_dodge(.8), width = .7) +
    geom_errorbar(data = utrindex_counts_mean, aes(x = genename, ymin = mean-sem, ymax = mean+sem, group = condition, color = condition), position = position_dodge(.8), width = .3) +
    scale_fill_manual(values = condition_col_list$condition) +
    scale_color_manual(values = condition_col_list$condition) +
    scale_y_continuous(expand = c(0,0))+
    #geom_point(data = gene_counts, aes(x = PAsite, y = counts, color = condition)) +
    theme_classic() +
    theme(axis.text.x = element_text(color = "black", angle = 30, hjust = 1, vjust = 1, size = 12),
          axis.text.y = element_text(color = "black", size = 12),
          axis.title.y = element_text(size = 14),
          axis.title.x = element_blank()) +
    #  xlab(paste(plot_gene, " isoforms")) +
    ylab("PPAU \n(percent of reads \nmapped to proximal PA site) (%)")
  
  print("got plot")
  
  ggsave(paste0(plot_gene, "_PPAU_bar.pdf"), p_utrindex_bar, width = 5, height = 4)
  return(utrindex_counts)
}


ppau_vs_logcpm <- function(used_condition = "N2", APA_comparison = "cfim1_vs_N2", logCPM, ppau_df){
  used_samples <- subset(sampleCondition, condition == used_condition)$sample
  
  used_expr <- data.frame(geneID = row.names(logCPM),
                          expr = rowMeans(logCPM[, used_samples]))
  # get APA summary 
  
  APAresult_sum_used <- subset(APAresult_sum_df, L1 == APA_comparison)
  
  cons <- unlist(strsplit(APA_comparison, "_vs_"))
  
  utrindex_mergedReps <- plot_proportion_df_mergeReps %>% 
    spread(condition, UTRindex) %>% 
    left_join(used_expr, by = c("Var1" = "geneID")) %>% 
    left_join(APAresult_sum_used, by = c("Var1" = "geneid")) %>% 
    mutate(label = ifelse(is.na(label), "notSig", as.character(label))) %>% 
    mutate(label = factor(label, levels = comparison_labels))
  
  utrindex_mergedReps[,"ppau"] <- utrindex_mergedReps[, cons[1]] - utrindex_mergedReps[, cons[2]]
  
  p <- ggplot(utrindex_mergedReps) +
    geom_point(aes(y= ppau, x = expr, color = label), alpha = .7) +
    geom_hline(yintercept = 0) +
    scale_color_manual(values = comparison_cols) +
    xlab(paste0("average expression (log2CPM) of ", used_condition, " samples")) +
    ylab(paste0("pPAU (", cons[1], "-", cons[2], ")")) +
    theme_bw()
  ggsave(paste0("logcpm_vs_ppau_", APA_comparison, "_", Sys.Date(), ".pdf"), width = 8, height = 6)
  
}