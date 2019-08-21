# Some useful functions
require(GenomicRanges)
require(GenomicFeatures)
require(GenomicAlignments)
require(Rsamtools)
require(tidyverse)
require(biomaRt)
require(RColorBrewer)

# FUNCTIONS:

#'  Brutal function to convert from byTranscript to byGene GRangesList
#' 
#' @param GRangesList A GRangesList object created by eg. \cdsBy() 
#' with ranges grouped by transcripts ID
#' @param txdb A Txdb object that stores transcripts metadata
#' @param reduce Logical value whether to use \reduce() function on gene models 
byGene <- function(GRangesList, txdb, reduce = T){
  cleaned <- unlist(GRangesList)
  
  # Tx to Gene ID mapping
  txNames <- names(cleaned)
  mapID <- suppressMessages(AnnotationDbi::select(txdb, 
                                                  keys = as.character(txNames), 
                                                  keytype = 'TXNAME', 
                                                  columns = 'GENEID'))
  # Match Tx Id to Gene ID
  toKeep <- mapID$GENEID[!is.na(mapID$GENEID)]
  cleaned <- split(cleaned[!is.na(mapID$GENEID)], toKeep)
  
  if (reduce){
    cleaned <- IRanges::reduce(cleaned)
  }
  
  return(cleaned)
}

#' Define models for genomic regions
#' 
#'  @param txdb A Txdb object that stores transcripts metadata
#'  @param by Character "gene" or "transcript", specifies IDs class by genomic 
#'   should be grouped, currently only "gene" is working
#'  @param mode Character "union" or "unique":
#'   'Union' : takes all regions of given class from all isoforms
#'    eg. class 'exons' will consists from all avaliable
#'    exons that are assigned to a gene
#'  'Uniqe' : only unique regions from a class are considered,
#'    eg. regions that can be assigned as intron and exon
#'    depending on isoform will be discarded, currently avaliable
#'    only when reduced = T
#'  @param use.names Logical whether transcripts/exon names should be included in the model
#'  @param regions Character specifying genomic regions, avaliable: 
#'  "genes" (mandatory), "transcripts", "exons", "introns", "cds", "utr3", "utr5"
#'  @param reduce Logical value whether to use \reduce() function on regions models
#'  
#'  @return List with genomic coordinates of selected regions
buildModels <- function(txdb, mode = "union", by = "gene", use.names = F,
                        regions = c("genes", "transcripts", "exons", "introns", 
                                    "cds", "utr3", "utr5"),
                        reduced = F){
  regions <- as.factor(regions)
  
  if ("genes" %in% regions){
    print("Selecting genes ...")
    genes <- genes(txdb)
    
    if (reduced){
      genes <- GenomicRanges::reduce(genes)
      
    } else {
      genes$ensembl <- splitvec(genes$gene_id, "[.]", 1)
      ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
      biomart <- getBM(attributes=c('ensembl_gene_id',
                                    'gene_biotype'), 
                       filters = "ensembl_gene_id", 
                       values = genes$ensembl,
                       mart = ensembl)
      genes$biotype <- biomart$gene_biotype[match(genes$ensembl,
                                                  biomart$ensembl_gene_id)]
    }
    
  } 
  
  if ("transcripts" %in% regions){
    print("Selecting transcripts ...")
    txs <- transcriptsBy(txdb, by = by)         
    names(txs) <- txs$tx_name
    
    if (reduced){
      txs <- GenomicRanges::reduce(unlist(txs))
    }
  } 
  if ("exons" %in% regions){
    print("Selecting exons ...")
    exons <- exonsBy(txdb, by = by, use.names = use.names) 
    
    if (reduced){
      exons <- GenomicRanges::reduce(unlist(exons))
    }
  }
  if ("introns" %in% regions){
    print("Selecting introns ...")
    introns <- intronsByTranscript(txdb, use.names = T)
    
    if (by == "gene"){
      introns <- byGene(introns, txdb)
    }
    
    if (reduced){
      introns <- GenomicRanges::reduce(unlist(introns))
    }
  }
  if ("cds" %in% regions){
    print("Selecting CDS ...")
    cds <- cdsBy(txdb, by = by, use.names = use.names)
    
    if (reduced){
      cds <- GenomicRanges::reduce(unlist(cds))
    }
  }
  
  if ("utr5" %in% regions){
    print("Selecting 5'UTR ...")
    utr5 <- fiveUTRsByTranscript(txdb, use.names = T)
    
    if (by == "gene"){
      utr5 <- byGene(utr5, txdb)
    }
    
    if (reduced){
      utr5 <- GenomicRanges::reduce(unlist(utr5))
    }
  } 
  if ("utr3" %in% regions){
    print("Selecting 3'UTR ...")
    utr3 <- threeUTRsByTranscript(txdb, use.names = T)
    
    if (by == "gene"){
      utr3 <- byGene(utr3, txdb)
    }
    
    if (reduced){
      utr3 <- GenomicRanges::reduce(unlist(utr3))
    }
  }
  
  allRegions <- mget(ls()[ls() %in% regions])
  
  print("Done!")                   
  return(allRegions)
}

#' Small function to convert \countOverlaps() output to data frame
#' @param overlaps output of \countOverlaps()
#' @param class character, name of a column with computed counts
makeTab <- function(overlaps, class){
  tab <- data.frame(name = names(overlaps),
                    class = as.numeric(overlaps),
                    stringsAsFactors = F)
  colnames(tab) <- c("name", class)
  return(tab)
}

#' Function to count reads by genomic regions build with \buildModels function
#' @param allRegions List with genomic regions models (output of \buildModels function)
#' @param GAlignment BAM file loaded as GAlignment object
#' @param summaryTab Vector of integers with minimum overlap for a read to be assigned
#' to a region, length and order should correspond to genomic regions in \allRegions
#' @param summarize Logical, if \TRUE, reads will be counted by regions only 
#'
#' @return Data frame with number of reads per region (and per gene if summarize = F)
countReadsByRegion <- function(allRegions, GAlignment, minOverlap = NULL,
                               summarize = F){
  regions <- names(allRegions)
  
  # Arguments check  
  if (!(c("genes") %in% regions)){
    print("Genes coordinates not present in annotations!")
  }
  
  if (is.null(minOverlap)){
    minOverlap <- rep(1, length(regions))
  }
  names(minOverlap) <- regions
  
  # Select reads with at least 1 nucleotide overlap with any gene
  AlignToGene <- countOverlaps(GAlignment, allRegions$genes, 
                               minoverlap = minOverlap["genes"],
                               ignore.strand = T)
  GAlignment <- GAlignment[!AlignToGene > 1]
  
  if (!summarize){
    tabsToMerge <- list()
    for (i in 1:length(regions)){
      print(paste("Counting", regions[i], "..."))
      
      object <- makeTab(countOverlaps(allRegions[[regions[i]]], GAlignment, 
                                      minoverlap = minOverlap[regions[i]], 
                                      ignore.strand = T), 
                        class = regions[i])
      objectName <- paste0(regions[i], "Reads")
      tabsToMerge[[i]] <- object
    }
    
    RegionsTab <- tabsToMerge %>% reduce(full_join, by = "name") %>%
      mutate_all(~replace(., is.na(.), 0)) %>%
      mutate(biotype = allRegions$genes$biotype[match(name, allRegions$genes$gene_id)])
  } else {
    intergenicN <- sum(AlignToGene == 0)
    genesN <- sum(AlignToGene > 0)
    
    RegionsTab <- data.frame(region = c("intergenic", "genes"),
                             readsN = c(intergenicN, genesN),
                             stringsAsFactors = F)
    regions <- regions[-(which(regions == "genes"))]
    
    for (i in 1:length(regions)){
      print(paste("Counting", regions[i], "..."))
      
      readsN <- sum(countOverlaps(GAlignment, allRegions[[regions[i]]], 
                                  minoverlap = minOverlap[regions[i]], 
                                  ignore.strand = T) > 0)
      toAppend <- data.frame(region = regions[i],
                             readsN = readsN,
                             stringsAsFactors = F)
      
      RegionsTab <- rbind(RegionsTab, toAppend)
    }
  }
  
  return(RegionsTab)
}

# # Vectrorized strsplit  
splitvec <- function(vector, split, select, merge = "_"){
  processed <- sapply(vector, function(x){
    separated <- unlist(strsplit(x, split = split))[select]
    if (length(separated) > 1){
      return(paste(separated, collapse = merge))
    } else
      return(separated)
  })
  processed <- unname(processed)
  return(processed)
}

# GRAPH PARAMETERS:

# # Colour palettes
divergingPal11 <- brewer.pal(11, "Spectral")
divergingPal28 <- c("grey",
                    brewer.pal(9, "YlOrRd")[c(3,5,7,9)],
                    brewer.pal(9, "YlGnBu")[c(2,3,4,5,6,7,8,9)],
                    brewer.pal(9, "YlGn")[c(8,7,5)],
                    brewer.pal(9, "PuRd")[c(3,5,6,7)],
                    brewer.pal(9, "BuPu")[c(9,8,7,6,5,4,3,2)])

divergingPal29 <- divergingPal_long[c(5,4,3,2,6,7,16:14,13:8,17:30)]
divergingPal26<- divergingPal_long[c(5,4,3,2,13:8,7,16:14,17:28)]

# # Ggplot themes
nature_barplot <- function(base_size = 11,
                           base_family = "",
                           base_line_size = base_size / 170,
                           base_rect_size = base_size / 170){
  theme_classic(base_size = base_size, 
                base_family = base_family,
                base_line_size = base_line_size) %+replace%
    theme(panel.border= element_blank(),
          axis.title.x = element_text(color="black", size=15, margin = margin(t=0, r=5, b=0, l=0)),
          axis.title.y = element_text(color="black", size=15, angle = 90, margin = margin(t=0, r=5, b=0, l=0)),
          axis.text.y = element_text(color="black", size=15, hjust=0.95,vjust=0.2),
          axis.text.x = element_text(color="black", size=15),
          #axis.line.x = element_line(color="black", size = 0.5),
          #axis.line.y = element_line(color="black", size = 0.5, hjust = 1),
          axis.ticks.y = element_blank(),
          #axis.ticks.x = element_blank(),
          legend.title = element_text(color="black", size=15),
          legend.text = element_text(color="black", size=15),
          legend.position = "right",
          strip.text = element_text(color="black", size=15, margin = margin(2,0,2,0, "mm")))
} 