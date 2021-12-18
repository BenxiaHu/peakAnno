#' @title peakAnno: a package annotating peaks from m6A-seq and ChIP-seq/ATAC-seq/Cut&Tag.
#' @description  This function peakAnno() in the peakAnno package.
#' @param GTF a GTF file used for peakAnno.
#' @param organism peaks called from which species.
#' @param up the start position of gene promoter used for peakAnno.
#' @param down the end position of gene promoter used for peakAnno.
#' @param peaktype peaks must include chr,start,end, and strand (+/-/.) used for peakAnno.
#' @param bedfile a peak file used for peakAnno.
#' @param outpath the path saving the output file.
#' @param outfile the file name for the output file.
#' @usage peakAnno(GTF,organism,up,down,peaktype,bedfile,outpath,outfile).
#' @details peakAnno annotate peaks and return the genomic feature associated with peaks.
#' @return A peak annotation file will be returned at default.
#' @note GTF files are available https://www.gencodegenes.org/.
#' @note up: upstream of TSS; default=2000.
#' @note down: downstream of TSS; default=0.
#' @note peaktype: default=m6A.  
#' @note outpath: pathway of output file.
#' @note outfile: output file name.
#' @author Benxia Hu & Chunhui Gao
#' @return One file will be returned at default.
#' @examples
#' # run the function
#' # peakAnno(GTF,organism,up,down,peaktype,bedfile,outpath,outfile)
#' @export
"peakAnno"
options(warn=-1)
peakAnno <- function(GTF=GTF,organism=organism,up=2000,down=0,peaktype='m6A',bedfile,outpath,outfile="m6A_anno") {
    message("make an annotation file containing transcriptid,gene symbol and gene type.")
    gtf <- as.character(GTF)
    GTF <- data.table::fread(gtf,sep="\t",header=F)
    GTF <- GTF[GTF$V3 == 'transcript',]
    df <- data.frame(V9=GTF[,9])
    gene_anno <- df %>% tidyr::separate("V9", c(NA,"transcript_id","gene_type","gene_name",NA,NA,NA,NA,NA,NA,NA,NA),sep=";")
    df <- data.frame(geneinfo=gene_anno[,3])
    temp <- df %>% separate("geneinfo", c("A","gene_name"),sep='\\"')
    gene_anno$gene_name <- temp[,2]
    df <- data.frame(geneinfo=gene_anno[,2])
    temp <- df %>% separate("geneinfo", c("A","gene_type"),sep='\\"')
    gene_anno$gene_type <- temp[,2]
    df <- data.frame(geneinfo=gene_anno[,1])
    temp <- df %>% separate("geneinfo", c("A","transcript_id"),sep='\\"')
    gene_anno$transcript_id <- temp[,2]

    message("make an annotation file")
    txdb <- makeTxDbFromGFF(gtf, organism = organism, format = "gtf")

    message("define promoter regions")
    PR <- promoters(txdb, upstream=as.numeric(up), downstream=as.numeric(down))
    names(PR) <- NULL
    PR <- as.data.frame(PR)[,c(1:3,5,7)]
    PR$typeid <- 'promoter'

    genomicFeature <- PR
    getfeature <- function(txdbid,featureid){
        if(featureid == 'exon'){
            getresult <- as.data.frame(exonsBy(txdbid, by = "tx",use.names=TRUE))[,c(3:5,7,2)]
        }else if(featureid == 'intron'){
            getresult <- as.data.frame(intronsByTranscript(txdbid, use.names=TRUE))[,c(3:5,7,2)]
        }else if(featureid == 'cds'){
            getresult <- as.data.frame(cdsBy(txdbid, by = "tx",use.names=TRUE))[,c(3:5,7,2)]
        }else if(featureid == '5utr'){
            getresult <- as.data.frame(fiveUTRsByTranscript(txdbid, use.names=TRUE))[,c(3:5,7,2)]
        }else if(featureid == '3utr'){
            getresult <- as.data.frame(threeUTRsByTranscript(txdbid, use.names=TRUE))[,c(3:5,7,2)]
        }else{
           message("the feature id provided is not correct!")
        }
     getresult$typeid <- featureid
     colnames(getresult)[5] <- 'tx_name'
     return(getresult)
    }
    message("define other genomic features")
    Txfeature <- c('exon','intron','5utr','3utr','cds')
    for(i in 1:length(Txfeature)){
    genomicFeature <- rbind(genomicFeature,getfeature(txdb,Txfeature[i]))
    }
    genomic_region <- left_join(genomicFeature, gene_anno, by = c("tx_name" = "transcript_id"))
    genomic_region <- na.omit(genomic_region)
    genomic_region <- genomic_region[,-5]
    genomic_region <- genomic_region %>% unite("metadata", "strand":"gene_name", remove = TRUE,sep = ":")
    genomic_region <- GRanges(seqnames=genomic_region[,1],IRanges(start=genomic_region[,2],end=genomic_region[,3]),metadata=genomic_region[,4])
    genomic_region <- unlist(reduce(split(genomic_region, genomic_region$metadata)))
    a <- data.frame(geneinfo=names(genomic_region))
    names(genomic_region) <- NULL
    genomic_region <- as.data.frame(genomic_region)[,1:3]
    metadata <- a %>% separate("geneinfo", c("strand","typeid","genetype","genename"),sep=":")
    genomic_region <- cbind(genomic_region, metadata)
    colnames(genomic_region)[1] <- c('chr')

    message("annotation is running")
    ####
    if(as.character(peaktype) =='m6A'){
        genomic_region <- genomic_region[genomic_region$typeid != "intron",]
    }else if(as.character(peaktype) %in% c('ChIP-seq','ATAC-seq','Cut&Tag')){
        genomic_region <- genomic_region
    }else{
      message("the peak type provided is not correct")
    }

    annorange <- GRanges(seqnames=genomic_region[,1], IRanges(genomic_region[,2], genomic_region[,3]),strand=genomic_region[,4])
    mcols(annorange) <- genomic_region[,c(4:7)]
    input <- read.table(bedfile,sep="\t",header=T)
    inputrange <- GRanges(seqnames=input[,1], IRanges(input[,2], input[,3]),strand=input[,4])
    olap <- findOverlaps(inputrange, annorange)
    result <- inputrange[queryHits(olap)]
    mcols(result) <- cbind(mcols(result), mcols(annorange[subjectHits(olap)]))
    result <- as.data.frame(result)

    PCG <- result[(result$genetype=="protein_coding") & (result$typeid != 'exon'),]
    Non <- result[(result$genetype !="protein_coding"),]
    result <- rbind(PCG,Non)
    result <- unique((result)[,c(1:3,5,6:ncol(result))])
    colnames(result) <- c('chr_peak','start_peak','end_peak','strand_peak',
                           'strand','typeid','genetype','genename')
    result <- result %>% unite("flag", c('chr_peak','start_peak','end_peak','strand_peak',
                                         'strand','genetype','genename'), remove = TRUE,sep = ":")
    result <- result %>% group_by("flag") %>% summarise(typeid=paste(typeid,collapse=';'))
    result <- as.data.frame(result)
    df <- data.frame(geneinfo=result$flag)
    metadata <- df %>% separate("geneinfo", c('chr_peak','start_peak','end_peak','strand_peak',
                                     'strand','genetype','genename'),sep=":")
    result <- data.frame(metadata,typeid=result$typeid)
    write_xlsx(result, paste0(outpath,"/",outfile,"_anno.xlsx"))
    message("annotation is done")
}
