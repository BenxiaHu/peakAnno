# !/usr/bin/env Rscript
# December 12,2021
# Version 0.1.0
#
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
#' @author Benxia Hu
#' @return One file will be returned at default.
#' @import utils
#' @import data.table
#' @import tidyr
#' @import GenomicFeatures
#' @import dplyr
#' @import stats
#' @import GenomicRanges
#' @import writexl
#' @examples
#' # run the function
#' # peakAnno(GTF,organism,up,down,peaktype,bedfile,outpath,outfile)
#' @export
"peakAnno"
options(warn=-1)
peakAnno <- function(GTF=GTF,organism=organism,up=2000,down=0,peaktype='m6A',bedfile,outpath,outfile="m6A_anno") {
    print("make an annotation file containing transcriptid,gene symbol and gene type.")
    #gtf="/data/bxhu/project/database/hg38/gencode.v38.annotation.gtf"
    gtf <- as.character(GTF)
    GTF <- fread(gtf,sep="\t",header=F)
    GTF <- GTF[GTF$V3 == 'transcript',]
    df <- data.frame(V9=GTF[,9])
    hg38_anno <- df %>% separate(V9, c(NA,"transcript_id","gene_type","gene_name",NA,NA,NA,NA,NA,NA,NA,NA),sep=";")
    df <- data.frame(X=hg38_anno[,3])
    temp <- df %>% separate(X, c("A","gene_name"),sep='\\"')
    hg38_anno$gene_name <- temp[,2]
    df <- data.frame(X=hg38_anno[,2])
    temp <- df %>% separate(X, c("A","gene_type"),sep='\\"')
    hg38_anno$gene_type <- temp[,2]
    df <- data.frame(X=hg38_anno[,1])
    temp <- df %>% separate(X, c("A","transcript_id"),sep='\\"')
    hg38_anno$transcript_id <- temp[,2]

    #save(file="/data/bxhu/project/database/hg38/hg38_anno.Rdata",hg38_anno)

    print("make an annotation file")
    txdb <- makeTxDbFromGFF(gtf, organism = organism, format = "gtf")

    print("define promoter regions")
    PR <- promoters(txdb, upstream=as.numeric(up), downstream=as.numeric(down))
    names(PR) <- NULL
    PR <- as.data.frame(PR)[,c(1:3,5,7)]
    PR$type <- 'promoter'

    genomicFeature <- PR
    getfeature <- function(txdbid,featureid){
        if(featureid == 'exon'){
            getresult <- as.data.frame(exonsBy(txdbid, by = "tx",use.names=TRUE))[,c(3:5,7,2)]
        #print(head(getresult))
        }else if(featureid == 'intron'){
            getresult <- as.data.frame(intronsByTranscript(txdbid, use.names=TRUE))[,c(3:5,7,2)]
        }else if(featureid == 'cds'){
            getresult <- as.data.frame(cdsBy(txdbid, by = "tx",use.names=TRUE))[,c(3:5,7,2)]
        }else if(featureid == '5utr'){
            getresult <- as.data.frame(fiveUTRsByTranscript(txdbid, use.names=TRUE))[,c(3:5,7,2)]
        }else if(featureid == '3utr'){
            getresult <- as.data.frame(threeUTRsByTranscript(txdbid, use.names=TRUE))[,c(3:5,7,2)]
        }else{
            print("the feature id provided is not correct!")
        }
     getresult$type <- featureid
     colnames(getresult)[5] <- 'tx_name'
     return(getresult)
    }
    print("define other genomic features")
    Txfeature <- c('exon','intron','5utr','3utr','cds')
    for(i in 1:length(Txfeature)){
    genomicFeature <- rbind(genomicFeature,getfeature(txdb,Txfeature[i]))
    }
    hg38_genomic_region <- left_join(genomicFeature, hg38_anno, by = c("tx_name" = "transcript_id"))
    hg38_genomic_region <- na.omit(hg38_genomic_region)
    hg38_genomic_region <- hg38_genomic_region[,-5]
    hg38_genomic_region <- hg38_genomic_region %>% unite("metadata", strand:gene_name, remove = TRUE,sep = ":")
    hg38_genomic_region <- GRanges(seqnames=hg38_genomic_region[,1],IRanges(start=hg38_genomic_region[,2],end=hg38_genomic_region[,3]),metadata=hg38_genomic_region[,4])
    hg38_genomic_region <- unlist(reduce(split(hg38_genomic_region, hg38_genomic_region$metadata)))
    a <- data.frame(X=names(hg38_genomic_region))
    names(hg38_genomic_region) <- NULL
    hg38_genomic_region <- as.data.frame(hg38_genomic_region)[,1:3]
    metadata <- a %>% separate(X, c("strand","type","genetype","genename"),sep=":")
    hg38_genomic_region <- cbind(hg38_genomic_region, metadata)
    colnames(hg38_genomic_region)[1] <- c('chr')
    #save(hg38_genomic_region,file="/data/bxhu/project/database/hg38/hg38_genomic_region.Rdata")

    print("annotation is running")
    #load("/data/bxhu/project/database/hg38/hg38_genomic_region.Rdata")

    ####
    if(as.character(peaktype) =='m6A'){
        hg38_genomic_region <- hg38_genomic_region[hg38_genomic_region$type != "intron",]
    }else if(as.character(peaktype) %in% c('ChIP-seq','ATAC-seq','Cut&Tag')){
        hg38_genomic_region <- hg38_genomic_region
    }else{
      print("the peak type provided is not correct")
    }

    annorange <- GRanges(seqnames=hg38_genomic_region[,1], IRanges(hg38_genomic_region[,2], hg38_genomic_region[,3]),strand=hg38_genomic_region[,4])
    mcols(annorange) <- hg38_genomic_region[,c(4:7)]
    input <- read.table(bedfile,sep="\t",header=T)
    inputrange <- GRanges(seqnames=input[,1], IRanges(input[,2], input[,3]),strand=input[,4])
    #mcols(inputrange) <- input[,c(4:ncol(input))]

    olap <- findOverlaps(inputrange, annorange)
    result <- inputrange[queryHits(olap)]
    mcols(result) <- cbind(mcols(result), mcols(annorange[subjectHits(olap)]))
    result <- as.data.frame(result)

    PCG <- result[(result$genetype=="protein_coding") & (result$type != 'exon'),]
    Non <- result[(result$genetype !="protein_coding"),]
    result <- rbind(PCG,Non)
    result <- unique((result)[,c(1:3,5,6:ncol(result))])
    colnames(result) <- c('chr_peak','start_peak','end_peak','strand_peak',
                           'strand','type','genetype','genename')
    result <- result %>% unite("flag", c('chr_peak','start_peak','end_peak','strand_peak',
                                         'strand','genetype','genename'), remove = TRUE,sep = ":")
    result <- result %>% group_by(flag) %>% summarise(type=paste(type,collapse=';'))
    result <- as.data.frame(result)
    df <- data.frame(X=result$flag)
    metadata <- df %>% separate(X, c('chr_peak','start_peak','end_peak','strand_peak',
                                     'strand','genetype','genename'),sep=":")
    result <- data.frame(metadata,type=result$type)
    write_xlsx(result, paste0(outpath,"/",outfile,"_anno.xlsx"))
    print("annotation is done")
}
