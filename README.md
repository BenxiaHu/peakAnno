# peakAnno

**peakAnno** is a R package to annotate genomic regions (peaks) called/defined from different High-throughput sequencing data, including ChIP-seq, ATAC-seq and MeRIP-seq. 

###install:  
devtools::install_github("https://github.com/BioEpi/peakAnno")

###usage:  
peakAnno(GTF,organism,up,down,peaktype,bedfile,outpath,outfile)  

***GTF***: GTF files are available from https://www.gencodegenes.org/.  
***organism***: the names of specie, such as homo sapiens and Mus musculus.  
***up***: the distance to the upstream of transcription start sites. default value is 2000.    
***down***: the distance to the downstream of transcription start sites. default value is 0.  
***peaktype***: what kind of sequencing data are used to call peaks.  default value is 'm6A'.  
***bedfile***: the name of peaks. The bed file must includes chromosome, start, end and strand (+/-/.). It is better to provide the absolute path of the file.  
***outpath***: the path is used to save the results.  
***outfile***: the file name is used to save the results.  