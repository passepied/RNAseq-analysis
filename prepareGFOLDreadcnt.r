##This R function is used to convert a RNAseq count matrix into several .read_cnt files, which can
##be used by GFOLD.

##To make this R function work, you will need:
##1. count matrix 
##2. condition is the name of condition for each column of count matrix
##3. conversion as the conversion table to look up the gene name of refseq genes, which is downloaded from UCSC genome browser https://genome.ucsc.edu/cgi-bin/hgTables?hgsid=389294497_dNrNdo8TdPUFtoHRCLpqeZkYo68n&clade=mammal&org=Human&db=hg19&hgta_group=genes&hgta_track=refGene&hgta_table=0&hgta_regionType=genome&position=chr1%3A148549986-148555668&hgta_outputType=primaryTable&hgta_outFileName=refgene_symbol.txt
##4. exon_length matrix with only one columns storing the exon length for each transcript

prepareGFOLDreadcnt<-function(data=data, condition=condition, conversion=conversion,exon_length=exon_length)
{
	#source("makeGFOLDdat.r")
	library(edgeR)
	y<-DGEList(counts=data,group=condition)

	y<-calcNormFactors(y)
	y1<-round(cpm(y))
	#Filter genes with cpm lower than 1.
	y1<-y1[rowSums(y1)>=1,] 
	colnames(y1)<-condition
	y1<-y1[intersect(x=rownames(y1),rownames(exon_length)),]
	for (i in 1:dim(y1)[2])
	{
		result<-makeGFOLDdat(count.dat=y1,ref.dat=conversion,exon.length=exon_length,i=i)
		write.table(result,file=sprintf("%s.read_cnt",colnames(y1)[i]),quote=F,sep="\t",col.names=F)
	}

}

