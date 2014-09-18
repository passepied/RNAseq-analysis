##Read count matrix transform into read_cnt for GFOLD

makeGFOLDdat<-function(count.dat, ref.dat, lib.size, exon.length,i)
{
	rpkm<-count.dat[,i]/exon.length[rownames(count.dat),1]/colSums(count.dat)[i]*1e9
	result<-cbind(as.vector(ref.dat[rownames(count.dat),13]),count.dat[,i],exon.length[rownames(count.dat),1],rpkm)	
	return(result)
}

