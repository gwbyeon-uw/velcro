library(data.table)
library(limma)
library(edgeR)
library(locfdr)
library(biomaRt)
library(GenomicFeatures)
library(Biostrings)
library(tidyverse)
library(wiggleplotr)

options(stringsAsFactors=F)

binCounts=data.frame(fread("./dat/winCounts.gz",sep="\t",header=T))
rownames(binCounts)=binCounts$id
binCounts.noin=binCounts[,-10] #Drop input sample column
binCounts.mat=as.matrix(binCounts.noin[,-(1:6)])

expr=DGEList(counts=binCounts.mat)
expr=calcNormFactors(expr,method="TMM") #TMM scaling factors
normFactors=expr$samples$norm.factors
names(normFactors)=rownames(expr$sample)

samples=c("ES","ES","ES","WT","WT","WT") #Sample labels
design=model.matrix(~0+samples)  
colnames(design)=unique(samples)
v=voom(expr[rowSums(expr$counts)>=30,],plot=T,design=design,span=0.1) #Voom; read count >=30

cont=makeContrasts(ES-WT,ES,WT,levels=design)
fit=lmFit(v,design)
fit2=contrasts.fit(fit,cont)
fit2.eb=eBayes(fit2)

pos=binCounts[rownames(fit),1:6] #Positions
rownames(pos)=pos$name

tt=topTable(fit2.eb,n=Inf,coef = 1,sort.by = "none",confint=T) #Summary table
lfdr=locfdr(tt$t,df=25,bre=150,mlests=c(-0.5,1.0)) #local FDR
tt$locfdr=lfdr$fdr
tt$bhfdr=p.adjust(tt$P.Value,method="BH") #Benjamini Hochberg
tt=cbind(tt,fit2.eb[rownames(tt),]$coefficients,binCounts[rownames(tt),c(1:4,6)])

#Let's add annotations
ensembl_mart = useMart("ENSEMBL_MART_ENSEMBL",host="http://jul2018.archive.ensembl.org")
ensembl_dataset = useDataset("mmusculus_gene_ensembl",mart=ensembl_mart)
annotm=getBM(c("mgi_symbol","ensembl_gene_id","ensembl_transcript_id"),mart=ensembl_dataset)
annotm.nodup=annotm[!duplicated(annotm$ensembl_transcript_id),] #ENST->ENSG,MGI
rownames(annotm.nodup)=annotm.nodup$ensembl_transcript_id
annotm.nodup.ensg=unique(annotm.nodup[,1:2])
rownames(annotm.nodup.ensg)=annotm.nodup.ensg$ensembl_gene_id #ENSG->MGI

#txdb=makeTxDbFromBiomart(biomart="ENSEMBL_MART_ENSEMBL",host="http://jul2018.archive.ensembl.org",dataset="mmusculus_gene_ensembl")
#saveDb(txdb, "./annotations/txdb.rds")
txdb = loadDb("./dat/txdb.rds") #load txdb
seqlevels(txdb)=seqlevels0(txdb)
seqlevels(txdb)=as.character(c(1:19,"X","Y","MT"))
seqlevels(txdb)[seqlevels(txdb)%in%as.character(c(1:19,"X","Y","MT"))]=c(paste("chr",seqlevels(txdb)[seqlevels(txdb)%in%as.character(c(1:19,"X","Y"))],sep=""),"chrM") #add chr to chromosome names

tx=as.data.frame(transcripts(txdb,columns=c("tx_name","gene_id","tx_id"))) #transcript annotation table
tx$gene_id=unlist(tx$gene_id)
rownames(tx)=tx$tx_name

utr5=fiveUTRsByTranscript(txdb) #5'UTRs
utr3=threeUTRsByTranscript(txdb) #3'UTRs
cds=cdsBy(txdb) #CDS
exons=exonsBy(txdb,by="gene") #Exons

#Find overlaps with annotated regions
pos.gr=GRanges(seqnames=Rle(pos[,1]),IRanges(start=pos[,2],end=pos[,3]),Rle(strand(pos[,6]))) #GRange object for all positions in data table
pos.gr.utr5=as.data.frame(findOverlaps(pos.gr,utr5))
pos.gr.utr3=as.data.frame(findOverlaps(pos.gr,utr3))
pos.gr.cds=as.data.frame(findOverlaps(pos.gr,cds))
pos.gr.exon=as.data.frame(findOverlaps(pos.gr,exons))
pos.gr.exon.nobound=as.data.frame(findOverlaps(pos.gr,exons,type="within"))
pos.gr.exon.nodup=pos.gr.exon[!duplicated(pos.gr.exon$queryHits),]
rownames(pos.gr.exon.nodup)=pos.gr.exon.nodup$queryHits
utr5overlap=unique(pos.gr.utr5$queryHits)
utr3overlap=unique(pos.gr.utr3$queryHits)
cdsoverlap=unique(pos.gr.cds$queryHits)
exonoverlap=unique(as.data.frame(pos.gr.exon)$queryHits)
exonoverlap.nobound=unique(as.data.frame(pos.gr.exon.nobound)$queryHits)

#Update summary table
tt$overlapscds=F
tt$overlaps5utr=F
tt$overlaps3utr=F
tt$overlapsexon=F
tt$overlapsexonbound=F
tt$overlapscds[cdsoverlap]=T
tt$overlaps3utr[utr3overlap]=T
tt$overlaps5utr[utr5overlap]=T
tt$overlapsexon[exonoverlap]=T
tt$overlapsexonbound[exonoverlap.nobound]=T
tt$nearest=NA

useid=as.character(1:nrow(tt))[as.character(1:nrow(tt))%in%rownames(pos.gr.exon.nodup)]
tt[as.numeric(useid),"nearest"]=names(exons)[pos.gr.exon.nodup[useid,"subjectHits"]]
tt$symbol=NA
tt$symbol[!is.na(tt$nearest)]=annotm.nodup.ensg[tt$nearest[!is.na(tt$nearest)],"mgi_symbol"]
tt$ensg=tt$nearest #Nearest ensemble gene id for each window

outTable=tt[,c("name","chr","start","end","strand","overlaps5utr","overlaps3utr","overlapscds","overlapsexon","ensg","symbol","ES","WT","logFC","CI.L","CI.R","t","P.Value","locfdr","bhfdr")]
outTable.order=outTable[order(outTable$P.Value),]
outTable.order.format=outTable.order[,-c(1)]
colnames(outTable.order.format)=c("Chromosome","Start","End","Strand","Overlaps 5'UTR","Overlaps 3'UTR","Overlaps CDS","Overlaps Exon", "Overlapping Ensembl Gene ID","Overlapping Ensembl Gene Symbol","ES","WT","Log2 Fold Change","Left 95% CI","Right 95% CI","t-statistic","p-value","local FDR","Benjamini-Hochberg FDR")
write.table(outTable.order.format,file="outTable.genome.tsv",sep="\t",row.names=F,col.names=T,quote=F)

#Plotting
txlist=as.data.frame(transcripts(txdb))$tx_name
exons.name = exonsBy(txdb, by = "tx", use.names = TRUE)
cds.name = cdsBy(txdb, by = "tx", use.names = TRUE)
seqlevels(exons.name)[seqlevels(exons.name)%in%as.character(c(1:19,"X","Y"))]=paste("chr",seqlevels(exons.name)[seqlevels(exons.name)%in%as.character(c(1:19,"X","Y"))],sep="")
seqlevels(cds.name)[seqlevels(cds.name)%in%as.character(c(1:19,"X","Y"))]=paste("chr",seqlevels(cds.name)[seqlevels(cds.name)%in%as.character(c(1:19,"X","Y"))],sep="")

selected_attributes = c("ensembl_transcript_id", "ensembl_gene_id", "external_gene_name", "strand", "gene_biotype", "transcript_biotype")
transcript_metadata=getBM(attributes = selected_attributes, mart=ensembl_dataset)
transcript_metadata = dplyr::rename(transcript_metadata, transcript_id = ensembl_transcript_id, gene_id = ensembl_gene_id, gene_name = external_gene_name)

#Counts data for plotting
sample_data=dplyr::data_frame(sample_id=gsub("\\.","-",colnames(binCounts)[-(1:6)]),condition=factor(c("ES","ES","ES","input","WT","WT","WT"),levels=c("ES","WT","input")),scaling_factor=c(1,1,1,1,1,1,1))
sample_data=dplyr::mutate(sample_data, track_id = condition, colour_group = condition)
sample_data_plus = sample_data %>% dplyr::mutate(bigWig =  paste0(getwd(),"/bw/",sample_id, ".plus.bw"))
sample_data_minus = sample_data %>% dplyr::mutate(bigWig =  paste0(getwd(),"/bw/",sample_id, ".minus.bw"))

utr5Lengths=transcriptLengths(txdb,with.utr5_len=T)
rownames(utr5Lengths)=utr5Lengths$tx_name

#Scaling factor for plotting
exprAll=DGEList(counts=binCounts[,-c(1:6)])
exprAll=calcNormFactors(exprAll,method="TMM")
normFactorsAll=(exprAll$samples$lib.size*exprAll$samples$norm.factors)
libScalingFactors=1/normFactorsAll*10^6

plotCov=function(plotsymbol) {
  selected_transcripts = transcript_metadata %>%dplyr::filter(gene_name == plotsymbol,transcript_biotype == "protein_coding")

  tx_ids = intersect(selected_transcripts$transcript_id, txlist)
  if(length(tx_ids)==0) {
    return(NA)
  }
  plotstrand=as.character(strand(exons.name[tx_ids[1]])[[1]])[1]
  colourpalette=c("#d82433","#183f6d","#000000")
  utr5lengths=utr5Lengths[tx_ids,"utr5_len"]
  
  inputScaling=max(libScalingFactors[-4])/libScalingFactors[4]/(max(binCounts[rownames(subset(tt,symbol==plotsymbol)),"input.TTGCGTAC"])/max(apply(binCounts.mat[rownames(subset(tt,symbol==plotsymbol)),,drop=F],2,max)))
  
  if(plotstrand=="+") {
    sample_data_use=sample_data_plus
  }
  else {
    sample_data_use=sample_data_minus
  }
  
  sample_data_use$scaling_factor[4]=1/inputScaling
  
  plotheights=c(1,1)
  plotregion=c(as.character(seqnames(range(GenomicRanges::reduce(unlist(exons.name[tx_ids]))))),start(range(GenomicRanges::reduce(unlist(exons.name[tx_ids]))))-100,end(range(GenomicRanges::reduce(unlist(exons.name[tx_ids]))))+100,as.character(strand(range(GenomicRanges::reduce(unlist(exons.name[tx_ids]))))))
  plotobj=list()
  plotobj[["cov"]]=plotCoverage(exons.name[tx_ids],cds.name[tx_ids] , transcript_metadata, sample_data_use, rescale_introns = T ,heights=plotheights,transcript_label=T,mean_only=F,alpha=0.3,fill_palette=colourpalette,coverage_type="both",plot_fraction=1,return_subplots_list=F,region_coords=as.numeric(plotregion[2:3]))
  plotobj[["region"]]=plotregion
  return(plotobj)
}

#Plot these genes
plotlist=unique(read.csv("./plotlist.txt",header=F)[,1])
trackplots=list()
for(plotsymbol in plotlist) {
  trackplots[[plotsymbol]]=plotCov(plotsymbol)
  plotregion=trackplots[[plotsymbol]][[2]]
  pdf(paste0("./",plotsymbol,"_",paste(plotregion,collapse="_"),".cov.pdf"),width=8,height=6)
  print(trackplots[[plotsymbol]][1])
  dev.off()
}