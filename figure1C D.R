
options(stringsAsFactors=F)
library(limma)
library(ggpubr)

inputFile="symbol.txt"       
cliFile="clinical.txt"         
gene="CXCR4"                    
setwd("C:\\Users\\86159\\Desktop\\TCGA\\Work\\09clinCor")    


rt=read.table(inputFile,sep="\t",header=T,check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]


group=sapply(strsplit(colnames(data),"\\-"),"[",4)
group=sapply(strsplit(group,""),"[",1)
group=gsub("2","1",group)
data=data[,group==0]


data=rbind(data,gene=data[gene,])
exp=as.matrix(t(data[c("gene",gene),]))
rownames(exp)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*","\\1\\-\\2\\-\\3",rownames(exp))
exp=avereps(exp)

cli=read.table(cliFile,sep="\t",header=T,check.names=F,row.names=1)

samSample=intersect(row.names(exp),row.names(cli))
exp=exp[samSample,]
cli=cli[samSample,]
rt=cbind(exp,cli)

for(clinical in colnames(rt[,3:ncol(rt)])){
	data=rt[c(gene,clinical)]
	colnames(data)=c("gene","clinical")
	data=data[(data[,"clinical"]!="unknow"),]
	group=levels(factor(data$clinical))
	data$clinical=factor(data$clinical, levels=group)
	comp=combn(group,2)
	my_comparisons=list()
    for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
	boxplot=ggboxplot(data, x="clinical", y="gene", color="clinical",
	          xlab=clinical,
	          ylab=paste(gene,"expression"),
	          legend.title=clinical,
	          add = "jitter")+ 
	stat_compare_means(comparisons = my_comparisons)
	pdf(file=paste0(clinical,".pdf"),width=5.5,height=5)
	print(boxplot)
	dev.off()
}



