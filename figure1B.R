
setwd("C:\\Users\\86159\\Desktop\\TCGA\\Work\\05pairedplot")                         
group = read.table("samplegroup.txt",header=F,sep="\t")                             
df = read.table("pairedInput.txt",row.names=1,header=T,sep="\t",check.names=F)       
m = match(group[,1],rownames(df))
df = df[m,]


Plot = function(data,group,outpdf){
	xfactors = as.factor(group[,2])
	xnumsample = as.numeric(xfactors)
	xaxis = levels(xfactors)
	link = group[,3]
	links = unique(group[,3])
	x1data = data[xnumsample==1]
	x2data = data[xnumsample==2]
  wilcoxP=wilcox.test(x1data,x2data)$p.value
  pvalue=signif(wilcoxP,4)
  if(pvalue<0.001){
     pvalue=signif(pvalue,4)
     pvalue=format(pvalue, scientific = TRUE)
  }else{
     pvalue=round(pvalue,3)
  }
	pdf(outpdf,width=6,height=5)
	par(las=1)
	plot(1,xlim=c(0.5,2.5),ylim=c(0,max(data)*1.2),type="n",xlab="",ylab="",xaxt="n")
	points(rep(1,length(x1data)),x1data,pch=16,cex=2,col="blue")
	points(rep(2,length(x2data)),x2data,pch=15,cex=2,col="red")
	axis(1,1:2,xaxis)
	for(i in links){
		w = which(link==i)
		x1 = xnumsample[w[1]]
		y1 = data[w[1]]
		x2 = xnumsample[w[2]]
		y2 = data[w[2]]
		segments(x1,y1,x2,y2)
	}
	par(xpd=T)
	arrows(1,max(data)*1.1,2,max(data)*1.1,angle=90,code=3,length=0.1)
	text(1.5,max(data)*1.1,paste("p =",pvalue),pos=3,cex=1)
	dev.off()
}

data = df[,1]
cell = colnames(df)[1]
cell = gsub(' ','_',cell)
outpdf = paste0(cell,".pdf")
Plot(data,group,outpdf)

