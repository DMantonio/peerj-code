

setwd("C:\\Users\\lexb4\\Desktop\\TCGAimmune\\15.survival")    

picDir="picture"                                              
dir.create(picDir)

library(survival)
rt=read.table("survival.txt",header=T,sep="\t",check.names=F)
rt$futime=rt$futime/365                                        
outTab=data.frame()

for(gene in colnames(rt[,4:ncol(rt)])){
  a=rt[,gene]<=median(rt[,gene])
  diff=survdiff(Surv(futime, fustat) ~a,data = rt)
  pValue=1-pchisq(diff$chisq,df=1)
  outTab=rbind(outTab,cbind(gene=gene,pvalue=pValue))
  pValue=round(pValue,3)
  #pValue=format(pValue, scientific = TRUE)

  fit <- survfit(Surv(futime, fustat) ~ a, data = rt)
  summary(fit)

  tiff(file=paste(picDir,"\\",gene,".survival.tiff",sep=""),
       width = 14,            
       height =14,            
       units ="cm",
       compression="lzw",
       bg="white",
       res=600)
  plot(fit, 
     lwd=2,
     col=c("red","blue"),
     xlab="Time (year)",
     mark.time=T,
     ylab="Survival rate",
     main=paste(gene,"(p=", pValue ,")",sep="") )
  legend("topright", 
       c("High","Low"), 
       lwd=2, 
       col=c("red","blue"))
  dev.off()
}
write.table(outTab,file="survival.xls",sep="\t",row.names=F,quote=F)
