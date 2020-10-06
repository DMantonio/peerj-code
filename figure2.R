

setwd("C:\\Users\\86159\\Desktop\\TCGA\\Work\\07survival")    
gene="RHBDF1"

library(survival)
library(survminer)
rt=read.table("survival.txt",header=T,sep="\t",check.names=F)
rt$futime=rt$futime/365                                       
outTab=data.frame()

a=rt[,gene]<=median(rt[,gene])
diff=survdiff(Surv(futime, fustat) ~a,data = rt)
pValue=1-pchisq(diff$chisq,df=1)
scoreKm=cbind(gene=gene,KM=pValue)
	outTab=rbind(outTab,scoreKm)
if(pValue<0.001){
		pValue="p<0.001"
	}else{
		pValue=paste0("p=",sprintf("%.03f",pValue))
	}

fit <- survfit(Surv(futime, fustat) ~ a, data = rt)


titleName=gene
	surPlot=ggsurvplot(fit, 
						data=rt,
						conf.int=TRUE,
						pval=pValue,
						pval.size=5,
						risk.table=T,
						#ncensor.plot = TRUE,
						legend.labs=c("high","low"),
						legend.title=titleName,
						xlab="Time(years)",
						break.time.by = 1,
						risk.table.title="",
						palette=c("red", "blue"),
						risk.table.height=.25)          
	pdf(file=paste0("sur.",gene,".pdf"), width = 6.5, height = 5.5,onefile = FALSE)
	print(surPlot)
	dev.off()
summary(fit)
