###Video source: http://study.163.com/provider/1026136977/index.htm?share=2&shareId=1026136977
######Video source: http://www.biowolf.cn/shop/
######生信自学网: http://www.biowolf.cn/
######合作邮箱：2749657388@qq.com
######答疑微信: 18520221056

#install.packages('survival')

setwd("C:\\Users\\86159\\Desktop\\TCGA\\Work\\12uniCox")
library(survival)
rt=read.table("coxInput.txt",header=T,sep="\t",check.names=F,row.names=1)
#rt[,"NME3"]=log2(rt[,"NME3"]+1)

outTab=data.frame()
for(i in colnames(rt[,3:ncol(rt)])){
 cox <- coxph(Surv(futime, fustat) ~ rt[,i], data = rt)
 coxSummary = summary(cox)
 coxP=coxSummary$coefficients[,"Pr(>|z|)"]
 outTab=rbind(outTab,
              cbind(id=i,
              HR=coxSummary$conf.int[,"exp(coef)"],
              HR.95L=coxSummary$conf.int[,"lower .95"],
              HR.95H=coxSummary$conf.int[,"upper .95"],
              pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
              )
}
write.table(outTab,file="uniCox.xls",sep="\t",row.names=F,quote=F)

###Video source: http://study.163.com/provider/1026136977/index.htm?share=2&shareId=1026136977
######Video source: http://www.biowolf.cn/shop/
######生信自学网: http://www.biowolf.cn/
######合作邮箱：2749657388@qq.com
######答疑微信: 18520221056