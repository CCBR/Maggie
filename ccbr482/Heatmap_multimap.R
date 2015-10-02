
setwd("/Users/maggiec/GitHub/Maggie/ccbr482/data/")
tab=read.table("EZH2_C.merge.sorted2.bam.alignstat")
tab=read.table("EZH2_RA.merge.sorted2.bam.alignstat")
tab=read.table("data/Input_merge2.sorted2.bam.alignstat")
tab=read.table("MYCN_C.sorted.bam.alignstat")
tab=read.table("MYCN_RA.sorted.bam.alignstat")
hg19=read.table("hg19.chrom.sizes")
tail(tab)
tab[1:10,]
x=vector()
for (i in 1:dim(tab)[1]){
x=rbind(x,unlist(strsplit(as.character(tab$V1[i]),",")))
}
x=cbind(x,tab$V3)



#testFunc <- function(a, b) a + b
str=function(a) unlist(strsplit(as.character(a),"chr"))[2]
#cbind(df,f = mapply(function(x,y) x+y, df$x, df$z) )
hg19=cbind(hg19,f = mapply(function(x) str(x), hg19$V1) )

x1=match(x[,1],hg19$f)
x2=match(x[,2],hg19$f)

y=cbind(x,hg19[x1,2],hg19[x2,2])
y=y[complete.cases(y),]
y=as.data.frame(y)
colnames(y)=c("chra","chrb","counts","lengtha","lengthb")
y[,3:5] <- apply(y[,3:5],2,as.numeric)
y=cbind(y,sum=mapply(function(x,y) x+y, y[,4],y[,5]))
y=cbind(y,perc=mapply(function(x,y) x/y*100, y[,3],y$sum))
class(y$chra) <- "character"
class(y$chrb) <- "character"
class(y$chra) <- "numeric"
class(y$chrb) <- "numeric"

chra=y$chra
chrb=y$chrb
y <- y[order(chra, chrb),]
#image(z$chra,z$chrb,z$perc)
plot.new()
  op <- par(pty = "s")
label=c(seq(1:22),"x","y")
image(xtabs(perc~chra+chrb, y),xaxt="n",yaxt="n",
           col = topo.colors(12))
axis( 1, at=seq(0,1,length.out=24), labels= as.character(label),cex.axis=0.8,las=2)
axis( 2, at=seq(0,1,length.out=24), labels= as.character(label),cex.axis=0.8,las=2)
require(fields)
image.plot(xtabs(perc~chra+chrb, y),add=TRUE,smallplot=c(.15, .17, .5, .85))
