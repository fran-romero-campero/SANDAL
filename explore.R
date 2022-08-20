plot.ld.sd.ll <- function(gene.id, gene.name, ld.ll.gene.expression,sd.ll.gene.expression,ld=T,sd=T,ll=T)
{
  ld.zt <- paste("ld",paste0("zt",sprintf(fmt = "%02d",seq(from=0,to=20,by=4))),sep="_")
  current.gene.expression.ld <- scale(ld.gene.expression[gene.id,c(paste(ld.zt,1,sep="_"),paste(ld.zt,2,sep="_"),paste(ld.zt,3,sep="_"))])
  
  sd.zt <- paste("sd",paste0("zt",sprintf(fmt = "%02d",seq(from=0,to=20,by=4))),sep="_")
  current.gene.expression.sd <- scale(sd.gene.expression[gene.id,c(paste(sd.zt,1,sep="_"),paste(sd.zt,2,sep="_"),paste(sd.zt,3,sep="_"))])
  
  nd.zt <- paste0("ZT",seq(from=0,to=21,by=3))
  current.gene.expression.nd <- scale(nd.gene.expression[gene.id,c(paste(nd.zt,1,sep="_"),paste(nd.zt,2,sep="_"),paste(nd.zt,3,sep="_"))])
  
  max.expr <- max(c(current.gene.expression.ld, current.gene.expression.nd, current.gene.expression.sd))
  min.expr <- min(c(current.gene.expression.ld, current.gene.expression.nd, current.gene.expression.sd))
  
  plot(x = -10,y= -10,axes=F,xlab="",ylab="",
       ylim=c(min.expr-2, max.expr),xlim=c(0,72),
       main=paste(gene.id, gene.name,sep=" - "),cex.main=2)
  
  if(ld)
    lines(x = seq(from=0,by=4,to=68),current.gene.expression.ld,type="o",lwd=3,col="blue")
  
  if(sd)
    lines(x = seq(from=0,by=4,to=68),current.gene.expression.sd,type="o",lwd=3,col="red")
  
  if(nd)
    lines(x = seq(from=0,by=3,to=69),current.gene.expression.nd,type="o",lwd=3,col="black")
  
  
  for(i in 0:2)
  {
    current.line <- 0.5
    
    if(ld)
    {
      polygon(x = c(24*i, 24*i+16, 24*i+16, 24*i),
              y=c(min.expr-current.line, min.expr-current.line, min.expr-(current.line+0.25), min.expr-(current.line +0.25)),
              lwd=2,border="blue")
      polygon(x = c(24*i+16,24*(i+1),24*(i+1),24*i+16),
              y=c(min.expr-current.line, min.expr-current.line, min.expr-(current.line+0.25), min.expr-(current.line +0.25)),
              lwd=2,border="blue",col="blue")
      current.line <- current.line + 0.5
    }
    
    if(nd)
    {
      polygon(x = c(24*i, 24*i+12, 24*i+12, 24*i),
              y=c(min.expr-current.line, min.expr-current.line, min.expr-(current.line+0.25), min.expr-(current.line +0.25)),
              lwd=2,border="black")
      polygon(x = c(24*i+12,24*(i+1),24*(i+1),24*i+12),
              y=c(min.expr-current.line, min.expr-current.line, min.expr-(current.line+0.25), min.expr-(current.line +0.25)),
              lwd=2,border="black",col="black")
      current.line <- current.line + 0.5
    }
    
    if(sd)
    {
      polygon(x = c(24*i, 24*i+8, 24*i+8, 24*i),
              y=c(min.expr-current.line, min.expr-current.line, min.expr-(current.line+0.25), min.expr-(current.line +0.25)),
              lwd=2,border="red")
      polygon(x = c(24*i+8,24*(i+1),24*(i+1),24*i+8),
              y=c(min.expr-current.line, min.expr-current.line, min.expr-(current.line+0.25), min.expr-(current.line +0.25)),
              lwd=2,border="red",col="red")
    }
  }
  
  return(0)  
}


grid.newpage()
draw.triple.venn(area1 = length(complete.ld.circadian.genes),area2 = length(complete.sd.circadian.genes),area3 = length(complete.nd.circadian.genes),n12 = length(intersect(complete.ld.circadian.genes,complete.sd.circadian.genes)),n23 = length(intersect(complete.sd.circadian.genes,complete.nd.circadian.genes)), n13 = length(intersect(complete.ld.circadian.genes,complete.nd.circadian.genes)), n123 = length(intersect(intersect(complete.sd.circadian.genes, complete.ld.circadian.genes),complete.nd.circadian.genes)),lwd = 3,category = c("LD","SD","ND"),euler.d = T,col = c("blue","red","darkgrey"),fill = c(" blue","red","darkgrey"),alpha = 0.3,cex = 2,cat.cex = 2)

set.seed(112)
data <- matrix(sample(1:30,15) , nrow=3)
colnames(data) <- c("A","B","C","D","E")
rownames(data) <- c("var1","var2","var3")


data


data <- matrix(c(length(rhythmic.genes.ld.12),length(setdiff(rhythmic.genes.ld,rhythmic.genes.ld.12)),length(non.rhythmic.ld.genes),
length(rhythmic.genes.nd.12),length(setdiff(rhythmic.genes.nd,rhythmic.genes.nd.12)),length(non.rhythmic.nd.genes),
length(rhythmic.genes.sd.12),length(setdiff(rhythmic.genes.sd,rhythmic.genes.sd.12)),length(non.rhythmic.sd.genes)),nrow=3)



length(rhythmic.genes.nd.12)/length(ostta.genes)
length(setdiff(rhythmic.genes.nd,rhythmic.genes.nd.12))/length(ostta.genes)
length(non.rhythmic.nd.genes)/length(ostta.genes)


# Get the stacked barplot
barplot(data, lwd=2, names.arg = c("LD","ND","SD"),cex.names = 2,cex.axis = 1.5,
        col=colors()[c(23,89,12)], 
        border="black",#"white", 
        space=0.25, 
        font.axis=2, 
        xlab="",ylim=c(0,8000))



data <- matrix(c(length(rhythmic.genes.ld.12),length(setdiff(rhythmic.genes.ld,rhythmic.genes.ld.12)),length(non.rhythmic.ld.genes),
                 length(rhythmic.genes.nd.12),length(setdiff(rhythmic.genes.nd,rhythmic.genes.nd.12)),length(non.rhythmic.nd.genes),
                 length(rhythmic.genes.sd.12),length(setdiff(rhythmic.genes.sd,rhythmic.genes.sd.12)),length(non.rhythmic.sd.genes)),nrow=3)

par(lwd=4)
barplot(data, lwd=4, names.arg = c("LD","ND","SD"),cex.names = 2,cex.axis = 1.5,
        col=colors()[c(23,89,12)], 
        border="black", 
        space=0.25, 
        font.axis=2, 
        xlab="",ylim=c(0,8000))
dev.off()


par(lwd=4)
barplot(data, lwd=4, names.arg = c("LD","ND","SD"),cex.names = 2,cex.axis = 1.5,
        col=c("firebrick3","red","white",
              "black","grey","white",
              "blue","lightblue","white"), 
        border="black", 
        space=0.25, 
        font.axis=2, 
        xlab="",ylim=c(0,8000))
dev.off()


colnames(ld.ll.gene.expression)
ld.zt %in% colnames(gene.expression)

c(paste(ld.zt,1,sep="_"),
  paste(ld.zt,2,sep="_"),
  paste(ld.zt,3,sep="_"),
  paste(ld.zt,4,sep="_"),
  paste(ld.zt,5,sep="_")) %in% colnames(gene.expression)
rownames(gene.expression)

plot.ld.ll <- function(gene.id, gene.name, gene.expression)
{
    ld.zt <- paste("ld",paste0("zt",sprintf(fmt = "%02d",seq(from=0,to=20,by=4))),sep="_")
    current.gene.expression.ld.ll <- gene.expression[gene.id,c(paste(ld.zt,1,sep="_"),
                                                               paste(ld.zt,2,sep="_"),
                                                               paste(ld.zt,3,sep="_"),
                                                               paste(ld.zt,4,sep="_"),
                                                               paste(ld.zt,5,sep="_"))]

    min.expression <- min(current.gene.expression.ld.ll)
    max.expression <- max(current.gene.expression.ld.ll)
    range.expression <- max.expression - min.expression
    
    expression.step <- floor(range.expression / 5)
    
    plot(current.gene.expression.ld.ll,type="o",lwd=3,col="blue",axes=F,xlab="",ylab="FPKM",
         ylim=c(min.expression-expression.step,max.expression),
         cex.lab=1.3,main=paste(gene.id, gene.name,sep=" - "),cex.main=2)
    axis(side=2,lwd=3)
    axis(side = 1,pos = min.expression - 1.1*expression.step, at = seq(from=1,to=30),
         labels = rep(paste("ZT",seq(from=0,to=20,by=4)),5),las=2,lwd=3)
    
    polygon(x = c(1,5,5,1),y=c(min.expression-expression.step/2,
                               min.expression-expression.step/2,
                               min.expression-expression.step,
                               min.expression-expression.step),lwd=2,border="blue")
    
    polygon(x = c(5,7,7,5),y=c(min.expression-expression.step/2,
                               min.expression-expression.step/2,
                               min.expression-expression.step,
                               min.expression-expression.step),lwd=2,border="blue",col="blue")
    
    polygon(x = c(7,11,11,7),y=c(min.expression-expression.step/2,
                                 min.expression-expression.step/2,
                                 min.expression-expression.step,
                                 min.expression-expression.step),lwd=2,border="blue")
    
    polygon(x = c(11,13,13,11),y=c(min.expression-expression.step/2,
                                   min.expression-expression.step/2,
                                   min.expression-expression.step,
                                   min.expression-expression.step),lwd=2,border="blue",col="blue")
    
    polygon(x = c(13,17,17,13),y=c(min.expression-expression.step/2,
                                   min.expression-expression.step/2,
                                   min.expression-expression.step,
                                   min.expression-expression.step),lwd=2,border="blue")
    
    polygon(x = c(17,19,19,17),y=c(min.expression-expression.step/2,
                                   min.expression-expression.step/2,
                                   min.expression-expression.step,
                                   min.expression-expression.step),lwd=2,border="blue",col="blue")
    
    polygon(x = c(19,23,23,19),y=c(min.expression-expression.step/2,
                                   min.expression-expression.step/2,
                                   min.expression-expression.step,
                                   min.expression-expression.step),lwd=2,border="blue")
    
    polygon(x = c(23,25,25,23),y=c(min.expression-expression.step/2,
                                   min.expression-expression.step/2,
                                   min.expression-expression.step,
                                   min.expression-expression.step),lwd=2,border="blue",col="lightblue")
    
    polygon(x = c(25,29,29,25),y=c(min.expression-expression.step/2,
                                   min.expression-expression.step/2,
                                   min.expression-expression.step,
                                   min.expression-expression.step),lwd=2,border="blue")
    polygon(x = c(29,30,30,29),y=c(min.expression-expression.step/2,
                                   min.expression-expression.step/2,
                                   min.expression-expression.step,
                                   min.expression-expression.step),lwd=2,border="blue",col="lightblue")
    
}

plot.ld.dd <- function(gene.id, gene.name, gene.expression)
{
  ld.zt <- paste("ld",paste0("zt",sprintf(fmt = "%02d",seq(from=0,to=20,by=4))),sep="_")
  current.gene.expression.ld.ll <- gene.expression[gene.id,c(paste(ld.zt,1,sep="_"),
                                                             paste(ld.zt,2,sep="_"),
                                                             paste(ld.zt,3,sep="_"),
                                                             paste(ld.zt,6,sep="_"),
                                                             paste(ld.zt,7,sep="_"))]
  
  min.expression <- min(current.gene.expression.ld.ll)
  max.expression <- max(current.gene.expression.ld.ll)
  range.expression <- max.expression - min.expression
  
  expression.step <- floor(range.expression / 5)
  
  plot(current.gene.expression.ld.ll,type="o",lwd=3,col="blue",axes=F,xlab="",ylab="FPKM",
       ylim=c(min.expression-expression.step,max.expression),
       cex.lab=1.3,main=paste(gene.id, gene.name,sep=" - "),cex.main=2)
  axis(side=2,lwd=3)
  axis(side = 1,pos = min.expression - 1.1*expression.step, at = seq(from=1,to=30),
       labels = rep(paste("ZT",seq(from=0,to=20,by=4)),5),las=2,lwd=3)
  
  polygon(x = c(1,5,5,1),y=c(min.expression-expression.step/2,
                             min.expression-expression.step/2,
                             min.expression-expression.step,
                             min.expression-expression.step),lwd=2,border="blue")
  
  polygon(x = c(5,7,7,5),y=c(min.expression-expression.step/2,
                             min.expression-expression.step/2,
                             min.expression-expression.step,
                             min.expression-expression.step),lwd=2,border="blue",col="blue")
  
  polygon(x = c(7,11,11,7),y=c(min.expression-expression.step/2,
                               min.expression-expression.step/2,
                               min.expression-expression.step,
                               min.expression-expression.step),lwd=2,border="blue")
  
  polygon(x = c(11,13,13,11),y=c(min.expression-expression.step/2,
                                 min.expression-expression.step/2,
                                 min.expression-expression.step,
                                 min.expression-expression.step),lwd=2,border="blue",col="blue")
  
  polygon(x = c(13,17,17,13),y=c(min.expression-expression.step/2,
                                 min.expression-expression.step/2,
                                 min.expression-expression.step,
                                 min.expression-expression.step),lwd=2,border="blue")
  
  polygon(x = c(17,19,19,17),y=c(min.expression-expression.step/2,
                                 min.expression-expression.step/2,
                                 min.expression-expression.step,
                                 min.expression-expression.step),lwd=2,border="blue",col="blue")
  
  polygon(x = c(19,23,23,19),y=c(min.expression-expression.step/2,
                                 min.expression-expression.step/2,
                                 min.expression-expression.step,
                                 min.expression-expression.step),lwd=2,border="blue",col="lightblue")
  
  polygon(x = c(23,25,25,23),y=c(min.expression-expression.step/2,
                                 min.expression-expression.step/2,
                                 min.expression-expression.step,
                                 min.expression-expression.step),lwd=2,border="blue",col="blue")
  
  polygon(x = c(25,29,29,25),y=c(min.expression-expression.step/2,
                                 min.expression-expression.step/2,
                                 min.expression-expression.step,
                                 min.expression-expression.step),lwd=2,border="blue",col="lightblue")
  polygon(x = c(29,30,30,29),y=c(min.expression-expression.step/2,
                                 min.expression-expression.step/2,
                                 min.expression-expression.step,
                                 min.expression-expression.step),lwd=2,border="blue",col="blue")
  
}



plot.ld.ll(gene.id = "ostta04g00450",gene.name = "MCM5", gene.expression = gene.expression)

i <- 1
par(mfrow=c(2,1))

current.gene <- complete.ld.dd.rhythmic.genes[i]
plot.ld.ll(gene.id = current.gene,gene.name = current.gene, gene.expression = gene.expression)
plot.ld.dd(gene.id = current.gene,gene.name = current.gene, gene.expression = gene.expression)
i <- i + 1

dev.off()


only.ld.ll.no.dd <- setdiff(complete.ld.ll.rhythmic.genes,rhythmic.ld.ll.dd)

i <- 1
par(mfrow=c(2,2))
par(mfrow=c(1,2))
gene.set <- setdiff(rhythmic.ld.ll.dd,rhythmic.sd.ll.dd)
gene.set <- intersect(rhythmic.genes.sd.12,intersect(rhythmic.genes.sd.ll,rhythmic.genes.sd.dd))
length(gene.set)

current.gene <- gene.set[i]
# plot.ld.ll(gene.id = current.gene,gene.name = current.gene, gene.expression = gene.expression)
# results.ld.ll[current.gene,]
# plot.ld.dd(gene.id = current.gene,gene.name = current.gene, gene.expression = gene.expression)
# results.ld.dd[current.gene,]
plot.sd.ll(gene.id = current.gene,gene.name = current.gene, gene.expression = gene.expression)
plot.sd.dd(gene.id = current.gene,gene.name = current.gene, gene.expression = gene.expression)
i <- i + 1

write.table(x=gene.set,file="sd_2_peaks_ll_dd_1_peak.tsv",quote = F,row.names = F,col.names = F)

dev.off()



length(complete.ld.dd.rhythmic.genes)
length(complete.ld.ll.rhythmic.genes)

length(non.rhythmic.ld.dd.genes)
length(non.rhythmic.ld.ll.genes)




rhythmic.ld.ll.dd <- intersect(complete.ld.dd.rhythmic.genes,complete.ld.ll.rhythmic.genes)
length(rhythmic.ld.ll.dd)
par(mfrow=c(2,1))

i <- 1
print(i)
current.gene <- rhythmic.ld.ll.dd[i]
plot.ld.ll(gene.id = current.gene,gene.name = current.gene, gene.expression = gene.expression)
plot.ld.dd(gene.id = current.gene,gene.name = current.gene, gene.expression = gene.expression)
i <- i + 1


head(results.ld.dd )

results.ld.dd["ostta01g04510",]


length(setdiff(complete.ld.ll.rhythmic.genes,complete.ld.dd.rhythmic.genes))
rhythmic.ld.ll <- setdiff(complete.ld.ll.rhythmic.genes,complete.ld.dd.rhythmic.genes)

i <- 1
print(i)
current.gene <- rhythmic.ld.ll[i]
plot.ld.ll(gene.id = current.gene,gene.name = current.gene, gene.expression = gene.expression)
plot.ld.dd(gene.id = current.gene,gene.name = current.gene, gene.expression = gene.expression)
i <- i + 1


length(setdiff(complete.ld.dd.rhythmic.genes,complete.ld.ll.rhythmic.genes))
rhythmic.ld.dd <- setdiff(complete.ld.dd.rhythmic.genes,complete.ld.ll.rhythmic.genes)

i <- 1
print(i)
current.gene <- rhythmic.ld.dd[i]
plot.ld.ll(gene.id = current.gene,gene.name = current.gene, gene.expression = gene.expression)
plot.ld.dd(gene.id = current.gene,gene.name = current.gene, gene.expression = gene.expression)
i <- i + 1


length(setdiff(non.rhythmic.ld.dd.genes, complete.ld.ll.rhythmic.genes))
rhythmic.ld <- setdiff(non.rhythmic.ld.dd.genes, complete.ld.ll.rhythmic.genes)

i <- 1
print(i)
current.gene <- rhythmic.ld[i]
plot.ld.ll(gene.id = current.gene,gene.name = current.gene, gene.expression = gene.expression)
plot.ld.dd(gene.id = current.gene,gene.name = current.gene, gene.expression = gene.expression)
i <- i + 1


data.ld.ll.dd <- matrix(c(length(rhythmic.ld.ll.dd),
                          length(rhythmic.ld.ll),
                          length(rhythmic.ld.dd),
                          length(rhythmic.ld),
                          length(non.rhythmic.ld.genes)),nrow=5)

sum(data.ld.ll.dd)

barplot(data.ld.ll.dd, lwd=2, names.arg = "",cex.names = 2,cex.axis = 1.5,
        col=colorRampPalette(c("blue", "white"))(20)[c(1,10,14,17,20)], 
        border="black", 
        font.axis=2, 
        xlab="",ylim=c(0,8000))



length(setdiff(non.rhythmic.ld.ll.genes, complete.ld.dd.rhythmic.genes))

blue.gradient <- colorRampPalette(c("blue", "white"))
colfunc(5)


length(circadian.ld)
non.rhythmic.ld.dd.genes
non.rhythmic.ld.ll.genes


rhythmic.sd.ll.dd <- intersect(complete.sd.dd.rhythmic.genes,complete.sd.ll.rhythmic.genes)
length(rhythmic.sd.ll.dd)

rhythmic.sd.ll <- setdiff(complete.sd.ll.rhythmic.genes,complete.sd.dd.rhythmic.genes)
length(rhythmic.sd.ll)

rhythmic.sd.dd <- setdiff(complete.sd.dd.rhythmic.genes,complete.sd.ll.rhythmic.genes)
length(rhythmic.sd.dd)

rhythmic.sd <- setdiff(non.rhythmic.sd.dd.genes, complete.sd.ll.rhythmic.genes)
length(rhythmic.sd)

length(non.rhythmic.sd.genes)


data.sd.ll.dd <- matrix(c(length(rhythmic.sd.ll.dd),
                          length(rhythmic.sd.ll),
                          length(rhythmic.sd.dd),
                          length(rhythmic.sd),
                          length(non.rhythmic.sd.genes)),nrow=5)


sum(data.sd.ll.dd)

barplot(data.sd.ll.dd, lwd=2, names.arg = "",cex.names = 2,cex.axis = 1.5,
        col=colorRampPalette(c("red", "white"))(20)[c(1,10,14,17,20)], 
        border="black", 
        font.axis=2, 
        xlab="",ylim=c(0,8000))


plot.sd.ll <- function(gene.id, gene.name, gene.expression)
{
  sd.zt <- paste("sd",paste0("zt",sprintf(fmt = "%02d",seq(from=0,to=20,by=4))),sep="_")
  current.gene.expression.sd.ll <- gene.expression[gene.id,c(paste(sd.zt,1,sep="_"),
                                                             paste(sd.zt,2,sep="_"),
                                                             paste(sd.zt,3,sep="_"),
                                                             paste(sd.zt,4,sep="_"),
                                                             paste(sd.zt,5,sep="_"))]
  
  min.expression <- min(current.gene.expression.sd.ll)
  max.expression <- max(current.gene.expression.sd.ll)
  range.expression <- max.expression - min.expression
  
  expression.step <- floor(range.expression / 5)
  
  plot(current.gene.expression.sd.ll,type="o",lwd=3,col="red",axes=F,xlab="",ylab="FPKM",
       ylim=c(min.expression-expression.step,max.expression),
       cex.lab=1.3,main=paste(gene.id, gene.name,sep=" - "),cex.main=2)
  axis(side=2,lwd=3)
  axis(side = 1,pos = min.expression - 1.1*expression.step, at = seq(from=1,to=30),
       labels = rep(paste("ZT",seq(from=0,to=20,by=4)),5),las=2,lwd=3)
  
  polygon(x = c(1,3,3,1),y=c(min.expression-expression.step/2,
                             min.expression-expression.step/2,
                             min.expression-expression.step,
                             min.expression-expression.step),lwd=2,border="red")
  
  polygon(x = c(3,7,7,3),y=c(min.expression-expression.step/2,
                             min.expression-expression.step/2,
                             min.expression-expression.step,
                             min.expression-expression.step),lwd=2,border="red",col="red")
  
  polygon(x = c(7,9,9,7),y=c(min.expression-expression.step/2,
                               min.expression-expression.step/2,
                               min.expression-expression.step,
                               min.expression-expression.step),lwd=2,border="red")
  
  polygon(x = c(9,13,13,9),y=c(min.expression-expression.step/2,
                                 min.expression-expression.step/2,
                                 min.expression-expression.step,
                                 min.expression-expression.step),lwd=2,border="red",col="red")
  
  polygon(x = c(13,15,15,13),y=c(min.expression-expression.step/2,
                                 min.expression-expression.step/2,
                                 min.expression-expression.step,
                                 min.expression-expression.step),lwd=2,border="red")
  
  polygon(x = c(15,19,19,15),y=c(min.expression-expression.step/2,
                                 min.expression-expression.step/2,
                                 min.expression-expression.step,
                                 min.expression-expression.step),lwd=2,border="red",col="red")
  
  polygon(x = c(19,21,21,19),y=c(min.expression-expression.step/2,
                                 min.expression-expression.step/2,
                                 min.expression-expression.step,
                                 min.expression-expression.step),lwd=2,border="red")
  
  polygon(x = c(21,25,25,21),y=c(min.expression-expression.step/2,
                                 min.expression-expression.step/2,
                                 min.expression-expression.step,
                                 min.expression-expression.step),lwd=2,border="red",col="salmon")
  
  polygon(x = c(25,27,27,25),y=c(min.expression-expression.step/2,
                                 min.expression-expression.step/2,
                                 min.expression-expression.step,
                                 min.expression-expression.step),lwd=2,border="red")
  polygon(x = c(27,30,30,27),y=c(min.expression-expression.step/2,
                                 min.expression-expression.step/2,
                                 min.expression-expression.step,
                                 min.expression-expression.step),lwd=2,border="red",col="salmon")
  
}



plot.sd.dd <- function(gene.id, gene.name, gene.expression)
{
  sd.zt <- paste("sd",paste0("zt",sprintf(fmt = "%02d",seq(from=0,to=20,by=4))),sep="_")
  current.gene.expression.sd.ll <- gene.expression[gene.id,c(paste(sd.zt,1,sep="_"),
                                                             paste(sd.zt,2,sep="_"),
                                                             paste(sd.zt,3,sep="_"),
                                                             paste(sd.zt,6,sep="_"),
                                                             paste(sd.zt,7,sep="_"))]
  
  min.expression <- min(current.gene.expression.sd.ll)
  max.expression <- max(current.gene.expression.sd.ll)
  range.expression <- max.expression - min.expression
  
  expression.step <- floor(range.expression / 5)
  
  plot(current.gene.expression.sd.ll,type="o",lwd=3,col="red",axes=F,xlab="",ylab="FPKM",
       ylim=c(min.expression-expression.step,max.expression),
       cex.lab=1.3,main=paste(gene.id, gene.name,sep=" - "),cex.main=2)
  axis(side=2,lwd=3)
  axis(side = 1,pos = min.expression - 1.1*expression.step, at = seq(from=1,to=30),
       labels = rep(paste("ZT",seq(from=0,to=20,by=4)),5),las=2,lwd=3)
  
  polygon(x = c(1,3,3,1),y=c(min.expression-expression.step/2,
                             min.expression-expression.step/2,
                             min.expression-expression.step,
                             min.expression-expression.step),lwd=2,border="red")
  
  polygon(x = c(3,7,7,3),y=c(min.expression-expression.step/2,
                             min.expression-expression.step/2,
                             min.expression-expression.step,
                             min.expression-expression.step),lwd=2,border="red",col="red")
  
  polygon(x = c(7,9,9,7),y=c(min.expression-expression.step/2,
                             min.expression-expression.step/2,
                             min.expression-expression.step,
                             min.expression-expression.step),lwd=2,border="red")
  
  polygon(x = c(9,13,13,9),y=c(min.expression-expression.step/2,
                               min.expression-expression.step/2,
                               min.expression-expression.step,
                               min.expression-expression.step),lwd=2,border="red",col="red")
  
  polygon(x = c(13,15,15,13),y=c(min.expression-expression.step/2,
                                 min.expression-expression.step/2,
                                 min.expression-expression.step,
                                 min.expression-expression.step),lwd=2,border="red")
  
  polygon(x = c(15,19,19,15),y=c(min.expression-expression.step/2,
                                 min.expression-expression.step/2,
                                 min.expression-expression.step,
                                 min.expression-expression.step),lwd=2,border="red",col="red")
  
  polygon(x = c(19,21,21,19),y=c(min.expression-expression.step/2,
                                 min.expression-expression.step/2,
                                 min.expression-expression.step,
                                 min.expression-expression.step),lwd=2,border="red",col="salmon")
  
  polygon(x = c(21,25,25,21),y=c(min.expression-expression.step/2,
                                 min.expression-expression.step/2,
                                 min.expression-expression.step,
                                 min.expression-expression.step),lwd=2,border="red",col="red")
  
  polygon(x = c(25,27,27,25),y=c(min.expression-expression.step/2,
                                 min.expression-expression.step/2,
                                 min.expression-expression.step,
                                 min.expression-expression.step),lwd=2,border="red",col="salmon")
  polygon(x = c(27,30,30,27),y=c(min.expression-expression.step/2,
                                 min.expression-expression.step/2,
                                 min.expression-expression.step,
                                 min.expression-expression.step),lwd=2,border="red",col="red")
  
}


par(mfrow=c(2,1))
i <- 1
print(i)
current.gene <- gene.i#sd.dark.repressed[i]
plot.sd.ll(gene.id = current.gene,gene.name = current.gene, gene.expression = gene.expression)
plot.sd.dd(gene.id = current.gene,gene.name = current.gene, gene.expression = gene.expression)
i <- i + 1







i <- 1
print(i)
current.gene <- sd.dark.repressed[i]
plot.sd.dd(gene.id = current.gene,gene.name = current.gene, gene.expression = gene.expression)
i <- i + 1


i <- 1
print(i)
current.gene <- sd.dark.activated[i]
plot.sd.dd(gene.id = current.gene,gene.name = current.gene, gene.expression = gene.expression)
i <- i + 1


i <- 1
print(i)
current.gene <- sd.light.repressed[i]
plot.sd.ll(gene.id = current.gene,gene.name = current.gene, gene.expression = gene.expression)
i <- i + 1


i <- 1
print(i)
current.gene <- sd.light.activated[i]
plot.sd.ll(gene.id = current.gene,gene.name = current.gene, gene.expression = gene.expression)
i <- i + 1


i <- 1
print(i)
current.gene <- ld.dark.repressed[i]
plot.ld.dd(gene.id = current.gene,gene.name = current.gene, gene.expression = gene.expression)
i <- i + 1


i <- 1
print(i)
current.gene <- ld.dark.activated[i]
plot.ld.dd(gene.id = current.gene,gene.name = current.gene, gene.expression = gene.expression)
i <- i + 1


length(ld.dark.repressed)
length(sd.dark.repressed)
ld.sd.dark.repressed <- intersect(ld.dark.repressed,sd.dark.repressed)
length(ld.sd.dark.repressed)
i <- 1
print(i)
current.gene <- ld.sd.dark.repressed[i]
plot.ld.dd(gene.id = current.gene,gene.name = current.gene, gene.expression = gene.expression)
plot.sd.dd(gene.id = current.gene,gene.name = current.gene, gene.expression = gene.expression)
i <- i + 1
a

length(sd.dark.repressed)
length(sd.light.activated)
sd.dark.repressed.light.activated <- intersect(sd.dark.repressed,sd.light.activated)
length(sd.dark.repressed.light.activated)
i <- 1
print(i)
current.gene <- sd.dark.repressed.light.activated[i]
plot.sd.ll(gene.id = current.gene,gene.name = current.gene, gene.expression = gene.expression)
plot.sd.dd(gene.id = current.gene,gene.name = current.gene, gene.expression = gene.expression)
i <- i + 1
a

length(sd.dark.activated)
length(sd.light.repressed)
sd.dark.activated.light.repressed <- intersect(sd.dark.activated,sd.light.repressed)
length(sd.dark.activated.light.repressed)
i <- 1
print(i)
current.gene <- sd.dark.activated.light.repressed[i]
plot.sd.ll(gene.id = current.gene,gene.name = current.gene, gene.expression = gene.expression)
plot.sd.dd(gene.id = current.gene,gene.name = current.gene, gene.expression = gene.expression)
i <- i + 1
a







i <- 1
print(i)
current.gene <- rhythmic.sd.ll[i]
plot.sd.ll(gene.id = current.gene,gene.name = current.gene, gene.expression = gene.expression)
plot.sd.dd(gene.id = current.gene,gene.name = current.gene, gene.expression = gene.expression)
i <- i + 1


i <- 1
print(i)
current.gene <- rhythmic.sd.dd[i]
plot.sd.ll(gene.id = current.gene,gene.name = current.gene, gene.expression = gene.expression)
plot.sd.dd(gene.id = current.gene,gene.name = current.gene, gene.expression = gene.expression)
i <- i + 1


i <- 1
print(i)
current.gene <- rhythmic.sd[i]
plot.sd.ll(gene.id = current.gene,gene.name = current.gene, gene.expression = gene.expression)
plot.sd.dd(gene.id = current.gene,gene.name = current.gene, gene.expression = gene.expression)
i <- i + 1




i <- 1
print(i)
current.gene <- rhythmic.genes.sd.12[i]
plot.sd.ll(gene.id = current.gene,gene.name = current.gene, gene.expression = gene.expression)
plot.sd.dd(gene.id = current.gene,gene.name = current.gene, gene.expression = gene.expression)
i <- i + 1



plot.photo.skoto.period <- function(gene.id, gene.name, gene.expression)
{
  sd.zt <- paste("sd",paste0("zt",sprintf(fmt = "%02d",seq(from=0,to=20,by=4))),sep="_")
  
  current.gene.expression.sd <- gene.expression[gene.id,c(paste(sd.zt,1,sep="_"),
                                                          paste(sd.zt,2,sep="_"),
                                                          paste(sd.zt,3,sep="_"))]
  
  current.gene.expression.sd.ll <- gene.expression[gene.id,c(paste(sd.zt,4,sep="_"),
                                                             paste(sd.zt,5,sep="_"))]
  
  current.gene.expression.sd.dd <- gene.expression[gene.id,c(paste(sd.zt,6,sep="_"),
                                                             paste(sd.zt,7,sep="_"))]
  
  
  
  min.expression <- min(c(current.gene.expression.sd,current.gene.expression.sd.ll,current.gene.expression.sd.dd))
  max.expression <- max(c(current.gene.expression.sd,current.gene.expression.sd.ll,current.gene.expression.sd.dd))
  range.expression <- max.expression - min.expression
  
  expression.step <- floor(range.expression / 5)
  
  plot(x=-20,y=-20,type="o",lwd=3,col="red",axes=F,xlab="",ylab="FPKM",
       ylim=c(min.expression-expression.step,max.expression),xlim=c(0,32),
       cex.lab=1.3,main=paste(gene.id, gene.name,sep=" - "),cex.main=2)
  
  lines(c(current.gene.expression.sd,current.gene.expression.sd.ll),type="o",lwd=3,col="salmon")
  lines(c(current.gene.expression.sd,current.gene.expression.sd.dd),type="o",lwd=3,col="red")
  
  axis(side=2,lwd=3)
  axis(side = 1,pos = min.expression - 1.1*expression.step, at = seq(from=1,to=30),
       labels = rep(paste("ZT",seq(from=0,to=20,by=4)),5),las=2,lwd=3)
  
  polygon(x = c(1,3,3,1),y=c(min.expression-expression.step/2,
                             min.expression-expression.step/2,
                             min.expression-expression.step,
                             min.expression-expression.step),lwd=2,border="red")
  
  polygon(x = c(3,7,7,3),y=c(min.expression-expression.step/2,
                             min.expression-expression.step/2,
                             min.expression-expression.step,
                             min.expression-expression.step),lwd=2,border="red",col="red")
  
  polygon(x = c(7,9,9,7),y=c(min.expression-expression.step/2,
                             min.expression-expression.step/2,
                             min.expression-expression.step,
                             min.expression-expression.step),lwd=2,border="red")
  
  polygon(x = c(9,13,13,9),y=c(min.expression-expression.step/2,
                               min.expression-expression.step/2,
                               min.expression-expression.step,
                               min.expression-expression.step),lwd=2,border="red",col="red")
  
  polygon(x = c(13,15,15,13),y=c(min.expression-expression.step/2,
                                 min.expression-expression.step/2,
                                 min.expression-expression.step,
                                 min.expression-expression.step),lwd=2,border="red")
  
  polygon(x = c(15,19,19,15),y=c(min.expression-expression.step/2,
                                 min.expression-expression.step/2,
                                 min.expression-expression.step,
                                 min.expression-expression.step),lwd=2,border="red",col="red")
  
  polygon(x = c(19,21,21,19),y=c(min.expression-expression.step/2,
                                 min.expression-expression.step/2,
                                 min.expression-expression.step,
                                 min.expression-expression.step),lwd=2,border="red")
  
  polygon(x = c(21,25,25,21),y=c(min.expression-expression.step/2,
                                 min.expression-expression.step/2,
                                 min.expression-expression.step,
                                 min.expression-expression.step),lwd=2,border="red",col="salmon")
  
  polygon(x = c(25,27,27,25),y=c(min.expression-expression.step/2,
                                 min.expression-expression.step/2,
                                 min.expression-expression.step,
                                 min.expression-expression.step),lwd=2,border="red")
  polygon(x = c(27,30,30,27),y=c(min.expression-expression.step/2,
                                 min.expression-expression.step/2,
                                 min.expression-expression.step,
                                 min.expression-expression.step),lwd=2,border="red",col="salmon")
  
}


plot.photo.skoto.period.1 <- function(gene.id, gene.name, gene.expression)
{
  sd.zt <- paste("sd",paste0("zt",sprintf(fmt = "%02d",seq(from=0,to=20,by=4))),sep="_")
  
  current.gene.expression.sd <- gene.expression[gene.id,c(paste(sd.zt,1,sep="_"),
                                                          paste(sd.zt,2,sep="_"),
                                                          paste(sd.zt,3,sep="_"))]
  
  current.gene.expression.sd <- current.gene.expression.sd / max(current.gene.expression.sd)
  
  current.gene.expression.sd.ll <- gene.expression[gene.id,c(paste(sd.zt,4,sep="_"),
                                                             paste(sd.zt,5,sep="_"))]
  
  
  current.gene.expression.sd.ll <- current.gene.expression.sd.ll / max(current.gene.expression.sd.ll)
  
  current.gene.expression.sd.dd <- gene.expression[gene.id,c(paste(sd.zt,6,sep="_"),
                                                             paste(sd.zt,7,sep="_"))]
  
  current.gene.expression.sd.dd <- current.gene.expression.sd.dd / max(current.gene.expression.sd.dd)
  
  min.expression <- min(c(current.gene.expression.sd,current.gene.expression.sd.ll,current.gene.expression.sd.dd))
  max.expression <- max(c(current.gene.expression.sd,current.gene.expression.sd.ll,current.gene.expression.sd.dd))
  range.expression <- max.expression - min.expression
  
  expression.step <- floor(range.expression / 5)
  
  plot(x=-20,y=-20,type="o",lwd=3,col="red",axes=F,xlab="",ylab="Normalized Expression",
       ylim=c(-0.2,1.1),xlim=c(0,32),
       cex.lab=1.3,main=paste(gene.id, gene.name,sep=" - "),cex.main=2)
  
  lines(c(current.gene.expression.sd,current.gene.expression.sd.ll),type="o",lwd=3,col="salmon")
  lines(current.gene.expression.sd,type="o",lwd=3,col="red")
  #lines(c(current.gene.expression.sd,current.gene.expression.sd.dd),type="o",lwd=3,col="red")
  
  axis(side=2,lwd=3,at = seq(from=0,to=1,by=0.2))
  axis(side = 1,pos = -0.2, at = seq(from=1,to=30),
       labels = rep(paste("ZT",seq(from=0,to=20,by=4)),5),las=2,lwd=3)
  
  polygon(x = c(1,3,3,1),y=c(0,
                             0,
                             -0.14,
                             -0.14),lwd=2,border="red")
  
  polygon(x = c(3,7,7,3),y=c(0,
                             0,
                             -0.14,
                             -0.14),lwd=2,border="red",col="red")
  
  polygon(x = c(7,9,9,7),y=c(0,
                             0,
                             -0.14,
                             -0.14),lwd=2,border="red")
  
  polygon(x = c(9,13,13,9),y=c(0,
                               0,
                               -0.14,
                               -0.14),lwd=2,border="red",col="red")
  
  polygon(x = c(13,15,15,13),y=c(0,
                                 0,
                                 -0.14,
                                 -0.14),lwd=2,border="red")
  
  polygon(x = c(15,19,19,15),y=c(0,
                                 0,
                                 -0.14,
                                 -0.14),lwd=2,border="red",col="red")
  
  polygon(x = c(19,21,21,19),y=c(0,0,-0.07,-0.07),lwd=2,border="red")
  #polygon(x = c(19,21,21,19),y=c(-0.07,-0.07,-0.14,-0.14),lwd=2,border="red",col="salmon")
  
  polygon(x = c(21,25,25,21),y=c(0,0,-0.07,-0.07),lwd=2,border="red")
  #polygon(x = c(21,25,25,21),y=c(-0.07,-0.07,-0.14,-0.14),lwd=2,border="red",col="salmon")
  
  polygon(x = c(25,27,27,25),y=c(0,0,-0.07,-0.07),lwd=2,border="red")
  #polygon(x = c(25,27,27,25),y=c(-0.07,-0.07,-0.14,-0.14),lwd=2,border="red",col="salmon")
  
  polygon(x = c(27,30,30,27),y=c(0,0,-0.07,-0.07),lwd=2,border="red")
  #polygon(x = c(27,30,30,27),y=c(-0.07,-0.07,-0.14,-0.14),lwd=2,border="red",col="salmon")
}


plot.photo.skoto.period.2 <- function(gene.id, gene.name, gene.expression)
{
  sd.zt <- paste("sd",paste0("zt",sprintf(fmt = "%02d",seq(from=0,to=20,by=4))),sep="_")
  
  current.gene.expression.sd <- gene.expression[gene.id,c(paste(sd.zt,1,sep="_"),
                                                          paste(sd.zt,2,sep="_"),
                                                          paste(sd.zt,3,sep="_"))]
  
  current.gene.expression.sd <- current.gene.expression.sd / max(current.gene.expression.sd)
  
  current.gene.expression.sd.ll <- gene.expression[gene.id,c(paste(sd.zt,4,sep="_"),
                                                             paste(sd.zt,5,sep="_"))]
  
  
  current.gene.expression.sd.ll <- current.gene.expression.sd.ll / max(current.gene.expression.sd.ll)
  
  current.gene.expression.sd.dd <- gene.expression[gene.id,c(paste(sd.zt,6,sep="_"),
                                                             paste(sd.zt,7,sep="_"))]
  
  current.gene.expression.sd.dd <- current.gene.expression.sd.dd / max(current.gene.expression.sd.dd)
  
  min.expression <- min(c(current.gene.expression.sd,current.gene.expression.sd.ll,current.gene.expression.sd.dd))
  max.expression <- max(c(current.gene.expression.sd,current.gene.expression.sd.ll,current.gene.expression.sd.dd))
  range.expression <- max.expression - min.expression
  
  expression.step <- floor(range.expression / 5)
  
  plot(x=-20,y=-20,type="o",lwd=3,col="red",axes=F,xlab="",ylab="Normalized Expression",
       ylim=c(-0.2,1.1),xlim=c(0,32),
       cex.lab=1.3,main=paste(gene.id, gene.name,sep=" - "),cex.main=2)
  
  #lines(c(current.gene.expression.sd,current.gene.expression.sd.ll),type="o",lwd=3,col="salmon")
  lines(c(current.gene.expression.sd,current.gene.expression.sd.dd),type="o",lwd=3,col="red")
  
  axis(side=2,lwd=3,at = seq(from=0,to=1,by=0.2))
  axis(side = 1,pos = -0.2, at = seq(from=1,to=30),
       labels = rep(paste("ZT",seq(from=0,to=20,by=4)),5),las=2,lwd=3)
  
  polygon(x = c(1,3,3,1),y=c(0,
                             0,
                             -0.14,
                             -0.14),lwd=2,border="red")
  
  polygon(x = c(3,7,7,3),y=c(0,
                             0,
                             -0.14,
                             -0.14),lwd=2,border="red",col="red")
  
  polygon(x = c(7,9,9,7),y=c(0,
                             0,
                             -0.14,
                             -0.14),lwd=2,border="red")
  
  polygon(x = c(9,13,13,9),y=c(0,
                               0,
                               -0.14,
                               -0.14),lwd=2,border="red",col="red")
  
  polygon(x = c(13,15,15,13),y=c(0,
                                 0,
                                 -0.14,
                                 -0.14),lwd=2,border="red")
  
  polygon(x = c(15,19,19,15),y=c(0,
                                 0,
                                 -0.14,
                                 -0.14),lwd=2,border="red",col="red")
  
  #polygon(x = c(19,21,21,19),y=c(0,0,-0.07,-0.07),lwd=2,border="red")
  polygon(x = c(19,21,21,19),y=c(-0.07,-0.07,-0.14,-0.14),lwd=2,border="red",col="salmon")
  
  #polygon(x = c(21,25,25,21),y=c(0,0,-0.07,-0.07),lwd=2,border="red")
  polygon(x = c(21,25,25,21),y=c(-0.07,-0.07,-0.14,-0.14),lwd=2,border="red",col="salmon")
  
  #polygon(x = c(25,27,27,25),y=c(0,0,-0.07,-0.07),lwd=2,border="red")
  polygon(x = c(25,27,27,25),y=c(-0.07,-0.07,-0.14,-0.14),lwd=2,border="red",col="salmon")
  
  #polygon(x = c(27,30,30,27),y=c(0,0,-0.07,-0.07),lwd=2,border="red")
  polygon(x = c(27,30,30,27),y=c(-0.07,-0.07,-0.14,-0.14),lwd=2,border="red",col="salmon")
}

plot.photo.skoto.period.3 <- function(gene.id, gene.name, gene.expression)
{
  sd.zt <- paste("sd",paste0("zt",sprintf(fmt = "%02d",seq(from=0,to=20,by=4))),sep="_")
  
  current.gene.expression.sd <- gene.expression[gene.id,c(paste(sd.zt,1,sep="_"),
                                                          paste(sd.zt,2,sep="_"),
                                                          paste(sd.zt,3,sep="_"))]
  
  current.gene.expression.sd <- current.gene.expression.sd / max(current.gene.expression.sd)
  
  current.gene.expression.sd.ll <- gene.expression[gene.id,c(paste(sd.zt,4,sep="_"),
                                                             paste(sd.zt,5,sep="_"))]
  
  
  current.gene.expression.sd.ll <- current.gene.expression.sd.ll / max(current.gene.expression.sd.ll)
  
  current.gene.expression.sd.dd <- gene.expression[gene.id,c(paste(sd.zt,6,sep="_"),
                                                             paste(sd.zt,7,sep="_"))]
  
  current.gene.expression.sd.dd <- current.gene.expression.sd.dd / max(current.gene.expression.sd.dd)
  
  min.expression <- min(c(current.gene.expression.sd,current.gene.expression.sd.ll,current.gene.expression.sd.dd))
  max.expression <- max(c(current.gene.expression.sd,current.gene.expression.sd.ll,current.gene.expression.sd.dd))
  range.expression <- max.expression - min.expression
  
  expression.step <- floor(range.expression / 5)
  
  plot(x=-20,y=-20,type="o",lwd=3,col="red",axes=F,xlab="",ylab="Normalized Expression",
       ylim=c(-0.2,1.1),xlim=c(0,32),
       cex.lab=1.3,main=paste(gene.id, gene.name,sep=" - "),cex.main=2)
  
  lines(c(current.gene.expression.sd,current.gene.expression.sd.ll),type="o",lwd=3,col="salmon")
  lines(c(current.gene.expression.sd,current.gene.expression.sd.dd),type="o",lwd=3,col="red")
  
  axis(side=2,lwd=3,at = seq(from=0,to=1,by=0.2))
  axis(side = 1,pos = -0.2, at = seq(from=1,to=30),
       labels = rep(paste("ZT",seq(from=0,to=20,by=4)),5),las=2,lwd=3)
  
  polygon(x = c(1,3,3,1),y=c(0,
                             0,
                             -0.14,
                             -0.14),lwd=2,border="red")
  
  polygon(x = c(3,7,7,3),y=c(0,
                             0,
                             -0.14,
                             -0.14),lwd=2,border="red",col="red")
  
  polygon(x = c(7,9,9,7),y=c(0,
                             0,
                             -0.14,
                             -0.14),lwd=2,border="red")
  
  polygon(x = c(9,13,13,9),y=c(0,
                               0,
                               -0.14,
                               -0.14),lwd=2,border="red",col="red")
  
  polygon(x = c(13,15,15,13),y=c(0,
                                 0,
                                 -0.14,
                                 -0.14),lwd=2,border="red")
  
  polygon(x = c(15,19,19,15),y=c(0,
                                 0,
                                 -0.14,
                                 -0.14),lwd=2,border="red",col="red")
  
  polygon(x = c(19,21,21,19),y=c(0,0,-0.07,-0.07),lwd=2,border="red")
  polygon(x = c(19,21,21,19),y=c(-0.07,-0.07,-0.14,-0.14),lwd=2,border="red",col="salmon")
  
  polygon(x = c(21,25,25,21),y=c(0,0,-0.07,-0.07),lwd=2,border="red")
  polygon(x = c(21,25,25,21),y=c(-0.07,-0.07,-0.14,-0.14),lwd=2,border="red",col="salmon")
  
  polygon(x = c(25,27,27,25),y=c(0,0,-0.07,-0.07),lwd=2,border="red")
  polygon(x = c(25,27,27,25),y=c(-0.07,-0.07,-0.14,-0.14),lwd=2,border="red",col="salmon")
  
  polygon(x = c(27,30,30,27),y=c(0,0,-0.07,-0.07),lwd=2,border="red")
  polygon(x = c(27,30,30,27),y=c(-0.07,-0.07,-0.14,-0.14),lwd=2,border="red",col="salmon")
}


plot.photo.skoto.period.1(gene.id = "ostta01g01700", gene.name="", gene.expression)
plot.photo.skoto.period.2(gene.id = "ostta01g01700", gene.name="", gene.expression)
plot.photo.skoto.period.3(gene.id = "ostta01g01700", gene.name="", gene.expression)


plot.ld.ll(gene.id = "ostta01g04430", gene.name="", gene.expression)


plot.photo.skoto.period.1(gene.id = "ostta01g04430", gene.name="", gene.expression)
plot.photo.skoto.period.2(gene.id = "ostta01g04430", gene.name="", gene.expression)
plot.photo.skoto.period.3(gene.id = "ostta01g04430", gene.name="", gene.expression)

plot.photo.skoto.period.1(gene.id = "ostta04g02250", gene.name="", gene.expression)
plot.photo.skoto.period.2(gene.id = "ostta04g02250", gene.name="", gene.expression)
plot.photo.skoto.period.3(gene.id = "ostta04g02250", gene.name="", gene.expression)


plot.photo.skoto.period(gene.id = "ostta01g02130", gene.name="", gene.expression)
plot.photo.skoto.period(gene.id = "ostta03g00470", gene.name="", gene.expression)
plot.photo.skoto.period(gene.id = "ostta03g01970", gene.name="", gene.expression)
plot.photo.skoto.period(gene.id = "ostta03g04470", gene.name="", gene.expression)
plot.photo.skoto.period(gene.id = "ostta04g02250", gene.name="", gene.expression)





plot.photo.skoto.period(gene.id = "ostta04g01990", gene.name="", gene.expression)


plot.photo.skoto.period(gene.id = "ostta04g01830", gene.name="", gene.expression)


plot.photo.skoto.period(gene.id = "ostta03g05120", gene.name="", gene.expression)

plot.photo.skoto.period(gene.id = "ostta03g05120", gene.name="", gene.expression)



plot.photo.skoto.period(gene.id = "ostta03g05090", gene.name="", gene.expression)



plot.photo.skoto.period(gene.id = "ostta03g04200", gene.name="", gene.expression)


plot.photo.skoto.period(gene.id = "ostta03g03260", gene.name="", gene.expression)


plot.photo.skoto.period(gene.id = "ostta03g02120", gene.name="", gene.expression)
plot.photo.skoto.period(gene.id = "ostta03g01550", gene.name="", gene.expression)


plot.photo.skoto.period(gene.id = "ostta03g00940", gene.name="", gene.expression)

plot.photo.skoto.period(gene.id = "ostta03g00110", gene.name="", gene.expression)


plot.photo.skoto.period(gene.id = "ostta02g04900", gene.name="", gene.expression)

plot.photo.skoto.period(gene.id = "ostta02g04550", gene.name="", gene.expression)


plot.photo.skoto.period(gene.id = "ostta02g04500", gene.name="", gene.expression)


plot.photo.skoto.period(gene.id = "ostta02g02930", gene.name="", gene.expression)
plot.photo.skoto.period(gene.id = "ostta02g02890", gene.name="", gene.expression)
plot.photo.skoto.period(gene.id = "ostta02g00930", gene.name="", gene.expression)



plot.photo.skoto.period(gene.id = "ostta01g03660", gene.name="", gene.expression)
plot.photo.skoto.period(gene.id = "ostta01g03630", gene.name="", gene.expression)
plot.photo.skoto.period(gene.id = "ostta01g01380", gene.name="", gene.expression)
plot.photo.skoto.period(gene.id = "ostta01g00470", gene.name="", gene.expression)
plot.photo.skoto.period(gene.id = "ostta01g01750", gene.name="", gene.expression)




plot.ld.ll <- function(gene.id, gene.name, gene.expression)#, mean.expression.ld, mean.expression.sd)
{
  #expression.ld <- mean.expression.ld[gene.id,]
  #expression.sd <- mean.expression.sd[gene.id,]
  
  # min.expression.ld.sd <- min(c(expression.ld,expression.sd))
  # max.expression.ld.sd <- max(c(expression.ld,expression.sd))
  # range.expression.ld.sd <- max.expression.ld.sd - min.expression.ld.sd
  # 
  # expression.step.ld.sd <- floor(range.expression.ld.sd / 10)
  # 
  # png(filename = paste(paste(c("mean",gene.name,gene.id),collapse="_"),".png",sep=""))
  # plot(mean.expression.ld[gene.id,],type="o",col="blue",lwd=3,axes=F,xlab="",ylab="FPKM",
  #      ylim=c(min.expression.ld.sd-2*expression.step.ld.sd,max.expression.ld.sd),cex.lab=1.3,main=paste(gene.id, gene.name,sep=" - "),cex.main=2)
  # lines(mean.expression.sd[gene.id,],type="o",col="red",lwd=3)
  # 
  # axis(side=2,lwd = 3)
  # axis(side = 1,pos = min.expression.ld.sd - 2.2 * expression.step.ld.sd, at = seq(from=1,to=6),
  #      labels = paste("ZT",seq(from=0,to=20,by=4)),las=2,lwd=3)
  # 
  # 
  # polygon(x = c(1,5,5,1),y=c(min.expression.ld.sd-expression.step.ld.sd/2,
  #                            min.expression.ld.sd-expression.step.ld.sd/2,
  #                            min.expression.ld.sd-expression.step.ld.sd,
  #                            min.expression.ld.sd-expression.step.ld.sd),lwd=2,border="blue")
  # 
  # polygon(x = c(5,6,6,5),y=c(min.expression.ld.sd-expression.step.ld.sd/2,
  #                            min.expression.ld.sd-expression.step.ld.sd/2,
  #                            min.expression.ld.sd-expression.step.ld.sd,
  #                            min.expression.ld.sd-expression.step.ld.sd),lwd=2,border="blue",col="blue")
  # 
  # polygon(x = c(1,3,3,1),y=c(min.expression.ld.sd-(3/2)*expression.step.ld.sd,
  #                            min.expression.ld.sd-(3/2)*expression.step.ld.sd,
  #                            min.expression.ld.sd-2*expression.step.ld.sd,
  #                            min.expression.ld.sd-2*expression.step.ld.sd),lwd=2,border="red")
  # 
  # polygon(x = c(3,6,6,3),y=c(min.expression.ld.sd-(3/2)*expression.step.ld.sd,
  #                            min.expression.ld.sd-(3/2)*expression.step.ld.sd,
  #                            min.expression.ld.sd-2*expression.step.ld.sd,
  #                            min.expression.ld.sd-2*expression.step.ld.sd),lwd=2,border="red",col="red")
  # dev.off()
  
  
  current.gene.expression.ld.ll <- gene.expression[gene.id,1:30]
  
  min.expression <- min(current.gene.expression.ld.ll)
  max.expression <- max(current.gene.expression.ld.ll)
  range.expression <- max.expression - min.expression
  
  expression.step <- floor(range.expression / 5)
  
  #png(filename = paste(paste(c("ld_five_days",gene.name,gene.id),collapse="_"),".png",sep=""))
  plot(current.gene.expression.ld.ll,type="o",lwd=3,col="blue",axes=F,xlab="",ylab="FPKM",
       ylim=c(min.expression-expression.step,max.expression),
       cex.lab=1.3,main=paste(gene.id, gene.name,sep=" - "),cex.main=2)
  axis(side=2,lwd=3)
  axis(side = 1,pos = min.expression - 1.1*expression.step, at = seq(from=1,to=30),
       labels = rep(paste("ZT",seq(from=0,to=20,by=4)),5),las=2,lwd=3)
  
  polygon(x = c(1,5,5,1),y=c(min.expression-expression.step/2,
                             min.expression-expression.step/2,
                             min.expression-expression.step,
                             min.expression-expression.step),lwd=2,border="blue")
  
  polygon(x = c(5,7,7,5),y=c(min.expression-expression.step/2,
                             min.expression-expression.step/2,
                             min.expression-expression.step,
                             min.expression-expression.step),lwd=2,border="blue",col="blue")
  
  polygon(x = c(7,11,11,7),y=c(min.expression-expression.step/2,
                               min.expression-expression.step/2,
                               min.expression-expression.step,
                               min.expression-expression.step),lwd=2,border="blue")
  
  polygon(x = c(11,13,13,11),y=c(min.expression-expression.step/2,
                                 min.expression-expression.step/2,
                                 min.expression-expression.step,
                                 min.expression-expression.step),lwd=2,border="blue",col="blue")
  
  polygon(x = c(13,17,17,13),y=c(min.expression-expression.step/2,
                                 min.expression-expression.step/2,
                                 min.expression-expression.step,
                                 min.expression-expression.step),lwd=2,border="blue")
  
  polygon(x = c(17,19,19,17),y=c(min.expression-expression.step/2,
                                 min.expression-expression.step/2,
                                 min.expression-expression.step,
                                 min.expression-expression.step),lwd=2,border="blue",col="blue")
  
  polygon(x = c(19,23,23,19),y=c(min.expression-expression.step/2,
                                 min.expression-expression.step/2,
                                 min.expression-expression.step,
                                 min.expression-expression.step),lwd=2,border="blue")
  
  polygon(x = c(23,25,25,23),y=c(min.expression-expression.step/2,
                                 min.expression-expression.step/2,
                                 min.expression-expression.step,
                                 min.expression-expression.step),lwd=2,border="blue",col="lightblue")
  
  polygon(x = c(25,29,29,25),y=c(min.expression-expression.step/2,
                                 min.expression-expression.step/2,
                                 min.expression-expression.step,
                                 min.expression-expression.step),lwd=2,border="blue")
  polygon(x = c(29,30,30,29),y=c(min.expression-expression.step/2,
                                 min.expression-expression.step/2,
                                 min.expression-expression.step,
                                 min.expression-expression.step),lwd=2,border="blue",col="lightblue")
  
  #dev.off()
  return(0)  
  
  
  
  sd.zt <- paste("sd",paste0("zt",sprintf(fmt = "%02d",seq(from=0,to=20,by=4))),sep="_")
  current.gene.expression.sd.ll <- gene.expression[gene.id,c(paste(sd.zt,1,sep="_"),paste(sd.zt,2,sep="_"),paste(sd.zt,3,sep="_"),paste(sd.zt,4,sep="_"),paste(sd.zt,5,sep="_"))]
  
  min.expression <- min(current.gene.expression.sd.ll)
  max.expression <- max(current.gene.expression.sd.ll)
  range.expression <- max.expression - min.expression
  
  expression.step <- floor(range.expression / 5)
  
  plot(current.gene.expression.sd.ll,type="o",lwd=3,col="red",axes=F,xlab="",ylab="FPKM",
       ylim=c(min.expression-expression.step,max.expression),
       cex.lab=1.3,main=paste(gene.id, gene.name,sep=" - "),cex.main=2)
  axis(side=2,lwd=3)
  axis(side = 1,pos = min.expression - 1.1*expression.step, at = seq(from=1,to=30),
       labels = rep(paste("ZT",seq(from=0,to=20,by=4)),5),las=2,lwd=3)
  
  polygon(x = c(1,3,3,1),y=c(min.expression-expression.step/2,
                             min.expression-expression.step/2,
                             min.expression-expression.step,
                             min.expression-expression.step),lwd=2,border="red")
  
  polygon(x = c(3,7,7,3),y=c(min.expression-expression.step/2,
                             min.expression-expression.step/2,
                             min.expression-expression.step,
                             min.expression-expression.step),lwd=2,border="red",col="red")
  
  polygon(x = c(7,9,9,7),y=c(min.expression-expression.step/2,
                             min.expression-expression.step/2,
                             min.expression-expression.step,
                             min.expression-expression.step),lwd=2,border="red")
  
  polygon(x = c(9,13,13,9),y=c(min.expression-expression.step/2,
                               min.expression-expression.step/2,
                               min.expression-expression.step,
                               min.expression-expression.step),lwd=2,border="red",col="red")
  
  polygon(x = c(13,15,15,13),y=c(min.expression-expression.step/2,
                                 min.expression-expression.step/2,
                                 min.expression-expression.step,
                                 min.expression-expression.step),lwd=2,border="red")
  
  polygon(x = c(15,19,19,15),y=c(min.expression-expression.step/2,
                                 min.expression-expression.step/2,
                                 min.expression-expression.step,
                                 min.expression-expression.step),lwd=2,border="red",col="red")
  
  polygon(x = c(19,21,21,19),y=c(min.expression-expression.step/2,
                                 min.expression-expression.step/2,
                                 min.expression-expression.step,
                                 min.expression-expression.step),lwd=2,border="red")
  
  polygon(x = c(21,25,25,21),y=c(min.expression-expression.step/2,
                                 min.expression-expression.step/2,
                                 min.expression-expression.step,
                                 min.expression-expression.step),lwd=2,border="red",col="lightsalmon")
  
  polygon(x = c(25,27,27,25),y=c(min.expression-expression.step/2,
                                 min.expression-expression.step/2,
                                 min.expression-expression.step,
                                 min.expression-expression.step),lwd=2,border="red")
  polygon(x = c(27,30,30,27),y=c(min.expression-expression.step/2,
                                 min.expression-expression.step/2,
                                 min.expression-expression.step,
                                 min.expression-expression.step),lwd=2,border="red",col="lightsalmon")
  
}

plot.ld.ll(gene.id = "ostta09g02190",gene.name = "ostta09g02190", gene.expression)

plot.ld.ll(gene.id = "ostta02g01510",gene.name = "HDAC", gene.expression)


colnames(gene.expression)

plot.ld <- function(gene.id, gene.name, gene.expression)#, mean.expression.ld, mean.expression.sd)
{

  current.gene.expression.ld.ll <- gene.expression[gene.id,1:18]#30]
  
  min.expression <- min(current.gene.expression.ld.ll)
  max.expression <- max(current.gene.expression.ld.ll)
  range.expression <- max.expression - min.expression
  
  expression.step <- floor(range.expression / 5)
  
#  png(filename = paste(paste(c("ld_three_days",gene.name,gene.id),collapse="_"),".png",sep=""),width = 1000)
  plot(current.gene.expression.ld.ll,type="o",lwd=5,col="blue",axes=F,xlab="",ylab="FPKM",
       ylim=c(min.expression-expression.step,max.expression),
       cex.lab=1.3,main=gene.id,cex.main=2)
  axis(side=2,lwd=3)
  axis(side = 1,pos = min.expression - 1.1*expression.step, at = seq(from=1,to=30),
       labels = rep(paste("ZT",seq(from=0,to=20,by=4)),5),las=2,lwd=3)
  
  polygon(x = c(1,5,5,1),y=c(min.expression-expression.step/2,
                             min.expression-expression.step/2,
                             min.expression-expression.step,
                             min.expression-expression.step),lwd=2,border="blue")
  
  polygon(x = c(5,7,7,5),y=c(min.expression-expression.step/2,
                             min.expression-expression.step/2,
                             min.expression-expression.step,
                             min.expression-expression.step),lwd=2,border="blue",col="blue")
  
  polygon(x = c(7,11,11,7),y=c(min.expression-expression.step/2,
                               min.expression-expression.step/2,
                               min.expression-expression.step,
                               min.expression-expression.step),lwd=2,border="blue")
  
  polygon(x = c(11,13,13,11),y=c(min.expression-expression.step/2,
                                 min.expression-expression.step/2,
                                 min.expression-expression.step,
                                 min.expression-expression.step),lwd=2,border="blue",col="blue")
  
  polygon(x = c(13,17,17,13),y=c(min.expression-expression.step/2,
                                 min.expression-expression.step/2,
                                 min.expression-expression.step,
                                 min.expression-expression.step),lwd=2,border="blue")
  
  polygon(x = c(17,19,19,17),y=c(min.expression-expression.step/2,
                                 min.expression-expression.step/2,
                                 min.expression-expression.step,
                                 min.expression-expression.step),lwd=2,border="blue",col="blue")
  
  polygon(x = c(19,23,23,19),y=c(min.expression-expression.step/2,
                                 min.expression-expression.step/2,
                                 min.expression-expression.step,
                                 min.expression-expression.step),lwd=2,border="blue")
  
  # polygon(x = c(23,25,25,23),y=c(min.expression-expression.step/2,
  #                                min.expression-expression.step/2,
  #                                min.expression-expression.step,
  #                                min.expression-expression.step),lwd=2,border="blue",col="lightblue")
  # 
  # polygon(x = c(25,29,29,25),y=c(min.expression-expression.step/2,
  #                                min.expression-expression.step/2,
  #                                min.expression-expression.step,
  #                                min.expression-expression.step),lwd=2,border="blue")
  # polygon(x = c(29,30,30,29),y=c(min.expression-expression.step/2,
  #                                min.expression-expression.step/2,
  #                                min.expression-expression.step,
  #                                min.expression-expression.step),lwd=2,border="blue",col="lightblue")
  
 # dev.off()
  return(0)  
  
  
  
  sd.zt <- paste("sd",paste0("zt",sprintf(fmt = "%02d",seq(from=0,to=20,by=4))),sep="_")
  current.gene.expression.sd.ll <- gene.expression[gene.id,c(paste(sd.zt,1,sep="_"),paste(sd.zt,2,sep="_"),paste(sd.zt,3,sep="_"),paste(sd.zt,4,sep="_"),paste(sd.zt,5,sep="_"))]
  
  min.expression <- min(current.gene.expression.sd.ll)
  max.expression <- max(current.gene.expression.sd.ll)
  range.expression <- max.expression - min.expression
  
  expression.step <- floor(range.expression / 5)
  
  plot(current.gene.expression.sd.ll,type="o",lwd=3,col="red",axes=F,xlab="",ylab="FPKM",
       ylim=c(min.expression-expression.step,max.expression),
       cex.lab=1.3,main=paste(gene.id, gene.name,sep=" - "),cex.main=2)
  axis(side=2,lwd=3)
  axis(side = 1,pos = min.expression - 1.1*expression.step, at = seq(from=1,to=30),
       labels = rep(paste("ZT",seq(from=0,to=20,by=4)),5),las=2,lwd=3)
  
  polygon(x = c(1,3,3,1),y=c(min.expression-expression.step/2,
                             min.expression-expression.step/2,
                             min.expression-expression.step,
                             min.expression-expression.step),lwd=2,border="red")
  
  polygon(x = c(3,7,7,3),y=c(min.expression-expression.step/2,
                             min.expression-expression.step/2,
                             min.expression-expression.step,
                             min.expression-expression.step),lwd=2,border="red",col="red")
  
  polygon(x = c(7,9,9,7),y=c(min.expression-expression.step/2,
                             min.expression-expression.step/2,
                             min.expression-expression.step,
                             min.expression-expression.step),lwd=2,border="red")
  
  polygon(x = c(9,13,13,9),y=c(min.expression-expression.step/2,
                               min.expression-expression.step/2,
                               min.expression-expression.step,
                               min.expression-expression.step),lwd=2,border="red",col="red")
  
  polygon(x = c(13,15,15,13),y=c(min.expression-expression.step/2,
                                 min.expression-expression.step/2,
                                 min.expression-expression.step,
                                 min.expression-expression.step),lwd=2,border="red")
  
  polygon(x = c(15,19,19,15),y=c(min.expression-expression.step/2,
                                 min.expression-expression.step/2,
                                 min.expression-expression.step,
                                 min.expression-expression.step),lwd=2,border="red",col="red")
  
  polygon(x = c(19,21,21,19),y=c(min.expression-expression.step/2,
                                 min.expression-expression.step/2,
                                 min.expression-expression.step,
                                 min.expression-expression.step),lwd=2,border="red")
  
  polygon(x = c(21,25,25,21),y=c(min.expression-expression.step/2,
                                 min.expression-expression.step/2,
                                 min.expression-expression.step,
                                 min.expression-expression.step),lwd=2,border="red",col="lightsalmon")
  
  polygon(x = c(25,27,27,25),y=c(min.expression-expression.step/2,
                                 min.expression-expression.step/2,
                                 min.expression-expression.step,
                                 min.expression-expression.step),lwd=2,border="red")
  polygon(x = c(27,30,30,27),y=c(min.expression-expression.step/2,
                                 min.expression-expression.step/2,
                                 min.expression-expression.step,
                                 min.expression-expression.step),lwd=2,border="red",col="lightsalmon")
  
}

plot.ld(gene.id = "ostta16g02620",gene.name = "ostta16g02620", gene.expression)
plot.ld(gene.id = "ostta08g00300",gene.name = "", gene.expression)
plot.ld(gene.id = "ostta05g01490",gene.name = "", gene.expression)
plot.ld(gene.id = "ostta11g01030",gene.name = "", gene.expression)
plot.ld(gene.id = "ostta18g01947",gene.name = "", gene.expression)
plot.ld(gene.id = "ostta01g05660",gene.name = "", gene.expression)
plot.ld(gene.id = "ostta11g00690",gene.name = "", gene.expression)
plot.ld(gene.id = "ostta06g02940",gene.name = "", gene.expression)
plot.ld(gene.id = "ostta05g00790",gene.name = "", gene.expression)
plot.ld(gene.id = "ostta06g03990",gene.name = "", gene.expression)
plot.ld(gene.id = "ostta07g04370",gene.name = "", gene.expression)
plot.ld(gene.id = "ostta06g00260",gene.name = "", gene.expression)
plot.ld(gene.id = "ostta02g04360",gene.name = "", gene.expression)
plot.ld(gene.id = "ostta03g00350",gene.name = "", gene.expression)

plot.ld(gene.id = "ostta06g02340",gene.name = "", gene.expression)


plot.sd <- function(gene.id, gene.name, gene.expression)#, mean.expression.ld, mean.expression.sd)
{
  
  sd.zt <- paste("sd",paste0("zt",sprintf(fmt = "%02d",seq(from=0,to=20,by=4))),sep="_")
  current.gene.expression.sd <- gene.expression[gene.id,c(paste(sd.zt,1,sep="_"),paste(sd.zt,2,sep="_"),paste(sd.zt,3,sep="_"))]
  
  min.expression <- min(current.gene.expression.sd)
  max.expression <- max(current.gene.expression.sd)
  range.expression <- max.expression - min.expression
  
  expression.step <- floor(range.expression / 5)

 # png(filename = paste(paste(c("sd_three_days",gene.name,gene.id),collapse="_"),".png",sep=""),width = 1000)
  plot(current.gene.expression.sd,type="o",lwd=5,col="red",axes=F,xlab="",ylab="FPKM",
       ylim=c(min.expression-expression.step,max.expression),
       cex.lab=1.3,main=gene.id,cex.main=2)
  axis(side=2,lwd=3)
  axis(side = 1,pos = min.expression - 1.1*expression.step, at = seq(from=1,to=30),
       labels = rep(paste("ZT",seq(from=0,to=20,by=4)),5),las=2,lwd=3)
  
  polygon(x = c(1,3,3,1),y=c(min.expression-expression.step/2,
                             min.expression-expression.step/2,
                             min.expression-expression.step,
                             min.expression-expression.step),lwd=2,border="red")
  
  polygon(x = c(3,7,7,3),y=c(min.expression-expression.step/2,
                             min.expression-expression.step/2,
                             min.expression-expression.step,
                             min.expression-expression.step),lwd=2,border="red",col="red")
  
  polygon(x = c(7,9,9,7),y=c(min.expression-expression.step/2,
                             min.expression-expression.step/2,
                             min.expression-expression.step,
                             min.expression-expression.step),lwd=2,border="red")
  
  polygon(x = c(9,13,13,9),y=c(min.expression-expression.step/2,
                               min.expression-expression.step/2,
                               min.expression-expression.step,
                               min.expression-expression.step),lwd=2,border="red",col="red")
  
  polygon(x = c(13,15,15,13),y=c(min.expression-expression.step/2,
                                 min.expression-expression.step/2,
                                 min.expression-expression.step,
                                 min.expression-expression.step),lwd=2,border="red")
  
  polygon(x = c(15,19,19,15),y=c(min.expression-expression.step/2,
                                 min.expression-expression.step/2,
                                 min.expression-expression.step,
                                 min.expression-expression.step),lwd=2,border="red",col="red")
#  dev.off()
  
}

plot.sd(gene.id = "ostta16g02620",gene.name = "", gene.expression)
plot.sd(gene.id = "ostta08g00300",gene.name = "", gene.expression)
plot.sd(gene.id = "ostta05g01490",gene.name = "", gene.expression)
plot.sd(gene.id = "ostta11g01030",gene.name = "", gene.expression)
plot.sd(gene.id = "ostta18g01947",gene.name = "", gene.expression)
plot.sd(gene.id = "ostta01g05660",gene.name = "", gene.expression)
plot.sd(gene.id = "ostta11g00690",gene.name = "", gene.expression)
plot.sd(gene.id = "ostta06g02940",gene.name = "", gene.expression)
plot.sd(gene.id = "ostta05g00790",gene.name = "", gene.expression)
plot.sd(gene.id = "ostta06g03990",gene.name = "", gene.expression)
plot.sd(gene.id = "ostta07g04370",gene.name = "", gene.expression)
plot.sd(gene.id = "ostta06g00260",gene.name = "", gene.expression)
plot.sd(gene.id = "ostta02g04360",gene.name = "", gene.expression)
plot.sd(gene.id = "ostta03g00350",gene.name = "", gene.expression)

plot.sd(gene.id = "ostta06g02340",gene.name = "", gene.expression)








res.sort <- sort(x = results.ld.ll$pVal,decreasing = F,index.return=T)

int.genes <- rownames(results.ld.ll)[res.sort$ix]

i <- 1
plot.ld.ll(gene.id = int.genes[i],gene.name = , gene.expression)
i <- i + 1

i <- i - 1

non.cycling.int.genes <- setdiff(complete.ld.rhythmic.genes,complete.ld.ll.rhythmic.genes)
length(int.genes)


i <- 1
plot.ld.ll(gene.id = non.cycling.int.genes[i],gene.name = "???", gene.expression)
i <- i + 1

## SPLS

sd.ld.genes <- intersect(complete.ld.rhythmic.genes,complete.sd.rhythmic.genes)
length(sd.ld.genes)

library(mixOmics)
library(ggplot2)
library(gridExtra)
library(grid)
library(lattice)

## Loading gene expression data that will be used as predictors (X)
tf.names <- read.table(file="transcription_factors_list.tsv",header=T,as.is=T)$gene_id
length(tf.names)

rhythmic.tfs <- intersect(tf.names,sd.ld.genes)
length(rhythmic.tfs)

ld.zt <- paste("ld",paste0("zt",sprintf(fmt = "%02d",seq(from=0,to=20,by=4))),sep="_")
sd.zt <- paste("sd",paste0("zt",sprintf(fmt = "%02d",seq(from=0,to=20,by=4))),sep="_")

gene.expression.ld.sd <- gene.expression[,c(paste(ld.zt,1,sep="_"),paste(ld.zt,2,sep="_"),paste(ld.zt,3,sep="_"),
                                            paste(sd.zt,1,sep="_"),paste(sd.zt,2,sep="_"),paste(sd.zt,3,sep="_"))]


colnames(gene.expression.ld.sd)


tfs.expression.sd.ld <- gene.expression.ld.sd[rhythmic.tfs,]
X <- t(tfs.expression.sd.ld) 
dim(X)


rhythmic.expression.ld.sd <- gene.expression.ld.sd[setdiff(sd.ld.genes,rhythmic.tfs),]

Y <- t(rhythmic.expression.ld.sd)
dim(Y)

myresult <- pls(X, Y, ncomp = 10, mode = "regression")
typeof(myresult)
attributes(myresult)

## We plot the explained variance by each component
tiff(filename = "explained_variance_10_component_replicates.tiff")
par(mfrow = c(1,2))
barplot(myresult$explained_variance$X, las = 2, main = "X",col=rainbow(5),ylim=c(0,0.7))
lines(x=c(0,12),y=c(0.1,0.1),col="black",lwd=4)
barplot(myresult$explained_variance$Y, las = 2, main = "Y",col=rainbow(5),ylim=c(0,0.7))
lines(x=c(0,12),y=c(0.1,0.1),col="black",lwd=4)
dev.off()

## Three components are selected and we apply PLS regression with ncomp=3
myresult <- pls(X, Y, ncomp = 3, mode = "regression", scale = TRUE) 

## Visualization of the individuals using the first two new components XY combined
tiff(filename = "individuals_visualization_on_the_new_components_joint_XY_with_replicates.tiff")
plotIndiv(object = myresult, comp = 1:2, rep.space = "XY-variate",
          ind.names = rownames(X),
          group = c(1,1,
                    2,2,
                    3,3,
                    4,4,
                    5,5,
                    6,6,
                    7,7,
                    8,8,
                    9,9,
                    10,10,
                    11,11,
                    12,12,
                    13,13,
                    14,14,
                    15,15,
                    16,16,
                    17,17,
                    18,18),
          # col = rep(rainbow(4),each=3), style = "ggplot2", ellipse = FALSE, ellipse.level = 0.95, centroid = FALSE,lwd=3,cex=6,xlim=c(-30,30)) 
          # col = rainbow(4), 
          col = rainbow(18),#c("black", "red","blue", "green"), 
          style = "ggplot2", ellipse = FALSE, ellipse.level = 0.95, centroid = FALSE,lwd=3,cex=6,xlim=c(-60,60)) 
dev.off()

## Individuals clustered according to their genotypes

tiff(filename = "individuals_visualization_on_the_new_components_separated_X_Y_with_replicates.tiff")
plotIndiv(object = myresult, comp = c(1,2),point.lwd = 4,cex = 5,
          ind.names = rownames(X),
          group = c(1,1,
                    2,2,
                    3,3,
                    4,4,
                    5,5,
                    6,6,
                    7,7,
                    8,8,
                    9,9,
                    10,10,
                    11,11,
                    12,12,
                    13,13,
                    14,14,
                    15,15,
                    16,16,
                    17,17,
                    18,18),
          # col = rep(rainbow(4),each=3), style = "ggplot2", ellipse = FALSE, 
          col = rainbow(18), style = "ggplot2", ellipse = FALSE, 
          ellipse.level = 0.95, centroid = FALSE,
          star = FALSE, legend = FALSE, abline = TRUE, alpha = 0)
dev.off()

## Individuals clustered according to their genotypes
tiff(filename = "correlation_circle_pls_with_replicates.tiff", 
     width = 5, height = 5, units = 'in', res = 200)
plotVar(object = myresult, comp = 1:2, cex = c(3, 0.2), col = c("forestgreen", "red3"))
dev.off()



## Cross validation of the model using LOO (Leave One Out)
myperfLoo = perf(myresult, validation = "loo", progressBar = TRUE)

myperfLoo$MSEP
myperfLoo$PRESS
myperfLoo$R2
myperfLoo$Q2
myperfLoo$Q2.total
myperfLoo$RSS

sum(myperfLoo$MSEP[,3])/length(sd.ld.genes)
## FC = 1.75 ---> 0.016

# tiff(filename = "MSEP_pls_with_replicates.tiff", 
#      width = 6, height = 6, units = "in", res = 150)
# plot(x = myperfLoo, 
#      criterion = "MSEP", 
#      xlab = "number of components",
#      ylab = NULL,
#      LimQ2 = 0.0975,
#      LimQ2.col = "none",type="l", lwd=3)
# dev.off()


## Applying Sparse PLS
mySPresult <- spls(X, Y, ncomp = 3, mode = 'regression', keepX = c(50, 50, 50), keepY = c(3000, 3000, 3000),
                   scale = TRUE)



mySPresult <- spls(X, Y, ncomp = 3, mode = 'regression', keepX = c(40, 40, 40), keepY = c(3000, 3000, 3000),
                   scale = TRUE)

mySPresult <- spls(X, Y, ncomp = 3, mode = 'regression', keepX = c(35, 35, 35), keepY = c(2800, 2800, 2800),
                   scale = TRUE)


## Visualization of the individuals using the first two new components XY combined
tiff(filename = "spls_individuals_visualization_on_the_new_components_joint_XY_with_replicates.tiff")
plotIndiv(object = mySPresult, comp = 1:2, rep.space = "XY-variate",
          ind.names = rownames(X),
          group = rep(1:18,each=2),
          # col = rep(rainbow(3),each=4), style = "ggplot2", ellipse = FALSE, ellipse.level = 0.95, centroid = FALSE,lwd=3,cex=6,xlim=c(-30,30)) 
          col = rainbow(18), style = "ggplot2", ellipse = FALSE, ellipse.level = 0.95, centroid = FALSE,lwd=3,cex=6,xlim=c(-30,30)) 
dev.off()


tiff(filename = "spls_correlation_circle_with_replicates.tiff", 
     width = 5, height = 5, units = 'in', res = 200)
plotVar(mySPresult, comp = 1:2, cex = c(4, 0.2), col = c("forestgreen", "red3"))
dev.off()

spls.loo <- perf(mySPresult, ncomp = 3, mode = 'regression', keepX = c(40,40,40), keepY = c(3000, 3000, 3000),validation = 'loo')      

spls.loo <- perf(mySPresult, ncomp = 3, mode = 'regression', keepX = c(35,35,35), keepY = c(2800, 2800, 2800),validation = 'loo')      


tiff(filename = "MSEP_spls_with_replicates.tiff")
palette(rainbow(3))
par(mfrow = c(1, 3))
for(i in 1:3){
  spls.rmsep <- sqrt(spls.loo$MSEP[i, ])
  pls.rmsep <- sqrt(myperfLoo$MSEP[i, ])
  matplot(cbind(spls.rmsep, pls.rmsep), type = 'l', col = i, ylab = 'MSEP', lwd = 2,
          xlab = 'dim', lty = c(1, 2), axes = FALSE)
  axis(1, 1:3, labels = 1:3)
  axis(2)
  title(main = paste(rownames(spls.loo$MSEP)[i]))
}
palette("default")
dev.off()

## Bipartite representation
color.edge <- colorRampPalette(c("red4", "white", "darkgreen"))
spl.th <- 0.75
spl.th <- 0.72

spl.th <- 0.7

spl.th <- 0.6
spl.th <- 0.5
tiff(filename = "bipartite_graph_with_replicates.tiff")#, height = 6, width = 6,
#units = "in", res = 150)
res <- network(mySPresult, comp = 1:2, cutoff = spl.th, shape.node = c("rectangle", "rectangle"),cex.node.name = 1,
               color.node = c("white", "coral1"), color.edge = color.edge(10),save = "png",name.save = "network2.png")
dev.off()

attributes(res)
res$gR
"ostta01g05650" %in% colnames(res$M)
"ostta13g01820" %in% rownames(res$M)

network.matrix <- res$M

max(network.matrix["ostta13g01820",])
network.matrix[,"ostta01g05650"]



library(igraph)

write.graph(graph = res$gR,file = "fisrt_try_ostta_network_toc1_5.gml",format = "gml")  

ostta.network <- read.graph(file = "fisrt_try_ostta_network.gml",format = "gml")

sort(degree(ostta.network),decreasing = T)[1:10]

ostta.graph <- as_edgelist(graph = ostta.network,names = T)

install.packages("linkcomm")
library(linkcomm)
lc <- getLinkCommunities(network = ostta.graph, hcmethod = "single")

print(lc)
save(lc,file = "link_communities_directed.rda")
plot(lc, type = "graph", layout = layout.fruchterman.reingold)
plot(lc, type = "graph", layout = "spencer.circle")
plot(lc, type = "graph", layout = "spencer.circle", shownodesin = 3)
plot(lc, type = "members")
plot(lc, type = "summary")
plot(lc, type = "dend")



cluster.files <- list.files(path = "LDpeak_trough_clusters/")
gene.cluster <- read.table(file="LDpeak_trough_clusters/cluster_peak_ZT0_trough_ZT12.txt",header = F,as.is = T)[[1]]
gene.cluster <- read.table(file="LDpeak_trough_clusters/cluster_peak_ZT0_trough_ZT16.txt",header = F,as.is = T)[[1]]
gene.cluster <- read.table(file="LDpeak_trough_clusters/cluster_peak_ZT0_trough_ZT20.txt",header = F,as.is = T)[[1]]

gene.cluster <- read.table(file="LDpeak_trough_clusters/cluster_peak_ZT12_trough_ZT0.txt",header = F,as.is = T)[[1]]


i <- 1
gene.cluster <- read.table(file=paste0("LDpeak_trough_clusters/",cluster.files[i]),header = F,as.is = T)[[1]]
length(gene.cluster)
i <- i + 1


ld.zt <- paste("ld",paste0("zt",sprintf(fmt = "%02d",seq(from=0,to=20,by=4))),sep="_")

cluster.ld.ll.expression <- gene.expression[gene.cluster,c(paste(ld.zt,1,sep="_"),paste(ld.zt,2,sep="_"),paste(ld.zt,3,sep="_"),paste(ld.zt,4,sep="_"),paste(ld.zt,5,sep="_"))]

cluster.ld.ll.scaled.expression <- apply(X = cluster.ld.ll.expression,MARGIN = 1,FUN = scale)
current.gene.expression.ld.ll <- 5*rowMeans(cluster.ld.ll.scaled.expression)

min.expression <- min(current.gene.expression.ld.ll)
max.expression <- max(current.gene.expression.ld.ll)
range.expression <- max.expression - min.expression

expression.step <- floor(range.expression / 5)

plot(current.gene.expression.ld.ll,type="o",lwd=6,col="blue",axes=F,xlab="",ylab="FPKM",
     ylim=c(min.expression-expression.step,max.expression),
     cex.lab=1.3,main=paste(gene.id, gene.name,sep=" - "),cex.main=2)

#axis(side=2,lwd=3,labels="")
axis(side = 1,pos = min.expression - 1.1*expression.step, at = seq(from=1,to=30),
     labels = rep(paste("ZT",seq(from=0,to=20,by=4)),5),las=2,lwd=3)

polygon(x = c(1,5,5,1),y=c(min.expression-expression.step/2,
                           min.expression-expression.step/2,
                           min.expression-expression.step,
                           min.expression-expression.step),lwd=2,border="blue")

polygon(x = c(5,7,7,5),y=c(min.expression-expression.step/2,
                           min.expression-expression.step/2,
                           min.expression-expression.step,
                           min.expression-expression.step),lwd=2,border="blue",col="blue")

polygon(x = c(7,11,11,7),y=c(min.expression-expression.step/2,
                             min.expression-expression.step/2,
                             min.expression-expression.step,
                             min.expression-expression.step),lwd=2,border="blue")

polygon(x = c(11,13,13,11),y=c(min.expression-expression.step/2,
                               min.expression-expression.step/2,
                               min.expression-expression.step,
                               min.expression-expression.step),lwd=2,border="blue",col="blue")

polygon(x = c(13,17,17,13),y=c(min.expression-expression.step/2,
                               min.expression-expression.step/2,
                               min.expression-expression.step,
                               min.expression-expression.step),lwd=2,border="blue")

polygon(x = c(17,19,19,17),y=c(min.expression-expression.step/2,
                               min.expression-expression.step/2,
                               min.expression-expression.step,
                               min.expression-expression.step),lwd=2,border="blue",col="blue")

polygon(x = c(19,23,23,19),y=c(min.expression-expression.step/2,
                               min.expression-expression.step/2,
                               min.expression-expression.step,
                               min.expression-expression.step),lwd=2,border="blue")

polygon(x = c(23,25,25,23),y=c(min.expression-expression.step/2,
                               min.expression-expression.step/2,
                               min.expression-expression.step,
                               min.expression-expression.step),lwd=2,border="blue",col="lightblue")

polygon(x = c(25,29,29,25),y=c(min.expression-expression.step/2,
                               min.expression-expression.step/2,
                               min.expression-expression.step,
                               min.expression-expression.step),lwd=2,border="blue")
polygon(x = c(29,30,30,29),y=c(min.expression-expression.step/2,
                               min.expression-expression.step/2,
                               min.expression-expression.step,
                               min.expression-expression.step),lwd=2,border="blue",col="lightblue")

i <- i - 1



sd.zt <- paste("sd",paste0("zt",sprintf(fmt = "%02d",seq(from=0,to=20,by=4))),sep="_")

cluster.sd.ll.expression <- gene.expression[gene.cluster,c(paste(sd.zt,1,sep="_"),paste(sd.zt,2,sep="_"),paste(sd.zt,3,sep="_"),paste(sd.zt,4,sep="_"),paste(sd.zt,5,sep="_"))]

cluster.sd.ll.scaled.expression <- apply(X = cluster.sd.ll.expression,MARGIN = 1,FUN = scale)
current.gene.expression.sd.ll <- 5*rowMeans(cluster.sd.ll.scaled.expression)

min.expression <- min(current.gene.expression.sd.ll)
max.expression <- max(current.gene.expression.sd.ll)
range.expression <- max.expression - min.expression

expression.step <- floor(range.expression / 5)

plot(current.gene.expression.sd.ll,type="o",lwd=6,col="red",axes=F,xlab="",ylab="FPKM",
     ylim=c(min.expression-expression.step,max.expression),
     cex.lab=1.3,main=paste(gene.id, gene.name,sep=" - "),cex.main=2)
#axis(side=2,lwd=3)
axis(side = 1,pos = min.expression - 1.1*expression.step, at = seq(from=1,to=30),
     labels = rep(paste("ZT",seq(from=0,to=20,by=4)),5),las=2,lwd=3)

polygon(x = c(1,3,3,1),y=c(min.expression-expression.step/2,
                           min.expression-expression.step/2,
                           min.expression-expression.step,
                           min.expression-expression.step),lwd=2,border="red")

polygon(x = c(3,7,7,3),y=c(min.expression-expression.step/2,
                           min.expression-expression.step/2,
                           min.expression-expression.step,
                           min.expression-expression.step),lwd=2,border="red",col="red")

polygon(x = c(7,9,9,7),y=c(min.expression-expression.step/2,
                           min.expression-expression.step/2,
                           min.expression-expression.step,
                           min.expression-expression.step),lwd=2,border="red")

polygon(x = c(9,13,13,9),y=c(min.expression-expression.step/2,
                             min.expression-expression.step/2,
                             min.expression-expression.step,
                             min.expression-expression.step),lwd=2,border="red",col="red")

polygon(x = c(13,15,15,13),y=c(min.expression-expression.step/2,
                               min.expression-expression.step/2,
                               min.expression-expression.step,
                               min.expression-expression.step),lwd=2,border="red")

polygon(x = c(15,19,19,15),y=c(min.expression-expression.step/2,
                               min.expression-expression.step/2,
                               min.expression-expression.step,
                               min.expression-expression.step),lwd=2,border="red",col="red")

polygon(x = c(19,21,21,19),y=c(min.expression-expression.step/2,
                               min.expression-expression.step/2,
                               min.expression-expression.step,
                               min.expression-expression.step),lwd=2,border="red")

polygon(x = c(21,25,25,21),y=c(min.expression-expression.step/2,
                               min.expression-expression.step/2,
                               min.expression-expression.step,
                               min.expression-expression.step),lwd=2,border="red",col="lightsalmon")

polygon(x = c(25,27,27,25),y=c(min.expression-expression.step/2,
                               min.expression-expression.step/2,
                               min.expression-expression.step,
                               min.expression-expression.step),lwd=2,border="red")
polygon(x = c(27,30,30,27),y=c(min.expression-expression.step/2,
                               min.expression-expression.step/2,
                               min.expression-expression.step,
                               min.expression-expression.step),lwd=2,border="red",col="lightsalmon")
                               
                               
                               

par(lwd=2)
boxplot(apply(X = gene.expression[intersect(complete.ld.rhythmic.genes,complete.sd.rhythmic.genes),],MARGIN = 1,FUN = max),apply(X = gene.expression[intersect(non.rhythmic.ld.genes,non.rhythmic.sd.genes),], MARGIN = 1,FUN = max),outline=F,col=colors()[c(89,12)],names=c("",""),ylab="",cex.lab=1.5,las=2,cex.axis=1.5)
dev.off()

par(lwd=2)
boxplot(apply(X = gene.expression[intersect(complete.ld.rhythmic.genes,complete.sd.rhythmic.genes),c(sd.zt.i,ld.zt.i)],MARGIN = 1,FUN = max),apply(X = gene.expression[intersect(non.rhythmic.ld.genes,non.rhythmic.sd.genes),c(sd.zt.i,ld.zt.i)], MARGIN = 1,FUN = max),outline=F,col=colors()[c(89,12)],names=c("",""),ylab="",cex.lab=1.5,las=2,cex.axis=1.5)
dev.off()

wilcox.test(apply(X = gene.expression[intersect(complete.ld.rhythmic.genes,complete.sd.rhythmic.genes),],MARGIN = 1,FUN = max),1.75*apply(X = gene.expression[unique(c(non.rhythmic.ld.genes,non.rhythmic.sd.genes)),], MARGIN = 1,FUN = max),alternative="greater")



boxplot(apply(X = gene.expression[intersect(complete.ld.rhythmic.genes,complete.sd.rhythmic.genes),],MARGIN = 1,FUN = max),apply(X = gene.expression[unique(c(non.rhythmic.ld.genes,non.rhythmic.sd.genes)),], MARGIN = 1,FUN = max),outline=F)


wilcox.test(apply(X = gene.expression[intersect(complete.ld.rhythmic.genes,complete.sd.rhythmic.genes),],MARGIN = 1,FUN = max),1.75*apply(X = gene.expression[unique(c(non.rhythmic.ld.genes,non.rhythmic.sd.genes)),], MARGIN = 1,FUN = max),alternative="greater")

wilcox.test(apply(X = gene.expression[intersect(complete.ld.rhythmic.genes,complete.sd.rhythmic.genes),c(sd.zt.i,ld.zt.i)],MARGIN = 1,FUN = max),3*apply(X = gene.expression[intersect(non.rhythmic.ld.genes,non.rhythmic.sd.genes),c(sd.zt.i,ld.zt.i)], MARGIN = 1,FUN = max),alternative="greater")


data.ld.ll <- matrix(c(length(complete.ld.ll.rhythmic.genes),length(setdiff(complete.ld.rhythmic.genes,complete.ld.ll.rhythmic.genes)),length(non.rhythmic.ld.genes)),nrow=3)
sum(c(length(complete.ld.ll.rhythmic.genes),length(setdiff(complete.ld.rhythmic.genes,complete.ld.ll.rhythmic.genes)),length(non.rhythmic.ld.genes)))

par(lwd=3)
barplot(data.ld.ll, lwd=2, names.arg = "",cex.names = 2,cex.axis = 1.5,
        col=c("blue","lightblue","white"), 
        border="black", 
        font.axis=2, 
        xlab="",ylim=c(0,8000))
dev.off()


data.sd.ll <- matrix(c(length(complete.sd.ll.rhythmic.genes),length(setdiff(complete.sd.rhythmic.genes,complete.sd.ll.rhythmic.genes)),length(non.rhythmic.sd.genes)),nrow=3)
sum(c(length(complete.sd.ll.rhythmic.genes),length(setdiff(complete.sd.rhythmic.genes,complete.sd.ll.rhythmic.genes)),length(non.rhythmic.sd.genes)))

par(lwd=3)
barplot(data.sd.ll, lwd=2, names.arg = "",cex.names = 2,cex.axis = 1.5,
        col=c("red","lightpink","white"), 
        border="black", 
        font.axis=2, 
        xlab="",ylim=c(0,8000))
dev.off()



plot(x = -10,y= -10,axes=F,xlab="",ylab="",
     ylim=c(min.expr-2, max.expr),xlim=c(0,72),
     main=paste(gene.id, gene.name,sep=" - "),cex.main=2)

lines(x = seq(from=0,by=4,to=68),current.gene.expression.ld,type="o",lwd=5,col="blue")

if(!scale)
{
  axis(side=2,lwd=3)
}

###------------------

gene.set <- setdiff(complete.ld.rhythmic.genes,complete.ld.ll.rhythmic.genes)

i <- 1
gene.id <- complete.ld.ll.rhythmic.genes[i] 
gene.id <- gene.set[i]
ld.zt <- paste("ld",paste0("zt",sprintf(fmt = "%02d",seq(from=0,to=20,by=4))),sep="_")
current.gene.expression.ld.ll <- gene.expression[gene.id,c(paste(ld.zt,1,sep="_"),paste(ld.zt,2,sep="_"),paste(ld.zt,3,sep="_"),paste(ld.zt,4,sep="_"),paste(ld.zt,5,sep="_"))]

min.expression <- min(current.gene.expression.ld.ll)
max.expression <- max(current.gene.expression.ld.ll)
range.expression <- max.expression - min.expression

expression.step <- floor(range.expression / 5)

plot(current.gene.expression.ld.ll,type="o",lwd=8,col="blue",axes=F,xlab="",ylab="",
     ylim=c(min.expression-expression.step,max.expression),
     cex.lab=1.3,main=gene.id,cex.main=2)
axis(side=2,lwd=4,cex.axis=1.5,las=1)
axis(side = 1,pos = min.expression - 1.1*expression.step, at = seq(from=1,to=30),
     labels = rep("",30),las=2,lwd=3,cex.axis=1.5)
     #labels = rep(paste("ZT",seq(from=0,to=20,by=4)),5),las=2,lwd=3,cex.axis=1.5)

polygon(x = c(1,5,5,1),y=c(min.expression-expression.step/2,
                           min.expression-expression.step/2,
                           min.expression-expression.step,
                           min.expression-expression.step),lwd=4,border="blue")

polygon(x = c(5,7,7,5),y=c(min.expression-expression.step/2,
                           min.expression-expression.step/2,
                           min.expression-expression.step,
                           min.expression-expression.step),lwd=4,border="blue",col="blue")

polygon(x = c(7,11,11,7),y=c(min.expression-expression.step/2,
                             min.expression-expression.step/2,
                             min.expression-expression.step,
                             min.expression-expression.step),lwd=4,border="blue")

polygon(x = c(11,13,13,11),y=c(min.expression-expression.step/2,
                               min.expression-expression.step/2,
                               min.expression-expression.step,
                               min.expression-expression.step),lwd=4,border="blue",col="blue")

polygon(x = c(13,17,17,13),y=c(min.expression-expression.step/2,
                               min.expression-expression.step/2,
                               min.expression-expression.step,
                               min.expression-expression.step),lwd=4,border="blue")

polygon(x = c(17,19,19,17),y=c(min.expression-expression.step/2,
                               min.expression-expression.step/2,
                               min.expression-expression.step,
                               min.expression-expression.step),lwd=4,border="blue",col="blue")

polygon(x = c(19,23,23,19),y=c(min.expression-expression.step/2,
                               min.expression-expression.step/2,
                               min.expression-expression.step,
                               min.expression-expression.step),lwd=4,border="blue")

polygon(x = c(23,25,25,23),y=c(min.expression-expression.step/2,
                               min.expression-expression.step/2,
                               min.expression-expression.step,
                               min.expression-expression.step),lwd=4,border="blue",col="lightblue")

polygon(x = c(25,29,29,25),y=c(min.expression-expression.step/2,
                               min.expression-expression.step/2,
                               min.expression-expression.step,
                               min.expression-expression.step),lwd=4,border="blue")
polygon(x = c(29,30,30,29),y=c(min.expression-expression.step/2,
                               min.expression-expression.step/2,
                               min.expression-expression.step,
                               min.expression-expression.step),lwd=4,border="blue",col="lightblue")

i <- i + 1
##---------------
i <- 1

gene.set <- setdiff(complete.sd.rhythmic.genes,complete.sd.ll.rhythmic.genes)


gene.id <- complete.sd.ll.rhythmic.genes[i]
gene.id <- gene.set[i]
current.gene.expression.sd.ll <- gene.expression[gene.id,31:60]

min.expression <- min(current.gene.expression.sd.ll)
max.expression <- max(current.gene.expression.sd.ll)
range.expression <- max.expression - min.expression

expression.step <- floor(range.expression / 5)

plot(current.gene.expression.sd.ll,type="o",lwd=8,col="red",axes=F,xlab="",ylab="",
     ylim=c(min.expression-expression.step,max.expression),
     cex.lab=1.3,main=gene.id,cex.main=2)
axis(side=2,lwd=4,cex.axis=1.5,las=1)
axis(side = 1,pos = min.expression - 1.1*expression.step, at = seq(from=1,to=30),
     labels = rep("",30),las=2,lwd=3)

polygon(x = c(1,3,3,1),y=c(min.expression-expression.step/2,
                           min.expression-expression.step/2,
                           min.expression-expression.step,
                           min.expression-expression.step),lwd=4,border="red")

polygon(x = c(3,7,7,3),y=c(min.expression-expression.step/2,
                           min.expression-expression.step/2,
                           min.expression-expression.step,
                           min.expression-expression.step),lwd=4,border="red",col="red")

polygon(x = c(7,9,9,7),y=c(min.expression-expression.step/2,
                           min.expression-expression.step/2,
                           min.expression-expression.step,
                           min.expression-expression.step),lwd=4,border="red")

polygon(x = c(9,13,13,9),y=c(min.expression-expression.step/2,
                             min.expression-expression.step/2,
                             min.expression-expression.step,
                             min.expression-expression.step),lwd=4,border="red",col="red")

polygon(x = c(13,15,15,13),y=c(min.expression-expression.step/2,
                               min.expression-expression.step/2,
                               min.expression-expression.step,
                               min.expression-expression.step),lwd=4,border="red")

polygon(x = c(15,19,19,15),y=c(min.expression-expression.step/2,
                               min.expression-expression.step/2,
                               min.expression-expression.step,
                               min.expression-expression.step),lwd=4,border="red",col="red")

polygon(x = c(19,21,21,19),y=c(min.expression-expression.step/2,
                               min.expression-expression.step/2,
                               min.expression-expression.step,
                               min.expression-expression.step),lwd=4,border="red")

polygon(x = c(21,25,25,21),y=c(min.expression-expression.step/2,
                               min.expression-expression.step/2,
                               min.expression-expression.step,
                               min.expression-expression.step),lwd=4,border="red",col="lightsalmon")

polygon(x = c(25,27,27,25),y=c(min.expression-expression.step/2,
                               min.expression-expression.step/2,
                               min.expression-expression.step,
                               min.expression-expression.step),lwd=4,border="red")
polygon(x = c(27,30,30,27),y=c(min.expression-expression.step/2,
                               min.expression-expression.step/2,
                               min.expression-expression.step,
                               min.expression-expression.step),lwd=4,border="red",col="lightsalmon")

i <- i + 1





##---------------

library(circacompare)

?circacompare

dim(ld.ll.gene.expression)

ld.zt <- paste("ld",paste0("zt",sprintf(fmt = "%02d",seq(from=0,to=20,by=4))),sep="_")

circacompare.ld.ll <- matrix(nrow=length(complete.ld.ll.rhythmic.genes),ncol=15)
rownames(circacompare.ld.ll) <- complete.ld.ll.rhythmic.genes

for(i in 1:length(complete.ld.ll.rhythmic.genes))
{
  gene.i <- complete.ld.ll.rhythmic.genes[i]
  
  ld.expression.i <- ld.ll.gene.expression[gene.i,c(paste(ld.zt,2,sep="_"),paste(ld.zt,3,sep="_"))]
  ll.expression.i <- ld.ll.gene.expression[gene.i,c(paste(ld.zt,4,sep="_"),paste(ld.zt,5,sep="_"))]
  
  time.points <- seq(from=0,by=4,length.out = 12)
  
  ld.ll.df <- data.frame(time=c(time.points,time.points),
                         measure=c(ld.expression.i, ll.expression.i),
                         group=c(rep("ld",12),rep("ll",12)))
  
  out.i <- circacompare(x = ld.ll.df, col_time = "time", col_group = "group", col_outcome = "measure",alpha_threshold = 1)
  circacompare.ld.ll[i,] <- out.i[[2]][,2]
 }


par(lwd=3)
boxplot(circacompare.ld.ll[,"ld amplitude estimate" ],circacompare.ld.ll[,"ll amplitude estimate" ],outline=F,col=c("blue","lightblue"),names=c("LD","LL"),cex.axis=2)
boxplot(circacompare.ld.ll[,"Amplitude difference estimate" ],outline=F,ylim=c(-90,80))
sum(circacompare.ld.ll[,"P-value for amplitude difference" ] < 0.05)

boxplot(circacompare.ld.ll[,"ld mesor estimate" ],circacompare.ld.ll[,"ll mesor estimate" ],outline=F,col=c("blue","lightblue"),names=c("LD","LL"),cex.axis=2)
boxplot(circacompare.ld.ll[,"Mesor difference estimate" ],outline=F,ylim=c(-90,80))
sum(circacompare.ld.ll[,"P-value for mesor difference" ] < 0.05)

sum(circacompare.ld.ll[,"P-value for difference in phase" ] < 0.05)


sd.zt <- paste("sd",paste0("zt",sprintf(fmt = "%02d",seq(from=0,to=20,by=4))),sep="_")

circacompare.sd.ll <- matrix(nrow=length(complete.sd.ll.rhythmic.genes),ncol=15)
rownames(circacompare.sd.ll) <- complete.sd.ll.rhythmic.genes
i <- 1
for(i in 1:length(complete.sd.ll.rhythmic.genes))
{
  gene.i <- complete.sd.ll.rhythmic.genes[i]
  
  sd.expression.i <- sd.ll.gene.expression[gene.i,c(paste(sd.zt,2,sep="_"),paste(sd.zt,3,sep="_"))]
  ll.expression.i <- sd.ll.gene.expression[gene.i,c(paste(sd.zt,4,sep="_"),paste(sd.zt,5,sep="_"))]
  
  time.points <- seq(from=0,by=4,length.out = 12)
  
  sd.ll.df <- data.frame(time=c(time.points,time.points),
                         measure=c(sd.expression.i, ll.expression.i),
                         group=c(rep("sd",12),rep("ll",12)))
  
  out.i <- circacompare(x = sd.ll.df, col_time = "time", col_group = "group", col_outcome = "measure",alpha_threshold = 1)
  out.i[[1]]
  plot.sd.ll(gene.id = gene.i,gene.name = "",gene.expression = gene.expression)
  out.i[[2]]
  circacompare.sd.ll[i,] <- out.i[[2]][,2]
}

colnames(circacompare.sd.ll) <- out.i[[2]][,1]

par(lwd=3)
boxplot(circacompare.sd.ll[,"sd amplitude estimate" ],circacompare.sd.ll[,"ll amplitude estimate" ],outline=F,col=c("red","lightsalmon"),names=c("SD","LL"),cex.axis=2)
boxplot(circacompare.sd.ll[,"Amplitude difference estimate" ],outline=F,ylim=c(-90,80))
sum(circacompare.sd.ll[,"P-value for amplitude difference" ] < 0.05)

boxplot(circacompare.sd.ll[,"sd mesor estimate" ],circacompare.sd.ll[,"ll mesor estimate" ],outline=F,col=c("red","lightsalmon"),names=c("SD","LL"),cex.axis=2)
boxplot(circacompare.ld.ll[,"Mesor difference estimate" ],outline=F,ylim=c(-90,80))
sum(circacompare.ld.ll[,"P-value for mesor difference" ] < 0.05)

sum(circacompare.ld.ll[,"P-value for difference in phase" ] < 0.05)


grid.newpage()
draw.pairwise.venn(area1 = length(intersect(complete.ld.sd.rhythmic.genes,complete.ld.ll.rhythmic.genes)),area2 = length(intersect(complete.ld.sd.rhythmic.genes,complete.sd.ll.rhythmic.genes)),cross.area = length(intersect(complete.ld.sd.rhythmic.genes,intersect(complete.ld.ll.rhythmic.genes,complete.sd.ll.rhythmic.genes))),lwd = 3,category = c("LD","SD"),euler.d = T,col = c("blue","red"),fill = c(" blue","red"),alpha = 0.3,cex = 2,cat.cex = 2)


grid.newpage()
draw.pairwise.venn(area1 = length(complete.ld.ll.rhythmic.genes),area2 = length(complete.sd.ll.rhythmic.genes),cross.area = length(intersect(complete.ld.ll.rhythmic.genes,complete.sd.ll.rhythmic.genes)),lwd = 3,category = c("LD","SD"),euler.d = T,col = c("blue","red"),fill = c(" blue","red"),alpha = 0.3,cex = 2,cat.cex = 2)



wave.form <- function(mesor, amplitude,period,phase,time=seq(from=0,to=48,by=0.01))
{
  y <- mesor + amplitude*cos(0.0174533*period*(time - phase))
  return(y)
}


## effect on synchronization
plot(x=0,y=0,col="white",ylim=c(20,80),xlim=c(0,120),xlab="",ylab="",axes=F)
time <- seq(from=0,to=72,by=0.01)
time.2 <- seq(from=72,to=120,by=0.01)
N <- 50
norm.random.1 <- rnorm(n = N,mean = 0,sd = 1)
norm.random.2 <- rnorm(n = N,mean = 0,sd = 6)
syn.waves <- matrix(data = 0,nrow = N,ncol = length(time))
asyn.waves <- matrix(data = 0,nrow = N,ncol = length(time.2))

for(i in 1:N)
{
  syn.waves[i,] <- wave.form(mesor = 50,amplitude = 20,period = 15,phase = 16+norm.random.1[i],time=time)
  asyn.waves[i,] <- wave.form(mesor = 50,amplitude = 20,period = 15,phase = 16+norm.random.2[i],time=time.2)
  lines(time,syn.waves[i,],type="l",lwd=1,col="azure2")
  lines(time.2,asyn.waves[i,],type="l",lwd=1,col="azure2")
  
}

lines(x=c(time,time.2),y=c(colMeans(syn.waves),colMeans(asyn.waves)),type="l",lwd=5,col="blue")

polygon(x = c(0,16,16,0),y=c(25,25,22,22),lwd=4,border="blue")
polygon(x = c(16,24,24,16),y=c(25,25,22,22),lwd=4,border="blue",col="blue")
polygon(x = c(24,40,40,24),y=c(25,25,22,22),lwd=4,border="blue")
polygon(x = c(40,48,48,40),y=c(25,25,22,22),lwd=4,border="blue",col="blue")
polygon(x = c(48,64,64,48),y=c(25,25,22,22),lwd=4,border="blue")
polygon(x = c(64,72,72,64),y=c(25,25,22,22),lwd=4,border="blue",col="blue")
polygon(x = c(72,88,88,72),y=c(25,25,22,22),lwd=4,border="blue")
polygon(x = c(88,96,96,88),y=c(25,25,22,22),lwd=4,border="blue",col="lightblue")
polygon(x = c(96,112,112,96),y=c(25,25,22,22),lwd=4,border="blue")
polygon(x = c(112,120,120,112),y=c(25,25,22,22),lwd=4,border="blue",col="lightblue")


## effect on amplitude
plot(x=0,y=0,col="white",ylim=c(20,80),xlim=c(0,120),xlab="",ylab="",axes=F)
time <- seq(from=0,to=72,by=0.01)
time.2 <- seq(from=72,to=120,by=0.01)
N <- 50
norm.random.1 <- - abs(rnorm(n = N,mean = 0,sd = 5))
norm.random.2 <- -abs(rnorm(n = N,mean = 0,sd = 15))
syn.waves <- matrix(data = 0,nrow = N,ncol = length(time))
asyn.waves <- matrix(data = 0,nrow = N,ncol = length(time.2))
for(i in 1:N)
{
  syn.waves[i,] <- wave.form(mesor = 50,amplitude = 20+norm.random.1[i],period = 15,phase = 16,time=time)
  asyn.waves[i,] <- wave.form(mesor = 50,amplitude = abs(15+norm.random.2[i]),period = 15,phase = 16,time=time.2)
  lines(time,syn.waves[i,],type="l",lwd=1,col="azure2")
  lines(time.2,asyn.waves[i,],type="l",lwd=1,col="azure2")
  
}

lines(x=c(time,time.2),y=c(colMeans(syn.waves),colMeans(asyn.waves)),type="l",lwd=5,col="blue")
polygon(x = c(0,16,16,0),y=c(25,25,22,22),lwd=4,border="blue")
polygon(x = c(16,24,24,16),y=c(25,25,22,22),lwd=4,border="blue",col="blue")
polygon(x = c(24,40,40,24),y=c(25,25,22,22),lwd=4,border="blue")
polygon(x = c(40,48,48,40),y=c(25,25,22,22),lwd=4,border="blue",col="blue")
polygon(x = c(48,64,64,48),y=c(25,25,22,22),lwd=4,border="blue")
polygon(x = c(64,72,72,64),y=c(25,25,22,22),lwd=4,border="blue",col="blue")
polygon(x = c(72,88,88,72),y=c(25,25,22,22),lwd=4,border="blue")
polygon(x = c(88,96,96,88),y=c(25,25,22,22),lwd=4,border="blue",col="lightblue")
polygon(x = c(96,112,112,96),y=c(25,25,22,22),lwd=4,border="blue")
polygon(x = c(112,120,120,112),y=c(25,25,22,22),lwd=4,border="blue",col="lightblue")


## effect on synchronization + amplitude
plot(x=0,y=0,col="white",ylim=c(20,80),xlim=c(0,120),xlab="",ylab="",axes=F)
time <- seq(from=0,to=72,by=0.01)
time.2 <- seq(from=72,to=120,by=0.01)
N <- 50
norm.random.1 <- - abs(rnorm(n = N,mean = 0,sd = 1))
norm.random.2 <- -abs(rnorm(n = N,mean = 0,sd = 8))
norm.random.3 <- rnorm(n = N,mean = 0,sd = 1)
norm.random.4 <- rnorm(n = N,mean = 0,sd = 4)

syn.waves <- matrix(data = 0,nrow = N,ncol = length(time))
asyn.waves <- matrix(data = 0,nrow = N,ncol = length(time.2))
for(i in 1:N)
{
  syn.waves[i,] <- wave.form(mesor = 50,amplitude = 20+norm.random.1[i],period = 15+norm.random.3[i],phase = 16,time=time)
  asyn.waves[i,] <- wave.form(mesor = 50,amplitude = abs(20+norm.random.2[i]),period = 15+norm.random.4[i],phase = 16,time=time.2)
  lines(time,syn.waves[i,],type="l",lwd=1,col="azure2")
  lines(time.2,asyn.waves[i,],type="l",lwd=1,col="azure2")
  
}

lines(x=c(time,time.2),y=c(colMeans(syn.waves),colMeans(asyn.waves)),type="l",lwd=5,col="blue")
polygon(x = c(0,16,16,0),y=c(25,25,22,22),lwd=4,border="blue")
polygon(x = c(16,24,24,16),y=c(25,25,22,22),lwd=4,border="blue",col="blue")
polygon(x = c(24,40,40,24),y=c(25,25,22,22),lwd=4,border="blue")
polygon(x = c(40,48,48,40),y=c(25,25,22,22),lwd=4,border="blue",col="blue")
polygon(x = c(48,64,64,48),y=c(25,25,22,22),lwd=4,border="blue")
polygon(x = c(64,72,72,64),y=c(25,25,22,22),lwd=4,border="blue",col="blue")
polygon(x = c(72,88,88,72),y=c(25,25,22,22),lwd=4,border="blue")
polygon(x = c(88,96,96,88),y=c(25,25,22,22),lwd=4,border="blue",col="lightblue")
polygon(x = c(96,112,112,96),y=c(25,25,22,22),lwd=4,border="blue")
polygon(x = c(112,120,120,112),y=c(25,25,22,22),lwd=4,border="blue",col="lightblue")








## Effect over amplitude
N <- 50
time <- seq(from=0,to=72,by=0.01)
norm.random <- rnorm(n = N,mean = 0,sd = 1)
syn.waves <- matrix(data = 0,nrow = N,ncol = length(time))
plot(time,wave.form(mesor = 50,amplitude = 20,period = 15,phase = 16,time=time),type="l",lwd=4,col="white",axes=F,xlab="",ylab="")
#colors <- rainbow(N)
for(i in 1:N)
{
  syn.waves[i,] <- wave.form(mesor = 50,amplitude = 20,period = 15,phase = 16+norm.random[i],time=time)
  lines(time,wave.form(mesor = 50,amplitude = 20,period = 15,phase = 16+norm.random[i],time=time),type="l",lwd=1,col="azure2")
}

lines(time,colMeans(syn.waves),type="l",lwd=5,col="blue")


N <- 50
time <- seq(from=0,to=48,by=0.01)
norm.random <- -abs(rnorm(n = N,mean = 0,sd = 10))
asyn.waves <- matrix(data = 0,nrow = N,ncol = length(time))
plot(time,wave.form(mesor = 50,amplitude = 20,period = 15,phase = 16,time=time),type="l",lwd=4,col="white",axes=F,xlab="",ylab="")
#colors <- rainbow(N)
for(i in 1:N)
{
  asyn.waves[i,] <- wave.form(mesor = 50,amplitude = 20+norm.random[i],period = 15,phase = 16,time=time)
  lines(time,wave.form(mesor = 50,amplitude = 20+norm.random[i],period = 15,phase = 16,time=time),type="l",lwd=1,col="azure2")
}

lines(time,colMeans(asyn.waves),type="l",lwd=5,col="blue")






lines(time,wave.form(mesor = 50,amplitude = 20,period = 15,phase = 16),type="l",lwd=4,col="blue")



expression.level <- read.table(file = "creinhardtii_expression_level.txt",header=T,as.is=T)

cre.expression <- as.matrix(expression.level[,2:ncol(expression.level)])
is.matrix(cre.expression)
rownames(cre.expression) <- expression.level$gene


results.cre <- rain(t(cre.expression), deltat=3, period=12, verbose=FALSE, nr.series=2)
sum(results.cre$pVal < 0.05)/nrow(cre.expression)
rhythmic.genes.cre.12 <- rownames(subset(results.cre, pVal < 0.05))
length(rhythmic.genes.cre.12)

## Example genes with two peaks

res.sd.12 <- sort.int(results.sd.12$pVal,decreasing = F,index.return = T)
res.nd.12 <- sort.int(results.nd.12$pVal,decreasing = F,index.return = T)


gene.set <- intersect(rownames(results.sd.12)[res.sd.12$ix[1:300]],
rownames(results.nd.12)[res.nd.12$ix[1:300]])

gene.set <- rownames(results.sd.12)[res.sd.12$ix[1:500]]

common.12.sd.nd <- intersect(rhythmic.genes.nd.12,rhythmic.genes.sd.12)
gene.set <- intersect(gene.set,common.12.sd.nd)
length(gene.set)

i <- 1
gene.id <- gene.set[i]
current.gene.expression.sd <- gene.expression[gene.id,31:49]

min.expression <- min(current.gene.expression.sd)
max.expression <- max(current.gene.expression.sd)
range.expression <- max.expression - min.expression

expression.step <- floor(range.expression / 5)

plot(current.gene.expression.sd,type="o",lwd=8,col="red",axes=F,xlab="",ylab="",
     ylim=c(min.expression-expression.step,max.expression),
     cex.lab=1.3,main=gene.id,cex.main=2)
axis(side=2,lwd=4,cex.axis=1.5,las=1)
axis(side = 1,pos = min.expression - 1.1*expression.step, at = seq(from=1,to=19),
     labels = rep("",19),las=2,lwd=3)

polygon(x = c(1,3,3,1),y=c(min.expression-expression.step/2,
                           min.expression-expression.step/2,
                           min.expression-expression.step,
                           min.expression-expression.step),lwd=4,border="red")

polygon(x = c(3,7,7,3),y=c(min.expression-expression.step/2,
                           min.expression-expression.step/2,
                           min.expression-expression.step,
                           min.expression-expression.step),lwd=4,border="red",col="red")

polygon(x = c(7,9,9,7),y=c(min.expression-expression.step/2,
                           min.expression-expression.step/2,
                           min.expression-expression.step,
                           min.expression-expression.step),lwd=4,border="red")

polygon(x = c(9,13,13,9),y=c(min.expression-expression.step/2,
                             min.expression-expression.step/2,
                             min.expression-expression.step,
                             min.expression-expression.step),lwd=4,border="red",col="red")

polygon(x = c(13,15,15,13),y=c(min.expression-expression.step/2,
                               min.expression-expression.step/2,
                               min.expression-expression.step,
                               min.expression-expression.step),lwd=4,border="red")

polygon(x = c(15,19,19,15),y=c(min.expression-expression.step/2,
                               min.expression-expression.step/2,
                               min.expression-expression.step,
                               min.expression-expression.step),lwd=4,border="red",col="red")

## ND expression
nd.zt <- paste("ZT",seq(from=0,to=21,by=3),sep="")



current.gene.expression.nd <- nd.gene.expression[gene.id,c(paste(nd.zt,1,sep="_"),paste(nd.zt,2,sep="_"),paste(nd.zt,3,sep="_"))]
gene.id <- "ostta03g03500"
current.gene.expression.nd <- nd.gene.expression["ostta03g03500",c(paste(nd.zt,1,sep="_"),paste(nd.zt,2,sep="_"),paste(nd.zt,3,sep="_"))]


min.expression <- min(current.gene.expression.nd)
max.expression <- max(current.gene.expression.nd)
range.expression <- max.expression - min.expression

expression.step <- ceiling(range.expression / 4)/2

plot(current.gene.expression.nd,type="o",lwd=8,col="black",axes=F,xlab="",ylab="",
     ylim=c(min.expression-expression.step,1.05*max.expression),
     cex.lab=1.3,main=gene.id,cex.main=2)
axis(side=2,lwd=4,cex.axis=1.5,las=1)
axis(side = 1,pos = min.expression - 1.1*expression.step, at = seq(from=1,to=24),
     labels = rep("",24),las=2,lwd=3)

polygon(x = c(1,5,5,1),y=c(min.expression-expression.step/2,
                           min.expression-expression.step/2,
                           min.expression-expression.step,
                           min.expression-expression.step),lwd=4,border="black")

polygon(x = c(5,9,9,5),y=c(min.expression-expression.step/2,
                           min.expression-expression.step/2,
                           min.expression-expression.step,
                           min.expression-expression.step),lwd=4,border="black",col="black")

polygon(x = c(9,13,13,9),y=c(min.expression-expression.step/2,
                           min.expression-expression.step/2,
                           min.expression-expression.step,
                           min.expression-expression.step),lwd=4,border="black")

polygon(x = c(13,17,17,13),y=c(min.expression-expression.step/2,
                             min.expression-expression.step/2,
                             min.expression-expression.step,
                             min.expression-expression.step),lwd=4,border="black",col="black")

polygon(x = c(17,21,21,17),y=c(min.expression-expression.step/2,
                               min.expression-expression.step/2,
                               min.expression-expression.step,
                               min.expression-expression.step),lwd=4,border="black")

polygon(x = c(21,24,24,21),y=c(min.expression-expression.step/2,
                               min.expression-expression.step/2,
                               min.expression-expression.step,
                               min.expression-expression.step),lwd=4,border="black",col="black")

## LD expression
ld.zt <- paste("ld",paste0("zt",sprintf(fmt = "%02d",seq(from=0,to=20,by=4))),sep="_")
current.gene.expression.ld <- gene.expression[gene.id,c(paste(ld.zt,1,sep="_"),paste(ld.zt,2,sep="_"),paste(ld.zt,3,sep="_"),"ld_zt00_4")]

min.expression <- min(current.gene.expression.ld)
max.expression <- max(current.gene.expression.ld)
range.expression <- max.expression - min.expression

expression.step <- floor(range.expression / 5)

plot(current.gene.expression.ld,type="o",lwd=8,col="blue",axes=F,xlab="",ylab="",
     ylim=c(min.expression-expression.step,max.expression),
     cex.lab=1.3,main=gene.id,cex.main=2)
axis(side=2,lwd=4,cex.axis=1.5,las=1)
axis(side = 1,pos = min.expression - 1.1*expression.step, at = seq(from=1,to=19),
     labels = rep("",19),las=2,lwd=3,cex.axis=1.5)

#labels = rep(paste("ZT",seq(from=0,to=20,by=4)),5),las=2,lwd=3,cex.axis=1.5)

polygon(x = c(1,5,5,1),y=c(min.expression-expression.step/2,
                           min.expression-expression.step/2,
                           min.expression-expression.step,
                           min.expression-expression.step),lwd=4,border="blue")

polygon(x = c(5,7,7,5),y=c(min.expression-expression.step/2,
                           min.expression-expression.step/2,
                           min.expression-expression.step,
                           min.expression-expression.step),lwd=4,border="blue",col="blue")

polygon(x = c(7,11,11,7),y=c(min.expression-expression.step/2,
                             min.expression-expression.step/2,
                             min.expression-expression.step,
                             min.expression-expression.step),lwd=4,border="blue")

polygon(x = c(11,13,13,11),y=c(min.expression-expression.step/2,
                               min.expression-expression.step/2,
                               min.expression-expression.step,
                               min.expression-expression.step),lwd=4,border="blue",col="blue")

polygon(x = c(13,17,17,13),y=c(min.expression-expression.step/2,
                               min.expression-expression.step/2,
                               min.expression-expression.step,
                               min.expression-expression.step),lwd=4,border="blue")

polygon(x = c(17,19,19,17),y=c(min.expression-expression.step/2,
                               min.expression-expression.step/2,
                               min.expression-expression.step,
                               min.expression-expression.step),lwd=4,border="blue",col="blue")

i <- i + 1

#ostta07g03850 This is a family of eukaryotic membrane proteins which incorporate serine into membranes and facilitate the synthesis of the serine-derived lipids phosphatidylserine and sphingolipid 
#ostta10g02350 Glycosyl transferase family 2


i <- 1
gene.id <- gene.set[i]
current.gene.expression.sd <- gene.expression[gene.id,31:60]

min.expression <- min(current.gene.expression.sd)
max.expression <- max(current.gene.expression.sd)
range.expression <- max.expression - min.expression

expression.step <- floor(range.expression / 5)

plot(current.gene.expression.sd,type="o",lwd=8,col="red",axes=F,xlab="",ylab="",
     ylim=c(min.expression-expression.step,max.expression),
     cex.lab=1.3,main=gene.id,cex.main=2)
axis(side=2,lwd=4,cex.axis=1.5,las=1)
axis(side = 1,pos = min.expression - 1.1*expression.step, at = seq(from=1,to=30),
     labels = rep("",30),las=2,lwd=3)

polygon(x = c(1,3,3,1),y=c(min.expression-expression.step/2,
                           min.expression-expression.step/2,
                           min.expression-expression.step,
                           min.expression-expression.step),lwd=4,border="red")

polygon(x = c(3,7,7,3),y=c(min.expression-expression.step/2,
                           min.expression-expression.step/2,
                           min.expression-expression.step,
                           min.expression-expression.step),lwd=4,border="red",col="red")

polygon(x = c(7,9,9,7),y=c(min.expression-expression.step/2,
                           min.expression-expression.step/2,
                           min.expression-expression.step,
                           min.expression-expression.step),lwd=4,border="red")

polygon(x = c(9,13,13,9),y=c(min.expression-expression.step/2,
                             min.expression-expression.step/2,
                             min.expression-expression.step,
                             min.expression-expression.step),lwd=4,border="red",col="red")

polygon(x = c(13,15,15,13),y=c(min.expression-expression.step/2,
                               min.expression-expression.step/2,
                               min.expression-expression.step,
                               min.expression-expression.step),lwd=4,border="red")

polygon(x = c(15,19,19,15),y=c(min.expression-expression.step/2,
                               min.expression-expression.step/2,
                               min.expression-expression.step,
                               min.expression-expression.step),lwd=4,border="red",col="red")

polygon(x = c(19,21,21,19),y=c(min.expression-expression.step/2,
                               min.expression-expression.step/2,
                               min.expression-expression.step,
                               min.expression-expression.step),lwd=4,border="red",col="white")

polygon(x = c(21,25,25,21),y=c(min.expression-expression.step/2,
                               min.expression-expression.step/2,
                               min.expression-expression.step,
                               min.expression-expression.step),lwd=4,border="red",col="salmon")

polygon(x = c(25,27,27,25),y=c(min.expression-expression.step/2,
                               min.expression-expression.step/2,
                               min.expression-expression.step,
                               min.expression-expression.step),lwd=4,border="red",col="white")

polygon(x = c(27,30,30,27),y=c(min.expression-expression.step/2,
                               min.expression-expression.step/2,
                               min.expression-expression.step,
                               min.expression-expression.step),lwd=4,border="red",col="salmon")


## LD expression
ld.zt <- paste("ld",paste0("zt",sprintf(fmt = "%02d",seq(from=0,to=20,by=4))),sep="_")
current.gene.expression.ld <- gene.expression[gene.id,c(paste(ld.zt,1,sep="_"),paste(ld.zt,2,sep="_"),paste(ld.zt,3,sep="_"),paste(ld.zt,4,sep="_"),paste(ld.zt,5,sep="_"))]

min.expression <- min(current.gene.expression.ld)
max.expression <- max(current.gene.expression.ld)
range.expression <- max.expression - min.expression

expression.step <- floor(range.expression / 5)

plot(current.gene.expression.ld,type="o",lwd=8,col="blue",axes=F,xlab="",ylab="",
     ylim=c(min.expression-expression.step,max.expression),
     cex.lab=1.3,main=gene.id,cex.main=2)
axis(side=2,lwd=4,cex.axis=1.5,las=1)
axis(side = 1,pos = min.expression - 1.1*expression.step, at = seq(from=1,to=30),
     labels = rep("",30),las=2,lwd=3,cex.axis=1.5)

#labels = rep(paste("ZT",seq(from=0,to=20,by=4)),5),las=2,lwd=3,cex.axis=1.5)

polygon(x = c(1,5,5,1),y=c(min.expression-expression.step/2,
                           min.expression-expression.step/2,
                           min.expression-expression.step,
                           min.expression-expression.step),lwd=4,border="blue")

polygon(x = c(5,7,7,5),y=c(min.expression-expression.step/2,
                           min.expression-expression.step/2,
                           min.expression-expression.step,
                           min.expression-expression.step),lwd=4,border="blue",col="blue")

polygon(x = c(7,11,11,7),y=c(min.expression-expression.step/2,
                             min.expression-expression.step/2,
                             min.expression-expression.step,
                             min.expression-expression.step),lwd=4,border="blue")

polygon(x = c(11,13,13,11),y=c(min.expression-expression.step/2,
                               min.expression-expression.step/2,
                               min.expression-expression.step,
                               min.expression-expression.step),lwd=4,border="blue",col="blue")

polygon(x = c(13,17,17,13),y=c(min.expression-expression.step/2,
                               min.expression-expression.step/2,
                               min.expression-expression.step,
                               min.expression-expression.step),lwd=4,border="blue")

polygon(x = c(17,19,19,17),y=c(min.expression-expression.step/2,
                               min.expression-expression.step/2,
                               min.expression-expression.step,
                               min.expression-expression.step),lwd=4,border="blue",col="blue")

polygon(x = c(19,23,23,19),y=c(min.expression-expression.step/2,
                               min.expression-expression.step/2,
                               min.expression-expression.step,
                               min.expression-expression.step),lwd=4,border="blue",col="white")

polygon(x = c(23,25,25,23),y=c(min.expression-expression.step/2,
                               min.expression-expression.step/2,
                               min.expression-expression.step,
                               min.expression-expression.step),lwd=4,border="blue",col="lightblue")

polygon(x = c(25,29,29,25),y=c(min.expression-expression.step/2,
                               min.expression-expression.step/2,
                               min.expression-expression.step,
                               min.expression-expression.step),lwd=4,border="blue",col="white")

polygon(x = c(29,30,30,29),y=c(min.expression-expression.step/2,
                               min.expression-expression.step/2,
                               min.expression-expression.step,
                               min.expression-expression.step),lwd=4,border="blue",col="lightblue")

i <- i + 1


i <- 1

gene.id <- rhythmic.genes.cre.12[i]

cre.zt <- paste0("ZT",sprintf(fmt = "%02d",seq(from=0,to=21,by=3)))

current.cre.gene.expression <- cre.expression[gene.id,c(paste(cre.zt,0,sep="_"),paste(cre.zt,1,sep="_"))]

min.expression <- min(current.cre.gene.expression)
max.expression <- max(current.cre.gene.expression)
range.expression <- max.expression - min.expression

expression.step <- ceiling(range.expression / 4)/2

plot(c(current.cre.gene.expression,current.cre.gene.expression[1]),type="o",lwd=8,col="black",axes=F,xlab="",ylab="",
     ylim=c(min.expression-expression.step,1.05*max.expression),
     cex.lab=1.3,main=gene.id,cex.main=2)
axis(side=2,lwd=4,cex.axis=1.5,las=1)
axis(side = 1,pos = min.expression - 1.1*expression.step, at = seq(from=1,to=17),
     labels = rep("",17),las=2,lwd=3)

polygon(x = c(1,5,5,1),y=c(min.expression-expression.step/2,
                           min.expression-expression.step/2,
                           min.expression-expression.step,
                           min.expression-expression.step),lwd=4,border="black")

polygon(x = c(5,9,9,5),y=c(min.expression-expression.step/2,
                           min.expression-expression.step/2,
                           min.expression-expression.step,
                           min.expression-expression.step),lwd=4,border="black",col="black")

polygon(x = c(9,13,13,9),y=c(min.expression-expression.step/2,
                             min.expression-expression.step/2,
                             min.expression-expression.step,
                             min.expression-expression.step),lwd=4,border="black")

polygon(x = c(13,17,17,13),y=c(min.expression-expression.step/2,
                               min.expression-expression.step/2,
                               min.expression-expression.step,
                               min.expression-expression.step),lwd=4,border="black",col="black")
print(gene.id)
i <- i + 1



genes.losing.rhythmicity.ld <- setdiff(complete.ld.rhythmic.genes, complete.ld.ll.rhythmic.genes)

length(intersect(rhythmic.genes.sd.12,genes.losing.rhythmicity.ld))
length(genes.losing.rhythmicity.ld)
length(rhythmic.genes.sd.12)

genes.losing.rhythmicity.sd <-setdiff(complete.sd.rhythmic.genes,complete.sd.ll.rhythmic.genes)

length(intersect(rhythmic.genes.sd.12,genes.losing.rhythmicity.sd))
length(genes.losing.rhythmicity.sd)
length(rhythmic.genes.sd.12)






plot.ld(gene.id="ostta16g02620", gene.name, gene.expression)


## Comparing LD and SD with circacompare
library(circacompare)
ld.zt <- paste("ld",paste0("zt",sprintf(fmt = "%02d",seq(from=0,to=20,by=4))),sep="_")
sd.zt <- paste("sd",paste0("zt",sprintf(fmt = "%02d",seq(from=0,to=20,by=4))),sep="_")

circacompare.ld.sd <- matrix(nrow=length(complete.ld.sd.rhythmic.genes),ncol=15)
rownames(circacompare.ld.sd) <- complete.ld.sd.rhythmic.genes

for(i in 1:length(complete.ld.sd.rhythmic.genes))
{
  gene.i <- complete.ld.sd.rhythmic.genes[i]
  
  ld.expression.i <- gene.expression[gene.i,c(paste(ld.zt,1,sep="_"),paste(ld.zt,2,sep="_"),paste(ld.zt,3,sep="_"))]
  sd.expression.i <- gene.expression[gene.i,c(paste(sd.zt,1,sep="_"),paste(sd.zt,2,sep="_"),paste(sd.zt,3,sep="_"))]
  
  time.points <- seq(from=0,by=4,length.out = 18)
  
  ld.sd.df <- data.frame(time=c(time.points,time.points),
                         measure=c(ld.expression.i, sd.expression.i),
                         group=c(rep("ld",18),rep("sd",18)))
  
  out.i <- circacompare(x = ld.sd.df, col_time = "time", col_group = "group", col_outcome = "measure",alpha_threshold = 1)
  circacompare.ld.sd[i,] <- out.i[[2]][,2]
}

colnames(circacompare.ld.sd) <- out.i[[2]][,1]

par(lwd=3)
boxplot(circacompare.ld.sd[,"ld amplitude estimate" ],circacompare.ld.sd[,"sd amplitude estimate" ],outline=F,col=c("blue","red"),names=c("LD","SD"),cex.axis=2)


sum(circacompare.ld.sd[,"Amplitude difference estimate" ] < 0)
sum(circacompare.ld.sd[,"Amplitude difference estimate" ] < 0) / length(complete.ld.sd.rhythmic.genes)
sum(circacompare.ld.sd[,"P-value for amplitude difference" ] < 0.05)
sum(circacompare.ld.sd[,"P-value for amplitude difference" ] < 0.05)/length(complete.ld.sd.rhythmic.genes)

sum(circacompare.ld.sd[,"Amplitude difference estimate" ] > 0)
sum(circacompare.ld.sd[,"Amplitude difference estimate" ] > 0) / length(complete.ld.sd.rhythmic.genes)
sum(circacompare.ld.sd[circacompare.ld.sd[,"Amplitude difference estimate" ] > 0,"P-value for amplitude difference" ] < 0.05)
write(x = names(which(circacompare.ld.sd[circacompare.ld.sd[,"Amplitude difference estimate" ] > 0,"P-value for amplitude difference" ] < 0.05)),file="genes_increasing_amplitude_in_SD.txt")

sum(circacompare.ld.sd[,"P-value for amplitude difference" ] < 0.05)/length(complete.ld.sd.rhythmic.genes)

par(lwd=3)
boxplot(circacompare.ld.sd[,"ld mesor estimate" ],circacompare.ld.sd[,"sd mesor estimate" ],outline=F,col=c("blue","red"),names=c("LD","SD"),cex.axis=2)
sum(circacompare.ld.sd[,"Mesor difference estimate" ] < 0) /length(complete.ld.sd.rhythmic.genes)
sum(circacompare.ld.sd[,"Mesor difference estimate" ] > 0) /length(complete.ld.sd.rhythmic.genes)
sum(circacompare.ld.sd[,"P-value for mesor difference" ] < 0.05)/length(complete.ld.sd.rhythmic.genes)


sum(circacompare.ld.sd[,"Phase difference estimate" ] < 0)
sum(circacompare.ld.sd[,"Phase difference estimate" ] < 0) / length(complete.ld.sd.rhythmic.genes)
sum(circacompare.ld.sd[,"P-value for difference in phase" ] < 0.05)
sum(circacompare.ld.sd[circacompare.ld.sd[,"Phase difference estimate" ] < 0,"P-value for difference in phase" ] < 0.05)/length(complete.ld.sd.rhythmic.genes)

par(lwd=3)
hist(circacompare.ld.sd[,"ld peak time" ],col="blue",
     ylim=c(0,1200),xlim=c(0,24),main="",axes=FALSE,xlab="",ylab="")
axis(side = 2,lwd=3,cex.axis=1.3,las=2) #530 450

par(lwd=3)
hist(circacompare.ld.sd[,"sd peak time" ],col="red",
     ylim=c(0,1200),xlim=c(0,24),main="",axes=FALSE,xlab="",ylab="")
axis(side = 2,lwd=3,cex.axis=1.3,las=2) #530 450


head(circacompare.ld.sd)

mean(circacompare.ld.sd[,"Phase difference estimate"])

sum(circacompare.ld.sd[,"P-value for difference in phase" ] < 0.05)/length(complete.ld.sd.rhythmic.genes)

boxplot(circacompare.ld.sd[,"ld peak time" ],circacompare.ld.sd[,"sd peak time" ],outline=F,col=c("blue","red"),names=c("LD","SD"),cex.axis=2)


hist(circacompare.ld.sd[,"sd peak time" ],col="red",ylim=c(0,1200),breaks=24)

box()


boxplot(circacompare.ld.sd[,"ld peak time" ],circacompare.ld.sd[,"sd peak time" ],outline=F,col=c("blue","red"),names=c("LD","SD"),cex.axis=2)
hist(circacompare.ld.sd[,"ld peak time" ],col=rgb(0,0,1,0.5),ylim=c(0,1200))
hist(circacompare.ld.sd[,"sd peak time" ],col=rgb(1,0,0,0.5),ylim=c(0,1200),add=T)


h2<-rnorm(1000,4)
h1<-rnorm(1000,6)

hist(h1, col=rgb(1,0,0,0.5),xlim=c(0,10), ylim=c(0,200), main="Overlapping Histogram", xlab="Variable")
hist(h2, col=rgb(0,0,1,0.5), add=T)
box()


sum(circacompare.ld.sd[,"Phase difference estimate" ] < 0) /length(complete.ld.sd.rhythmic.genes)
sum(circacompare.ld.sd[,"Phase difference estimate" ] > 0) /length(complete.ld.sd.rhythmic.genes)
sum(circacompare.ld.sd[,"P-value for difference in phase" ] < 0.05)/length(complete.ld.sd.rhythmic.genes)


sum(circacompare.ld.sd[,"Phase difference estimate" ] < 2 & circacompare.ld.sd[,"Phase difference estimate" ] > - 2) /length(complete.ld.sd.rhythmic.genes)
sum(circacompare.ld.sd[,"Phase difference estimate" ] < 2 & circacompare.ld.sd[,"Phase difference estimate" ] > - 2)

non.responsive.genes.phase <- names(which(circacompare.ld.sd[,"Phase difference estimate" ] < 2 & circacompare.ld.sd[,"Phase difference estimate" ] > - 2))
# non.responsive.genes.mesor <- names(which(circacompare.ld.sd[,"Mesor difference estimate" ] < 5 & circacompare.ld.sd[,"Mesor difference estimate" ] > - 5))
# non.responsive.genes.amplitude <- names(which(circacompare.ld.sd[,"Amplitude difference estimate" ] < 5 & circacompare.ld.sd[,"Amplitude difference estimate" ] > - 5))

non.responsive.genes.mesor <- names(which((abs(circacompare.ld.sd[, "ld mesor estimate"] - circacompare.ld.sd[, "sd mesor estimate"]) / circacompare.ld.sd[, "ld mesor estimate"]) < 0.1))
non.responsive.genes.amplitude <-names(which((abs(circacompare.ld.sd[, "ld amplitude estimate"] - circacompare.ld.sd[, "sd amplitude estimate"]) / circacompare.ld.sd[, "ld mesor estimate"]) < 0.1))


non.responsive.genes <- intersect(intersect(non.responsive.genes.mesor,non.responsive.genes.phase),non.responsive.genes.amplitude)
length(non.responsive.genes)

i <- 1

gene.i <- non.responsive.genes[i]

ld.expression.i <- gene.expression[gene.i,c(paste(ld.zt,1,sep="_"),paste(ld.zt,2,sep="_"),paste(ld.zt,3,sep="_"))]
sd.expression.i <- gene.expression[gene.i,c(paste(sd.zt,1,sep="_"),paste(sd.zt,2,sep="_"),paste(sd.zt,3,sep="_"))]

time.points <- seq(from=0,by=4,length.out = 18)

ld.sd.df <- data.frame(time=c(time.points,time.points),
                       measure=c(ld.expression.i, sd.expression.i),
                       group=c(rep("ld",18),rep("sd",18)))

circacompare(x = ld.sd.df, col_time = "time", col_group = "group", col_outcome = "measure",alpha_threshold = 1)

i <- i + 1



getwd()
write.table(x = non.responsive.genes,file = "non_responsive_genes.txt",quote = F,row.names = F,col.names = F)



plot.ld.sd <- function(gene.id, gene.name, gene.expression)
{

  ld.zt <- paste("ld",paste0("zt",sprintf(fmt = "%02d",seq(from=0,to=20,by=4))),sep="_")
  current.gene.expression.ld <- gene.expression[gene.id,c(paste(ld.zt,1,sep="_"),paste(ld.zt,2,sep="_"),paste(ld.zt,3,sep="_"))]
  
  sd.zt <- paste("sd",paste0("zt",sprintf(fmt = "%02d",seq(from=0,to=20,by=4))),sep="_")
  current.gene.expression.sd <- gene.expression[gene.id,c(paste(sd.zt,1,sep="_"),paste(sd.zt,2,sep="_"),paste(sd.zt,3,sep="_"))]

  max.expr <- max(c(current.gene.expression.ld, current.gene.expression.sd))
  min.expr <- min(c(current.gene.expression.ld, current.gene.expression.sd))
  
  range.expr <- (max.expr - min.expr)
  
  plot(x = -10,y= -10,axes=F,xlab="",ylab="",
       ylim=c(min.expr-0.2*range.expr, max.expr+0.2*range.expr),xlim=c(0,72),
       main=paste(gene.id, gene.name,sep=" - "),cex.main=2)
  lines(x = seq(from=0,by=4,to=68),current.gene.expression.ld,type="o",lwd=5,col="blue")
  lines(x = seq(from=0,by=4,to=68),current.gene.expression.sd,type="o",lwd=5,col="red")
  axis(side=2,lwd=3,las=2,cex.axis=1.5)

  for(i in 0:2)
  {
    current.line <- (max.expr - min.expr)/20
    
    polygon(x = c(24*i, 24*i+16, 24*i+16, 24*i),
            y=c(min.expr-current.line, min.expr-current.line, 
                min.expr-2.5*(current.line), min.expr-2.5*(current.line)),
            lwd=3,border="blue")
    polygon(x = c(24*i+16,24*(i+1),24*(i+1),24*i+16),
            y=c(min.expr-current.line, min.expr-current.line, 
                min.expr-2.5*(current.line), min.expr-2.5*(current.line)),
            lwd=3,border="blue",col="blue")
    
    polygon(x = c(24*i, 24*i+8, 24*i+8, 24*i),
            y=c(min.expr-3*current.line, min.expr-3*current.line, 
                min.expr-4.5*(current.line), min.expr-4.5*(current.line)),
            lwd=3,border="red")
    polygon(x = c(24*i+8,24*(i+1),24*(i+1),24*i+8),
            y=c(min.expr-3*current.line, min.expr-3*current.line, 
                min.expr-4.5*(current.line), min.expr-4.5*(current.line)),
            lwd=3,border="red",col="red")
  }
  
  return(0)  
}

genes.amplitude.phase <- intersect(genes.decrease.amplitude.effect.ld.sd,genes.phase.effect.ld.sd)

i <- 1
plot.ld.sd(gene.id = genes.amplitude.phase[i],gene.name = "MCM6", gene.expression)
print(genes.amplitude.phase[i])
print(i)
i <- i + 1

plot.ld.sd(gene.id = "ostta14g01410",gene.name = "EZH2", gene.expression)


plot.ld.sd(gene.id = "ostta03g04220",gene.name = "HAT", gene.expression)

plot.ld.ll(gene.id = "ostta03g04220",gene.name = "HAT", gene.expression)
plot.ld.dd(gene.id = "ostta03g04220",gene.name = "HAT", gene.expression)


plot.sd.ll(gene.id = "ostta03g04220",gene.name = "HAT", gene.expression)
plot.sd.dd(gene.id = "ostta03g04220",gene.name = "HAT", gene.expression)


plot.ld.sd(gene.id = "ostta14g01410",gene.name = "EZH2", gene.expression)

plot.ld.ll(gene.id = "ostta14g01410",gene.name = "EZH2", gene.expression)
plot.ld.dd(gene.id = "ostta14g01410",gene.name = "EZH2", gene.expression)


plot.sd.ll(gene.id = "ostta14g01410",gene.name = "EZH2", gene.expression)
plot.sd.dd(gene.id = "ostta14g01410",gene.name = "EZH2", gene.expression)


plot.ld.sd(gene.id = "ostta02g01510",gene.name = "HDAC", gene.expression)




plot.ld.sd <- function(gene.id, gene.name, gene.expression)
{
  
  ld.zt <- paste("ld",paste0("zt",sprintf(fmt = "%02d",seq(from=0,to=20,by=4))),sep="_")
  current.gene.expression.ld <- gene.expression[gene.id,c(paste(ld.zt,1,sep="_"),paste(ld.zt,2,sep="_"),paste(ld.zt,3,sep="_"))]
  
  sd.zt <- paste("sd",paste0("zt",sprintf(fmt = "%02d",seq(from=0,to=20,by=4))),sep="_")
  current.gene.expression.sd <- gene.expression[gene.id,c(paste(sd.zt,1,sep="_"),paste(sd.zt,2,sep="_"),paste(sd.zt,3,sep="_"))]
  
  max.expr <- max(c(current.gene.expression.ld, current.gene.expression.sd))
  min.expr <- min(c(current.gene.expression.ld, current.gene.expression.sd))
  
  range.expr <- (max.expr - min.expr)
  
  plot(x = -10,y= -10,axes=F,xlab="",ylab="",
       ylim=c(min.expr-0.2*range.expr, max.expr+0.2*range.expr),xlim=c(0,72),
       main=paste(gene.id, gene.name,sep=" - "),cex.main=2)
  lines(x = seq(from=0,by=4,to=68),current.gene.expression.ld,type="o",lwd=5,col="blue")
  lines(x = seq(from=0,by=4,to=68),current.gene.expression.sd,type="o",lwd=5,col="red")
  axis(side=2,lwd=3,las=2,cex.axis=1.5)
  
  for(i in 0:2)
  {
    current.line <- (max.expr - min.expr)/20
    
    polygon(x = c(24*i, 24*i+16, 24*i+16, 24*i),
            y=c(min.expr-current.line, min.expr-current.line, 
                min.expr-2.5*(current.line), min.expr-2.5*(current.line)),
            lwd=3,border="blue")
    polygon(x = c(24*i+16,24*(i+1),24*(i+1),24*i+16),
            y=c(min.expr-current.line, min.expr-current.line, 
                min.expr-2.5*(current.line), min.expr-2.5*(current.line)),
            lwd=3,border="blue",col="blue")
    
    polygon(x = c(24*i, 24*i+8, 24*i+8, 24*i),
            y=c(min.expr-3*current.line, min.expr-3*current.line, 
                min.expr-4.5*(current.line), min.expr-4.5*(current.line)),
            lwd=3,border="red")
    polygon(x = c(24*i+8,24*(i+1),24*(i+1),24*i+8),
            y=c(min.expr-3*current.line, min.expr-3*current.line, 
                min.expr-4.5*(current.line), min.expr-4.5*(current.line)),
            lwd=3,border="red",col="red")
  }
  
  return(0)  
}

## UN PICO
#ostta01g00280
#ostta01g00330
#ostta01g00360
#ostta01g01220
#ostta01g01420
#ostta01g01560 muy bueno
#ostta01g01790 muy bueno
#ostta01g02250 muy bueno
#ostta01g03390 muy bueno
#ostta01g03960 muy bueno
#ostta01g05210 muy bueno
#ostta01g05750 muy bueno
#ostta01g06150 muy bueno CYCB
plot.ld.sd(gene.id="ostta01g06150", gene.name="CYCB", gene.expression)
#ostta01g06380 muy bueno Regulator of chromosome condensation, RCC1
plot.ld.sd(gene.id="ostta01g06380", gene.name="RCC1", gene.expression)
plot.ld.sd(gene.id="ostta01g01790", gene.name="RCC1", gene.expression)
#ostta02g02090 muy bueno

plot.ld.sd(gene.id="ostta01g06150", gene.name="CYCB", gene.expression)
plot.ld.sd(gene.id="ostta01g00790", gene.name="ADS1", gene.expression)

plot.ld.sd(gene.id="ostta09g02190", gene.name="HisKinLOV", gene.expression)

plot.sd.ll(gene.id="ostta09g02190", gene.name="ostta09g02190", gene.expression)

#ostta01g01700 doble pico
#ostta01g02130 doble pico
#ostta01g02270 doble pico
#ostta01g02210 doble pico
#ostta01g02390 doble pico
#ostta01g03060 doble pico
#ostta01g05010 doble pico
#ostta01g05220 doble pico
# 

## DOOS PICOS
plot.ld.sd(gene.id="ostta01g01700", gene.name="ostta01g01700", gene.expression)
plot.ld.ll(gene.id="ostta01g01700", gene.name="ostta01g01700", gene.expression)
plot.ld.dd(gene.id="ostta01g01700", gene.name="ostta01g01700", gene.expression)
plot.sd.ll(gene.id="ostta01g01700", gene.name="ostta01g01700", gene.expression)
plot.sd.dd(gene.id="ostta01g01700", gene.name="ostta01g01700", gene.expression)

current.gene <- "ostta04g04000"
plot.ld.sd(gene.id=current.gene, gene.name=current.gene, gene.expression)
plot.ld.ll(gene.id=current.gene, gene.name=current.gene, gene.expression)
plot.ld.dd(gene.id=current.gene, gene.name=current.gene, gene.expression)
plot.sd.ll(gene.id=current.gene, gene.name=current.gene, gene.expression)
plot.sd.dd(gene.id=current.gene, gene.name=current.gene, gene.expression)

current.gene <- "ostta04g02740"
plot.ld.sd(gene.id=current.gene, gene.name=current.gene, gene.expression)
plot.ld.ll(gene.id=current.gene, gene.name=current.gene, gene.expression)
plot.ld.dd(gene.id=current.gene, gene.name=current.gene, gene.expression)
plot.sd.ll(gene.id=current.gene, gene.name=current.gene, gene.expression)
plot.sd.dd(gene.id=current.gene, gene.name=current.gene, gene.expression)


plot.ld.sd(gene.id="ostta01g02130", gene.name="ostta01g02130", gene.expression)
plot.ld.sd(gene.id="ostta01g03060", gene.name="ostta01g03060", gene.expression)
plot.ld.sd(gene.id="ostta01g05010", gene.name="ostta01g05010", gene.expression)
plot.ld.sd(gene.id="ostta01g05220", gene.name="ostta01g05220", gene.expression)


plot.ld.sd(gene.id="ostta01g02390", gene.name="ostta01g02390", gene.expression)
plot.ld.sd(gene.id="ostta01g02210", gene.name="ostta01g02210", gene.expression)
plot.ld.sd(gene.id="ostta01g02270", gene.name="ostta01g02270", gene.expression)

## Clustering LD
cluster.circadian.genes <- vector(mode="character", length=length(complete.ld.sd.rhythmic.genes))

classify.genes <- function(rain.results, gene.set)
{
  circadian.rain.results <- rain.results[gene.set,]
  
  gene.peaks <- circadian.rain.results$phase - 4
  gene.troughs <- (gene.peaks + circadian.rain.results$peak.shape) %% 24
  
  return(data.frame(gene=gene.set,peak=gene.peaks,trough=gene.troughs,stringsAsFactors = F))
}

ld.gene.classification <- classify.genes(rain.results = results.ld, gene.set = complete.ld.sd.rhythmic.genes)

zt0.ld <- subset(ld.gene.classification, peak == 0)$gene
zt4.ld <- subset(ld.gene.classification, peak == 4)$gene
zt8.ld <- subset(ld.gene.classification, peak == 8)$gene
zt12.ld <- subset(ld.gene.classification, peak == 12)$gene
zt16.ld <- subset(ld.gene.classification, peak == 16)$gene
zt20.ld <- subset(ld.gene.classification, peak == 20)$gene

ld.mean.gene.expression <- matrix(ncol = 6,data = c(rowMeans(gene.expression[,c("ld_zt00_1","ld_zt00_2","ld_zt00_3")]),
                                                    rowMeans(gene.expression[,c("ld_zt04_1","ld_zt04_2","ld_zt04_3")]),
                                                    rowMeans(gene.expression[,c("ld_zt08_1","ld_zt08_2","ld_zt08_3")]),
                                                    rowMeans(gene.expression[,c("ld_zt12_1","ld_zt12_2","ld_zt12_3")]),
                                                    rowMeans(gene.expression[,c("ld_zt16_1","ld_zt16_2","ld_zt16_3")]),
                                                    rowMeans(gene.expression[,c("ld_zt20_1","ld_zt20_2","ld_zt20_3")])))
rownames(ld.mean.gene.expression) <- rownames(ld.gene.expression)
colnames(ld.mean.gene.expression) <- c("zt00","zt04","zt08","zt12","zt16","zt20")

ld.mean.expression.sorted <- ld.mean.gene.expression[c(zt0.ld,zt4.ld,zt8.ld,
                                                       zt12.ld,zt16.ld,zt20.ld),]

scale.and.interpolate.mean.expression <- function(mean.expression)
{
  time.points <- seq(from=0,to=24,4)
  sorted.genes <- rownames(mean.expression)
  
  scaled.gene.expression <- matrix(nrow=nrow(mean.expression),ncol=49)
  rownames(scaled.gene.expression) <- sorted.genes
  colnames(scaled.gene.expression) <- paste("ZT",seq(from=0,to=24,by=0.5),sep="_")
  
  i <- 1
  for(i in 1:length(sorted.genes))
  {
    current.gene <- sorted.genes[i]
    
    expression.current.gene.i <- mean.expression[current.gene,]
    #expression.current.gene.i <- short.day.mean.expression[current.gene,]
    
    expression.current.gene.i <- c(expression.current.gene.i, expression.current.gene.i[1])
    #scaled.gene.expression[i,] <- approx(x=time.points,expression.current.gene.i,xout=0:24)[[2]]
    
    scaled.gene.expression[i,] <- approx(x=time.points,expression.current.gene.i,xout=seq(from=0,to=24,by=0.5))[[2]]
    scaled.gene.expression[i,] <- scaled.gene.expression[i,]/max(scaled.gene.expression[i,])
  }
  
  return(scaled.gene.expression)
}

scaled.ld.mean.expression.sorted <- scale.and.interpolate.mean.expression(mean.expression = ld.mean.expression.sorted)

library(gplots)
colfunc <- colorRampPalette(c("black","blue","yellow"))
nrow(scaled.ld.mean.expression.sorted)
png(filename = "heatmap_LD.png",width = 400)
heatmap.2(scaled.ld.mean.expression.sorted,Colv =FALSE,trace=c("none"),dendrogram = "none",col=colfunc(100),density.info=c("none"),key=FALSE,margins = c(2,8),cexRow = 1.2,Rowv=FALSE,labCol = c(""),labRow = c(""))
dev.off()

ind.1 <- seq(from=1,by=3,to=nrow(scaled.ld.mean.expression.sorted))
length(ind.1)
ind.2 <- seq(from=2,by=3,to=nrow(scaled.ld.mean.expression.sorted))
length(ind.2)
ind.3 <- seq(from=3,by=3,to=nrow(scaled.ld.mean.expression.sorted))
length(ind.3)

ld.cir.heatmap <- (scaled.ld.mean.expression.sorted[ind.1[1:1757],] +
  scaled.ld.mean.expression.sorted[ind.2,] +
  scaled.ld.mean.expression.sorted[ind.3,])/3

circos.par("track.height" = 0.3)
col_fun_g1 = colorRamp2(c(0,0.5,1), c("black", "blue", "yellow"))
circos.heatmap(t(ld.cir.heatmap), #data in a matrix
               col = col_fun_g1, #color function defined previosly
               show.sector.labels = F, #labels of track (col names) or levels
               cluster = F) #F to avoid cluster organization 
circos.clear()

  
sd.cir.heatmap <- (scaled.sd.mean.expression.sorted[ind.1[1:1757],] +
                     scaled.sd.mean.expression.sorted[ind.2,] +
                     scaled.sd.mean.expression.sorted[ind.3,])/3

circos.par("track.height" = 0.3)
col_fun_g1 = colorRamp2(c(0,0.5,1), c("black", "blue", "yellow"))
circos.heatmap(t(sd.cir.heatmap), #data in a matrix
               col = col_fun_g1, #color function defined previosly
               show.sector.labels = F, #labels of track (col names) or levels
               cluster = F) #F to avoid cluster organization 
circos.clear()





circos.par("track.height" = 0.3)
col_fun_g1 = colorRamp2(c(0,0.5,1), c("black", "blue", "yellow"))
circos.heatmap(t(scaled.ld.mean.expression.sorted), #data in a matrix
               col = col_fun_g1, #color function defined previosly
               show.sector.labels = F, #labels of track (col names) or levels
               cluster = F) #F to avoid cluster organization 
circos.clear()





sd.gene.classification <- classify.genes(rain.results = results.sd, gene.set = complete.ld.sd.rhythmic.genes)

zt0.sd <- subset(sd.gene.classification, peak == 0)$gene
zt4.sd <- subset(sd.gene.classification, peak == 4)$gene
zt8.sd <- subset(sd.gene.classification, peak == 8)$gene
zt12.sd <- subset(sd.gene.classification, peak == 12)$gene
zt16.sd <- subset(sd.gene.classification, peak == 16)$gene
zt20.sd <- subset(sd.gene.classification, peak == 20)$gene

sd.mean.gene.expression <- matrix(ncol = 6,data = c(rowMeans(gene.expression[,c("sd_zt00_1","sd_zt00_2","sd_zt00_3")]),
                                                    rowMeans(gene.expression[,c("sd_zt04_1","sd_zt04_2","sd_zt04_3")]),
                                                    rowMeans(gene.expression[,c("sd_zt08_1","sd_zt08_2","sd_zt08_3")]),
                                                    rowMeans(gene.expression[,c("sd_zt12_1","sd_zt12_2","sd_zt12_3")]),
                                                    rowMeans(gene.expression[,c("sd_zt16_1","sd_zt16_2","sd_zt16_3")]),
                                                    rowMeans(gene.expression[,c("sd_zt20_1","sd_zt20_2","sd_zt20_3")])))
rownames(sd.mean.gene.expression) <- rownames(sd.gene.expression)
colnames(sd.mean.gene.expression) <- c("zt00","zt04","zt08","zt12","zt16","zt20")

sd.mean.expression.sorted <- sd.mean.gene.expression[c(zt0.sd,zt4.sd,zt8.sd,
                                                       zt12.sd,zt16.sd,zt20.sd),]
scaled.sd.mean.expression.sorted <- scale.and.interpolate.mean.expression(mean.expression = sd.mean.expression.sorted)

png(filename = "heatmap_SD.png",width = 400)
heatmap.2(scaled.sd.mean.expression.sorted,Colv =FALSE,trace=c("none"),dendrogram = "none",col=colfunc(100),density.info=c("none"),key=FALSE,margins = c(2,8),cexRow = 1.2,Rowv=FALSE,labCol = c(""),labRow = c(""))
dev.off()


plot.ld <- function(gene.id, gene.name, gene.expression)
{
  
  ld.zt <- paste("ld",paste0("zt",sprintf(fmt = "%02d",seq(from=0,to=20,by=4))),sep="_")
  current.gene.expression.ld <- as.numeric(gene.expression[gene.id,c(paste(ld.zt,1,sep="_"),paste(ld.zt,2,sep="_"),paste(ld.zt,3,sep="_"))])
  
  #sd.zt <- paste("sd",paste0("zt",sprintf(fmt = "%02d",seq(from=0,to=20,by=4))),sep="_")
  #current.gene.expression.sd <- gene.expression[gene.id,c(paste(sd.zt,1,sep="_"),paste(sd.zt,2,sep="_"),paste(sd.zt,3,sep="_"))]
  
  max.expr <- max(c(current.gene.expression.ld))#, current.gene.expression.sd))
  min.expr <- min(c(current.gene.expression.ld))#, current.gene.expression.sd))
  
  range.expr <- (max.expr - min.expr)
  
  plot(x = -10,y= -10,axes=F,xlab="",ylab="",
       ylim=c(min.expr-0.2*range.expr, max.expr+0.2*range.expr),xlim=c(0,72),
       main=paste(gene.id, gene.name,sep=" - "),cex.main=2)
  lines(x = seq(from=0,by=4,to=68),current.gene.expression.ld,type="o",lwd=5,col="blue")
  #lines(x = seq(from=0,by=4,to=68),current.gene.expression.sd,type="o",lwd=5,col="red")
  axis(side=2,lwd=3,las=2,cex.axis=1.5)
  
  for(i in 0:2)
  {
    current.line <- (max.expr - min.expr)/20
    
    polygon(x = c(24*i, 24*i+16, 24*i+16, 24*i),
            y=c(min.expr-current.line, min.expr-current.line, 
                min.expr-2.5*(current.line), min.expr-2.5*(current.line)),
            lwd=3,border="blue")
    polygon(x = c(24*i+16,24*(i+1),24*(i+1),24*i+16),
            y=c(min.expr-current.line, min.expr-current.line, 
                min.expr-2.5*(current.line), min.expr-2.5*(current.line)),
            lwd=3,border="blue",col="blue")
    
    # polygon(x = c(24*i, 24*i+8, 24*i+8, 24*i),
    #         y=c(min.expr-3*current.line, min.expr-3*current.line, 
    #             min.expr-4.5*(current.line), min.expr-4.5*(current.line)),
    #         lwd=3,border="red")
    # polygon(x = c(24*i+8,24*(i+1),24*(i+1),24*i+8),
    #         y=c(min.expr-3*current.line, min.expr-3*current.line, 
    #             min.expr-4.5*(current.line), min.expr-4.5*(current.line)),
    #         lwd=3,border="red",col="red")
  }
  
  return(0)  
}

plot.ld(gene.id = "ostta14g01410",gene.name = "EZH2", gene.expression)
plot.ld(gene.id = "ostta11g01290",gene.name = "TX", gene.expression)

plot.ld(gene.id = "ostta03g01090",gene.name = "TX", gene.expression)



plot.ld.sd(gene.id = "ostta05g01490",gene.name = "ostta05g01490",gene.expression = gene.expression )
plot.ld.ll(gene.id = "ostta05g01490",gene.name = "ostta05g01490",gene.expression = gene.expression )
plot.ld.dd(gene.id = "ostta05g01490",gene.name = "ostta05g01490",gene.expression = gene.expression )
plot.sd.ll(gene.id = "ostta05g01490",gene.name = "ostta05g01490",gene.expression = gene.expression )
plot.sd.dd(gene.id = "ostta05g01490",gene.name = "ostta05g01490",gene.expression = gene.expression )

plot.ld.sd(gene.id = "ostta08g00300",gene.name = "ostta08g00300",gene.expression = gene.expression )
plot.ld.ll(gene.id = "ostta08g00300",gene.name = "ostta08g00300",gene.expression = gene.expression )
plot.ld.dd(gene.id = "ostta08g00300",gene.name = "ostta08g00300",gene.expression = gene.expression )
plot.sd.ll(gene.id = "ostta08g00300",gene.name = "ostta08g00300",gene.expression = gene.expression )
plot.sd.dd(gene.id = "ostta08g00300",gene.name = "ostta08g00300",gene.expression = gene.expression )


plot.ld.sd(gene.id = "ostta09g02720",gene.name = "ostta09g02720",gene.expression = gene.expression )
plot.ld.sd(gene.id = "ostta09g02720",gene.name = "ostta09g02720",gene.expression = gene.expression )
plot.ld.sd(gene.id = "ostta11g01030",gene.name = "ostta11g01030",gene.expression = gene.expression )
plot.ld.sd(gene.id = "ostta16g02620",gene.name = "ostta16g02620",gene.expression = gene.expression )
plot.ld.sd(gene.id = "ostta18g01947",gene.name = "ostta18g01947",gene.expression = gene.expression )


plot.ld.sd(gene.id = "ostta06g02940",gene.name = "GBSS",gene.expression = gene.expression )
plot.sd.dd(gene.id = "ostta06g02940",gene.name = "GBSS",gene.expression = gene.expression )
plot.sd.ll(gene.id = "ostta06g02940",gene.name = "GBSS",gene.expression = gene.expression )


genes.helen <- c("ostta05g01490", "ostta08g00300","ostta09g02720", "ostta11g01030","ostta16g02620","ostta18g01947")

helen.gene.expression <- gene.expression[genes.helen,]
write.table(x = helen.gene.expression,file = "helen_gene_expression.tsv",row.names = T,col.names = T,sep = "\t",quote=F)


plot.ld.sd(gene.id = "ostta02g03030",gene.name = "CSP",gene.expression = gene.expression )
plot.ld.ll(gene.id = "ostta02g03030",gene.name = "CSP",gene.expression = gene.expression )

plot.ld.sd(gene.id = "ostta13g00710",gene.name = "HSP20",gene.expression = gene.expression )
plot.ld.ll(gene.id = "ostta13g00710",gene.name = "HSP20",gene.expression = gene.expression )

plot.ld.sd(gene.id = "ostta07g03750",gene.name = "JMJC",gene.expression = gene.expression )
plot.ld.ll(gene.id = "ostta07g03750",gene.name = "JMJC",gene.expression = gene.expression )

plot.ld.sd(gene.id = "ostta01g00260",gene.name = "TR1",gene.expression = gene.expression )
plot.ld.ll(gene.id = "ostta01g00260",gene.name = "TR1",gene.expression = gene.expression )

plot.ld.sd(gene.id = "ostta01g06130",gene.name = "TR2",gene.expression = gene.expression )
plot.ld.ll(gene.id = "ostta01g06130",gene.name = "TR2",gene.expression = gene.expression )

plot.ld.sd(gene.id = "ostta19g00045",gene.name = "C5Met",gene.expression = gene.expression )
plot.ld.ll(gene.id = "ostta19g00045",gene.name = "C5Met",gene.expression = gene.expression )


plot.ld.sd(gene.id = "ostta11g02575",gene.name = "Int",gene.expression = gene.expression )
plot.ld.ll(gene.id = "ostta11g02575",gene.name = "Int",gene.expression = gene.expression )

plot.ld.sd(gene.id = "ostta11g02720",gene.name = "?",gene.expression = gene.expression )
plot.ld.ll(gene.id = "ostta11g02720",gene.name = "?",gene.expression = gene.expression )

plot.ld.sd(gene.id = "ostta20g00470",gene.name = "?",gene.expression = gene.expression )
plot.ld.ll(gene.id = "ostta20g00470",gene.name = "?",gene.expression = gene.expression )


plot.ld.sd(gene.id = "ostta11g02490",gene.name = "WRKY",gene.expression = gene.expression )
plot.ld.ll(gene.id = "ostta11g02490",gene.name = "WRKY",gene.expression = gene.expression )



library(circacompare)
ld.zt <- paste("ld",paste0("zt",sprintf(fmt = "%02d",seq(from=0,to=20,by=4))),sep="_")

circacompare.ld.ll <- matrix(nrow=length(complete.ld.ll.rhythmic.genes),ncol=15)
rownames(circacompare.ld.ll) <- complete.ld.ll.rhythmic.genes

for(i in 1:length(complete.ld.ll.rhythmic.genes))
{
  print(i)
  gene.i <- complete.ld.ll.rhythmic.genes[i]
  print(gene.i)
  
  ld.expression.i <- gene.expression[gene.i,c(paste(ld.zt,2,sep="_"),paste(ld.zt,3,sep="_"))]
  ll.expression.i <- gene.expression[gene.i,c(paste(ld.zt,4,sep="_"),paste(ld.zt,5,sep="_"))]
  
  time.points <- seq(from=0,by=4,length.out = 12)
  
  ld.ll.df <- data.frame(time=c(time.points,time.points),
                         measure=c(ld.expression.i, ll.expression.i),
                         group=c(rep("ld",12),rep("ll",12)))
  
  out.i <- circacompare(x = ld.ll.df, col_time = "time", col_group = "group", col_outcome = "measure",alpha_threshold = 1)
  circacompare.ld.ll[i,] <- out.i[[2]][,2]
}

colnames(circacompare.ld.ll) <- out.i[[2]][,1]

par(lwd=3)
boxplot(circacompare.ld.ll[,"ld amplitude estimate" ],circacompare.ld.ll[,"ll amplitude estimate" ],outline=F,col=c("blue","lightblue"),names=c("LD","LL"),cex.axis=2)
sum(circacompare.ld.ll[,"Amplitude difference estimate" ] < 0)
sum(circacompare.ld.ll[,"Amplitude difference estimate" ] < 0) / length(complete.ld.ll.rhythmic.genes)
sum(circacompare.ld.ll[,"P-value for amplitude difference" ] < 0.05)
sum(circacompare.ld.ll[,"P-value for amplitude difference" ] < 0.05)/length(complete.ld.ll.rhythmic.genes)

boxplot(circacompare.ld.ll[,"ld mesor estimate" ],circacompare.ld.ll[,"ll mesor estimate" ],outline=F,col=c("blue","lightblue"),names=c("LD","LL"),cex.axis=2)
sum(circacompare.ld.ll[,"Mesor difference estimate" ] < 0) /length(complete.ld.ll.rhythmic.genes)
sum(circacompare.ld.ll[,"Mesor difference estimate" ] > 0) /length(complete.ld.ll.rhythmic.genes)
sum(circacompare.ld.ll[,"P-value for mesor difference" ] < 0.05)/length(complete.ld.ll.rhythmic.genes)

boxplot(circacompare.ld.ll[,"Phase difference estimate"],outline=F,col=c("lightblue"),cex.axis=2)

sum(circacompare.ld.ll[,] > 0)
sum(circacompare.ld.ll[,"P-value for difference in phase" ] < 0.05)

summary(circacompare.ld.ll[,"Phase difference estimate"])

png(filename = "LD_LL_DD.png",width = 500,height = 500)
grid.newpage()
draw.pairwise.venn(area1 = length(complete.ld.ll.rhythmic.genes),
                   area2 = length(complete.ld.dd.rhythmic.genes),
                   cross.area = length(intersect(complete.ld.ll.rhythmic.genes,complete.ld.dd.rhythmic.genes)),
                   lwd = 6,alpha = 0.7,fill = c("lightblue","darkblue"),cex = 3,category = c("LD/LL","LD/DD"),cat.cex = 3,cat.pos = c(-10,10))
dev.off()

i <- 1                   
current.gene <- ld.dark.repressed[i]
plot.ld.dd(gene.id = current.gene,gene.name = current.gene, gene.expression = gene.expression)
i <- i + 1


i <- 1                   
current.gene <- ld.light.repressed.dark.activated[i]
plot.ld.dd(gene.id = current.gene,gene.name = current.gene, gene.expression = gene.expression)
plot.ld.ll(gene.id = current.gene,gene.name = current.gene, gene.expression = gene.expression)
i <- i + 1


## rhythmicity circacompare
## LD LL
nrow(gene.expression)
number.genes
library(circacompare)
ld.zt <- paste("ld",paste0("zt",sprintf(fmt = "%02d",seq(from=0,to=20,by=4))),sep="_")

genes <- rownames(gene.expression)

circacompare.ld.ll <- matrix(nrow=number.genes,ncol=15)
rownames(circacompare.ld.ll) <- genes

time.points.ld <- seq(from=0,by=4,length.out = 18)
time.points.ll <- seq(from=0,by=4,length.out = 12)

for(i in 1:length(genes))
{
  print(i)
  gene.i <- genes[i]
  
  ld.expression.i <- gene.expression[gene.i,c(paste(ld.zt,1,sep="_"),paste(ld.zt,2,sep="_"),paste(ld.zt,3,sep="_"))]
  ll.expression.i <- gene.expression[gene.i,c(paste(ld.zt,4,sep="_"),paste(ld.zt,5,sep="_"))]
  
  ld.ll.df <- data.frame(time=c(time.points.ld,time.points.ll),
                         measure=c(ld.expression.i, ll.expression.i),
                         group=c(rep("ld",18),rep("ll",12)))
  
  out.i <- circacompare(x = ld.ll.df, col_time = "time", col_group = "group", col_outcome = "measure",alpha_threshold = 1)
  
  if(!is.null(out.i))
  {
    circacompare.ld.ll[i,] <- out.i[[2]][,2]
  }
}

rhythmic.ld <- genes[circacompare.ld.ll[,1] < 0.05]
rhythmic.ld.ll <- intersect(rhythmic.ld, genes[circacompare.ld.ll[,2] < 0.05])
length(rhythmic.ld)
length(rhythmic.ld.ll)

## LD DD
circacompare.ld.dd <- matrix(nrow=number.genes,ncol=15)
rownames(circacompare.ld.dd) <- genes

time.points.ld <- seq(from=0,by=4,length.out = 18)
time.points.dd <- seq(from=0,by=4,length.out = 12)

for(i in 1:length(genes))
{
  print(i)
  gene.i <- genes[i]
  
  ld.expression.i <- gene.expression[gene.i,c(paste(ld.zt,1,sep="_"),paste(ld.zt,2,sep="_"),paste(ld.zt,3,sep="_"))]
  dd.expression.i <- gene.expression[gene.i,c(paste(ld.zt,6,sep="_"),paste(ld.zt,7,sep="_"))]
  
  ld.dd.df <- data.frame(time=c(time.points.ld,time.points.dd),
                         measure=c(ld.expression.i, dd.expression.i),
                         group=c(rep("ld",18),rep("dd",12)))
  
  out.i <- circacompare(x = ld.dd.df, col_time = "time", col_group = "group", col_outcome = "measure",alpha_threshold = 1)
  
  if(!is.null(out.i))
  {
    circacompare.ld.dd[i,] <- out.i[[2]][,2]
  }
}

rhythmic.ld <- genes[circacompare.ld.dd[,2] < 0.05]
rhythmic.ld.dd <- intersect(rhythmic.ld, genes[circacompare.ld.dd[,1] < 0.05])
length(rhythmic.ld)
length(rhythmic.ld.dd)

rhythmic.ld.ll.dd <- intersect(rhythmic.ld.ll, rhythmic.ld.dd)
length(rhythmic.ld.ll.dd)

## SD LL
sd.zt <- paste("sd",paste0("zt",sprintf(fmt = "%02d",seq(from=0,to=20,by=4))),sep="_")

circacompare.sd.ll <- matrix(nrow=number.genes,ncol=15)
rownames(circacompare.ld.ll) <- genes

time.points.ld <- seq(from=0,by=4,length.out = 18)
time.points.ll <- seq(from=0,by=4,length.out = 12)

for(i in 1:length(genes))
{
  print(i)
  gene.i <- genes[i]
  
  sd.expression.i <- gene.expression[gene.i,c(paste(sd.zt,1,sep="_"),paste(sd.zt,2,sep="_"),paste(sd.zt,3,sep="_"))]
  ll.expression.i <- gene.expression[gene.i,c(paste(sd.zt,4,sep="_"),paste(sd.zt,5,sep="_"))]
  
  sd.ll.df <- data.frame(time=c(time.points.ld,time.points.ll),
                         measure=c(sd.expression.i, ll.expression.i),
                         group=c(rep("ld",18),rep("ll",12)))
  
  out.i <- circacompare(x = sd.ll.df, col_time = "time", col_group = "group", col_outcome = "measure",alpha_threshold = 1)
  
  if(!is.null(out.i))
  {
    circacompare.sd.ll[i,] <- out.i[[2]][,2]
  }
}

rhythmic.sd <- genes[circacompare.sd.ll[,1] < 0.05]
rhythmic.sd.ll <- intersect(rhythmic.sd, genes[circacompare.sd.ll[,2] < 0.05])
length(rhythmic.sd)
length(rhythmic.sd.ll)

## SD DD
circacompare.sd.dd <- matrix(nrow=number.genes,ncol=15)
rownames(circacompare.sd.dd) <- genes

time.points.sd <- seq(from=0,by=4,length.out = 18)
time.points.dd <- seq(from=0,by=4,length.out = 12)

for(i in 1:length(genes))
{
  print(i)
  gene.i <- genes[i]
  
  sd.expression.i <- gene.expression[gene.i,c(paste(sd.zt,1,sep="_"),paste(sd.zt,2,sep="_"),paste(sd.zt,3,sep="_"))]
  dd.expression.i <- gene.expression[gene.i,c(paste(sd.zt,6,sep="_"),paste(sd.zt,7,sep="_"))]
  
  sd.dd.df <- data.frame(time=c(time.points.sd,time.points.dd),
                         measure=c(sd.expression.i, dd.expression.i),
                         group=c(rep("sd",18),rep("dd",12)))
  
  out.i <- circacompare(x = sd.dd.df, col_time = "time", col_group = "group", col_outcome = "measure",alpha_threshold = 1)
  
  if(!is.null(out.i))
  {
    circacompare.sd.dd[i,] <- out.i[[2]][,2]
  }
}

rhythmic.sd <- genes[circacompare.sd.dd[,2] < 0.05]
rhythmic.sd.dd <- intersect(rhythmic.sd, genes[circacompare.sd.dd[,1] < 0.05])
length(rhythmic.sd)
length(rhythmic.sd.dd)

rhythmic.sd.ll.dd <- intersect(rhythmic.sd.ll, rhythmic.sd.dd)
length(rhythmic.sd.ll.dd)

rhythmic.sd.ld.ll.dd <- intersect(rhythmic.sd.ll.dd, rhythmic.ld.ll.dd)
length(rhythmic.sd.ld.ll.dd)
length(rhythmic.sd.ld.ll.dd)/number.genes

length(intersect(rhythmic.sd.ll,rhythmic.ld.ll))

write(x = genes.two.peaks.sd.one.peak.ll.dd,file = "genes_two_peaks_sd_one_peak_ll_dd.txt")

## Two peaks become one single peak
i <- 57
gene.i <- genes.two.peaks.sd.one.peak.ll.dd[i]
#gene.i <- "ostta04g05380"
#gene.i <- "ostta01g03910"
#gene.i <- "ostta04g05380"
#gene.i <- "ostta05g02930"
#gene.i <- "ostta05g04830"
#gene.i <- "ostta08g00330"
#gene.i <- "ostta09g01760"
#gene.i <- "ostta17g02120"
#gene.i <- "ostta06g02940"
sd.ll.expression.i <- gene.expression[gene.i,c(paste(sd.zt,4,sep="_"),
                                               paste(sd.zt,5,sep="_"))]
sd.dd.expression.i <- gene.expression[gene.i,c(paste(sd.zt,6,sep="_"),
                                               paste(sd.zt,7,sep="_"))]

time.points <- seq(from=0,by=4,length.out = 12)

sd.ll.dd.df <- data.frame(time=c(time.points,time.points),
                          measure=c(sd.ll.expression.i, sd.dd.expression.i),
                          group=c(rep("ll",12),rep("dd",12)))

out.i <- circacompare(x = sd.ll.dd.df, col_time = "time", 
                      col_group = "group", col_outcome = "measure",alpha_threshold = 1)


current.gene <- gene.i#sd.dark.repressed[i]
#current.gene <- "ostta04g05380"
#current.gene <- "ostta01g03910"
#current.gene <- "ostta04g05380"
#current.gene <- "ostta05g02930"
#current.gene <- "ostta05g04830"
#current.gene <- "ostta08g00330"
#current.gene <- "ostta09g01760"
#current.gene <- "ostta17g02120"
#current.gene <- "ostta06g02940"
plot.sd.ll(gene.id = current.gene,gene.name = current.gene, 
           gene.expression = gene.expression)
lines(c(out.i$fit$m$fitted()[1:6],
        out.i$fit$m$fitted()[1:12],
        out.i$fit$m$fitted()[1:12]),lty=2,lwd=3,col="salmon")

plot.sd.dd(gene.id = current.gene,gene.name = current.gene, 
           gene.expression = gene.expression)
lines(c(out.i$fit$m$fitted()[13:18],
        out.i$fit$m$fitted()[13:24],
        out.i$fit$m$fitted()[13:24]),lty=2,lwd=3,col="brown3")


plot.sd(gene.id = current.gene,gene.name = current.gene,gene.expression)
lines(c(out.i$fit$m$fitted()[1:6],
        out.i$fit$m$fitted()[1:12]),lty=2,lwd=3,col="salmon")
lines(c(out.i$fit$m$fitted()[13:18],
        out.i$fit$m$fitted()[13:24]),lty=2,lwd=3,col="brown3")

# lines(c(out.i$fit$m$fitted()[1:6],
#         out.i$fit$m$fitted()[1:12])+c(out.i$fit$m$fitted()[13:18],
#         out.i$fit$m$fitted()[13:24]),lty=2,lwd=3,col="magenta")



ld.ll.expression.i <- gene.expression[gene.i,c(paste(ld.zt,4,sep="_"),
                                               paste(ld.zt,5,sep="_"))]
ld.dd.expression.i <- gene.expression[gene.i,c(paste(ld.zt,6,sep="_"),
                                               paste(ld.zt,7,sep="_"))]

time.points <- seq(from=0,by=4,length.out = 12)

ld.ll.dd.df <- data.frame(time=c(time.points,time.points),
                          measure=c(ld.ll.expression.i, ld.dd.expression.i),
                          group=c(rep("ll",12),rep("dd",12)))

out.i <- circacompare(x = ld.ll.dd.df, col_time = "time", 
                      col_group = "group", col_outcome = "measure",alpha_threshold = 1)

current.gene <- gene.i#sd.dark.repressed[i]
plot.ld.ll(gene.id = current.gene,gene.name = current.gene, gene.expression = gene.expression)
lines(c(out.i$fit$m$fitted()[1:6],
        out.i$fit$m$fitted()[1:12],
        out.i$fit$m$fitted()[1:12]),lty=2,lwd=3,col="lightblue")


plot.ld.dd(gene.id = current.gene,gene.name = current.gene, gene.expression = gene.expression)
lines(c(out.i$fit$m$fitted()[13:18],
        out.i$fit$m$fitted()[13:24],
        out.i$fit$m$fitted()[13:24]),lty=2,lwd=3,col="blue4")


plot.ld(gene.id = current.gene,gene.name = current.gene,gene.expression)
lines(c(out.i$fit$m$fitted()[1:6],
        out.i$fit$m$fitted()[1:12]),lty=2,lwd=3,col="lightblue")
lines(c(out.i$fit$m$fitted()[13:18],
        out.i$fit$m$fitted()[13:24]),lty=2,lwd=3,col="blue4")

# lines(c(out.i$fit$m$fitted()[1:6],
#         out.i$fit$m$fitted()[1:12])+c(out.i$fit$m$fitted()[13:18],
#         out.i$fit$m$fitted()[13:24]),lty=2,lwd=3,col="black")



plot(x=1,y=1)
i <- i + 1

#ostta01g01410
#ostta01g02230
#ostta01g02910
#ostta01g03450
#ostta01g03910 se ve bien ld
#ostta03g01980
#ostta03g02050
#ostta03g03980
#ostta03g04470
#ostta04g00430
#ostta04g00580
#ostta04g01910
#ostta04g02530
#ostta04g02640
#ostta04g02740
#ostta04g03390
#ostta04g04000
#ostta04g03730
#ostta04g04350
#ostta04g04900
#ostta04g05380 -----------
#ostta05g01610
#ostta05g01730
#ostta05g02930 ---------
#ostta05g03730
#ostta05g04830 ---------
#ostta06g01570
#ostta07g00740
#ostta07g02570
#ostta07g04510
#ostta08g00330 ---------- PSI2 phytochrome signalling 2
#ostta09g01760 --
#ostta10g01980
#ostta10g03130
#ostta10g03480
#ostta12g01690
#ostta14g01350
#ostta14g01470
#ostta15g00340
#ostta17g02120 --------

two.peaks.genes <- read.table(file = "two_peaks.txt",header = F,as.is = T)[[1]]

subset(gene.annotation, geneID %in% two.peaks.genes)

ostta2atha <- read.table(file = "annotation/ostta2atha.tsv",header = F)

write(x = subset(ostta2atha, V1 %in% two.peaks.genes)$V2,file = "atha.txt")


i <- 11

which(genes.two.peaks.sd.one.peak.ll.dd == "ostta05g01730")


## Dynamic change from LD to SD
wave <- function(k, alpha, phi, x)
{
  res <- (k  + alpha * cos((2*pi/24) * (x - phi))) 
  res[res < 0] <- 0
  return(res)
}

gene.id <- "ostta08g00710"
gene.id <- "ostta13g02440"
gene.id <- "ostta10g00920"
gene.id <- "ostta06g02940"
ind <- 1
gene.id <- complete.ld.rhythmic.genes[ind]

ld.zt <- paste("ld",paste0("zt",sprintf(fmt = "%02d",seq(from=0,to=20,by=4))),sep="_")
current.gene.expression.ld <- gene.expression[gene.id,c(paste(ld.zt,1,sep="_"),paste(ld.zt,2,sep="_"),paste(ld.zt,3,sep="_"))]

sd.zt <- paste("sd",paste0("zt",sprintf(fmt = "%02d",seq(from=0,to=20,by=4))),sep="_")
current.gene.expression.sd <- gene.expression[gene.id,c(paste(sd.zt,1,sep="_"),paste(sd.zt,2,sep="_"),paste(sd.zt,3,sep="_"))]


max.expr <- ceiling(max(c(unlist(current.gene.expression.ld), unlist(current.gene.expression.sd))))
min.expr <- floor(min(c(unlist(current.gene.expression.ld), unlist(current.gene.expression.sd))))
range.expr <- max.expr - min.expr

width.rectangule <- floor(range.expr / 10)

time.points <- seq(from=0,by=4,length.out = 18)

ld.sd.df <- data.frame(time=c(time.points,time.points),
                        measure=c(unlist(current.gene.expression.ld), 
                                  unlist(current.gene.expression.sd)),
                        group=c(rep("ld",18),rep("sd",18)))

out.i <- circacompare(x = ld.sd.df, col_time = "time", 
                      col_group = "group", col_outcome = "measure",alpha_threshold = 1)

circa.values <- out.i$summary$value
names(circa.values) <- out.i$summary$parameter

mesor.ld <- circa.values["ld mesor estimate"]
mesor.sd <- circa.values["sd mesor estimate"]

amplitude.ld <- circa.values["ld amplitude estimate"]
amplitude.sd <- circa.values["sd amplitude estimate"]
  
phase.ld <- circa.values["ld peak time hours"]
phase.sd <- circa.values["sd peak time hours"]

transitions.sd.ld <- matrix(nrow = 48,ncol=18)
mesors <- seq(from=mesor.sd,to=mesor.ld,length.out = 48)
amplitudes <- seq(from=amplitude.sd,to=amplitude.ld,length.out = 48)

if(phase.sd <= phase.ld)
{
  phases <- seq(from=phase.sd,to=phase.ld,length.out = 48)
} else
{
  until.dawn <- 24 - phase.sd
  
  phases <- seq(from=0,to=phase.ld+until.dawn,length.out = 48)-until.dawn
  phases[phases < 0] <- phases[phases < 0] + 24
}

for(i in 1:48)
{
  transitions.sd.ld[i,] <- wave(k = mesors[i], alpha = amplitudes[i], 
                                phi = phases[i],x = time.points)   
}

months <- c("January", "February", "March", "April","May", "June", "July", 
            "August", "September", "October", "November", "December")

months <- rep(months,each=8)
colors <- colorRampPalette(c("red","blue"))(48)
dusk <- seq(from=8,to=16,length.out = 48)
k <- 1
for(j in 1:3)
{
  for(i in 1:48)
  {
    plot(time.points,current.gene.expression.ld,type="o",col="blue",lwd=0.25,pch=19,cex=0.5,
         main=paste(c(months[k], " (Day ",floor(dusk[i]), " hours : Night ",
                      24-floor(dusk[i]), " hours)"),collapse=""),
         ylim=c(min.expr - 2 * width.rectangule,max.expr),axes=F,xlab="",ylab="")
    lines(time.points,current.gene.expression.sd,type="o",col="red",lwd=0.25,pch=19,cex=0.5)
    
    for(l in 0:2)
    {
      polygon(x = c(24*l, 24*l+dusk[i], 24*l+dusk[i], 24*l),
              y=c(min.expr-1.5*width.rectangule, min.expr-1.5*width.rectangule, 
                  min.expr-0.5*width.rectangule, min.expr-0.5*width.rectangule),
              lwd=2,border=colors[i])
      polygon(x = c(24*l+dusk[i],24*(l+1),24*(l+1),24*l+dusk[i]),
              y=c(min.expr-1.5*width.rectangule, min.expr-1.5*width.rectangule, 
                  min.expr-0.5*width.rectangule, min.expr-0.5*width.rectangule),
              lwd=2,border=colors[i],col=colors[i])
    }
    
    lines(time.points,transitions.sd.ld[i,],type="l",col=colors[i],lwd=2)
    k <- k + 1
    Sys.sleep(0.25)
    print(i)
  }
  
  for(i in 48:1)
  {
    # plot(time.points,current.gene.expression.ld,type="o",col="blue",lwd=0.25,pch=19,cex=0.5,
    #      main=months[k],ylim=c(min.expr - 2 * width.rectangule,max.expr))
    # lines(time.points,current.gene.expression.sd,type="o",col="red",lwd=0.25,pch=19,cex=0.5)
    
    plot(time.points,current.gene.expression.ld,type="o",col="blue",lwd=0.25,pch=19,cex=0.5,
         main=paste(c(months[k], " (Day ",floor(dusk[i]), " hours : Night ",
                      24-floor(dusk[i]), " hours)"),collapse=""),
         ylim=c(min.expr - 2 * width.rectangule,max.expr),axes=F,xlab="",ylab="")
    lines(time.points,current.gene.expression.sd,type="o",col="red",lwd=0.25,pch=19,cex=0.5)
    
    
    for(l in 0:2)
    {
      polygon(x = c(24*l, 24*l+dusk[i], 24*l+dusk[i], 24*l),
              y=c(min.expr-1.5*width.rectangule, min.expr-1.5*width.rectangule, 
                  min.expr-0.5*width.rectangule, min.expr-0.5*width.rectangule),
              lwd=2,border=colors[i])
      polygon(x = c(24*l+dusk[i],24*(l+1),24*(l+1),24*l+dusk[i]),
              y=c(min.expr-1.5*width.rectangule, min.expr-1.5*width.rectangule, 
                  min.expr-0.5*width.rectangule, min.expr-0.5*width.rectangule),
              lwd=2,border=colors[i],col=colors[i])
    }
    
    lines(time.points,transitions.sd.ld[i,],type="l",col=colors[i],lwd=2)
    k <- k + 1
    Sys.sleep(0.25)
    print(i)
  }
  k <- 1
}

ind <- ind + 1


gene.id <- "ostta04g05380"
gene.id <- "ostta08g00330"
gene.id <- "ostta17g02120"
gene.id <- "ostta05g04830"
gene.id <- "ostta05g02930"
gene.id <- "ostta04g05380"
gene.id <- "ostta04g02740" #
gene.id <- "ostta04g04000"

gene.id <- "ostta06g02940"
#gene.id <- complete.ld.rhythmic.genes[ind]

ld.zt <- paste("ld",paste0("zt",sprintf(fmt = "%02d",seq(from=0,to=20,by=4))),sep="_")
current.gene.expression.ld <- gene.expression[gene.id,c(paste(ld.zt,1,sep="_"),paste(ld.zt,2,sep="_"),paste(ld.zt,3,sep="_"))]
current.gene.expression.ld.ll <- gene.expression[gene.id,c(paste(ld.zt,4,sep="_"),paste(ld.zt,5,sep="_"))]
current.gene.expression.ld.dd <- gene.expression[gene.id,c(paste(ld.zt,6,sep="_"),paste(ld.zt,7,sep="_"))]

sd.zt <- paste("sd",paste0("zt",sprintf(fmt = "%02d",seq(from=0,to=20,by=4))),sep="_")
current.gene.expression.sd <- gene.expression[gene.id,c(paste(sd.zt,1,sep="_"),paste(sd.zt,2,sep="_"),paste(sd.zt,3,sep="_"))]
current.gene.expression.sd.ll <- gene.expression[gene.id,c(paste(sd.zt,4,sep="_"),paste(sd.zt,5,sep="_"))]
current.gene.expression.sd.dd <- gene.expression[gene.id,c(paste(sd.zt,6,sep="_"),paste(sd.zt,7,sep="_"))]

time.points <- seq(from=0,by=4,length.out = 12)

sd.ll.dd.df <- data.frame(time=c(time.points,time.points),
                          measure=c(current.gene.expression.sd.ll, 
                                    current.gene.expression.sd.dd),
                          group=c(rep("sd_ll",12),rep("sd_dd",12)))

sd.out.i <- circacompare(x = sd.ll.dd.df, col_time = "time", 
                         col_group = "group", 
                         col_outcome = "measure",alpha_threshold = 1)
circa.values.sd <- sd.out.i$summary$value
names(circa.values.sd) <- sd.out.i$summary$parameter
  
mesor.sd.ll <- circa.values.sd["sd_ll mesor estimate"]
mesor.sd.dd <- circa.values.sd["sd_dd mesor estimate"]

amplitude.sd.ll <- circa.values.sd["sd_ll amplitude estimate"]
amplitude.sd.dd <- circa.values.sd["sd_dd amplitude estimate"]

phase.sd.ll <- circa.values.sd["sd_ll peak time hours"]
phase.sd.dd <- circa.values.sd["sd_dd peak time hours"]

wave.2 <- function(beta, phi, x)
{
 res <- (1  + beta * cos((2*pi/24) * (x - phi))) 
 res[res < 0] <- 0
 return(res)
}

two.waves.combination.2 <- function(x,beta1,beta2,phi1,phi2)
{
 wave1 <- wave.2(beta1,phi1,x)
 wave2 <- wave.2(beta2,phi2,x)
 two.waves <- matrix(data = c(wave1,wave2),ncol=2)
 res <- apply(X = two.waves,MARGIN = 1,FUN = max)
 
 return(res)
}

two.waves.combination.3 <- function(x,beta1,beta2,phi,offset)
{
 phi1 <- phi
 phi2 <- (phi + offset) %% 24
 wave1 <- wave.2(beta1,phi1,x)
 wave2 <- wave.2(beta2,phi2,x)
 two.waves <- matrix(data = c(wave1,wave2),ncol=2)
 res <- apply(X = two.waves,MARGIN = 1,FUN = max)
 
 return(res)
}



time.points <- seq(from=0,by=4,length.out = 18)
gene.expression.df.2 <- data.frame(x=time.points,y=current.gene.expression.sd/mean(current.gene.expression.sd))



m <- nls( y ~ two.waves.combination.2(x,beta1,beta2,phi1,phi2), data = gene.expression.df.2,
          start = list(beta1=amplitude.sd.ll/mesor.sd.ll,
                       beta2=amplitude.sd.dd/mesor.sd.dd,
                       phi1=phase.sd.ll,
                       phi2=phase.sd.dd),
          lower = c(0.1,0.1,0,0), upper = c(2,2,24,24),algorithm = "port", 
          control=stats::nls.control(warnOnly = T),
          #control=stats::nls.control(maxiter = 100, minFactor = 1/10000,warnOnly = T),
          trace = T)

summariy.m <- summary(m)
estimates <- summariy.m$coefficients[,"Estimate"]
beta1.estimate <- estimates[["beta1"]]
beta2.estimate <- estimates[["beta2"]]
phi1.estimate <- estimates[["phi1"]]
phi2.estimate <- estimates[["phi2"]]

dos.ondas.2 <- two.waves.combination.2(x=seq(from=0,to=68,by=0.01),
                                       beta1=beta1.estimate,
                                       beta2=beta2.estimate,
                                       phi1=phi1.estimate,
                                       phi2=phi2.estimate)


plot(time.points,current.gene.expression.sd/mean(current.gene.expression.sd),type="o",col="red",lwd=0.25,pch=19,cex=0.5)#,
lines(seq(from=0,to=68,by=0.01),dos.ondas.2,type="l",col="blue",lwd=2,pch=19,cex=0.5)


wave1 <- wave.2(beta=beta1.estimate,phi=phi1.estimate,x=seq(from=0,to=68,by=0.01))
wave2 <- wave.2(beta=beta2.estimate,phi=phi2.estimate,x=seq(from=0,to=68,by=0.01))

time.points <- seq(from=0,by=4,length.out = 18)
plot(time.points,current.gene.expression.sd,type="o",col="red",lwd=0.25,pch=19,cex=0.5)#,
lines(seq(from=0,to=68,by=0.01),mean(current.gene.expression.sd)*dos.ondas.2,type="l",col="blue",lwd=2,pch=19,cex=0.5)

lines(seq(from=0,to=68,by=0.01),mean(current.gene.expression.sd)*wave1,type="l",col="blue",lwd=2,pch=19,cex=0.5,lty=2)
lines(seq(from=0,to=68,by=0.01),mean(current.gene.expression.sd)*wave2,type="l",col="blue",lwd=2,pch=19,cex=0.5,lty=2)


m.2 <- nls( y ~ two.waves.combination.3(x,beta1,beta2,phi,offset), data = gene.expression.df.2,
          start = list(beta1=amplitude.sd.dd/mesor.sd.dd,
                       beta2=amplitude.sd.ll/mesor.sd.ll,
                       phi=phase.sd.dd,
                       offset=11),#phase.sd.ll-phase.sd.dd),
          lower = c(0.1,0.1,0,11), upper = c(2,2,24,13),algorithm = "port", 
          control=stats::nls.control(maxiter = 100, minFactor = 1/10000,warnOnly = T),
          trace = T)


summariy.m.2 <- summary(m.2)
estimates.2 <- summariy.m.2$coefficients[,"Estimate"]
beta1.estimate.2 <- estimates.2[["beta1"]]
beta2.estimate.2 <- estimates.2[["beta2"]]
phi1.estimate.2 <- estimates.2[["phi"]]
phi2.estimate.2 <- ((phi1.estimate.2 + estimates.2[["offset"]]) %% 24)

dos.ondas.2 <- two.waves.combination.2(x=seq(from=0,to=68,by=0.01),
                                       beta1=beta1.estimate.2,
                                       beta2=beta2.estimate.2,
                                       phi1=phi1.estimate.2,
                                       phi2=phi2.estimate.2)


plot(time.points,current.gene.expression.sd/mean(current.gene.expression.sd),type="o",col="red",lwd=0.25,pch=19,cex=0.5)#,
lines(seq(from=0,to=68,by=0.01),dos.ondas.2,type="l",col="blue",lwd=2,pch=19,cex=0.5)


wave1 <- wave.2(beta=beta1.estimate.2,phi=phi1.estimate.2,x=seq(from=0,to=68,by=0.01))
wave2 <- wave.2(beta=beta2.estimate.2,phi=phi2.estimate.2,x=seq(from=0,to=68,by=0.01))

time.points <- seq(from=0,by=4,length.out = 18)
plot(time.points,current.gene.expression.sd,type="o",col="red",lwd=0.25,pch=19,cex=0.5)#,
lines(seq(from=0,to=68,by=0.01),mean(current.gene.expression.sd)*dos.ondas.2,type="l",col="red",lwd=2,pch=19,cex=0.5)

lines(seq(from=0,to=68,by=0.01),mean(current.gene.expression.sd)*wave1,type="l",col="red",lwd=2,pch=19,cex=0.5,lty=2)
lines(seq(from=0,to=68,by=0.01),mean(current.gene.expression.sd)*wave2,ltype="l",col="red",lwd=2,pch=19,cex=0.5,lty=2)

plot(time.points,current.gene.expression.ld,type="o",col="blue",lwd=0.25,pch=19,cex=0.5)#,







time.points <- seq(from=0,by=4,length.out = 18)
wave.sd.ll <- wave(k = mesor.sd.ll, alpha = amplitude.sd.ll, 
     phi = phase.sd.ll+9,x = seq(from=0,to=68,by=0.01))
wave.sd.dd <- wave(k = mesor.sd.dd, alpha = amplitude.sd.dd, 
                   phi = phase.sd.dd,x = seq(from=0,to=68,by=0.01))

combine.waves <- function(wave1,wave2)
{
  final.wave <- wave1
  final.wave[wave1 < wave2] <- wave2[wave1 < wave2]
  return(final.wave)
}

wave.sd.ll.dd <- combine.waves(wave1 = wave.sd.ll, wave2 = wave.sd.dd)

plot(time.points,current.gene.expression.sd,type="o",col="red",lwd=0.25,pch=19,cex=0.5)#,
     # main=paste(c(months[k], " (Day ",floor(dusk[i]), " hours : Night ",
     #              24-floor(dusk[i]), " hours)"),collapse=""),
     # ylim=c(min.expr - 2 * width.rectangule,max.expr),axes=F,xlab="",ylab="")
lines(seq(from=0,to=68,by=0.01),wave.sd.ll,type="l",col="red",lwd=0.25,pch=19,cex=0.5,lty=2)
lines(seq(from=0,to=68,by=0.01),wave.sd.dd,type="l",col="red",lwd=0.25,pch=19,cex=0.5,lty=2)
lines(seq(from=0,to=68,by=0.01),wave.sd.ll.dd,type="l",col="red",lwd=2,pch=19,cex=0.5)




max.expr <- ceiling(max(c(unlist(current.gene.expression.ld), unlist(current.gene.expression.sd))))
min.expr <- floor(min(c(unlist(current.gene.expression.ld), unlist(current.gene.expression.sd))))
range.expr <- max.expr - min.expr

width.rectangule <- floor(range.expr / 10)

time.points <- seq(from=0,by=4,length.out = 18)

ld.sd.df <- data.frame(time=c(time.points,time.points),
                       measure=c(unlist(current.gene.expression.ld), 
                                 unlist(current.gene.expression.sd)),
                       group=c(rep("ld",18),rep("sd",18)))

wave <- function(k, alpha, phi, x)
{
  res <- (k  + alpha * cos((2*pi/24) * (x - phi))) 
  res[res < 0] <- 0
  return(res)
}

two.waves.combination <- function(x,k1,k2,alpha1,alpha2,phi1,phi2)
{
  wave1 <- wave(k1,alpha1,phi1,x)
  wave2 <- wave(k2,alpha2,phi2,x)
  two.waves <- matrix(data = c(wave1,wave2),ncol=2)
  res <- apply(X = two.waves,MARGIN = 1,FUN = max)
  
  
  # final.wave <- wave1
  # final.wave[wave1 < wave2] <- wave2[wave1 < wave2]
  # return(final.wave)
  
  return(res)
}







gene.expression.df <- data.frame(x=time.points,y=current.gene.expression.sd)

m <- nls( y ~ two.waves.combination(x,k1,k2,alpha1,alpha2,phi1,phi2), data = gene.expression.df,
          start = list(k1=mesor.sd.ll,k2=mesor.sd.dd,alpha1=amplitude.sd.ll,alpha2=amplitude.sd.dd,phi1=phase.sd.ll+6,phi2=phase.sd.dd),trace = T)

dos.ondas <- two.waves.combination(x=seq(from=0,to=68,by=0.01),k1=mesor.sd.ll,k2=mesor.sd.dd,alpha1=amplitude.sd.ll,alpha2=amplitude.sd.dd,phi1=phase.sd.ll+6,phi2=phase.sd.dd)
lines(seq(from=0,to=68,by=0.01),dos.ondas,type="l",col="black",lwd=2,pch=19,cex=0.5)


wave.sd.ll <- wave(k = mesor.sd.ll, alpha = amplitude.sd.ll, 
                   phi = phase.sd.ll-6,x = seq(from=0,to=68,by=0.01))
wave.sd.dd <- wave(k = mesor.sd.dd, alpha = amplitude.sd.dd, 
                   phi = phase.sd.dd,x = seq(from=0,to=68,by=0.01))



wave.2 <- function(beta, phi, x)
{
  res <- (1  + beta * cos((2*pi/24) * (x - phi))) 
  res[res < 0] <- 0
  return(res)
}

two.waves.combination.2 <- function(x,beta1,beta2,phi1,phi2)
{
  wave1 <- wave.2(beta1,phi1,x)
  wave2 <- wave.2(beta2,phi2,x)
  two.waves <- matrix(data = c(wave1,wave2),ncol=2)
  res <- apply(X = two.waves,MARGIN = 1,FUN = max)

  return(res)
}

gene.expression.df.2 <- data.frame(x=time.points,y=current.gene.expression.sd/mean(current.gene.expression.sd))

warnOnly = TRUE
m <- nls( y ~ two.waves.combination.2(x,beta1,beta2,phi1,phi2), data = gene.expression.df.2,
          start = list(beta1=amplitude.sd.ll/mesor.sd.ll,
                       beta2=amplitude.sd.dd/mesor.sd.dd,
                       phi1=phase.sd.ll+6,
                       phi2=phase.sd.dd),
          lower = c(0.1,0.1,0,0), upper = c(0.9,0.9,24,24),algorithm = "port", 
          control=stats::nls.control(maxiter = 100, minFactor = 1/10000,warnOnly = T),
          trace = T)

summariy.m <- summary(m)
estimates <- summariy.m$coefficients[,"Estimate"]
beta1.estimate <- estimates[["beta1"]]
beta2.estimate <- estimates[["beta2"]]
phi1.estimate <- estimates[["phi1"]]
phi2.estimate <- estimates[["phi2"]]

dos.ondas.2 <- two.waves.combination.2(x=seq(from=0,to=68,by=0.01),
                                       beta1=beta1.estimate,
                                       beta2=beta2.estimate,
                                       phi1=phi1.estimate,
                                       phi2=phi2.estimate)


plot(time.points,current.gene.expression.sd/mean(current.gene.expression.sd),type="o",col="red",lwd=0.25,pch=19,cex=0.5)#,
lines(seq(from=0,to=68,by=0.01),dos.ondas.2,type="l",col="blue",lwd=2,pch=19,cex=0.5)


wave1 <- wave.2(beta=beta1.estimate,phi=phi1.estimate,x=seq(from=0,to=68,by=0.01))
wave2 <- wave.2(beta=beta2.estimate,phi=phi2.estimate,x=seq(from=0,to=68,by=0.01))

plot(time.points,current.gene.expression.sd,type="o",col="red",lwd=0.25,pch=19,cex=0.5)#,
lines(seq(from=0,to=68,by=0.01),mean(current.gene.expression.sd)*dos.ondas.2,type="l",col="blue",lwd=2,pch=19,cex=0.5)

lines(seq(from=0,to=68,by=0.01),mean(current.gene.expression.sd)*wave1,ltype="l",col="blue",lwd=2,pch=19,cex=0.5,lty=2)
lines(seq(from=0,to=68,by=0.01),mean(current.gene.expression.sd)*wave2,ltype="l",col="blue",lwd=2,pch=19,cex=0.5,lty=2)







m <- nls( y ~ two.waves.combination.2(x,beta1,beta2,phi1,phi2), data = gene.expression.df.2,
          start = list(beta1=0.303490,beta2=0.620566,
                       phi1=19.5612,phi2=11.1297),
          lower = c(0.3,0.3,0,0), upper = c(1,1,24,24),algorithm = "port", 
          control=stats::nls.control(maxiter = 100, minFactor = 1/10000),
          trace = T)


dos.ondas.2 <- two.waves.combination.2(x=seq(from=0,to=68,by=0.01),
                                       beta1=0.213682,
                                       beta2=0.618108,phi1=19.3875,phi2=10.8411)

dos.ondas.2 <- two.waves.combination.2(x=seq(from=0,to=68,by=0.01),
                                       beta1=0.303490,
                                       beta2=0.620567,phi1=19.5612,phi2=11.1297)


plot(time.points,current.gene.expression.sd/mean(current.gene.expression.sd),type="o",col="red",lwd=0.25,pch=19,cex=0.5)#,
lines(seq(from=0,to=68,by=0.01),dos.ondas.2,type="l",col="blue",lwd=2,pch=19,cex=0.5)


library(nlstools)
library(minpack.lm)

m <- nlsLM( y ~ two.waves.combination.2(x,beta1,beta2,phi1,phi2), data = gene.expression.df.2,
            start = list(beta1=amplitude.sd.ll/mesor.sd.ll,beta2=amplitude.sd.dd/mesor.sd.dd,
                         phi1=phase.sd.ll,phi2=phase.sd.dd),
          lower = c(0.05,0.05,0,0), upper = c(1,1,24,24),
          algorithm = "port", control=stats::nls.control(maxiter = 10000, minFactor = 1/1000000),
          trace = T)

summary(m)
dos.ondas.2 <- two.waves.combination.2(x=seq(from=0,to=68,by=0.01),
                                       beta1=0.1700,
                                       beta2=0.4821,
                                       phi1=16.3252,
                                       phi2=11.3710)


plot(time.points,current.gene.expression.sd/mean(current.gene.expression.sd),type="o",col="red",lwd=0.25,pch=19,cex=0.5)#,
lines(seq(from=0,to=68,by=0.01),dos.ondas.2,type="l",col="blue",lwd=2,pch=19,cex=0.5)


gene.id <- "ostta05g02930"
sd.zt <- paste("sd",paste0("zt",sprintf(fmt = "%02d",seq(from=0,to=20,by=4))),sep="_")
current.gene.expression.sd <- gene.expression[gene.id,c(paste(sd.zt,1,sep="_"),paste(sd.zt,2,sep="_"),paste(sd.zt,3,sep="_"))]

gene.expression.df.2 <- data.frame(x=time.points,y=current.gene.expression.sd/mean(current.gene.expression.sd))



m <- nlsLM( y ~ two.waves.combination.2(x,beta1,beta2,phi1,phi2), data = gene.expression.df.2,
            start = list(beta1=0.128812+0.3*0.128812,beta2=0.243531+0.3*0.243531,
                         phi1=9.80799+0.3*9.80799,phi2=7.75815+0.3*7.75815),
            lower = c(0.05,0.05,0,0), upper = c(1,1,24,24),algorithm = "port", control=stats::nls.control(maxiter = 100, minFactor = 1/10000),
            trace = T)


#-----------------------------------------
## DOOS PICOOS
#-----------------------------------------

wave.2 <- function(beta, phi, x)
{
 res <- (1  + beta * cos((2*pi/24) * (x - phi))) 
 res[res < 0] <- 0
 return(res)
}

two.waves.combination.3 <- function(x,beta1,beta2,phi,offset)
{
 phi1 <- phi
 phi2 <- (phi + offset) %% 24
 wave1 <- wave.2(beta1,phi1,x)
 wave2 <- wave.2(beta2,phi2,x)
 two.waves <- matrix(data = c(wave1,wave2),ncol=2)
 res <- apply(X = two.waves,MARGIN = 1,FUN = max)
 
 return(res)
}

current.gene <- "ostta04g02740"
plot.ld.sd(gene.id=current.gene, gene.name=current.gene, gene.expression)
plot.ld.ll(gene.id=current.gene, gene.name=current.gene, gene.expression)
plot.ld.dd(gene.id=current.gene, gene.name=current.gene, gene.expression)
plot.sd.ll(gene.id=current.gene, gene.name=current.gene, gene.expression)
plot.sd.dd(gene.id=current.gene, gene.name=current.gene, gene.expression)



#ostta01g01700 doble pico
#ostta01g02130 doble pico
#ostta01g02270 doble pico
#ostta01g02210 doble pico
#ostta01g02390 doble pico
#ostta01g03060 doble pico
#ostta01g05010 doble pico
#ostta01g05220 doble pico

current.gene <- "ostta04g04000" ## PWD
current.gene <- "ostta06g02940" ## GBSS
current.gene <- "ostta01g01700"
current.gene <- "ostta01g02130"
current.gene <- "ostta01g02270"
current.gene <- "ostta01g02210"
current.gene <- "ostta01g02390"
current.gene <- "ostta01g03060"
current.gene <- "ostta01g05010"
plot.ld.sd(gene.id=current.gene, gene.name="OsttaGBSS", gene.expression)
plot.ld.ll(gene.id=current.gene, gene.name="OsttaGBSS", gene.expression)
plot.ld.dd(gene.id=current.gene, gene.name="OsttaGBSS", gene.expression)
plot.sd.ll(gene.id=current.gene, gene.name="OsttaGBSS", gene.expression)
plot.sd.dd(gene.id=current.gene, gene.name="OsttaGBSS", gene.expression)

gene.id <- current.gene

ld.zt <- paste("ld",paste0("zt",sprintf(fmt = "%02d",seq(from=0,to=20,by=4))),sep="_")
current.gene.expression.ld <- gene.expression[gene.id,c(paste(ld.zt,1,sep="_"),paste(ld.zt,2,sep="_"),paste(ld.zt,3,sep="_"))]
current.gene.expression.ld.ll <- gene.expression[gene.id,c(paste(ld.zt,4,sep="_"),paste(ld.zt,5,sep="_"))]
current.gene.expression.ld.dd <- gene.expression[gene.id,c(paste(ld.zt,6,sep="_"),paste(ld.zt,7,sep="_"))]

sd.zt <- paste("sd",paste0("zt",sprintf(fmt = "%02d",seq(from=0,to=20,by=4))),sep="_")
current.gene.expression.sd <- gene.expression[gene.id,c(paste(sd.zt,1,sep="_"),paste(sd.zt,2,sep="_"),paste(sd.zt,3,sep="_"))]
current.gene.expression.sd.ll <- gene.expression[gene.id,c(paste(sd.zt,4,sep="_"),paste(sd.zt,5,sep="_"))]
current.gene.expression.sd.dd <- gene.expression[gene.id,c(paste(sd.zt,6,sep="_"),paste(sd.zt,7,sep="_"))]

time.points <- seq(from=0,by=4,length.out = 12)

sd.ll.dd.df <- data.frame(time=c(time.points,time.points),
                          measure=c(current.gene.expression.sd.ll, 
                                    current.gene.expression.sd.dd),
                          group=c(rep("sd_ll",12),rep("sd_dd",12)))

sd.out.i <- circacompare(x = sd.ll.dd.df, col_time = "time", 
                         col_group = "group", 
                         col_outcome = "measure",alpha_threshold = 1)
circa.values.sd <- sd.out.i$summary$value
names(circa.values.sd) <- sd.out.i$summary$parameter

mesor.sd.ll <- circa.values.sd["sd_ll mesor estimate"]
mesor.sd.dd <- circa.values.sd["sd_dd mesor estimate"]

amplitude.sd.ll <- circa.values.sd["sd_ll amplitude estimate"]
amplitude.sd.dd <- circa.values.sd["sd_dd amplitude estimate"]

phase.sd.ll <- circa.values.sd["sd_ll peak time hours"]
phase.sd.dd <- circa.values.sd["sd_dd peak time hours"]

time.points <- seq(from=0,by=4,length.out = 18)
gene.expression.df.2 <- data.frame(x=time.points,
                                   y=current.gene.expression.sd/mean(current.gene.expression.sd))

m.2 <- nls( y ~ two.waves.combination.3(x,beta1,beta2,phi,offset), data = gene.expression.df.2,
            start = list(beta1=amplitude.sd.dd/mesor.sd.dd,
                         beta2=amplitude.sd.ll/mesor.sd.ll,
                         phi=phase.sd.dd,
                         offset=12),#phase.sd.ll-phase.sd.dd),
            lower = c(0.02,0.02,0,8), upper = c(2,2,24,16),algorithm = "port", 
            control=stats::nls.control(maxiter = 100, minFactor = 1/10000,warnOnly = T),
            trace = T)

summariy.m.2 <- summary(m.2)
estimates.2 <- summariy.m.2$coefficients[,"Estimate"]
beta1.estimate.2 <- estimates.2[["beta1"]]
beta2.estimate.2 <- estimates.2[["beta2"]]
phi1.estimate.2 <- estimates.2[["phi"]]
phi2.estimate.2 <- ((phi1.estimate.2 + estimates.2[["offset"]]) %% 24)


dos.ondas.2 <- two.waves.combination.2(x=seq(from=0,to=68,by=0.01),
                                       beta1=beta1.estimate.2,
                                       beta2=beta2.estimate.2,
                                       phi1=phi1.estimate.2,
                                       phi2=phi2.estimate.2)

wave1 <- wave.2(beta=beta1.estimate.2,phi=phi1.estimate.2,x=seq(from=0,to=68,by=0.01))
wave2 <- wave.2(beta=beta2.estimate.2,phi=phi2.estimate.2,x=seq(from=0,to=68,by=0.01))

time.points <- seq(from=0,by=4,length.out = 18)
plot(time.points,current.gene.expression.sd,type="o",col="red",lwd=0.25,pch=19,cex=0.5)#,
lines(seq(from=0,to=68,by=0.01),mean(current.gene.expression.sd)*dos.ondas.2,type="l",col="red",lwd=2,pch=19,cex=0.5)
lines(seq(from=0,to=68,by=0.01),mean(current.gene.expression.sd)*wave1,type="l",col="red",lwd=2,pch=19,cex=0.5,lty=2)
lines(seq(from=0,to=68,by=0.01),mean(current.gene.expression.sd)*wave2,type="l",col="red",lwd=2,pch=19,cex=0.5,lty=2)

ld.zt <- paste("ld",paste0("zt",sprintf(fmt = "%02d",seq(from=0,to=20,by=4))),sep="_")
current.gene.expression.ld <- gene.expression[gene.id,c(paste(ld.zt,1,sep="_"),paste(ld.zt,2,sep="_"),paste(ld.zt,3,sep="_"))]

sd.zt <- paste("sd",paste0("zt",sprintf(fmt = "%02d",seq(from=0,to=20,by=4))),sep="_")
current.gene.expression.sd <- gene.expression[gene.id,c(paste(sd.zt,1,sep="_"),paste(sd.zt,2,sep="_"),paste(sd.zt,3,sep="_"))]


max.expr <- ceiling(max(c(unlist(current.gene.expression.ld), unlist(current.gene.expression.sd))))
min.expr <- floor(min(c(unlist(current.gene.expression.ld), unlist(current.gene.expression.sd))))
range.expr <- max.expr - min.expr

width.rectangule <- floor(range.expr / 10)

time.points <- seq(from=0,by=4,length.out = 18)

ld.sd.df <- data.frame(time=c(time.points,time.points),
                       measure=c(unlist(current.gene.expression.ld), 
                                 unlist(current.gene.expression.sd)),
                       group=c(rep("ld",18),rep("sd",18)))

out.i <- circacompare(x = ld.sd.df, col_time = "time", 
                      col_group = "group", col_outcome = "measure",alpha_threshold = 1)

circa.values <- out.i$summary$value
names(circa.values) <- out.i$summary$parameter

mesor.ld <- circa.values["ld mesor estimate"]
amplitude.ld <- circa.values["ld amplitude estimate"]
phase.ld <- circa.values["ld peak time hours"]


mesor.sd <- circa.values["sd mesor estimate"]
amplitude.sd <- circa.values["sd amplitude estimate"]
phase.sd <- circa.values["sd peak time hours"]

phase.sd.1 <- phi1.estimate.2
phase.sd.2 <- phi2.estimate.2
# phase.sd.1 <- phi2.estimate.2
# phase.sd.2 <- phi1.estimate.2 


phase.sd <- phase.sd.2

transitions.sd.ld <- matrix(nrow = 48,ncol=18)
fixed.sd.ld <- matrix(nrow = 48,ncol=18)
mesors <- seq(from=mesor.sd,to=mesor.ld,length.out = 48)
amplitudes <- seq(from=amplitude.sd,to=amplitude.ld,length.out = 48)

if(phase.sd <= phase.ld)
{
 phases <- seq(from=phase.sd,to=phase.ld,length.out = 48)
} else
{
 until.dawn <- 24 - phase.sd
 
 phases <- seq(from=0,to=phase.ld+until.dawn,length.out = 48)-until.dawn
 phases[phases < 0] <- phases[phases < 0] + 24
}

for(i in 1:48)
{
 transitions.sd.ld[i,] <- wave(k = mesors[i], alpha = amplitudes[i], 
                               phi = phases[i],x = time.points)   
 fixed.sd.ld[i,] <- wave(k = mesor.sd, alpha = amplitude.sd, 
                         phi = phase.sd.1,x = time.points)
}

months <- c("January", "February", "March", "April","May", "June", "July", 
            "August", "September", "October", "November", "December")

months <- rep(months,each=8)
colors <- colorRampPalette(c("red","blue"))(48)
dusk <- seq(from=8,to=16,length.out = 48)
k <- 1
for(j in 1:2)#3)
{
 for(i in 1:48)
 {
  plot(time.points,current.gene.expression.ld,type="o",col="blue",lwd=0.25,pch=19,cex=0.5,
       main=paste(c(months[k], " (Day ",floor(dusk[i]), " hours : Night ",
                    24-floor(dusk[i]), " hours)"),collapse=""),
       ylim=c(min.expr - 2 * width.rectangule,max.expr),axes=F,xlab="",ylab="")
  lines(time.points,current.gene.expression.sd,type="o",col="red",lwd=0.25,pch=19,cex=0.5)
  
  for(l in 0:2)
  {
   polygon(x = c(24*l, 24*l+dusk[i], 24*l+dusk[i], 24*l),
           y=c(min.expr-1.5*width.rectangule, min.expr-1.5*width.rectangule, 
               min.expr-0.5*width.rectangule, min.expr-0.5*width.rectangule),
           lwd=2,border=colors[i])
   polygon(x = c(24*l+dusk[i],24*(l+1),24*(l+1),24*l+dusk[i]),
           y=c(min.expr-1.5*width.rectangule, min.expr-1.5*width.rectangule, 
               min.expr-0.5*width.rectangule, min.expr-0.5*width.rectangule),
           lwd=2,border=colors[i],col=colors[i])
  }
  
  lines(time.points,transitions.sd.ld[i,],type="l",col=colors[i],lwd=1.5,lty=2)
  lines(time.points,fixed.sd.ld[i,],type="l",col=colors[i],lwd=1.5,lty=3)
  combined.wave.i <- apply(X = matrix(data = c(transitions.sd.ld[i,],fixed.sd.ld[i,]),ncol=2),MARGIN = 1,FUN = max)
  lines(time.points,combined.wave.i,type="l",col=colors[i],lwd=2,lty=1)
  
  k <- k + 1
  Sys.sleep(0.25)
  print(i)
 }
 
 for(i in 48:1)
 {
  # plot(time.points,current.gene.expression.ld,type="o",col="blue",lwd=0.25,pch=19,cex=0.5,
  #      main=months[k],ylim=c(min.expr - 2 * width.rectangule,max.expr))
  # lines(time.points,current.gene.expression.sd,type="o",col="red",lwd=0.25,pch=19,cex=0.5)
  
  plot(time.points,current.gene.expression.ld,type="o",col="blue",lwd=0.25,pch=19,cex=0.5,
       main=paste(c(months[k], " (Day ",floor(dusk[i]), " hours : Night ",
                    24-floor(dusk[i]), " hours)"),collapse=""),
       ylim=c(min.expr - 2 * width.rectangule,max.expr),axes=F,xlab="",ylab="")
  lines(time.points,current.gene.expression.sd,type="o",col="red",lwd=0.25,pch=19,cex=0.5)
  
  
  for(l in 0:2)
  {
   polygon(x = c(24*l, 24*l+dusk[i], 24*l+dusk[i], 24*l),
           y=c(min.expr-1.5*width.rectangule, min.expr-1.5*width.rectangule, 
               min.expr-0.5*width.rectangule, min.expr-0.5*width.rectangule),
           lwd=2,border=colors[i])
   polygon(x = c(24*l+dusk[i],24*(l+1),24*(l+1),24*l+dusk[i]),
           y=c(min.expr-1.5*width.rectangule, min.expr-1.5*width.rectangule, 
               min.expr-0.5*width.rectangule, min.expr-0.5*width.rectangule),
           lwd=2,border=colors[i],col=colors[i])
  }
  
  lines(time.points,transitions.sd.ld[i,],type="l",col=colors[i],lwd=1.5,lty=2)
  lines(time.points,fixed.sd.ld[i,],type="l",col=colors[i],lwd=1.5,lty=3)
  combined.wave.i <- apply(X = matrix(data = c(transitions.sd.ld[i,],fixed.sd.ld[i,]),ncol=2),MARGIN = 1,FUN = max)
  lines(time.points,combined.wave.i,type="l",col=colors[i],lwd=2,lty=1)

  k <- k + 1
  Sys.sleep(0.25)
  print(i)
 }
 k <- 1
}







