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
plot.ld.dd(gene.id = "ostta04g00450",gene.name = "MCM5", gene.expression = gene.expression)


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

plot.ld <- function(gene.id, gene.name, gene.expression)
{

  current.gene.expression.ld.ll <- gene.expression[gene.id,1:18]
  
  min.expression <- min(current.gene.expression.ld.ll)
  max.expression <- max(current.gene.expression.ld.ll)
  range.expression <- max.expression - min.expression
  
  expression.step <- floor(range.expression / 5)
  
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
  return(0)  
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

plot.ld.sd(gene.id = "ostta01g05650",gene.name = "ATG8", gene.expression)



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
  axis(side = 1,pos = min.expr - 1.1*range.expr, at = seq(from=1,to=30),
       labels = rep(paste("ZT",seq(from=0,to=20,by=4)),5),las=2,lwd=3)
  
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

## Dyname change from LD to SD evaluated with ND
wave <- function(k, alpha, phi, x)
{
 res <- (k  + alpha * cos((2*pi/24) * (x - phi))) 
 res[res < 0] <- 0
 return(res)
}

gene.id <- "ostta08g00710"

rhythmic.ld.nd.sd.single.peak <- intersect(intersect(setdiff(rhythmic.genes.ld,rhythmic.genes.ld.12),
                    setdiff(rhythmic.genes.nd,rhythmic.genes.nd.12)),
          setdiff(rhythmic.genes.sd,rhythmic.genes.sd.12))

library(circacompare)
length(rhythmic.ld.nd.sd.single.peak)
time.points <- seq(from=0,by=4,length.out = 18)
nd.time.points <- seq(from=0,to=69,by=3)

ld.zt <- paste("ld",paste0("zt",sprintf(fmt = "%02d",seq(from=0,to=20,by=4))),sep="_")
nd.zt <- paste0("ZT",seq(from=0,to=21,by=3))
sd.zt <- paste("sd",paste0("zt",sprintf(fmt = "%02d",seq(from=0,to=20,by=4))),sep="_")

model.deviations <- vector(mode = "numeric",length = length(rhythmic.ld.nd.sd.single.peak))
for(i in 1:length(rhythmic.ld.nd.sd.single.peak))
{
 print(i)
 gene.id <- rhythmic.ld.nd.sd.single.peak[i]
 
 current.gene.expression.ld <- gene.expression[gene.id,c(paste(ld.zt,1,sep="_"),paste(ld.zt,2,sep="_"),paste(ld.zt,3,sep="_"))]
 current.gene.expression.nd <- nd.gene.expression[gene.id,c(paste(nd.zt,1,sep="_"),paste(nd.zt,2,sep="_"),paste(nd.zt,3,sep="_"))]
 current.gene.expression.sd <- gene.expression[gene.id,c(paste(sd.zt,1,sep="_"),paste(sd.zt,2,sep="_"),paste(sd.zt,3,sep="_"))]
 
 ld.sd.df <- data.frame(time=c(time.points,time.points),
                        measure=c(unlist(current.gene.expression.ld), 
                                  unlist(current.gene.expression.sd)),
                        group=c(rep("ld",18),rep("sd",18)))
 out.ld.sd <- circacompare(x = ld.sd.df, col_time = "time", 
                           col_group = "group", col_outcome = "measure",
                           alpha_threshold = 1)
 circa.values <- out.ld.sd$summary$value
 names(circa.values) <- out.ld.sd$summary$parameter
 
 phase.ld <- circa.values["ld peak time hours"]
 phase.sd <- circa.values["sd peak time hours"]
 
 
 ld.nd.df <- data.frame(time=c(time.points,nd.time.points),
                        measure=c(unlist(current.gene.expression.ld),
                                  unlist(current.gene.expression.nd)),
                        group=c(rep("ld",18),rep("nd",24)))
 out.ld.nd <- circacompare(x = ld.nd.df, col_time = "time", 
                           col_group = "group", col_outcome = "measure",
                           alpha_threshold = 1)
 
 circa.values <- out.ld.nd$summary$value
 names(circa.values) <- out.ld.nd$summary$parameter
 
 phase.nd <- circa.values["nd peak time hours"]
 
 distancia.ld.sd <- abs(phase.ld - phase.sd)

 if(distancia.ld.sd > 12)
 {
  distancia.ld.sd <- 24 - distancia.ld.sd
 }
 
 if(phase.ld > phase.sd)
 {
  predicted.phase.nd <- phase.sd + distancia.ld.sd/2
 } else
 {
  predicted.phase.nd <- phase.ld + distancia.ld.sd/2
 }
 
 if (predicted.phase.nd > 24)
 {
  predicted.phase.nd <- predicted.phase.nd - 24
 }
 
 model.deviations[i] <- abs(phase.nd - predicted.phase.nd)
}

hist(model.deviations,breaks=0:24)
hist(model.deviations,breaks=seq(from=0,to=24,by=1))

sum(model.deviations < 1)/length(model.deviations)
sum(model.deviations < 2)/length(model.deviations)
sum(model.deviations < 3)/length(model.deviations)
sum(model.deviations < 4)/length(model.deviations)
sum(model.deviations < 5)/length(model.deviations)



nd.gene.expression
colnames(nd.gene.expression)
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
    
    
    #lines(nd.time.points,current.gene.expression.nd,type="o",col="red",lwd=0.25,pch=19,cex=0.5)   
    
     
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

two.waves.combination.3 <- function(x,beta1,beta2,phi,offset) #con la nueva wave
{
 phi1 <- phi
 phi2 <- (phi + offset) %% 24
 wave1 <- wave.form(mesor = 1, amplitude = beta1, period = 24, phase = phi1, time = x)
 wave2 <- wave.form(mesor = 1, amplitude = beta2, period = 24, phase = phi2, time = x)
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
current.gene <- "ostta18g01250"
plot.ld.sd(gene.id=current.gene, gene.name="GBSS", gene.expression)
plot.expression.profile(gene.id ="ostta06g02940", gene.name = "GBSS", 
                        gene.expression,entrainment="SD",
                        free.running="LL",param=F)
plot.expression.profile(gene.id ="ostta06g02940", gene.name = "GBSS", 
                        gene.expression,entrainment="SD",
                        free.running="DD",param=F)

plot.ld.ll(gene.id=current.gene, gene.name="OsttaGBSS", gene.expression)
plot.ld.dd(gene.id=current.gene, gene.name="OsttaGBSS", gene.expression)
plot.sd.ll(gene.id=current.gene, gene.name="OsttaGBSS", gene.expression)
plot.sd.dd(gene.id=current.gene, gene.name="OsttaGBSS", gene.expression)







write.table(x = subset(ostta2atha, ostta %in% genes.two.peaks.sd.one.peak.ll.dd), file = "/home/fran/tmp/good_candidates_2_peaks.tsv",sep = "\t",row.names = F,col.names = F,quote=F)

head(ostta2atha)
      
      
      

gene.id <- "ostta07g01610" #FtsZ
gene.id <- "ostta14g01160" #HMG
gene.id <- "ostta04g04000" #PWD
gene.id <- "ostta18g01250" #FNR
gene.id <- "ostta04g02740" #PRK
gene.id <- "ostta04g00070" #Psb33
gene.id <- "ostta15g02670" #PSAD
gene.id <- "ostta06g02940" #GBSSok
gene.id <- "ostta07g02570" #SMC6A

gene.name <- "GBSS"


bimodal.gene.model <- function(gene.id, gene.name)
{
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
 time.points.smooth <- seq(from=0,to=68,by=1)
 gene.expression.df.2 <- data.frame(x=time.points,
                                    y=current.gene.expression.sd/mean(current.gene.expression.sd))
 
 m.2 <- nls( y ~ two.waves.combination(x,beta1,beta2,phi,offset), data = gene.expression.df.2,
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
 
 dos.ondas.2 <- two.waves.combination(x=seq(from=0,to=68,by=0.01), 
                                      beta1=beta1.estimate.2,
                                      beta2=beta2.estimate.2,
                                      phi=phi1.estimate.2,
                                      offset=phi2.estimate.2 - phi1.estimate.2)
 
 wave1 <- wave.form(mesor = 1, amplitude = beta1.estimate.2, phase = phi1.estimate.2,period = 24, time = seq(from=0,to=68,by=0.01))
 wave2 <- wave.form(mesor = 1, amplitude = beta2.estimate.2, phase = phi2.estimate.2,period = 24, time = seq(from=0,to=68,by=0.01))
 
 time.points <- seq(from=0,by=4,length.out = 18)

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
 amplitude.sd.1 <- beta2.estimate.2*mesor.sd
 phase.sd.2 <- phi2.estimate.2
 # phase.sd.1 <- phi2.estimate.2
 # phase.sd.2 <- phi1.estimate.2
 
 phase.sd <- phase.sd.2
 
 # transitions.sd.ld <- matrix(nrow = 48,ncol=18)
 # fixed.sd.ld <- matrix(nrow = 48,ncol=18)
 transitions.sd.ld <- matrix(nrow = 48,ncol=length(time.points.smooth))
 fixed.sd.ld <- matrix(nrow = 48,ncol=length(time.points.smooth))
 
 mesors <- seq(from=mesor.sd,to=mesor.ld,length.out = 48)
 amplitudes <- seq(from=amplitude.sd.1,to=amplitude.ld,length.out = 48)
 
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
  # transitions.sd.ld[i,] <- wave.form(mesor = mesors[i], amplitude = amplitudes[i], 
  #                                    phase = phases[i],time = time.points,period = 24)   
  # fixed.sd.ld[i,] <- wave.form(mesor = mesor.sd, amplitude = amplitude.sd, 
  #                              phase = phase.sd.1,time = time.points,period = 24)
  
  transitions.sd.ld[i,] <- wave.form(mesor = mesors[i], amplitude = amplitudes[i], 
                                     phase = phases[i],time = time.points.smooth,period = 24)   
  fixed.sd.ld[i,] <- wave.form(mesor = mesor.sd, amplitude = amplitude.sd, 
                               phase = phase.sd.1,time = time.points.smooth,period = 24)
  
 }
 
 months <- c("January", "February", "March", "April","May", "June", "July", 
             "August", "September", "October", "November", "December")
 
 months <- rep(months,each=8)
 colors <- colorRampPalette(c("red","blue"))(48)
 dusk <- seq(from=8,to=16,length.out = 48)
 k <- 1
 
 image.number <- 1
 image.dir <- paste0("gif_images/",gene.id)
 dir.create(path = image.dir)
 
 for(i in 1:48)
 {
  png(filename = paste0(image.dir,"/",paste0(paste(gene.id,sprintf(fmt = "%03d",image.number),sep="_"),".png")))
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
  
  text(x = 32,y=min.expr-2.2*width.rectangule,labels=paste(gene.id,gene.name,sep=" - "),cex=1.2,adj=0.5)
  
  
  # lines(time.points,transitions.sd.ld[i,],type="l",col=colors[i],lwd=1.5,lty=2)
  # lines(time.points,fixed.sd.ld[i,],type="l",col=colors[i],lwd=1.5,lty=3)
  # combined.wave.i <- apply(X = matrix(data = c(transitions.sd.ld[i,],fixed.sd.ld[i,]),ncol=2),MARGIN = 1,FUN = max)
  # lines(time.points,combined.wave.i,type="l",col=colors[i],lwd=2,lty=1)

  lines(time.points.smooth,transitions.sd.ld[i,],type="l",col=colors[i],lwd=1.5,lty=2)
  lines(time.points.smooth,fixed.sd.ld[i,],type="l",col=colors[i],lwd=1.5,lty=3)
  combined.wave.i <- apply(X = matrix(data = c(transitions.sd.ld[i,],fixed.sd.ld[i,]),ncol=2),MARGIN = 1,FUN = max)
  lines(time.points.smooth,combined.wave.i,type="l",col=colors[i],lwd=2,lty=1)
  
  dev.off()
  image.number <- image.number + 1
  k <- k + 1
  # Sys.sleep(0.25)
  # print(i)
 }
 
 for(i in 48:1)
 {
  png(filename = paste0(image.dir,"/",paste0(paste(gene.id,sprintf(fmt = "%03d",image.number),sep="_"),".png")))
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
  
  # lines(time.points,transitions.sd.ld[i,],type="l",col=colors[i],lwd=1.5,lty=2)
  # lines(time.points,fixed.sd.ld[i,],type="l",col=colors[i],lwd=1.5,lty=3)
  # combined.wave.i <- apply(X = matrix(data = c(transitions.sd.ld[i,],fixed.sd.ld[i,]),ncol=2),MARGIN = 1,FUN = max)
  # lines(time.points,combined.wave.i,type="l",col=colors[i],lwd=2,lty=1)
  
  lines(time.points.smooth,transitions.sd.ld[i,],type="l",col=colors[i],lwd=1.5,lty=2)
  lines(time.points.smooth,fixed.sd.ld[i,],type="l",col=colors[i],lwd=1.5,lty=3)
  combined.wave.i <- apply(X = matrix(data = c(transitions.sd.ld[i,],fixed.sd.ld[i,]),ncol=2),MARGIN = 1,FUN = max)
  lines(time.points.smooth,combined.wave.i,type="l",col=colors[i],lwd=2,lty=1)
  
  dev.off()
  image.number <- image.number + 1
  k <- k + 1
  # Sys.sleep(0.25)
  # print(i)
 }
 
 imgs <- list.files(image.dir, full.names = TRUE)
 img.list <- lapply(imgs, image_read)
 img.joined <- image_join(img.list)
 img.animated <- image_animate(img.joined, fps = 5)
 image_write(image = img.animated,
             path = paste0(c("gif_images/",paste(gene.id,gene.name,sep="_"),".gif"),collapse=""))
 unlink(x = image.dir,recursive = TRUE)
}

bimodal.gene.model(gene.id = "ostta06g02940",gene.name = "GBSS")
bimodal.gene.model(gene.id = "ostta04g04000",gene.name = "PWD")



bimodal.gene.model(gene.id = "ostta07g02570") #SMC6A
bimodal.gene.model(gene.id = "ostta15g02670") #PSAD

i <- 1

bimodal.gene.model(gene.id = genes.two.peaks.sd.one.peak.ll.dd[i])
i <- i + 1
i <- 3
i <- 10
i <- 11
i <- 27
i <- 34  # good
i <- 35 # good
i <- 49 # PWD
i <- 50
i <- 56
i <- 68

gene.id %in% rhythmic.genes.nd.12




i <- 23

gene.i <- bona.fide.circadian[i]

plot.ld.ll(gene.id = gene.i,gene.name = gene.i, gene.expression)
plot.ld.dd(gene.id = gene.i,gene.name = gene.i, gene.expression)
plot.sd.ll.param(gene.id=gene.i, gene.name=gene.i, gene.expression,param=F)
plot.sd.ll.param(gene.id=gene.i, gene.name=gene.i, gene.expression,param=T)
plot.sd.dd(gene.id = gene.i,gene.name = gene.i, gene.expression)
i <- i + 1
print(i)

#19 21 23 25 26

wave.form <- function(mesor, amplitude,period,phase,time=seq(from=0,to=72,by=0.01))
{
 y <- mesor + amplitude*cos(((2*pi)/period)*(time - phase))
 return(y)
}

plot.expression.profile <- function(gene.id, gene.name, gene.expression,
                             entrainment,free.running,param=F)
{
 if (entrainment == "SD")
 {
  zts <- paste("sd",paste0("zt",sprintf(fmt = "%02d",seq(from=0,to=20,by=4))),sep="_")
  entrainment.col <- "red"
  free.running.col <- "lightsalmon"
 } else if (entrainment == "LD")
 {
  zts <- paste("ld",paste0("zt",sprintf(fmt = "%02d",seq(from=0,to=20,by=4))),sep="_")
  entrainment.col <- "blue"
  free.running.col <- "lightblue"
 } else
 {
  print("Unkwnon entrainment. Possible entrainments SD and LD")
 }
 
 if (free.running == "LL")
 {
  free.running.indeces <- c(4,5)
 } else if (free.running == "DD")
 {
  free.running.indeces <- c(6,7)
 } else
 {
  print("Unknown free running. Possible options LL and DD")
 }
 
 current.gene.expression.entrainment <- gene.expression[gene.id,c(paste(zts,1,sep="_"),
                                                                  paste(zts,2,sep="_"),
                                                                  paste(zts,3,sep="_"))]
 current.gene.expression.free.running <- gene.expression[gene.id,c(paste(zts,free.running.indeces[1],sep="_"),
                                                                   paste(zts,free.running.indeces[2],sep="_"))]
 current.gene.expression <- c(current.gene.expression.entrainment,current.gene.expression.free.running)
 min.expression <- min(current.gene.expression)
 max.expression <- max(current.gene.expression)
 range.expression <- max.expression - min.expression
 
 expression.step <- floor(range.expression / 5)
 time.points.entrainment <- seq(from=0,by=4,length.out = 18)
 time.points.free.running <- seq(from=0,by=4,length.out = 12)
 
 
 if(param)
 {
  circacompare.data <- data.frame(time=c(time.points.entrainment,time.points.free.running),
                               measure=c(current.gene.expression.entrainment,
                                         current.gene.expression.free.running),
                               group=c(rep("Entrainment",18),rep("Free_running",12)))
  
  result.i <- circacompare(x = circacompare.data, 
                          col_time = "time", 
                          col_group = "group", 
                          col_outcome = "measure",
                          alpha_threshold = 1)
  result.i.values <- result.i$summary$value
  names(result.i.values) <- result.i$summary$parameter
  time.entrainment <- seq(from=0,to=72,by=0.01)
  time.free.running <- seq(from=72,to=120,by=0.01)
  line.entrainment <- wave.form(mesor = result.i.values["Entrainment mesor estimate"],
            amplitude = result.i.values["Entrainment amplitude estimate"],
            period = result.i.values["Shared period estimate"],phase = result.i.values["Entrainment peak time hours"],
            time = time.entrainment)
  line.free.running <- wave.form(mesor = result.i.values["Free_running mesor estimate"],
                       amplitude = result.i.values["Free_running amplitude estimate"],
                       period = result.i.values["Shared period estimate"],
                       phase = result.i.values["Free_running peak time hours"],
                       time = time.free.running)
  plot(x = c(time.points.entrainment,time.points.free.running+max(time.points.sd)+4),#type="o",
       y = current.gene.expression,lwd=1,col=entrainment.col,axes=F,xlab="",ylab="FPKM",
       ylim=c(min.expression-expression.step,max.expression),
       cex.lab=1.3,main=paste(gene.id, gene.name,sep=" - "),cex.main=2,pch=20,cex=0.8)
  lines(x = c(time.entrainment,time.free.running), 
        y = c(line.entrainment,line.free.running), 
        type="l",
        col=entrainment.col,lwd=3)
 } else
 {
  plot(x= c(time.points.entrainment, 72 +time.points.free.running),
       y = current.gene.expression,type="o",lwd=3,col=entrainment.col,axes=F,xlab="",ylab="FPKM",
       ylim=c(min.expression-expression.step,max.expression),
       cex.lab=1.3,main=paste(gene.id, gene.name,sep=" - "),cex.main=2)
 }
 
 axis(side=2,lwd=3)
 axis(side = 1,pos = min.expression - 1.1*expression.step, at = seq(from=0,to=116,by=4),
      labels = rep(paste("ZT",seq(from=0,to=20,by=4)),5),las=2,lwd=3)
 if(entrainment == "SD")
 {
  polygon(x = c(0,8,8,0),y=c(min.expression-expression.step/2,
                             min.expression-expression.step/2,
                             min.expression-expression.step,
                             min.expression-expression.step),lwd=2,border="red")
  polygon(x = c(8,24,24,8),y=c(min.expression-expression.step/2,
                               min.expression-expression.step/2,
                               min.expression-expression.step,
                               min.expression-expression.step),lwd=2,border="red",col="red")
  
  polygon(x = c(24,32,32,24),y=c(min.expression-expression.step/2,
                                 min.expression-expression.step/2,
                                 min.expression-expression.step,
                                 min.expression-expression.step),lwd=2,border="red")
  polygon(x = c(32,48,48,32),y=c(min.expression-expression.step/2,
                                 min.expression-expression.step/2,
                                 min.expression-expression.step,
                                 min.expression-expression.step),lwd=2,border="red",col="red")
  
  polygon(x = c(48,56,56,48),y=c(min.expression-expression.step/2,
                                 min.expression-expression.step/2,
                                 min.expression-expression.step,
                                 min.expression-expression.step),lwd=2,border="red")
  polygon(x = c(56,72,72,56),y=c(min.expression-expression.step/2,
                                 min.expression-expression.step/2,
                                 min.expression-expression.step,
                                 min.expression-expression.step),lwd=2,border="red",col="red")
  if(free.running == "LL")
  {
   polygon(x = c(72,80,80,72),y=c(min.expression-expression.step/2,
                                  min.expression-expression.step/2,
                                  min.expression-expression.step,
                                  min.expression-expression.step),lwd=2,border="red")
   polygon(x = c(80,96,96,80),y=c(min.expression-expression.step/2,
                                  min.expression-expression.step/2,
                                  min.expression-expression.step,
                                  min.expression-expression.step),lwd=2,border="red",col="salmon")
   
   polygon(x = c(96,104,104,96),y=c(min.expression-expression.step/2,
                                    min.expression-expression.step/2,
                                    min.expression-expression.step,
                                    min.expression-expression.step),lwd=2,border="red")
   polygon(x = c(104,120,120,104),y=c(min.expression-expression.step/2,
                                      min.expression-expression.step/2,
                                      min.expression-expression.step,
                                      min.expression-expression.step),lwd=2,border="red",col="salmon")
  } else if (free.running == "DD")
  {
   polygon(x = c(72,80,80,72),y=c(min.expression-expression.step/2,
                                  min.expression-expression.step/2,
                                  min.expression-expression.step,
                                  min.expression-expression.step),lwd=2,border="red",col="salmon")
   polygon(x = c(80,96,96,80),y=c(min.expression-expression.step/2,
                                  min.expression-expression.step/2,
                                  min.expression-expression.step,
                                  min.expression-expression.step),lwd=2,border="red",col="red")
   
   polygon(x = c(96,104,104,96),y=c(min.expression-expression.step/2,
                                    min.expression-expression.step/2,
                                    min.expression-expression.step,
                                    min.expression-expression.step),lwd=2,border="red",col="salmon")
   polygon(x = c(104,120,120,104),y=c(min.expression-expression.step/2,
                                      min.expression-expression.step/2,
                                      min.expression-expression.step,
                                      min.expression-expression.step),lwd=2,border="red",col="red")
  }
 } else if (entrainment == "LD")
 {
  polygon(x = c(0,16,16,0),y=c(min.expression-expression.step/2,
                             min.expression-expression.step/2,
                             min.expression-expression.step,
                             min.expression-expression.step),lwd=2,border="blue")
  polygon(x = c(16,24,24,16),y=c(min.expression-expression.step/2,
                               min.expression-expression.step/2,
                               min.expression-expression.step,
                               min.expression-expression.step),lwd=2,border="blue",col="blue")
  
  polygon(x = c(24,40,40,24),y=c(min.expression-expression.step/2,
                                 min.expression-expression.step/2,
                                 min.expression-expression.step,
                                 min.expression-expression.step),lwd=2,border="blue")
  polygon(x = c(40,48,48,40),y=c(min.expression-expression.step/2,
                                 min.expression-expression.step/2,
                                 min.expression-expression.step,
                                 min.expression-expression.step),lwd=2,border="blue",col="blue")
  
  polygon(x = c(48,64,64,48),y=c(min.expression-expression.step/2,
                                 min.expression-expression.step/2,
                                 min.expression-expression.step,
                                 min.expression-expression.step),lwd=2,border="blue")
  polygon(x = c(64,72,72,64),y=c(min.expression-expression.step/2,
                                 min.expression-expression.step/2,
                                 min.expression-expression.step,
                                 min.expression-expression.step),lwd=2,border="blue",col="blue")
  
  if(free.running == "LL")
  {
   polygon(x = c(72,88,88,72),y=c(min.expression-expression.step/2,
                                  min.expression-expression.step/2,
                                  min.expression-expression.step,
                                  min.expression-expression.step),lwd=2,border="blue")
   polygon(x = c(88,96,96,88),y=c(min.expression-expression.step/2,
                                  min.expression-expression.step/2,
                                  min.expression-expression.step,
                                  min.expression-expression.step),lwd=2,border="blue",col="lightblue")
   
   polygon(x = c(96,112,112,96),y=c(min.expression-expression.step/2,
                                    min.expression-expression.step/2,
                                    min.expression-expression.step,
                                    min.expression-expression.step),lwd=2,border="blue")
   polygon(x = c(112,120,120,112),y=c(min.expression-expression.step/2,
                                      min.expression-expression.step/2,
                                      min.expression-expression.step,
                                      min.expression-expression.step),lwd=2,border="blue",col="lightblue")
  } else if (free.running == "DD")
  {
   polygon(x = c(72,88,88,72),y=c(min.expression-expression.step/2,
                                  min.expression-expression.step/2,
                                  min.expression-expression.step,
                                  min.expression-expression.step),lwd=2,border="blue",col="lightblue")
   polygon(x = c(88,96,96,88),y=c(min.expression-expression.step/2,
                                  min.expression-expression.step/2,
                                  min.expression-expression.step,
                                  min.expression-expression.step),lwd=2,border="blue",col="blue")
   
   polygon(x = c(96,112,112,96),y=c(min.expression-expression.step/2,
                                    min.expression-expression.step/2,
                                    min.expression-expression.step,
                                    min.expression-expression.step),lwd=2,border="blue", col="lightblue")
   polygon(x = c(112,120,120,112),y=c(min.expression-expression.step/2,
                                      min.expression-expression.step/2,
                                      min.expression-expression.step,
                                      min.expression-expression.step),lwd=2,border="blue",col="blue")
  }
 }
}


gene.i <- "ostta03g00350" # fructose 1,6-bisphosphate phosphatase
gene.i <- "ostta10g01980" # BETA CAROTENOID HYDROXYLASE 1 (BCH1) bueno para phase swift
gene.i <- "ostta04g01820" # ELIP2
gene.i <- "ostta17g00100" # GUN4
gene.i <- "ostta16g02310" # LHCB4.3
gene.i <- "ostta03g01980" # LIGHT-HARVESTING-LIKE 3:2 (LIL3:2)

gene.i <- "ostta04g02740" # PHOSPHORIBULOKINASE (PRK)
gene.i <- "ostta03g05500" # SBPASE

gene.i <- "ostta15g02670" # PSAD
gene.i <- "ostta04g02510" # RUBISCO ACTIVASE (RCA)



 
gene.i <- "ostta07g01610" # FTSZ2-1 chloroplast division #bona.fide.circadian[i]
gene.i <- "ostta08g01100" # CHLG chlorophyll synthesis


gene.i <- "ostta02g03330" # CAB
gene.i <- "ostta03g05620" # PHOTOLYASE 1 (PHR1)

gene.i <- "ostta02g03860" # PSAE
gene.i <- "ostta04g01790" # PSAF

plot.expression.profile(gene.id=gene.i, gene.name=gene.i, gene.expression,entrainment="LD",free.running="LL",param=T)
plot.expression.profile(gene.id=gene.i, gene.name=gene.i, gene.expression,entrainment="LD",free.running="LL",param=F)

plot.expression.profile(gene.id=gene.i, gene.name=gene.i, gene.expression,entrainment="LD",free.running="DD",param=T)
plot.expression.profile(gene.id=gene.i, gene.name=gene.i, gene.expression,entrainment="LD",free.running="DD",param=F)


plot.expression.profile(gene.id=gene.i, gene.name=gene.i, gene.expression,entrainment="SD",free.running="LL",param=T)
plot.expression.profile(gene.id=gene.i, gene.name=gene.i, gene.expression,entrainment="SD",free.running="LL",param=F)

plot.expression.profile(gene.id=gene.i, gene.name=gene.i, gene.expression,entrainment="SD",free.running="DD",param=T)
plot.expression.profile(gene.id=gene.i, gene.name=gene.i, gene.expression,entrainment="SD",free.running="DD",param=F)

i <- i + 1
print(i)

#19 21 23 25 26

gene.ld.dd <- read.table(file="gene_sets/LD_and_DD_rhythmic_genes.tsv",header = F,as.is = T)[[1]]
gene.ld.dd.atha <- subset(ostta2atha, ostta %in% gene.ld.dd)[[2]] 
write.table(x = subset(ostta2atha, ostta %in% gene.ld.dd), 
            file="/home/fran/tmp/LD_SD_atha_genes.tsv",row.names = F,quote=F)


library(org.At.tair.db)
enrich.go.ld.dd.atha <- enrichGO(gene = gene.ld.dd.atha,
                                              OrgDb = org.At.tair.db,
                                              ont = "BP",
                                              pAdjustMethod = "BH",
                                              pvalueCutoff  = 0.05,
                                              keyType = "TAIR")
write.table(x = as.data.frame(enrich.go.ld.dd.atha),file="/home/fran/tmp/enrichment_atha_ld_dd.tsv",quote = F,sep = "\t",row.names = F,col.names = F)


## SPLS
sd.ld.genes <- intersect(complete.ld.rhythmic.genes,complete.sd.rhythmic.genes)
length(sd.ld.genes)

"ostta06g02340" %in%sd.ld.genes

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

## Predictors are gene expression
rhythmic.expression.ld.sd <- gene.expression.ld.sd[rhythmic.tfs,]
X <- t(rhythmic.expression.ld.sd) 
dim(X)
rownames(X)
colnames(X)
## Response variables are percentage of cells in the different cell cycle phases
cell.cycle.data <- read.table(file="physiological_data/cell_cycle_data_spls.csv",header = T,sep = ",",as.is=T)
dim(cell.cycle.data)

Y <- cell.cycle.data[,3:5]
dim(Y)
rownames(X)
colnames(Y)

myresult <- pls(X, Y, ncomp = 10, mode = "regression")
typeof(myresult)
attributes(myresult)

## We plot the explained variance by each component
tiff(filename = "explained_variance_10_component_replicates.tiff")
par(mfrow = c(1,2))
barplot(myresult$prop_expl_var$X, las = 2, main = "X",col=rainbow(5),ylim=c(0,0.7))
lines(x=c(0,12),y=c(0.1,0.1),col="black",lwd=4)
barplot(myresult$prop_expl_var$Y, las = 2, main = "Y",col=rainbow(5),ylim=c(0,0.7))
lines(x=c(0,12),y=c(0.25,0.25),col="black",lwd=4)
dev.off()

## Three components are selected and we apply PLS regression with ncomp=3
myresult <- pls(X, Y, ncomp = 3, mode = "regression", scale = TRUE) 

## Visualization of the individuals using the first two new components XY combined
# tiff(filename = "individuals_visualization_on_the_new_components_joint_XY_with_replicates.tiff")
# plotIndiv(object = myresult, comp = 1:2, rep.space = "XY-variate",
#           ind.names = rownames(X),
#           group = c(1,2,3,4,5,6,
#                     1,2,3,4,5,6,
#                     1,2,3,4,5,6,
#                     7,8,9,10,11,12,
#                     7,8,9,10,11,12,
#                     7,8,9,10,11,12),
#           # col = rep(rainbow(4),each=3), style = "ggplot2", ellipse = FALSE, ellipse.level = 0.95, centroid = FALSE,lwd=3,cex=6,xlim=c(-30,30)) 
#           # col = rainbow(4), 
#           col = rainbow(12),#c("black", "red","blue", "green"), 
#           style = "ggplot2", ellipse = FALSE, ellipse.level = 0.95, centroid = FALSE,lwd=3,cex=6,xlim=c(-60,60)) 
# dev.off()

## Individuals clustered according to their genotypes

tiff(filename = "individuals_visualization_on_the_new_components_separated_X_Y_with_replicates.tiff")
plotIndiv(object = myresult, comp = 1:2,point.lwd = 0.5,cex = 4,
          ind.names = rownames(X),
          group = c(1,2,3,4,5,6,
                    1,2,3,4,5,6,
                    1,2,3,4,5,6,
                    7,8,9,10,11,12,
                    7,8,9,10,11,12,
                    7,8,9,10,11,12),
          # col = rep(rainbow(4),each=3), style = "ggplot2", ellipse = FALSE, 
          col = rainbow(12), style = "ggplot2", ellipse = F, 
          ellipse.level = 0.95, centroid = FALSE,
          star = FALSE, legend = FALSE, abline = TRUE, alpha = 0)
dev.off()

## Individuals clustered according to their genotypes
tiff(filename = "correlation_circle_pls_with_replicates.tiff", 
     width = 5, height = 5, units = 'in', res = 200)
plotVar(object = myresult, comp = 1:2, cex = c(3, 4), col = c("forestgreen", "red3"))
dev.off()



## Cross validation of the model using LOO (Leave One Out)
myperfLoo = perf(myresult, validation = "loo", progressBar = TRUE)

myperfLoo$measures$MSEP
myperfLoo$PRESS
myperfLoo$R2
myperfLoo$Q2
myperfLoo$Q2.total
myperfLoo$RSS

sum(myperfLoo$measures$MSEP$summary[,3])/length(sd.ld.genes)
## FC = 1.75 ---> 0.016

tiff(filename = "MSEP_pls_with_replicates.tiff",
     width = 6, height = 6, units = "in", res = 150)
plot(x = myperfLoo,
     criterion = "MSEP",
     xlab = "number of components",
     ylab = NULL,
     LimQ2 = 0.0975,
     LimQ2.col = "none",type="line", lwd=3)
dev.off()

## Applying Sparse PLS
mySPresult <- spls(X, Y, ncomp = 3, mode = 'regression', 
                   keepX = c(40, 40, 40), 
                   keepY = c(3, 3, 3),
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
plotVar(mySPresult, comp = 1:2, cex = c(3, 5), col = c("forestgreen", "red3"))
dev.off()

## Bipartite representation
color.edge <- colorRampPalette(c("red4", "white", "darkgreen"))
spl.th <- 0.5
tiff(filename = "bipartite_graph_with_replicates.tiff")#, height = 6, width = 6,
#units = "in", res = 150)
res <- network(mySPresult, comp = 1:2, cutoff = spl.th, shape.node = c("rectangle", "circle"),cex.node.name = 0.9,
               color.node = c("white", "coral1"), color.edge = color.edge(10),save = "jpeg",name.save = "network_cca1")
dev.off()


## Proteomics ##

zt0 <- read.table(file="/home/fran/Nextcloud2/Microalgas/analysis_minotaur/swath/ld_1_2_3/fran_analysis_3/zt00.tsv",header = T,as.is=T,sep="\t")
zt0 <- zt0[,c(1,6:11)]
head(zt0)
nrow(zt0)
zt0.3 <- read.table(file="/home/fran/Nextcloud2/Microalgas/analysis_minotaur/swath/ld_1_2_3/fran_analysis_3/zt00_3.tsv",header = T,as.is=T,sep="\t")
head(zt0.3)
nrow(zt0.3)

zt4 <- read.table(file="/home/fran/Nextcloud2/Microalgas/analysis_minotaur/swath/ld_1_2_3/fran_analysis_3/zt04.tsv",header = T,as.is=T,sep="\t")
zt4 <- zt4[,c(1,6:11)]
head(zt4)
zt4.3 <- read.table(file="/home/fran/Nextcloud2/Microalgas/analysis_minotaur/swath/ld_1_2_3/fran_analysis_3/zt04_3.tsv",header = T,as.is=T,sep="\t")
head(zt4.3)

zt8 <- read.table(file="/home/fran/Nextcloud2/Microalgas/analysis_minotaur/swath/ld_1_2_3/fran_analysis_3/zt08.tsv",header = T,as.is=T,sep="\t")
zt8 <- zt8[,c(1,6:11)]
head(zt8)
zt8.3 <- read.table(file="/home/fran/Nextcloud2/Microalgas/analysis_minotaur/swath/ld_1_2_3/fran_analysis_3/zt08_3.tsv",header = T,as.is=T,sep="\t")
head(zt8.3)

zt12 <- read.table(file="/home/fran/Nextcloud2/Microalgas/analysis_minotaur/swath/ld_1_2_3/fran_analysis_3/zt12.tsv",header = T,as.is=T,sep="\t")
zt12 <- zt12[,c(1,6:11)]
head(zt12)
zt12.3 <- read.table(file="/home/fran/Nextcloud2/Microalgas/analysis_minotaur/swath/ld_1_2_3/fran_analysis_3/zt12_3.tsv",header = T,as.is=T,sep="\t")
head(zt12.3)

zt16 <- read.table(file="/home/fran/Nextcloud2/Microalgas/analysis_minotaur/swath/ld_1_2_3/fran_analysis_3/zt16.tsv",header = T,as.is=T,sep="\t")
zt16 <- zt16[,c(1,6:11)]
head(zt16)
zt16.3 <- read.table(file="/home/fran/Nextcloud2/Microalgas/analysis_minotaur/swath/ld_1_2_3/fran_analysis_3/zt16_3.tsv",header = T,as.is=T,sep="\t")
head(zt16.3)

zt20 <- read.table(file="/home/fran/Nextcloud2/Microalgas/analysis_minotaur/swath/ld_1_2_3/fran_analysis_3/zt20.tsv",header = T,as.is=T,sep="\t")
zt20 <- zt20[,c(1,6:11)]
head(zt20)
zt20.3 <- read.table(file="/home/fran/Nextcloud2/Microalgas/analysis_minotaur/swath/ld_1_2_3/fran_analysis_3/zt20_3.tsv",header = T,as.is=T,sep="\t")
head(zt20.3)


total.proteins <- unique(c(zt0$Peak.Name,zt0.3$Peak.Name,
                           zt4$Peak.Name,zt4.3$Peak.Name, 
                           zt8$Peak.Name,zt8.3$Peak.Name,
                           zt12$Peak.Name,zt12.3$Peak.Name,
                           zt16$Peak.Name,zt16.3$Peak.Name,
                           zt20$Peak.Name,zt20.3$Peak.Name))
length(total.proteins)

prot.zt0.1 <- zt0$SWATH_O.tauri_25.11.20_1.1
names(prot.zt0.1) <- zt0$Peak.Name
prot.zt0.2 <- zt0$SWATH_O.tauri_25.11.20_1.2
names(prot.zt0.2) <- zt0$Peak.Name
prot.zt0.3 <- zt0$SWATH_O.tauri_25.11.20_1.3
names(prot.zt0.3) <- zt0$Peak.Name
prot.zt0.4 <- zt0$Swath.O..tauri_01122020_1.1.prima
names(prot.zt0.4) <- zt0$Peak.Name
prot.zt0.5 <- zt0$Swath.O..tauri_01122020_1.2.prima
names(prot.zt0.5) <- zt0$Peak.Name
prot.zt0.6 <- zt0$Swath.O..tauri_01122020_1.3.prima
names(prot.zt0.6) <- zt0$Peak.Name
prot.zt0.7 <- zt0.3$SWATH_O.Tauri_Replica3_12.04.21_ZT20_1
names(prot.zt0.7) <- zt0.3$Peak.Name
prot.zt0.8 <- zt0.3$SWATH_O.Tauri_Replica3_12.04.21_ZT20_2
names(prot.zt0.8) <- zt0.3$Peak.Name
prot.zt0.9 <- zt0.3$SWATH_O.Tauri_Replica3_12.04.21_ZT20_3
names(prot.zt0.9) <- zt0.3$Peak.Name


prot.zt4.1 <- zt4$SWATH_O.tauri_25.11.20_2.1
names(prot.zt4.1) <- zt4$Peak.Name
prot.zt4.2 <- zt4$SWATH_O.tauri_25.11.20_2.2
names(prot.zt4.2) <- zt4$Peak.Name
prot.zt4.3 <- zt4$SWATH_O.tauri_25.11.20_2.3
names(prot.zt4.3) <- zt4$Peak.Name
prot.zt4.4 <- zt4$Swath.O..tauri_01122020_2.1prima
names(prot.zt4.4) <- zt4$Peak.Name
prot.zt4.5 <- zt4$Swath.O..tauri_01122020_2.2prima
names(prot.zt4.5) <- zt4$Peak.Name
prot.zt4.6 <- zt4$Swath.O..tauri_01122020_2.3prima
names(prot.zt4.6) <- zt4$Peak.Name
prot.zt4.7 <- zt4.3$SWATH_O.Tauri_Replica3_12.04.21_ZT4_1
names(prot.zt4.7) <- zt4.3$Peak.Name
prot.zt4.8 <- zt4.3$SWATH_O.Tauri_Replica3_12.04.21_ZT4_2
names(prot.zt4.8) <- zt4.3$Peak.Name
prot.zt4.9 <- zt4.3$SWATH_O.Tauri_Replica3_12.04.21_ZT4_3
names(prot.zt4.9) <- zt4.3$Peak.Name


prot.zt8.1 <- zt8$SWATH_O.tauri_25.11.20_3.1
names(prot.zt8.1) <- zt8$Peak.Name
prot.zt8.2 <- zt8$SWATH_O.tauri_25.11.20_3.2
names(prot.zt8.2) <- zt8$Peak.Name
prot.zt8.3 <- zt8$SWATH_O.tauri_25.11.20_3.3
names(prot.zt8.3) <- zt8$Peak.Name
prot.zt8.4 <- zt8$Swath.O.tauri.cont02122020_3.1prima
names(prot.zt8.4) <- zt8$Peak.Name
prot.zt8.5 <- zt8$Swath.O.tauri.cont02122020_3.2prima
names(prot.zt8.5) <- zt8$Peak.Name
prot.zt8.6 <- zt8$Swath.O.tauri.cont02122020_3.3prima
names(prot.zt8.6) <- zt8$Peak.Name
prot.zt8.7 <- zt8.3$SWATH_O.Tauri_Replica3_12.04.21_ZT8_1
names(prot.zt8.7) <- zt8.3$Peak.Name
prot.zt8.8 <- zt8.3$SWATH_O.Tauri_Replica3_12.04.21_ZT8_2
names(prot.zt8.8) <- zt8.3$Peak.Name
prot.zt8.9 <- zt8.3$SWATH_O.Tauri_Replica3_12.04.21_ZT8_3
names(prot.zt8.9) <- zt8.3$Peak.Name


prot.zt12.1 <- zt12$SWATH_O.tauri_25.11.20_4.1
names(prot.zt12.1) <- zt12$Peak.Name
prot.zt12.2 <- zt12$SWATH_O.tauri_25.11.20_4.2
names(prot.zt12.2) <- zt12$Peak.Name
prot.zt12.3 <- zt12$SWATH_O.tauri_25.11.20_4.3
names(prot.zt12.3) <- zt12$Peak.Name
prot.zt12.4 <- zt12$Swath.O.tauri.cont02122020_4.1prima
names(prot.zt12.4) <- zt12$Peak.Name
prot.zt12.5 <- zt12$Swath.O.tauri.cont02122020_4.2prima
names(prot.zt12.5) <- zt12$Peak.Name
prot.zt12.6 <- zt12$Swath.O.tauri.cont02122020_4.3prima
names(prot.zt12.6) <- zt12$Peak.Name
prot.zt12.7 <- zt12.3$SWATH_O.Tauri_Replica3_12.04.21_ZT12_1
names(prot.zt12.7) <- zt12.3$Peak.Name
prot.zt12.8 <- zt12.3$SWATH_O.Tauri_Replica3_12.04.21_ZT12_2
names(prot.zt12.8) <- zt12.3$Peak.Name
prot.zt12.9 <- zt12.3$SWATH_O.Tauri_Replica3_12.04.21_ZT12_3
names(prot.zt12.9) <- zt12.3$Peak.Name


prot.zt16.1 <- zt16$SWATH_O.tauri_25.11.20_5.1
names(prot.zt16.1) <- zt16$Peak.Name
prot.zt16.2 <- zt16$SWATH_O.tauri_25.11.20_5.2
names(prot.zt16.2) <- zt16$Peak.Name
prot.zt16.3 <- zt16$SWATH_O.tauri_25.11.20_5.3
names(prot.zt16.3) <- zt16$Peak.Name
prot.zt16.4 <- zt20$Swath.O.tauri.cont02122020_6.1prima
names(prot.zt16.4) <- zt20$Peak.Name
prot.zt16.5 <- zt20$SWATH_O.tauri_25.11.20_6.2
names(prot.zt16.5) <- zt20$Peak.Name
prot.zt16.6 <- zt20$Swath.O.tauri.cont02122020_6.3prima
names(prot.zt16.6) <- zt20$Peak.Name
prot.zt16.7 <- zt16.3$SWATH_O.Tauri_Replica3_12.04.21_ZT16_1
names(prot.zt16.7) <- zt16.3$Peak.Name
prot.zt16.8 <- zt16.3$SWATH_O.Tauri_Replica3_12.04.21_ZT16_2
names(prot.zt16.8) <- zt16.3$Peak.Name
prot.zt16.9 <- zt16.3$SWATH_O.Tauri_Replica3_12.04.21_ZT16_3
names(prot.zt16.9) <- zt16.3$Peak.Name


prot.zt20.1 <- zt20$SWATH_O.tauri_25.11.20_6.1
names(prot.zt20.1) <- zt20$Peak.Name
prot.zt20.2 <- zt20$SWATH_O.tauri_25.11.20_6.2
names(prot.zt20.2) <- zt20$Peak.Name
prot.zt20.3 <- zt20$SWATH_O.tauri_25.11.20_6.3
names(prot.zt20.3) <- zt20$Peak.Name
prot.zt20.4 <- zt16$Swath.O.tauri.cont02122020_5.1prima
names(prot.zt20.4) <- zt16$Peak.Name
prot.zt20.5 <- zt16$Swath.O.tauri.cont02122020_5.2prima
names(prot.zt20.5) <- zt16$Peak.Name
prot.zt20.6 <- zt16$Swath.O.tauri.cont02122020_5.3prima
names(prot.zt20.6) <- zt16$Peak.Name
prot.zt20.7 <- zt20.3$SWATH_O.Tauri_Replica3_12.04.21_ZT20_1
names(prot.zt20.7) <- zt20.3$Peak.Name
prot.zt20.8 <- zt20.3$SWATH_O.Tauri_Replica3_12.04.21_ZT20_2
names(prot.zt20.8) <- zt20.3$Peak.Name
prot.zt20.9 <- zt20.3$SWATH_O.Tauri_Replica3_12.04.21_ZT20_3
names(prot.zt20.9) <- zt20.3$Peak.Name




swath.raw.data <- matrix(NA,nrow = length(total.proteins),ncol=54)
rownames(swath.raw.data) <- total.proteins
colnames(swath.raw.data) <- c(paste("zt0",1:9,sep="_"),paste("zt4",1:9,sep="_"),
                              paste("zt8",1:9,sep="_"),paste("zt12",1:9,sep="_"),
                              paste("zt16",1:9,sep="_"),paste("zt20",1:9,sep="_"))

swath.raw.data[names(prot.zt0.1),"zt0_1"] <- prot.zt0.1
swath.raw.data[names(prot.zt0.2),"zt0_2"] <- prot.zt0.2
swath.raw.data[names(prot.zt0.3),"zt0_3"] <- prot.zt0.3
swath.raw.data[names(prot.zt0.4),"zt0_4"] <- prot.zt0.4
swath.raw.data[names(prot.zt0.5),"zt0_5"] <- prot.zt0.5
swath.raw.data[names(prot.zt0.6),"zt0_6"] <- prot.zt0.6
swath.raw.data[names(prot.zt0.7),"zt0_7"] <- prot.zt0.7
swath.raw.data[names(prot.zt0.8),"zt0_8"] <- prot.zt0.8
swath.raw.data[names(prot.zt0.9),"zt0_9"] <- prot.zt0.9

swath.raw.data[names(prot.zt4.1),"zt4_1"] <- prot.zt4.1
swath.raw.data[names(prot.zt4.2),"zt4_2"] <- prot.zt4.2
swath.raw.data[names(prot.zt4.3),"zt4_3"] <- prot.zt4.3
swath.raw.data[names(prot.zt4.4),"zt4_4"] <- prot.zt4.4
swath.raw.data[names(prot.zt4.5),"zt4_5"] <- prot.zt4.5
swath.raw.data[names(prot.zt4.6),"zt4_6"] <- prot.zt4.6
swath.raw.data[names(prot.zt4.7),"zt4_7"] <- prot.zt4.7
swath.raw.data[names(prot.zt4.8),"zt4_8"] <- prot.zt4.8
swath.raw.data[names(prot.zt4.9),"zt4_9"] <- prot.zt4.9

swath.raw.data[names(prot.zt8.1),"zt8_1"] <- prot.zt8.1
swath.raw.data[names(prot.zt8.2),"zt8_2"] <- prot.zt8.2
swath.raw.data[names(prot.zt8.3),"zt8_3"] <- prot.zt8.3
swath.raw.data[names(prot.zt8.4),"zt8_4"] <- prot.zt8.4
swath.raw.data[names(prot.zt8.5),"zt8_5"] <- prot.zt8.5
swath.raw.data[names(prot.zt8.6),"zt8_6"] <- prot.zt8.6
swath.raw.data[names(prot.zt8.7),"zt8_7"] <- prot.zt8.7
swath.raw.data[names(prot.zt8.8),"zt8_8"] <- prot.zt8.8
swath.raw.data[names(prot.zt8.9),"zt8_9"] <- prot.zt8.9

swath.raw.data[names(prot.zt12.1),"zt12_1"] <- prot.zt12.1
swath.raw.data[names(prot.zt12.2),"zt12_2"] <- prot.zt12.2
swath.raw.data[names(prot.zt12.3),"zt12_3"] <- prot.zt12.3
swath.raw.data[names(prot.zt12.4),"zt12_4"] <- prot.zt12.4
swath.raw.data[names(prot.zt12.5),"zt12_5"] <- prot.zt12.5
swath.raw.data[names(prot.zt12.6),"zt12_6"] <- prot.zt12.6
swath.raw.data[names(prot.zt12.7),"zt12_7"] <- prot.zt12.7
swath.raw.data[names(prot.zt12.8),"zt12_8"] <- prot.zt12.8
swath.raw.data[names(prot.zt12.9),"zt12_9"] <- prot.zt12.9

swath.raw.data[names(prot.zt16.1),"zt16_1"] <- prot.zt16.1
swath.raw.data[names(prot.zt16.2),"zt16_2"] <- prot.zt16.2
swath.raw.data[names(prot.zt16.3),"zt16_3"] <- prot.zt16.3
swath.raw.data[names(prot.zt16.4),"zt16_4"] <- prot.zt16.4
swath.raw.data[names(prot.zt16.5),"zt16_5"] <- prot.zt16.5
swath.raw.data[names(prot.zt16.6),"zt16_6"] <- prot.zt16.6
swath.raw.data[names(prot.zt16.7),"zt16_7"] <- prot.zt16.7
swath.raw.data[names(prot.zt16.8),"zt16_8"] <- prot.zt16.8
swath.raw.data[names(prot.zt16.9),"zt16_9"] <- prot.zt16.9

swath.raw.data[names(prot.zt20.1),"zt20_1"] <- prot.zt20.1
swath.raw.data[names(prot.zt20.2),"zt20_2"] <- prot.zt20.2
swath.raw.data[names(prot.zt20.3),"zt20_3"] <- prot.zt20.3
swath.raw.data[names(prot.zt20.4),"zt20_4"] <- prot.zt20.4
swath.raw.data[names(prot.zt20.5),"zt20_5"] <- prot.zt20.5
swath.raw.data[names(prot.zt20.6),"zt20_6"] <- prot.zt20.6
swath.raw.data[names(prot.zt20.7),"zt20_7"] <- prot.zt20.7
swath.raw.data[names(prot.zt20.8),"zt20_8"] <- prot.zt20.8
swath.raw.data[names(prot.zt20.9),"zt20_9"] <- prot.zt20.9

png(filename = "boxplot_before_normalization.png",width = 1000)
boxplot(swath.raw.data,outline=F,las=2,col=rep(rainbow(6),each=9),
        main="Before Normalization",cex.main=2)
dev.off()

swath.raw.data.frame <- data.frame(rownames(swath.raw.data),swath.raw.data)
head(swath.raw.data.frame)
colnames(swath.raw.data.frame)[1] <- "ProtID"

sorted.prots <- sort.int(x = swath.raw.data.frame$ProtID,index.return = T)

swath.raw.data.frame <- swath.raw.data.frame[sorted.prots$ix,]
head(swath.raw.data.frame)
swath.raw.data.frame <- swath.raw.data.frame[-(1:4),]
head(swath.raw.data.frame)

rev(rownames(swath.raw.data.frame))

swath.raw.data.frame <- swath.raw.data.frame[-nrow(swath.raw.data.frame),]
write.table(x = swath.raw.data.frame,file = "LD_swath_raw_data.tsv",
            quote = F,row.names = F,
            sep = "\t")

for(i in 1:nrow(ld.normalized.proteomic.data))
{
 if(sum(is.na(ld.normalized.proteomic.data[i,])) != 0)
 {
  na.points <- colnames(ld.normalized.proteomic.data)[which(is.na(ld.normalized.proteomic.data[i,]))]
  
  zts <- unique(sapply(X = strsplit(na.points,split="_"), FUN = function(x){ return(x[1])}))
  
  for(j in 1:length(zts))
  {
   current.time.point <- as.numeric(strsplit(zts[j],split="t")[[1]][2])
   
   if(current.time.point > 0)
   {
    previous.time.point <- paste("zt",current.time.point-4,sep="")
   } else
   {
    previous.time.point <- "zt20"
   }
   
   if(current.time.point < 20)
   {
    next.time.point <- paste("zt",current.time.point+4,sep="")
   } else
   {
    next.time.point <- "zt0"
   }
   
   ld.normalized.proteomic.data[i,paste(zts[j],1:9,sep="_")] <-  mean(as.numeric(c(ld.normalized.proteomic.data[i,paste(previous.time.point,1:9,sep="_")],
   ld.normalized.proteomic.data[i,paste(next.time.point,1:9,sep="_")])))
  }
 }
}


## SD proteomic data
zt0.1 <- read.table(file="/home/fran/Nextcloud2/Microalgas/analysis_minotaur/swath/sd_1_2_3/zt00_1.tsv",header = T,as.is=T,sep="\t")
zt0.1 <- zt0.1[,c(1,6:8)]
zt0.2 <- read.table(file="/home/fran/Nextcloud2/Microalgas/analysis_minotaur/swath/sd_1_2_3/zt00_2.tsv",header = T,as.is=T,sep="\t")
zt0.2 <- zt0.2[,c(1,6:8)]
zt0.3 <- read.table(file="/home/fran/Nextcloud2/Microalgas/analysis_minotaur/swath/sd_1_2_3/zt00_3.tsv",header = T,as.is=T,sep="\t")
zt0.3 <- zt0.3[,c(1,6:8)]

zt4.1 <- read.table(file="/home/fran/Nextcloud2/Microalgas/analysis_minotaur/swath/sd_1_2_3/zt04_1.tsv",header = T,as.is=T,sep="\t")
zt4.1 <- zt4.1[,c(1,6:8)]
zt4.2 <- read.table(file="/home/fran/Nextcloud2/Microalgas/analysis_minotaur/swath/sd_1_2_3/zt04_2.tsv",header = T,as.is=T,sep="\t")
zt4.2 <- zt4.2[,c(1,6:8)]
zt4.3 <- read.table(file="/home/fran/Nextcloud2/Microalgas/analysis_minotaur/swath/sd_1_2_3/zt04_3.tsv",header = T,as.is=T,sep="\t")
zt4.3 <- zt4.3[,c(1,6:8)]

zt8.1 <- read.table(file="/home/fran/Nextcloud2/Microalgas/analysis_minotaur/swath/sd_1_2_3/zt08_1.tsv",header = T,as.is=T,sep="\t")
zt8.1 <- zt8.1[,c(1,6:8)]
zt8.2 <- read.table(file="/home/fran/Nextcloud2/Microalgas/analysis_minotaur/swath/sd_1_2_3/zt08_2.tsv",header = T,as.is=T,sep="\t")
zt8.2 <- zt8.2[,c(1,6:8)]
zt8.3 <- read.table(file="/home/fran/Nextcloud2/Microalgas/analysis_minotaur/swath/sd_1_2_3/zt08_3.tsv",header = T,as.is=T,sep="\t")
zt8.3 <- zt8.3[,c(1,6:8)]

zt12.1 <- read.table(file="/home/fran/Nextcloud2/Microalgas/analysis_minotaur/swath/sd_1_2_3/zt12_1.tsv",header = T,as.is=T,sep="\t")
zt12.1 <- zt12.1[,c(1,6:8)]
zt12.2 <- read.table(file="/home/fran/Nextcloud2/Microalgas/analysis_minotaur/swath/sd_1_2_3/zt12_2.tsv",header = T,as.is=T,sep="\t")
zt12.2 <- zt12.2[,c(1,6:8)]
zt12.3 <- read.table(file="/home/fran/Nextcloud2/Microalgas/analysis_minotaur/swath/sd_1_2_3/zt12_3.tsv",header = T,as.is=T,sep="\t")
zt12.3 <- zt12.3[,c(1,6:8)]

zt16.1 <- read.table(file="/home/fran/Nextcloud2/Microalgas/analysis_minotaur/swath/sd_1_2_3/zt16_1.tsv",header = T,as.is=T,sep="\t")
zt16.1 <- zt16.1[,c(1,6:8)]
zt16.2 <- read.table(file="/home/fran/Nextcloud2/Microalgas/analysis_minotaur/swath/sd_1_2_3/zt16_2.tsv",header = T,as.is=T,sep="\t")
zt16.2 <- zt16.2[,c(1,6:8)]
zt16.3 <- read.table(file="/home/fran/Nextcloud2/Microalgas/analysis_minotaur/swath/sd_1_2_3/zt16_3.tsv",header = T,as.is=T,sep="\t")
zt16.3 <- zt16.3[,c(1,6:8)]

zt20.1 <- read.table(file="/home/fran/Nextcloud2/Microalgas/analysis_minotaur/swath/sd_1_2_3/zt20_1.tsv",header = T,as.is=T,sep="\t")
zt20.1 <- zt20.1[,c(1,6:8)]
zt20.2 <- read.table(file="/home/fran/Nextcloud2/Microalgas/analysis_minotaur/swath/sd_1_2_3/zt20_2.tsv",header = T,as.is=T,sep="\t")
zt20.2 <- zt20.2[,c(1,6:8)]
zt20.3 <- read.table(file="/home/fran/Nextcloud2/Microalgas/analysis_minotaur/swath/sd_1_2_3/zt20_3.tsv",header = T,as.is=T,sep="\t")
zt20.3 <- zt20.3[,c(1,6:8)]

total.proteins <- unique(c(zt0.1$Peak.Name,zt0.2$Peak.Name,zt0.3$Peak.Name,
                           zt4.1$Peak.Name,zt4.2$Peak.Name,zt4.3$Peak.Name, 
                           zt8.1$Peak.Name,zt8.2$Peak.Name,zt8.3$Peak.Name,
                           zt12.1$Peak.Name,zt12.2$Peak.Name,zt12.3$Peak.Name,
                           zt16.1$Peak.Name,zt16.2$Peak.Name,zt16.3$Peak.Name,
                           zt20.1$Peak.Name,zt20.2$Peak.Name,zt20.3$Peak.Name))
length(total.proteins)

prot.zt0.1 <- zt0.1$SWATH_O.tauri_C.corto_03.05.21_ZT0_1
names(prot.zt0.1) <- zt0.1$Peak.Name
prot.zt0.2 <- zt0.1$SWATH_O.tauri_C.corto_03.05.21_ZT0_2
names(prot.zt0.2) <- zt0.1$Peak.Name
prot.zt0.3 <- zt0.1$SWATH_O.tauri_C.corto_03.05.21_ZT0_3
names(prot.zt0.3) <- zt0.1$Peak.Name
prot.zt0.4 <- zt0.2$SWATH_O.tauri_Ciclo_Corto_Replica_2_10.05.21_ZT0_1
names(prot.zt0.4) <- zt0.2$Peak.Name
prot.zt0.5 <- zt0.2$SWATH_O.tauri_Ciclo_Corto_Replica_2_10.05.21_ZT0_2
names(prot.zt0.5) <- zt0.2$Peak.Name
prot.zt0.6 <- zt0.2$SWATH_O.tauri_Ciclo_Corto_Replica_2_10.05.21_ZT0_3
names(prot.zt0.6) <- zt0.2$Peak.Name
prot.zt0.7 <- zt0.3$SWATH_O.Tauri_Replica3_12.04.21_ZT0_1
names(prot.zt0.7) <- zt0.3$Peak.Name
prot.zt0.8 <- zt0.3$SWATH_O.Tauri_Replica3_12.04.21_ZT0_2
names(prot.zt0.8) <- zt0.3$Peak.Name
prot.zt0.9 <- zt0.3$SWATH_O.Tauri_Replica3_12.04.21_ZT0_3
names(prot.zt0.9) <- zt0.3$Peak.Name

prot.zt4.1 <- zt4.1$SWATH_O.tauri_C.corto_03.05.21_ZT4_1
names(prot.zt4.1) <- zt4.1$Peak.Name
prot.zt4.2 <- zt4.1$SWATH_O.tauri_C.corto_03.05.21_ZT4_2
names(prot.zt4.2) <- zt4.1$Peak.Name
prot.zt4.3 <- zt4.1$SWATH_O.tauri_C.corto_03.05.21_ZT4_3
names(prot.zt4.3) <- zt4.1$Peak.Name
prot.zt4.4 <- zt4.2$SWATH_O.tauri_Ciclo_Corto_Replica_2_10.05.21_ZT4_1
names(prot.zt4.4) <- zt4.2$Peak.Name
prot.zt4.5 <- zt4.2$SWATH_O.tauri_Ciclo_Corto_Replica_2_10.05.21_ZT4_2
names(prot.zt4.5) <- zt4.2$Peak.Name
prot.zt4.6 <- zt4.2$SWATH_O.tauri_Ciclo_Corto_Replica_2_10.05.21_ZT4_3
names(prot.zt4.6) <- zt4.2$Peak.Name
prot.zt4.7 <- zt4.3$SWATH_O.tauri_Ciclo_Corto_Replica_3_17.05.21_ZT4_1
names(prot.zt4.7) <- zt4.3$Peak.Name
prot.zt4.8 <- zt4.3$SWATH_O.tauri_Ciclo_Corto_Replica_3_17.05.21_ZT4_2
names(prot.zt4.8) <- zt4.3$Peak.Name
prot.zt4.9 <- zt4.3$SWATH_O.tauri_Ciclo_Corto_Replica_3_17.05.21_ZT4_3
names(prot.zt4.9) <- zt4.3$Peak.Name

prot.zt8.1 <- zt8.1$SWATH_O.tauri_C.corto_03.05.21_ZT8_1
names(prot.zt8.1) <- zt8.1$Peak.Name
prot.zt8.2 <- zt8.1$SWATH_O.tauri_C.corto_03.05.21_ZT8_2
names(prot.zt8.2) <- zt8.1$Peak.Name
prot.zt8.3 <- zt8.1$SWATH_O.tauri_C.corto_03.05.21_ZT8_3
names(prot.zt8.3) <- zt8.1$Peak.Name
prot.zt8.4 <- zt8.2$SWATH_O.tauri_Ciclo_Corto_Replica_2_10.05.21_ZT8_1
names(prot.zt8.4) <- zt8.2$Peak.Name
prot.zt8.5 <- zt8.2$SWATH_O.tauri_Ciclo_Corto_Replica_2_10.05.21_ZT8_2
names(prot.zt8.5) <- zt8.2$Peak.Name
prot.zt8.6 <- zt8.2$SWATH_O.tauri_Ciclo_Corto_Replica_2_10.05.21_ZT8_3
names(prot.zt8.6) <- zt8.2$Peak.Name
prot.zt8.7 <- zt8.3$SWATH_O.Tauri_Replica3_12.04.21_ZT8_1
names(prot.zt8.7) <- zt8.3$Peak.Name
prot.zt8.8 <- zt8.3$SWATH_O.Tauri_Replica3_12.04.21_ZT8_2
names(prot.zt8.8) <- zt8.3$Peak.Name
prot.zt8.9 <- zt8.3$SWATH_O.Tauri_Replica3_12.04.21_ZT8_3
names(prot.zt8.9) <- zt8.3$Peak.Name

prot.zt12.1 <- zt12.1$SWATH_O.tauri_C.corto_03.05.21_ZT12_1
names(prot.zt12.1) <- zt12.1$Peak.Name
prot.zt12.2 <- zt12.1$SWATH_O.tauri_C.corto_03.05.21_ZT12_2
names(prot.zt12.2) <- zt12.1$Peak.Name
prot.zt12.3 <- zt12.1$SWATH_O.tauri_C.corto_03.05.21_ZT12_3
names(prot.zt12.3) <- zt12.1$Peak.Name
prot.zt12.4 <- zt12.2$SWATH_O.tauri_Ciclo_Corto_Replica_2_10.05.21_ZT12_1
names(prot.zt12.4) <- zt12.2$Peak.Name
prot.zt12.5 <- zt12.2$SWATH_O.tauri_Ciclo_Corto_Replica_2_10.05.21_ZT12_2
names(prot.zt12.5) <- zt12.2$Peak.Name
prot.zt12.6 <- zt12.2$SWATH_O.tauri_Ciclo_Corto_Replica_2_10.05.21_ZT12_3
names(prot.zt12.6) <- zt12.2$Peak.Name
prot.zt12.7 <- zt12.3$SWATH_O.Tauri_Replica3_12.04.21_ZT12_1
names(prot.zt12.7) <- zt12.3$Peak.Name
prot.zt12.8 <- zt12.3$SWATH_O.Tauri_Replica3_12.04.21_ZT12_2
names(prot.zt12.8) <- zt12.3$Peak.Name
prot.zt12.9 <- zt12.3$SWATH_O.Tauri_Replica3_12.04.21_ZT12_3
names(prot.zt12.9) <- zt12.3$Peak.Name

prot.zt16.1 <- zt16.1$SWATH_O.tauri_C.corto_03.05.21_ZT16_1
names(prot.zt16.1) <- zt16.1$Peak.Name
prot.zt16.2 <- zt16.1$SWATH_O.tauri_C.corto_03.05.21_ZT16_2
names(prot.zt16.2) <- zt16.1$Peak.Name
prot.zt16.3 <- zt16.1$SWATH_O.tauri_C.corto_03.05.21_ZT16_3
names(prot.zt16.3) <- zt16.1$Peak.Name
prot.zt16.4 <- zt16.2$SWATH_O.tauri_Ciclo_Corto_Replica_2_10.05.21_ZT16_1
names(prot.zt16.4) <- zt16.2$Peak.Name
prot.zt16.5 <- zt16.2$SWATH_O.tauri_Ciclo_Corto_Replica_2_10.05.21_ZT16_2
names(prot.zt16.5) <- zt16.2$Peak.Name
prot.zt16.6 <- zt16.2$SWATH_O.tauri_Ciclo_Corto_Replica_2_10.05.21_ZT16_3
names(prot.zt16.6) <- zt16.2$Peak.Name
prot.zt16.7 <- zt16.3$SWATH_O.Tauri_Replica3_12.04.21_ZT16_1
names(prot.zt16.7) <- zt16.3$Peak.Name
prot.zt16.8 <- zt16.3$SWATH_O.Tauri_Replica3_12.04.21_ZT16_2
names(prot.zt16.8) <- zt16.3$Peak.Name
prot.zt16.9 <- zt16.3$SWATH_O.Tauri_Replica3_12.04.21_ZT16_3
names(prot.zt16.9) <- zt16.3$Peak.Name

prot.zt20.1 <- zt20.1$SWATH_O.tauri_C.corto_03.05.21_ZT20_1
names(prot.zt20.1) <- zt20.1$Peak.Name
prot.zt20.2 <- zt20.1$SWATH_O.tauri_C.corto_03.05.21_ZT20_2
names(prot.zt20.2) <- zt20.1$Peak.Name
prot.zt20.3 <- zt20.1$SWATH_O.tauri_C.corto_03.05.21_ZT20_3
names(prot.zt20.3) <- zt20.1$Peak.Name
prot.zt20.4 <- zt20.2$SWATH_O.tauri_Ciclo_Corto_Replica_2_10.05.21_ZT20_1
names(prot.zt20.4) <- zt20.2$Peak.Name
prot.zt20.5 <- zt20.2$SWATH_O.tauri_Ciclo_Corto_Replica_2_10.05.21_ZT20_2
names(prot.zt20.5) <- zt20.2$Peak.Name
prot.zt20.6 <- zt20.2$SWATH_O.tauri_Ciclo_Corto_Replica_2_10.05.21_ZT20_3
names(prot.zt20.6) <- zt20.2$Peak.Name
prot.zt20.7 <- zt20.3$SWATH_O.tauri_Ciclo_Corto_Replica_3_17.05.21_ZT20_1
names(prot.zt20.7) <- zt20.3$Peak.Name
prot.zt20.8 <- zt20.3$SWATH_O.tauri_Ciclo_Corto_Replica_3_17.05.21_ZT20_2
names(prot.zt20.8) <- zt20.3$Peak.Name
prot.zt20.9 <- zt20.3$SWATH_O.tauri_Ciclo_Corto_Replica_3_17.05.21_ZT20_3
names(prot.zt20.9) <- zt20.3$Peak.Name

swath.raw.data <- matrix(NA,nrow = length(total.proteins),ncol=54)
rownames(swath.raw.data) <- total.proteins
colnames(swath.raw.data) <- c(paste("zt0",1:9,sep="_"),paste("zt4",1:9,sep="_"),
                              paste("zt8",1:9,sep="_"),paste("zt12",1:9,sep="_"),
                              paste("zt16",1:9,sep="_"),paste("zt20",1:9,sep="_"))

swath.raw.data[names(prot.zt0.1),"zt0_1"] <- prot.zt0.1
swath.raw.data[names(prot.zt0.2),"zt0_2"] <- prot.zt0.2
swath.raw.data[names(prot.zt0.3),"zt0_3"] <- prot.zt0.3
swath.raw.data[names(prot.zt0.4),"zt0_4"] <- prot.zt0.4
swath.raw.data[names(prot.zt0.5),"zt0_5"] <- prot.zt0.5
swath.raw.data[names(prot.zt0.6),"zt0_6"] <- prot.zt0.6
swath.raw.data[names(prot.zt0.7),"zt0_7"] <- prot.zt0.7
swath.raw.data[names(prot.zt0.8),"zt0_8"] <- prot.zt0.8
swath.raw.data[names(prot.zt0.9),"zt0_9"] <- prot.zt0.9

swath.raw.data[names(prot.zt4.1),"zt4_1"] <- prot.zt4.1
swath.raw.data[names(prot.zt4.2),"zt4_2"] <- prot.zt4.2
swath.raw.data[names(prot.zt4.3),"zt4_3"] <- prot.zt4.3
swath.raw.data[names(prot.zt4.4),"zt4_4"] <- prot.zt4.4
swath.raw.data[names(prot.zt4.5),"zt4_5"] <- prot.zt4.5
swath.raw.data[names(prot.zt4.6),"zt4_6"] <- prot.zt4.6
swath.raw.data[names(prot.zt4.7),"zt4_7"] <- prot.zt4.7
swath.raw.data[names(prot.zt4.8),"zt4_8"] <- prot.zt4.8
swath.raw.data[names(prot.zt4.9),"zt4_9"] <- prot.zt4.9

swath.raw.data[names(prot.zt8.1),"zt8_1"] <- prot.zt8.1
swath.raw.data[names(prot.zt8.2),"zt8_2"] <- prot.zt8.2
swath.raw.data[names(prot.zt8.3),"zt8_3"] <- prot.zt8.3
swath.raw.data[names(prot.zt8.4),"zt8_4"] <- prot.zt8.4
swath.raw.data[names(prot.zt8.5),"zt8_5"] <- prot.zt8.5
swath.raw.data[names(prot.zt8.6),"zt8_6"] <- prot.zt8.6
swath.raw.data[names(prot.zt8.7),"zt8_7"] <- prot.zt8.7
swath.raw.data[names(prot.zt8.8),"zt8_8"] <- prot.zt8.8
swath.raw.data[names(prot.zt8.9),"zt8_9"] <- prot.zt8.9

swath.raw.data[names(prot.zt12.1),"zt12_1"] <- prot.zt12.1
swath.raw.data[names(prot.zt12.2),"zt12_2"] <- prot.zt12.2
swath.raw.data[names(prot.zt12.3),"zt12_3"] <- prot.zt12.3
swath.raw.data[names(prot.zt12.4),"zt12_4"] <- prot.zt12.4
swath.raw.data[names(prot.zt12.5),"zt12_5"] <- prot.zt12.5
swath.raw.data[names(prot.zt12.6),"zt12_6"] <- prot.zt12.6
swath.raw.data[names(prot.zt12.7),"zt12_7"] <- prot.zt12.7
swath.raw.data[names(prot.zt12.8),"zt12_8"] <- prot.zt12.8
swath.raw.data[names(prot.zt12.9),"zt12_9"] <- prot.zt12.9

swath.raw.data[names(prot.zt16.1),"zt16_1"] <- prot.zt16.1
swath.raw.data[names(prot.zt16.2),"zt16_2"] <- prot.zt16.2
swath.raw.data[names(prot.zt16.3),"zt16_3"] <- prot.zt16.3
swath.raw.data[names(prot.zt16.4),"zt16_4"] <- prot.zt16.4
swath.raw.data[names(prot.zt16.5),"zt16_5"] <- prot.zt16.5
swath.raw.data[names(prot.zt16.6),"zt16_6"] <- prot.zt16.6
swath.raw.data[names(prot.zt16.7),"zt16_7"] <- prot.zt16.7
swath.raw.data[names(prot.zt16.8),"zt16_8"] <- prot.zt16.8
swath.raw.data[names(prot.zt16.9),"zt16_9"] <- prot.zt16.9

swath.raw.data[names(prot.zt20.1),"zt20_1"] <- prot.zt20.1
swath.raw.data[names(prot.zt20.2),"zt20_2"] <- prot.zt20.2
swath.raw.data[names(prot.zt20.3),"zt20_3"] <- prot.zt20.3
swath.raw.data[names(prot.zt20.4),"zt20_4"] <- prot.zt20.4
swath.raw.data[names(prot.zt20.5),"zt20_5"] <- prot.zt20.5
swath.raw.data[names(prot.zt20.6),"zt20_6"] <- prot.zt20.6
swath.raw.data[names(prot.zt20.7),"zt20_7"] <- prot.zt20.7
swath.raw.data[names(prot.zt20.8),"zt20_8"] <- prot.zt20.8
swath.raw.data[names(prot.zt20.9),"zt20_9"] <- prot.zt20.9

boxplot(swath.raw.data,outline=F,las=2,col=rep(rainbow(6),each=9),
        main="Before Normalization",cex.main=2)

swath.raw.data.frame <- data.frame(rownames(swath.raw.data),swath.raw.data)
head(swath.raw.data.frame)
colnames(swath.raw.data.frame)[1] <- "ProtID"
write.table(x = swath.raw.data.frame,file = "swath_proteomic_data/SD_swath_raw_data.tsv",
            quote = F,row.names = F,
            sep = "\t")

ld.protein.abundance.1 <- 2^ld.protein.abundance

sd.ld <- apply(X = ld.protein.abundance.1,MARGIN = 1,FUN = sd)
names(sort(sd.ld)[1:10])

#"ostta03g01680" "ostta04g05460" "ostta02g04910" "ostta12g02510" "ostta06g00135"
#[6] "ostta16g02045" "ostta04g03880" "ostta03g00120" "ostta15g01560" "ostta07g01600"

current.gene <- "ostta07g01600"
plot.ld.gene.prot(gene.id=current.gene, gene.name=current.gene, gene.expression=gene.expression, protein.data=ld.protein.abundance)
plot.sd.gene.prot(gene.id=current.gene, gene.name=current.gene, gene.expression=gene.expression, protein.data=sd.protein.abundance)



genes <- c("ostta18g01570", "ostta01g06150", "ostta18g01420")
gene.names <- c("CycA", "CycB", "CycD")
entrainment <- "SD"
free.running <- "DD"

plot.multiple.genes.expression.profile <- function(genes, gene.names, gene.expression,
                                    entrainment,free.running,param=F)
{
 if (entrainment == "SD")
 {
  zts <- paste("sd",paste0("zt",sprintf(fmt = "%02d",seq(from=0,to=20,by=4))),sep="_")
  entrainment.col <- "red"
  free.running.col <- "lightsalmon"
 } else if (entrainment == "LD")
 {
  zts <- paste("ld",paste0("zt",sprintf(fmt = "%02d",seq(from=0,to=20,by=4))),sep="_")
  entrainment.col <- "blue"
  free.running.col <- "lightblue"
 } else
 {
  print("Unkwnon entrainment. Possible entrainments SD and LD")
 }
 
 if (free.running == "LL")
 {
  free.running.indeces <- c(4,5)
 } else if (free.running == "DD")
 {
  free.running.indeces <- c(6,7)
 } else
 {
  print("Unknown free running. Possible options LL and DD")
 }
 
 current.genes.expression.entrainment <- gene.expression[genes,c(paste(zts,1,sep="_"),
                                                                  paste(zts,2,sep="_"),
                                                                  paste(zts,3,sep="_"))]
 current.genes.expression.free.running <- gene.expression[genes,c(paste(zts,free.running.indeces[1],sep="_"),
                                                                   paste(zts,free.running.indeces[2],sep="_"))]
 current.genes.expression <- cbind(current.genes.expression.entrainment,current.genes.expression.free.running)
 min.expression <- min(current.genes.expression)
 max.expression <- max(current.genes.expression)
 range.expression <- max.expression - min.expression
 
 expression.step <- floor(range.expression / 5)
 time.points.entrainment <- seq(from=0,by=4,length.out = 18)
 time.points.free.running <- seq(from=0,by=4,length.out = 12)
 
 
 if(param)
 {
  circacompare.data <- data.frame(time=c(time.points.entrainment,time.points.free.running),
                                  measure=c(current.gene.expression.entrainment,
                                            current.gene.expression.free.running),
                                  group=c(rep("Entrainment",18),rep("Free_running",12)))
  
  result.i <- circacompare(x = circacompare.data, 
                           col_time = "time", 
                           col_group = "group", 
                           col_outcome = "measure",
                           alpha_threshold = 1)
  result.i.values <- result.i$summary$value
  names(result.i.values) <- result.i$summary$parameter
  time.entrainment <- seq(from=0,to=72,by=0.01)
  time.free.running <- seq(from=72,to=120,by=0.01)
  line.entrainment <- wave.form(mesor = result.i.values["Entrainment mesor estimate"],
                                amplitude = result.i.values["Entrainment amplitude estimate"],
                                period = result.i.values["Shared period estimate"],phase = result.i.values["Entrainment peak time hours"],
                                time = time.entrainment)
  line.free.running <- wave.form(mesor = result.i.values["Free_running mesor estimate"],
                                 amplitude = result.i.values["Free_running amplitude estimate"],
                                 period = result.i.values["Shared period estimate"],
                                 phase = result.i.values["Free_running peak time hours"],
                                 time = time.free.running)
  plot(x = c(time.points.entrainment,time.points.free.running+max(time.points.sd)+4),#type="o",
       y = current.gene.expression,lwd=1,col=entrainment.col,axes=F,xlab="",ylab="FPKM",
       ylim=c(min.expression-expression.step,max.expression),
       cex.lab=1.3,main=paste(gene.id, gene.name,sep=" - "),cex.main=2,pch=20,cex=0.8)
  lines(x = c(time.entrainment,time.free.running), 
        y = c(line.entrainment,line.free.running), 
        type="l",
        col=entrainment.col,lwd=3)
 } else
 {
  plot(x= c(time.points.entrainment, 72 +time.points.free.running),
       y = current.genes.expression[1,],type="l",lwd=3,col=entrainment.col,axes=F,xlab="",ylab="FPKM",
       ylim=c(min.expression-expression.step,max.expression),
       cex.lab=1.3,main=paste(gene.id, gene.name,sep=" - "),cex.main=2,pch=0)
  if(length(genes) > 1)
  {
   for(i in 2:length(genes))
   {
    lines(x = c(time.points.entrainment, 72 +time.points.free.running), 
          y = current.genes.expression[i,],type="l",lwd=3,col=entrainment.col, lty=i, pch=(i-1) )
   }
   
  }
  
   
   
 }
 
 axis(side=2,lwd=3)
 axis(side = 1,pos = min.expression - 1.1*expression.step, at = seq(from=0,to=116,by=4),
      labels = rep(paste("ZT",seq(from=0,to=20,by=4)),5),las=2,lwd=3)
 if(entrainment == "SD")
 {
  polygon(x = c(0,8,8,0),y=c(min.expression-expression.step/2,
                             min.expression-expression.step/2,
                             min.expression-expression.step,
                             min.expression-expression.step),lwd=2,border="red")
  polygon(x = c(8,24,24,8),y=c(min.expression-expression.step/2,
                               min.expression-expression.step/2,
                               min.expression-expression.step,
                               min.expression-expression.step),lwd=2,border="red",col="red")
  
  polygon(x = c(24,32,32,24),y=c(min.expression-expression.step/2,
                                 min.expression-expression.step/2,
                                 min.expression-expression.step,
                                 min.expression-expression.step),lwd=2,border="red")
  polygon(x = c(32,48,48,32),y=c(min.expression-expression.step/2,
                                 min.expression-expression.step/2,
                                 min.expression-expression.step,
                                 min.expression-expression.step),lwd=2,border="red",col="red")
  
  polygon(x = c(48,56,56,48),y=c(min.expression-expression.step/2,
                                 min.expression-expression.step/2,
                                 min.expression-expression.step,
                                 min.expression-expression.step),lwd=2,border="red")
  polygon(x = c(56,72,72,56),y=c(min.expression-expression.step/2,
                                 min.expression-expression.step/2,
                                 min.expression-expression.step,
                                 min.expression-expression.step),lwd=2,border="red",col="red")
  if(free.running == "LL")
  {
   polygon(x = c(72,80,80,72),y=c(min.expression-expression.step/2,
                                  min.expression-expression.step/2,
                                  min.expression-expression.step,
                                  min.expression-expression.step),lwd=2,border="red")
   polygon(x = c(80,96,96,80),y=c(min.expression-expression.step/2,
                                  min.expression-expression.step/2,
                                  min.expression-expression.step,
                                  min.expression-expression.step),lwd=2,border="red",col="salmon")
   
   polygon(x = c(96,104,104,96),y=c(min.expression-expression.step/2,
                                    min.expression-expression.step/2,
                                    min.expression-expression.step,
                                    min.expression-expression.step),lwd=2,border="red")
   polygon(x = c(104,120,120,104),y=c(min.expression-expression.step/2,
                                      min.expression-expression.step/2,
                                      min.expression-expression.step,
                                      min.expression-expression.step),lwd=2,border="red",col="salmon")
  } else if (free.running == "DD")
  {
   polygon(x = c(72,80,80,72),y=c(min.expression-expression.step/2,
                                  min.expression-expression.step/2,
                                  min.expression-expression.step,
                                  min.expression-expression.step),lwd=2,border="red",col="salmon")
   polygon(x = c(80,96,96,80),y=c(min.expression-expression.step/2,
                                  min.expression-expression.step/2,
                                  min.expression-expression.step,
                                  min.expression-expression.step),lwd=2,border="red",col="red")
   
   polygon(x = c(96,104,104,96),y=c(min.expression-expression.step/2,
                                    min.expression-expression.step/2,
                                    min.expression-expression.step,
                                    min.expression-expression.step),lwd=2,border="red",col="salmon")
   polygon(x = c(104,120,120,104),y=c(min.expression-expression.step/2,
                                      min.expression-expression.step/2,
                                      min.expression-expression.step,
                                      min.expression-expression.step),lwd=2,border="red",col="red")
  }
 } else if (entrainment == "LD")
 {
  polygon(x = c(0,16,16,0),y=c(min.expression-expression.step/2,
                               min.expression-expression.step/2,
                               min.expression-expression.step,
                               min.expression-expression.step),lwd=2,border="blue")
  polygon(x = c(16,24,24,16),y=c(min.expression-expression.step/2,
                                 min.expression-expression.step/2,
                                 min.expression-expression.step,
                                 min.expression-expression.step),lwd=2,border="blue",col="blue")
  
  polygon(x = c(24,40,40,24),y=c(min.expression-expression.step/2,
                                 min.expression-expression.step/2,
                                 min.expression-expression.step,
                                 min.expression-expression.step),lwd=2,border="blue")
  polygon(x = c(40,48,48,40),y=c(min.expression-expression.step/2,
                                 min.expression-expression.step/2,
                                 min.expression-expression.step,
                                 min.expression-expression.step),lwd=2,border="blue",col="blue")
  
  polygon(x = c(48,64,64,48),y=c(min.expression-expression.step/2,
                                 min.expression-expression.step/2,
                                 min.expression-expression.step,
                                 min.expression-expression.step),lwd=2,border="blue")
  polygon(x = c(64,72,72,64),y=c(min.expression-expression.step/2,
                                 min.expression-expression.step/2,
                                 min.expression-expression.step,
                                 min.expression-expression.step),lwd=2,border="blue",col="blue")
  
  if(free.running == "LL")
  {
   polygon(x = c(72,88,88,72),y=c(min.expression-expression.step/2,
                                  min.expression-expression.step/2,
                                  min.expression-expression.step,
                                  min.expression-expression.step),lwd=2,border="blue")
   polygon(x = c(88,96,96,88),y=c(min.expression-expression.step/2,
                                  min.expression-expression.step/2,
                                  min.expression-expression.step,
                                  min.expression-expression.step),lwd=2,border="blue",col="lightblue")
   
   polygon(x = c(96,112,112,96),y=c(min.expression-expression.step/2,
                                    min.expression-expression.step/2,
                                    min.expression-expression.step,
                                    min.expression-expression.step),lwd=2,border="blue")
   polygon(x = c(112,120,120,112),y=c(min.expression-expression.step/2,
                                      min.expression-expression.step/2,
                                      min.expression-expression.step,
                                      min.expression-expression.step),lwd=2,border="blue",col="lightblue")
  } else if (free.running == "DD")
  {
   polygon(x = c(72,88,88,72),y=c(min.expression-expression.step/2,
                                  min.expression-expression.step/2,
                                  min.expression-expression.step,
                                  min.expression-expression.step),lwd=2,border="blue",col="lightblue")
   polygon(x = c(88,96,96,88),y=c(min.expression-expression.step/2,
                                  min.expression-expression.step/2,
                                  min.expression-expression.step,
                                  min.expression-expression.step),lwd=2,border="blue",col="blue")
   
   polygon(x = c(96,112,112,96),y=c(min.expression-expression.step/2,
                                    min.expression-expression.step/2,
                                    min.expression-expression.step,
                                    min.expression-expression.step),lwd=2,border="blue", col="lightblue")
   polygon(x = c(112,120,120,112),y=c(min.expression-expression.step/2,
                                      min.expression-expression.step/2,
                                      min.expression-expression.step,
                                      min.expression-expression.step),lwd=2,border="blue",col="blue")
  }
 }
}



genes.scaled.profile <- scaled.S.ld.mean.gene.expression
proteins.scaled.profile <- scaled.S.ld.mean.protein.abundance
cell.cycle.scaled.profile <- scaled.S.percentage.mean.LD

genes.splines <- spline(x=seq(from=0,to=24,by=4),y=c(genes.scaled.profile,genes.scaled.profile[1]),n = 100)
proteins.splines <- spline(x=seq(from=0,to=24,by=4),y=c(proteins.scaled.profile,proteins.scaled.profile[1]),n = 100)
cell.cycle.splines <- spline(x=seq(from=0,to=24,by=4),y=c(cell.cycle.scaled.profile,cell.cycle.scaled.profile[1]),n = 100)


if(min(genes.splines$y) < 0)
{
 genes.splines$y <- genes.splines$y - min(genes.splines$y)
}

if(min(proteins.splines$y) < 0)
{
 proteins.splines$y <- proteins.splines$y - min(proteins.splines$y)
}

if(min(cell.cycle.splines$y) < 0)
{
 cell.cycle.splines$y <- cell.cycle.splines$y - min(cell.cycle.splines$y)
}


plot(x = genes.splines$x,y = genes.splines$y,type="l",ylim=c(-4,12),axes=F,ylab="",xlab="",lwd=3)
lines(x = genes.splines$x, y = -genes.splines$y,lwd=3)
polygon(x=c(genes.splines$x,rev(genes.splines$x)),y = c(genes.splines$y,rev(-genes.splines$y)),col="lightblue")

lines(x = proteins.splines$x,y = 5 + proteins.splines$y,lwd=3)
lines(x = proteins.splines$x,y = 5 - proteins.splines$y,lwd=3)
polygon(x=c(proteins.splines$x,rev(proteins.splines$x)),y = c(5 + proteins.splines$y,rev(5 - proteins.splines$y)),col="blue")

lines(x = cell.cycle.splines$x,y = 10 + cell.cycle.splines$y,lwd=3)
lines(x = cell.cycle.splines$x,y = 10 - cell.cycle.splines$y,lwd=3)
polygon(x=c(cell.cycle.splines$x,rev(cell.cycle.splines$x)),y = c(10 + cell.cycle.splines$y,rev(10 - cell.cycle.splines$y)),col="darkblue")

polygon(x = c(0,0,16,16),y=c(-4,-3,-3,-4),border="blue",lwd=2)
polygon(x = c(16,16,24,24),y=c(-4,-3,-3,-4),lwd=2,col="blue",border="blue") #510 420


y = 0
y =  5
y = 10


help("smooth.spline")
help(approx)
res <- spline(x=seq(from=0,to=24,by=4),y =  c(scaled.S.ld.mean.gene.expression,scaled.S.ld.mean.gene.expression[1]),n = 100)

res$x
res$y

plot(x=seq(from=0,to=24,by=4), c(scaled.S.ld.mean.gene.expression,scaled.S.ld.mean.gene.expression[1]),type="l") 
lines(res$x,res$y)

if(min(res$y) < 0)
{
 res$y <- res$y - min(res$y)
}

plot(x=res$x,y = res$y, type="l",ylim=c(-2,2),lwd=3)
lines(x=res$x,y = -res$y,lwd=3)




lines(x=seq(from=0,to=24,by=4),c(scaled.S.ld.mean.protein.abundance,scaled.S.ld.mean.protein.abundance[1]))

lines(x=seq(from=0,to=24,by=4),c(scaled.S.percentage.mean.LD,scaled.S.percentage.mean.LD[1]))
scaled.S.percentage.mean.LD

library(circlize)

color.legend <- colorRamp2(breaks = c(0,50,100),colors = c("black","blue","yellow"))
color.legend(0:100)

plot(x=0:100,y=rep(1,101),col=color.legend(0:100),pch=15,cex=10,axes=F,xlab="",ylab="")
    


plot(x = seq(from=0,to=24,by=4),y = c(ld.cyca.expression,ld.cyca.expression[1]),type="o",lwd=3,pch=16,cex=1.5,axes=F,xlab="",ylab="")
lines(x = seq(from=0,to=24,by=4),y = c(ld.cycb.expression,ld.cycb.expression[1]),type="o",lwd=3,pch=16,cex=1.5)
lines(x = seq(from=0,to=24,by=4),y = c(ld.cycd.expression,ld.cycd.expression[1]),type="o",lwd=3,pch=16,cex=1.5)

min(sd.cycb.expression)
max(sd.cyca.expression)
plot(x = seq(from=0,to=24,by=4),y = c(sd.cyca.expression,sd.cyca.expression[1]),type="o",lwd=3,pch=16,cex=1.5,axes=F,xlab="",ylab="",ylim=c(0,1300))
lines(x = seq(from=0,to=24,by=4),y = c(sd.cycb.expression,sd.cycb.expression[1]),type="o",lwd=3,pch=16,cex=1.5)
lines(x = seq(from=0,to=24,by=4),y = c(sd.cycd.expression,sd.cycd.expression[1]),type="o",lwd=3,pch=16,cex=1.5)

plot(G1.percentage,type="o")

plot(G2.percentage,type="o")

plot(S.percentage,type="o")

points.to.remove <- 1:6
#LD 0.008629813 SD 0.010169254       Phase difference estimate -2.145421445
#          P-value for difference in phase  0.260179413
points.to.remove <- 7:12
#LD 5.752909e-05 SD 1.100559e-02
#Phase difference estimate -2.201612e+00
# P-value for difference in phase  1.576913e-01
points.to.remove <- 13:18
#LD 0.003019613 SD 0.004623511
#Phase difference estimate -2.150013434
#P-value for difference in phase  0.263884834
points.to.remove <- 19:24

cell.cycle.time.points <- seq(from=0,by=4,length.out = 18)
circacompare.data <- data.frame(time=c(cell.cycle.time.points,
                                       cell.cycle.time.points),
                                measure=c(G2.percentage.LD[1:18],
                                          G2.percentage.SD[1:18]),
                                group=c(rep("LD",18),rep("SD",18)))

result.i <- circacompare(x = circacompare.data, 
                         col_time = "time", 
                         col_group = "group", 
                         col_outcome = "measure",
                         alpha_threshold = 1)
result.i.values <- result.i$summary$value
names(result.i.values) <- result.i$summary$parameter
time.ld.sd <- seq(from=0,to=72,by=0.01)
line.ld <- wave.form(mesor = result.i.values["LD mesor estimate"],
                              amplitude = result.i.values["LD amplitude estimate"],
                              period = result.i.values["Shared period estimate"],
                              phase = result.i.values["LD peak time hours"],
                              time = time.ld.sd)
line.sd <- wave.form(mesor = result.i.values["SD mesor estimate"],
                               amplitude = result.i.values["SD amplitude estimate"],
                               period = result.i.values["Shared period estimate"],
                               phase = result.i.values["SD peak time hours"],
                               time = time.ld.sd)
plot(x = time.ld.sd,type="l",
     y = line.ld,lwd=3,col="blue",axes=F,xlab="",ylab="",
     #ylim=c(min.expression-expression.step,max.expression),
     cex.lab=1.3,main="",cex.main=2,pch=20,cex=0.8,ylim=c(5,25))
lines(x = time.ld.sd, y = line.sd,col="red", 
      type="l",lwd=3)
points(x = seq(from=0,by=4,length.out=18),G2.percentage.LD[1:18],pch=19,col="blue",cex=0.5)
points(x = seq(from=0,by=4,length.out=18),G2.percentage.SD[1:18],pch=19,col="red",cex=0.5)
axis(side=2,lwd=3,at = seq(from=10,to=25,by=5),cex.axis=1.2)

polygon(x = c(0,8,8,0),y=c(5, 5, 6.5, 6.5),lwd=2,border="red")
polygon(x = c(8,24,24,8),y=c(5, 5, 6.5, 6.5),lwd=2,border="red",col="red")
polygon(x = c(24,32,32,24),y=c(5, 5, 6.5, 6.5),lwd=2,border="red")
polygon(x = c(32,48,48,32),y=c(5, 5, 6.5, 6.5),lwd=2,border="red",col="red")
polygon(x = c(48,56,56,48),y=c(5, 5, 6.5, 6.5),lwd=2,border="red")
polygon(x = c(56,72,72,56),y=c(5, 5, 6.5, 6.5),lwd=2,border="red",col="red")

polygon(x = c(0,16,16,0),y=c(7, 7, 8.5, 8.5),lwd=2,border="blue")
polygon(x = c(16,24,24,16),y=c(7, 7, 8.5, 8.5),lwd=2,border="blue",col="blue")
polygon(x = c(24,40,40,24),y=c(7, 7, 8.5, 8.5),lwd=2,border="blue")
polygon(x = c(40,48,48,40),y=c(7, 7, 8.5, 8.5),lwd=2,border="blue",col="blue")
polygon(x = c(48,64,64,48),y=c(7, 7, 8.5, 8.5),lwd=2,border="blue")
polygon(x = c(64,72,72,64),y=c(7, 7, 8.5, 8.5),lwd=2,border="blue",col="blue")
#590 420




ld.cyca.expression <- ld.mean.gene.expression["ostta18g01570",]
ld.cycb.expression <- ld.mean.gene.expression["ostta01g06150",]
ld.cycd.expression <- ld.mean.gene.expression["ostta18g01420",]

plot(x = seq(from=0,to=24,by=4),y = c(ld.cyca.expression,ld.cyca.expression[1]),
     type="o",lwd=3,axes=F,xlab="",ylab="", ylim=c(-50,500),
     cex.lab=1.3,main="",cex.main=2,pch=20,cex=1.5,col="gold")
lines(x = seq(from=0,to=24,by=4),y = c(ld.cycb.expression,ld.cycb.expression[1]),
      lwd=3, col="purple",type="o",pch=20,cex=1.5)
lines(x = seq(from=0,to=24,by=4),y = c(ld.cycd.expression,ld.cycd.expression[1]),
      lwd=3, col="darkorange2",type="o", pch=20, cex=1.5)

axis(side=2,lwd=3,at = seq(from=0,to=500,by=100),cex.axis=1.2,las=2)

polygon(x = c(0,16,16,0),y=c(-50, -50, -20, -20),lwd=2,border="blue")
polygon(x = c(16,24,24,16),y=c(-50, -50, -20, -20),lwd=2,border="blue",col="blue") #490 460





plot(x = seq(from=0,to=24,by=4),y = c(sd.cyca.expression,sd.cyca.expression[1]),
     type="o",lwd=3,axes=F,xlab="",ylab="", ylim=c(-150,1400),
     cex.lab=1.3,main="",cex.main=2,pch=20,cex=1.5,col="gold")
lines(x = seq(from=0,to=24,by=4),y = c(sd.cycb.expression,sd.cycb.expression[1]),
      lwd=3, col="purple",type="o",pch=20,cex=1.5)
lines(x = seq(from=0,to=24,by=4),y = c(sd.cycd.expression,sd.cycd.expression[1]),
      lwd=3, col="darkorange2",type="o", pch=20, cex=1.5)

axis(side=2,lwd=3,at = seq(from=0,to=1400,by=200),cex.axis=1.2,las=2)

polygon(x = c(0,8,8,0),y=c(-150, -150, -50, -50),lwd=2,border="red")
polygon(x = c(8,24,24,8),y=c(-150, -150, -50, -50),lwd=2,border="red",col="red") #490 460


ld.cdka.expression <- scale(ld.mean.gene.expression["ostta04g00110",],center = T,scale = T)[,1]
ld.cdka.abundance <-scale(ld.mean.protein.abundance["ostta04g00110",],center = T,scale = T)[,1]

plot(x = seq(from=0,to=24,by=4),y = c(ld.cdka.expression,ld.cdka.expression[1]),
     type="o",lwd=3,axes=F,xlab="",ylab="", ylim=c(-3,2),
     cex.lab=1.3,main="",cex.main=2,pch=20,cex=1.5,col="gold",lty=2)
lines(x = seq(from=0,to=24,by=4),y = c(ld.cdka.abundance,ld.cdka.abundance[1]),
      lwd=3, col="gold",type="o",pch=20,cex=1.5)

axis(side=2,lwd=3,at = seq(from=-2,to=2,by=1),labels=rep("",5),cex.axis=1.2,las=2)

polygon(x = c(0,16,16,0),y=c(-2.5, -2.5, -2.1, -2.1),lwd=2,border="blue")
polygon(x = c(16,24,24,16),y=c(-2.5, -2.5, -2.1, -2.1),lwd=2,border="blue",col="blue") #490 460


sd.cdka.expression <- scale(sd.mean.gene.expression["ostta04g00110",],center = T,scale = T)[,1]
sd.cdka.abundance <- scale(sd.mean.protein.abundance["ostta04g00110",],center = T, scale = T)[,1]

plot(x = seq(from=0,to=24,by=4),y = c(sd.cdka.expression,sd.cdka.expression[1]),
     type="o",lwd=3,axes=F,xlab="",ylab="", ylim=c(-3,2),
     cex.lab=1.3,main="",cex.main=2,pch=20,cex=1.5,col="gold",lty=2)
lines(x = seq(from=0,to=24,by=4),y = c(sd.cdka.abundance,sd.cdka.abundance[1]),
      lwd=3, col="gold",type="o",pch=20,cex=1.5)

axis(side=2,lwd=3,at = seq(from=-2,to=2,by=1),labels=rep("",5),cex.axis=1.2,las=2)

polygon(x = c(0,8,8,0),y=c(-2.5, -2.5, -2.1, -2.1),lwd=2,border="red")
polygon(x = c(8,24,24,8),y=c(-2.5, -2.5, -2.1, -2.1),lwd=2,border="red",col="red") #490 460

ld.cdkb.expression <- scale(ld.mean.gene.expression["ostta15g00670",],center = T,scale = T)[,1]
ld.cdkb.abundance <- scale(ld.mean.protein.abundance["ostta15g00670",],center = T,scale = T)[,1]

plot(x = seq(from=0,to=24,by=4),y = c(ld.cdkb.expression,ld.cdkb.expression[1]),
     type="o",lwd=3,axes=F,xlab="",ylab="", ylim=c(-3,2),
     cex.lab=1.3,main="",cex.main=2,pch=20,cex=1.5,col="purple",lty=2)
lines(x = seq(from=0,to=24,by=4),y = c(ld.cdkb.abundance,ld.cdkb.abundance[1]),
      lwd=3, col="purple",type="o",pch=20,cex=1.5)

axis(side=2,lwd=3,at = seq(from=-2,to=2,by=1),labels=rep("",5),cex.axis=1.2,las=2)

polygon(x = c(0,16,16,0),y=c(-2.5, -2.5, -2.1, -2.1),lwd=2,border="blue")
polygon(x = c(16,24,24,16),y=c(-2.5, -2.5, -2.1, -2.1),lwd=2,border="blue",col="blue") #490 460


sd.cdkb.expression <- scale(sd.mean.gene.expression["ostta15g00670",],center = T,scale = T)[,1]
#sd.cdkb.abundance <- sd.mean.protein.abundance["ostta15g00670",]

plot(x = seq(from=0,to=24,by=4),y = c(sd.cdkb.expression,sd.cdkb.expression[1]),
     type="o",lwd=3,axes=F,xlab="",ylab="", ylim=c(-3,2),
     cex.lab=1.3,main="",cex.main=2,pch=20,cex=1.5,col="purple",lty=2)
#lines(x = seq(from=0,to=24,by=4),y = c(sd.cdka.abundance,sd.cdka.abundance[1]),
#      lwd=3, col="gold",type="o",pch=20,cex=1.5)

axis(side=2,lwd=3,at = seq(from=-2,to=2,by=1),labels=rep("",5),cex.axis=1.2,las=2)

polygon(x = c(0,8,8,0),y=c(-2.5, -2.5, -2.1, -2.1),lwd=2,border="red")
polygon(x = c(8,24,24,8),y=c(-2.5, -2.5, -2.1, -2.1),lwd=2,border="red",col="red") #490 460


ld.cdks.expression <- ld.mean.gene.expression["ostta07g01870",]
ld.cdks.abundance <- ld.mean.protein.abundance["ostta07g01870",]

sd.cdks.expression <- sd.mean.gene.expression["ostta07g01870",]
sd.cdks.abundance <- sd.mean.protein.abundance["ostta07g01870",]




plot(x = seq(from=0,by=4,to=68), unlist(ld.carotenoids[1,]),type="o",lwd=3,col="blue",ylim=c(4,12),main=rownames(ld.carotenoids)[1],axes=F,xlab="",ylab="",cex.main=2)
lines(x = seq(from=0,by=4,to=68), unlist(sd.carotenoids[1,]),type="o",lwd=3,col="red")
axis(side = 2,lwd=3,cex.axis=1.2,las=2)
polygon(x =  c(0,0,16,16),y=c(5,4.5,4.5,5),border="blue",lwd=2)
polygon(x = c(16,16,24,24),y=c(5,4.5,4.5,5),lwd=2,col="blue",border="blue")
polygon(x =  c(24,24,40,40),y=c(5,4.5,4.5,5),border="blue",lwd=2)
polygon(x = c(40,40,48,48),y=c(5,4.5,4.5,5),lwd=2,col="blue",border="blue")
polygon(x =  c(48,48,64,64),y=c(5,4.5,4.5,5),border="blue",lwd=2)
polygon(x = c(64,64,72,72),y=c(5,4.5,4.5,5),lwd=2,col="blue",border="blue")

polygon(x = c(0,0,8,8),y=c(4.25,3.75,3.75,4.25),border="red",lwd=2)
polygon(x = c(8,8,24,24),y=c(4.25,3.75,3.75,4.25),border="red",col="red",lwd=2)
polygon(x = c(24,24,32,32),y=c(4.25,3.75,3.75,4.25),border="red",lwd=2)
polygon(x = c(32,32,48,48),y=c(4.25,3.75,3.75,4.25),border="red",col="red",lwd=2)
polygon(x = c(48,48,56,56),y=c(4.25,3.75,3.75,4.25),border="red",lwd=2)
polygon(x = c(56,56,72,72),y=c(4.25,3.75,3.75,4.25),border="red",col="red",lwd=2)

plot(x = seq(from=0,by=4,to=68), unlist(ld.carotenoids[2,]),type="o",lwd=3,col="blue",ylim=c(9,18),main=rownames(ld.carotenoids)[2],axes=F,xlab="",ylab="")
lines(x = seq(from=0,by=4,to=68), unlist(sd.carotenoids[2,]),type="o",lwd=3,col="red")
axis(side = 2,lwd=3,cex.axis=1.2,las=2)
polygon(x =  c(0,0,16,16),y=c(10,9.5,9.5,10),border="blue",lwd=2)
polygon(x = c(16,16,24,24),y=c(10,9.5,9.5,10),lwd=2,col="blue",border="blue")
polygon(x =  c(24,24,40,40),y=c(10,9.5,9.5,10),border="blue",lwd=2)
polygon(x = c(40,40,48,48),y=c(10,9.5,9.5,10),lwd=2,col="blue",border="blue")
polygon(x =  c(48,48,64,64),y=c(10,9.5,9.5,10),border="blue",lwd=2)
polygon(x = c(64,64,72,72),y=c(10,9.5,9.5,10),lwd=2,col="blue",border="blue")

polygon(x =  c(0,0,8,8),y=c(9.25,8.75,8.75,9.25),border="red",lwd=2)
polygon(x = c(8,8,24,24),y=c(9.25,8.75,8.75,9.25),lwd=2,col="red",border="red")
polygon(x =  c(24,24,32,32),y=c(9.25,8.75,8.75,9.25),border="red",lwd=2)
polygon(x = c(32,32,48,48),y=c(9.25,8.75,8.75,9.25),lwd=2,col="red",border="red")
polygon(x =  c(48,48,56,56),y=c(9.25,8.75,8.75,9.25),border="red",lwd=2)
polygon(x = c(56,56,72,72),y=c(9.25,8.75,8.75,9.25),lwd=2,col="red",border="red")



plot(x = seq(from=0,by=4,to=68), unlist(ld.carotenoids[3,]),type="o",lwd=3,col="blue",ylim=c(21,34),main=rownames(ld.carotenoids)[3],axes=F,xlab="",ylab="")
lines(x = seq(from=0,by=4,to=68), unlist(sd.carotenoids[3,]),type="o",lwd=3,col="red")
axis(side = 2,lwd=3,cex.axis=1.2,las=2)
polygon(x =  c(0,0,16,16),y=c(23.5,22.5,22.5,23.5),border="blue",lwd=2)
polygon(x = c(16,16,24,24),y=c(23.5,22.5,22.5,23.5),lwd=2,col="blue",border="blue")
polygon(x =  c(24,24,40,40),y=c(23.5,22.5,22.5,23.5),border="blue",lwd=2)
polygon(x = c(40,40,48,48),y=c(23.5,22.5,22.5,23.5),lwd=2,col="blue",border="blue")
polygon(x =  c(48,48,64,64),y=c(23.5,22.5,22.5,23.5),border="blue",lwd=2)
polygon(x = c(64,64,72,72),y=c(23.5,22.5,22.5,23.5),lwd=2,col="blue",border="blue")

polygon(x =  c(0,0,8,8),y=c(22,21,21,22),border="red",lwd=2)
polygon(x = c(8,8,24,24),y=c(22,21,21,22),lwd=2,col="red",border="red")
polygon(x =  c(24,24,32,32),y=c(22,21,21,22),border="red",lwd=2)
polygon(x = c(32,32,48,48),y=c(22,21,21,22),lwd=2,col="red",border="red")
polygon(x =  c(48,48,56,56),y=c(22,21,21,22),border="red",lwd=2)
polygon(x = c(56,56,72,72),y=c(22,21,21,22),lwd=2,col="red",border="red")


plot(x = seq(from=0,by=4,to=68), unlist(ld.carotenoids[4,]),type="o",lwd=3,col="blue",ylim=c(-2.5,22),main=rownames(ld.carotenoids)[4],axes=F,xlab="",ylab="")
lines(x = seq(from=0,by=4,to=68), unlist(sd.carotenoids[4,]),type="o",lwd=3,col="red")
axis(side = 2,lwd=3,cex.axis=1.2,las=2)
polygon(x =  c(0,0,16,16),y=c(1,-0.5,-0.5,1),border="blue",lwd=2)
polygon(x = c(16,16,24,24),y=c(1,-0.5,-0.5,1),lwd=2,col="blue",border="blue")
polygon(x =  c(24,24,40,40),y=c(1,-0.5,-0.5,1),border="blue",lwd=2)
polygon(x = c(40,40,48,48),y=c(1,-0.5,-0.5,1),lwd=2,col="blue",border="blue")
polygon(x =  c(48,48,64,64),y=c(1,-0.5,-0.5,1),border="blue",lwd=2)
polygon(x = c(64,64,72,72),y=c(1,-0.5,-0.5,1),lwd=2,col="blue",border="blue")

polygon(x =  c(0,0,8,8),y=c(-1,-2.5,-2.5,-1),border="red",lwd=2)
polygon(x = c(8,8,24,24),y=c(-1,-2.5,-2.5,-1),lwd=2,col="red",border="red")
polygon(x =  c(24,24,32,32),y=c(-1,-2.5,-2.5,-1),border="red",lwd=2)
polygon(x = c(32,32,48,48),y=c(-1,-2.5,-2.5,-1),lwd=2,col="red",border="red")
polygon(x =  c(48,48,56,56),y=c(-1,-2.5,-2.5,-1),border="red",lwd=2)
polygon(x = c(56,56,72,72),y=c(-1,-2.5,-2.5,-1),lwd=2,col="red",border="red")




plot(x = seq(from=0,by=4,to=68), unlist(ld.carotenoids[5,]),type="o",lwd=3,col="blue",ylim=c(2,8),main=rownames(ld.carotenoids)[5],axes=F,xlab="",ylab="")
lines(x = seq(from=0,by=4,to=68), unlist(sd.carotenoids[5,]),type="o",lwd=3,col="red")
axis(side = 2,lwd=3,cex.axis=1.2,las=2)
polygon(x =  c(0,0,16,16),y=c(3.5,3,3,3.5),border="blue",lwd=2)
polygon(x = c(16,16,24,24),y=c(3.5,3,3,3.5),lwd=2,col="blue",border="blue")
polygon(x =  c(24,24,40,40),y=c(3.5,3,3,3.5),border="blue",lwd=2)
polygon(x = c(40,40,48,48),y=c(3.5,3,3,3.5),lwd=2,col="blue",border="blue")
polygon(x =  c(48,48,64,64),y=c(3.5,3,3,3.5),border="blue",lwd=2)
polygon(x = c(64,64,72,72),y=c(3.5,3,3,3.5),lwd=2,col="blue",border="blue")

polygon(x =  c(0,0,8,8),y=c(2.75,2.25,2.25,2.75),border="red",lwd=2)
polygon(x = c(8,8,24,24),y=c(2.75,2.25,2.25,2.75),lwd=2,col="red",border="red")
polygon(x =  c(24,24,32,32),y=c(2.75,2.25,2.25,2.75),border="red",lwd=2)
polygon(x = c(32,32,48,48),y=c(2.75,2.25,2.25,2.75),lwd=2,col="red",border="red")
polygon(x =  c(48,48,56,56),y=c(2.75,2.25,2.25,2.75),border="red",lwd=2)
polygon(x = c(56,56,72,72),y=c(2.75,2.25,2.25,2.75),lwd=2,col="red",border="red")



plot(x = seq(from=0,by=4,to=68), unlist(ld.carotenoids[6,]),type="o",lwd=3,col="blue",ylim=c(-4,22),main=rownames(ld.carotenoids)[6],axes=F,xlab="",ylab="")
lines(x = seq(from=0,by=4,to=68), unlist(sd.carotenoids[6,]),type="o",lwd=3,col="red")
axis(side = 2,lwd=3,cex.axis=1.2,las=2)
polygon(x =  c(0,0,16,16),y=c(0,-1.5,-1.5,0),border="blue",lwd=2)
polygon(x = c(16,16,24,24),y=c(0,-1.5,-1.5,0),lwd=2,col="blue",border="blue")
polygon(x =  c(24,24,40,40),y=c(0,-1.5,-1.5,0),border="blue",lwd=2)
polygon(x = c(40,40,48,48),y=c(0,-1.5,-1.5,0),lwd=2,col="blue",border="blue")
polygon(x =  c(48,48,64,64),y=c(0,-1.5,-1.5,0),border="blue",lwd=2)
polygon(x = c(64,64,72,72),y=c(0,-1.5,-1.5,0),lwd=2,col="blue",border="blue")

polygon(x =  c(0,0,8,8),y=c(-2,-3.5,-3.5,-2),border="red",lwd=2)
polygon(x = c(8,8,24,24),y=c(-2,-3.5,-3.5,-2),lwd=2,col="red",border="red")
polygon(x =  c(24,24,32,32),y=c(-2,-3.5,-3.5,-2),border="red",lwd=2)
polygon(x = c(32,32,48,48),y=c(-2,-3.5,-3.5,-2),lwd=2,col="red",border="red")
polygon(x =  c(48,48,56,56),y=c(-2,-3.5,-3.5,-2),border="red",lwd=2)
polygon(x = c(56,56,72,72),y=c(-2,-3.5,-3.5,-2),lwd=2,col="red",border="red")


plot(x = seq(from=0,by=4,to=68), unlist(ld.carotenoids[7,]),type="o",lwd=3,col="blue",ylim=c(-4,13),main=rownames(ld.carotenoids)[7],axes=F,xlab="",ylab="")
lines(x = seq(from=0,by=4,to=68), unlist(sd.carotenoids[7,]),type="o",lwd=3,col="red")
axis(side = 2,lwd=3,cex.axis=1.2,las=2)
polygon(x =  c(0,0,16,16),y=c(0,-1.5,-1.5,0),border="blue",lwd=2)
polygon(x = c(16,16,24,24),y=c(0,-1.5,-1.5,0),lwd=2,col="blue",border="blue")
polygon(x =  c(24,24,40,40),y=c(0,-1.5,-1.5,0),border="blue",lwd=2)
polygon(x = c(40,40,48,48),y=c(0,-1.5,-1.5,0),lwd=2,col="blue",border="blue")
polygon(x =  c(48,48,64,64),y=c(0,-1.5,-1.5,0),border="blue",lwd=2)
polygon(x = c(64,64,72,72),y=c(0,-1.5,-1.5,0),lwd=2,col="blue",border="blue")
polygon(x =  c(0,0,8,8),y=c(-2,-3.5,-3.5,-2),border="red",lwd=2)
polygon(x = c(8,8,24,24),y=c(-2,-3.5,-3.5,-2),lwd=2,col="red",border="red")
polygon(x =  c(24,24,32,32),y=c(-2,-3.5,-3.5,-2),border="red",lwd=2)
polygon(x = c(32,32,48,48),y=c(-2,-3.5,-3.5,-2),lwd=2,col="red",border="red")
polygon(x =  c(48,48,56,56),y=c(-2,-3.5,-3.5,-2),border="red",lwd=2)
polygon(x = c(56,56,72,72),y=c(-2,-3.5,-3.5,-2),lwd=2,col="red",border="red")

plot(x = seq(from=0,by=4,to=68), unlist(ld.carotenoids[8,]),type="o",lwd=3,col="blue",ylim=c(-1,15),main=rownames(ld.carotenoids)[8],axes=F,xlab="",ylab="")
lines(x = seq(from=0,by=4,to=68), unlist(sd.carotenoids[8,]),type="o",lwd=3,col="red")
axis(side = 2,lwd=3,cex.axis=1.2,las=2)
polygon(x =  c(0,0,16,16),y=c(2,1,1,2),border="blue",lwd=2)
polygon(x = c(16,16,24,24),y=c(2,1,1,2),lwd=2,col="blue",border="blue")
polygon(x =  c(24,24,40,40),y=c(2,1,1,2),border="blue",lwd=2)
polygon(x = c(40,40,48,48),y=c(2,1,1,2),lwd=2,col="blue",border="blue")
polygon(x =  c(48,48,64,64),y=c(2,1,1,2),border="blue",lwd=2)
polygon(x = c(64,64,72,72),y=c(2,1,1,2),lwd=2,col="blue",border="blue")
polygon(x =  c(0,0,8,8),y=c(0.5,-0.5,-0.5,0.5),border="red",lwd=2)
polygon(x = c(8,8,24,24),y=c(0.5,-0.5,-0.5,0.5),lwd=2,col="red",border="red")
polygon(x =  c(24,24,32,32),y=c(0.5,-0.5,-0.5,0.5),border="red",lwd=2)
polygon(x = c(32,32,48,48),y=c(0.5,-0.5,-0.5,0.5),lwd=2,col="red",border="red")
polygon(x =  c(48,48,56,56),y=c(0.5,-0.5,-0.5,0.5),border="red",lwd=2)
polygon(x = c(56,56,72,72),y=c(0.5,-0.5,-0.5,0.5),lwd=2,col="red",border="red")

plot(x = seq(from=0,by=4,to=68), unlist(ld.carotenoids[9,]),type="o",lwd=3,col="blue",ylim=c(-1,3),main=rownames(ld.carotenoids)[9],axes=F,xlab="",ylab="")
lines(x = seq(from=0,by=4,to=68), unlist(sd.carotenoids[9,]),type="o",lwd=3,col="red")
axis(side = 2,lwd=3,cex.axis=1.2,las=2)
polygon(x =  c(0,0,16,16),y=c(0,-0.2,-0.2,0),border="blue",lwd=2)
polygon(x = c(16,16,24,24),y=c(0,-0.2,-0.2,0),lwd=2,col="blue",border="blue")
polygon(x =  c(24,24,40,40),y=c(0,-0.2,-0.2,0),border="blue",lwd=2)
polygon(x = c(40,40,48,48),y=c(0,-0.2,-0.2,0),lwd=2,col="blue",border="blue")
polygon(x =  c(48,48,64,64),y=c(0,-0.2,-0.2,0),border="blue",lwd=2)
polygon(x = c(64,64,72,72),y=c(0,-0.2,-0.2,0),lwd=2,col="blue",border="blue")
polygon(x =  c(0,0,8,8),y=c(-0.3,-0.5,-0.5,-0.3),border="red",lwd=2)
polygon(x = c(8,8,24,24),y=c(-0.3,-0.5,-0.5,-0.3),lwd=2,col="red",border="red")
polygon(x =  c(24,24,32,32),y=c(-0.3,-0.5,-0.5,-0.3),border="red",lwd=2)
polygon(x = c(32,32,48,48),y=c(-0.3,-0.5,-0.5,-0.3),lwd=2,col="red",border="red")
polygon(x =  c(48,48,56,56),y=c(-0.3,-0.5,-0.5,-0.3),border="red",lwd=2)
polygon(x = c(56,56,72,72),y=c(-0.3,-0.5,-0.5,-0.3),lwd=2,col="red",border="red")


plot(x = seq(from=0,by=4,to=68), unlist(ld.carotenoids[10,]),type="o",lwd=3,col="blue",ylim=c(0,12),main=rownames(ld.carotenoids)[10],axes=F,xlab="",ylab="")
lines(x = seq(from=0,by=4,to=68), unlist(sd.carotenoids[10,]),type="o",lwd=3,col="red")
axis(side = 2,lwd=3,cex.axis=1.2,las=2)
polygon(x =  c(0,0,16,16),y=c(3,2,2,3),border="blue",lwd=2)
polygon(x = c(16,16,24,24),y=c(3,2,2,3),lwd=2,col="blue",border="blue")
polygon(x =  c(24,24,40,40),y=c(3,2,2,3),border="blue",lwd=2)
polygon(x = c(40,40,48,48),y=c(3,2,2,3),lwd=2,col="blue",border="blue")
polygon(x =  c(48,48,64,64),y=c(3,2,2,3),border="blue",lwd=2)
polygon(x = c(64,64,72,72),y=c(3,2,2,3),lwd=2,col="blue",border="blue")
polygon(x =  c(0,0,8,8),y=c(1.5,0.5,0.5,1.5),border="red",lwd=2)
polygon(x = c(8,8,24,24),y=c(1.5,0.5,0.5,1.5),lwd=2,col="red",border="red")
polygon(x =  c(24,24,32,32),y=c(1.5,0.5,0.5,1.5),border="red",lwd=2)
polygon(x = c(32,32,48,48),y=c(1.5,0.5,0.5,1.5),lwd=2,col="red",border="red")
polygon(x =  c(48,48,56,56),y=c(1.5,0.5,0.5,1.5),border="red",lwd=2)
polygon(x = c(56,56,72,72),y=c(1.5,0.5,0.5,1.5),lwd=2,col="red",border="red")

plot(x = seq(from=0,by=4,to=68), unlist(ld.carotenoids[11,]),type="o",lwd=3,col="blue",ylim=c(-3,8),main=rownames(ld.carotenoids)[11],axes=F,xlab="",ylab="")
lines(x = seq(from=0,by=4,to=68), unlist(sd.carotenoids[11,]),type="o",lwd=3,col="red")
axis(side = 2,lwd=3,cex.axis=1.2,las=2)
polygon(x =  c(0,0,16,16),y=c(-0.5,-1,-1,-0.5),border="blue",lwd=2)
polygon(x = c(16,16,24,24),y=c(-0.5,-1,-1,-0.5),lwd=2,col="blue",border="blue")
polygon(x =  c(24,24,40,40),y=c(-0.5,-1,-1,-0.5),border="blue",lwd=2)
polygon(x = c(40,40,48,48),y=c(-0.5,-1,-1,-0.5),lwd=2,col="blue",border="blue")
polygon(x =  c(48,48,64,64),y=c(-0.5,-1,-1,-0.5),border="blue",lwd=2)
polygon(x = c(64,64,72,72),y=c(-0.5,-1,-1,-0.5),lwd=2,col="blue",border="blue")
polygon(x =  c(0,0,8,8),y=c(-1.2,-1.7,-1.7,-1.2),border="red",lwd=2)
polygon(x = c(8,8,24,24),y=c(-1.2,-1.7,-1.7,-1.2),lwd=2,col="red",border="red")
polygon(x =  c(24,24,32,32),y=c(-1.2,-1.7,-1.7,-1.2),border="red",lwd=2)
polygon(x = c(32,32,48,48),y=c(-1.2,-1.7,-1.7,-1.2),lwd=2,col="red",border="red")
polygon(x =  c(48,48,56,56),y=c(-1.2,-1.7,-1.7,-1.2),border="red",lwd=2)
polygon(x = c(56,56,72,72),y=c(-1.2,-1.7,-1.7,-1.2),lwd=2,col="red",border="red")






gene <- "ostta16g00660"
ld.current.gene.expression <- scale(ld.mean.gene.expression[gene,],center = T,scale = T)[,1]
ld.current.protein.abundance <- scale(ld.mean.protein.abundance[gene,],center = T,scale = T)[,1]


plot(x = seq(from=0,to=24,by=4),y = c(ld.current.gene.expression,ld.current.gene.expression[1]),
     type="o",lwd=3,axes=F,xlab="",ylab="", ylim=c(-3,2),
     cex.lab=1.3,main="",cex.main=2,pch=20,cex=1.5,col="lightblue",lty=2)

arrows(x0 = seq(from=0,to=24,by=4),y0 =  c(ld.current.gene.expression,ld.current.gene.expression[1]) - 
        c(scaled.sem.ld.gene.expression[gene,],scaled.sem.ld.gene.expression[gene,1]),
       x1 = seq(from=0,to=24,by=4), y1 = c(ld.current.gene.expression,ld.current.gene.expression[1]) + 
        c(scaled.sem.ld.gene.expression[gene,],scaled.sem.ld.gene.expression[gene,1]),
       length = 0.05,angle = 90,code = 3,col = "lightblue",lwd=3)

lines(x = seq(from=0,to=24,by=4),y = c(ld.current.protein.abundance,ld.current.protein.abundance[1]),
      lwd=3, col="blue",type="o",pch=20,cex=1.5)
arrows(x0 = seq(from=0,to=24,by=4),y0 =  c(ld.current.protein.abundance,ld.current.protein.abundance[1]) - 
        c(scaled.sem.ld.protein.abundance[gene,],scaled.sem.ld.protein.abundance[gene,1]),
       x1 = seq(from=0,to=24,by=4), y1 = c(ld.current.protein.abundance,ld.current.protein.abundance[1]) + 
        c(scaled.sem.ld.protein.abundance[gene,],scaled.sem.ld.protein.abundance[gene,1]),length = 0.05,angle = 90,code = 3,col = "blue",lwd=3)
axis(side=2,lwd=3,at = seq(from=-2,to=2,by=1),labels=rep("",5),cex.axis=1.2,las=2)

polygon(x = c(0,16,16,0),y=c(-3, -3, -2.6, -2.6),lwd=2,border="blue")
polygon(x = c(16,24,24,16),y=c(-3, -3, -2.6, -2.6),lwd=2,border="blue",col="blue") #490 460



gene <- "ostta16g00670"
ld.current.gene.expression <- scale(ld.mean.gene.expression[gene,],center = T,scale = T)[,1]
ld.current.protein.abundance <- scale(ld.mean.protein.abundance[gene,],center = T,scale = T)[,1]


plot(x = seq(from=0,to=24,by=4),y = c(ld.current.gene.expression,ld.current.gene.expression[1]),
     type="o",lwd=3,axes=F,xlab="",ylab="", ylim=c(-3,2),
     cex.lab=1.3,main="",cex.main=2,pch=20,cex=1.5,col="lightblue",lty=2)

arrows(x0 = seq(from=0,to=24,by=4),y0 =  c(ld.current.gene.expression,ld.current.gene.expression[1]) - 
        c(scaled.sem.ld.gene.expression[gene,],scaled.sem.ld.gene.expression[gene,1]),
       x1 = seq(from=0,to=24,by=4), y1 = c(ld.current.gene.expression,ld.current.gene.expression[1]) + 
        c(scaled.sem.ld.gene.expression[gene,],scaled.sem.ld.gene.expression[gene,1]),
       length = 0.05,angle = 90,code = 3,col = "lightblue",lwd=3)

lines(x = seq(from=0,to=24,by=4),y = c(ld.current.protein.abundance,ld.current.protein.abundance[1]),
      lwd=3, col="blue",type="o",pch=20,cex=1.5)
arrows(x0 = seq(from=0,to=24,by=4),y0 =  c(ld.current.protein.abundance,ld.current.protein.abundance[1]) - 
        c(scaled.sem.ld.protein.abundance[gene,],scaled.sem.ld.protein.abundance[gene,1]),
       x1 = seq(from=0,to=24,by=4), y1 = c(ld.current.protein.abundance,ld.current.protein.abundance[1]) + 
        c(scaled.sem.ld.protein.abundance[gene,],scaled.sem.ld.protein.abundance[gene,1]),length = 0.05,angle = 90,code = 3,col = "blue",lwd=3)
axis(side=2,lwd=3,at = seq(from=-2,to=2,by=1),labels=rep("",5),cex.axis=1.2,las=2)

polygon(x = c(0,16,16,0),y=c(-3, -3, -2.6, -2.6),lwd=2,border="blue")
polygon(x = c(16,24,24,16),y=c(-3, -3, -2.6, -2.6),lwd=2,border="blue",col="blue") #490 460




gene <- "ostta16g00660"
sd.current.gene.expression <- scale(sd.mean.gene.expression[gene,],center = T,scale = T)[,1]
sd.current.protein.abundance <- scale(sd.mean.protein.abundance[gene,],center = T,scale = T)[,1]


plot(x = seq(from=0,to=24,by=4),y = c(sd.current.gene.expression,sd.current.gene.expression[1]),
     type="o",lwd=3,axes=F,xlab="",ylab="", ylim=c(-3,2),
     cex.lab=1.3,main="",cex.main=2,pch=20,cex=1.5,col="salmon",lty=2)

arrows(x0 = seq(from=0,to=24,by=4),y0 =  c(sd.current.gene.expression,sd.current.gene.expression[1]) - 
        c(scaled.sem.sd.gene.expression[gene,],scaled.sem.sd.gene.expression[gene,1]),
       x1 = seq(from=0,to=24,by=4), y1 = c(sd.current.gene.expression,sd.current.gene.expression[1]) + 
        c(scaled.sem.sd.gene.expression[gene,],scaled.sem.sd.gene.expression[gene,1]),
       length = 0.05,angle = 90,code = 3,col = "salmon",lwd=3)

lines(x = seq(from=0,to=24,by=4),y = c(sd.current.protein.abundance,sd.current.protein.abundance[1]),
      lwd=3, col="red",type="o",pch=20,cex=1.5)
arrows(x0 = seq(from=0,to=24,by=4),y0 =  c(sd.current.protein.abundance,sd.current.protein.abundance[1]) - 
        c(scaled.sem.sd.protein.abundance[gene,],scaled.sem.sd.protein.abundance[gene,1]),
       x1 = seq(from=0,to=24,by=4), y1 = c(sd.current.protein.abundance,sd.current.protein.abundance[1]) + 
        c(scaled.sem.sd.protein.abundance[gene,],scaled.sem.sd.protein.abundance[gene,1]),length = 0.05,angle = 90,code = 3,col = "red",lwd=3)
axis(side=2,lwd=3,at = seq(from=-2,to=2,by=1),labels=rep("",5),cex.axis=1.2,las=2)

polygon(x = c(0,8,8,0),y=c(-3, -3, -2.6, -2.6),lwd=2,border="red")
polygon(x = c(8,24,24,8),y=c(-3, -3, -2.6, -2.6),lwd=2,border="red",col="red") #490 460



gene <- "ostta16g00670"
sd.current.gene.expression <- scale(sd.mean.gene.expression[gene,],center = T,scale = T)[,1]
sd.current.protein.abundance <- scale(sd.mean.protein.abundance[gene,],center = T,scale = T)[,1]


plot(x = seq(from=0,to=24,by=4),y = c(sd.current.gene.expression,sd.current.gene.expression[1]),
     type="o",lwd=3,axes=F,xlab="",ylab="", ylim=c(-3,2),
     cex.lab=1.3,main="",cex.main=2,pch=20,cex=1.5,col="salmon",lty=2)

arrows(x0 = seq(from=0,to=24,by=4),y0 =  c(sd.current.gene.expression,sd.current.gene.expression[1]) - 
        c(scaled.sem.sd.gene.expression[gene,],scaled.sem.sd.gene.expression[gene,1]),
       x1 = seq(from=0,to=24,by=4), y1 = c(sd.current.gene.expression,sd.current.gene.expression[1]) + 
        c(scaled.sem.sd.gene.expression[gene,],scaled.sem.sd.gene.expression[gene,1]),
       length = 0.05,angle = 90,code = 3,col = "salmon",lwd=3)

lines(x = seq(from=0,to=24,by=4),y = c(sd.current.protein.abundance,sd.current.protein.abundance[1]),
      lwd=3, col="red",type="o",pch=20,cex=1.5)
arrows(x0 = seq(from=0,to=24,by=4),y0 =  c(sd.current.protein.abundance,sd.current.protein.abundance[1]) - 
        c(scaled.sem.sd.protein.abundance[gene,],scaled.sem.sd.protein.abundance[gene,1]),
       x1 = seq(from=0,to=24,by=4), y1 = c(sd.current.protein.abundance,sd.current.protein.abundance[1]) + 
        c(scaled.sem.sd.protein.abundance[gene,],scaled.sem.sd.protein.abundance[gene,1]),length = 0.05,angle = 90,code = 3,col = "red",lwd=3)
axis(side=2,lwd=3,at = seq(from=-2,to=2,by=1),labels=rep("",5),cex.axis=1.2,las=2)

polygon(x = c(0,8,8,0),y=c(-3, -3, -2.6, -2.6),lwd=2,border="red")
polygon(x = c(8,24,24,8),y=c(-3, -3, -2.6, -2.6),lwd=2,border="red",col="red") #490 460



ld.mean.zeaxanthin <- (ld.carotenoids["zeaxanthin",paste(paste("ZT",seq(from=0,to=20,by=4),sep=""),"1",sep="_")] +
                       ld.carotenoids["zeaxanthin",paste(paste("ZT",seq(from=0,to=20,by=4),sep=""),"2",sep="_")] +
                       ld.carotenoids["zeaxanthin",paste(paste("ZT",seq(from=0,to=20,by=4),sep=""),"3",sep="_")]) / 3
ld.mean.zeaxanthin.splines <- spline(x=seq(from=0,to=24,by=4),y=c(ld.mean.zeaxanthin, ld.mean.zeaxanthin[1]),n = 100)
col_fun_zeaxanthin = colorRamp2(c(min(ld.mean.zeaxanthin.splines$y),max(ld.mean.zeaxanthin.splines$y)), 
                        c("black", "orangered"))
cols <- col_fun_zeaxanthin(x = ld.mean.zeaxanthin.splines$y)
plot(x=1,y=1,col="white",xlim=c(-1,25),ylim=c(-10,10),axes=F,xlab="",ylab="")
for(i in 1:99)
{
 polygon(x = c(ld.mean.zeaxanthin.splines$x[i],ld.mean.zeaxanthin.splines$x[i+1],
               ld.mean.zeaxanthin.splines$x[i+1],ld.mean.zeaxanthin.splines$x[i]),
         y=c(-2.5, -2.5, 2.5, 2.5),lwd=2,border=cols[i],col=cols[i]) 
}


ld.mean.violaxanthin <- (ld.carotenoids["violaxanthin",paste(paste("ZT",seq(from=0,to=20,by=4),sep=""),"1",sep="_")] +
                        ld.carotenoids["violaxanthin",paste(paste("ZT",seq(from=0,to=20,by=4),sep=""),"2",sep="_")] +
                        ld.carotenoids["violaxanthin",paste(paste("ZT",seq(from=0,to=20,by=4),sep=""),"3",sep="_")]) / 3
ld.mean.violaxanthin.splines <- spline(x=seq(from=0,to=24,by=4),y=c(ld.mean.violaxanthin, ld.mean.violaxanthin[1]),n = 100)
col_fun_violaxanthin = colorRamp2(c(min(ld.mean.violaxanthin.splines$y),max(ld.mean.violaxanthin.splines$y)), 
                                c("black", "orange"))
cols <- col_fun_violaxanthin(x = ld.mean.violaxanthin.splines$y)
plot(x=1,y=1,col="white",xlim=c(-1,25),ylim=c(-10,10),axes=F,xlab="",ylab="")
for(i in 1:99)
{
 polygon(x = c(ld.mean.violaxanthin.splines$x[i],ld.mean.violaxanthin.splines$x[i+1],
               ld.mean.violaxanthin.splines$x[i+1],ld.mean.violaxanthin.splines$x[i]),
         y=c(-2.5, -2.5, 2.5, 2.5),lwd=2,border=cols[i],col=cols[i]) 
}


ld.mean.antheraxanthin <- (ld.carotenoids["antheraxanthin",paste(paste("ZT",seq(from=0,to=20,by=4),sep=""),"1",sep="_")] +
                          ld.carotenoids["antheraxanthin",paste(paste("ZT",seq(from=0,to=20,by=4),sep=""),"2",sep="_")] +
                          ld.carotenoids["antheraxanthin",paste(paste("ZT",seq(from=0,to=20,by=4),sep=""),"3",sep="_")]) / 3
ld.mean.antheraxanthin.splines <- spline(x=seq(from=0,to=24,by=4),y=c(ld.mean.antheraxanthin, 
                                                                      ld.mean.antheraxanthin[1]),n = 100)
col_fun_antheraxanthin = colorRamp2(c(min(ld.mean.antheraxanthin.splines$y),
                                      max(ld.mean.antheraxanthin.splines$y)), 
                                  c("black", "yellow"))
cols <- col_fun_antheraxanthin(x = ld.mean.antheraxanthin.splines$y)
plot(x=1,y=1,col="white",xlim=c(-1,25),ylim=c(-10,10),axes=F,xlab="",ylab="")
for(i in 1:99)
{
 polygon(x = c(ld.mean.antheraxanthin.splines$x[i],ld.mean.antheraxanthin.splines$x[i+1],
               ld.mean.antheraxanthin.splines$x[i+1],ld.mean.antheraxanthin.splines$x[i]),
         y=c(-2.5, -2.5, 2.5, 2.5),lwd=2,border=cols[i],col=cols[i]) 
}


#------------------
sd.mean.zeaxanthin <- (sd.carotenoids["zeaxanthin",paste(paste("ZT",seq(from=0,to=20,by=4),sep=""),"1",sep="_")] +
                        sd.carotenoids["zeaxanthin",paste(paste("ZT",seq(from=0,to=20,by=4),sep=""),"2",sep="_")] +
                        sd.carotenoids["zeaxanthin",paste(paste("ZT",seq(from=0,to=20,by=4),sep=""),"3",sep="_")]) / 3
sd.mean.zeaxanthin.splines <- spline(x=seq(from=0,to=24,by=4),y=c(sd.mean.zeaxanthin, sd.mean.zeaxanthin[1]),n = 100)
col_fun_zeaxanthin = colorRamp2(c(min(sd.mean.zeaxanthin.splines$y),max(sd.mean.zeaxanthin.splines$y)), 
                                c("black", "orangered"))
cols <- col_fun_zeaxanthin(x = sd.mean.zeaxanthin.splines$y)
plot(x=1,y=1,col="white",xlim=c(-1,25),ylim=c(-10,10),axes=F,xlab="",ylab="")
for(i in 1:99)
{
 polygon(x = c(sd.mean.zeaxanthin.splines$x[i],sd.mean.zeaxanthin.splines$x[i+1],
               sd.mean.zeaxanthin.splines$x[i+1],sd.mean.zeaxanthin.splines$x[i]),
         y=c(-2.5, -2.5, 2.5, 2.5),lwd=2,border=cols[i],col=cols[i]) 
}


sd.mean.violaxanthin <- (sd.carotenoids["violaxanthin",paste(paste("ZT",seq(from=0,to=20,by=4),sep=""),"1",sep="_")] +
                          sd.carotenoids["violaxanthin",paste(paste("ZT",seq(from=0,to=20,by=4),sep=""),"2",sep="_")] +
                          sd.carotenoids["violaxanthin",paste(paste("ZT",seq(from=0,to=20,by=4),sep=""),"3",sep="_")]) / 3
sd.mean.violaxanthin.splines <- spline(x=seq(from=0,to=24,by=4),y=c(sd.mean.violaxanthin, sd.mean.violaxanthin[1]),n = 100)
col_fun_violaxanthin = colorRamp2(c(min(sd.mean.violaxanthin.splines$y),max(sd.mean.violaxanthin.splines$y)), 
                                  c("black", "orange"))
cols <- col_fun_violaxanthin(x = sd.mean.violaxanthin.splines$y)
plot(x=1,y=1,col="white",xlim=c(-1,25),ylim=c(-10,10),axes=F,xlab="",ylab="")
for(i in 1:99)
{
 polygon(x = c(sd.mean.violaxanthin.splines$x[i],sd.mean.violaxanthin.splines$x[i+1],
               sd.mean.violaxanthin.splines$x[i+1],sd.mean.violaxanthin.splines$x[i]),
         y=c(-2.5, -2.5, 2.5, 2.5),lwd=2,border=cols[i],col=cols[i]) 
}


sd.mean.antheraxanthin <- (sd.carotenoids["antheraxanthin",paste(paste("ZT",seq(from=0,to=20,by=4),sep=""),"1",sep="_")] +
                            sd.carotenoids["antheraxanthin",paste(paste("ZT",seq(from=0,to=20,by=4),sep=""),"2",sep="_")] +
                            sd.carotenoids["antheraxanthin",paste(paste("ZT",seq(from=0,to=20,by=4),sep=""),"3",sep="_")]) / 3
sd.mean.antheraxanthin.splines <- spline(x=seq(from=0,to=24,by=4),y=c(sd.mean.antheraxanthin, 
                                                                      sd.mean.antheraxanthin[1]),n = 100)
col_fun_antheraxanthin = colorRamp2(c(min(sd.mean.antheraxanthin.splines$y),
                                      max(sd.mean.antheraxanthin.splines$y)), 
                                    c("black", "yellow"))
cols <- col_fun_antheraxanthin(x = sd.mean.antheraxanthin.splines$y)
plot(x=1,y=1,col="white",xlim=c(-1,25),ylim=c(-10,10),axes=F,xlab="",ylab="")
for(i in 1:99)
{
 polygon(x = c(sd.mean.antheraxanthin.splines$x[i],sd.mean.antheraxanthin.splines$x[i+1],
               sd.mean.antheraxanthin.splines$x[i+1],sd.mean.antheraxanthin.splines$x[i]),
         y=c(-2.5, -2.5, 2.5, 2.5),lwd=2,border=cols[i],col=cols[i]) 
}



gene.id <- "ostta16g00670"

heatmap.gene.protein <- function(gene.id)
{
 current.gene.mean.sd.gene.expression <- (scaled.sd.gene.expression[gene.id,seq(from=1,to=16,by=3)] + scaled.sd.gene.expression[gene.id,seq(from=2,to=17,by=3)] + scaled.sd.gene.expression[gene.id,seq(from=3,to=18,by=3)])/3
 current.gene.mean.sd.gene.expression <- c(current.gene.mean.sd.gene.expression, current.gene.mean.sd.gene.expression[1])
 spline.current.gene.mean.sd.gene.expression <- spline(x=seq(from=0,to=24,by=4),y=current.gene.mean.sd.gene.expression,n = 100)
 col_fun_gene = colorRamp2(c(min(spline.current.gene.mean.sd.gene.expression$y),
                             (min(spline.current.gene.mean.sd.gene.expression$y) + max(spline.current.gene.mean.sd.gene.expression$y)) /2,
                             max(spline.current.gene.mean.sd.gene.expression$y)), 
                           c("black","blue","yellow"))
 cols.gene <- col_fun_gene(x = spline.current.gene.mean.sd.gene.expression$y)
 
 par(mar=c(0,0,0,0))
 plot(x=1,y=1,col="white",xlim=c(-1,25),ylim=c(0,17),axes=F,xlab="",ylab="")
 
 polygon(x = c(0,8,8,0),y=c(2, 2, 0, 0),lwd=2,border="red")
 polygon(x = c(8,24,24,8),y=c(2, 2, 0, 0),lwd=2,border="red",col="red") #490 460
 
 for(i in 1:99)
 {
  polygon(x = c(spline.current.gene.mean.sd.gene.expression$x[i],spline.current.gene.mean.sd.gene.expression$x[i+1],
                spline.current.gene.mean.sd.gene.expression$x[i+1],spline.current.gene.mean.sd.gene.expression$x[i]),
          y=c(2.5, 2.5, 5, 5),lwd=2,border=cols.gene[i],col=cols.gene[i]) 
 }
 
 
 if(gene.id %in% rownames(scaled.sd.protein.abundance))
 {
  current.protein.mean.sd.abundance <- (scaled.sd.protein.abundance[gene.id,seq(from=1,to=16,by=3)] + scaled.sd.protein.abundance[gene.id,seq(from=2,to=17,by=3)] + scaled.sd.protein.abundance[gene.id,seq(from=3,to=18,by=3)])/3
  current.protein.mean.sd.abundance <- c(current.protein.mean.sd.abundance, current.protein.mean.sd.abundance[1])
  spline.current.protein.mean.sd.abundance <- spline(x=seq(from=0,to=24,by=4),y=current.protein.mean.sd.abundance,n = 100)
  col_fun_protein = colorRamp2(c(min(spline.current.protein.mean.sd.abundance$y),
                                 (min(spline.current.protein.mean.sd.abundance$y) + max(spline.current.protein.mean.sd.abundance$y)) /2,
                                 max(spline.current.protein.mean.sd.abundance$y)), 
                               c("black","blue","yellow"))
  cols.protein <- col_fun_protein(x = spline.current.protein.mean.sd.abundance$y)
  
  for(i in 1:99)
  {
   polygon(x = c(spline.current.gene.mean.sd.gene.expression$x[i],spline.current.gene.mean.sd.gene.expression$x[i+1],
                 spline.current.gene.mean.sd.gene.expression$x[i+1],spline.current.gene.mean.sd.gene.expression$x[i]),
           y=c(5.25, 5.25, 7.75, 7.75),lwd=2,border=cols.protein[i],col=cols.protein[i]) 
  }
 }
 
 current.gene.mean.ld.gene.expression <- (scaled.ld.gene.expression[gene.id,seq(from=1,to=16,by=3)] + scaled.ld.gene.expression[gene.id,seq(from=2,to=17,by=3)] + scaled.ld.gene.expression[gene.id,seq(from=3,to=18,by=3)])/3
 current.gene.mean.ld.gene.expression <- c(current.gene.mean.ld.gene.expression, current.gene.mean.ld.gene.expression[1])
 spline.current.gene.mean.ld.gene.expression <- spline(x=seq(from=0,to=24,by=4),y=current.gene.mean.ld.gene.expression,n = 100)
 col_fun_gene = colorRamp2(c(min(spline.current.gene.mean.ld.gene.expression$y),
                             (min(spline.current.gene.mean.ld.gene.expression$y) + max(spline.current.gene.mean.ld.gene.expression$y)) /2,
                             max(spline.current.gene.mean.ld.gene.expression$y)), 
                           c("black","blue","yellow"))
 cols.gene <- col_fun_gene(x = spline.current.gene.mean.ld.gene.expression$y)
 
 
 polygon(x = c(0,16,16,0),y=c(10.5, 10.5, 8.5, 8.5),lwd=2,border="blue")
 polygon(x = c(16,24,24,16),y=c(10.5, 10.5, 8.5, 8.5),lwd=2,border="blue",col="blue") #490 460
 
 for(i in 1:99)
 {
  polygon(x = c(spline.current.gene.mean.ld.gene.expression$x[i],spline.current.gene.mean.ld.gene.expression$x[i+1],
                spline.current.gene.mean.ld.gene.expression$x[i+1],spline.current.gene.mean.ld.gene.expression$x[i]),
          y=c(11, 11, 13.5, 13.5),lwd=2,border=cols.gene[i],col=cols.gene[i]) 
 }
 
 
 if(gene.id %in% rownames(scaled.ld.protein.abundance))
 {
  current.protein.mean.ld.abundance <- (scaled.ld.protein.abundance[gene.id,seq(from=1,to=16,by=3)] + scaled.ld.protein.abundance[gene.id,seq(from=2,to=17,by=3)] + scaled.ld.protein.abundance[gene.id,seq(from=3,to=18,by=3)])/3
  current.protein.mean.ld.abundance <- c(current.protein.mean.ld.abundance, current.protein.mean.ld.abundance[1])
  spline.current.protein.mean.ld.abundance <- spline(x=seq(from=0,to=24,by=4),y=current.protein.mean.ld.abundance,n = 100)
  col_fun_protein = colorRamp2(c(min(spline.current.protein.mean.ld.abundance$y),
                                 (min(spline.current.protein.mean.ld.abundance$y) + max(spline.current.protein.mean.ld.abundance$y)) /2,
                                 max(spline.current.protein.mean.ld.abundance$y)), 
                               c("black","blue","yellow"))
  cols.protein <- col_fun_protein(x = spline.current.protein.mean.ld.abundance$y)
  
  for(i in 1:99)
  {
   polygon(x = c(spline.current.gene.mean.ld.gene.expression$x[i],spline.current.gene.mean.ld.gene.expression$x[i+1],
                 spline.current.gene.mean.ld.gene.expression$x[i+1],spline.current.gene.mean.ld.gene.expression$x[i]),
           y=c(13.75, 13.75, 16.25, 16.25),lwd=2,border=cols.protein[i],col=cols.protein[i]) 
  }
 }
}

png(filename = "DXS_ostta02g01730.png", width = 410, height = 280)
heatmap.gene.protein(gene.id = "ostta02g01730")
dev.off()

png(filename = "DXR_ostta04g03270.png", width = 410, height = 280)
heatmap.gene.protein(gene.id = "ostta04g03270")
dev.off()

png(filename = "CMS_ostta07g03570.png", width = 410, height = 280)
heatmap.gene.protein(gene.id = "ostta07g03570")
dev.off()


png(filename = "CMK_ostta16g01910.png", width = 410, height = 280)
heatmap.gene.protein(gene.id = "ostta16g01910")
dev.off()

png(filename = "MCS_ostta11g00050.png", width = 410, height = 280)
heatmap.gene.protein(gene.id = "ostta11g00050")
dev.off()


png(filename = "HDS_ostta09g00960.png", width = 410, height = 280)
heatmap.gene.protein(gene.id = "ostta09g00960")
dev.off()

png(filename = "HDR_ostta08g01180.png", width = 410, height = 280)
heatmap.gene.protein(gene.id = "ostta08g01180")
dev.off()

png(filename = "GPPS_ostta04g03170.png", width = 410, height = 280)
heatmap.gene.protein(gene.id = "ostta04g03170")
dev.off()

png(filename = "GGPPS_ostta12g01280.png", width = 410, height = 280)
heatmap.gene.protein(gene.id = "ostta12g01280")
dev.off()

png(filename = "PSY_ostta05g03530.png", width = 410, height = 280)
heatmap.gene.protein(gene.id = "ostta05g03530")
dev.off()

png(filename = "PDS_ostta16g00590.png", width = 410, height = 280)
heatmap.gene.protein(gene.id = "ostta16g00590")
dev.off()

png(filename = "PDS_ostta10g02860.png", width = 410, height = 280)
heatmap.gene.protein(gene.id = "ostta10g02860")
dev.off()

png(filename = "ZDS_ostta16g00590.png", width = 410, height = 280)
heatmap.gene.protein(gene.id = "ostta16g00590")
dev.off()

png(filename = "LCYbe_ostta14g00700.png", width = 410, height = 280)
heatmap.gene.protein(gene.id = "ostta14g00700")
dev.off()

png(filename = "CYP97_ostta01g05550.png", width = 410, height = 280)
heatmap.gene.protein(gene.id = "ostta01g05550")
dev.off()


png(filename = "CYP97_ostta15g00680.png", width = 410, height = 280)
heatmap.gene.protein(gene.id = "ostta15g00680")
dev.off()


png(filename = "CHYb_ostta03g05610.png", width = 410, height = 280)
heatmap.gene.protein(gene.id = "ostta03g05610")
dev.off()

png(filename = "CHYb_ostta18g01970.png", width = 410, height = 280)
heatmap.gene.protein(gene.id = "ostta18g01970")
dev.off()

png(filename = "ECH_ostta09g02540.png", width = 410, height = 280)
heatmap.gene.protein(gene.id = "ostta09g02540")
dev.off()


png(filename = "VDE_ostta16g00660.png", width = 410, height = 280)
heatmap.gene.protein(gene.id = "ostta16g00660")
dev.off()

png(filename = "ZEP_ostta16g00670.png", width = 410, height = 280)
heatmap.gene.protein(gene.id = "ostta16g00670")
dev.off()

png(filename = "ZEP_ostta02g02500.png", width = 410, height = 280)
heatmap.gene.protein(gene.id = "ostta02g02500")
dev.off()


png(filename = "NXS_ostta18g01970.png", width = 410, height = 280)
heatmap.gene.protein(gene.id = "ostta18g01970")
dev.off()

png(filename = "NXS_ostta03g05610.png", width = 410, height = 280)
heatmap.gene.protein(gene.id = "ostta03g05610")
dev.off()

png(filename = "CHYb_ostta13g02440.png", width = 410, height = 280)
heatmap.gene.protein(gene.id = "ostta13g02440")
dev.off()

head(ld.mean.gene.expression)

ld.mean.gene.correlation <- sort(cor(t(ld.mean.gene.expression))["ostta15g02860",],decreasing = T)
write.table(x = data.frame(names(ld.mean.gene.correlation),ld.mean.gene.correlation),file = "LD_gene_correlation.tsv",quote = F,sep = "\t",row.names = F,col.names = F)

sd.mean.gene.correlation <- sort(cor(t(sd.mean.gene.expression))["ostta15g02860",],decreasing = T)
write.table(x = data.frame(names(sd.mean.gene.correlation),sd.mean.gene.correlation),file = "SD_gene_correlation.tsv",quote = F,sep = "\t",row.names = F,col.names = F)

ld.mean.protein.correlation <- sort(cor(t(ld.mean.protein.abundance))["ostta15g02860",],decreasing = T)
write.table(x = data.frame(names(ld.mean.protein.correlation),ld.mean.protein.correlation),file = "LD_protein_correlation.tsv",quote = F,sep = "\t",row.names = F,col.names = F)

sd.mean.protein.correlation <- sort(cor(t(sd.mean.protein.abundance))["ostta15g02860",],decreasing = T)
write.table(x = data.frame(names(sd.mean.protein.correlation),sd.mean.protein.correlation),file = "SD_protein_correlation.tsv",quote = F,sep = "\t",row.names = F,col.names = F)



carotenoid <- "neoxanthin"
carotenoid.abundance <- sd.carotenoids

heatmap.carotenoid <- function(carotenoid, carotenoid.abundance)
{
 mean.carotenoid <- (carotenoid.abundance[carotenoid,paste(paste("ZT",seq(from=0,to=20,by=4),sep=""),"1",sep="_")] +
                      carotenoid.abundance[carotenoid,paste(paste("ZT",seq(from=0,to=20,by=4),sep=""),"2",sep="_")] +
                      carotenoid.abundance[carotenoid,paste(paste("ZT",seq(from=0,to=20,by=4),sep=""),"3",sep="_")]) / 3
 
 mean.carotenoid.splines <- spline(x=seq(from=0,to=24,by=4),y=c(mean.carotenoid, mean.carotenoid[1]),n = 100)
 col_fun_carotenoid = colorRamp2(c(min(mean.carotenoid.splines$y),
                                   (min(mean.carotenoid.splines$y) + max(mean.carotenoid.splines$y))/2,
                                   max(mean.carotenoid.splines$y)), 
                                 c("black", "blue", "yellow"))
 cols <- col_fun_carotenoid(x = mean.carotenoid.splines$y)
 par(mar=c(0,0,0,0))
 plot(x=1,y=1,col="white",xlim=c(-1,25),ylim=c(-2,2),axes=F,xlab="",ylab="")
 for(i in 1:99)
 {
  polygon(x = c(mean.carotenoid.splines$x[i],mean.carotenoid.splines$x[i+1],
                mean.carotenoid.splines$x[i+1],mean.carotenoid.splines$x[i]),
          y=c(-2, -2, 2, 2),lwd=2,border=cols[i],col=cols[i]) 
 }
}

rownames(ld.carotenoids)

png(filename = "LD_a_carotene.png", width = 410, height = 35)
heatmap.carotenoid(carotenoid="a-carotene", carotenoid.abundance=ld.carotenoids)
dev.off()

png(filename = "SD_a_carotene.png", width = 410, height = 35)
heatmap.carotenoid(carotenoid="a-carotene", carotenoid.abundance=sd.carotenoids)
dev.off()

png(filename = "LD_lutein.png", width = 410, height = 35)
heatmap.carotenoid(carotenoid="lutein", carotenoid.abundance=ld.carotenoids)
dev.off()

png(filename = "LD_prasinoxanthin.png", width = 410, height = 35)
heatmap.carotenoid(carotenoid="prasinoxanthin", carotenoid.abundance=ld.carotenoids)
dev.off()

png(filename = "SD_prasinoxanthin.png", width = 410, height = 35)
heatmap.carotenoid(carotenoid="prasinoxanthin", carotenoid.abundance=sd.carotenoids)
dev.off()

png(filename = "LD_dihydrolutein.png", width = 410, height = 35)
heatmap.carotenoid(carotenoid="dihydrolutein", carotenoid.abundance=ld.carotenoids)
dev.off()

png(filename = "SD_dihydrolutein.png", width = 410, height = 35)
heatmap.carotenoid(carotenoid="dihydrolutein", carotenoid.abundance=sd.carotenoids)
dev.off()

png(filename = "LD_micromonal.png", width = 410, height = 35)
heatmap.carotenoid(carotenoid="micromonal", carotenoid.abundance=ld.carotenoids)
dev.off()

png(filename = "SD_micromonal.png", width = 410, height = 35)
heatmap.carotenoid(carotenoid="micromonal", carotenoid.abundance=sd.carotenoids)
dev.off()

png(filename = "LD_uriolide.png", width = 410, height = 35)
heatmap.carotenoid(carotenoid="uriolide", carotenoid.abundance=ld.carotenoids)
dev.off()

png(filename = "SD_uriolide.png", width = 410, height = 35)
heatmap.carotenoid(carotenoid="uriolide", carotenoid.abundance=sd.carotenoids)
dev.off()

png(filename = "LD_b_carotene.png", width = 410, height = 35)
heatmap.carotenoid(carotenoid="b-carotene", carotenoid.abundance=ld.carotenoids)
dev.off()

png(filename = "SD_b_carotene.png", width = 410, height = 35)
heatmap.carotenoid(carotenoid="b-carotene", carotenoid.abundance=sd.carotenoids)
dev.off()

png(filename = "LD_zeaxanthin.png", width = 410, height = 35)
heatmap.carotenoid(carotenoid="zeaxanthin", carotenoid.abundance=ld.carotenoids)
dev.off()

png(filename = "SD_zeaxanthin.png", width = 410, height = 35)
heatmap.carotenoid(carotenoid="zeaxanthin", carotenoid.abundance=sd.carotenoids)
dev.off()

png(filename = "LD_antheraxanthin.png", width = 410, height = 35)
heatmap.carotenoid(carotenoid="antheraxanthin", carotenoid.abundance=ld.carotenoids)
dev.off()

png(filename = "SD_antheraxanthin.png", width = 410, height = 35)
heatmap.carotenoid(carotenoid="antheraxanthin", carotenoid.abundance=sd.carotenoids)
dev.off()

png(filename = "LD_violaxanthin.png", width = 410, height = 35)
heatmap.carotenoid(carotenoid="violaxanthin", carotenoid.abundance=ld.carotenoids)
dev.off()

png(filename = "SD_violaxanthin.png", width = 410, height = 35)
heatmap.carotenoid(carotenoid="violaxanthin", carotenoid.abundance=sd.carotenoids)
dev.off()

png(filename = "LD_neoxanthin.png", width = 410, height = 35)
heatmap.carotenoid(carotenoid="neoxanthin", carotenoid.abundance=ld.carotenoids)
dev.off()

png(filename = "SD_neoxanthin.png", width = 410, height = 35)
heatmap.carotenoid(carotenoid="neoxanthin", carotenoid.abundance=sd.carotenoids)
dev.off()





png(filename = "PPH1_ostta02g02330.png", width = 410, height = 280)
heatmap.gene.protein(gene.id = "ostta02g02330")
dev.off()

png(filename = "LHCSR_ostta05g01780.png", width = 410, height = 280)
heatmap.gene.protein(gene.id = "ostta05g01780")
dev.off()

png(filename = "PsbS_ostta06g04600.png", width = 410, height = 280)
heatmap.gene.protein(gene.id = "ostta06g04600")
dev.off()

png(filename = "PsbO_ostta14g00150.png", width = 410, height = 280)
heatmap.gene.protein(gene.id = "ostta14g00150")
dev.off()

png(filename = "PsbP_ostta14g02630.png", width = 410, height = 280)
heatmap.gene.protein(gene.id = "ostta14g02630")
dev.off()

png(filename = "PsbPL_ostta07g02340.png", width = 410, height = 280)
heatmap.gene.protein(gene.id = "ostta07g02340")
dev.off()

png(filename = "PsbQ1_ostta16g01620.png", width = 410, height = 280)
heatmap.gene.protein(gene.id = "ostta16g01620")
dev.off()

png(filename = "PsbR_ostta05g04560.png", width = 410, height = 280)
heatmap.gene.protein(gene.id = "ostta05g04560")
dev.off()

png(filename = "PsbW_ostta02g02320.png", width = 410, height = 280)
heatmap.gene.protein(gene.id = "ostta02g02320")
dev.off()

png(filename = "PsbX_ostta02g02560.png", width = 410, height = 280)
heatmap.gene.protein(gene.id = "ostta02g02560")
dev.off()

png(filename = "PsbY_ostta12g00520.png", width = 410, height = 280)
heatmap.gene.protein(gene.id = "ostta12g00520")
dev.off()

png(filename = "Psb27A_ostta08g03340.png", width = 410, height = 280)
heatmap.gene.protein(gene.id = "ostta08g03340")
dev.off()

png(filename = "Psb27B_ostta12g01970.png", width = 410, height = 280)
heatmap.gene.protein(gene.id = "ostta12g01970")
dev.off()

png(filename = "Psb27C_ostta04g04065.png", width = 410, height = 280)
heatmap.gene.protein(gene.id = "ostta04g04065")
dev.off()

png(filename = "Psb28_ostta16g00650.png", width = 410, height = 280)
heatmap.gene.protein(gene.id = "ostta16g00650")
dev.off()

png(filename = "LHCB4_ostta16g00650.png", width = 410, height = 280)
heatmap.gene.protein(gene.id = "ostta16g00650")
dev.off()

png(filename = "LHCB5_ostta14g00065.png", width = 410, height = 280)
heatmap.gene.protein(gene.id = "ostta14g00065")
dev.off()


png(filename = "LHC14_ostta06g04600.png", width = 410, height = 280)
heatmap.gene.protein(gene.id = "ostta06g04600")
dev.off()

png(filename = "LHC15_ostta11g00990.png", width = 410, height = 280)
heatmap.gene.protein(gene.id = "ostta11g00990")
dev.off()

png(filename = "PTOX_ostta16g00930.png", width = 410, height = 280)
heatmap.gene.protein(gene.id = "ostta16g00930")
dev.off()


png(filename = "PetC_ostta01g06610.png", width = 410, height = 280)
heatmap.gene.protein(gene.id = "ostta01g06610")
dev.off()

png(filename = "PetM_ostta04g04150.png", width = 410, height = 280)
heatmap.gene.protein(gene.id = "ostta04g04150") #no
dev.off()

png(filename = "PetN_ostta12g00550.png", width = 410, height = 280)
heatmap.gene.protein(gene.id = "ostta12g00550") #no
dev.off()

png(filename = "PetE_ostta09g04240.png", width = 410, height = 280)
heatmap.gene.protein(gene.id = "ostta09g04240")
dev.off()

png(filename = "PetF_ostta17g00310.png", width = 410, height = 280)
heatmap.gene.protein(gene.id = "ostta17g00310")
dev.off()



png(filename = "PsaD_ostta15g02670.png", width = 410, height = 280)
heatmap.gene.protein(gene.id = "ostta15g02670")
dev.off()

png(filename = "PsaE_ostta02g03860.png", width = 410, height = 280)
heatmap.gene.protein(gene.id = "ostta02g03860")
dev.off()

png(filename = "PsaF_ostta04g01790.png", width = 410, height = 280)
heatmap.gene.protein(gene.id = "ostta04g01790")
dev.off()

png(filename = "PsaG_ostta05g00800.png", width = 410, height = 280)
heatmap.gene.protein(gene.id = "ostta05g00800")
dev.off()

png(filename = "PsaH_ostta15g00990.png", width = 410, height = 280)
heatmap.gene.protein(gene.id = "ostta15g00990")
dev.off()

png(filename = "PsaI_ostta17g00570.png", width = 410, height = 280)
heatmap.gene.protein(gene.id = "ostta17g00570") #no
dev.off()

png(filename = "PsaL_ostta02g00580.png", width = 410, height = 280)
heatmap.gene.protein(gene.id = "ostta02g00580") 
dev.off()

png(filename = "PsaN_ostta06g00250.png", width = 410, height = 280)
heatmap.gene.protein(gene.id = "ostta06g00250")
dev.off()

png(filename = "PsaO_ostta02g05340.png", width = 410, height = 280)
heatmap.gene.protein(gene.id = "ostta02g05340") #no
dev.off()


png(filename = "PGR5_ostta13g00960.png", width = 410, height = 280)
heatmap.gene.protein(gene.id = "ostta13g00960")
dev.off()

png(filename = "PGRL1_ostta11g00330.png", width = 410, height = 280)
heatmap.gene.protein(gene.id = "ostta11g00330")
dev.off()

png(filename = "Fd_ostta17g00310.png", width = 410, height = 280)
heatmap.gene.protein(gene.id = "ostta17g00310")
dev.off()


png(filename = "FNR_ostta18g01250.png", width = 410, height = 280)
heatmap.gene.protein(gene.id = "ostta18g01250")
dev.off()


#ATP synthase


#nitrate assimilation
png(filename = "NRT2_ostta10g00950.png", width = 410, height = 280)
heatmap.gene.protein(gene.id = "ostta10g00950")
dev.off()


png(filename = "NRT3_ostta10g00940.png", width = 410, height = 280)
heatmap.gene.protein(gene.id = "ostta10g00940")
dev.off()


png(filename = "NR_ostta10g00920.png", width = 410, height = 280)
heatmap.gene.protein(gene.id = "ostta10g00920")
dev.off()

png(filename = "NIR_ostta10g00930.png", width = 410, height = 280)
heatmap.gene.protein(gene.id = "ostta10g00930")
dev.off()

png(filename = "GS_ostta01g05020.png", width = 410, height = 280)
heatmap.gene.protein(gene.id = "ostta01g05020")
dev.off()


png(filename = "GOGAT_ostta14g01900.png", width = 410, height = 280)
heatmap.gene.protein(gene.id = "ostta14g01900")
dev.off()



png(filename = "ATPB_ostta01g00830.png", width = 410, height = 280)
heatmap.gene.protein(gene.id = "ostta01g00830")
dev.off()


png(filename = "ATPC1_ostta09g01080.png", width = 410, height = 280)
heatmap.gene.protein(gene.id = "ostta09g01080")
dev.off()

png(filename = "ATPC2_ostta01g03770.png", width = 410, height = 280)
heatmap.gene.protein(gene.id = "ostta01g03770")
dev.off()

png(filename = "ATPD_ostta07g01350.png", width = 410, height = 280)
heatmap.gene.protein(gene.id = "ostta07g01350")
dev.off()

png(filename = "ATPG_ostta09g01080.png", width = 410, height = 280)
heatmap.gene.protein(gene.id = "ostta09g01080")
dev.off()

## Calvin cycle
png(filename = "RBCS_ostta18g01890.png", width = 410, height = 280)
heatmap.gene.protein(gene.id = "ostta18g01890")
dev.off()

png(filename = "PGK_ostta06g00700.png", width = 410, height = 280)
heatmap.gene.protein(gene.id = "ostta06g00700")
dev.off()


png(filename = "GAPDHA_ostta01g01560.png", width = 410, height = 280)
heatmap.gene.protein(gene.id = "ostta01g01560")
dev.off()


png(filename = "GAPDHB_ostta10g03420.png", width = 410, height = 280)
heatmap.gene.protein(gene.id = "ostta10g03420")
dev.off()


png(filename = "TPI_ostta09g00120.png", width = 410, height = 280)
heatmap.gene.protein(gene.id = "ostta09g00120")
dev.off()

png(filename = "FBA1_ostta01g03040.png", width = 410, height = 280)
heatmap.gene.protein(gene.id = "ostta01g03040")
dev.off()

png(filename = "FBA2_ostta03g00660.png", width = 410, height = 280)
heatmap.gene.protein(gene.id = "ostta03g00660")
dev.off()

png(filename = "FBPase_ostta03g00350.png", width = 410, height = 280)
heatmap.gene.protein(gene.id = "ostta03g00350")
dev.off()

png(filename = "FBPase_ostta14g01010.png", width = 410, height = 280)
heatmap.gene.protein(gene.id = "ostta14g01010")
dev.off()



png(filename = "TRK_ostta07g04370.png", width = 410, height = 280)
heatmap.gene.protein(gene.id = "ostta07g04370")
dev.off()


png(filename = "SBPase_ostta03g05500.png", width = 410, height = 280)
heatmap.gene.protein(gene.id = "ostta03g05500")
dev.off()

png(filename = "RPI1_ostta05g04130.png", width = 410, height = 280)
heatmap.gene.protein(gene.id = "ostta05g04130")
dev.off()

png(filename = "RPI2_ostta17g02290.png", width = 410, height = 280)
heatmap.gene.protein(gene.id = "ostta17g02290")
dev.off()


png(filename = "RPE_ostta06g01650.png", width = 410, height = 280)
heatmap.gene.protein(gene.id = "ostta06g01650")
dev.off()

png(filename = "PRK_ostta04g02740.png", width = 410, height = 280)
heatmap.gene.protein(gene.id = "ostta04g02740")
dev.off()

png(filename = "PGI_ostta11g02780.png", width = 410, height = 280)
heatmap.gene.protein(gene.id = "ostta11g02780")
dev.off()

png(filename = "PGM_ostta06g02930.png", width = 410, height = 280)
heatmap.gene.protein(gene.id = "ostta06g02930")
dev.off()

png(filename = "PGM_ostta06g02930.png", width = 410, height = 280)
heatmap.gene.protein(gene.id = "ostta06g02930")
dev.off()


png(filename = "PGM_ostta15g02520.png", width = 410, height = 280)
heatmap.gene.protein(gene.id = "ostta15g02520")
dev.off()

png(filename = "APL_ostta07g03440.png", width = 410, height = 280)
heatmap.gene.protein(gene.id = "ostta07g03440")
dev.off()

png(filename = "APS_ostta07g03070.png", width = 410, height = 280)
heatmap.gene.protein(gene.id = "ostta07g03070")
dev.off()

png(filename = "GBSS_ostta06g02940.png", width = 410, height = 280)
heatmap.gene.protein(gene.id = "ostta06g02940")
dev.off()

png(filename = "SS1_ostta13g01180.png", width = 410, height = 280)
heatmap.gene.protein(gene.id = "ostta13g01180")
dev.off()

png(filename = "SS2_ostta16g02480.png", width = 410, height = 280)
heatmap.gene.protein(gene.id = "ostta16g02480")
dev.off()

png(filename = "SS3_ostta16g01560.png", width = 410, height = 280)
heatmap.gene.protein(gene.id = "ostta16g01560")
dev.off()


png(filename = "SBE_ostta03g00870.png", width = 410, height = 280)
heatmap.gene.protein(gene.id = "ostta03g00870")
dev.off()

png(filename = "SBE_ostta04g03940.png", width = 410, height = 280)
heatmap.gene.protein(gene.id = "ostta04g03940")
dev.off()

png(filename = "ISA_ostta12g00320.png", width = 410, height = 280)
heatmap.gene.protein(gene.id = "ostta12g00320")
dev.off()

png(filename = "AMY_ostta10g00260.png", width = 410, height = 280)
heatmap.gene.protein(gene.id = "ostta10g00260")
dev.off()






gene <- "ostta10g00260"
ld.current.gene.expression <- scale(ld.mean.gene.expression[gene,],center = T,scale = T)[,1]
ld.current.protein.abundance <- scale(ld.mean.protein.abundance[gene,],center = T,scale = T)[,1]

plot(x = seq(from=0,to=24,by=4),y = c(ld.current.gene.expression,ld.current.gene.expression[1]),
     type="o",lwd=3,axes=F,xlab="",ylab="", ylim=c(-3,2),
     cex.lab=1.3,main="",cex.main=2,pch=20,cex=1.5,col="lightblue",lty=2)

arrows(x0 = seq(from=0,to=24,by=4),y0 =  c(ld.current.gene.expression,ld.current.gene.expression[1]) - 
        c(scaled.sem.ld.gene.expression[gene,],scaled.sem.ld.gene.expression[gene,1]),
       x1 = seq(from=0,to=24,by=4), y1 = c(ld.current.gene.expression,ld.current.gene.expression[1]) + 
        c(scaled.sem.ld.gene.expression[gene,],scaled.sem.ld.gene.expression[gene,1]),
       length = 0.05,angle = 90,code = 3,col = "lightblue",lwd=3)

lines(x = seq(from=0,to=24,by=4),y = c(ld.current.protein.abundance,ld.current.protein.abundance[1]),
      lwd=3, col="blue",type="o",pch=20,cex=1.5)
arrows(x0 = seq(from=0,to=24,by=4),y0 =  c(ld.current.protein.abundance,ld.current.protein.abundance[1]) - 
        c(scaled.sem.ld.protein.abundance[gene,],scaled.sem.ld.protein.abundance[gene,1]),
       x1 = seq(from=0,to=24,by=4), y1 = c(ld.current.protein.abundance,ld.current.protein.abundance[1]) + 
        c(scaled.sem.ld.protein.abundance[gene,],scaled.sem.ld.protein.abundance[gene,1]),length = 0.05,angle = 90,code = 3,col = "blue",lwd=3)
axis(side=2,lwd=3,at = seq(from=-2,to=2,by=1),labels=rep("",5),cex.axis=1.2,las=2)

polygon(x = c(0,16,16,0),y=c(-3, -3, -2.6, -2.6),lwd=2,border="blue")
polygon(x = c(16,24,24,16),y=c(-3, -3, -2.6, -2.6),lwd=2,border="blue",col="blue") #490 460


sd.current.gene.expression <- scale(sd.mean.gene.expression[gene,],center = T,scale = T)[,1]
sd.current.protein.abundance <- scale(sd.mean.protein.abundance[gene,],center = T,scale = T)[,1]

plot(x = seq(from=0,to=24,by=4),y = c(sd.current.gene.expression,sd.current.gene.expression[1]),
     type="o",lwd=3,axes=F,xlab="",ylab="", ylim=c(-3,2),
     cex.lab=1.3,main="",cex.main=2,pch=20,cex=1.5,col="salmon",lty=2)

arrows(x0 = seq(from=0,to=24,by=4),y0 =  c(sd.current.gene.expression,sd.current.gene.expression[1]) - 
        c(scaled.sem.sd.gene.expression[gene,],scaled.sem.sd.gene.expression[gene,1]),
       x1 = seq(from=0,to=24,by=4), y1 = c(sd.current.gene.expression,sd.current.gene.expression[1]) + 
        c(scaled.sem.sd.gene.expression[gene,],scaled.sem.sd.gene.expression[gene,1]),
       length = 0.05,angle = 90,code = 3,col = "salmon",lwd=3)

lines(x = seq(from=0,to=24,by=4),y = c(sd.current.protein.abundance,sd.current.protein.abundance[1]),
      lwd=3, col="red",type="o",pch=20,cex=1.5)
arrows(x0 = seq(from=0,to=24,by=4),y0 =  c(sd.current.protein.abundance,sd.current.protein.abundance[1]) - 
        c(scaled.sem.sd.protein.abundance[gene,],scaled.sem.sd.protein.abundance[gene,1]),
       x1 = seq(from=0,to=24,by=4), y1 = c(sd.current.protein.abundance,sd.current.protein.abundance[1]) + 
        c(scaled.sem.sd.protein.abundance[gene,],scaled.sem.sd.protein.abundance[gene,1]),length = 0.05,angle = 90,code = 3,col = "red",lwd=3)
axis(side=2,lwd=3,at = seq(from=-2,to=2,by=1),labels=rep("",5),cex.axis=1.2,las=2)

polygon(x = c(0,8,8,0),y=c(-3, -3, -2.6, -2.6),lwd=2,border="red")
polygon(x = c(8,24,24,8),y=c(-3, -3, -2.6, -2.6),lwd=2,border="red",col="red") #490 460





starch.data <- read.table(file="physiological_data/starch/starch.tsv")

sorted.zts <- c(paste("ZT0",1:2,sep="_"),paste("ZT4",1:2,sep="_"),paste("ZT8",1:2,sep="_"),
                paste("ZT12",1:2,sep="_"),paste("ZT16",1:2,sep="_"),paste("ZT20",1:2,sep="_"))

starch.content.ld <- starch.data$V2[starch.data$V3 == "LD"]
names(starch.content.ld) <- starch.data$V1[starch.data$V3 == "LD"]
rain.ld.starch <- rain(starch.content.ld[sorted.zts], deltat=4, period=24, verbose=FALSE, nr.series=2)

starch.content.sd <- starch.data$V2[starch.data$V3 == "SD"]
names(starch.content.sd) <- starch.data$V1[starch.data$V3 == "SD"]
rain.sd.starch <- rain(starch.content.sd[sorted.zts], deltat=4, period=24, verbose=FALSE, nr.series=2)

library(circacompare)
starch.time.points <- seq(from=0,by=4,length.out = 12)
circacompare.data <- data.frame(time=c(starch.time.points,
                                       starch.time.points),
                                measure=c(starch.content.ld,
                                          starch.content.sd),
                                group=c(rep("LD",12),rep("SD",12)))

result.starch <- circacompare(x = circacompare.data, 
                              col_time = "time", 
                              col_group = "group", 
                              col_outcome = "measure",
                              alpha_threshold = 1)
result.starch.values <- result.starch$summary$value
names(result.starch.values) <- result.starch$summary$parameter

# SD starch content was significantly lower than under LD conditions with 
# a p-value of 1.7 x 10-5. SD starch content oscillates around 29.04 % whereas under
# LD conditions 33.45 %. 
# 
# A significant phase shift of almost 5 hours is detected as a response to photoperiod
# shortening under SD with respect to LD with a p-value of 7.14 x 10-5. Under LD conditions
# phase or the time of highest starch content was detected at around ZT4 whereas under LD
# condtions it was estimated to be at around ZT9. 


plot(x = seq(from=0,by=4,length.out=12),starch.content.ld,type="o",col="blue",lwd=4,ylim=c(19,45),axes=F,xlab="",ylab="")
axis(side=2,lwd=3,at = seq(from=25,to=45,by=5),cex.axis=1.2,las=2)

polygon(x = c(0,16,16,0),y=c(22, 22, 24, 24),lwd=2,border="blue")
polygon(x = c(16,24,24,16),y=c(22, 22, 24, 24),lwd=2,border="blue",col="blue") 
polygon(x = c(24,40,40,24),y=c(22, 22, 24, 24),lwd=2,border="blue")
polygon(x = c(40,48,48,40),y=c(22, 22, 24, 24),lwd=2,border="blue",col="blue") 

lines(x = seq(from=0,by=4,length.out=12),y=starch.content.sd, type="o",col="red",lwd=4)

polygon(x = c(0,8,8,0),y=c(19, 19, 21, 21),lwd=2,border="red")
polygon(x = c(8,24,24,8),y=c(19, 19, 21, 21),lwd=2,border="red",col="red") 
polygon(x = c(24,32,32,24),y=c(19, 19, 21, 21),lwd=2,border="red")
polygon(x = c(32,48,48,32),y=c(19, 19, 21, 21),lwd=2,border="red",col="red") #380 375




sd.mean.starch <- (starch.content.sd[1:6] + starch.content.sd[7:12])/2
sd.mean.starch.splines <- spline(x=seq(from=0,to=24,by=4),y=c(sd.mean.starch, sd.mean.starch[1]),n = 100)
col_fun_sd_starch = colorRamp2(c(min(sd.mean.starch.splines$y),max(sd.mean.starch.splines$y)), 
                                c("black", "white"))
cols <- col_fun_sd_starch(x = sd.mean.starch.splines$y)
plot(x=1,y=1,col="white",xlim=c(-1,25),ylim=c(-10,10),axes=F,xlab="",ylab="")
for(i in 1:99)
{
 polygon(x = c(sd.mean.starch.splines$x[i],sd.mean.starch.splines$x[i+1],
               sd.mean.starch.splines$x[i+1],sd.mean.starch.splines$x[i]),
         y=c(-2.5, -2.5, 2.5, 2.5),lwd=2,border=cols[i],col=cols[i]) 
}




ld.mean.starch <- (starch.content.ld[1:6] + starch.content.ld[7:12])/2
ld.mean.starch.splines <- spline(x=seq(from=0,to=24,by=4),y=c(ld.mean.starch, ld.mean.starch[1]),n = 100)
col_fun_ld_starch = colorRamp2(c(min(ld.mean.starch.splines$y),max(ld.mean.starch.splines$y)), 
                               c("black", "white"))
cols <- col_fun_ld_starch(x = ld.mean.starch.splines$y)
plot(x=1,y=1,col="white",xlim=c(-1,25),ylim=c(-10,10),axes=F,xlab="",ylab="")
for(i in 1:99)
{
 polygon(x = c(ld.mean.starch.splines$x[i],ld.mean.starch.splines$x[i+1],
               ld.mean.starch.splines$x[i+1],ld.mean.starch.splines$x[i]),
         y=c(-2.5, -2.5, 2.5, 2.5),lwd=2,border=cols[i],col=cols[i]) 
}







png(filename = "CAH_ostta01g05830.png", width = 410, height = 280)
heatmap.gene.protein(gene.id = "ostta01g05830")
dev.off()

png(filename = "PEPC_ostta12g02960.png", width = 410, height = 280)
heatmap.gene.protein(gene.id = "ostta12g02960")
dev.off()





lineplot.ld.gene.protein <- function(gene)
{
 ld.current.gene.expression <- scale(ld.mean.gene.expression[gene,],center = T,scale = T)[,1]
 ld.current.protein.abundance <- scale(ld.mean.protein.abundance[gene,],center = T,scale = T)[,1]
 
 plot(x = seq(from=0,to=24,by=4),y = c(ld.current.gene.expression,ld.current.gene.expression[1]),
      type="o",lwd=3,axes=F,xlab="",ylab="", ylim=c(-3,2.2),
      cex.lab=1.3,main="",cex.main=2,pch=20,cex=1.5,col="lightblue",lty=2)
 
 arrows(x0 = seq(from=0,to=24,by=4),y0 =  c(ld.current.gene.expression,ld.current.gene.expression[1]) - 
         c(scaled.sem.ld.gene.expression[gene,],scaled.sem.ld.gene.expression[gene,1]),
        x1 = seq(from=0,to=24,by=4), y1 = c(ld.current.gene.expression,ld.current.gene.expression[1]) + 
         c(scaled.sem.ld.gene.expression[gene,],scaled.sem.ld.gene.expression[gene,1]),
        length = 0.05,angle = 90,code = 3,col = "lightblue",lwd=3)
 
 lines(x = seq(from=0,to=24,by=4),y = c(ld.current.protein.abundance,ld.current.protein.abundance[1]),
       lwd=3, col="blue",type="o",pch=20,cex=1.5)
 arrows(x0 = seq(from=0,to=24,by=4),y0 =  c(ld.current.protein.abundance,ld.current.protein.abundance[1]) - 
         c(scaled.sem.ld.protein.abundance[gene,],scaled.sem.ld.protein.abundance[gene,1]),
        x1 = seq(from=0,to=24,by=4), y1 = c(ld.current.protein.abundance,ld.current.protein.abundance[1]) + 
         c(scaled.sem.ld.protein.abundance[gene,],scaled.sem.ld.protein.abundance[gene,1]),length = 0.05,angle = 90,code = 3,col = "blue",lwd=3)
 axis(side=2,lwd=3,at = seq(from=-2,to=2,by=1),labels=rep("",5),cex.axis=1.2,las=2)
 
 polygon(x = c(0,16,16,0),y=c(-3, -3, -2.6, -2.6),lwd=2,border="blue")
 polygon(x = c(16,24,24,16),y=c(-3, -3, -2.6, -2.6),lwd=2,border="blue",col="blue") #490 460
}

lineplot.sd.gene.protein <- function(gene)
{
 sd.current.gene.expression <- scale(sd.mean.gene.expression[gene,],center = T,scale = T)[,1]
 sd.current.protein.abundance <- scale(sd.mean.protein.abundance[gene,],center = T,scale = T)[,1]
 
 plot(x = seq(from=0,to=24,by=4),y = c(sd.current.gene.expression,sd.current.gene.expression[1]),
      type="o",lwd=3,axes=F,xlab="",ylab="", ylim=c(-3,2.2),
      cex.lab=1.3,main="",cex.main=2,pch=20,cex=1.5,col="salmon",lty=2)
 
 arrows(x0 = seq(from=0,to=24,by=4),y0 =  c(sd.current.gene.expression,sd.current.gene.expression[1]) - 
         c(scaled.sem.sd.gene.expression[gene,],scaled.sem.sd.gene.expression[gene,1]),
        x1 = seq(from=0,to=24,by=4), y1 = c(sd.current.gene.expression,sd.current.gene.expression[1]) + 
         c(scaled.sem.sd.gene.expression[gene,],scaled.sem.sd.gene.expression[gene,1]),
        length = 0.05,angle = 90,code = 3,col = "salmon",lwd=3)
 
 lines(x = seq(from=0,to=24,by=4),y = c(sd.current.protein.abundance,sd.current.protein.abundance[1]),
       lwd=3, col="red",type="o",pch=20,cex=1.5)
 arrows(x0 = seq(from=0,to=24,by=4),y0 =  c(sd.current.protein.abundance,sd.current.protein.abundance[1]) - 
         c(scaled.sem.sd.protein.abundance[gene,],scaled.sem.sd.protein.abundance[gene,1]),
        x1 = seq(from=0,to=24,by=4), y1 = c(sd.current.protein.abundance,sd.current.protein.abundance[1]) + 
         c(scaled.sem.sd.protein.abundance[gene,],scaled.sem.sd.protein.abundance[gene,1]),length = 0.05,angle = 90,code = 3,col = "red",lwd=3)
 axis(side=2,lwd=3,at = seq(from=-2,to=2,by=1),labels=rep("",5),cex.axis=1.2,las=2)
 
 polygon(x = c(0,8,8,0),y=c(-3, -3, -2.6, -2.6),lwd=2,border="red")
 polygon(x = c(8,24,24,8),y=c(-3, -3, -2.6, -2.6),lwd=2,border="red",col="red") #490 460
}



dev.off()
lineplot.ld.gene.protein(gene = "ostta10g00260")
lineplot.sd.gene.protein(gene = "ostta10g00260")


lineplot.ld.gene.protein(gene = "ostta10g00920")
lineplot.sd.gene.protein(gene = "ostta10g00920")


enzymatic.activity <- read.table(file="physiological_data/nitrate_assimilation/enzymatic_activities.csv",header=T,sep=",")


mean.enzymatic.activity <- (enzymatic.activity$NR[enzymatic.activity$Condition == "SD"][1:6] + enzymatic.activity$NR[enzymatic.activity$Condition == "SD"][7:12])/2
mean.enzymatic.activity <- (enzymatic.activity$GS[enzymatic.activity$Condition == "SD"][1:6] + enzymatic.activity$GS[enzymatic.activity$Condition == "SD"][7:12])/2

mean.enzymatic.activity.splines <- spline(x=seq(from=0,to=24,by=4),y=c(mean.enzymatic.activity, mean.enzymatic.activity[1]),n = 100)
col_fun_enzymatic_activity = colorRamp2(c(min(mean.enzymatic.activity.splines$y),max(mean.enzymatic.activity.splines$y)), 
                               c("black", "yellow"))
cols <- col_fun_enzymatic_activity(x = mean.enzymatic.activity.splines$y)
plot(x=1,y=1,col="white",xlim=c(-1,25),ylim=c(-10,10),axes=F,xlab="",ylab="")
for(i in 1:99)
{
 polygon(x = c(mean.enzymatic.activity.splines$x[i],mean.enzymatic.activity.splines$x[i+1],
               mean.enzymatic.activity.splines$x[i+1],mean.enzymatic.activity.splines$x[i]),
         y=c(-2.5, -2.5, 2.5, 2.5),lwd=2,border=cols[i],col=cols[i]) 
}



plot((enzymatic.activity$NR[enzymatic.activity$Condition == "LD"][1:6] + enzymatic.activity$NR[enzymatic.activity$Condition == "LD"][7:12])/2,type="o")
plot((enzymatic.activity$NR[enzymatic.activity$Condition == "SD"][1:6] + enzymatic.activity$NR[enzymatic.activity$Condition == "SD"][7:12])/2,type="o")

plot((enzymatic.activity$GS[enzymatic.activity$Condition == "LD"][1:6] + enzymatic.activity$GS[enzymatic.activity$Condition == "LD"][7:12])/2,type="o")
plot((enzymatic.activity$GS[enzymatic.activity$Condition == "SD"][1:6] + enzymatic.activity$GS[enzymatic.activity$Condition == "SD"][7:12])/2,type="o")


lineplot.ld.gene.protein(gene = "ostta01g05020")
lineplot.sd.gene.protein(gene = "ostta01g05020")

# Cytosolic ribsomes
lineplot.ld.gene.protein(gene = "ostta01g00485")
lineplot.sd.gene.protein(gene = "ostta01g00485")
heatmap.gene.protein(gene.id = "ostta01g00485")

lineplot.ld.gene.protein(gene = "ostta15g01160")
lineplot.sd.gene.protein(gene = "ostta15g01160")
heatmap.gene.protein(gene.id = "ostta15g01160")

lineplot.ld.gene.protein(gene = "ostta09g04010")
lineplot.sd.gene.protein(gene = "ostta09g04010")
heatmap.gene.protein(gene.id = "ostta09g04010")

lineplot.ld.gene.protein(gene = "ostta03g01600")
lineplot.sd.gene.protein(gene = "ostta03g01600")
heatmap.gene.protein(gene.id = "ostta03g01600")

lineplot.ld.gene.protein(gene = "ostta04g01860")
lineplot.sd.gene.protein(gene = "ostta04g01860")
heatmap.gene.protein(gene.id = "ostta04g01860")

lineplot.ld.gene.protein(gene = "ostta04g01760")
lineplot.sd.gene.protein(gene = "ostta04g01760")
heatmap.gene.protein(gene.id = "ostta04g01760")

##Plastid ribosomes
lineplot.ld.gene.protein(gene = "ostta15g02200")
lineplot.sd.gene.protein(gene = "ostta15g02200")
heatmap.gene.protein(gene.id = "ostta15g02200")


ld.proteins <- read.table(file = "physiological_data/proteins/proteins_ld.tsv",header = T)
ld.proteins.mean <- (ld.proteins$Mean[1:6] + ld.proteins$Mean[7:12] + ld.proteins$Mean[13:18])/3

sd.proteins <- read.table(file = "physiological_data/proteins/proteins_sd.tsv",header = T)
sd.proteins.mean <- (sd.proteins$Mean[1:6] + sd.proteins$Mean[7:12] + sd.proteins$Mean[13:18])/3

plot(ld.proteins.mean,type="o")
plot(sd.proteins.mean,type="o")

nrow(sd.proteins)/6
plot(sd.proteins$Mean,type="o")


plot(x = seq(from=0,by=4,length.out=12),starch.content.ld,type="o",col="blue",lwd=4,ylim=c(19,45),axes=F,xlab="",ylab="")
axis(side=2,lwd=3,at = seq(from=25,to=45,by=5),cex.axis=1.2,las=2)

polygon(x = c(0,16,16,0),y=c(22, 22, 24, 24),lwd=2,border="blue")
polygon(x = c(16,24,24,16),y=c(22, 22, 24, 24),lwd=2,border="blue",col="blue") 
polygon(x = c(24,40,40,24),y=c(22, 22, 24, 24),lwd=2,border="blue")
polygon(x = c(40,48,48,40),y=c(22, 22, 24, 24),lwd=2,border="blue",col="blue") 
