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

plot.ld.ll <- function(gene.id, gene.name, gene.expression,)
{
    ld.zt <- paste("ld",paste0("zt",sprintf(fmt = "%02d",seq(from=0,to=20,by=4))),sep="_")
    current.gene.expression.ld.ll <- gene.expression[gene.id,c(paste(ld.zt,1,sep="_"),paste(ld.zt,2,sep="_"),paste(ld.zt,3,sep="_"),paste(ld.zt,4,sep="_"),paste(ld.zt,5,sep="_"))]

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

plot.ld.ll(gene.id = "ostta04g00450",gene.name = "MCM5", ld.ll.gene.expression)




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
  
  png(filename = paste(paste(c("ld_five_days",gene.name,gene.id),collapse="_"),".png",sep=""))
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
  
  dev.off()
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


res.sort <- sort(x = results.ld.ll$pVal,decreasing = F,index.return=T)

int.genes <- rownames(results.ld.ll)[res.sort$ix]

i <- 1
plot.ld.ll(gene.id = int.genes[i],gene.name = "???", gene.expression)
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
