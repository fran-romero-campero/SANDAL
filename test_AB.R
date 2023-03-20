######Comparacion con circacompare de ondas LD y SD proteomica global


ld.normalized.proteomic.data <- read.table(file = "swath_proteomic_data/LD_normalization/Quantile-normalized.txt",header = T,sep="\t",as.is=T)
head(ld.normalized.proteomic.data)
sum(is.na(ld.normalized.proteomic.data))
for(i in 1:nrow(ld.normalized.proteomic.data))
{
  if(sum(is.na(ld.normalized.proteomic.data[i,])) != 0)
  { 
    na.points <- colnames(ld.normalized.proteomic.data)[which(is.na(ld.normalized.proteomic.data[i,]))]
    for(j in 1:length(na.points))
    {
      zt <- strsplit(na.points[j],split="_")[[1]][1]  
      
      imputed.value <- mean(as.numeric(ld.normalized.proteomic.data[i,paste(zt,1:9,sep="_")]),na.rm = T)
      if(!is.nan(imputed.value))
      {
        ld.normalized.proteomic.data[i,na.points[j]] <- imputed.value
      }
    }
  }
}

sum(is.na(ld.normalized.proteomic.data))

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

sum(is.na(ld.normalized.proteomic.data))

ld.normalized.proteomic.data[is.na(ld.normalized.proteomic.data)] <- min(ld.normalized.proteomic.data[,2:55],na.rm = T)
head(ld.normalized.proteomic.data)

ld.protein.ids <- ld.normalized.proteomic.data$ProtID

ld.normalized.proteomic.data <- as.matrix(ld.normalized.proteomic.data[,2:55])
head(ld.normalized.proteomic.data)
rownames(ld.normalized.proteomic.data) <- ld.protein.ids

############################SD lectura e imputacion
sd.normalized.proteomic.data <- read.table(file = "swath_proteomic_data/SD_normalization/Quantile-normalized.txt",header = T,sep="\t",as.is=T)
head(sd.normalized.proteomic.data)
sum(is.na(sd.normalized.proteomic.data))
for(i in 1:nrow(sd.normalized.proteomic.data))
{
  if(sum(is.na(sd.normalized.proteomic.data[i,])) != 0)
  { 
    na.points <- colnames(sd.normalized.proteomic.data)[which(is.na(sd.normalized.proteomic.data[i,]))]
    for(j in 1:length(na.points))
    {
      zt <- strsplit(na.points[j],split="_")[[1]][1]  
      
      imputed.value <- mean(as.numeric(sd.normalized.proteomic.data[i,paste(zt,1:9,sep="_")]),na.rm = T)
      if(!is.nan(imputed.value))
      {
        sd.normalized.proteomic.data[i,na.points[j]] <- imputed.value
      }
    }
  }
}

sum(is.na(sd.normalized.proteomic.data))

for(i in 1:nrow(sd.normalized.proteomic.data))
{
  if(sum(is.na(sd.normalized.proteomic.data[i,])) != 0)
  {
    na.points <- colnames(sd.normalized.proteomic.data)[which(is.na(sd.normalized.proteomic.data[i,]))]
    
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
      
      sd.normalized.proteomic.data[i,paste(zts[j],1:9,sep="_")] <-  mean(as.numeric(c(sd.normalized.proteomic.data[i,paste(previous.time.point,1:9,sep="_")],
                                                                                      sd.normalized.proteomic.data[i,paste(next.time.point,1:9,sep="_")])))
    }
  }
}

sum(is.na(sd.normalized.proteomic.data))

sd.normalized.proteomic.data[is.na(sd.normalized.proteomic.data)] <- min(sd.normalized.proteomic.data[,2:55],na.rm = T)
head(sd.normalized.proteomic.data)

sd.protein.ids <- sd.normalized.proteomic.data$ProtID
sd.normalized.proteomic.data <- as.matrix(sd.normalized.proteomic.data[,2:55])
head(sd.normalized.proteomic.data)
rownames(sd.normalized.proteomic.data) <- sd.protein.ids

#####################seleccionar rhythmicos
rhythmic.sd<-read.table(file = "protein_sets/rhythmic_proteins_sd.txt")
rhythmic.sd <- rhythmic.sd$V1
rhythmic.ld<-read.table(file = "protein_sets/rhythmic_proteins_ld.txt")
rhythmic.ld <- rhythmic.ld$V1

ld.sd.rhythmic.proteins <- intersect(rhythmic.ld, rhythmic.sd)
length(rhythmic.ld)
length(rhythmic.sd)
length(ld.sd.rhythmic.proteins)



########################circacomprare

circacompare.ld.sd.prot <- matrix(nrow=length(ld.sd.rhythmic.proteins),ncol=15)
rownames(circacompare.ld.sd.prot) <- ld.sd.rhythmic.proteins

for(i in 1:length(ld.sd.rhythmic.proteins))
{
  gene.id <- ld.sd.rhythmic.proteins[i]
  
  zts <- rep(paste("zt", seq(from=0,by=4,to=20), sep = ""),9)
  order <- paste(zts,rep(1:9, each=6), sep = "_")
  ld.expression.i <- ld.normalized.proteomic.data[gene.id,order]
  sd.expression.i <- sd.normalized.proteomic.data[gene.id,order]
  
  time.points <- seq(from=0,by=4,length.out = 54)
  
  ld.sd.df <- data.frame(time=c(time.points,time.points),
                         measure=c(ld.expression.i, sd.expression.i),
                         group=c(rep("protein_ld",54),rep("protein_sd",54)))
  
  out.i <- circacompare(x = ld.sd.df, col_time = "time", col_group = "group", col_outcome = "measure",alpha_threshold = 1)
  circacompare.ld.sd.prot[i,] <- out.i[[2]][,2]
}

colnames(circacompare.ld.sd.prot) <- out.i[[2]][,1]

sum(circacompare.ld.sd.prot[,"Amplitude difference estimate" ] < 0)
sum(circacompare.ld.sd.prot[,"Amplitude difference estimate" ] < 0) / length(ld.sd.rhythmic.proteins)
sum(circacompare.ld.sd.prot[circacompare.ld.sd.prot[,"Amplitude difference estimate" ] < 0,"P-value for amplitude difference" ] < 0.05)
genes.decrease.amplitude.effect.ld.sd <- names(which(circacompare.ld.sd.prot[circacompare.ld.sd.prot[,"Amplitude difference estimate" ] < 0,"P-value for amplitude difference" ] < 0.05)) 
sum(circacompare.ld.sd.prot[circacompare.ld.sd.prot[,"Amplitude difference estimate" ] < 0,"P-value for amplitude difference" ] < 0.05)/length(ld.sd.rhythmic.proteins)
87/190


par(lwd=3)
boxplot(circacompare.ld.sd.prot[,"protein_ld amplitude estimate" ],circacompare.ld.sd.prot[,"protein_sd amplitude estimate" ],outline=F,col=c("blue","red"),names=c("LD","SD"),cex.axis=2)
wilcox.test(x = circacompare.ld.sd.prot[,"protein_ld amplitude estimate" ], y = circacompare.ld.sd.prot[,"protein_sd amplitude estimate" ], alternative = "greater")$p.value
#0.12

sum(circacompare.ld.sd.prot[,"Mesor difference estimate" ] < 0)
sum(circacompare.ld.sd.prot[,"Mesor difference estimate" ] < 0) / length(ld.sd.rhythmic.proteins)
sum(circacompare.ld.sd.prot[circacompare.ld.sd.prot[,"Mesor difference estimate" ] < 0,"P-value for mesor difference" ] < 0.05)
sum(circacompare.ld.sd.prot[circacompare.ld.sd.prot[,"Mesor difference estimate" ] < 0,"P-value for mesor difference" ] < 0.05)/length(ld.sd.rhythmic.proteins)

sum(circacompare.ld.sd.prot[,"Mesor difference estimate" ] > 0)
sum(circacompare.ld.sd.prot[,"Mesor difference estimate" ] > 0) / length(ld.sd.rhythmic.proteins)

sum(circacompare.ld.sd.prot[circacompare.ld.sd.prot[,"Mesor difference estimate" ] > 0,"P-value for mesor difference" ] < 0.05)

sum(circacompare.ld.sd.prot[circacompare.ld.sd.prot[,"Mesor difference estimate" ] > 0,"P-value for mesor difference" ] < 0.05)/length(ld.sd.rhythmic.proteins)

#No clear effect over the mesor was detected as 30.96% of the rhythmic genes significantly increase the mean value around which they oscillate and 19.48% significantly decrease it under short day conditions when compared to long day conditions.  

par(lwd=3)
boxplot(circacompare.ld.sd.prot[,"protein_ld mesor estimate" ],circacompare.ld.sd.prot[,"protein_sd mesor estimate" ],outline=F,col=c("blue","red"),names=c("LD","SD"),cex.axis=2)

sum(circacompare.ld.sd.prot[,"Phase difference estimate" ] < 0)
sum(circacompare.ld.sd.prot[,"Phase difference estimate" ] < 0) / length(ld.sd.rhythmic.proteins)
sum(circacompare.ld.sd.prot[circacompare.ld.sd.prot[,"Phase difference estimate" ] < 0,"P-value for difference in phase" ] < 0.05)
genes.phase.effect.ld.sd <- names(which(circacompare.ld.sd.prot[circacompare.ld.sd.prot[,"Phase difference estimate" ] < 0,"P-value for difference in phase" ] < 0.05))
sum(circacompare.ld.sd.prot[circacompare.ld.sd.prot[,"Phase difference estimate" ] < 0,"P-value for difference in phase" ] < 0.05)/length(ld.sd.rhythmic.proteins)
109/179
