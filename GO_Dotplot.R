library(data.table)
# GO Breakdown ------------------------------------------------------------
#PHB specific
GO_PHB_specific_BP_UP <- fread("~/Aman/GO/2024-06-04_PHB-specific-Up_BP.txt",sep="\t",header=T,stringsAsFactors = F)
GO_PHB_specific_BP_UP[,`:=`("ONT"=rep("BP",GO_PHB_specific_BP_UP[,.N]),"COMP"=rep("PHB_specific_UP",GO_PHB_specific_BP_UP[,.N]))]
GO_PHB_specific_BP_DOWN <- fread("~/Aman/GO/2024-06-04_PHB-specific-Down_BP.txt",sep="\t",header=T,stringsAsFactors = F)
GO_PHB_specific_BP_DOWN[,`:=`("ONT"=rep("BP",GO_PHB_specific_BP_DOWN[,.N]),"COMP"=rep("PHB_specific_DOWN",GO_PHB_specific_BP_DOWN[,.N]))]
GO_PHB_specific_CC_UP <- fread("~/Aman/GO/2024-06-04_PHB-specific-Up_CC.txt",sep="\t",header=T,stringsAsFactors = F)
GO_PHB_specific_CC_UP[,`:=`("ONT"=rep("CC",GO_PHB_specific_CC_UP[,.N]),"COMP"=rep("PHB_specific_UP",GO_PHB_specific_CC_UP[,.N]))]
GO_PHB_specific_CC_DOWN <- fread("~/Aman/GO/2024-06-04_PHB-specific-Down_CC.txt",sep="\t",header=T,stringsAsFactors = F)
GO_PHB_specific_CC_DOWN[,`:=`("ONT"=rep("CC",GO_PHB_specific_CC_DOWN[,.N]),"COMP"=rep("PHB_specific_DOWN",GO_PHB_specific_CC_DOWN[,.N]))]
GO_PHB_specific_MF_UP <- fread("~/Aman/GO/2024-06-04_PHB-specific-Up_MF.txt",sep="\t",header=T,stringsAsFactors = F)
GO_PHB_specific_MF_UP[,`:=`("ONT"=rep("MF",GO_PHB_specific_MF_UP[,.N]),"COMP"=rep("PHB_specific_UP",GO_PHB_specific_MF_UP[,.N]))]
GO_PHB_specific_MF_DOWN <- fread("~/Aman/GO/2024-06-04_PHB-specific-Down_MF.txt",sep="\t",header=T,stringsAsFactors = F)
GO_PHB_specific_MF_DOWN[,`:=`("ONT"=rep("MF",GO_PHB_specific_MF_DOWN[,.N]),"COMP"=rep("PHB_specific_DOWN",GO_PHB_specific_MF_DOWN[,.N]))]

#CNA specific
GO_CNA_specific_BP_UP <- fread("~/Aman/GO/2024-06-04_CNA-specific-Up_BP.txt",sep="\t",header=T,stringsAsFactors = F)
GO_CNA_specific_BP_UP[,`:=`("ONT"=rep("BP",GO_CNA_specific_BP_UP[,.N]),"COMP"=rep("CNA_specific_UP",GO_CNA_specific_BP_UP[,.N]))]
GO_CNA_specific_BP_DOWN <- fread("~/Aman/GO/2024-06-04_CNA-specific-Down_BP.txt",sep="\t",header=T,stringsAsFactors = F)
GO_CNA_specific_BP_DOWN[,`:=`("ONT"=rep("BP",GO_CNA_specific_BP_DOWN[,.N]),"COMP"=rep("CNA_specific_DOWN",GO_CNA_specific_BP_DOWN[,.N]))]
GO_CNA_specific_CC_UP <- fread("~/Aman/GO/2024-06-04_CNA-specific-Up_CC.txt",sep="\t",header=T,stringsAsFactors = F)
GO_CNA_specific_CC_UP[,`:=`("ONT"=rep("CC",GO_CNA_specific_CC_UP[,.N]),"COMP"=rep("CNA_specific_UP",GO_CNA_specific_CC_UP[,.N]))]
GO_CNA_specific_CC_DOWN <- fread("~/Aman/GO/2024-06-04_CNA-specific-Down_CC.txt",sep="\t",header=T,stringsAsFactors = F)
GO_CNA_specific_CC_DOWN[,`:=`("ONT"=rep("CC",GO_CNA_specific_CC_DOWN[,.N]),"COMP"=rep("CNA_specific_DOWN",GO_CNA_specific_CC_DOWN[,.N]))]
GO_CNA_specific_MF_UP <- fread("~/Aman/GO/2024-06-04_CNA-specific-Up_MF.txt",sep="\t",header=T,stringsAsFactors = F)
GO_CNA_specific_MF_UP[,`:=`("ONT"=rep("MF",GO_CNA_specific_MF_UP[,.N]),"COMP"=rep("CNA_specific_UP",GO_CNA_specific_MF_UP[,.N]))]
GO_CNA_specific_MF_DOWN <- fread("~/Aman/GO/2024-06-04_CNA-specific-Down_MF.txt",sep="\t",header=T,stringsAsFactors = F)
GO_CNA_specific_MF_DOWN[,`:=`("ONT"=rep("MF",GO_CNA_specific_MF_DOWN[,.N]),"COMP"=rep("CNA_specific_DOWN",GO_CNA_specific_MF_DOWN[,.N]))]

#Mutually-regulated_concordant
GO_Mutually_regulated_concordant_BP_UP <- fread("~/Aman/GO/2024-06-04_Mutually-regulated_concordant-Up_BP.txt",sep="\t",header=T,stringsAsFactors = F)
GO_Mutually_regulated_concordant_BP_UP[,`:=`("ONT"=rep("BP",GO_Mutually_regulated_concordant_BP_UP[,.N]),"COMP"=rep("Mutually-regulated_concordant_UP",GO_Mutually_regulated_concordant_BP_UP[,.N]))]
GO_Mutually_regulated_concordant_BP_DOWN <- fread("~/Aman/GO/2024-06-04_Mutually-regulated_concordant-Down_BP.txt",sep="\t",header=T,stringsAsFactors = F)
GO_Mutually_regulated_concordant_BP_DOWN[,`:=`("ONT"=rep("BP",GO_Mutually_regulated_concordant_BP_DOWN[,.N]),"COMP"=rep("Mutually-regulated_concordant_DOWN",GO_Mutually_regulated_concordant_BP_DOWN[,.N]))]
GO_Mutually_regulated_concordant_CC_UP <- fread("~/Aman/GO/2024-06-04_Mutually-regulated_concordant-Up_CC.txt",sep="\t",header=T,stringsAsFactors = F)
GO_Mutually_regulated_concordant_CC_UP[,`:=`("ONT"=rep("CC",GO_Mutually_regulated_concordant_CC_UP[,.N]),"COMP"=rep("Mutually-regulated_concordant_UP",GO_Mutually_regulated_concordant_CC_UP[,.N]))]
GO_Mutually_regulated_concordant_CC_DOWN <- fread("~/Aman/GO/2024-06-04_Mutually-regulated_concordant-Down_CC.txt",sep="\t",header=T,stringsAsFactors = F)
GO_Mutually_regulated_concordant_CC_DOWN[,`:=`("ONT"=rep("CC",GO_Mutually_regulated_concordant_CC_DOWN[,.N]),"COMP"=rep("Mutually-regulated_concordant_DOWN",GO_Mutually_regulated_concordant_CC_DOWN[,.N]))]
GO_Mutually_regulated_concordant_MF_UP <- fread("~/Aman/GO/2024-06-04_Mutually-regulated_concordant-Up_MF.txt",sep="\t",header=T,stringsAsFactors = F)
GO_Mutually_regulated_concordant_MF_UP[,`:=`("ONT"=rep("MF",GO_Mutually_regulated_concordant_MF_UP[,.N]),"COMP"=rep("Mutually-regulated_concordant_UP",GO_Mutually_regulated_concordant_MF_UP[,.N]))]
GO_Mutually_regulated_concordant_MF_DOWN <- fread("~/Aman/GO/2024-06-04_Mutually-regulated_concordant-Down_MF.txt",sep="\t",header=T,stringsAsFactors = F)
GO_Mutually_regulated_concordant_MF_DOWN[,`:=`("ONT"=rep("MF",GO_Mutually_regulated_concordant_MF_DOWN[,.N]),"COMP"=rep("Mutually-regulated_concordant_DOWN",GO_Mutually_regulated_concordant_MF_DOWN[,.N]))]


Combined_GO <- rbind(GO_PHB_specific_BP_UP,GO_PHB_specific_BP_DOWN,GO_PHB_specific_CC_UP,GO_PHB_specific_CC_DOWN,GO_PHB_specific_MF_UP,GO_PHB_specific_MF_DOWN,
                     GO_CNA_specific_BP_UP,GO_CNA_specific_BP_DOWN,GO_CNA_specific_CC_UP,GO_CNA_specific_CC_DOWN,GO_CNA_specific_MF_UP,GO_CNA_specific_MF_DOWN,
                     GO_Mutually_regulated_concordant_BP_UP,GO_Mutually_regulated_concordant_BP_DOWN,GO_Mutually_regulated_concordant_CC_UP,
                     GO_Mutually_regulated_concordant_CC_DOWN,GO_Mutually_regulated_concordant_MF_UP,GO_Mutually_regulated_concordant_MF_DOWN)

Combined_GO <- Combined_GO[,.(GO.ID,Genes_in_Term)]
#Combined_GO <- Combined_GO[-which(duplicated(Combined_GO[,1])),] 
setkeyv(Combined_GO,"GO.ID")






#FOCUSED TERMS LIST
FOCUS_TERM <- fread("~/Aman/GO/Selected by Aman/Aman_GO_Focus_terms.txt",sep="\t",header=T,stringsAsFactors = F)
setkeyv(FOCUS_TERM,"GO.ID")


# Dotplots for GO ---------------------------------------------------------

Sizemap <- function(FE,breaks=c(-3,-2,-1,0,1,2,3),cex.val = seq(1/3,8/3,1/3)){
  if(is.na(FE)) return(0)
  logFE <- log(FE,2)
  for(i in 1:length(breaks)){
    #print(i)
    if(i==1){
      if(logFE <= breaks[i+1]){
        #print(cex.val[i])
        return(cex.val[i])
      }
    }else if(i == length(breaks)){
      if(logFE >= breaks[i]){
        #print(cex.val[i])
        return(cex.val[i])
      }
      
      
    } else {
      if((logFE > breaks[i]&logFE <= breaks[i+1])){
        #print(cex.val[i])
        return(cex.val[i])
      }
      
    }
    
    
  }
  
  
  
} 

#Sizemap(5,breaks = 0:3,cex.val = 1:5)

Signmap <- function(FE,col=c("black","white")){
  if(is.na(FE)) return(0)
  logFE <- log(FE,2)
  if(logFE<0) return(col[2])
  if(logFE>0) return(col[1])
  if(logFE==0) return(rgb(0,0,0,0))
  
} 

library(RColorBrewer)
blueramp <- colorRampPalette(colors = c("white","darkblue"))

Colmap <- function(FDR,breaks=seq(0,10,1),color.pallete) {
  if(is.na(FDR)) return(0)
  logFDR <- -1*log(FDR,10)
  cols <- color.pallete(length(breaks))
  for(i in 1:length(breaks)){
    if(i == length(breaks)){
      if(logFDR >= breaks[i]){
        return(cols[i])
      }
    } else {
      if((logFDR >= breaks[i]&logFDR < breaks[i+1])){
        return(cols[i])
      }
      
    }
    
  }
  
} 



#grey FDR
png("GO_dot_plot_Grey_FDR.png",width = 9,height = 15,units = "in",res = 300,bg = rgb(0,0,0,0))

layout(matrix(c(rep(1,12),2,3,4),ncol=5,nrow=3))
#main plot area, dot plot
par(mar=c(11.1,22.1,0.1,0.5))
cols <- unique(Combined_GO[,COMP])[c(2,1,6,5,4,3)]
rows <- FOCUS_TERM[order(`Order on figure`,decreasing = T),1]#New_list[rev(c(6,14,16,5,15,12,8,9,3,4,1,2,11,13,7,10)),]
number.cols <- length(cols)  
number.rows <- nrow(rows)
row.pos <- seq(0.5,number.rows*2.5,2.5)
col.pos <- seq(0.5,number.cols*0.5,0.5)

plot(0,0,col=rgb(0,0,0,0),xlim=c(0,max(col.pos)+0.5),ylim=c(-10,max(row.pos)),axes=F,xlab="",ylab="")
for(i in 1:length(row.pos)){
  lines(x = c(0,max(col.pos)+0.5),y=rep(row.pos[i],2),lty=2,col="grey60")
}
for(i in 1:length(col.pos)){
  lines(x = rep(col.pos[i],2), y = c(0,max(row.pos)), lty = 2, col = "grey60")
}

greyramp <- colorRampPalette(colors = c("white","grey20"))

for(row in 1:length(unlist(rows[,1]))){
  for(col in 1:length(cols)){
    print(row)
    print(col)
    print(rows[row])
    print(cols[col])
    print(Combined_GO[.(cols[col],rows[row]),c(FE,two.sideFDR)])
    color <- Colmap(FDR = Combined_GO[.(cols[col],rows[row]),two.sideFDR],breaks = -1*log(c(1,0.05,0.005,0.0005,0.00005,0.000005),10),color.pallete = greyramp)
    size <- Sizemap(FE = abs(Combined_GO[.(cols[col],rows[row]),FE]),breaks = seq(0,3,1),cex.val = c(1:5))
    signFE <- Signmap(FE = Combined_GO[.(cols[col],rows[row]),FE],col = c("black","darkorange2"))
    print(color)
    print(size)
    print(signFE)
    points(x=col.pos[col],y=row.pos[row],pch=21,col=signFE,bg=color,cex=size)
    
  }
}
text(x = col.pos,y=-1.5,labels = gsub(cols,pattern = "_\\w+_",replacement = " "),
     srt=45,adj=1,cex=1.3,font=2)
mtext(2,at = row.pos,text = unlist(FOCUS_TERM[rows[,1],2]),las=2,cex=0.8,font=2)
par(mar=c(0.5,0.5,0.5,0.5))
#FC key
plot(0,0,col=rgb(0,0,0,0),axes=F,xlab="",ylab="",xlim=c(0,2),ylim=c(0,15))
breaks=c(0,1,2,3)
cex.val = seq(1:4)
for(steps in 1:(length(breaks)+2)){
  points(x=0.25,y=(steps*1),pch=21,col="black",bg="grey50",cex=cex.val[steps])
}
text(x = rep(1.25,4),y = seq(1,4,1),labels =  c("0 <= FE < 1","1 <= FE < 2","2 <= FE < 3","FE >= 3     "),font=2,cex=1.3,)
text(x=1,y=5.5,font=2, labels = bquote(Log[2]~Fold~Enrichment),cex=1.3)

color.pallete <- greyramp(length(col.breaks+1))
plot(0,0,col=rgb(0,0,0,0),axes=F,xlab="",ylab="",ylim=c(0,length(col.breaks)/2+1.5),xlim=c(0,2))
text(x=0.65,y=4,font=2, labels = "FDR",cex=1.3)

for(i in seq(1,(length(col.breaks)),1)){
  print(i)
  polygon(x = c(0,0.5,0.5,0,0),y=c(0.5*i,0.5*i,0.5*i+0.5,0.5*i+0.5,0.5*i),col = color.pallete[i],border = "black")
}

text(x = rep(0.75,length(col.breaks)),y = seq(0.5,length(col.breaks)/2+0.5,0.5),
     labels =  c("1","5x10^-2","5x10^-3","5x10^-4","5x10^-5","5x10^-6","<=5x10^-7"),
     las=2,font=2,cex=1.3,adj = c(0,0.5))


dev.off()

