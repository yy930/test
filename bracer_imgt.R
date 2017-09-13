library(data.table)
library(stringr)
options(stringsAsFactors = FALSE)
# #get input for imgt
# bracer1<-read.table(file = "/Users/yingy_adm/Desktop/get_ig_genes/full_bracerseq_1.txt",header=T,sep = "\t")
# bracer2<-read.table(file = "/Users/yingy_adm/Desktop/get_ig_genes/full_bracerseq_2.txt",header=T,sep = "\t")
# bracer3<-read.table(file = "/Users/yingy_adm/Desktop/get_ig_genes/full_bracerseq_3.txt",header=T,sep = "\t")
# bracer4<-read.table(file = "/Users/yingy_adm/Desktop/get_ig_genes/full_bracerseq_4.txt",header=T,sep = "\t")
# bracer5<-read.table(file = "/Users/yingy_adm/Desktop/get_ig_genes/full_bracerseq_5.txt",header=T,sep = "\t")
# bracer<-rbind(bracer1,bracer2,bracer3,bracer4,bracer5)
# bracer$ID<-paste(bracer$CellID,bracer$ReceptorID,bracer$TPM,sep="|")
# bracer<-bracer[,c("ID","full_sequence")]
# bracer$ID<-paste(">",bracer$ID,sep = "")
# write.table(file="/Users/yingy_adm/Desktop/bracer/bracer.txt",bracer, sep="\t", col.names=F)
celltype <- read.delim("~/Documents/sc_count/celltype_all.txt")
bracer_imgtout<-read.table(file = "/Users/yingy_adm/Desktop/bracer/bracer_imgt.txt",header=T,sep = "\t")
flag<-str_split_fixed(bracer_imgtout$ID,"[|]",3)
bracer_imgtout$CellID<-flag[,1]
bracer_imgtout$ReceptorID<-flag[,2]
bracer_imgtout$TPM<-flag[,3]
bracer_imgtout<-bracer_imgtout[,c("CellID","ReceptorID","TPM","V.GENE.and.allele","J.GENE.and.allele","CDR3.IMGT","CDR3.IMGT..AA.")]
bracer_imgtout$V.GENE.and.allele<-str_split_fixed(bracer_imgtout$V.GENE.and.allele,"[*]",2)[,1]
bracer_imgtout$J.GENE.and.allele<-str_split_fixed(bracer_imgtout$J.GENE.and.allele,"[*]",2)[,1]
bracer_imgtout<-bracer_imgtout[bracer_imgtout$TPM>1,]
colnames(bracer_imgtout)<-c("CellID","ReceptorID","TPM","V.GENE","J.GENE","CDR3nt","CDR3aa")

bracer1<-bracer_imgtout[grepl("^[^1234]",bracer_imgtout$CellID),]
bracer2<-bracer_imgtout[grepl("^1",bracer_imgtout$CellID),]
bracer3<-bracer_imgtout[grepl("^2",bracer_imgtout$CellID),]
bracer4<-bracer_imgtout[grepl("^3",bracer_imgtout$CellID),]
bracer5<-bracer_imgtout[grepl("^4",bracer_imgtout$CellID),]
bracerlist<-list(bracer1,bracer2,bracer3,bracer4,bracer5)

#aggretate BCRs in bracer by cellid&receptorid
for (i in 1:length(bracerlist)){
  #bracerlist[[i]]$flag<-paste(bracerlist[[i]]$CellID,bracerlist[[i]]$ReceptorID,bracerlist[[i]]$full_sequence,sep='|')
  bracerlist[[i]]$flag<-paste(bracerlist[[i]]$CellID,bracerlist[[i]]$ReceptorID,bracerlist[[i]]$CDR3nt,bracerlist[[i]]$CDR3aa,sep='|')
  bracerlist[[i]]$TPM<-as.numeric(bracerlist[[i]]$TPM)
  agg_bracer<-aggregate(TPM~flag,data=bracerlist[[i]],FUN=sum)
  agg_bracer$CellID<-str_split_fixed(agg_bracer$flag, "[|]", 4)[,1]
  #agg_bracer$ReceptorID_seq<-str_split_fixed(agg_bracer$flag, "[|]", 2)[,2]
  agg_bracer$ReceptorID<-str_split_fixed(agg_bracer$flag, "[|]", 4)[,2]#########'3'for 2
  agg_bracer$CDR3nt<-str_split_fixed(agg_bracer$flag, "[|]", 4)[,3]
  agg_bracer$CDR3aa<-str_split_fixed(agg_bracer$flag, "[|]", 4)[,4]
  bracerlist[[i]]<-agg_bracer[ , -which(names(agg_bracer)=="flag" )]
}

bracer1<-bracerlist[[1]]
bracer2<-bracerlist[[2]]
bracer3<-bracerlist[[3]]
bracer4<-bracerlist[[4]]
bracer5<-bracerlist[[5]]
##################try CDR3nt#####no correction#####
for (i in 1:length(bracerlist)){
  bracerlist[[i]]$flag<-paste(bracerlist[[i]]$CellID,bracerlist[[i]]$CDR3nt,sep = "|")
  bracerlist[[i]]$TPM<-as.numeric(bracerlist[[i]]$TPM)
  agg_bracer<-aggregate(TPM~flag,data=bracerlist[[i]],FUN=sum)
  collapsed_id<-aggregate(ReceptorID~flag,data=bracerlist[[i]],FUN=c)
  agg_bracer<-merge(agg_bracer, collapsed_id, by="flag")
  agg_bracer$CellID<-str_split_fixed(agg_bracer$flag, "[|]", 2)[,1]
  agg_bracer$CDR3nt<-str_split_fixed(agg_bracer$flag, "[|]", 2)[,2]
  bracerlist[[i]]<-agg_bracer[ , -which(names(agg_bracer)=="flag" )]
}

common_nt<-matrix(nrow = 5,ncol = 5)
for (i in 1:5){
  for (j in 1:5){
    common_nt[i,j]=length(intersect(bracerlist[[i]]$CDR3nt,bracerlist[[j]]$CDR3nt))
  }
}

#     [,1] [,2] [,3] [,4] [,5]
# [1,]  161   29   30   34   32
# [2,]   29  115   88   35   30
# [3,]   30   88  421   40   36
# [4,]   34   35   40  202  115
# [5,]   32   30   36  115  360

ntbcr1<-as.data.frame(table(bracerlist[[1]]$CDR3nt))#ReceptorID_seq
ntbcr2<-as.data.frame(table(bracerlist[[2]]$CDR3nt))
ntbcr3<-as.data.frame(table(bracerlist[[3]]$CDR3nt))
ntbcr4<-as.data.frame(table(bracerlist[[4]]$CDR3nt))
ntbcr5<-as.data.frame(table(bracerlist[[5]]$CDR3nt))
colnames(ntbcr1)[1]=colnames(ntbcr2)[1]=colnames(ntbcr3)[1]=colnames(ntbcr4)[1]=colnames(ntbcr5)[1]="CDR3nt"
ntbcrlist<-list(ntbcr1,ntbcr2,ntbcr3,ntbcr4,ntbcr5)

#get cells in which each BCR has been detected
for (plateID in 1:5){
  i=0
  for (nt in ntbcrlist[[plateID]]$CDR3nt){
    i=i+1
    ind_cell <- bracerlist[[plateID]]$CDR3nt==nt
    cells <- bracerlist[[plateID]]$CellID[ind_cell]
    cells<-paste0(cells,collapse = ",")
    ntbcrlist[[plateID]]$cells[i] <- cells
    tpms <- paste(bracerlist[[plateID]]$TPM[ind_cell],collapse = ",")
    ntbcrlist[[plateID]]$TPMs[i] <- tpms
    receptor<-paste(unique(bracerlist[[plateID]]$ReceptorID[ind_cell]),collapse = ",")
    ntbcrlist[[plateID]]$Receptors[i] <- receptor
  }
  rownames(ntbcrlist[[plateID]])<-ntbcrlist[[plateID]]$ReceptorID
}

ntbcr1 <- data.frame(lapply(ntbcrlist[[1]], as.character), stringsAsFactors=FALSE)
ntbcr2 <- data.frame(lapply(ntbcrlist[[2]], as.character), stringsAsFactors=FALSE)
ntbcr3 <- data.frame(lapply(ntbcrlist[[3]], as.character), stringsAsFactors=FALSE)
ntbcr4 <- data.frame(lapply(ntbcrlist[[4]], as.character), stringsAsFactors=FALSE)
ntbcr5 <- data.frame(lapply(ntbcrlist[[5]], as.character), stringsAsFactors=FALSE)

ntbcr1<-ntbcr1[ntbcr1$Freq>1,]#10 TCRs
ntbcr23<-merge(ntbcr2,ntbcr3,by=c("CDR3nt","CDR3nt"),all = FALSE,suffixes=c("_plate2","_plate3"))
ntbcr23$sumfreq<-as.numeric(ntbcr23$Freq_plate2)+as.numeric(ntbcr23$Freq_plate3)
ntbcr23$cells<-paste(ntbcr23$cells_plate2,ntbcr23$cells_plate3,sep = ",")
ntbcr23$TPMs<-paste(ntbcr23$TPMs_plate2,ntbcr23$TPMs_plate3,sep = ",")
ntbcr45<-merge(ntbcr4,ntbcr5,by=c("CDR3nt","CDR3nt"),all = FALSE,suffixes=c("_plate4","_plate5"))
ntbcr45$sumfreq<-as.numeric(ntbcr45$Freq_plate4)+as.numeric(ntbcr45$Freq_plate5)
ntbcr45$cells<-paste(ntbcr45$cells_plate4,ntbcr45$cells_plate5,sep = ",")
ntbcr45$TPMs<-paste(ntbcr45$TPMs_plate4,ntbcr45$TPMs_plate5,sep = ",")

#identity sources
ntbcr1_source<-find_source_plate(ntbcr1,type = "P")
ntbcr23_source<-find_source_plate(ntbcr23,type = "P")
ntbcr45_source<-find_source_plate(ntbcr45,type = "P")

#create sourcelist of their highest expressed a,b vgenes to check if a cell is source for v gene
ntbracer<-rbind(bracerlist[[1]],bracerlist[[2]],bracerlist[[3]],
               bracerlist[[4]],bracerlist[[5]])
if (exists("ntbcrsourcelist")) {rm(ntbcrsourcelist)}
for (cell in celltype$well){
  h_id=NA
  l_id=NA
  k_id=NA
  max_h_id=NA
  max_l_id=NA
  id<-which(ntbracer$CellID==cell)
  h_id<-grep("^IGHV",ntbracer$ReceptorID[id])
  l_id<-grep("^IGLV",ntbracer$ReceptorID[id])
  k_id<-grep("^IGKV",ntbracer$ReceptorID[id])
  if(length(h_id)!=0)
  {max_h_id<-which.max(ntbracer$TPM[id][h_id])}
  if(length(l_id)!=0|length(k_id)!=0)
  {max_l_id<-which.max(c(ntbracer$TPM[id][l_id],ntbracer$TPM[id][k_id]))}
  if (length(max_h_id)>0|(length(max_l_id)>0))
  {
    record<-c(cell,ntbracer$CDR3nt[id][h_id][max_h_id],ntbracer$CDR3nt[id][c(l_id,k_id)][max_l_id])
    ifelse(exists("ntbcrsourcelist"),ntbcrsourcelist<-rbind(ntbcrsourcelist,record),ntbcrsourcelist<-record)
  }
}
rownames(ntbcrsourcelist)<-ntbcrsourcelist[,1]

#apply check_source function to filter potential source
ntbcr1_source$filted_source<-NA
for (i in 1:length(ntbcr1_source$CDR3nt)){
  if(ntbcr1_source$source[i]!="")
  {ntbcr1_source$filted_source[i]<-check_source(ntbcr1_source$CDR3nt[i],ntbcr1_source$source[i],use="ntbcrsourcelist")}
}

ntbcr23_source$filted_source<-NA
for (i in 1:length(ntbcr23_source$CDR3nt)){
  if(ntbcr23_source$source[i]!="")
  {ntbcr23_source$filted_source[i]<-check_source(ntbcr23_source$CDR3nt[i],ntbcr23_source$source[i],use="ntbcrsourcelist")}
}

ntbcr45_source$filted_source<-NA
for (i in 1:length(ntbcr45_source$CDR3nt)){
  if(ntbcr45_source$source[i]!="")
  {ntbcr45_source$filted_source[i]<-check_source(ntbcr45_source$CDR3nt[i],ntbcr45_source$source[i],use="ntbcrsourcelist")}
}

array_nt1<-get_array(ntbcr1_source,n_plate="F",plate=1,use="nt")#use c(1,2,3,4,6,8)
array_nt23<-get_array(ntbcr23_source,n_plate="T",plate=23,use="nt")
array_nt45<-get_array(ntbcr45_source,n_plate="T",plate=45,use="nt")

###################################

collapse_xnt_mismatch<-function(bcr,num_nt){
  bcr$V<-str_split_fixed(bcr$ReceptorID, "_", 3)[,1]
  bcr$N<-str_split_fixed(bcr$ReceptorID, "_", 3)[,2]
  bcr$J<-str_split_fixed(bcr$ReceptorID, "_", 3)[,3]
  bcr$V_J<-paste(bcr$V,bcr$J,sep='_')
  #bcr$topTPM<-unlist(lapply(str_split(bcr$TPM,","),function(x) max(as.numeric(x))))
  for (i in unique(bcr$V_J)){
    ids<-which(bcr$V_J==i)
    if (length(ids)>1){
      #      pairp=combn(bcr$ReceptorID[ids],2)
      #ind_maxtpm<-which.max(bcr[ids,]$topTPM)
      ind_maxtpm<-which.max(bcr[ids,]$TPM)
      #      ref<-bcr[ids,]$N[ind_maxtpm]
      #ref<-bcr[ind_maxtpm,]$full_sequence??????
      ref<-bcr[ids,]$CDR3nt[ind_maxtpm]
      refaa<-bcr[ids,]$CDR3aa[ind_maxtpm]
      ids <- ids[-ind_maxtpm]
      for (j in ids)
      {
        #      if (adist(bcr$N[j],ref)<=num_nt) 
        if (adist(bcr$CDR3nt[j],ref)<=num_nt)
        {
          #           bcr$N[j]<-ref #change the nt seq
          bcr$CDR3nt[j]<-ref
          bcr$CDR3aa[j]<-refaa
        }
      }  
    }
  }
  bcr<-bcr[,-which(names(bcr)%in%c("V","N","J","V_J"))]
  return(bcr)}

collapsed_bracer1<-collapse_xnt_mismatch(bracer1,num_nt=2)
collapsed_bracer2<-collapse_xnt_mismatch(bracer2,num_nt=2)
collapsed_bracer3<-collapse_xnt_mismatch(bracer3,num_nt=2)
collapsed_bracer4<-collapse_xnt_mismatch(bracer4,num_nt=2)
collapsed_bracer5<-collapse_xnt_mismatch(bracer5,num_nt=2)

collapsed_bracerlist<-list(collapsed_bracer1,collapsed_bracer2,collapsed_bracer3,collapsed_bracer4,collapsed_bracer5)

for (i in 1:length(collapsed_bracerlist)){
  collapsed_bracerlist[[i]]$flag<-paste(collapsed_bracerlist[[i]]$CellID,collapsed_bracerlist[[i]]$CDR3aa,sep = "|")
  collapsed_bracerlist[[i]]$TPM<-as.numeric(collapsed_bracerlist[[i]]$TPM)
  agg_bracer<-aggregate(TPM~flag,data=collapsed_bracerlist[[i]],FUN=sum)
  collapsed_id<-aggregate(ReceptorID~flag,data=collapsed_bracerlist[[i]],FUN=c)
  agg_bracer<-merge(agg_bracer, collapsed_id, by="flag")
  agg_bracer$CellID<-str_split_fixed(agg_bracer$flag, "[|]", 2)[,1]
  agg_bracer$CDR3aa<-str_split_fixed(agg_bracer$flag, "[|]", 2)[,2]
  collapsed_bracerlist[[i]]<-agg_bracer[ , -which(names(agg_bracer)=="flag" )]
}

common_collapsed_bcr<-matrix(nrow = 5,ncol = 5)
for (i in 1:5){
  for (j in 1:5){
    common_collapsed_bcr[i,j]=length(intersect(collapsed_bracerlist[[i]]$CDR3aa,collapsed_bracerlist[[j]]$CDR3aa))
    #common_collapsed_bcr[i,j]=length(intersect(collapsed_bracerlist[[i]]$CDR3nt,collapsed_bracerlist[[j]]$CDR3nt))
  }
}

#      [,1] [,2] [,3] [,4] [,5]
# [1,]  124   22   24   26   24
# [2,]   22   85   68   30   25
# [3,]   24   68  198   35   29
# [4,]   26   30   35  155   92
# [5,]   24   25   29   92  214

#freq for each BCR in plate
sbcr1<-as.data.frame(table(collapsed_bracerlist[[1]]$CDR3aa))#ReceptorID_seq
sbcr2<-as.data.frame(table(collapsed_bracerlist[[2]]$CDR3aa))
sbcr3<-as.data.frame(table(collapsed_bracerlist[[3]]$CDR3aa))
sbcr4<-as.data.frame(table(collapsed_bracerlist[[4]]$CDR3aa))
sbcr5<-as.data.frame(table(collapsed_bracerlist[[5]]$CDR3aa))
colnames(sbcr1)[1]=colnames(sbcr2)[1]=colnames(sbcr3)[1]=colnames(sbcr4)[1]=colnames(sbcr5)[1]="CDR3aa"
sbcrlist<-list(sbcr1,sbcr2,sbcr3,sbcr4,sbcr5)

#get cells in which each BCR has been detected
for (plateID in 1:5){
  i=0
  for (aa in sbcrlist[[plateID]]$CDR3aa){
    i=i+1
    ind_cell <- collapsed_bracerlist[[plateID]]$CDR3aa==aa
    cells <- collapsed_bracerlist[[plateID]]$CellID[ind_cell]
    cells<-paste0(cells,collapse = ",")
    sbcrlist[[plateID]]$cells[i] <- cells
    tpms <- paste(collapsed_bracerlist[[plateID]]$TPM[ind_cell],collapse = ",")
    sbcrlist[[plateID]]$TPMs[i] <- tpms
    receptor<-paste(unique(collapsed_bracerlist[[plateID]]$ReceptorID[ind_cell]),collapse = ",")
    sbcrlist[[plateID]]$Receptors[i] <- receptor
  }
  #rownames(bcrlist[[plateID]])<-bcrlist[[plateID]]$ReceptorID
}

sbcr1 <- data.frame(lapply(sbcrlist[[1]], as.character), stringsAsFactors=FALSE)
sbcr2 <- data.frame(lapply(sbcrlist[[2]], as.character), stringsAsFactors=FALSE)
sbcr3 <- data.frame(lapply(sbcrlist[[3]], as.character), stringsAsFactors=FALSE)
sbcr4 <- data.frame(lapply(sbcrlist[[4]], as.character), stringsAsFactors=FALSE)
sbcr5 <- data.frame(lapply(sbcrlist[[5]], as.character), stringsAsFactors=FALSE)

sbcr1<-sbcr1[sbcr1$Freq>1,]#10 TCRs
sbcr23<-merge(sbcr2,sbcr3,by=c("CDR3aa","CDR3aa"),all = FALSE,suffixes=c("_plate2","_plate3"))
sbcr23$sumfreq<-as.numeric(sbcr23$Freq_plate2)+as.numeric(sbcr23$Freq_plate3)
sbcr23$cells<-paste(sbcr23$cells_plate2,sbcr23$cells_plate3,sep = ",")
sbcr23$TPMs<-paste(sbcr23$TPMs_plate2,sbcr23$TPMs_plate3,sep = ",")
sbcr45<-merge(sbcr4,sbcr5,by=c("CDR3aa","CDR3aa"),all = FALSE,suffixes=c("_plate4","_plate5"))
sbcr45$sumfreq<-as.numeric(sbcr45$Freq_plate4)+as.numeric(sbcr45$Freq_plate5)
sbcr45$cells<-paste(sbcr45$cells_plate4,sbcr45$cells_plate5,sep = ",")
sbcr45$TPMs<-paste(sbcr45$TPMs_plate4,sbcr45$TPMs_plate5,sep = ",")

#identity sources
sbcr1_source<-find_source_plate(sbcr1,type = "P")
sbcr23_source<-find_source_plate(sbcr23,type = "P")
sbcr45_source<-find_source_plate(sbcr45,type = "P")

#create sourcelist of their highest expressed a,b vgenes to check if a cell is source for v gene
sbracer<-rbind(collapsed_bracerlist[[1]],collapsed_bracerlist[[2]],collapsed_bracerlist[[3]],
                collapsed_bracerlist[[4]],collapsed_bracerlist[[5]])
if (exists("sbcrsourcelist")) {rm(sbcrsourcelist)}
for (cell in celltype$well){
  h_id=NA
  l_id=NA
  k_id=NA
  max_h_id=NA
  max_l_id=NA
  id<-which(sbracer$CellID==cell)
  h_id<-grep("^IGHV",sbracer$ReceptorID[id])
  l_id<-grep("^IGLV",sbracer$ReceptorID[id])
  k_id<-grep("^IGKV",sbracer$ReceptorID[id])
  if(length(h_id)!=0)
  {max_h_id<-which.max(sbracer$TPM[id][h_id])}
  if(length(l_id)!=0|length(k_id)!=0)
  {max_l_id<-which.max(c(sbracer$TPM[id][l_id],sbracer$TPM[id][k_id]))}
  if (length(max_h_id)>0|(length(max_l_id)>0))
  {
    record<-c(cell,sbracer$CDR3aa[id][h_id][max_h_id],sbracer$CDR3aa[id][c(l_id,k_id)][max_l_id])
    ifelse(exists("sbcrsourcelist"),sbcrsourcelist<-rbind(sbcrsourcelist,record),sbcrsourcelist<-record)
  }
}
rownames(sbcrsourcelist)<-sbcrsourcelist[,1]

#apply check_source function to filter potential source
sbcr1_source$filted_source<-NA
for (i in 1:length(sbcr1_source$CDR3aa)){
  if(sbcr1_source$source[i]!="")
  {sbcr1_source$filted_source[i]<-check_source(sbcr1_source$CDR3aa[i],sbcr1_source$source[i],use="sbcrsourcelist")}
}

sbcr23_source$filted_source<-NA
for (i in 1:length(sbcr23_source$CDR3aa)){
  if(sbcr23_source$source[i]!="")
  {sbcr23_source$filted_source[i]<-check_source(sbcr23_source$CDR3aa[i],sbcr23_source$source[i],use="sbcrsourcelist")}
}

sbcr45_source$filted_source<-NA
for (i in 1:length(sbcr45_source$CDR3aa)){
  if(sbcr45_source$source[i]!="")
  {sbcr45_source$filted_source[i]<-check_source(sbcr45_source$CDR3aa[i],sbcr45_source$source[i],use="sbcrsourcelist")}
}

array_aa1<-get_array(sbcr1_source,n_plate="F",plate=1,use="aa")#use c(1,2,3,4,6,8)
array_aa23<-get_array(sbcr23_source,n_plate="T",plate=23,use="aa")
array_aa45<-get_array(sbcr45_source,n_plate="T",plate=45,use="aa")
array_bcr1<-get_array(sbcr1_source,n_plate="F",plate=1,use="bcr")#use c(5,6,7,8)
array_bcr23<-get_array(bcr23_source,n_plate="T",plate=23,use="bcr")
array_bcr45<-get_array(bcr45_source,n_plate="T",plate=45,use="bcr")

