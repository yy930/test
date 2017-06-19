library(data.table)
#function to rename cells
getfile <- function(taskid){
  letter = c("A","B","C","D","E","F","G","H")
  num = c("1","2","3","4","5","6","7","8","9","10","11","12")
  if (as.numeric(taskid)<97)
  {if (as.numeric(taskid) %% 12 == 0)
  { l=letter[as.numeric(taskid)/12]
  n="12"} 
    else{
      l = letter[as.numeric(taskid)/12+1]
      n = num[as.numeric(taskid)%%12]}
    return (paste("1",l,n,sep=""))
  }
  if (as.numeric(taskid)>96){
    smalln=as.numeric(taskid)-96
    if (smalln%%12 == 0){
      l = letter[smalln/12]
      n="12"
    }
    else{
      l = letter[smalln/12+1]
      n = num[smalln%%12]     
    }
    return (paste("2",l,n,sep=""))
  }
}

getfilelist<-function(input_file_list){
  out_file_list=vector('character')
  for (infile in input_file_list){
    out_file_list<-c(out_file_list,getfile(infile))
  }
  return (out_file_list)   
}

options(stringsAsFactors = FALSE)
celltype <- read.delim("~/Documents/sc_count/celltype_all.txt")
# #read data, subset, rename colomn
transtpm_sep<-read.table(file = "/Users/yingy_adm/Desktop/tracer-2017/transtpm_sep.txt",row.names=1,header = T, sep = " ",check.names =F)
transtpm_nov<-read.table(file = "/Users/yingy_adm/Desktop/tracer-2017/transtpm_nov.txt",row.names=1,header = T, sep = " ",check.names =F)
transtpm_feb<-read.table(file = "/Users/yingy_adm/Desktop/tracer-2017/transtpm_feb.txt",row.names=1,header = T, sep = " ",check.names =F)
# transtpm_sep<-transtpm_sep[grepl("^TRAV|^TRBV|^TRDV|^TRGV",rownames(transtpm_sep)),]
# transtpm_nov<-transtpm_nov[grepl("^TRAV|^TRBV|^TRDV|^TRGV",rownames(transtpm_nov)),]
# transtpm_feb<-transtpm_feb[grepl("^TRAV|^TRBV|^TRDV|^TRGV",rownames(transtpm_feb)),]
# colnames(transtpm_sep)<-celltype$well[1:96]
# colnames(transtpm_nov)<-celltype$well[97:288]
# colnames(transtpm_feb)<-celltype$well[289:480]

transtpm_sep$TranscriptID<-rownames(transtpm_sep)
transtpm_sep<-melt(transtpm_sep,id="TranscriptID")
transtpm_sep<-transtpm_sep[transtpm_sep$value>0,]

transtpm_nov$TranscriptID<-rownames(transtpm_nov)
transtpm_nov<-melt(transtpm_nov,id="TranscriptID")
transtpm_nov<-transtpm_nov[transtpm_nov$value>0,]

transtpm_feb$TranscriptID<-rownames(transtpm_feb)
transtpm_feb<-melt(transtpm_feb,id="TranscriptID")
transtpm_feb<-transtpm_feb[transtpm_feb$value>0,]

colnames(transtpm_sep)=colnames(transtpm_nov)=colnames(transtpm_feb)=c("TranscriptID","CellID","TPM")

v1<-as.data.frame(table(transtpm_sep$TranscriptID))
vplate2<-transtpm_nov[grepl("^1",transtpm_nov$CellID),]
vplate2$CellID<-gsub("^1","",vplate2$CellID)
v2<-as.data.frame(table(vplate2$TranscriptID))
vplate3<-transtpm_nov[grepl("^2",transtpm_nov$CellID),]
vplate3$CellID<-gsub("^2","",vplate3$CellID)
v3<-as.data.frame(table(vplate3$TranscriptID))
vplate4<-transtpm_feb[grepl("^3",transtpm_feb$CellID),]
vplate4$CellID<-gsub("^3","",vplate4$CellID)
v4<-as.data.frame(table(vplate4$TranscriptID))
vplate5<-transtpm_feb[grepl("^4",transtpm_feb$CellID),]
vplate5$CellID<-gsub("^4","",vplate5$CellID)
v5<-as.data.frame(table(vplate5$TranscriptID))
colnames(v1)[1]=colnames(v2)[1]=colnames(v3)[1]=colnames(v4)[1]=colnames(v5)[1]="V_GENE"
vlist<-list(v1,v2,v3,v4,v5)
vplatelist<-list(transtpm_sep,vplate2,vplate3,vplate4,vplate5)

#get cells in which each TCR has been detected
for (plateID in 1:5){
  i=0
  for (gene in vlist[[plateID]]$V_GENE){
    i=i+1
    ind_cell <- vplatelist[[plateID]]$TranscriptID==gene
    cells <- vplatelist[[plateID]]$CellID[ind_cell]
    cells<-paste0(cells,collapse = ",")
    vlist[[plateID]]$cells[i] <- cells
    tpms <- paste(vplatelist[[plateID]]$TPM[ind_cell],collapse = ",")
    vlist[[plateID]]$TPMs[i] <- tpms
  }
}

v1 <- data.frame(lapply(vlist[[1]], as.character), stringsAsFactors=FALSE)
v2 <- data.frame(lapply(vlist[[2]], as.character), stringsAsFactors=FALSE)
v3 <- data.frame(lapply(vlist[[3]], as.character), stringsAsFactors=FALSE)
v4 <- data.frame(lapply(vlist[[4]], as.character), stringsAsFactors=FALSE)
v5 <- data.frame(lapply(vlist[[5]], as.character), stringsAsFactors=FALSE)

#function for a receptor's expression in plate 
expr_v_mat<-function(plateID,geneID){
  cells<-vlist[[plateID]][geneID,3]#example
  tpms<-vlist[[plateID]][geneID,4]
  mat_plate<-matrix(nrow=8,ncol=12,dimnames=list(c("A","B","C","D","E","F","G","H"),as.character(1:12)))
  vec_cell<-strsplit(cells,split = ",")[[1]]
  vec_cellrow<-substr(vec_cell,1,1)
  vec_cellcol<-substring(vec_cell,2)
  vec_tpm <-as.numeric(strsplit(tpms,split = ",")[[1]])
  #cell_tpm_mat<-data.frame(vec_cellrow,vec_cellcol,vec_tpm)
  for (i in (1:length(vec_cell))){
    row_id=which(rownames(mat_plate)==vec_cellrow[i])
    col_id=which(colnames(mat_plate)==vec_cellcol[i])
    mat_plate[row_id,col_id]<-vec_tpm[i]
  }
  mat_plate<-as.numeric(mat_plate)
  mat<-matrix(log2(as.numeric(mat_plate)+1),8,12)
  rownames(mat)<-c("A","B","C","D","E","F","G","H")
  colnames(mat)<-as.character(1:12)
  return(mat)}

v23<-merge(v2,v3,by=c("V_GENE","V_GENE"),all = FALSE,suffixes=c("_plate2","_plate3"))
v23$sumfreq<-as.numeric(v23$Freq_plate2)+as.numeric(v23$Freq_plate3)
v23<-v23[order(-v23[,8]),]

v45<-merge(v4,v5,by=c("V_GENE","V_GENE"),all = FALSE,suffixes=c("_plate4","_plate5"))
v45$sumfreq<-as.numeric(v45$Freq_plate4)+as.numeric(v45$Freq_plate5)
v45<-v45[order(-v45[,8]),]

vplotlist1<-list()
v1<-v1[v1$Freq>1,]
v1<-v1[order(-as.numeric(v1$Freq)),]
#for (i in 1:dim(v1)[1]){
for (i in 1:12){
  receptor=v1$V_GENE[i]
  plateID=1
  geneID=which(vlist[[1]]$V_GENE==receptor)
  mat<-expr_v_mat(plateID,geneID)
  vplotlist1[[i]]<-visheatmap(mat,.title=paste(receptor," in Plate",plateID),.labs = c("", ""),.legend = "log2(TPM)")}
pdf("/Users/yingy_adm/Desktop/v_heatmap1_celltype.pdf",height=65,width=8)
multiplot(vplotlist1[[1]], vplotlist1[[2]], vplotlist1[[3]], vplotlist1[[4]],vplotlist1[[5]],
          vplotlist1[[6]],vplotlist1[[7]],vplotlist1[[8]],vplotlist1[[9]],vplotlist1[[10]],
          vplotlist1[[11]],vplotlist1[[12]])
dev.off()

vplotlist23<-list()
for (i in 1:12){########!
  receptor=v23$V_GENE[i]
  plateID1=2
  geneID=which(v2$V_GENE==receptor)
  mat1<-expr_v_mat(plateID1,geneID)
  plateID2=3
  geneID=which(v3$V_GENE==receptor)
  mat2<-expr_v_mat(plateID2,geneID)
  rownames(mat2)<-paste(rownames(mat2),"*")
  mat<-rbind(mat1,mat2)
  vplotlist23[[i]]<-visheatmap(mat,.title=paste(receptor," in Plate",plateID1," + Plate",plateID2,"*"),.labs = c("", ""),.legend = "log2(TPM)")}
pdf("/Users/yingy_adm/Desktop/v_heatmap23_celltype.pdf",height=110,width=8)
multiplot(vplotlist23[[1]], vplotlist23[[2]], vplotlist23[[3]], vplotlist23[[4]],vplotlist23[[5]],
          vplotlist23[[6]],vplotlist23[[7]],vplotlist23[[8]],vplotlist23[[9]],vplotlist23[[10]],
          vplotlist23[[11]],vplotlist23[[12]])
dev.off()  

vplotlist45<-list()
for (i in 1:12){########!
  receptor=v45$V_GENE[i]
  plateID1=4
  geneID=which(v4$V_GENE==receptor)
  mat1<-expr_v_mat(plateID1,geneID)
  plateID2=5
  geneID=which(v5$V_GENE==receptor)
  mat2<-expr_v_mat(plateID2,geneID)
  rownames(mat2)<-paste(rownames(mat2),"*")
  mat<-rbind(mat1,mat2)
  vplotlist45[[i]]<-visheatmap(mat,.title=paste(receptor," in Plate",plateID1," + Plate",plateID2,"*"),.labs = c("", ""),.legend = "log2(TPM)")}
pdf("/Users/yingy_adm/Desktop/v_heatmap45_celltype.pdf",height=110,width=8)
multiplot(vplotlist45[[1]], vplotlist45[[2]], vplotlist45[[3]], vplotlist45[[4]],vplotlist45[[5]],
          vplotlist45[[6]],vplotlist45[[7]],vplotlist45[[8]],vplotlist45[[9]],vplotlist45[[10]],
          vplotlist45[[11]],vplotlist45[[12]])
dev.off()  
#################################################

tracer<-read.table(file = "/Users/yingy_adm/Desktop/tracer-2017/output.txt",header = T, sep = "\t")
tracer$CellID<-gsub("_","",tracer$CellID)
plate1<-subset(tracer,CellID %in% celltype$well[1:96])
plate2<-subset(tracer,CellID %in% as.character(1:96))
plate2$CellID<-gsub("^1","",getfilelist(plate2$CellID))
plate3<-subset(tracer,CellID %in% as.character(97:192))
plate3$CellID<-gsub("^2","",getfilelist(plate3$CellID))
plate4<-subset(tracer,CellID %in% celltype$well[97:192])
plate4$CellID<-gsub("^1","",plate4$CellID)
plate5<-subset(tracer,CellID %in% celltype$well[193:288])
plate5$CellID<-gsub("^2","",plate5$CellID)
platelist<-list(plate1,plate2,plate3,plate4,plate5)

#freq for each TCR in plate
tcr1<-as.data.frame(table(plate1$ReceptorID))
tcr2<-as.data.frame(table(plate2$ReceptorID))
tcr3<-as.data.frame(table(plate3$ReceptorID))
tcr4<-as.data.frame(table(plate4$ReceptorID))
tcr5<-as.data.frame(table(plate5$ReceptorID))
colnames(tcr1)[1]=colnames(tcr2)[1]=colnames(tcr3)[1]=colnames(tcr4)[1]=colnames(tcr5)[1]="TCR"
tcrlist<-list(tcr1,tcr2,tcr3,tcr4,tcr5)

#get cells in which each TCR has been detected
for (plateID in 1:5){
  i=0
  for (receptor in tcrlist[[plateID]]$TCR){
    i=i+1
    ind_cell <- platelist[[plateID]]$ReceptorID==receptor
    cells <- platelist[[plateID]]$CellID[ind_cell]
    cells<-paste0(cells,collapse = ",")
    tcrlist[[plateID]]$cells[i] <- cells
    tpms <- paste(platelist[[plateID]]$TPM[ind_cell],collapse = ",")
    tcrlist[[plateID]]$TPMs[i] <- tpms
  }
}

tcr1 <- data.frame(lapply(tcrlist[[1]], as.character), stringsAsFactors=FALSE)
tcr2 <- data.frame(lapply(tcrlist[[2]], as.character), stringsAsFactors=FALSE)
tcr3 <- data.frame(lapply(tcrlist[[3]], as.character), stringsAsFactors=FALSE)
tcr4 <- data.frame(lapply(tcrlist[[4]], as.character), stringsAsFactors=FALSE)
tcr5 <- data.frame(lapply(tcrlist[[5]], as.character), stringsAsFactors=FALSE)

#number of detected receptors and V genes in each cell
#table(strsplit(paste(tcrlist[[2]]$cells,collapse = ","),","))
#number of detected receptors and V genes in all plates
length(unique(c(tcr1$TCR,tcr2$TCR,tcr3$TCR,tcr4$TCR,tcr5$TCR)))#173
length(unique(c(v1$V_GENE,v2$V_GENE,v3$V_GENE,v4$V_GENE,v5$V_GENE)))#183

library(tcR)
library(scater)
#function for a receptor's expression in plate 
expr_tcr_mat<-function(plateID,tcrID){
  cells<-tcrlist[[plateID]][tcrID,3]#example
  tpms<-tcrlist[[plateID]][tcrID,4]
  mat_plate<-matrix(nrow=8,ncol=12,dimnames=list(c("A","B","C","D","E","F","G","H"),as.character(1:12)))
  vec_cell<-strsplit(cells,split = ",")[[1]]
  vec_cellrow<-substr(vec_cell,1,1)
  vec_cellcol<-substring(vec_cell,2)
  vec_tpm <-as.numeric(strsplit(tpms,split = ",")[[1]])
  #cell_tpm_mat<-data.frame(vec_cellrow,vec_cellcol,vec_tpm)
  for (i in (1:length(vec_cell))){
    row_id=which(rownames(mat_plate)==vec_cellrow[i])
    col_id=which(colnames(mat_plate)==vec_cellcol[i])
    mat_plate[row_id,col_id]<-vec_tpm[i]
  }
  mat_plate<-as.numeric(mat_plate)
  mat<-matrix(log2(as.numeric(mat_plate)+1),8,12)
  rownames(mat)<-c("A","B","C","D","E","F","G","H")
  colnames(mat)<-as.character(1:12)
  return(mat)}

tcr23<-merge(tcr2,tcr3,by=c("TCR","TCR"),all = FALSE,suffixes=c("_plate2","_plate3"))
tcr23$sumfreq<-as.numeric(tcr23$Freq_plate2)+as.numeric(tcr23$Freq_plate3)
tcr23<-tcr23[order(-tcr23[,8]),]

tcr45<-merge(tcr4,tcr5,by=c("TCR","TCR"),all = FALSE,suffixes=c("_plate4","_plate5"))
tcr45$sumfreq<-as.numeric(tcr45$Freq_plate4)+as.numeric(tcr45$Freq_plate5)
tcr45<-tcr45[order(-tcr45[,8]),]

plotlist1<-list()
tcr1<-tcr1[tcr1$Freq>1,]
for (i in 1:dim(tcr1)[1]){
  receptor=tcr1$TCR[i]
  plateID=1
  tcrID=which(tcrlist[[1]]$TCR==receptor)
  mat<-expr_tcr_mat(plateID,tcrID)
  plotlist1[[i]]<-visheatmap(mat,.title=paste(receptor," in Plate",plateID),.labs = c("", ""),.legend = "log2(TPM)")}
pdf("/Users/yingy_adm/Desktop/tcr_heatmap1_celltype.pdf",height=60,width=8)
multiplot(plotlist1[[1]], plotlist1[[2]], plotlist1[[3]], plotlist1[[4]],plotlist1[[5]],
          plotlist1[[6]],plotlist1[[7]],plotlist1[[8]],plotlist1[[9]],plotlist1[[10]])
dev.off()

plotlist23<-list()
for (i in 1:dim(tcr23)[1]){
  receptor=tcr23$TCR[i]
  plateID1=2
  tcrID=which(tcr2$TCR==receptor)
  mat1<-expr_tcr_mat(plateID1,tcrID)
  plateID2=3
  tcrID=which(tcr3$TCR==receptor)
  mat2<-expr_tcr_mat(plateID2,tcrID)
  rownames(mat2)<-paste(rownames(mat2),"*")
  mat<-rbind(mat1,mat2)
  plotlist23[[i]]<-visheatmap(mat,.title=paste(receptor," in Plate",plateID1," + Plate",plateID2,"*"),.labs = c("", ""),.legend = "log2(TPM)")}
pdf("/Users/yingy_adm/Desktop/tcr_heatmap23_celltype.pdf",height=110,width=8)
multiplot(plotlist23[[1]], plotlist23[[2]], plotlist23[[3]], plotlist23[[4]],plotlist23[[5]],
          plotlist23[[6]],plotlist23[[7]],plotlist23[[8]],plotlist23[[9]],plotlist23[[10]],
          plotlist23[[11]],plotlist23[[12]])
dev.off()  

plotlist45<-list()
for (i in 1:dim(tcr45)[1]){
  receptor=tcr45$TCR[i]
  plateID1=4
  tcrID=which(tcr4$TCR==receptor)
  mat1<-expr_tcr_mat(plateID1,tcrID)
  plateID2=5
  tcrID=which(tcr5$TCR==receptor)
  mat2<-expr_tcr_mat(plateID2,tcrID)
  rownames(mat2)<-paste(rownames(mat2),"*")
  mat<-rbind(mat1,mat2)
  plotlist45[[i]]<-visheatmap(mat,.title=paste(receptor," in Plate",plateID1," + Plate",plateID2,"*"),.labs = c("", ""),.legend = "log2(TPM)")}
pdf("/Users/yingy_adm/Desktop/tcr_heatmap45_celltype.pdf",height=100,width=10)
multiplot(plotlist45[[1]], plotlist45[[2]], plotlist45[[3]], plotlist45[[4]],plotlist45[[5]],plotlist45[[6]],plotlist45[[7]],plotlist45[[8]],plotlist45[[9]],plotlist45[[10]])
dev.off()  

#test if common TCRs from different plates
intersect(tcrlist[[1]]$TCR,tcrlist[[2]]$TCR)#1
intersect(tcrlist[[2]]$TCR,tcrlist[[3]]$TCR)#18
intersect(tcrlist[[4]]$TCR,tcrlist[[5]]$TCR)#15
intersect(tcrlist[[1]]$TCR,tcrlist[[3]]$TCR)#1
intersect(tcrlist[[4]]$TCR,tcrlist[[3]]$TCR)#1

visheatmap<-function(.data, .title = "Number of shared clonotypes", .labs = c("Sample", "Sample"),
                     .legend = "Shared clonotypes", .na.value = NA, .text = T, .scientific = FALSE,
                     .signif.digits = 4, .size.text = 4, .no.legend = F, .no.labs = F) 
{
  if (has.class(.data, "data.frame")) {
    names <- .data[, 1]
    .data <- as.matrix(.data[, -1])
    row.names(.data) <- names
  }
  if (is.null(colnames(.data))) {
    colnames(.data) <- paste0("C", 1:ncol(.data))
  }
  if (is.null(row.names(.data))) {
    row.names(.data) <- paste0("C", 1:nrow(.data))
  }
  .data[is.na(.data)] <- .na.value
  tmp <- as.data.frame(.data)
  tmp$name <- row.names(.data)
  m <- melt(tmp, id.var = c("name"))
  m[, 1] <- factor(m[, 1], levels = rev(rownames(.data)))
  m[, 2] <- factor(m[, 2], levels = colnames(.data))
  .cg <- .colourblind.gradient(min(m$value), max(m$value))
  #m$label <- format(m$value, scientific = .scientific, digits = .signif.digits)
  m$label <- as.vector(matrix(celltype$cell.type[289:480],nrow = 16, ncol=12, byrow=T))#!!!!!!
  p <- ggplot(m, aes(x = variable, y = name, fill = value))
  p <- p + geom_tile(aes(fill = value), colour = "white")
  if (.text) {
    p <- p + geom_text(aes(fill = value, label = label), size = .size.text)#!!!!!!!!!!
  }
  p <- p + .cg
  p <- p + ggtitle(.title) + guides(fill = guide_colourbar(title = .legend)) + 
    xlab(.labs[1]) + ylab(.labs[2]) + coord_fixed() + theme_linedraw() + 
    theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(expand = c(0, 
                                                                                0)) + scale_y_discrete(expand = c(0, 0))
  if (.no.legend) {
    p <- p + theme(legend.position = "none")
  }
  if (.no.labs) {
    p <- p + theme(axis.title.x = element_blank(), axis.title.y = element_blank())
  }
  p
}

.colourblind.gradient <- function (.min = NA, .max = NA, .colour = F) {
  #   cs <- c("#FFFFD9", "#41B6C4", "#225EA8")
  #   cs <- c("#FFFFBB", "#41B6C4", "#225EA8")
  #   cs <- c("#FFBB00", "#41B6C4", "#225EA8") <- old version
  #   cs <- c("#FF4B20", "#FFB433", "#C6EDEC", "#85CFFF", "#0348A6")
  # cs <- c("#FF4B20", "#FFB433", "#C6FDEC", "#7AC5FF", "#0348A6")
  # scale_fill_gradientn(guide='colourbar', colours=c("#0072B2", "#EEEEEE", "#D55E00")
  
  cs <- c(c("#0072B2", "#EEEEEE", "#D55E00"))
  
  if (!is.na(.min)) {
    if (.colour) {
      scale_colour_gradientn(limits = c(.min, .max), guide='colorbar', colours = cs, na.value = 'grey60')
    } else {
      scale_fill_gradientn(limits = c(.min, .max), guide='colorbar', colours = cs, na.value = 'grey60')
    }
  } else {
    if (.colour) {
      scale_colour_gradientn(colours = cs, na.value = 'grey60')
    } else {
      scale_fill_gradientn(colours = cs, na.value = 'grey60')
    }
  }
}

common_tcr<-matrix(nrow = 5,ncol = 5)
for (i in 1:5){
  for (j in 1:5){
    common_tcr[i,j]=length(intersect(tcrlist[[i]]$TCR,tcrlist[[j]]$TCR))
  }
}

common_v<-matrix(nrow = 5,ncol = 5)
for (i in 1:5){
  for (j in 1:5){
    common_v[i,j]=length(intersect(vlist[[i]]$V_GENE,vlist[[j]]$V_GENE))
  }
}