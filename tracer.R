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

#read data, subset, rename colomn
options(stringsAsFactors = FALSE)
celltype <- read.delim("~/Documents/sc_count/celltype_all.txt")
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

#number of detected receptors in each cell
#table(strsplit(paste(tcrlist[[2]]$cells,collapse = ","),","))

library(tcR)
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
  m$label <- as.vector(matrix(celltype$cell.type[289:480],nrow = 16, ncol=12, byrow=T))
  p <- ggplot(m, aes(x = variable, y = name, fill = value))
  p <- p + geom_tile(aes(fill = value), colour = "white")
  if (.text) {
    p <- p + geom_text(aes(fill = value, label = label), size = .size.text)
#    p <- p + geom_label(data=subset(mtcars, wt > 4 | mpg > 25),
                        #geom_label(aes(fill = factor(cyl)), colour = "white", fontface = "bold")
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

