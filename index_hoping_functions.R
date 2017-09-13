#heatmap for plate 1
visheatmap1<-function(.data, .title = "Number of shared clonotypes", .labs = c("Sample", "Sample"),
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
  m$label <- as.vector(matrix(celltype$cell.type[1:96],nrow = 8, ncol=12, byrow=T))#!!!!!!
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

#heatmap for plate 2+3
visheatmap23<-function(.data, .title = "Number of shared clonotypes", .labs = c("Sample", "Sample"),
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
  m$label <- as.vector(matrix(celltype$cell.type[97:288],nrow = 16, ncol=12, byrow=T))#!!!!!!
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

#heatmap for plate 4+5
visheatmap45<-function(.data, .title = "Number of shared clonotypes", .labs = c("Sample", "Sample"),
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

#function to rename files
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

#function to get a list of cells with converted names using function "getfile" 
getfilelist<-function(input_file_list){
  out_file_list=vector('character')
  for (infile in input_file_list){
    out_file_list<-c(out_file_list,getfile(infile))
  }
  return (out_file_list)   
}

#function to aggregate V_GENEs
agg_v_gene<-function(v_cell_tpm){
  v_cell_tpm$V_GENE<-str_split_fixed(v_cell_tpm$V_GENE,pattern ="[*]",n=2)[,1]
  v_cell_tpm$flag<-paste(v_cell_tpm$V_GENE,v_cell_tpm$CellID,sep='|')
  agg_v<-aggregate(TPM~flag,data=v_cell_tpm,FUN=sum)
  agg_v$V_GENE<-str_split_fixed(agg_v$flag,pattern ="[|]",n=2)[,1]
  agg_v$CellID<-str_split_fixed(agg_v$flag,pattern ="[|]",n=2)[,2]
  agg_v<-agg_v[,c("V_GENE","CellID","TPM")]
  return(agg_v)
}

#function to get expression in plate use={"tcr","v","bcr","vj","aa"}
expr_mat<-function(plateID,geneID,use){
  if (use=="v") {
    cells<-vlist[[plateID]][geneID,3]
    tpms<-vlist[[plateID]][geneID,4]
  }
  if (use=="tcr"){
    cells<-tcrlist[[plateID]][tcrID,3]
    tpms<-tcrlist[[plateID]][tcrID,4]
  }
  if (use=="vj"){
    cells<-vjlist[[plateID]][geneID,3]
    tpms<-vjlist[[plateID]][geneID,4]
  }
  if (use=="bcr") {
    cells<-bcrlist[[plateID]][geneID,3]
    tpms<-bcrlist[[plateID]][geneID,4]
  }
  if (use=="aa") {
    cells<-sbcrlist[[plateID]][geneID,3]
    tpms<-sbcrlist[[plateID]][geneID,4]
  }
  mat_plate<-matrix(nrow=8,ncol=12,dimnames=list(c("A","B","C","D","E","F","G","H"),as.character(1:12)))
  vec_cell<-strsplit(cells,split = ",")[[1]]
  vec_cellrow<-gsub("([0-9]+)", "",vec_cell)
  vec_cellcol<-gsub(".*[A-Z]", "", vec_cell)
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
  return (mat)
}

expr_matrix<-function(vec_cell,vec_ratio,plate=1){
  mat_plate<-matrix(nrow=8,ncol=12,dimnames=list(c("A","B","C","D","E","F","G","H"),as.character(1:12)))
  vec_cellrow<-gsub("([0-9]+)$", "", vec_cell)
  vec_cellcol<-gsub("[^\\d]+", "", vec_cell, perl=TRUE)
  for (i in (1:length(vec_cell))){
    row_id=which(rownames(mat_plate)==vec_cellrow[i])
    col_id=which(colnames(mat_plate)==vec_cellcol[i])
    mat_plate[row_id,col_id]<-vec_ratio[i]
  }
  mat_plate[mat_plate==1]<-0
#  mat_plate[is.na(mat_plate)|]<-0
  return (mat_plate)
}

expr_2matrix<-function(vec_cell,vec_ratio,plate=23){
  ifelse(plate==23,
    mat_plate<-matrix(nrow=16,ncol=12,dimnames=list(c("1A","1B","1C","1D","1E","1F","1G","1H","2A","2B","2C","2D","2E","2F","2G","2H"),as.character(1:12))),
    mat_plate<-matrix(nrow=16,ncol=12,dimnames=list(c("3A","3B","3C","3D","3E","3F","3G","3H","4A","4B","4C","4D","4E","4F","4G","4H"),as.character(1:12))))
  vec_cellrow<-gsub("([0-9]+)$", "", vec_cell)
  vec_cellcol<-gsub("[^\\d]+", "", gsub("^([0-9]+)", "",vec_cell), perl=TRUE)
  for (i in (1:length(vec_cell))){
    row_id=which(rownames(mat_plate)==vec_cellrow[i])
    col_id=which(colnames(mat_plate)==vec_cellcol[i])
    mat_plate[row_id,col_id]<-vec_ratio[i]
  }
  mat_plate[mat_plate==1]<-0
  return (mat_plate)
}

#function to get input variable name
get_input_name <- function(v1) {
  deparse(substitute(v1))
}

addstar<-function(x){a=paste("*",x,sep = "")
return(a)}
cat<-function(x){return (paste(x,collapse = ","))}

#function to identify sources!crosspos,celltype=T,tpm>2^9)
# find_source_1plate<-function(tcr){
#   tcr$source<-""
#   cell_list<-str_split(tcr$cells,",")
#   tpm_list<-str_split(tcr$TPMs,",")
#   for (tcrind in 1:length(cell_list)){
#     rows<-substr(cell_list[[tcrind]],1,1)
#     columns<-substring(cell_list[[tcrind]],2)
#     freqtable_row<-table(rows)>1
#     common_rows<-names(freqtable_row[freqtable_row=="TRUE"])
#     freqtable_col<-table(columns)>1
#     common_columns<-names(freqtable_col[freqtable_col=="TRUE"])
#     if (length(common_columns)==0|length(common_rows)==0){
#       maxcellind<-which.max(tpm_list[[tcrind]])
#         cell_type<-celltype$cell.type[celltype$well==paste(rows[maxcellind],columns[maxcellind],sep = "")]##!!
#       if (cell_type=="T"){
#         tcr$source[tcrind]<-paste(tcr$source[tcrind],paste(rows[maxcellind],columns[maxcellind],sep = ""),sep = ",")
#       }
#     }else{
#       # for (cellind in 1:length(cell_list[[tcrind]])){
#       #     cell_type<-celltype$cell.type[celltype$well==paste(rows[cellind],columns[cellind],sep = "")]##!!
#       #   tpm<-tpm_list[[tcrind]][cellind]
#       #   if (rows[cellind]%in%common_rows && columns[cellind]%in%common_columns &&
#       #       cell_type=="T" && 
#       #       tpm%in%as.character(sort(as.numeric(tpm_list[[tcrind]]), decreasing=TRUE)[1:(length(cell_list[[tcrind]])*0.3)]))
#       #   {
#       #     tcr$source[tcrind]<-paste(tcr$source[tcrind],paste(common_rows,common_columns,sep = ""),sep = ",")
#       #   }
#       # }
#       cross_pos<-sort(apply(expand.grid(common_rows,common_columns), 1, paste, collapse = "", sep = ""))
#       cross_cells<-cell_list[[tcrind]][cell_list[[tcrind]]%in%cross_pos]
#       source_cells<-cross_cells[as.numeric(tpm_list[[tcrind]][match(cross_cells,cell_list[[tcrind]])])>2^6]#512
#       source_cells<-source_cells[celltype$cell.type[match(source_cells,celltype$well)]=="T"]
#       tcr$source[tcrind]<-paste(tcr$source[tcrind],paste(source_cells, collapse = ','),sep = ",")
#     }
#   }
#   tcr$source<-gsub("^,","",tcr$source)
#   return(tcr)
# }

#function to identify sources
find_source_plate<-function(tcr,type="T")#type="P"
{tcr$source<-""
 for (tcrind in 1:length(tcr$cells)){
      source_cells=""
      cell_vec<-unlist(str_split(tcr$cells[tcrind],","))
      tpm_vec<-as.numeric(unlist(str_split(tcr$TPMs[tcrind],",")))
      rows<-gsub("([0-9]+)$", "",cell_vec)
      freqtable_row<-table(rows)>1
      common_rows<-names(freqtable_row[freqtable_row=="TRUE"])
      columns<-gsub("[^\\d]+", "",gsub("^([0-9]+)","",cell_vec),perl=TRUE)
      freqtable_col<-table(columns)>1
      common_columns<-names(freqtable_col[freqtable_col=="TRUE"])
      if (length(common_columns)+length(common_rows)>0){#number of common lines
        if (length(common_columns)==0|length(common_rows)==0){#cells in the same column(row)
           maxcellind<-which.max(tpm_vec)
           tpm_source<-as.numeric(tpm_vec[maxcellind])
           well<-paste(rows[maxcellind],columns[maxcellind],sep = "")
           cell_type<-celltype$cell.type[celltype$well==well]##!!
           if (type=="T"){
              if (cell_type=="T" && tpm_source>2^5 && 
                  all(tpm_source>5*tpm_vec[!cell_vec %in% c(well,"A1","1A1","2A1","3A1","4A1")]))
                 {tcr$source[tcrind]<-paste(tcr$source[tcrind],well,sep = ",")}}
           else if(type=="P"){
              if (cell_type=="P" && tpm_source>2^5 && 
                  all(tpm_source>5*tpm_vec[!cell_vec %in% c(well,"A1","1A1","2A1","3A1","4A1")]))
                  {tcr$source[tcrind]<-paste(tcr$source[tcrind],well,sep = ",")}}
            }else{#cells in multiple columns(rows)
                   cross_pos<-sort(apply(expand.grid(common_rows,common_columns), 1, paste, collapse = "", sep = ""))
                   cross_cells<-cell_vec[cell_vec%in%cross_pos]
                   source_cells<-cross_cells[as.numeric(tpm_vec[match(cross_cells,cell_vec)])>2^5]#512
                   if (length(source_cells)>0){
                     scells<-source_cells
                     source_tpms<-tpm_vec[match(source_cells,cell_vec)]
                     if(exists("cells_in_line_sources")){rm(cells_in_line_sources)}
                     for (i in 1:length(source_cells)){
                         ifelse(length(source_cells)==1,scell<-source_cells,scell<-source_cells[i])
                         ifelse(length(source_cells)==1,source_tpm<-as.numeric(source_tpms),source_tpm<-as.numeric(source_tpms[i]))
                         source_row<-gsub("([0-9]+)$", "",scell)
                         source_column<-gsub("[^\\d]+", "",gsub("^([0-9]+)","",scell),perl=TRUE)
                         is_cells_in_line<-(rows==source_row|columns==source_column)
                         cells_in_line<-cell_vec[is_cells_in_line]
                         ifelse(exists("cells_in_line_sources"),cells_in_line_sources<-c(cells_in_line_sources,cells_in_line),cells_in_line_sources<-cells_in_line)
                         tpm_in_line<-as.numeric(tpm_vec[match(cells_in_line,cell_vec)])
                      if (!all(source_tpm>5*tpm_in_line[!cells_in_line %in% c(cells_in_line_sources,"A1","1A1","2A1","3A1","4A1")]) | source_tpm<2^5){
                             scells<-scells[scells != scell]#remove scell from source_cells
                             cells_in_line_sources[!cells_in_line_sources%in%cells_in_line]
                           }
                       }
                     ifelse(type=="T",source_cells<-scells[celltype$cell.type[match(scells,celltype$well)]=="T"],
                           source_cells<-scells[celltype$cell.type[match(scells,celltype$well)]=="P"])        
                     }
                 }
             } 
           tcr$source[tcrind]<-paste(source_cells,collapse = ",")
         }
       return(tcr)
     }

#add a list of source_cells for each V_genes in v1 
append_subset_v<-function(tcr_source,v){
  v_gene_vec<-str_split_fixed(tcr_source$TCR,"_",2)[,1]#[nchar(tcr1_source$source)!=0]
  tcr_source$V_GENE<-v_gene_vec
  sv<-v[str_split_fixed(v$V_GENE,"_",2)[,1] %in% v_gene_vec,]
  sv$V_GENE<-str_split_fixed(sv$V_GENE,"_",2)[,1]
  sv<-merge(sv,aggregate(source~V_GENE,data = tcr_source,FUN = c),by="V_GENE")
  sv$source<-lapply(sv$source,function(x){paste(x,collapse = ",")})
  return(sv)
}

get_potential_source1<-function(tcr_source,sv)
{potential_source<-rep(NA,dim(sv)[1])
v_vec<-unique(str_split_fixed(tcr_source$TCR,"_",2)[,1])
v_vec[1]<-"TRAV38-2_DV8"########only plate1
for (i in 1:length(v_vec)){
  acell<-sv$cells[which(sv$V_GENE==v_vec[i])]
  acell<-unlist(str_split(acell,","))
  is_t<-celltype$cell.type[match(acell,celltype$well)]=="T"
  atpm<-sv$TPMs[which(sv$V_GENE==v_vec[i])]
  atpm<-as.numeric(unlist(str_split(atpm,",")))
  is_hightpm<-atpm>2^5
  rows<-gsub("([0-9]+)$", "", acell)
  columns<-gsub("[^\\d]+", "", acell, perl=TRUE)
  is_spread<-unlist(lapply(rows,function(x) length(x[rows==x])>1))&unlist(lapply(columns,function(x) length(x[columns==x])>1))
  potential_source[i]<-paste(acell[is_t&is_hightpm&is_spread],collapse = ",")
}
return(potential_source)
}

get_potential_source23<-function(tcr_source,sv)
{potential_source<-rep(NA,dim(sv)[1])
v_vec<-unique(str_split_fixed(tcr_source$TCR,"_",2)[,1])
for (i in 1:length(v_vec)){
  acell<-sv$cells[which(sv$V_GENE==v_vec[i])]
  acell<-unlist(str_split(acell,","))
#  acell[!grepl("[*]",acell)]<-paste("1",acell[!grepl("[*]",acell)],sep = "")
#  acell<-gsub("[*]","2",acell)
  is_t<-celltype$cell.type[match(acell,celltype$well)]=="T"
  atpm<-sv$TPMs[which(sv$V_GENE==v_vec[i])]
  atpm<-as.numeric(unlist(str_split(atpm,",")))
  is_hightpm<-atpm>2^5
  rows<-gsub("([0-9]+)$", "", acell)
  columns<-gsub("[^\\d]+", "", acell, perl=TRUE)
  is_spread<-unlist(lapply(rows,function(x) length(x[rows==x])>1))&unlist(lapply(columns,function(x) length(x[columns==x])>1))
  potential_source[i]<-paste(acell[is_t&is_hightpm&is_spread],collapse = ",")
}
return(potential_source)
}

get_potential_source45<-function(tcr_source,sv)
{potential_source<-rep(NA,dim(sv)[1])
v_vec<-unique(str_split_fixed(tcr_source$TCR,"_",2)[,1])
for (i in 1:length(v_vec)){
  acell<-sv$cells[which(sv$V_GENE==v_vec[i])]
  acell<-unlist(str_split(acell,","))
#  acell[!grepl("[*]",acell)]<-paste("3",acell[!grepl("[*]",acell)],sep = "")
#  acell<-gsub("[*]","4",acell)
  is_t<-celltype$cell.type[match(acell,celltype$well)]=="T"
  atpm<-sv$TPMs[which(sv$V_GENE==v_vec[i])]
  atpm<-as.numeric(unlist(str_split(atpm,",")))
  is_hightpm<-atpm>2^5
  rows<-gsub("([0-9]+)$", "", acell)
  columns<-gsub("[^\\d]+", "", acell, perl=TRUE)
  is_spread<-unlist(lapply(rows,function(x) length(x[rows==x])>1))&unlist(lapply(columns,function(x) length(x[columns==x])>1))
  potential_source[i]<-paste(acell[is_t&is_hightpm&is_spread],collapse = ",")
}
return(potential_source)
}

#to check if a cell is source for v gene
check_source<-function(gene,cell_list,use="vsourcelist")#use=tcrsourcelist,vsourcelist,bcrsourcelist,sbcrsourcelist,ntbcrsourcelist
{cell_list<-unlist(strsplit(cell_list,split=","))
n=length(cell_list)
ifs<-rep(NA,n)
for (i in 1:n){
  if (!is.na(cell_list[1])){
    scells<-NA
    if(use=="vsourcelist"){ifs[i]<-gene%in%sourcelist[cell_list[i],]}
    else if(use=="tcrsourcelist"){
      ifelse(grepl("^TRA",gene),ab<-2,ab<-3)
      ifs[i]<-gene%in%tcrsourcelist[cell_list[i],]|is.na(tcrsourcelist[cell_list[i],ab])
      }
    else if(use=="bcrsourcelist"){ifs[i]<-gene%in%bcrsourcelist[cell_list[i],]}
    else if(use=="sbcrsourcelist"){ifs[i]<-gene%in%sbcrsourcelist[cell_list[i],]}
    else if(use=="ntbcrsourcelist"){ifs[i]<-gene%in%ntbcrsourcelist[cell_list[i],]}
    scells<-paste(cell_list[ifs],collapse = ",") 
  }
}
return (scells)
}

#filter non_source
filter_non_source<-function(gene,cell_list,use)#use=tcrsourcelist,vsourcelist,bcrsourcelist,sbcrsourcelist,ntbcrsourcelist
{#cell_list<-unlist(strsplit(cell_list,split=","))
n=length(cell_list)
ifs<-rep(NA,n)
for (i in 1:n){
  if (!is.na(cell_list[1])){
    scells<-NA
    if(use=="v"){ifs[i]<-!gene%in%sourcelist[cell_list[i],]}
    else if(use=="tcr"){
      ifelse(grepl("^TRA",gene),ab<-2,ab<-3)
      ifs[i]<-!gene%in%tcrsourcelist[cell_list[i],]
    }
    else if(use=="bcr"){!ifs[i]<-gene%in%bcrsourcelist[cell_list[i],]}
    else if(use=="aa"){!ifs[i]<-gene%in%sbcrsourcelist[cell_list[i],]}
    else if(use=="nt"){!ifs[i]<-gene%in%ntbcrsourcelist[cell_list[i],]}
    scells<-paste(cell_list[ifs],collapse = ",") 
  }
}
return (scells)
}

#get expression array from v gene$filtered_source
get_array<-function(sv,n_plate="T"){
  source_cells<-unlist(str_split(paste(sv$filted_source[sv$filted_source!=""],collapse = ","),","))
  ifelse(n_plate=="T",vexpr<-array(0,dim = c(16,12,length(source_cells))),vexpr<-array(0,dim = c(8,12,length(source_cells))))
  k=1
  for (i in 1:length(sv$filted_source)){#loop through all v genes in sv
    if(!is.na(sv$filted_source[i])&&sv$filted_source[i]!=""){
      cell_vec<-unlist(str_split(sv$cells[[i]],","))
      tpm_vec<-unlist(str_split(sv$TPMs[[i]],","))
      source_cell_vec<-unlist(str_split(sv$filted_source[[i]],","))
      source_cell_vec<-gsub("^[1,3]","",source_cell_vec)
      source_cell_vec<-gsub("^[2,4]","*",source_cell_vec)
      source_tpm_vec<-tpm_vec[cell_vec%in%source_cell_vec]
      rows<-gsub("([0-9]+)$", "",cell_vec)
      columns<-gsub("[^\\d]+", "",cell_vec,perl=TRUE)
      if (length(source_cell_vec)==1){#only one source cell
        id<-which(cell_vec==source_cell_vec)
        tpm_source<-tpm_vec[id]
        l<-rows%in%(gsub("([0-9]+)$", "",source_cell_vec))|columns%in%gsub("[^\\d]+", "",source_cell_vec,perl=TRUE)
        cells_line<-cell_vec[l]#cells in the same row|column with the source
        allrowcells<-paste(gsub("([0-9]+)$", "",source_cell_vec),1:12,sep="")
        ifelse(n_plate=="F",
               allcolumncells<-paste(c("A","B","C","D","E","F","G","H"),gsub("[^\\d]+", "",source_cell_vec,perl=TRUE),sep=""),
               allcolumncells<-paste(c("A","B","C","D","E","F","G","H","*A","*B","*C","*D","*E","*F","*G","*H"),gsub("[^\\d]+", "",source_cell_vec,perl=TRUE),sep=""))
        tpm_vec_line<-c(tpm_vec[l],rep(0,length(setdiff(union(allcolumncells,allrowcells),cells_line))))
        cells_line<-c(cells_line,setdiff(union(allcolumncells,allrowcells),cells_line))
        ratios<-as.numeric(tpm_vec_line)/as.numeric(tpm_source)
        ifelse(n_plate=="F",vexpr[,,k]<-expr_matrix(cells_line,ratios),vexpr[,,k]<-expr_2matrix(cells_line,ratios))
        k=k+1
      }
      else{#multiple sources
        for (j in 1:length(source_cell_vec)){
          tpm_vec<-unlist(str_split(sv$TPMs[[i]],","))
          id<-which(cell_vec==source_cell_vec[j])
          tpm_source<-tpm_vec[id]
          is_cells_in_row<-rows==(gsub("([0-9]+)$", "",source_cell_vec[j]))
          is_cells_in_col<-columns==gsub("[^\\d]+", "",source_cell_vec[j],perl=TRUE)
          is_cells_in_line<-is_cells_in_row|is_cells_in_col
          cell_in_line<-cell_vec[is_cells_in_line]#cells in the same row|column with the source:source_cell_vec[j]
          tpm_in_line<-tpm_vec[is_cells_in_line]
          for (t in cell_in_line) {
            if (t!=source_cell_vec[j]){
              t_row<-gsub("([0-9]+)$", "",t)
              t_col<-regmatches(t, regexpr("\\d$", t))
              is_sources_in_row<-gsub("([0-9]+)$", "",source_cell_vec)==t_row
              is_sources_in_col<-gsub("[^\\d]+", "",source_cell_vec, perl=TRUE)==t_col
              if ((sum(is_sources_in_row)+sum(is_sources_in_col))>1)
              {
                perc_tpm<-as.numeric(tpm_source)/sum(as.numeric(source_tpm_vec[is_sources_in_row|is_sources_in_col]))
                tpm_in_line[match(t,cell_in_line)]<-as.numeric(tpm_in_line[match(t,cell_in_line)])*perc_tpm
              }
            }
          }
          allrowcells<-paste(gsub("([0-9]+)$", "",source_cell_vec[j]),1:12,sep="")
          ifelse(n_plate=="F",
                 allcolumncells<-paste(c("A","B","C","D","E","F","G","H"),gsub("[^\\d]+", "",source_cell_vec[j],perl=TRUE),sep=""),
                 allcolumncells<-paste(c("A","B","C","D","E","F","G","H","*A","*B","*C","*D","*E","*F","*G","*H"),gsub("[^\\d]+", "",source_cell_vec[j],perl=TRUE),sep=""))
          app_tpm_in_line<-c(tpm_in_line,rep(0,length(setdiff(union(allcolumncells,allrowcells),cell_in_line))))
          app_cell_in_line<-c(cell_in_line,setdiff(union(allrowcells,allcolumncells),cell_in_line))
          ratios<-as.numeric(app_tpm_in_line)/as.numeric(tpm_source)
          ratios[which(app_cell_in_line==source_cell_vec[j])]<-0
          ifelse(n_plate=="F",vexpr[,,k]<-expr_matrix(app_cell_in_line,ratios),vexpr[,,k]<-expr_2matrix(app_cell_in_line,ratios))
          k=k+1
        }
      }
    }
  }
  return(vexpr)
}

#get expression array from tcr_source$source
get_array_tcr<-function(tcr_source,n_plate="T")
{source_cells<-unlist(str_split(paste(tcr_source$source[tcr_source$source!=""],collapse = ","),","))
ifelse(n_plate=="T",tcr_array<-array(0,dim = c(16,12,length(source_cells))),tcr_array<-array(0,dim = c(8,12,length(source_cells))))
k=1
for (i in 1:length(tcr_source$source)){
  cell_vec<-unlist(str_split(tcr_source$cells[[i]],","))
  tpm_vec<-unlist(str_split(tcr_source$TPMs[[i]],","))
  source_cell_vec<-unlist(str_split(tcr_source$source[[i]],","))
  source_cell_vec<-gsub("^[1,3]","",source_cell_vec)
  source_cell_vec<-gsub("^[2,4]","*",source_cell_vec)
  source_tpm_vec<-tpm_vec[cell_vec%in%source_cell_vec]
  rows<-gsub("([0-9]+)$", "",cell_vec)
  columns<-gsub("[^\\d]+", "",cell_vec,perl=TRUE)
  if (source_cell_vec!=""&&length(source_cell_vec)==1){#only one source cell
    id<-which(cell_vec==source_cell_vec)
    tpm_source<-tpm_vec[id]
    l<-rows%in%(gsub("([0-9]+)$", "",source_cell_vec))|columns%in%gsub("[^\\d]+", "",source_cell_vec,perl=TRUE)
    cells_line<-cell_vec[l]#cells in the same row|column with the source
    allrowcells<-paste(gsub("([0-9]+)$", "",source_cell_vec),1:12,sep="")
    ifelse(n_plate=="F",
           allcolumncells<-paste(c("A","B","C","D","E","F","G","H"),gsub("[^\\d]+", "",source_cell_vec,perl=TRUE),sep=""),
           allcolumncells<-paste(c("A","B","C","D","E","F","G","H","*A","*B","*C","*D","*E","*F","*G","*H"),gsub("[^\\d]+", "",source_cell_vec,perl=TRUE),sep=""))
    tpm_vec_line<-c(tpm_vec[l],rep(0,length(setdiff(union(allcolumncells,allrowcells),cells_line))))
    cells_line<-c(cells_line,setdiff(union(allcolumncells,allrowcells),cells_line))
    ratios<-as.numeric(tpm_vec_line)/as.numeric(tpm_source)
    ifelse(n_plate=="F",tcr_array[,,k]<-expr_matrix(cells_line,ratios),tcr_array[,,k]<-expr_2matrix(cells_line,ratios))
    k=k+1
  }else if(source_cell_vec!=""&&length(source_cell_vec)>1){#multiple source cell
    for (j in 1:length(source_cell_vec)){
      tpm_vec<-unlist(str_split(tcr_source$TPMs[[i]],","))
      id<-which(cell_vec==source_cell_vec[j])
      tpm_source<-tpm_vec[id]
      is_cells_in_row<-rows==(gsub("([0-9]+)$", "",source_cell_vec[j]))
      is_cells_in_col<-columns==gsub("[^\\d]+", "",source_cell_vec[j],perl=TRUE)
      is_cells_in_line<-is_cells_in_row|is_cells_in_col
      cell_in_line<-cell_vec[is_cells_in_line]#cells in the same row|column with the source:source_cell_vec[j]
      tpm_in_line<-tpm_vec[is_cells_in_line]
      for (t in cell_in_line) {
        if (t!=source_cell_vec[j]){
          t_row<-gsub("([0-9]+)$", "",t)
          t_col<-regmatches(t, regexpr("\\d$", t))
          is_sources_in_row<-gsub("([0-9]+)$", "",source_cell_vec)==t_row
          is_sources_in_col<-gsub("[^\\d]+", "",source_cell_vec)==t_col
          if ((sum(is_sources_in_row)+sum(is_sources_in_col))>1)
          {
            perc_tpm<-as.numeric(tpm_source)/sum(as.numeric(source_tpm_vec[is_sources_in_row|is_sources_in_col]))
            tpm_in_line[match(t,cell_in_line)]<-as.numeric(tpm_in_line[match(t,cell_in_line)])*perc_tpm
          }
        }
      }
      allrowcells<-paste(gsub("([0-9]+)$", "",source_cell_vec[j]),1:12,sep="")
      ifelse(n_plate=="F",
             allcolumncells<-paste(c("A","B","C","D","E","F","G","H"),gsub("[^\\d]+", "",source_cell_vec[j],perl=TRUE),sep=""),
             allcolumncells<-paste(c("A","B","C","D","E","F","G","H","*A","*B","*C","*D","*E","*F","*G","*H"),gsub("[^\\d]+", "",source_cell_vec[j],perl=TRUE),sep=""))
      app_tpm_in_line<-c(tpm_in_line,rep(0,length(setdiff(union(allcolumncells,allrowcells),cell_in_line))))
      app_cell_in_line<-c(cell_in_line,setdiff(union(allrowcells,allcolumncells),cell_in_line))
      ratios<-as.numeric(app_tpm_in_line)/as.numeric(tpm_source)
      ratios[which(app_cell_in_line==source_cell_vec[j])]<-0
      ifelse(n_plate=="F",tcr_array[,,k]<-expr_matrix(app_cell_in_line,ratios),tcr_array[,,k]<-expr_2matrix(app_cell_in_line,ratios))
      k=k+1
    }
  }
}
return(tcr_array)
}

collect_rows<-function(array,array_use=""){
  ifelse(array_use=="",n<-dim(array)[3],n<-length(array_use))
  target_rows<-matrix(NA,nrow = n,ncol = 12)
  row_freq<-rep(0,dim(array)[1])
  row_ids<-NA
  for (i in 1:n){
    row_idx<-which(rowSums(array[,,i])>=0)
    row_freq[row_idx]=row_freq[row_idx]+1
    target_rows[i,]<-array[,,i][row_idx,]
    ifelse(is.na(row_ids),row_ids<-row_idx,row_ids<-c(row_ids,row_idx))
  }
  return(list(target_rows,row_freq,row_ids))
}

collect_cols<-function(array,array_use=""){
  ifelse(array_use=="",n<-dim(array)[3],n<-length(array_use))
  k<-dim(array)[1]
  target_columns<-matrix(NA,nrow = k,ncol = n)
  col_freq<-rep(0,12)
  col_ids<-NA
  for (i in 1:n){
    col_idx<-which(colSums(array[,,i])>=0)
    col_freq[col_idx]=col_freq[col_idx]+1
    target_columns[,i]<-array[,,i][,col_idx]
    ifelse(is.na(col_ids),col_ids<-col_idx,col_ids<-c(col_ids,col_idx))
  }
  return(list(target_columns,col_freq,col_ids))
}

#correct the sequences such that allow one nt mismatch if BCRs with the same V_J
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
# 
#   bcr$ReceptorID<-paste(bcr$V,bcr$N,bcr$J,sep = "_")
#   bcr$Freq<-as.numeric(bcr$Freq)
#   # v1 <- aggregate(Freq~ReceptorID,data=bcr,FUN=sum)
#   # v2 <- aggregate(cells~ReceptorID,data=bcr,FUN=c)
#   # v3 <- aggregate(TPMs~ReceptorID,data=bcr,FUN=c)
#   # bcr <- merge(v1, v2, by="ReceptorID")
#   # bcr <- merge(bcr, v3, by="ReceptorID")
#   v1 <- aggregate(Freq~full_sequence,data=bcr,FUN=sum)
#   v2 <- aggregate(cells~full_sequence,data=bcr,FUN=c)
#   v3 <- aggregate(TPMs~full_sequence,data=bcr,FUN=c)
#   bcr <- merge(v1, v2, by="full_sequence")
#   bcr <- merge(bcr, v3, by="full_sequence")
#   for (i in 1:dim(bcr)[1]){
#     if (length(unlist(bcr$cells[i]))>1){
#       cells<-unlist(str_split(paste(unlist(bcr$cells[i]),collapse =","),","))
#       if(length(unique(cells))<length(cells)){
#         tpms<-as.numeric(unlist(str_split(paste(unlist(bcr$TPMs[i]),collapse =","),",")))
#         df<-data.frame(cells,tpms)
#         adf<-aggregate(tpms~cells,data = df, FUN=sum)
#         bcr$cells[i]<-paste((adf$cells),collapse =",")
#         bcr$TPMs[i]<-paste((adf$tpms),collapse =",")
#       }else{
#         bcr$cells[i]<-paste(unlist(bcr$cells[i]),collapse =",")
#         bcr$TPMs[i]<-paste(unlist(bcr$TPMs[i]),collapse =",")
#       }
#     }
#     else
#     { bcr$cells[i]<-unlist(bcr$cells[i])
#       bcr$TPMs[i]<-unlist(bcr$TPMs[i])
#     }
#   }
#   bcr$cells<-unlist(bcr$cells)
#   bcr$TPMs<-unlist(bcr$TPMs)
#   return(bcr)
# }

merge_seq<-function(bracer,bracer_seq){
  bracer<-bracer[order(bracer[,1], bracer[,2]), ]
  bracer_seq<-bracer_seq[order(bracer_seq[,1], bracer_seq[,2]), ]
  bracer<-cbind(bracer,bracer_seq$Sequence)
  colnames(bracer)[4]<-"sequence"
  return(bracer)
}
