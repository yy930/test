#n_plate="T"/"F";plate=1/23/45;use="v"/"tcr"/"bcr"/"aa/nt";
get_array<-function(v,n_plate="T",plate=23,use="v"){
  #ifelse(use%in%c("v","bcr","aa","nt"),source_cells_list<-str_split(paste(v$filted_source),","),source_cells_list<-str_split(paste(v$source),","))
  source_cells_list<-str_split(paste(v$filted_source),",")
  is_short_source<-unlist(lapply(source_cells_list, function(x) length(x)))<3
  #ifelse(use%in%c("v","bcr","aa","nt"),is_non_emp<-v$filted_source!="",is_non_emp<-v$source!="")
  is_non_emp<-v$filted_source!=""
  sourceid_to_use<-which(is_short_source&is_non_emp)
  #source_cells<-unlist(str_split(paste(v$filted_source[v$filted_source!=""],collapse = ","),","))
  source_cells<-unlist(source_cells_list[sourceid_to_use])
  ifelse(n_plate=="T",vexpr<-array(0,dim = c(16,12,length(source_cells))),vexpr<-array(0,dim = c(8,12,length(source_cells))))
  k=1
  v_tcr_vec<-rep(NA,length(source_cells))
  s_vec<-rep(NA,length(source_cells))
  for (i in sourceid_to_use){#loop through all v genes in v
      if (use=="v") {v_tcr<-v$V_GENE[i]}
      else if (use=="tcr") {v_tcr<-v$TCR[i]}
      else if (use=="bcr") {v_tcr<-v$ReceptorID[i]}
      else if (use=="aa") {v_tcr<-v$CDR3aa[i]}
      else if (use=="nt") {v_tcr<-v$CDR3nt[i]}
      cell_vec<-unlist(str_split(v$cells[[i]],","))
      tpm_vec<-unlist(str_split(v$TPMs[[i]],","))
#      ifelse(use%in%c("v","bcr","aa","nt"),source_cell_vec<-unlist(str_split(v$filted_source[[i]],",")),source_cell_vec<-unlist(str_split(v$source[[i]],",")))
      source_cell_vec<-unlist(str_split(v$filted_source[[i]],","))
#      source_cell_vec<-gsub("^[1,3]","",source_cell_vec)
#      source_cell_vec<-gsub("^[2,4]","*",source_cell_vec)
      source_tpm_vec<-tpm_vec[match(source_cell_vec,cell_vec)]
      rows<-gsub("([0-9]+)$", "",cell_vec)
      columns<-gsub("[^\\d]+", "",gsub("^([0-9]+)","",cell_vec),perl=TRUE)
      if (length(source_cell_vec)==1){#only one source cell
        source_col<-gsub("[^\\d]+", "",gsub("^([0-9]+)","",source_cell_vec),perl=TRUE)
        l<-rows%in%gsub("([0-9]+)$", "",source_cell_vec)|columns%in%source_col
        if (sum(l==FALSE)<6 && length(l)>length(source_cell_vec)+3 && sum(l)>(length(l)-sum(l)+1))
          {
          id<-which(cell_vec==source_cell_vec)
          tpm_source<-tpm_vec[id]
          cells_line<-cell_vec[l]#cells in the same row|column with the source
          tpms_line<-tpm_vec[l]#tpms in the same row|column with the source
          if (all(as.numeric(tpm_source)>as.numeric(tpms_line[-which(cells_line==source_cell_vec)])*5)){
                allrowcells<-paste(gsub("([0-9]+)$", "",source_cell_vec),1:12,sep="")
          ifelse(n_plate=="F",
               allcolumncells<-paste(c("A","B","C","D","E","F","G","H"),source_col,sep=""),
                 ifelse(plate==23,allcolumncells<-paste(c("1A","1B","1C","1D","1E","1F","1G","1H","2A","2B","2C","2D","2E","2F","2G","2H"),source_col,sep=""),
                 allcolumncells<-paste(c("3A","3B","3C","3D","3E","3F","3G","3H","4A","4B","4C","4D","4E","4F","4G","4H"),source_col,sep="")))
          tpm_vec_line<-c(tpm_vec[l],rep(0,length(setdiff(union(allcolumncells,allrowcells),cells_line))))
          cells_line<-c(cells_line,setdiff(union(allcolumncells,allrowcells),cells_line))
          ratios<-as.numeric(tpm_vec_line)/as.numeric(tpm_source)
          ifelse(n_plate=="F",vexpr[,,k]<-expr_matrix(cells_line,ratios,plate),vexpr[,,k]<-expr_2matrix(cells_line,ratios,plate))
          s_vec[k]<-source_cell_vec
          v_tcr_vec[k]<-v_tcr
          k=k+1      
          }
        }
      }
      else{#multiple sources
        for (j in 1:length(source_cell_vec)){
          #tpm_vec<-unlist(str_split(v$TPMs[[i]],","))
          id<-which(cell_vec%in%source_cell_vec)#[j])
          tpm_source<-tpm_vec[id][j]
          is_cells_in_rows<-rows%in%(gsub("([0-9]+)$", "",source_cell_vec))#[j]))
          is_cells_in_row<-rows%in%(gsub("([0-9]+)$", "",source_cell_vec[j]))
          source_cols<-gsub("[^\\d]+", "",gsub("^([0-9]+)","",source_cell_vec),perl=TRUE)
          source_col<-gsub("[^\\d]+", "",gsub("^([0-9]+)","",source_cell_vec[j]),perl=TRUE)
          is_cells_in_cols<-columns%in%source_cols
          is_cells_in_col<-columns%in%source_col
          is_cells_in_lines<-is_cells_in_rows|is_cells_in_cols
          is_cells_in_line<-is_cells_in_row|is_cells_in_col
          if (sum(is_cells_in_lines==FALSE)<6 && length(is_cells_in_lines)>length(source_cell_vec)+3 && 
              sum(is_cells_in_lines)>(length(is_cells_in_lines)-sum(is_cells_in_lines))){
          cell_in_line<-cell_vec[is_cells_in_line]#cells in the same row|column with the source:source_cell_vec[j]
          tpm_in_line<-tpm_vec[is_cells_in_line]
          for (t in cell_in_line) {
            if (t!=source_cell_vec[j]){
              t_row<-gsub("([0-9]+)$", "",t)
              t_col<-gsub("[^\\d]+", "",gsub("^([0-9]+)","",t),perl=TRUE)
              is_sources_in_row<-gsub("([0-9]+)$", "",source_cell_vec)==t_row
              is_sources_in_col<-source_cols==t_col
              perc_tpm<-as.numeric(tpm_source)/sum(as.numeric(source_tpm_vec[is_sources_in_row|is_sources_in_col]))
              if (sum(is_sources_in_row)>1)#2 sources in same row
              {if (t %in% cell_vec[is_cells_in_row])
                #tpm_in_line[match(t,cell_vec[is_cells_in_row])]<-as.numeric(tpm_in_line[match(t,cell_vec[is_cells_in_row])])*perc_tpm
                {tpm_in_line[match(t,cell_in_line)]<-as.numeric(tpm_in_line[match(t,cell_in_line)])*perc_tpm}
              }
              else if (sum(is_sources_in_col)>1)#2 sources in same column
              {if (t %in% cell_vec[is_cells_in_col])
                {tpm_in_line[match(t,cell_in_line)]<-as.numeric(tpm_in_line[match(t,cell_in_line)])*perc_tpm}
              }
              else#2 sources in different column/row
              {
                if(t%in%cell_vec[is_cells_in_rows&is_cells_in_cols])
                  {tpm_in_line[match(t,cell_in_line)]<-as.numeric(tpm_in_line[match(t,cell_in_line)])*perc_tpm}
              }
            }
          }
          if (all(as.numeric(tpm_source)>as.numeric(tpm_in_line[-which(cell_in_line==source_cell_vec[j])])*5))
          {allrowcells<-paste(gsub("([0-9]+)$", "",source_cell_vec[j]),1:12,sep="")###############
          ifelse(n_plate=="F",
                 allcolumncells<-paste(c("A","B","C","D","E","F","G","H"),source_col,sep=""),
                 ifelse(plate==23,allcolumncells<-paste(c("1A","1B","1C","1D","1E","1F","1G","1H","2A","2B","2C","2D","2E","2F","2G","2H"),source_col,sep=""),
                        allcolumncells<-paste(c("3A","3B","3C","3D","3E","3F","3G","3H","4A","4B","4C","4D","4E","4F","4G","4H"),source_col,sep="")))
          app_tpm_in_line<-c(tpm_in_line,rep(0,length(setdiff(union(allcolumncells,allrowcells),cell_in_line))))
          app_cell_in_line<-c(cell_in_line,setdiff(union(allrowcells,allcolumncells),cell_in_line))
          ratios<-as.numeric(app_tpm_in_line)/as.numeric(tpm_source)
          ratios[which(app_cell_in_line==source_cell_vec[j])]<-0
          ifelse(n_plate=="F",vexpr[,,k]<-expr_matrix(app_cell_in_line,ratios,plate),vexpr[,,k]<-expr_2matrix(app_cell_in_line,ratios,plate))
          s_vec[k]<-source_cell_vec[j]
          v_tcr_vec[k]<-v_tcr
          k=k+1}
          }
        }
      }
    }
  return(list(vexpr[,,1:(k-1)],v_tcr_vec[1:(k-1)],s_vec[1:(k-1)]))
}