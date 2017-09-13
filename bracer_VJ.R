#aggretate V_J in bracer
vjplatelist<-list()
for (i in 1:length(bracerlist)){
  vjplatelist[[i]]<-bracerlist[[i]][bracerlist[[i]]$TPM>=1,]
  V_N_J<-str_split_fixed(vjplatelist[[i]]$ReceptorID,"_",3)
  vjplatelist[[i]]$flag<-paste(vjplatelist[[i]]$CellID,V_N_J[,1],V_N_J[,3],sep='|')
  agg_bracer<-aggregate(TPM~flag,data=vjplatelist[[i]],FUN=sum)
  agg_bracer$CellID<-str_split_fixed(agg_bracer$flag, "[|]", 3)[,1]
  agg_bracer$V_J<-str_split_fixed(agg_bracer$flag, "[|]", 2)[,2]
  vjplatelist[[i]]<-agg_bracer[ , -which(names(agg_bracer)=="flag" )]
}

#freq for each vj in plate
vj1<-as.data.frame(table(vjplatelist[[1]]$V_J))
vj2<-as.data.frame(table(vjplatelist[[2]]$V_J))
vj3<-as.data.frame(table(vjplatelist[[3]]$V_J))
vj4<-as.data.frame(table(vjplatelist[[4]]$V_J))
vj5<-as.data.frame(table(vjplatelist[[5]]$V_J))
colnames(vj1)[1]=colnames(vj2)[1]=colnames(vj3)[1]=colnames(vj4)[1]=colnames(vj5)[1]="V_J"
vjlist<-list(vj1,vj2,vj3,vj4,vj5)

#get cells in which each vj has been detected
for (plateID in 1:5){
  i=0
  for (vj in vjlist[[plateID]]$V_J){
    i=i+1
    ind_cell <- vjplatelist[[plateID]]$V_J==vj
    cells <- vjplatelist[[plateID]]$CellID[ind_cell]
    cells<-paste0(cells,collapse = ",")
    vjlist[[plateID]]$cells[i] <- cells
    tpms <- paste(vjplatelist[[plateID]]$TPM[ind_cell],collapse = ",")
    vjlist[[plateID]]$TPMs[i] <- tpms
  }
}

vj1 <- data.frame(lapply(vjlist[[1]], as.character), stringsAsFactors=FALSE)
vj2 <- data.frame(lapply(vjlist[[2]], as.character), stringsAsFactors=FALSE)
vj3 <- data.frame(lapply(vjlist[[3]], as.character), stringsAsFactors=FALSE)
vj4 <- data.frame(lapply(vjlist[[4]], as.character), stringsAsFactors=FALSE)
vj5 <- data.frame(lapply(vjlist[[5]], as.character), stringsAsFactors=FALSE)

vj23<-merge(vj2,vj3,by=c("V_J","V_J"),all = FALSE,suffixes=c("_plate2","_plate3"))
vj23$sumfreq<-as.numeric(vj23$Freq_plate2)+as.numeric(vj23$Freq_plate3)

vj45<-merge(vj4,vj5,by=c("V_J","V_J"),all = FALSE,suffixes=c("_plate4","_plate5"))
vj45$sumfreq<-as.numeric(vj45$Freq_plate4)+as.numeric(vj45$Freq_plate5)

