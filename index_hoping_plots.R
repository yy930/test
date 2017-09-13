##############PLOT_V##############PLOT_V##############PLOT_V##############PLOT_V####
#plot for all v genes in plate1
vplotlist1<-list()
#length(which(as.numeric(v1$Freq)>3))
#n=0
for (i in 1:dim(v1)[1]){
#for (i in which(as.numeric(v1$Freq)>3)){
  receptor=v1$V_GENE[i]
  plateID=1
  geneID=which(vlist[[1]]$V_GENE==receptor)
  mat<-expr_mat(plateID,geneID,use = "v")
#  n=n+1
  vplotlist1[[i]]<-visheatmap1(mat,.title=paste(receptor," in Plate",plateID),.labs = c("", ""),.legend = "log2(TPM)")}
pdf("/Users/yingy_adm/Desktop/agg_v_heatmap1_celltype.pdf",height=115,width=8)
multiplot(vplotlist1[[1]], vplotlist1[[2]], vplotlist1[[3]], vplotlist1[[4]],vplotlist1[[5]],
          vplotlist1[[6]],vplotlist1[[7]],vplotlist1[[8]],vplotlist1[[9]],vplotlist1[[10]],
          vplotlist1[[11]],vplotlist1[[12]], vplotlist1[[13]], vplotlist1[[14]],vplotlist1[[15]],
          vplotlist1[[16]],vplotlist1[[17]],vplotlist1[[18]],vplotlist1[[19]],vplotlist1[[20]])
#          vplotlist1[[21]], vplotlist1[[22]], vplotlist1[[23]])
dev.off()

vplotlist23<-list()
for (i in 1:55){########!
  receptor=v23$V_GENE[i]
  plateID1=2
  geneID=which(v2$V_GENE==receptor)
  mat1<-expr_mat(plateID1,geneID,use = "v")
  plateID2=3
  geneID=which(v3$V_GENE==receptor)
  mat2<-expr_mat(plateID2,geneID,use = "v")
  rownames(mat2)<-paste(rownames(mat2),"*")
  mat<-rbind(mat1,mat2)
  vplotlist23[[i]]<-visheatmap23(mat,.title=paste(receptor," in Plate",plateID1," + Plate",plateID2,"*"),.labs = c("", ""),.legend = "log2(TPM)")}
pdf("/Users/yingy_adm/Desktop/agg_v_heatmap23_celltype.pdf",height=450,width=8)
multiplot(vplotlist23[[1]], vplotlist23[[2]], vplotlist23[[3]], vplotlist23[[4]],vplotlist23[[5]],
          vplotlist23[[6]],vplotlist23[[7]],vplotlist23[[8]],vplotlist23[[9]],vplotlist23[[10]],
          vplotlist23[[11]],vplotlist23[[12]],vplotlist23[[13]], vplotlist23[[14]],vplotlist23[[15]],
          vplotlist23[[16]],vplotlist23[[17]],vplotlist23[[18]],vplotlist23[[19]],vplotlist23[[20]],
          vplotlist23[[21]], vplotlist23[[22]], vplotlist23[[23]], vplotlist23[[24]],vplotlist23[[25]],
          vplotlist23[[26]],vplotlist23[[27]],vplotlist23[[28]],vplotlist23[[29]],vplotlist23[[30]],
          vplotlist23[[31]], vplotlist23[[32]], vplotlist23[[33]], vplotlist23[[34]],vplotlist23[[35]],
          vplotlist23[[36]],vplotlist23[[37]],vplotlist23[[38]],vplotlist23[[39]],vplotlist23[[40]],
          vplotlist23[[41]], vplotlist23[[42]], vplotlist23[[43]], vplotlist23[[44]],vplotlist23[[45]],
          vplotlist23[[46]],vplotlist23[[47]],vplotlist23[[48]],vplotlist23[[49]],vplotlist23[[50]],
          vplotlist23[[51]], vplotlist23[[52]], vplotlist23[[53]], vplotlist23[[54]],vplotlist23[[55]])
dev.off()  

vplotlist45<-list()
for (i in 1:27){########!
  receptor=v45$V_GENE[i]
  plateID1=4
  geneID=which(v4$V_GENE==receptor)
  mat1<-expr_mat(plateID1,geneID,use="v")
  plateID2=5
  geneID=which(v5$V_GENE==receptor)
  mat2<-expr_mat(plateID2,geneID,use="v")
  rownames(mat2)<-paste(rownames(mat2),"*")
  mat<-rbind(mat1,mat2)
  vplotlist45[[i]]<-visheatmap45(mat,.title=paste(receptor," in Plate",plateID1," + Plate",plateID2,"*"),.labs = c("", ""),.legend = "log2(TPM)")}
pdf("/Users/yingy_adm/Desktop/agg_v_heatmap45_celltype.pdf",height=220,width=8)
multiplot(vplotlist45[[1]], vplotlist45[[2]], vplotlist45[[3]], vplotlist45[[4]],vplotlist45[[5]],
          vplotlist45[[6]],vplotlist45[[7]],vplotlist45[[8]],vplotlist45[[9]],vplotlist45[[10]],
          vplotlist45[[11]],vplotlist45[[12]],vplotlist45[[13]], vplotlist45[[14]],vplotlist45[[15]],
          vplotlist45[[16]],vplotlist45[[17]],vplotlist45[[18]],vplotlist45[[19]],vplotlist45[[20]],
          vplotlist45[[21]], vplotlist45[[22]], vplotlist45[[23]], vplotlist45[[24]],vplotlist45[[25]],
          vplotlist45[[26]],vplotlist45[[27]])
dev.off()  

##############PLOT_TCR##############PLOT_TCR##############PLOT_TCR##############PLOT_TCR####
plotlist1<-list()
for (i in 1:dim(tcr1)[1]){
  receptor=tcr1$TCR[i]
  plateID=1
  tcrID=which(tcrlist[[1]]$TCR==receptor)
  mat<-expr_mat(plateID,tcrID,use = "tcr")
  plotlist1[[i]]<-visheatmap1(mat,.title=paste(receptor," in Plate",plateID),.labs = c("", ""),.legend = "log2(TPM)")}
pdf("/Users/yingy_adm/Desktop/tcr_heatmap1_celltype.pdf",height=60,width=8)
multiplot(plotlist1[[1]], plotlist1[[2]], plotlist1[[3]], plotlist1[[4]],plotlist1[[5]],
          plotlist1[[6]],plotlist1[[7]],plotlist1[[8]],plotlist1[[9]],plotlist1[[10]])
dev.off()

plotlist23<-list()
for (i in 1:dim(tcr23)[1]){
  receptor=tcr23$TCR[i]
  plateID1=2
  tcrID=which(tcr2$TCR==receptor)
  mat1<-expr_mat(plateID1,tcrID,use="tcr")
  plateID2=3
  tcrID=which(tcr3$TCR==receptor)
  mat2<-expr_mat(plateID2,tcrID,use="tcr")
  rownames(mat2)<-paste(rownames(mat2),"*")
  mat<-rbind(mat1,mat2)
  plotlist23[[i]]<-visheatmap23(mat,.title=paste(receptor," in Plate",plateID1," + Plate",plateID2,"*"),.labs = c("", ""),.legend = "log2(TPM)")}
pdf("/Users/yingy_adm/Desktop/tcr_heatmap23_celltype.pdf",height=172,width=8)
multiplot(plotlist23[[1]], plotlist23[[2]], plotlist23[[3]], plotlist23[[4]],plotlist23[[5]],
          plotlist23[[6]],plotlist23[[7]],plotlist23[[8]],plotlist23[[9]],plotlist23[[10]],
          plotlist23[[11]],plotlist23[[12]],plotlist23[[13]], plotlist23[[14]],plotlist23[[15]],
          plotlist23[[16]],plotlist23[[17]],plotlist23[[18]],plotlist23[[19]])
dev.off()  

plotlist45<-list()
for (i in 1:dim(tcr45)[1]){
  receptor=tcr45$TCR[i]
  plateID1=4
  tcrID=which(tcr4$TCR==receptor)
  mat1<-expr_mat(plateID1,tcrID,use="tcr")
  plateID2=5
  tcrID=which(tcr5$TCR==receptor)
  mat2<-expr_mat(plateID2,tcrID,use="tcr")
  rownames(mat2)<-paste(rownames(mat2),"*")
  mat<-rbind(mat1,mat2)
  plotlist45[[i]]<-visheatmap45(mat,.title=paste(receptor," in Plate",plateID1," + Plate",plateID2,"*"),.labs = c("", ""),.legend = "log2(TPM)")}
pdf("/Users/yingy_adm/Desktop/tcr_heatmap45_celltype.pdf",height=140,width=8)
multiplot(plotlist45[[1]], plotlist45[[2]], plotlist45[[3]], plotlist45[[4]],plotlist45[[5]],plotlist45[[6]],plotlist45[[7]],plotlist45[[8]],
          plotlist45[[9]],plotlist45[[10]],plotlist45[[11]], plotlist45[[12]], plotlist45[[13]], plotlist45[[14]],plotlist45[[15]])
dev.off()  

columnplot<-list()
columnplot[[1]]<-barplot(apply(array_sv1, 2, function(x) (sum(x[!is.na(x)])))/apply(array_sv1, 2, function(x) sum(!is.na(x))))
columnplot[[2]]<-barplot(apply(array_sv23, 2, function(x) (sum(x[!is.na(x)])))/apply(array_sv23, 2, function(x) sum(!is.na(x))))
columnplot[[3]]<-barplot(apply(array_sv45, 2, function(x) (sum(x[!is.na(x)])))/apply(array_sv45, 2, function(x) sum(!is.na(x))))
pdf("/Users/yingy_adm/Desktop/columnwise.pdf",height=25,width=8)
multiplot(columnplot[[1]], columnplot[[2]], columnplot[[3]])
dev.off()

####PLOT_V_J_BCR##############PLOT_V_J_BCR##############PLOT_V_J_BCR##############PLOT_V_J_BCR####
#plot for all V_J_BCR genes in plate1
vjplotlist1<-list()
n=1
for (i in which(as.numeric(vj1$Freq)>4)){
  vj=vj1$V_J[i]
  plateID=1
  geneID=which(vjlist[[1]]$V_J==vj)
  mat<-expr_mat(plateID,geneID,use = "vj")
  vjplotlist1[[n]]<-visheatmap1(mat,.title=paste(vj," in Plate",plateID),.labs = c("", ""),.legend = "log2(TPM)")
  n=n+1}
pdf("/Users/yingy_adm/Desktop/vj_plate1.pdf",height=240,width=8)
multiplot(vjplotlist1[[1]],vjplotlist1[[2]],vjplotlist1[[3]],vjplotlist1[[4]],vjplotlist1[[5]],
          vjplotlist1[[6]],vjplotlist1[[7]],vjplotlist1[[8]],vjplotlist1[[9]],vjplotlist1[[10]],
          vjplotlist1[[11]],vjplotlist1[[12]],vjplotlist1[[13]],vjplotlist1[[14]],vjplotlist1[[15]],
          vjplotlist1[[16]],vjplotlist1[[17]],vjplotlist1[[18]],vjplotlist1[[19]],vjplotlist1[[20]],
          vjplotlist1[[21]],vjplotlist1[[22]],vjplotlist1[[23]],vjplotlist1[[24]],vjplotlist1[[25]],
          vjplotlist1[[26]],vjplotlist1[[27]],vjplotlist1[[28]],vjplotlist1[[29]],vjplotlist1[[30]],
          vjplotlist1[[31]],vjplotlist1[[32]],vjplotlist1[[33]],vjplotlist1[[34]],vjplotlist1[[35]],
          vjplotlist1[[36]],vjplotlist1[[37]],vjplotlist1[[38]],vjplotlist1[[39]],vjplotlist1[[40]],
          vjplotlist1[[41]],vjplotlist1[[42]],vjplotlist1[[43]],vjplotlist1[[44]],vjplotlist1[[45]],
          vjplotlist1[[46]])
dev.off()

vjplotlist23<-list()
n=1
for (i in which(vj23$sumfreq>4)){########!
  receptor=vj23$V_J[i]
  plateID1=2
  geneID=which(vj2$V_J==receptor)
  mat1<-expr_mat(plateID1,geneID,use = "vj")
  plateID2=3
  geneID=which(vj3$V_J==receptor)
  mat2<-expr_mat(plateID2,geneID,use = "vj")
  rownames(mat2)<-paste(rownames(mat2),"*")
  mat<-rbind(mat1,mat2)
  vjplotlist23[[n]]<-visheatmap23(mat,.title=paste(receptor," in Plate",plateID1," + Plate",plateID2,"*"),.labs = c("", ""),.legend = "log2(TPM)")
  n=n+1}
pdf("/Users/yingy_adm/Desktop/vj_plate23.pdf",height=570,width=8)
multiplot(vjplotlist23[[1]], vjplotlist23[[2]], vjplotlist23[[3]], vjplotlist23[[4]],vjplotlist23[[5]],
          vjplotlist23[[6]],vjplotlist23[[7]],vjplotlist23[[8]],vjplotlist23[[9]],vjplotlist23[[10]],
          vjplotlist23[[11]],vjplotlist23[[12]],vjplotlist23[[13]], vjplotlist23[[14]],vjplotlist23[[15]],
          vjplotlist23[[16]],vjplotlist23[[17]],vjplotlist23[[18]],vjplotlist23[[19]],vjplotlist23[[20]],
          vjplotlist23[[21]], vjplotlist23[[22]], vjplotlist23[[23]], vjplotlist23[[24]],vjplotlist23[[25]],
          vjplotlist23[[26]],vjplotlist23[[27]],vjplotlist23[[28]],vjplotlist23[[29]],vjplotlist23[[30]],
          vjplotlist23[[31]], vjplotlist23[[32]], vjplotlist23[[33]], vjplotlist23[[34]],vjplotlist23[[35]],
          vjplotlist23[[36]],vjplotlist23[[37]],vjplotlist23[[38]],vjplotlist23[[39]],vjplotlist23[[40]],
          vjplotlist23[[41]], vjplotlist23[[42]], vjplotlist23[[43]], vjplotlist23[[44]],vjplotlist23[[45]],
          vjplotlist23[[46]],vjplotlist23[[47]],vjplotlist23[[48]],vjplotlist23[[49]],vjplotlist23[[50]],
          vjplotlist23[[51]], vjplotlist23[[52]], vjplotlist23[[53]], vjplotlist23[[54]],vjplotlist23[[55]],
          vjplotlist23[[56]],vjplotlist23[[57]])
dev.off()  

vjplotlist45<-list()
n=1
for (i in which(vj45$sumfreq>4)){########!
  receptor=vj45$V_J[i]
  plateID1=4
  geneID=which(vj4$V_J==receptor)
  mat1<-expr_mat(plateID1,geneID,use="vj")
  plateID2=5
  geneID=which(vj5$V_J==receptor)
  mat2<-expr_mat(plateID2,geneID,use="vj")
  rownames(mat2)<-paste(rownames(mat2),"*")
  mat<-rbind(mat1,mat2)
  vjplotlist45[[n]]<-visheatmap45(mat,.title=paste(receptor," in Plate",plateID1," + Plate",plateID2,"*"),.labs = c("", ""),.legend = "log2(TPM)")
  n<-n+1}
pdf("/Users/yingy_adm/Desktop/vj_plate45.pdf",height=570,width=8)
multiplot(vjplotlist45[[1]], vjplotlist45[[2]], vjplotlist45[[3]], vjplotlist45[[4]],vjplotlist45[[5]],
          vjplotlist45[[6]],vjplotlist45[[7]],vjplotlist45[[8]],vjplotlist45[[9]],vjplotlist45[[10]],
          vjplotlist45[[11]],vjplotlist45[[12]],vjplotlist45[[13]], vjplotlist45[[14]],vjplotlist45[[15]],
          vjplotlist45[[16]],vjplotlist45[[17]],vjplotlist45[[18]],vjplotlist45[[19]],vjplotlist45[[20]],
          vjplotlist45[[21]], vjplotlist45[[22]], vjplotlist45[[23]], vjplotlist45[[24]],vjplotlist45[[25]],
          vjplotlist45[[26]],vjplotlist45[[27]],vjplotlist45[[28]],vjplotlist45[[29]],vjplotlist45[[30]],
          vjplotlist45[[31]], vjplotlist45[[32]], vjplotlist45[[33]], vjplotlist45[[34]],vjplotlist45[[35]],
          vjplotlist45[[36]],vjplotlist45[[37]],vjplotlist45[[38]],vjplotlist45[[39]],vjplotlist45[[40]],
          vjplotlist45[[41]], vjplotlist45[[42]], vjplotlist45[[43]], vjplotlist45[[44]],vjplotlist45[[45]],
          vjplotlist45[[46]],vjplotlist45[[47]],vjplotlist45[[48]],vjplotlist45[[49]],vjplotlist45[[50]],
          vjplotlist45[[51]], vjplotlist45[[52]], vjplotlist45[[53]], vjplotlist45[[54]],vjplotlist45[[55]])
dev.off()  

#######BCR_PLOT########
bcrplotlist1<-list()
j=0
for (i in which(as.numeric(bcr1$Freq)>4)){
  receptor=bcr1$ReceptorID[i]
  plateID=1
  #tcrID=which(tcrlist[[1]]$TCR==receptor)
  mat<-expr_mat(plateID,receptor,use = "bcr")
  j=j+1
  bcrplotlist1[[j]]<-visheatmap1(mat,.title=paste(receptor," in Plate",plateID),.labs = c("", ""),.legend = "log2(TPM)")}
pdf("/Users/yingy_adm/Desktop/bcr_heatmap1.pdf",height=360,width=8)
multiplot(bcrplotlist1[[1]], bcrplotlist1[[2]], bcrplotlist1[[3]],bcrplotlist1[[4]],bcrplotlist1[[5]],
          bcrplotlist1[[6]], bcrplotlist1[[7]],bcrplotlist1[[8]],bcrplotlist1[[9]],bcrplotlist1[[10]],
          bcrplotlist1[[11]],bcrplotlist1[[12]],bcrplotlist1[[13]],bcrplotlist1[[14]],bcrplotlist1[[15]],
          bcrplotlist1[[16]],bcrplotlist1[[17]],bcrplotlist1[[18]],bcrplotlist1[[19]],bcrplotlist1[[20]],
          bcrplotlist1[[21]],bcrplotlist1[[22]],bcrplotlist1[[23]],bcrplotlist1[[24]],bcrplotlist1[[25]],
          bcrplotlist1[[26]],bcrplotlist1[[27]],bcrplotlist1[[28]],bcrplotlist1[[29]],bcrplotlist1[[30]],
          bcrplotlist1[[31]],bcrplotlist1[[32]], bcrplotlist1[[33]], bcrplotlist1[[34]],bcrplotlist1[[35]],
          bcrplotlist1[[36]],bcrplotlist1[[37]],bcrplotlist1[[38]],bcrplotlist1[[39]],bcrplotlist1[[40]],
          bcrplotlist1[[41]], bcrplotlist1[[42]], bcrplotlist1[[43]], bcrplotlist1[[44]],bcrplotlist1[[45]],
          bcrplotlist1[[46]],bcrplotlist1[[47]],bcrplotlist1[[48]],bcrplotlist1[[49]],bcrplotlist1[[50]],
          bcrplotlist1[[51]], bcrplotlist1[[52]], bcrplotlist1[[53]], bcrplotlist1[[54]],bcrplotlist1[[55]],
          bcrplotlist1[[56]],bcrplotlist1[[57]],bcrplotlist1[[58]],bcrplotlist1[[59]],bcrplotlist1[[60]],
          bcrplotlist1[[61]])
dev.off()

bcrplotlist23<-list()
j=0
for (i in which(as.numeric(bcr23$sumfreq)>4)){
  receptor=bcr23$ReceptorID[i]
  plateID1=2
  mat1<-expr_mat(plateID,receptor,use = "bcr")
  j=j+1
  plateID2=3
  mat2<-expr_mat(plateID2,receptor,use="bcr")
  rownames(mat2)<-paste(rownames(mat2),"*")
  mat<-rbind(mat1,mat2)
  bcrplotlist23[[j]]<-visheatmap23(mat,.title=paste(receptor," in Plate",plateID1," + Plate",plateID2,"*"),.labs = c("", ""),.legend = "log2(TPM)")}
pdf("/Users/yingy_adm/Desktop/bcr_heatmap23.pdf",height=615,width=8)
multiplot(bcrplotlist23[[1]], bcrplotlist23[[2]], bcrplotlist23[[3]], bcrplotlist23[[4]],bcrplotlist23[[5]],
          bcrplotlist23[[6]],bcrplotlist23[[7]],bcrplotlist23[[8]],bcrplotlist23[[9]],bcrplotlist23[[10]],
          bcrplotlist23[[11]],bcrplotlist23[[12]],bcrplotlist23[[13]], bcrplotlist23[[14]],bcrplotlist23[[15]],
          bcrplotlist23[[16]],bcrplotlist23[[17]],bcrplotlist23[[18]],bcrplotlist23[[19]],bcrplotlist23[[20]],
          bcrplotlist23[[21]], bcrplotlist23[[22]], bcrplotlist23[[23]], bcrplotlist23[[24]],bcrplotlist23[[25]],
          bcrplotlist23[[26]],bcrplotlist23[[27]],bcrplotlist23[[28]],bcrplotlist23[[29]],bcrplotlist23[[30]],
          bcrplotlist23[[31]], bcrplotlist23[[32]], bcrplotlist23[[33]], bcrplotlist23[[34]],bcrplotlist23[[35]],
          bcrplotlist23[[36]],bcrplotlist23[[37]],bcrplotlist23[[38]],bcrplotlist23[[39]],bcrplotlist23[[40]],
          bcrplotlist23[[41]], bcrplotlist23[[42]], bcrplotlist23[[43]], bcrplotlist23[[44]],bcrplotlist23[[45]],
          bcrplotlist23[[46]],bcrplotlist23[[47]],bcrplotlist23[[48]],bcrplotlist23[[49]],bcrplotlist23[[50]],
          bcrplotlist23[[51]], bcrplotlist23[[52]], bcrplotlist23[[53]], bcrplotlist23[[54]],bcrplotlist23[[55]],
          bcrplotlist23[[56]],bcrplotlist23[[57]],bcrplotlist23[[58]],bcrplotlist23[[59]],bcrplotlist23[[60]],
          bcrplotlist23[[61]], bcrplotlist23[[62]], bcrplotlist23[[63]], bcrplotlist23[[64]],bcrplotlist23[[65]],
          bcrplotlist23[[66]],bcrplotlist23[[67]],bcrplotlist23[[68]])
dev.off()  

bcrplotlist45<-list()
j=0
for (i in which(as.numeric(bcr45$sumfreq)>4)){
  receptor=bcr45$ReceptorID[i]
  plateID1=4
  mat1<-expr_mat(plateID1,receptor,use = "bcr")
  j=j+1
  plateID2=5
  mat2<-expr_mat(plateID2,receptor,use="bcr")
  rownames(mat2)<-paste(rownames(mat2),"*")
  mat<-rbind(mat1,mat2)
  bcrplotlist45[[j]]<-visheatmap45(mat,.title=paste(receptor," in Plate",plateID1," + Plate",plateID2,"*"),.labs = c("", ""),.legend = "log2(TPM)")}
pdf("/Users/yingy_adm/Desktop/bcr_heatmap45.pdf",height=923,width=8)
multiplot(bcrplotlist45[[1]], bcrplotlist45[[2]], bcrplotlist45[[3]], bcrplotlist45[[4]],bcrplotlist45[[5]],
          bcrplotlist45[[6]],bcrplotlist45[[7]],bcrplotlist45[[8]],bcrplotlist45[[9]],bcrplotlist45[[10]],
          bcrplotlist45[[11]],bcrplotlist45[[12]],bcrplotlist45[[13]], bcrplotlist45[[14]],bcrplotlist45[[15]],
          bcrplotlist45[[16]],bcrplotlist45[[17]],bcrplotlist45[[18]],bcrplotlist45[[19]],bcrplotlist45[[20]],
          bcrplotlist45[[21]], bcrplotlist45[[22]], bcrplotlist45[[23]], bcrplotlist45[[24]],bcrplotlist45[[25]],
          bcrplotlist45[[26]],bcrplotlist45[[27]],bcrplotlist45[[28]],bcrplotlist45[[29]],bcrplotlist45[[30]],
          bcrplotlist45[[31]], bcrplotlist45[[32]], bcrplotlist45[[33]], bcrplotlist45[[34]],bcrplotlist45[[35]],
          bcrplotlist45[[36]],bcrplotlist45[[37]],bcrplotlist45[[38]],bcrplotlist45[[39]],bcrplotlist45[[40]],
          bcrplotlist45[[41]], bcrplotlist45[[42]], bcrplotlist45[[43]], bcrplotlist45[[44]],bcrplotlist45[[45]],
          bcrplotlist45[[46]],bcrplotlist45[[47]],bcrplotlist45[[48]],bcrplotlist45[[49]],bcrplotlist45[[50]],
          bcrplotlist45[[51]], bcrplotlist45[[52]], bcrplotlist45[[53]], bcrplotlist45[[54]],bcrplotlist45[[55]],
          bcrplotlist45[[56]],bcrplotlist45[[57]],bcrplotlist45[[58]],bcrplotlist45[[59]],bcrplotlist45[[60]],
          bcrplotlist45[[61]], bcrplotlist45[[62]], bcrplotlist45[[63]], bcrplotlist45[[64]],bcrplotlist45[[65]],
          bcrplotlist45[[66]],bcrplotlist45[[67]],bcrplotlist45[[68]],bcrplotlist45[[69]],bcrplotlist45[[70]],
          bcrplotlist45[[71]], bcrplotlist45[[72]], bcrplotlist45[[73]], bcrplotlist45[[74]],bcrplotlist45[[75]],
          bcrplotlist45[[76]],bcrplotlist45[[77]],bcrplotlist45[[78]],bcrplotlist45[[79]],bcrplotlist45[[80]],
          bcrplotlist45[[81]], bcrplotlist45[[82]], bcrplotlist45[[83]], bcrplotlist45[[84]],bcrplotlist45[[85]],
          bcrplotlist45[[86]],bcrplotlist45[[87]],bcrplotlist45[[88]],bcrplotlist45[[89]],bcrplotlist45[[90]],
          bcrplotlist45[[91]], bcrplotlist45[[92]], bcrplotlist45[[93]], bcrplotlist45[[94]],bcrplotlist45[[95]],
          bcrplotlist45[[96]],bcrplotlist45[[97]],bcrplotlist45[[98]],bcrplotlist45[[99]],bcrplotlist45[[100]],
          bcrplotlist45[[101]], bcrplotlist45[[102]])
dev.off()  

columnplot<-list()
columnplot[[1]]<-barplot(apply(array_sv1, 2, function(x) (sum(x[!is.na(x)])))/apply(array_sv1, 2, function(x) sum(!is.na(x))))
columnplot[[2]]<-barplot(apply(array_sv23, 2, function(x) (sum(x[!is.na(x)])))/apply(array_sv23, 2, function(x) sum(!is.na(x))))
columnplot[[3]]<-barplot(apply(array_sv45, 2, function(x) (sum(x[!is.na(x)])))/apply(array_sv45, 2, function(x) sum(!is.na(x))))
pdf("/Users/yingy_adm/Desktop/columnwise.pdf",height=25,width=8)
multiplot(columnplot[[1]], columnplot[[2]], columnplot[[3]])
dev.off()
##################################################################
#plot for all V_J_BCR genes in plate1
aaplotlist1<-list()
n=1
for (i in which(as.numeric(sbcr1$Freq)>4)){
  aa=sbcr1$CDR3aa[i]
  plateID=1
  geneID=which(sbcrlist[[1]]$CDR3aa==aa)
  mat<-expr_mat(plateID,geneID,use = "aa")
  aaplotlist1[[n]]<-visheatmap1(mat,.title=paste(aa," in Plate",plateID),.labs = c("", ""),.legend = "log2(TPM)")
  n=n+1}
pdf("/Users/yingy_adm/Desktop/aa_plate1.pdf",height=300,width=8)
multiplot(aaplotlist1[[1]],aaplotlist1[[2]],aaplotlist1[[3]],aaplotlist1[[4]],aaplotlist1[[5]],
          aaplotlist1[[6]],aaplotlist1[[7]],aaplotlist1[[8]],aaplotlist1[[9]],aaplotlist1[[10]],
          aaplotlist1[[11]],aaplotlist1[[12]],aaplotlist1[[13]],aaplotlist1[[14]],aaplotlist1[[15]],
          aaplotlist1[[16]],aaplotlist1[[17]],aaplotlist1[[18]],aaplotlist1[[19]],aaplotlist1[[20]],
          aaplotlist1[[21]],aaplotlist1[[22]],aaplotlist1[[23]],aaplotlist1[[24]],aaplotlist1[[25]],
          aaplotlist1[[26]],aaplotlist1[[27]],aaplotlist1[[28]],aaplotlist1[[29]],aaplotlist1[[30]],
          aaplotlist1[[31]],aaplotlist1[[32]],aaplotlist1[[33]],aaplotlist1[[34]],aaplotlist1[[35]],
          aaplotlist1[[36]],aaplotlist1[[37]],aaplotlist1[[38]],aaplotlist1[[39]],aaplotlist1[[40]],
          aaplotlist1[[41]],aaplotlist1[[42]],aaplotlist1[[43]],aaplotlist1[[44]],aaplotlist1[[45]],
          aaplotlist1[[46]],aaplotlist1[[47]],aaplotlist1[[48]],aaplotlist1[[49]],aaplotlist1[[50]],
          aaplotlist1[[51]], aaplotlist1[[52]], aaplotlist1[[53]], aaplotlist1[[54]],aaplotlist1[[55]],
          aaplotlist1[[56]],aaplotlist1[[57]],aaplotlist1[[58]],aaplotlist1[[59]],aaplotlist1[[60]])
dev.off()

aaplotlist23<-list()
for (i in which(as.numeric(sbcr23$sumfreq)>4)){
  aa=sbcr23$CDR3aa[i]
  plateID1=2
  geneID=which(sbcrlist[[2]]$CDR3aa==aa)
  mat1<-expr_mat(plateID1,geneID,use = "aa")
  plateID2=3
  geneID=which(sbcrlist[[3]]$CDR3aa==aa)
  mat2<-expr_mat(plateID2,geneID,use = "aa")
  rownames(mat2)<-paste(rownames(mat2),"*")
  mat<-rbind(mat1,mat2)
  aaplotlist23[[i]]<-visheatmap23(mat,.title=paste(aa," in Plate",plateID1," + Plate",plateID2,"*"),.labs = c("", ""),.legend = "log2(TPM)")}
pdf("/Users/yingy_adm/Desktop/aa_plate23.pdf",height=570,width=8)
multiplot(aaplotlist23[[1]],aaplotlist23[[2]],aaplotlist23[[3]],aaplotlist23[[4]],aaplotlist23[[5]],
          aaplotlist23[[6]],aaplotlist23[[7]],aaplotlist23[[8]],aaplotlist23[[9]],aaplotlist23[[10]],
          aaplotlist23[[11]],aaplotlist23[[12]],aaplotlist23[[13]],aaplotlist23[[14]],aaplotlist23[[15]],
          aaplotlist23[[16]],aaplotlist23[[17]],aaplotlist23[[18]],aaplotlist23[[19]],aaplotlist23[[20]],
          aaplotlist23[[21]],aaplotlist23[[22]],aaplotlist23[[23]],aaplotlist23[[24]],aaplotlist23[[25]],
          aaplotlist23[[26]],aaplotlist23[[27]],aaplotlist23[[28]],aaplotlist23[[29]],aaplotlist23[[30]],
          aaplotlist23[[31]],aaplotlist23[[32]],aaplotlist23[[33]],aaplotlist23[[34]],aaplotlist23[[35]],
          aaplotlist23[[36]],aaplotlist23[[37]],aaplotlist23[[38]],aaplotlist23[[39]],aaplotlist23[[40]],
          aaplotlist23[[41]],aaplotlist23[[42]],aaplotlist23[[43]],aaplotlist23[[44]],aaplotlist23[[45]],
          aaplotlist23[[46]],aaplotlist23[[47]],aaplotlist23[[48]],aaplotlist23[[49]],aaplotlist23[[50]],
          aaplotlist23[[51]],aaplotlist23[[52]],aaplotlist23[[53]],aaplotlist23[[54]],aaplotlist23[[55]],
          aaplotlist23[[56]],aaplotlist23[[57]],aaplotlist23[[58]],aaplotlist23[[59]],aaplotlist23[[60]],
          aaplotlist23[[61]],aaplotlist23[[62]],aaplotlist23[[63]])
dev.off()

aaplotlist45<-list()
for (i in which(as.numeric(sbcr45$sumfreq)>4)){
  aa=sbcr45$CDR3aa[i]
  plateID1=4
  geneID=which(sbcrlist[[4]]$CDR3aa==aa)
  mat1<-expr_mat(plateID1,geneID,use = "aa")
  plateID2=5
  geneID=which(sbcrlist[[5]]$CDR3aa==aa)
  mat2<-expr_mat(plateID2,geneID,use = "aa")
  rownames(mat2)<-paste(rownames(mat2),"*")
  mat<-rbind(mat1,mat2)
  aaplotlist45[[i]]<-visheatmap45(mat,.title=paste(aa," in Plate",plateID1," + Plate",plateID2,"*"),.labs = c("", ""),.legend = "log2(TPM)")}
pdf("/Users/yingy_adm/Desktop/aa_plate45.pdf",height=750,width=8)
multiplot(aaplotlist45[[1]],aaplotlist45[[2]],aaplotlist45[[3]],aaplotlist45[[4]],aaplotlist45[[5]],
          aaplotlist45[[6]],aaplotlist45[[7]],aaplotlist45[[8]],aaplotlist45[[9]],aaplotlist45[[10]],
          aaplotlist45[[11]],aaplotlist45[[12]],aaplotlist45[[13]],aaplotlist45[[14]],aaplotlist45[[15]],
          aaplotlist45[[16]],aaplotlist45[[17]],aaplotlist45[[18]],aaplotlist45[[19]],aaplotlist45[[20]],
          aaplotlist45[[21]],aaplotlist45[[22]],aaplotlist45[[23]],aaplotlist45[[24]],aaplotlist45[[25]],
          aaplotlist45[[26]],aaplotlist45[[27]],aaplotlist45[[28]],aaplotlist45[[29]],aaplotlist45[[30]],
          aaplotlist45[[31]],aaplotlist45[[32]],aaplotlist45[[33]],aaplotlist45[[34]],aaplotlist45[[35]],
          aaplotlist45[[36]],aaplotlist45[[37]],aaplotlist45[[38]],aaplotlist45[[39]],aaplotlist45[[40]],
          aaplotlist45[[41]],aaplotlist45[[42]],aaplotlist45[[43]],aaplotlist45[[44]],aaplotlist45[[45]],
          aaplotlist45[[46]],aaplotlist45[[47]],aaplotlist45[[48]],aaplotlist45[[49]],aaplotlist45[[50]],
          aaplotlist45[[51]],aaplotlist45[[52]],aaplotlist45[[53]],aaplotlist45[[54]],aaplotlist45[[55]],
          aaplotlist45[[56]],aaplotlist45[[57]],aaplotlist45[[58]],aaplotlist45[[59]],aaplotlist45[[60]],
          aaplotlist45[[61]],aaplotlist45[[62]],aaplotlist45[[63]],aaplotlist45[[64]],aaplotlist45[[65]],
          aaplotlist45[[66]],aaplotlist45[[67]],aaplotlist45[[68]],aaplotlist45[[69]],aaplotlist45[[70]],
          aaplotlist45[[71]],aaplotlist45[[72]],aaplotlist45[[73]],aaplotlist45[[74]],aaplotlist45[[75]],
          aaplotlist45[[76]],aaplotlist45[[77]],aaplotlist45[[78]],aaplotlist45[[79]],aaplotlist45[[80]],
          aaplotlist45[[81]],aaplotlist45[[82]],aaplotlist45[[83]])
dev.off()


