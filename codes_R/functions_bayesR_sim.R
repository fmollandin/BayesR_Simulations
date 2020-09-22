# Load bayesR output files, preproccessing


## medium: is there any medium QTL ? wo= are the large qtl in the data? highfirst= how are ranked the simulated QTL in the simul_ped file?
sim_files2<-function(index1,index2,medium=T,wo=F,wo_part=F,highfirst=F,var_tot=100,h2=0.8,ratio_QTL=1){
  
  info_QTL_path=paste("info_QTL_",index1,sep="")
  simperf_path=paste("simperf_",index1,sep="")
  
  if(wo==T){
    index1=paste(index1,"wo",sep="")
    index2=paste(index2,"wo",sep="")
  }
  
  if(wo_part==T){
    index1=paste(index1,"wo",sep="")
  }
  
  scenario=index2
  frq_app_path=paste("sim_",index2,"_app.frq",sep="")
  hyp_app_path=paste("sim_",index2,"_app.hyp",sep="")
  model_app_path=paste("sim_",index2,"_app.model",sep="")
  param_app_path=paste("sim_",index2,"_app.param",sep="")
  gv_app_path=paste("sim_",index2,"_app.gv",sep="")
  gv_test_path=paste("sim_",index2,"_test.gv",sep="")
  fam_app_path=paste("sim_",index1,"_app.fam",sep="")
  bim_app_path=paste("sim_",index1,"_app.bim",sep="")
  fam_test_path=paste("sim_",index1,"_test.fam",sep="")
  bim_test_path=paste("sim_",index1,"_app.bim",sep="")
  
  info_QTL=read.table(info_QTL_path,header = F)
  colnames(info_QTL)=c("nQTL","Chr","Pos_Chr","Pos_All","Var","Frq","Effet")
  sim_perf=read.table(simperf_path,header=F)
  colnames(sim_perf)=c("id","TV")
  
  fam_app<-read.table(fam_app_path,header=F)
  colnames(fam_app)=c('id_famille','id','id_pere','id_mere','sexe','pheno')
  bim_app<-read.table(bim_app_path,header=F)
  colnames(bim_app)=c('chr','id_SNP','cm','bp','allele1','allele2')
  fam_test<-read.table(fam_test_path,header=F)
  colnames(fam_test)=c('id_famille','id','id_pere','id_mere','sexe','pheno')
  bim_test<-read.table(bim_test_path,header=F)
  colnames(bim_test)=c('chr','id_SNP','cm','bp','allele1','allele2')
  
  frq_app<-read.table(frq_app_path,header=F)
  colnames(frq_app)="frequence"
  rownames(frq_app)=bim_app$id_SNP
  
  hyp_app<-read.table(hyp_app_path,header=T,fill=TRUE)
  
  model_app<-read.table(model_app_path,header=F)
  model_app=spread(model_app,key=V1,value=V2)
  
  gv_app<-read.table(gv_app_path,header=F)
  rownames(gv_app)=fam_app$id
  colnames(gv_app)="GV"
  gv_app=data.frame(id=as.numeric(rownames(gv_app)),GV=gv_app[,1])
  
  param_app<-read.table(param_app_path,header=T)
  rownames(param_app)=rownames(frq_app)
  
  gv_test<-read.table(gv_test_path,header=F)
  rownames(gv_test)=fam_test$id
  colnames(gv_test)="GV"
  gv_test=data.frame(id=as.numeric(fam_test$id),gv_test)
  
  if(length(which(duplicated(info_QTL$Pos_All)))!=0){
    info_QTL=info_QTL[-which(duplicated(info_QTL$Pos_All)),]
  }
  
  SNP_info<-read.table("num_map",header=T)
  QTL_info_updated<-left_join(info_QTL,SNP_info,by=c("Pos_All"="rang"))
  QTL_info_updated=QTL_info_updated[,-8]
  
  
  nlow=length(which(QTL_info_updated$Var==names(table(QTL_info_updated$Var))[1]))
  nmedium=length(which(QTL_info_updated$Var==names(table(QTL_info_updated$Var))[2]))
  nhigh=length(which(QTL_info_updated$Var==names(table(QTL_info_updated$Var))[3]))
  
  if (highfirst==F){
    QTL_info_updated=cbind(QTL_info_updated,classe=c(rep(2,nlow),rep(3,nmedium),rep(4,nhigh)))
  } else {
    QTL_info_updated=cbind(QTL_info_updated,classe=c(rep(4,nhigh),rep(3,nmedium),rep(2,nlow)))
  }
  TV_app=merge(gv_app,sim_perf,by="id")
  corr_app=cor(TV_app$GV,TV_app$TV) 
  
  TV_test=merge(gv_test,sim_perf,by="id")
  corr_test=cor(TV_test$GV,TV_test$TV)
  
  
  param_updated=cbind(param_app,id=rownames(param_app))
  param_updated=merge(param_updated,QTL_info_updated[,c(7,8)],all.x=T)
  param_updated=left_join(param_updated, SNP_info)
  param_updated=param_updated[order(param_updated$rang),]
  rownames(param_updated)=param_updated$id
  
  pos_QTL= param_updated[,"rang"]
  
  param_updated=param_updated[,c("PIP1","PIP2","PIP3","PIP4","beta","variance",'Effet')]
  
  QTL_info_updated=QTL_info_updated[order(QTL_info_updated$Pos_All),]
  
  index_high=QTL_info_updated[which(QTL_info_updated$Var==names(table(QTL_info_updated$Var))[3]),]
  QTL_high=QTL_info_updated[which(QTL_info_updated$Var==names(table(QTL_info_updated$Var))[3]),]$id
  QTL_medium=QTL_info_updated[which(QTL_info_updated$Var==names(table(QTL_info_updated$Var))[2]),]$id
  QTL_low=QTL_info_updated[which(QTL_info_updated$Var==names(table(QTL_info_updated$Var))[1]),]$id
  QTL_null=setdiff(bim_app$id_SNP,QTL_info_updated$id)
  
  
  if(medium==F){
    area_array=as.vector(unlist(sapply(filter(QTL_info_updated,classe==3)$Pos_All,function(x) {(x-7):(x+7)})))
    QTL_ref=rep(QTL_medium,each=15)
    QTL_high=QTL_medium
  } else {
    area_array=as.vector(unlist(sapply(filter(QTL_info_updated,classe==4)$Pos_All,function(x) {(x-7):(x+7)})))
    QTL_ref=rep(QTL_high,each=15)
  }
  
  
  QTL_area=filter(SNP_info,rang %in% area_array)$id
  QTL_area=data.frame("id"=as.vector(QTL_area),"Ref"=as.vector(QTL_ref))
  
  narea=length(QTL_area$id)
  
  QTL_null=setdiff(QTL_null,QTL_area$id)
  
  map<-apply(param_updated[,1:4],1,which.max)
  map_nsnp<-apply(cbind(param_updated[,1],rowSums(param_updated[,2:4])),1,which.max)
  null_nsnp=which(map_nsnp==1)
  high_map=which(map==4)
  medium_map=which(map==3)
  low_map=which(map==2)
  null_map=which(map==1)
  
  param_updated$Effet[which(is.na(param_updated$Effet))]=0
  param_updated=cbind(param_updated,True_Cat=rep("Null",nrow(param_updated)),MAP=rep("Null",nrow(param_updated)),True_nsnp=rep("Null",nrow(param_updated)),MAP_nsnp=rep("Nsnp",nrow(param_updated)))
  param_updated$True_Cat=as.character(param_updated$True_Cat)
  param_updated$MAP=as.character(param_updated$MAP)
  param_updated$True_nsnp=as.character(param_updated$True_nsnp)
  param_updated$MAP_nsnp=as.character(param_updated$MAP_nsnp)
  
  param_updated$True_Cat[which(rownames(param_updated) %in% QTL_low)]="Low"
  param_updated$True_Cat[which(rownames(param_updated) %in% QTL_area$id)]="Area"
  param_updated$True_Cat[which(rownames(param_updated) %in% QTL_medium)]="Medium"
  param_updated$True_Cat[which(rownames(param_updated) %in% QTL_high)]="High"
  
  
  param_updated$MAP[low_map]="Low"
  param_updated$MAP[medium_map]="Medium"
  param_updated$MAP[high_map]="High"
  param_updated$MAP_nsnp[null_nsnp]="Null"
  param_updated$True_nsnp[which(rownames(param_updated) %in% QTL_high)]="Nsnp"
  param_updated$True_nsnp[which(rownames(param_updated) %in% QTL_area$id)]="Nsnp"
  
  param_updated=cbind(param_updated,Frq=frq_app)
  
  sumvar<-SumVar_Int(param_updated[,paste("PIP",1:4,sep="")])
  param_updated=cbind(param_updated,SumVar=sumvar)
  
  
  if(medium==F){
    param_updated$True_Cat[which(param_updated$True_Cat=="Medium")]="High"
  }
  
  param_updated=cbind(param_updated,"Position"=pos_QTL)
  
  param_updated=cbind(param_updated,V_i=(param_updated$beta^2+param_updated$variance))
  param_updated=cbind(param_updated,Part_Vi=param_updated$V_i/sum(param_updated$V_i))
  
  temp=cbind(id=rownames(param_updated),param_updated)
  param_QTL=inner_join(QTL_area,temp)
  rm(temp)
  
  if(highfirst==T){
    id=rownames(param_updated)[which(param_updated$True_Cat=="High")]
    QTL=filter(param_updated,True_Cat=="High")
    pi3<-info_QTL$Var[1]/(var_tot*ratio_QTL*h2)
    pi3_QTL<-pi3*2*QTL$frequence*(1-QTL$frequence)
    pi3_QTL=data.frame(id=id,pi3_True=pi3_QTL)
  } else {
    id=rownames(param_updated)[which(param_updated$True_Cat=="High")]
    QTL=filter(param_updated,True_Cat=="High")
    pi3<-info_QTL$Var[nrow(info_QTL)]/(var_tot*ratio_QTL*h2)
    pi3_QTL<-pi3*2*QTL$frequence*(1-QTL$frequence)
    pi3_QTL=data.frame(id=id,pi3_True=pi3_QTL)
  }
  
  
  param_ordered=cbind(param_updated,beta_abs=abs(param_updated$beta))
  param_ordered=param_ordered[order(abs(param_updated$beta),decreasing  =TRUE),]
  
  out=list(scenario,info_QTL,sim_perf,frq_app,hyp_app,model_app,gv_app,param_updated,param_ordered,gv_test,corr_app,corr_test,TV_test,TV_app,param_QTL,narea,pi3_QTL)
  names(out)=c("scenario","info_QTL","sim_perf","frq_app","hyp_app","model_app","gv_app","param_app","param_ordered","gv_test","cor_app","cor_test","TV_test","TV_app","param_QTL","narea","pi3_QTL")
  return(out)
}

## load multiple file numbered
sim_multi_files<-function(ind1,ind2,medium=F,wo=F,wo_part=F,highfirst=F){
  out=list()
  for(i in 1:length(ind1)){
    out[[i]]=sim_files2(ind1[i],ind2[i],medium,wo,wo_part,highfirst)
  }
  return(out)
}

## compare results of different scenarios in a dataframe
sim_comp<-function(liste){
  nlist=length(liste)
  tab_comp=data.frame(matrix(0,16,nlist))
  colnames(tab_comp)=sapply(liste,function(x) {x$'scenario'})
  rownames(tab_comp)=c("Cor_Training","Cor_Validation",paste("Nk",1:4,sep=""),"Nsnp",paste("Eff_MAP",1:4,sep=""),'MAP_Nsnp',"Ratio_TP_Nsnp","Ratio_TP_High","Sensitivity_Nsnp","Sensitivity_High")
  tab_comp[1,]=sapply(liste,function(x) {x$'cor_app'})
  tab_comp[2,]=sapply(liste,function(x) {x$'cor_test'})
  tab_comp[3,]=sapply(liste,function(x) {x$'model_app'$'Nk1'})
  tab_comp[4,]=sapply(liste,function(x) {x$'model_app'$'Nk2'})
  tab_comp[5,]=sapply(liste,function(x) {x$'model_app'$'Nk3'})
  tab_comp[6,]=sapply(liste,function(x) {x$'model_app'$'Nk4'})
  tab_comp[7,]=sapply(liste,function(x) {x$'model_app'$'Nsnp'})
  tab_comp[8,]=sapply(liste,function(x) {length(which(x$param_app$MAP=='Null'))})
  tab_comp[9,]=sapply(liste,function(x) {length(which(x$param_app$MAP=='Low'))})
  tab_comp[10,]=sapply(liste,function(x) {length(which(x$param_app$MAP=='Medium'))})
  tab_comp[11,]=sapply(liste,function(x) {length(which(x$param_app$MAP=='High'))})
  tab_comp[12,]=sapply(liste,function(x) {length(which(x$param_app$MAP_nsnp=='Nsnp'))})
  tab_comp[13,]=sapply(liste,function(x) {length(which(x$param_app$MAP_nsnp=='Nsnp' & x$param_app$True_nsnp=='Nsnp'))/length(which(x$param_app$MAP_nsnp=='Nsnp'))})
  tab_comp[13,which(is.na(tab_comp[13,]))]=0   # avoid na
  tab_comp[14,]=sapply(liste,function(x) {length(which(x$param_app$MAP=='High' & x$param_app$True_nsnp=="Nsnp"))/length(which(x$param_app$MAP=='High'))})
  tab_comp[14,which(is.na(tab_comp[14,]))]=0   # avoid na
  tab_comp[15,]=sapply(liste,function(x) {length(which(x$param_app$MAP_nsnp=='Nsnp' & x$param_app$True_nsnp=="Nsnp"))/length(which(x$param_app$True_nsnp=="Nsnp"))})
  tab_comp[16,]=sapply(liste,function(x) {length(which(x$param_app$MAP=='High' & x$param_app$True_nsnp=="Nsnp"))/length(which(x$param_app$True_nsnp=="Nsnp"))}) 
  return(tab_comp)
}

