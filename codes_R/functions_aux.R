hyp_hist<-function(data_hyp,hyper,bin=30){
  ggplot(data_hyp,aes_string(x=hyper))+
    geom_histogram(aes(y=..density..),color="red",fill="white",bins=bin)+
    geom_density(alpha=.2,fill="#FF6666",color="red")
}

hyp_plot<-function(data_hyp,hyper,data_mean=NA){
  if(!is.na(data_mean)){
    ggplot(data_hyp,aes_string(x="Replicate",y=hyper,colour=hyper))+
      geom_line()+
      geom_hline(data=data_mean,aes_string(yintercept=hyper))
  }else{
    ggplot(data_hyp,aes_string(x="Replicate",y=hyper,colour=hyper))+
      geom_line()
  }
}

hyp_hist_comp<-function(data_hyp,hyper,bin=30){
  a=paste0(hyper,1:4)
  temp<-gather(data_hyp[,a], key, Effectifs) %>%  select(-key)
  temp=cbind(temp,"Groupe"=rep(1:4,each=nrow(data_hyp)))
  temp[,2]=as.factor(temp[,2])
  ggplot(temp, aes(x=Effectifs, color=Groupe,fill=Groupe)) +
    geom_histogram(alpha=0.3, position="identity",bins=bin)+
    theme(legend.position="top")+
    xlab(hyper)
}

effet_chr<-function(data,chromosome,seuil){
  temp=data[which(data$chr==chromosome),]
  ggplot(temp,aes(x=bp,y=effet))+      
    geom_point()+
    geom_text(aes(label=ifelse(abs(effet)>seuil,rownames(temp),'')))+
    ggtitle(paste0("Effet estimé des SNPs du chromosome ",chromosome))
}  


top_SNP<-function(data_SNP,n){
  temp=data_SNP[1:n,]
  ggplot(temp,aes(x=as.factor(n:1),y=effet,label=rownames(temp)))+
    geom_point(color="red")+
    theme(axis.text.x = element_blank(),axis.ticks.x=element_blank())+
    xlab("SNP")+
    geom_text(size=3)
}

min_SNP<-function(data_SNP,n){
  m=nrow(data_SNP)
  temp=data_SNP[(m-n+1):m,]
  ggplot(temp,aes(x=as.factor(1:n),y=effet,label=rownames(temp)))+
    geom_point(color="red")+
    theme(axis.text.x = element_blank(),axis.ticks.x=element_blank())+
    xlab("SNP")+
    geom_text(size=3)
}

sim_files<-function(index,perf_index){
  frq_app_path=paste("sim_",index,"_app.frq",sep="")
  hyp_app_path=paste("sim_",index,"_app.hyp",sep="")
  model_app_path=paste("sim_",index,"_app.model",sep="")
  param_app_path=paste("sim_",index,"_app.param",sep="")
  gv_app_path=paste("sim_",index,"_app.gv",sep="")
  gv_test_path=paste("sim_",index,"_test.gv",sep="")
  fam_app_path=paste("sim_",index,"_app.fam",sep="")
  bim_app_path=paste("sim_",index,"_app.bim",sep="")
  fam_test_path=paste("sim_",index,"_test.fam",sep="")
  bim_test_path=paste("sim_",index,"_app.bim",sep="")
 
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
  
  perf=read.table(perf_index)
  colnames(perf)=c("id","sexe","Fiabilite_fille","Pheno","Set")
  
  param_app<-read.table(param_app_path,header=T)
  rownames(param_app)=rownames(frq_app)
  
  param_app=cbind(param_app,V_i=param_app$variance+param_app$beta^2)
  
  
  sumvar<-SumVar_Int(param_app[,paste("PIP",1:4,sep="")])
  param_app=cbind(param_app,SumInt=sumvar)
  
  
  map<-apply(param_app[,1:4],1,which.max)
  map_nsnp<-apply(cbind(param_app[,1],rowSums(param_app[,2:4])),1,which.max)
  null_nsnp=which(map_nsnp==1)
  high_map=which(map==4)
  medium_map=which(map==3)
  low_map=which(map==2)
  null_map=which(map==1)
  
  param_app=cbind(param_app,MAP=rep("Null",nrow(param_app)),MAP_nsnp=rep("Nsnp",nrow(param_app)))
  param_app$MAP=as.character(param_app$MAP)
  param_app$MAP_nsnp=as.character(param_app$MAP_nsnp)
  
  param_app$MAP[low_map]="Low"
  param_app$MAP[medium_map]="Medium"
  param_app$MAP[high_map]="High"
  param_app$MAP_nsnp[null_nsnp]="Null"
  
  gv_test<-read.table(gv_test_path,header=F)
  rownames(gv_test)=fam_test$id
  colnames(gv_test)="GV"
  gv_test=data.frame(id=as.numeric(fam_test$id),gv_test)
  
  TV_app=merge(gv_app,perf,by="id")
  corr_app=cor(TV_app$GV,TV_app$Pheno) 
  
  TV_test=merge(gv_test,perf,by="id")
  corr_test=cor(TV_test$GV,TV_test$Pheno)
  
  out=list(fam_app,bim_app,fam_test,bim_test,frq_app,hyp_app,model_app,gv_app,param_app,gv_test,corr_test,corr_app,TV_app,TV_test)
  names(out)=c("fam_app","bim_app","fam_test","bim_test","frq_app","hyp_app","model_app","gv_app","param_app","gv_test","cor_test","cor_app","TV_app","TV_test")
  return(out)
}


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
  
  
  #index_QTL=which( rownames(param_app) %in% QTL_info_updated$id)
 
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
 # param_updated$True_Effect[index_QTL]=arrange(QTL_info_updated,Pos_All)$Effet
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
  
  
  sumvar_MAF<-SumVar_Int(PIP=param_updated[,paste("PIP",1:4,sep="")],frq=param_updated$frequence)
  param_updated=cbind(param_updated,SumVar_MAF=sumvar_MAF)
  
  param_updated$post_var=param_updated$beta^2*2*param_updated$frequence*(1-param_updated$frequence)
  model_app$post_var_tot=sum(param_updated$post_var)
  param_updated$post_part_var=param_updated$post_var/model_app$post_var_tot
  
  
  if(medium==F){
    param_updated$True_Cat[which(param_updated$True_Cat=="Medium")]="High"
  }
  
  param_updated=cbind(param_updated,"Position"=pos_QTL)
  
  param_updated=cbind(param_updated,V_i=(param_updated$beta^2+param_updated$variance))
  
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

sim_multi_files<-function(ind1,ind2,medium=F,wo=F,wo_part=F,highfirst=F){
  out=list()
  for(i in 1:length(ind1)){
    out[[i]]=sim_files2(ind1[i],ind2[i],medium,wo,wo_part,highfirst)
  }
  return(out)
}

sim_mean_comp<-function(liste,wo=F){
  comp_mean=data.frame(sapply(liste,function(x) {rowMeans(x)}))
  if (wo==F){
    colnames(comp_mean)=sapply(liste,function(x) {str_sub(colnames(x)[1],end=nchar(colnames(x)[1])-1)})
  } else {
    colnames(comp_mean)=sapply(liste,function(x) {str_sub(colnames(x)[1],end=nchar(colnames(x)[1])-3)})  
  }
  return(comp_mean)
}

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
  tab_comp[13,which(is.na(tab_comp[13,]))]=0   # sinon pb de calcul de moyenne
  tab_comp[14,]=sapply(liste,function(x) {length(which(x$param_app$MAP=='High' & x$param_app$True_nsnp=="Nsnp"))/length(which(x$param_app$MAP=='High'))})
  tab_comp[14,which(is.na(tab_comp[14,]))]=0 
  tab_comp[15,]=sapply(liste,function(x) {length(which(x$param_app$MAP_nsnp=='Nsnp' & x$param_app$True_nsnp=="Nsnp"))/length(which(x$param_app$True_nsnp=="Nsnp"))})
  tab_comp[16,]=sapply(liste,function(x) {length(which(x$param_app$MAP=='High' & x$param_app$True_nsnp=="Nsnp"))/length(which(x$param_app$True_nsnp=="Nsnp"))}) #pb avec AREA
  return(tab_comp)
}



AUC_plot<-function(liste,nmin,nmax){
  AUC_comp=matrix(0,length(nmin:nmax),length(liste))
  rownames(AUC_comp)=nmin:nmax
  colnames(AUC_comp)=sapply(liste,function(x) {x$'scenario'})
  j=1
  for (i in nmin:nmax){
    AUC_comp[j,]=sapply(liste,function(x) {roc(x$param_ordered[1:i,"True_nsnp"],x$param_ordered[1:i,"beta_abs"])$auc[[1]]})
    j=j+1
  }
  AUC_melted=melt(AUC_comp,1)
  colnames(AUC_melted)=c("Top","Scenario","AUC")
  GG<-ggplot(AUC_melted,aes(x=Top,y=AUC,group=Scenario,color=Scenario)) +geom_line()
  return(GG)
}

intersect_top_SNP<-function(liste,n){
  intersect_top<-matrix(0,length(liste),length(liste))
  rownames(intersect_top)=lapply(liste,function(x){x$scenario})
  colnames(intersect_top)=lapply(liste,function(x){x$scenario})
  for (i in 1:length(liste)){
    for (j in 1:i){
      a=length(intersect(rownames(liste[[i]]$param_ordered)[1:n],rownames(liste[[j]]$param_ordered)[1:n]))
      intersect_top[i,j]=a
      intersect_top[j,i]=a
     }
  }
  return(intersect_top)
}

heatmap_matrix<-function(matrice){
  melted_matrix=melt(matrice)
  GG<-ggplot(melted_matrix,aes(x = Var1,y=Var2,fill=value))+
      geom_tile()
}

df_PartVar_QTL<-function(liste,gpin=c(0,0.0001,0.001,0.01),Cat="High"){
  df_PartVar=data.frame()
  j=0
  for (i in liste){
    j=j+1
    temp=sapply(i,function(x){as.matrix(filter(x$param_app,True_Cat==Cat)[,paste("PIP",1:4,sep="")])%*%gpin})
    df_temp=data.frame(SumInt=as.vector(temp),Scenario=rep(names(liste)[j],length(temp)))
    df_PartVar=rbind(df_PartVar,df_temp)
    rm(temp)
    rm(df_temp)
  }
  return(df_PartVar)
}

SumVar_Int<-function(PIP,nint=15,gpin=c(0,0.0001,0.001,0.01),frq=NA){
  n=nrow(PIP)
  
  if (is.na(frq)){
    part_var=as.matrix(PIP) %*% gpin  
  } else {
    part_var=as.matrix(PIP) %*% gpin * (2*frq*(1-frq))
  }
  
  part_int=rep(0,n-2*trunc(nint))
  
  for(i in (1+trunc(nint/2)):(n-trunc(nint/2))){
    part_int[i-trunc(nint/2)]=sum(part_var[(i-trunc(nint/2)):(i+trunc(nint/2))])
  }
  
  part_int=c(rep(min(part_int),trunc(nint/2)),part_int,rep(min(part_int),trunc(nint/2)))
  return(part_int)
}

SumVi=function(out,nint=15){
  n=nrow(out$param_app)
  SumV_i=rep(0,n-2*trunc(nint))
  
  for(i in (1+trunc(nint/2)):(n-trunc(nint/2))){
    SumV_i[i-trunc(nint/2)]=sum(out$param_app$V_i[(i-trunc(nint/2)):(i+trunc(nint/2))])
  }
  
  SumV_i=c(rep(min(SumV_i),trunc(nint/2)),SumV_i,rep(min(SumV_i),trunc(nint/2)))
  return(SumV_i)
}

QTL_Int<-function(liste,seuil=0.001){
  nQTL=sapply(liste,function(x) {length(uqnique(filter(x$param_QTL,SumVar>seuil)$Ref))})
  nSNP=sapply(liste,function(x) {nrow(filter(x$param_app,SumVar>seuil))})
  return(list("nQTL"=nQTL,"nSNP"=nSNP))
}

QTL_comp<-function(liste,seuil=0.001){
  df_QTL_comp=data.frame(Scenario=rep(names(liste),each=2*length(liste[[1]])),Effectif=unlist(sapply(liste,function(x){QTL_Int(x,seuil)})),Categorie=rep(rep(rownames(sapply(liste,function(x){QTL_Int(x)})),each=length(liste[[1]])),length(liste)))
  return(df_QTL_comp)
}

AUC_top<-function(liste,n=50){
  AUC_top<-sapply(liste,function(x) {roc(response=x$param_ordered[1:n,c("True_nsnp")],predictor=x$param_ordered[1:n,c("beta_abs")],levels=c("Nsnp","Null"))$auc[[1]]})
  return(AUC_top)
}

AUC_comp<-function(liste,namesl,n=50){
  df_QTL_comp=data.frame(Scenario=rep(namesl,each=length(liste[[1]])),AUC=as.vector(unlist(sapply(liste,function(x){AUC_top(x,n)}))))
  return(df_QTL_comp)
}

nqtl<-function(out,type){
  nqtl<-nrow(filter(out$param_app,True_Cat==type))
  return(nqtl)
}


df_max_int<-function(liste){
  df_max=data.frame()
  for (i in 1:length(liste)){
    temp=as.vector(sapply(liste[[i]],function(x) {aggregate(x$param_QTL$V_i,list(x$param_QTL$Ref),max)$x}))
    temp=data.frame(Vi=temp,Scenario=names(liste)[[i]])
    df_max=rbind(df_max,temp)
  }
  return(df_max)
}



max_int_var<-function(out,nQTL){
  max_int=c()
  for (i in 1:nQTL){
    max_int=c(max_int,max(filter(out$param_QTL,Ref==names(table(out$param_QTL$Ref))[i])$variance))
  }
  return(max_int)
}

max_int_sum<-function(out,nQTL,pondere=F){
  max_int=c()
  if (pondere==F){
    for (i in 1:nQTL){
      max_int=c(max_int,max(filter(out$param_QTL,Ref==names(table(out$param_QTL$Ref))[i])$SumVar))
    }
  } else {
    for (i in 1:nQTL){
      max_int=c(max_int,max(filter(out$param_QTL,Ref==names(table(out$param_QTL$Ref))[i])$SumVar_MAF))
    }
  }
  return(max_int)
}

which_max_int<-function(out,nQTL=5){
  which_max=c()
  for (i in 1:nQTL){
    temp=filter(out$param_QTL,Ref==names(table(out$param_QTL$Ref))[i])
    which_max=c(which_max,temp$id[which.max(abs(temp$beta))])
  }
  return(which_max)
}

narea_detected<-function(out){
  nMAP<-length(unique(filter(out$param_QTL,MAP_nsnp=="Nsnp")$Ref))
  return(nMAP)
}

df_QTL_detected<-function(liste,area=T){
  if (area==T){
    df_QTL=data.frame()
    for (i in 1:length(liste)){
      temp<-data.frame(QTL=unlist(lapply(liste[[i]],function(x) {unique(x$param_QTL$Ref)})))
      temp=cbind.data.frame(temp,Detected=rep("No",nrow(temp)),Scenario=names(liste)[[i]],stringsAsFactors = FALSE)
      temp_QTL=unlist(lapply(liste[[i]],function(x) {unique(filter(x$param_QTL,MAP_nsnp=="Nsnp")$Ref)}))
      temp[which(temp$QTL %in% temp_QTL),"Detected"]="Yes"
      df_QTL=rbind(df_QTL,temp)
      rm(temp)
    }
  } else {
    df_QTL=data.frame()
    for (i in 1:length(liste)){
      temp<-data.frame(QTL=unlist(lapply(liste[[i]],function(x) {unique(x$param_QTL$Ref)})))
      temp_MAF=unlist(lapply(liste[[i]],function(x) {filter(x$param_app,True_Cat=="High")$frequence}))
      temp=cbind.data.frame(temp,MAF=1-temp_MAF,Detected=rep("No",nrow(temp)),Scenario=names(liste)[[i]],stringsAsFactors = FALSE)
      temp_QTL=unlist(lapply(liste[[i]],function(x) {unique(filter(filter(x$param_QTL,MAP_nsnp=="Nsnp"),True_Cat=="High")$Ref)}))
      temp$Detected[temp_QTL]="Yes"
      df_QTL=rbind(df_QTL,temp)
      rm(temp)
    }
  }
  return(df_QTL)
}


df_QTL_detected<-function(liste,area=T){
  if (area==T){
    df_QTL=data.frame()
    for (i in 1:length(liste)){
      temp<-data.frame(QTL=unlist(lapply(liste[[i]],function(x) {unique(x$param_QTL$Ref)})))
      temp=cbind.data.frame(temp,Detected=rep("No",nrow(temp)),Scenario=names(liste)[[i]],stringsAsFactors = FALSE)
      temp_QTL=unlist(lapply(liste[[i]],function(x) {unique(filter(x$param_QTL,MAP_nsnp=="Nsnp")$Ref)}))
      temp$Detected[temp_QTL]="Yes"
      df_QTL=rbind(df_QTL,temp)
      rm(temp)
    }
  } else {
    df_QTL=data.frame()
    for (i in 1:length(liste)){
      temp<-data.frame(QTL=unlist(lapply(liste[[i]],function(x) {unique(x$param_QTL$Ref)})))
      temp_MAF=unlist(lapply(liste[[i]],function(x) {filter(x$param_app,True_Cat=="High")$frequence}))
      temp=cbind.data.frame(temp,MAF=1-temp_MAF,Detected=rep("No",nrow(temp)),Scenario=names(liste)[[i]],stringsAsFactors = FALSE)
      temp_QTL=unlist(lapply(liste[[i]],function(x) {unique(filter(filter(x$param_QTL,MAP_nsnp=="Nsnp"),True_Cat=="High")$Ref)}))
      temp$Detected[temp_QTL]="Yes"
      df_QTL=rbind(df_QTL,temp)
      rm(temp)
    }
  }
  return(df_QTL)
}

df_MAP_dist<-function(liste,nQTL=50){
  df_comp=data.frame()
  for (i in 1:length(liste)){
    temp_ind_map=data.frame(table(sapply(liste[[i]], function(x) {filter(x$param_QTL,True_Cat=="High")$MAP}))/nQTL)
    temp_ind_map=cbind(temp_ind_map,Scale=rep("QTL",nrow(temp_ind_map)),MAP=rep("MAP",nrow(temp_ind_map)))
    
    temp_ind_nsnp=data.frame(table(sapply(liste[[i]], function(x) {filter(x$param_QTL,True_Cat=="High")$MAP_nsnp}))/nQTL)
    temp_ind_nsnp=cbind(temp_ind_nsnp,Scale=rep("QTL",nrow(temp_ind_nsnp)),MAP=rep("MAP_nsnp",nrow(temp_ind_nsnp)))
    
    temp_comp=rbind(temp_ind_map,temp_ind_nsnp)
    temp_comp=cbind(temp_comp,Scenario=names(liste)[[i]])
    
    df_comp=rbind(df_comp,temp_comp)
  }
  return(df_comp)
}

area_stat<-function(out,nQTL=5){
  ntemp=narea_detected(out)
  joint_temp=left_join(filter(out$param_app,MAP_nsnp=="Nsnp"),out$param_QTL)
  MAP_unique=length(which(is.na(joint_temp$Ref)))+ntemp
  if (MAP_unique!=0){
    sensibility<-ntemp/MAP_unique
  } else {
    sensibility=0
  }
  sensitivity<-ntemp/nQTL
  return(c(sensibility,sensitivity))
}

stat_andrea<-function(out,nint=15,high=F){
  if (high==F){
    n=nrow(out$param_app)
    ratio_int=rep(0,n-2*trunc(nint))
    
    for(i in (1+trunc(nint/2)):(n-trunc(nint/2))){
      ratio_int[i-trunc(nint/2)]=out$param_app$SumVar[i]/mean(out$param_app$SumVar[(i-trunc(nint/2)):(i+trunc(nint/2))])
    }
    ratio_int=c(rep(min(ratio_int),trunc(nint/2)),ratio_int,rep(min(ratio_int),trunc(nint/2)))
      
  } else {
    forts<-which(out$param_app$True_Cat=="High")
    ratio_int=c()
    for (i in forts){
      ratio_int=cbind(ratio_int,out$param_app$SumVar[i]/mean(out$param_app$SumVar[(i-trunc(nint/2)):(i+trunc(nint/2))]))
    }
  }
  return(ratio_int)
}  

df_andrea_QTL<-function(liste,nint=15){
  df_andrea=data.frame()
  for (i in 1:length(liste)){
    temp=sapply(liste[[i]], function(x) {stat_andrea(x,nint=nint,high=T)})
    temp=cbind.data.frame(Scenario=rep(names(liste)[[i]],length(temp)),Pic=as.vector(temp))
    df_andrea=rbind.data.frame(df_andrea,temp)
    rm(temp)
  }
  return(df_andrea)
}


detection_andrea<-function(out,nint=15,seuil=5){
  QTL_andrea<-rownames(out$param_app[which(stat_andrea(out,nint)>seuil),])
  return(QTL_andrea)
}

df_PIP<-function(out,cat="High",FUN="moyenne"){
  if (FUN=="moyenne"){
    df_PIP1<-sapply(out,function(x) {colMeans(filter(x$param_app,True_Cat==cat)[,1:4])})
  }
  if (FUN=="mediane"){
    df_PIP1<-sapply(out,function (x) {apply(filter(x$param_app,True_Cat==cat)[,1:4],2,median)})
  }
  colnames(df_PIP1)=as.character(1:length(out))
  df_PIP1=gather(as.data.frame(df_PIP1),key="Simulation")
  df_PIP1=cbind(df_PIP1,PIP=rep(paste("PIP",1:4,sep=""),length(out)))
  return(df_PIP1)
}

df_PIP_list<-function(liste,cat1="High",FUN1="moyenne"){
  df_PIP_all<-data.frame()
  if (FUN1=="moyenne") {
    for (i in 1:length(liste)){
    temp=data.frame(df_PIP(liste[[i]],cat1,FUN=FUN1))
    temp=aggregate(temp$value,list(temp$PIP),mean)
    colnames(temp)=c("PIP","Moyenne")
    temp=cbind(temp,Scenario=rep(names(liste)[[i]],nrow(temp)))
    df_PIP_all=rbind(df_PIP_all,temp)
    rm(temp)
    }
  }
  if (FUN1=="mediane"){
    for (i in 1:length(liste)){
      temp=data.frame(df_PIP(liste[[i]],cat1,FUN=FUN1))
      temp=aggregate(temp$value,list(temp$PIP),median)
      colnames(temp)=c("PIP","Mediane")
      temp=cbind(temp,Scenario=rep(names(liste)[[i]],nrow(temp)))
      df_PIP_all=rbind(df_PIP_all,temp)
      rm(temp)
    }
  }
  return(df_PIP_all)
}

df_PIP_v2<-function(out){
  df_PIP1<-sapply(out,function(x) {colMeans((x$param_QTL %>% group_by(Ref)  %>% slice(which.min(PIP1))   %>% as.data.frame)[,paste("PIP",1:4,sep='')])})
  df_PIP1=gather(data.frame(t(df_PIP1)),key=PIP,value=Moyenne)
  df_PIP1=aggregate(df_PIP1$Moyenne,list(df_PIP1$PIP),mean)
  colnames(df_PIP1)=c("PIP","Moyenne")
  return(df_PIP1)
}

df_PIP_all_v2<-function(liste){
  df_PIP_all=data.frame()
  for (i in 1:length(liste)){
    temp=df_PIP_v2(liste[[i]])
    temp=cbind(temp,Scenario=rep(names(liste)[[i]],nrow(temp)))
    df_PIP_all=rbind(df_PIP_all,temp)
    rm(temp)
  }
  return(df_PIP_all)
}


df_PIP_med<-function(liste){
  df_PIP_all=matrix(ncol=4,nrow=length(liste))
  df_PIP_all=data.frame(df_PIP_all)
  colnames(df_PIP_all)=c("PIP1","PIP2","PIP3","PIP4")
  for (i in 1:length(liste)){
    for (j in 1:4){
      df_PIP_all[i,j]=median(sapply(liste[[i]],function(x) {matrix(filter(x$param_app,True_Cat=="High")[,j])}))
    }
  }
  df_PIP_all=gather(df_PIP_all,key="PIP",value="Mediane")
  df_PIP_all=cbind(df_PIP_all,Scenario=rep(names(liste),4))
  return(df_PIP_all)
}


best_int<-function(out){
  best_int=c()
  for (i in 1:length(unique(out$param_QTL$Ref))){
    temp=filter(out$param_QTL,Ref==unique(out$param_QTL$Ref)[i])
    best_int=c(best_int,temp$id[which.min(temp$PIP1)])
    rm(temp)
  }
  return(best_int)
}

SumInt_MAF<-function(out){
  SumInt2<-out$param_app$SumVar/(2*(1-out$param_app$frequence)*out$param_app$frequence)
  return(SumInt2)
}


pi3_MAF<-function(out,h2,ratio_QTL=1,var_tot=100){
  QTL=filter(out$param_app,True_Cat=="High")
  pi3<-out$info_QTL$Var[1]/(var_tot*ratio_QTL*h2)
  pi3_SNP<-pi3*2*QTL$frequence*(1-QTL$frequence)
  return(pi3_SNP)
}

df_pi3_merge<-function(out){
  df_pi3<-merge(out$param_QTL,out$pi3_QTL,by.x="Ref",by.y="id")
  return(df_pi3)
}

df_pi3_PIP<-function(liste){
  df_pi3_PIP<-do.call("rbind",lapply(out_R,function(x) {filter(df_pi3_merge(x),True_Cat=="High")[,c(paste("PIP",1:4,sep=""),"pi3_True")]}))
  return(df_pi3_PIP)
}

mean_pi3_true<-function(liste){
  pi3_mean=as.vector(sapply(liste_out,function(y) {mean(as.vector(sapply(y,function(x) {x$pi3_QTL$pi3_True})))}))
  scenario=names(liste)
  df_pi3_mean=data.frame(Scenario=scenario,pi3_mean=pi3_mean)
  return(df_pi3_mean)
}

ratio_SumInt<-function(liste,from="AC",Categorie="High"){
  base=which(names(liste)==from)
  df_ratio=data.frame()
  for (i in which(names(liste)!=from)){
    temp=data.frame(Ratio=mean(unlist(sapply(1:length(liste[[1]]),function(x) {filter(liste[[base]][[x]]$param_app,True_Cat==Categorie)$SumVar_MAF/filter(liste[[i]][[x]]$param_app,True_Cat==Categorie)$SumVar_MAF}))),Base=from,Scenario=names(liste)[[i]])
    df_ratio=rbind(df_ratio,temp)
    rm(temp)
  }
  return(df_ratio)
}

df_post_var<-function(liste){
  df_var=data.frame()
  for (i in 1:length(liste)){
    temp=lapply(liste[[i]],function(x) {filter(x$param_app,True_Cat=="High")$post_var})
    temp=data.frame(Post_Var=unlist(temp),Variance=rep(liste[[i]][[1]]$info_QTL$Var[1],length(temp)))
    df_var=rbind(temp,df_var)
  }
  return(df_var)
}

df_part_var<-function(liste){
  df_var=data.frame()
  for (i in 1:length(liste)){
    temp=lapply(liste[[i]],function(x) {filter(x$param_app,True_Cat=="High")$post_part_var})
    temp=data.frame(Part_Var=unlist(temp),Variance=rep(liste[[i]][[1]]$info_QTL$Var[1]/80,length(temp)))
    df_var=rbind(temp,df_var)
  }
  return(df_var)
}

df_part_high<-function(liste){
  df_var=data.frame()
  for (i in 1:length(liste)){
    temp=lapply(liste[[i]],function(x) {sum(filter(x$param_app,True_Cat=="High")$post_part_var)})
    temp=data.frame(Part_High=unlist(temp),Variance=rep(5*liste[[i]][[1]]$info_QTL$Var[1]/80,length(temp)),Repetition=as.character(1:10))
    df_var=rbind(temp,df_var)
  }
  return(df_var)
}

Sum_PostVar<-function(liste){
  df_Sum<-data.frame()
  for (i in 1:length(liste)){
    temp=lapply(liste[[i]], function(x) {aggregate(x$param_QTL$post_var,list(x$param_QTL$Ref),sum)$x})
    temp=data.frame(SumVar=unlist(temp),Variance=liste[[i]][[1]]$info_QTL$Var[1])
    df_Sum=rbind(temp,df_Sum)
    rm(temp)
  }
  return(df_Sum)
}

Max_PostVar<-function(liste){
  df_Max<-data.frame()
  for (i in 1:length(liste)){
    temp=lapply(liste[[i]], function(x) {aggregate(x$param_QTL$post_var,list(x$param_QTL$Ref),max)$x})
    temp=data.frame(MaxVar=unlist(temp),Variance=liste[[i]][[1]]$info_QTL$Var[1])
    df_Max=rbind(temp,df_Max)
    rm(temp)
  }
  return(df_Max)
}

PostVar_input<-function(liste){
  df_Var=data.frame()
  for (i in 1:length(liste)){
    temp=lapply(liste[[i]],function(x) {filter(x$param_QTL,True_Cat=="High")$Effet^2*filter(x$param_QTL,True_Cat=="High")$frequence*(1-filter(x$param_QTL,True_Cat=="High")$frequence)*2})
    temp=data.frame(VarEffect=unlist(temp),Variance=liste[[i]][[1]]$info_QTL$Var[1])
    df_Var=rbind(temp,df_Var)
  }
  return(df_Var)
}

df_PostVar_comp<-function(liste,liste_wo,var,nQTL=5,var_snp=80,k=1){
  df_true=data.frame(VarSum=c(var$Var*nQTL,var_snp-var$Var*nQTL),SNP=rep(c("High","Low"),each=length(var$Var)),VarInd=rep(var$Var,2))
  df_true=cbind(df_true,Data=rep("True",nrow(df_true)))
  df_PostVar=df_true
  for (i in 1:length(liste)){
    temp_high=mean(sapply(liste[[i]],function(x) {sum(filter(x$param_app,True_Cat=="High")$post_var)}))*k
    temp_area=mean(sapply(liste[[i]],function(x) {sum(filter(x$param_app,True_Cat=="Area")$post_var)}))*k
    temp_low=mean(sapply(liste[[i]],function(x) {sum(filter(x$param_app,True_Cat=="Low")$post_var)}))*k
    temp_null=mean(sapply(liste[[i]],function(x) {sum(filter(x$param_app,True_Cat=="Null")$post_var)}))*k
    df_temp=data.frame(VarSum=c(temp_high,temp_area,temp_low,temp_null),SNP=c("High","Area","Low","Null"),VarInd=rep(round(liste[[i]][[1]]$info_QTL$Var[[1]],4),digits=3))
    df_temp=cbind(df_temp,Data=rep("With",nrow(df_temp)))
    df_PostVar=rbind(df_PostVar,df_temp)
    rm(df_temp)
  }
  for (i in 1:length(liste_wo)){
    temp_high=mean(sapply(liste_wo[[i]],function(x) {sum(filter(x$param_app,True_Cat=="High")$post_var)}))*k
    temp_area=mean(sapply(liste_wo[[i]],function(x) {sum(filter(x$param_app,True_Cat=="Area")$post_var)}))*k
    temp_low=mean(sapply(liste_wo[[i]],function(x) {sum(filter(x$param_app,True_Cat=="Low")$post_var)}))*k
    temp_null=mean(sapply(liste_wo[[i]],function(x) {sum(filter(x$param_app,True_Cat=="Null")$post_var)}))*k
    df_temp=data.frame(VarSum=c(temp_high,temp_area,temp_low,temp_null),SNP=c("High","Area","Low","Null"),VarInd=rep(round(liste_wo[[i]][[1]]$info_QTL$Var[[1]],4),digits=3))
    df_temp=cbind(df_temp,Data=rep("WO",nrow(df_temp)))
    df_PostVar=rbind(df_PostVar,df_temp)
    rm(df_temp)
  }
  return(df_PostVar)
}

df_Vi_comp<-function(liste,var,nQTL=5,var_snp=80){
  df_true=data.frame(Vi=c(var$Var*nQTL,var_snp-var$Var*nQTL),SNP=rep(c("High","Low"),each=length(var$Var)),VarInd=rep(var$Var,2))
  df_true=cbind(df_true,Data=rep("True",nrow(df_true)))
  df_Vi=df_true
  for (i in 1:length(liste)){
    temp_high=mean(sapply(liste[[i]],function(x) {sum(filter(x$param_app,True_Cat=="High")$V_i)}))
    temp_area=mean(sapply(liste[[i]],function(x) {sum(filter(x$param_app,True_Cat=="Area")$V_i)}))
    temp_low=mean(sapply(liste[[i]],function(x) {sum(filter(x$param_app,True_Cat=="Low")$V_i)}))
    temp_null=mean(sapply(liste[[i]],function(x) {sum(filter(x$param_app,True_Cat=="Null")$V_i)}))
    df_temp=data.frame(Vi=c(temp_high,temp_area,temp_low,temp_null),SNP=c("High","Area","Low","Null"),VarInd=rep(round(liste[[i]][[1]]$info_QTL$Var[[1]],4),digits=3))
    df_temp=cbind(df_temp,Data=rep("With",nrow(df_temp)))
    df_Vi=rbind(df_Vi,df_temp)
    rm(df_temp)
  }
  return(df_Vi)
}


matrix_QTL_detected<-function(liste){
  df_QTL=data.frame(matrix(0,ncol=length(liste),50),stringsAsFactors = F)
  colnames(df_QTL)=names(liste)
  for (i in 1:length(liste)){
    if (i==1){
      rownames(df_QTL)=unlist(lapply(liste[[i]],function(x) {unique(x$param_QTL$Ref)}))
    }
    temp_QTL=unlist(lapply(liste[[i]],function(x) {unique(filter(x$param_QTL,MAP_nsnp=="Nsnp")$Ref)}))
    df_QTL[as.character(temp_QTL),i]=1
    rm(temp_QTL)
  }
  return(df_QTL)
}

matrix_QTL_integration<-function(liste){
  df_QTL=data.frame(matrix(0,ncol=length(liste),50),stringsAsFactors = F)
  rownames(df_QTL)=unlist(lapply(liste[[1]],function(x) {unique(x$param_QTL$Ref)}))
  colnames(df_QTL)=names(liste)
  for (i in 1:length(liste)){
    df_QTL[,i]=do.call(rbind,lapply(liste[[i]],function(x) {aggregate(1-x$param_QTL$PIP1,list(x$param_QTL$Ref),max)[unique(x$param_QTL$Ref),]}))$x
  }
  return(df_QTL)
}



max_LD<-function(liste){
  df_max_LD=do.call(rbind,lapply(liste,function(x) {aggregate(x$R2,list(x$SNP_A),max)}))
  colnames(df_max_LD)=c("QTL","Dprime")
  return(df_max_LD)
}

df_cor_mean<-function(liste_res){
  temp=sapply(liste_res,function(x) {rowMeans(x)["Cor_Validation"]})
  names(temp)=c()
  df_cor_mean=data.frame(Correlation=temp,Scenario=names(liste_res))
  return(df_cor_mean)
}

df_cor<-function(liste_res){
  df_cor_res=data.frame()
  for (i in 1:length(liste_res)){
    temp=data.frame(Correlation=t(liste_res[[i]][2,]),Scenario=names(liste_res)[i])
    colnames(temp)[1]="Correlation"
    df_cor_res=rbind(df_cor_res,temp)
  }
  return(df_cor_res)
}

df_model_va<-function(liste){
  temp=gather(data.frame(sapply(liste,function(y) {sapply(y,function(x) {x$model_app$Va})})),key=Scenario,value=Variance)
  temp=cbind(temp,Variable="Va")
  temp2=gather(data.frame(sapply(liste,function(y) {sapply(y,function(x) {sum(x$param_app$V_i)})})),key=Scenario,value=Variance)
  temp2=cbind(temp2,Variable="Sum_Vi")
  temp3=gather(data.frame(sapply(liste,function(y) {sapply(y,function(x) {sum(x$param_app$variance)})})),key=Scenario,value=Variance)
  temp3=cbind(temp3,Variable="Var_u")
  df_va=rbind(temp,temp2,temp3)
  return(df_va)
}

df_top_Vi_out<-function(out,n=10){
  top=head(arrange(out$param_app, desc(V_i)),n)$Position
  QTL=unique(filter(out$param_QTL,Position %in% head(arrange(out$param_app, desc(V_i)),n)$Position)$Ref)
  return(QTL)
}

df_top_ViMAF_out<-function(out,n=10){
  top=head(arrange(out$param_app, desc(V_i)),n)$Position
  QTL=unique(filter(out$param_QTL,Position %in% head(arrange(out$param_app, desc(VarMAF_Part)),n)$Position)$Ref)
  return(QTL)
}

df_top_Ui_out<-function(out,n=10){
  top=head(arrange(out$param_app, desc(abs(beta))),n)$Position
  QTL=unique(filter(out$param_QTL,Position %in% head(arrange(out$param_app, desc(abs(beta))),n)$Position)$Ref)
  return(QTL)
}

df_top_SumVar_out<-function(out,n=150){
  top=head(arrange(out$param_app, desc(abs(SumVar))),n)$Position
  QTL=unique(filter(out$param_QTL,Position %in% head(arrange(out$param_app, desc(abs(beta))),n)$Position)$Ref)
  return(QTL)
}

df_top_Vi<-function(liste,n=10){
  df_QTL=data.frame(matrix(0,ncol=length(liste),50),stringsAsFactors = F)
  colnames(df_QTL)=names(liste)
  for (i in 1:length(liste)){
    if (i==1){
      rownames(df_QTL)=unlist(lapply(liste[[i]],function(x) {unique(x$param_QTL$Ref)}))
    }
    temp_QTL=unlist(lapply(liste[[i]],function(x) {df_top_Vi_out(x,n)}))
    df_QTL[as.character(temp_QTL),i]=1
    rm(temp_QTL)
  }
  return(df_QTL)
}



df_top_Ui<-function(liste,n=10){
  df_QTL=data.frame(matrix(0,ncol=length(liste),50),stringsAsFactors = F)
  colnames(df_QTL)=names(liste)
  for (i in 1:length(liste)){
    if (i==1){
      rownames(df_QTL)=unlist(lapply(liste[[i]],function(x) {unique(x$param_QTL$Ref)}))
    }
    temp_QTL=unlist(lapply(liste[[i]],function(x) {df_top_Ui_out(x,n)}))
    df_QTL[as.character(temp_QTL),i]=1
    rm(temp_QTL)
  }
  return(df_QTL)
}

df_top_SumVar<-function(liste,n=150){
  df_QTL=data.frame(matrix(0,ncol=length(liste),50),stringsAsFactors = F)
  colnames(df_QTL)=names(liste)
  for (i in 1:length(liste)){
    if (i==1){
      rownames(df_QTL)=unlist(lapply(liste[[i]],function(x) {unique(x$param_QTL$Ref)}))
    }
    temp_QTL=unlist(lapply(liste[[i]],function(x) {df_top_SumVar_out(x,n)}))
    df_QTL[as.character(temp_QTL),i]=1
    rm(temp_QTL)
  }
  return(df_QTL)
}



df_comp_detected<-function(liste,nQTL=50){
  df_QTL=data.frame()
  for (i in 1:length(liste)){
    temp=data.frame(Detected=rep(0,nQTL),Method=rep("Both",nQTL))
    rownames(temp)=unlist(lapply(liste[[i]],function(x) {unique(x$param_QTL$Ref)}))
    temp1=unlist(lapply(liste[[i]],function(x) {df_top_Vi_out(x)}))
    temp2=unlist(lapply(liste[[i]],function(x) {unique(filter(x$param_QTL,MAP_nsnp=="Nsnp")$Ref)}))
    joint_QTL=intersect(temp1,temp2)
    temp1=setdiff(temp1,joint_QTL)
    temp2=setdiff(temp2,joint_QTL)
    
    temp$Method=as.character(temp$Method)
    
    temp[joint_QTL,"Detected"]=1
    temp[temp1,"Detected"]=1
    temp[temp2,"Detected"]=1
    
    temp[temp1,"Method"]="Vi"
    temp[temp2,"Method"]="MAP"
    
    temp3=cbind(temp,Scenario=names(liste)[[i]])
    df_QTL=rbind(df_QTL,temp3)
    rm(temp1)
    rm(temp2)
    rm(temp3)
  }
  return(df_QTL)
}

df_comp_detected2<-function(liste,nQTL=50){
  df_QTL=data.frame()
  for (i in 1:length(liste)){
    temp=data.frame(Detected=rep(0,nQTL),Method=rep("Both",nQTL))
    rownames(temp)=unlist(lapply(liste[[i]],function(x) {unique(x$param_QTL$Ref)}))
    temp1=unlist(lapply(liste[[i]],function(x) {df_top_Vi_out(x)}))
    temp2=unlist(lapply(liste[[i]],function(x) {df_top_SumVar_out(x,100)}))
    joint_QTL=intersect(temp1,temp2)
    temp1=setdiff(temp1,joint_QTL)
    temp2=setdiff(temp2,joint_QTL)
    
    temp$Method=as.character(temp$Method)
    
    temp[joint_QTL,"Detected"]=1
    temp[temp1,"Detected"]=1
    temp[temp2,"Detected"]=1
    
    temp[temp1,"Method"]="Vi"
    temp[temp2,"Method"]="SumVar"
    
    temp3=cbind(temp,Scenario=names(liste)[[i]])
    df_QTL=rbind(df_QTL,temp3)
    rm(temp1)
    rm(temp2)
    rm(temp3)
  }
  return(df_QTL)
}

df_comp_detected_v2<-function(liste,nQTL=50){
  df_QTL=data.frame()
  for (i in 1:length(liste)){
    temp=data.frame(Detected=rep(0,nQTL),Method=rep("MAP+Vi+SumVar",nQTL))
    rownames(temp)=unlist(lapply(liste[[i]],function(x) {unique(x$param_QTL$Ref)}))
 
    temp1=unlist(lapply(liste[[i]],function(x) {unique(filter(x$param_QTL,MAP_nsnp=="Nsnp")$Ref)}))
    temp2=unlist(lapply(liste[[i]],function(x) {df_top_Vi_out(x)}))
    temp3=unlist(lapply(liste[[i]],function(x) {df_top_SumVar_out(x)}))

    temp1=as.character(temp1)
    temp2=as.character(temp2)
    temp3=as.character(temp3)
    
    diff1=setdiff(temp2,temp1)
    diff2=setdiff(temp3,temp2)
    
    temp$Method=as.character(temp$Method)
    
    temp[temp3,"Detected"]=1
    
    temp[diff1,"Method"]="Vi+SumVar"
    temp[diff2,"Method"]="SumVar"
    
    temp3=cbind(temp,Scenario=names(liste)[[i]])
    df_QTL=rbind(df_QTL,temp3)
    rm(temp1)
    rm(temp2)
    rm(temp3)
  }
  return(df_QTL)
}


df_MAP_summary<-function(liste){
  df_MAP=data.frame()
  for (i in 1:length(liste)){
    temp_high=sum(unlist(lapply(liste[[i]],function(y) {length(which(y$param_app$MAP=="High"))})))
    temp_medium=sum(unlist(lapply(liste[[i]],function(y) {length(which(y$param_app$MAP=="Medium"))})))
    temp_low=sum(unlist(lapply(liste[[i]],function(y) {length(which(y$param_app$MAP=="Low"))})))
    df_temp=data.frame(MAP=c(temp_high,temp_medium,temp_low),Classe=c("High","Medium",'Low'))
    df_temp=cbind(df_temp,Scenario=names(liste)[i])
    df_MAP=rbind(df_MAP,df_temp)
  }
  return(df_MAP)
}

df_QTL_summary<-function(liste){
  df_MAP=data.frame()
  for (i in 1:length(liste)){
    temp_high=sum(unlist(lapply(liste[[i]],function(y) {length(unique((filter(y$param_QTL,MAP=="High"))$Ref))})))
    temp_medium=sum(unlist(lapply(liste[[i]],function(y) {length(unique((filter(y$param_QTL,MAP=="Medium"))$Ref))})))
    temp_low=sum(unlist(lapply(liste[[i]],function(y) {length(unique((filter(y$param_QTL,MAP=="Low"))$Ref))})))
    df_temp=data.frame(MAP=c(temp_high,temp_medium,temp_low),Classe=c("High","Medium",'Low'))
    df_temp=cbind(df_temp,Scenario=names(liste)[i])
    df_MAP=rbind(df_MAP,df_temp)
  }
  return(df_MAP)
}

df_Var_MAF<-function(liste){
  for (i in 1:length(liste)){
    for (j in 1:length(liste[[i]])){
      liste[[i]][[j]]$param_app=cbind(liste[[i]][[j]]$param_app,VarMAF=liste[[i]][[j]]$param_app$V_i*2*(1-liste[[i]][[j]]$param_app$frequence)*liste[[i]][[j]]$param_app$frequence)
      liste[[i]][[j]]$param_QTL=cbind(liste[[i]][[j]]$param_QTL,VarMAF=liste[[i]][[j]]$param_QTL$V_i*2*(1-liste[[i]][[j]]$param_QTL$frequence)*liste[[i]][[j]]$param_QTL$frequence)
    }
  }
  return(liste)
}


df_PartVar_MAF<-function(liste){
  for (i in 1:length(liste)){
    for (j in 1:length(liste[[i]])){
      liste[[i]][[j]]$param_app=cbind(liste[[i]][[j]]$param_app,VarMAF_Part=liste[[i]][[j]]$param_app$VarMAF/sum(liste[[i]][[j]]$param_app$VarMAF))
      liste[[i]][[j]]$param_QTL=cbind(liste[[i]][[j]]$param_QTL,VarMAF_Part=liste[[i]][[j]]$param_QTL$VarMAF/sum(liste[[i]][[j]]$param_QTL$VarMAF))
      }
  }
  return(liste)
}

df_Vi=function(liste){
  df_Var=data.frame()
  for (i in 1:length(liste)){
    temp=unlist(lapply(liste[[i]],function(x) {aggregate(x$param_QTL$V_i,list(x$param_QTL$Ref),max)$x}))
    temp2=data.frame(Vi=temp,Scenario=names(liste)[i])
    df_Var=rbind(df_Var,temp2)
  }
  return(df_Var)
}

Classif_SumVar<-function(out){
  classif_SumVar=rep(0,nrow(out$param_app))
  classif_SumVar[which(out$param_app$SumVar<10^(-4))]="Null"
  classif_SumVar[which(out$param_app$SumVar>10^(-4) & out$param_app$SumVar<10^(-3))]="Low"
  classif_SumVar[which(out$param_app$SumVar>10^(-3) & out$param_app$SumVar<10^(-2))]="Medium"
  classif_SumVar[which(out$param_app$SumVar>10^(-2))]="High"
  classif_SumVar2=rep(0,nrow(out$param_QTL))
  classif_SumVar2[which(out$param_QTL$SumVar<10^(-4))]="Null"
  classif_SumVar2[which(out$param_QTL$SumVar>10^(-4) & out$param_QTL$SumVar<10^(-3))]="Low"
  classif_SumVar2[which(out$param_QTL$SumVar>10^(-3) & out$param_QTL$SumVar<10^(-2))]="Medium"
  classif_SumVar2[which(out$param_QTL$SumVar>10^(-2))]="High"
  return(list(classif_SumVar,classif_SumVar2))
}

Classif_liste<-function(liste){
  for (i in 1:length(liste)){
    for (j in 1:length(liste[[i]])){
      res=Classif_SumVar(liste[[i]][[j]])
      liste[[i]][[j]]$param_app=cbind(liste[[i]][[j]]$param_app,Classif_SumVar=res[[1]])
      liste[[i]][[j]]$param_QTL=cbind(liste[[i]][[j]]$param_QTL,Classif_SumVar=res[[2]])
    }
  }
  return(liste)
}
