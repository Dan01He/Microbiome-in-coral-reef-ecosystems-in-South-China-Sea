 library(picante)
 library(iCAMP)
Bphy=read.tree('final.HEBAC_full.tre')
Fphy=read.tree('final.HEFUN_full.tre')

niche_names=c('sedi','wt')
K_names=c('b','e')
Region_names=c('Hainan','ZhongXisha','Nansha')

 for (i in niche_names) {
 for (j in K_names) {
 if (i=='sedi') {envs=get(paste0('envs_',i,'_',j)); env=envs[c(7:11)]; treat=envs[c('Space','Region2','Field')]}
    else { envs=get(paste0('envs_',i,'_',j)); env=envs[c(7:16)]; treat=envs[c('Space','Region2','Field')]}

save.wd=paste0('E:\\HE\\iCAMP\\',i,'_',j)
 if(!dir.exists(save.wd)){dir.create(save.wd)} 
setwd(save.wd)
 if(j=='b') {tree=get(paste0(toupper(j),'phy')); 
	comm=BAC_COM[match(rownames(envs),rownames(BAC_COM)),]; comm=na.omit(comm); comm=comm[,colSums(comm)>1]
	clas=BAC.tax} else {tree=get(paste0(toupper(j),'phy')); 
	comm=microEUK_COM[match(rownames(envs),rownames(microEUK_COM)),]; comm=na.omit(comm); comm=comm[,colSums(comm)>1]
	clas=microEUK.tax
	}
	
	sampid.check=match.name(rn.list=list(comm=comm,treat=treat,env=env))
	treat=sampid.check$treat
	comm=sampid.check$comm
	comm=comm[,colSums(comm)>0,drop=FALSE] 
	env=sampid.check$env 
	
	spid.check=match.name(cn.list=list(comm=comm),rn.list=list(clas=clas),tree.list=list(tree=tree))
	comm=spid.check$comm
	clas=spid.check$clas
	tree=spid.check$tree
 save(comm,tree,clas,treat,env, file=paste0(i,'_',j,'.icamp.RData'))
 
 for (k in Region_names) {
 save.wd=paste0('E:\\HE\\iCAMP\\',i,'_',j,'\\',k)
 if(!dir.exists(save.wd)){dir.create(save.wd)} 
 setwd(save.wd)
  comm=as.data.frame(comm)
  comm2=comm[treat$Region2==k,] # use comm2 instead of comm, for the latter will be re-used for the next k level, and cause error!
  comm2=comm2[,colSums(comm2)>0] 
  env2=env[treat$Region2==k,]
  treat2=subset(treat, Region2==k)
 
	spid.check=match.name(cn.list=list(comm=comm2),rn.list=list(clas=clas),tree.list=list(tree=tree))
	clas2=spid.check$clas
	tree2=spid.check$tree
 save(comm2,tree2,clas2,treat2,env2, file=paste0(i,'_',j,'_',k,'.Foricamp.RData'))
}

}
}


