library(iCAMP)
library(picante)
library(ape)

t0=Sys.time() # to calculate time cost

# 1 # set folder paths and file names, please change according to the folder paths and file names in your computer.
# the folder saving the input files
 #wd=getwd()
 # wd=paste0(curr.wd,'/','Example')
 # if(!dir.exists(wd)){dir.create(wd)}
setwd(wd)
# the folder to save the output. please change to a new folder even if you are just testing the example data.
 save.wd=paste0(wd,'/','Outputs4')
 if(!dir.exists(save.wd)){dir.create(save.wd)}


# 2 # key parameter setting
pfs=unlist(strsplit(wd,split='/'))
pfs
prefix=paste0(pfs[4], '_', pfs[5])  ### prefix of the output file names. usually use a project ID.
prefix

rand.time=1000  # randomization time, 1000 is usually enough. For example test, you may set as 100 or less to save time.
nworker=2 # nworker is thread number for parallel computing, which depends on the CPU core number of your computer.
memory.G=1 # to set the memory size as you need (but should be less than the available space in your hard disk), so that calculation of large tree will not be limited by physical memory. unit is Gb.

# 3 # load R packages and data
if(exists('comm2')) {
comm=comm2
tree=tree2
clas=clas2
treat=treat2
env=env2
}

comm=sqrt(comm)
# 4 # match sample IDs in OTU table and treatment information table
sampid.check=match.name(rn.list=list(comm=comm,treat=treat,env=env))
# sampid.check=match.name(rn.list=list(comm=comm,treat=treat)) # if you do not have env.file
# for the example data, the output should be "All match very well".
# for your data files, if you have not matched their IDs, the unmatched samples will be removed.
treat=sampid.check$treat
comm=sampid.check$comm
comm=comm[,colSums(comm)>0,drop=FALSE] # if some unmatched samples were removed, some OTUs may become ghosts, then you may use this line to remove them if necessary.
env=sampid.check$env # skip this if you do not have env.file

# 5 # match OTU IDs in OTU table and tree file
spid.check=match.name(cn.list=list(comm=comm),rn.list=list(clas=clas),tree.list=list(tree=tree))
# for the example data, the output should be "All match very well".
# for your data files, if you have not matched the IDs before, the unmatched OTUs will be removed.
comm=spid.check$comm
clas=spid.check$clas
tree=spid.check$tree

# 6 # calculate pairwise phylogenetic distance matrix.
# since microbial community data usually has a large number of species (OTUs or ASVs), we use "big.matrix" in R package "bigmemory" to handle the large phylogenetic distance matrix. 

setwd(save.wd)
if(!file.exists("pd.desc")) 
{
  pd.big=iCAMP::pdist.big(tree = tree, wd=save.wd, nworker = nworker, memory.G = memory.G)

}else{
  # if you already calculated the phylogenetic distance matrix in a previous run
  pd.big=list()
  pd.big$tip.label=read.csv(paste0(save.wd,"/pd.taxon.name.csv"),row.names = 1,stringsAsFactors = FALSE)[,1]
  pd.big$pd.wd=save.wd
  pd.big$pd.file="pd.desc"
  pd.big$pd.name.file="pd.taxon.name.csv"
}


####################
# 7 # Phylogenetic Signal Test
####################
# 7.1 # env factor transformation and standardization
################
envin=env
env=envin[,1:ncol(envin)]
lnx<-function(v)
{
  minv=min(v)
  if(minv<=0){minv=min(density(v)$x);v=v-minv}
  log(v)
}
envln=sapply(1:ncol(env),function(i){lnx(env[,i])})
colnames(envln)=colnames(env)
rownames(envln)=rownames(env)

if(sum(colnames(envln)%in%'pH')>0) envln[,which(colnames(envln)=="pH")]=env[,"pH"]



################
# 7.2 # niche difference
################
envi="TOC" # test representative factors
setwd(save.wd)
library(iCAMP)
envim=envln[,envi,drop=FALSE]
nichedi=iCAMP::dniche(env = envim,comm = comm,method = "niche.value",
                      nworker = nworker, out.dist=FALSE,bigmemo=F,
                      nd.wd=save.wd)



################
# 7.3 # root the tree
################
if(!ape::is.rooted(tree))
{
  tree.rt=iCAMP::midpoint.root.big(tree = tree, pd.desc = pd.big$pd.file,
                                   pd.spname = pd.big$tip.label,pd.wd = pd.big$pd.wd,
                                   nworker = nworker)
  tree=tree.rt$tree
}



####################
# 8 # Other approach: QPEN (quantifying community assembly processes based on entire-community null model analysis)
####################
8.1 # QPEN calculation
qpout=iCAMP::qpen.cm(comm=comm,pd=pd.big$pd.file,pd.big.wd=pd.big$pd.wd,
                     pd.big.spname=pd.big$tip.label,ab.weight=TRUE,
                     rand.time=rand.time, nworker=nworker,project=prefix,
                     wd=save.wd, save.bNTIRC=TRUE,memory.G=1)
qpout.uw=iCAMP::qpen.cm(comm=comm,pd=pd.big$pd.file,pd.big.wd=pd.big$pd.wd,
                     pd.big.spname=pd.big$tip.label,ab.weight=FALSE,
                     rand.time=rand.time, nworker=nworker,project=prefix,
                     wd=save.wd, save.bNTIRC=TRUE,memory.G=1)
if(dim(comm)[2]>2000) {
qpout.tpuw=iCAMP::qpen.cm(comm=comm[,order(colSums(comm),decreasing=T)[1:2000]],pd=pd.big$pd.file,pd.big.wd=pd.big$pd.wd,
                     pd.big.spname=pd.big$tip.label,ab.weight=FALSE,
                     rand.time=rand.time, nworker=nworker,project=prefix,
                     wd=save.wd, save.bNTIRC=TRUE,memory.G=1)

qpout.tp=iCAMP::qpen.cm(comm=comm[,order(colSums(comm),decreasing=T)[1:2000]],pd=pd.big$pd.file,pd.big.wd=pd.big$pd.wd,
                     pd.big.spname=pd.big$tip.label,ab.weight=TRUE,
                     rand.time=rand.time, nworker=nworker,project=prefix,
                     wd=save.wd, save.bNTIRC=TRUE,memory.G=1) } else if (dim(comm)[2]>1000)
 {
qpout.tpuw=iCAMP::qpen.cm(comm=comm[,order(colSums(comm),decreasing=T)[1:1000]],pd=pd.big$pd.file,pd.big.wd=pd.big$pd.wd,
                     pd.big.spname=pd.big$tip.label,ab.weight=FALSE,
                     rand.time=rand.time, nworker=nworker,project=prefix,
                     wd=save.wd, save.bNTIRC=TRUE,memory.G=1)

qpout.tp=iCAMP::qpen.cm(comm=comm[,order(colSums(comm),decreasing=T)[1:1000]],pd=pd.big$pd.file,pd.big.wd=pd.big$pd.wd,
                     pd.big.spname=pd.big$tip.label,ab.weight=TRUE,
                     rand.time=rand.time, nworker=nworker,project=prefix,
                     wd=save.wd, save.bNTIRC=TRUE,memory.G=1) } else
 {
qpout.tpuw=qpout.uw
qpout.tp=qpout
}

# 8.2 # significance test
qptest=iCAMP::qpen.test(qpen.result = qpout,treat = treat,rand.time = rand.time,
                        between.group = TRUE,out.detail=TRUE,silent=FALSE)
qptest.uw=iCAMP::qpen.test(qpen.result = qpout.uw,treat = treat,rand.time = rand.time,
                        between.group = TRUE,out.detail=TRUE,silent=FALSE)
qptest.tp=iCAMP::qpen.test(qpen.result = qpout.tp,treat = treat,rand.time = rand.time,
                        between.group = TRUE,out.detail=TRUE,silent=FALSE)
qptest.tpuw=iCAMP::qpen.test(qpen.result = qpout.tpuw,treat = treat,rand.time = rand.time,
                        between.group = TRUE,out.detail=TRUE,silent=FALSE)

save(qptest, qptest.uw, qptest.tp, qptest.tpuw, file = paste0(prefix,".cM.QPEN.boot.detail.rda"))


(t1=format(Sys.time()-t0)) 
# End #

Noprefix=c('aa','bb','prefix','wd','wds','region','Noprefix')
for( i in ls()) {
if ((sum(Noprefix %in% i)==0)&(!grepl('Hainan',i)&!grepl('ZhongXisha',i)&!grepl('Nansha',i))){
assign(paste0(prefix,'_',i),get(i))
rm(list=i)
}
}

save.image(file=paste0(prefix,'_','icamp_qpprefix.RData'))
rm(prefix)

#Real End #

