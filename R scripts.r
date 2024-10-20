require(ggplot2)
require(vegan)
require(reshape2)

### Fig.1
 
 envs_wt=read.csv(file='Env.seawater.csv',row.names=1,head=T)
 envs_sedi=read.csv(file='Env.sediment.csv',row.names=1,head=T)
 s.factor=read.csv(file='Metadata.csv.csv',row.names=1,head=T)
 
 envs_wt.m=melt(envs_wt, measure.vars=c(7:16))
 ggbarplot(data=envs_wt.m, x = 'Region2', y = "value",fill='Space', add = "mean_se", position = position_dodge( .9),facet.by='variable',scale='free')+theme(axis.text.x = element_text(angle = 15, hjust = 1))
 ggboxplot(data=envs_wt.m, x = 'Region2', y = "value",fill='Space',facet.by='variable',scale='free')

 envs_sedi.m=melt(envs_sedi, measure.vars=c(7:11))
 ggbarplot(data=envs_sedi.m, x = 'Region2', y = "value",fill='Space', add = "mean_se", position = position_dodge( .9),facet.by='variable',scale='free')+theme(axis.text.x = element_text(angle = 15, hjust = 1))
 ggboxplot(data=envs_sedi.m, x = 'Region2', y = "value",fill='Space',facet.by='variable',scale='free')


### Fig.2 
# Fig. 2a, b
microEUKVGalpha=data.frame(Shannon=diversity(microEUK_COM),InvSimp=diversity(microEUK_COM, "simpson"),Richness=specnumber(microEUK_COM),Pielou_Evenness=diversity(microEUK_COM)/log(specnumber(microEUK_COM)))
BACVGalpha=data.frame(Shannon=diversity(BAC_COM),InvSimp=diversity(BAC_COM, "simpson"),Richness=specnumber(BAC_COM),Pielou_Evenness=diversity(BAC_COM)/log(specnumber(BAC_COM)))

s.factor_b=na.omit(s.factor[rownames(BAC_COM),])
s.factor_e=na.omit(s.factor[rownames(microEUK_COM),])

envs_sedi_b=na.omit(envs_sedi[rownames(BAC_COM),])
envs_sedi_e=na.omit(envs_sedi[rownames(microEUK_COM),])
envs_wt_b=na.omit(envs_wt[rownames(BAC_COM),])
envs_wt_e=na.omit(envs_wt[rownames(microEUK_COM),])

BACVGalphaf=cbind(BACVGalpha,s.factor_b)
microEUKVGalphaf=cbind(microEUKVGalpha,s.factor_e)

library(reshape2)
BAC.alphaf.m=melt(BACVGalphaf,id=c(5:11))
microEUK.alphaf.m=melt(microEUKVGalphaf,id=c(5:11))
require(ggplot2)
ggplot(data=subset(BAC.alphaf.m,variable=='Richness'), aes(y=value,x=Region2))+geom_boxplot(aes(fill=Niche),width=0.5)+facet_wrap(~variable,scales="free_y")+theme_bw()+theme(line=element_blank())
ggplot(data=subset(microEUK.alphaf.m,variable=='Richness'), aes(y=value,x=Region2))+geom_boxplot(aes(fill=Niche),width=0.5)+facet_wrap(~variable,scales="free_y")+theme_bw()+theme(line=element_blank())

# Fig. 2c,d
bCOMTAX=read.csv('bCOMTAX.csv', row.names=1, head=T)
eCOMTAX1=read.csv('eCOMTAX.csv', row.names=1, head=T)
mid=base::intersect(colnames(bCOMTAX),colnames(eCOMTAX))
metaCOMTAX=rbind(bCOMTAX[,mid],eCOMTAX[,mid])
metaCOMTAX[metaCOMTAX$Genus=='Vibrio',]-> VbCOMTAX

summaryBy(.~Genus,data=VbCOMTAX[c( 1:(dim(VbCOMTAX)[2]-7),(dim(VbCOMTAX)[2]-1) )],FUN=sum)->vbrio.sum; rownames(vbrio.sum)=vbrio.sum[,1]; vbrio.sum=as.data.frame(t(vbrio.sum[-1]));
 rownames(vbrio.sum)=gsub('.sum','',rownames(vbrio.sum));
vbrio.taxonra=data.frame(Vibrio=round(vbrio.sum/34258*100,5), s.factor_b)
vbrio.mlt=melt(vbrio.taxonra,id=c((dim(vbrio.taxonra)[2]-7):(dim(vbrio.taxonra)[2])))
require(ggplot2)
ggplot(vbrio.mlt, aes(x=Region2, y=Vibrio)) +geom_boxplot(aes(fill=Niche))+facet_wrap(~Niche)+theme_bw() # may change



### Fig.3 
# Fig.3a 
COM=BAC_COM
bray=vegdist(sqrt(COM)); mds.bray=metaMDS(bray)
scr.bMDS=as.data.frame(cbind(scores(mds.bray),s.factor_b))
ggplot(data = scr.bMDS, mapping = aes(x = NMDS1, y = NMDS2))+geom_point(aes(colour = Region2, shape= Niche),size=3)+theme_bw()+
annotate("text",x=max(scr.bMDS$NMDS1)*0.89, y=max(scr.bMDS$NMDS2)*0.99, size=2.5,alpha=0.9, label=paste("Stress:",round(mds.bray$stress,2)),fontface='italic')

# Fig.3c
COM=microEUK_COM
bray=vegdist(sqrt(COM)); mds.bray=metaMDS(bray)
scr.bMDS=as.data.frame(cbind(scores(mds.bray),s.factor_e))
ggplot(data = scr.bMDS, mapping = aes(x = NMDS1, y = NMDS2))+geom_point(aes(colour = Region2, shape= Niche),size=3)+theme_bw()+
annotate("text",x=max(scr.bMDS$NMDS1)*0.89, y=max(scr.bMDS$NMDS2)*0.99, size=2.5,alpha=0.9, label=paste("Stress:",round(mds.bray$stress,2)),fontface='italic')+
annotate("text",x=0.5*(max(scr.bMDS$NMDS1)+min(scr.bMDS$NMDS1)),y=1.05*max(scr.bMDS$NMDS2), label="Micro-Eukaryote NMDS",size=4,color="darkblue")

# Fig. 3b, d
Niche1=c('sedi','water')
Region1=c('Hainan','ZhongXisha','Nansha')
for (i in Niche1) {
for (j in Region1) {
mid_b=as.vector(vegdist(sqrt(BAC_COM[s.factor_b$Niche==i&s.factor_b$Region2==j,]))); n=length(mid_b)
bray_b=data.frame(Bray_curtis=mid_b,Region2=rep(j,n),Niche=rep(i,n),Kingdom=rep('Bacteria',n))
mid_f=as.vector(vegdist(sqrt(microEUK_COM[s.factor_e$Niche==i&s.factor_e$Region2==j,]))); n=length(mid_f)
bray_f=data.frame(Bray_curtis=mid_f,Region2=rep(j,n),Niche=rep(i,n),Kingdom=rep('Eukaryotes',n))

assign(paste0('BAC_',i,'_bray_',j),bray_b)
assign(paste0('microEUK_',i,'_bray_',j),bray_f)
if(i=='sedi'&j=='Hainan'){mid_ball=bray_b} else {mid_ball=rbind(mid_ball,bray_b)}
if(i=='sedi'&j=='Hainan'){mid_fall=bray_f} else {mid_fall=rbind(mid_fall,bray_f)}
}
}
BACBray_niche_region2=mid_ball; rm(mid_ball)
microEUKBray_niche_region2=mid_fall; rm(mid_fall)
BACBray_niche_region2$Region2=factor(BACBray_niche_region2$Region2, levels=c('Hainan','ZhongXisha','Nansha'))
microEUKBray_niche_region2$Region2=factor(microEUKBray_niche_region2$Region2, levels=c('Hainan','ZhongXisha','Nansha'))
Bray_niche_region2=rbind(BACBray_niche_region2, microEUKBray_niche_region2)

library(ggpubr)
library(ggplot2)
ggplot(Bray_niche_region2, aes(x=Region2, y=Bray_curtis,fill=Niche)) +geom_bar(stat='summary', fun='mean',position=position_dodge( .9))+stat_summary(fun.data = 'mean_se', geom = "errorbar", colour = "black",width = 0.15,position = position_dodge( .9))+theme_bw()+facet_wrap(~Kingdom)


### Fig. 4
# Fig.4 a-d

# db-RDA analysis
ENVs=envs_sedi_b[c(7:11)]; COM=BAC_COM[match(rownames(ENVs),rownames(BAC_COM)),]; COM=na.omit(COM) 
  # ENVs=envs_sedi_e[c(7:11)]; COM=microEUK_COM[match(rownames(ENVs),rownames(microEUK_COM)),]; COM=na.omit(COM) ## uncomment this for sediment microeukaryotes
   # ENVs=envs_wt_b[c(7:16)]; COM=BAC_COM[match(rownames(ENVs),rownames(BAC_COM)),]; COM=na.omit(COM) ## uncomment this for water bacteria
    # ENVs=envs_wt_e[c(7:16)]; COM=microEUK_COM[match(rownames(ENVs),rownames(microEUK_COM)),]; COM=na.omit(COM) ## uncomment this for water microeukaryotes

COM=COM[,colSums(COM)>1]
dim(COM); COM=sqrt(COM)

capT=capscale(COM~.,data=as.data.frame(scale(ENVs)),dist="bray",add=T)
(cap_adj=(RsquareAdj(capT))$adj.r.squared)

captT=capture.output(capR2step <- ordiR2step(capscale(COM~1, as.data.frame(scale(ENVs)),dist="bray",add=T), scope = formula(capT), R2scope = cap_adj, direction = 'forward', permutations = 9999))

vif.cca(capT)
vif.cca(capR2step)

cap.scaling2 <- summary(capR2step, scaling = 2)
cap.site <- data.frame(cap.scaling2$sites)[1:2]
cap.env <- data.frame(cap.scaling2$biplot)[1:2]
cap.env$group=rownames(cap.env)

bacsedi_cap.site=cap.site; bacsedi_cap.env=cap.env
  # euksedi_cap.site=cap.site;  euksedi_cap.env=cap.env ## uncomment this for sediment microeukaryotes
   # bacwt_cap.site=cap.site;  bacwt_cap.env=cap.env ## uncomment this for water bacteria
    # eukwt_cap.site=cap.site;  eukwt_cap.env=cap.env ## uncomment this for water microeukaryotes


# plot db-RDA
for (i in c('sedi','wt')) {
for (j in c('bac','euk')) {
cap.site=get(paste0(j,i,'_cap.site'))
cap.env=get(paste0(j,i,'_cap.env'))

if (i == 'sedi') { 
 if (j == 'bacteria') { cap.site=cbind(cap.site,envs_sedi_b[rownames(cap.site),c(1:6,14)]); } else
  {cap.site=cbind(cap.site,envs_sedi_e[rownames(cap.site),c(1:6,14)])} 
 } else
 {
 if (j == 'bacteria') { cap.site=cbind(cap.site,envs_wt_b[rownames(cap.site),c(1:6,19)]) } else
  {cap.site=cbind(cap.site,envs_wt_e[rownames(cap.site),c(1:6,19)]) } 
 }

library(ggplot2)
color20<-c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000')

p.cap <- ggplot(data=cap.site, aes(CAP1, CAP2)) +
geom_point(aes(color =Region2, shape=Space),size=2,alpha=0.9) +
scale_color_manual(values = color20) +scale_shape_manual(values = c(16:20))+
theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), legend.title = element_blank(), legend.key = element_rect(fill = 'transparent')) + 
geom_vline(xintercept = 0, color = 'gray', linewidth = 0.5) + 
geom_hline(yintercept = 0, color = 'gray', linewidth = 0.5) +
labs(x=paste('CAP1',round(cap_adj*capT$CCA$eig[1]/sum(capT$CCA$eig)*100,2),'%'),y=paste('CAP2',round(cap_adj*capT$CCA$eig[2]/sum(capT$CCA$eig)*100,2),'%'))+
geom_segment(data = cap.env, aes(x = 0,y = 0, xend = 2*CAP1,yend = 2*CAP2), arrow = arrow(length = unit(0.2, 'cm')), size = 0.3, color = 'blue') +
geom_text(data = cap.env, aes(x=CAP1 * 2.5, y=CAP2 * 2.5, label = group), color = 'blue', size = 3)

assign(paste0(j,i,'_caplot'), p.cap)
pdf(paste0(j,i,'_db-RDA.pdf'),width=5.5,height=4)
print(p.cap)
dev.off()
}
}


# Fig 4. e-h
ENVs=envs_wt_e; ENVs.p=envs_wt_e[c(7:16)]
for (i in c('water','sedi')) {
for (j in c('b','e')) {
if (i=='water') {ENVs=envs_wt_b; ENVs.p=envs_wt_b[c(7:16)]} else
 if (j=='b') {ENVs=envs_sedi_b; ENVs.p=envs_sedi_b[c(7:11)] } else
{ENVs=envs_sedi_e; ENVs.p=envs_sedi_e[c(7:11)]}
ENVs.pln=sapply(1:ncol(ENVs.p),function(i){log2(ENVs.p[,i]+1e-20)})
colnames(ENVs.pln)=colnames(ENVs.p)
rownames(ENVs.pln)=rownames(ENVs.p)
if(sum(colnames(ENVs.pln)%in%'pH')>0) ENVs.pln[,which(colnames(ENVs.pln)=="pH")]=ENVs.p[,"pH"]
if (j=='b') {COM=BAC_COM[match(rownames(ENVs),rownames(BAC_COM)),]; COM=na.omit(COM); COM=sqrt(COM[,colSums(COM)>1])} else
{COM=microEUK_COM[match(rownames(ENVs),rownames(microEUK_COM)),]; COM=na.omit(COM); COM=sqrt(COM[,colSums(COM)>1])}

VP1=varpart(vegdist(COM),  as.data.frame(scale(ENVs.p)), ~PCNM1+PCNM2,~Region2, ~Space,  data=ENVs); 
VP2=varpart(vegdist(COM),  as.data.frame(scale(ENVs.p)), ~PCNM1+PCNM2, ~Region2, data=ENVs); 

assign(paste0('VP1.',i,'_',j),VP1)
assign(paste0('VP2.',i,'_',j),VP2)
pdf(paste0('VP1.',i,'_',j,'.pdf'))
plot( Xnames=c("Envrion.","PCNM",'Region','Space'),VP1)
#plot( Xnames=c("Envrion.","PCNM",'Region'),VP2)
dev.off()
}
}


### Fig. 5
source('H:/R/scripts_functions/HEsubsetcommenv.files.r')

## each niche-region qptest
wds=c( 'E:/HE/iCAMP/wt_b/', 'E:/HE/iCAMP/wt_e/', 'E:/HE/iCAMP/sedi_b/', 'E:/HE/iCAMP/sedi_e/')
region=c('Hainan','ZhongXisha','Nansha')
for (aa in wds) { # use aa to avoid confusing with for loops in the next source file.
for (bb in region) {  # use bb to avoid confusing with for loops in the next source file.
wd=paste0(aa,bb)
setwd(wd)
prefix=paste0(unlist(strsplit(wd,split='/'))[4], '_', unlist(strsplit(wd,split='/'))[5])
load(paste0(prefix,'.Foricamp.RData'))
Noprefix=c('aa','bb','prefix','wd','wds','region','Noprefix')
source('E:/HE/iCAMP/icamp.neat.simp3.r')
}
}
save.image('icamp_qp.All.RData')

files=list.files()
Rdfiles=files[grepl('RData',files)]
for (i in Rdfiles) load(i)

 niche_names=c('sedi','wt')
 K_names=c('b','e')
 Region_names=c('Hainan','ZhongXisha','Nansha')
 qpouts=NULL
 qpouts2=NULL
 qpouts3=NULL
 for (i in niche_names) {
 for (j in K_names) {
 for (k in Region_names) {
 grp=paste0(i,'_',j,'_',k)
 outs=(get(paste0(grp,'_qpout')))$ratio
 outs.tp=(get(paste0(grp,'_qpout.tp')))$ratio
 outs.uw=(get(paste0(grp,'_qpout.uw')))$ratio
 outs.tpuw=(get(paste0(grp,'_qpout.tpuw')))$ratio
 outs$grp=grp
 outs.tp$grp=grp
 outs.uw$grp=grp
  outs.tpuw$grp=grp

 if(i=='sedi'&j=='b'&k=='Hainan') {
 qpouts=outs
 qpouts.tp=outs.tp
 qpouts.uw=outs.uw
 qpouts.tpuw=outs.tpuw
 }   else  {
 qpouts=rbind(qpouts, outs)
 qpouts.tp=rbind(qpouts.tp, outs.tp)
 qpouts.uw=rbind(qpouts.uw, outs.uw)
 qpouts.tpuw=rbind(qpouts.tpuw, outs.tpuw)
 }
 }
 }
 }

for( i in c('qpouts','qpouts.tp','qpouts.uw','qpouts.tpuw')) {
mid=get(i)

mid$Style[grepl('wt',mid$grp)]='water'
mid$Style[grepl('sedi',mid$grp)]='sediment'
mid$Region2[grepl('Hainan',mid$grp)]='Hainan'
mid$Region2[grepl('ZhongXisha',mid$grp)]='ZhongXisha'
mid$Region2[grepl('Nansha',mid$grp)]='Nansha'
mid$Kingdom[grepl('_b',mid$grp)]='Bacteria'
mid$Kingdom[grepl('_e',mid$grp)]='Eukaryota'
mid$Region2=factor(mid$Region2,levels=c('Hainan','ZhongXisha','Nansha'))
assign(i, mid)
} 
# qpout.m=reshape2::melt(qpouts,id=c(6:10))
 qpout.m=reshape2::melt(qpouts.tp,id=c(6:10))
  # qpout.m=reshape2::melt(qpouts.uw,id=c(6:10))
   # qpout.m=reshape2::melt(qpouts.tpuw,id=c(6:10))

qpout.m$Region2=factor(qpout.m$Region2,levels=c('Hainan','ZhongXisha','Nansha'))
qpout.m$Style=factor(qpout.m$Style,levels=c('water','sediment'))
qpout.m$value=as.numeric(qpout.m$value)

# plot relative importance
ggplot(subset(qpout.m), aes(x = Region2, y = value, fill = variable)) +
  eom_bar(width = 1, stat = "identity", position="stack")+ facet_grid(Style~Kingdom)


### Fig. 6

library(pacman)
p_load(igraph)
p_load(RMThreshold)
p_load(WGCNA)
p_load(multtest)
p_load(vegan)
p_load(info.centrality)
allowWGCNAThreads(nThreads =3 )

meta_sampleID=base::intersect(rownames(microEUK_COM),rownames(BAC_COM))
meta_COM=cbind(BAC_COM[meta_sampleID,],microEUK_COM[meta_sampleID,])
META.tax=rbind(BAC.tax,microEUK.tax)

# (1.1) choose BAC COM
ENVs=envs_sedi_b[c(7:11)]; 
 # ENVs=envs_wt_b[c(7:16)];  ## to calculate for water networks, uncomment this line

bCOM=BAC_COM[match(rownames(ENVs),rownames(BAC_COM)),]; bCOM=na.omit(bCOM)
 bCOM=sqrt(bCOM[,order(colSums(bCOM),decreasing=T)[1:1000]])

# (1.2) choose microEUK COM
 ENVs=envs_sedi_e[c(7:11)]
 # ENVs=envs_wt_e[c(7:16)] ## to calculate for water networks, uncomment this line
eCOM=microEUK_COM[match(rownames(ENVs),rownames(microEUK_COM)),]; eCOM=na.omit(eCOM)
eCOM=sqrt(eCOM[,order(colSums(eCOM),decreasing=T)[1:1000]])


# (2) COM ASV correlation
metanet_splID=base::intersect(rownames(bCOM),rownames(eCOM))
netcom=cbind(bCOM[metanet_splID,],eCOM[metanet_splID,]) 
netcom=netcom[,colSums(netcom)>0]
meta_otutax=rbind(BAC.tax,microEUK.tax)
dim(netcom); dim(meta_otutax)

meta_sedi_netcom_Hainan=netcom[s.factor[rownames(netcom),]$Region2=='Hainan',]
meta_sedi_netcom_ZhongXisha=netcom[s.factor[rownames(netcom),]$Region2=='ZhongXisha',]
meta_sedi_netcom_Nansha=netcom[s.factor[rownames(netcom),]$Region2=='Nansha',]

## to calculate for water networks, uncomment these 3 lines
 # meta_wt_netcom_Hainan=netcom[s.factor[rownames(netcom),]$Region2=='Hainan',]
 # meta_wt_netcom_ZhongXisha=netcom[s.factor[rownames(netcom),]$Region2=='ZhongXisha',]
 # meta_wt_netcom_Nansha=netcom[s.factor[rownames(netcom),]$Region2=='Nansha',]


## correlation and adjust p
 (corAndPvalue(netcom, method="pearson"))$p ->rawp
 (corAndPvalue(netcom, method="pearson"))$cor ->corr

col.corr=mat2cols(corr); sum(col.corr$value>0.75);sum(col.corr$value<(-0.75))
col.rawp=mat2cols(rawp)
 is.numeric(col.rawp$value)
 mt.adjp=mt.rawp2adjp(col.rawp$value)   # mt.adjp is a list
 mt.index=mt.adjp$index # index from the smallest p value


# (3) RMThreshold AND NETWORK DEVOLUTION
isSymmetric(corr); corr=na.omit(corr)
 diag(corr) <- 0; 
 corr=abs(corr)  
 rm.matrix.validation(as.matrix(corr))


## random matrix threshold
#res=rm.get.threshold(corr)
 res=rm.get.threshold(corr,interval=c(0.75,0.95),unfold.method = "spline",interactive = F, discard.zeros=T)
 (thre=mean(c(res$sse.chosen, res$p.ks.chosen)))
 # (thre = mean(res$chosen.thresholds))
 cleaned.matrix <- rm.denoise.mat(corr, threshold = thre) ## for wt 0.874;  sedi, 0.866

 
 cleaned.matrix <- rm.discard.zeros(cleaned.matrix) ;diag(cleaned.matrix)=0
 sum(cleaned.matrix!=0)
 range(abs(cleaned.matrix[cleaned.matrix!=0]))
 
 diag_forND=abs(cleaned.matrix) 

 diag(diag_forND)=0
 isSymmetric.matrix(diag_forND)
 write.csv(diag_forND, file="diag_forND.csv")


## use python , paste the following commented lines in python console !
# python -m pip install pandas 
# python
# from pandas import pandas as pd
# import os
# os.chdir('H:/HEwtsedi')
# import ND3 # put ND3.py in the $RWD

# df=pd.read_table("diag_forND.csv",sep=',',header=0,index_col=0) 
# ND_df=ND3.ND(df)
# ND_df.to_csv('ND_out.csv')

ND_out=read.csv('ND_out.csv',head=T,row.names=1)
diag(ND_out)=1
cleaned.matrix2 <- rm.denoise.mat(as.matrix(ND_out), threshold = thre) 
cleaned.matrix2 <- rm.discard.zeros(cleaned.matrix2) ; diag(cleaned.matrix2)=0
 sum(cleaned.matrix2!=0)
 range(abs(cleaned.matrix2[cleaned.matrix2!=0]))

mat2cols(cleaned.matrix2)->col.ND_out; dim(col.ND_out) 
subset(col.ND_out,value!=0)->col.ND_out ; dim(col.ND_out)

edge_pair=paste(col.ND_out[,1],col.ND_out[,2])
ND.edges=col.corr[paste(col.corr[,1],col.corr[,2])%in%edge_pair|
paste(col.corr[,2],col.corr[,1])%in%edge_pair,]
sum(ND.edges$value<0); sum(ND.edges$value>0)

(RMT_ND.net <- graph_from_data_frame(ND.edges[c(1,2,3)],directed=F))

RMT_ND.nettax=meta_otutax[match(as_ids(V(RMT_ND.net)),rownames(meta_otutax)),]; RMT_ND.nettax=na.omit(RMT_ND.nettax)
RMT_ND.nettax=data.frame(ASVid=rownames(RMT_ND.nettax),RMT_ND.nettax)
ND.edges$Relationship=NA
ND.edges$Relationship[ND.edges$value>0]='Positive'
ND.edges$Relationship[ND.edges$value<0]='Negative'

graph_from_data_frame(ND.edges,directed=F,vertices=RMT_ND.nettax)-> RMT_ND.net
table(edge.attributes(RMT_ND.net)$Relationship)
 
RMT_ND.net -> meta_sedi.net
  # RMT_ND.net -> meta_wt.net   ## to calculate for water networks, uncomment this line

base::intersect(as_ids(E(meta_wt.net)),as_ids(E(meta_sedi.net)))

# (4)plot networks in igraph 

taxon32=c('Gammaproteobacteria','Bacteroidota','Planctomycetota','Latescibacterota','Alphaproteobacteria','Bacteria_unclassified','Verrucomicrobiota','Acidobacteriota','Patescibacteria','Actinobacteriota','Desulfobacterota','Spirochaetota','Chloroflexi','Myxococcota','Bdellovibrionota','Minor_bacteria_groups','Eukaryota_unclassified','Dinoflagellata','Diatomea','Ciliophora','Protalveolata','Chlorophyta_ph','Ochrophyta_ph','Peronosporomycetes','Cercozoa','Cryptophyceae_ph','Platyhelminthes','MAST−3','Holozoa_ph','Tunicata','Ascomycota','Minor_eukaryote_groups')


taxon32=factor(taxon32,levels=taxon32)
color_32=c('#C1E8E9','#1F91FF','#FF915F','#B29FFF','#FFCD00','#00E97E','#0DC7C8','#CC90CB','#739C2F','#EEDFA7','#5912F4','#C71585','#CCCCCC','#00293F','#E2F100','#666666','#87CEFA','#0915EC','#A52A2A','#00FF00','#FFA500','#A020F0','#008B00','#FFC0CB','#8B5A00','#8B8B00','#8299EA','#C6FF3F','#FF3BA1','#E62617','#006A66','#2B2B2B')

remove_isolates <- function(g) {
  g.degree <- degree(g)
  g.isolates <- V(g)[which(g.degree < 1)]
  return(delete_vertices(g, g.isolates))
}

## plot igraph for different niche networks in loops
subs.trt=c('meta_wt','meta_sedi')
for (i in subs.trt)
{library(igraph)
g=get(paste0(i,'.net'))
g=remove_isolates(g)
linecol_pal <- c("tomato1","skyblue")
edge_cat=as.numeric(ifelse(edge.attributes(g)$Relationship=='Negative',1,2))

sqrt(colSums(meta_COM[,vertex.attributes(g)$name]))->vertex.attributes(g)$Size ### important for plot
vertex.attributes(g)$taxon=factor(meta_otutax[vertex.attributes(g)$name,]$taxon,levels=levels(taxon32))
vertex.attributes(g)$taxon[is.na(vertex.attributes(g)$taxon)]=ifelse(vertex.attributes(g)$Kingdom[is.na(vertex.attributes(g)$taxon)]=='Bacteria', 'Minor_bacteria_groups', 'Minor_eukaryote_groups')

node_cat=as.numeric(vertex.attributes(g)$taxon)
V(g)$color <-color_32[node_cat] #

assign(paste0(i,'.net'),g)

V(g)$size=plotrix::rescale(c(vertex.attributes(g)$Size, 2, 300 ),c(1,7))[1:vcount(g)] # rescale size
V(g)$frame.color=NA; V(g)$label=NA
E(g)$color= c("tomato1","grey10")[edge_cat]; E(g)$width=0.01; E(g)$curved=T
E(g)$name=as_ids(E(g))
g <- set_graph_attr(g, "layout", layout_with_fr(g))
# g <- set_graph_attr(g, "layout", layout_nicely(g))
edge.attr=as.data.frame(edge.attributes(g))
edge.attr <- tidyr::separate(edge.attr, 'name', c("Source", "Target"),  "\\|", F, F) %>% dplyr::select(c(7,8),everything())  ## for Gephi
edge.attr$intertype=NA 
edge.attr$intertype[grepl('BAC_ASV',edge.attr$Source)&grepl('BAC_ASV',edge.attr$Target)]='Intra_Bacteria'
edge.attr$intertype[grepl('FUN_ASV',edge.attr$Source)&grepl('FUN_ASV',edge.attr$Target)]='Intra_Eukaryotes'
edge.attr$intertype[(grepl('BAC_ASV',edge.attr$Source)&grepl('FUN_ASV',edge.attr$Target))|(grepl('FUN_ASV',edge.attr$Source)&grepl('BAC_ASV',edge.attr$Target))]='Inter_Kingdoms'
edge.attr$intertype -> edge.attributes(g)$intertype

node_attr <- as.data.frame(vertex.attributes(g)) %>% dplyr::rename(Id=name, Ori_Size=Size)
write.csv(edge.attr, paste0(i,'.net.edges.csv'),row.names=F)
write.table(node_attr , paste0(i,'.net.nodes.tsv'), row.names=FALSE, sep='\t')
write.graph(g, paste0(i,'.net.edgelist'),"edgelist")
pdf(file=paste0(i,'.net.pdf'),width=4,height=4)
plot.igraph(g, main=paste0(i,'.net'))
dev.off()

library(qgraph)
pdf(file=paste0(i,'.qgraph.net.pdf'),width=4,height=4)
e <- read.table(file=paste0(i,'.net.edgelist'))
#l <- qgraph.layout.fruchtermanreingold(e,vcount=vcount(g))
l <- qgraph.layout.fruchtermanreingold(e,vcount=vcount(g),
      area=12*(vcount(g)^2),repulse.rad=(vcount(g)^4)) ## repulse.rad can change
plot.igraph(g,layout=l,main=paste0(i,'.net'))
dev.off()
}

## plot network legend
pdf(file='meta_net.legend.pdf')
plot.new()
legend("topright",legend=c("Negative","Positive"), col= c("tomato1","grey10"), lty=1, lwd=3, seg.len=1, bty="n",
title="Relationship")
legend("bottomleft",legend=levels(taxon32), col=color_32, pch=19, pt.cex=2, bty="n",title="Taxonomic Groups")
dev.off()


## (5)calculate region-net basic traits
region.trt=paste_lp(subs.trt,region)

region.trt.net.basic_traits=NULL;
for (i in region.trt)
{
g=get(paste0(i,'_net'))
mean_distance(g)->g1
mean(degree(g))->g2
centr_betw(g)$centralization->g3
centr_clo(g)$centralization->g4
edge_density(g)->g5
assortativity_degree(g)->g6
transitivity(g, type="global")->g7
max(membership(cluster_fast_greedy(g)))->g8
modularity(cluster_fast_greedy(g))->g9
vcount(g)->g10
ecount(g)->g11
diameter(g)->g12
g13<-round(sum(edge.attributes(g)$Relationship=='Positive')/
length(edge.attributes(g)$Relationship)*100,2)
g14<-info.centrality.network(g)
g15<-network.efficiency(g)

mid=vector(length=15) ## mid should be classifed before the loop
for (j in 1:15) {
mid[j]=get(paste0('g',j))
}

names(mid)=c('mean_distance',
'mean_degree',
'centr_betw',
'centr_clo',
'edge_density',
'assortativity_degree',
'transitivity',
'cluster_number',
'modularity',
'vcount',
'ecount',
'diameter',
'PPE',
'infocentrality', 
'efficiency')

rm(list=paste0('g',1:15))
assign(paste0(i,'net.basic_traits'),mid)

if(i==1) region.trt.net.basic_traits=get(paste0(i,'net.basic_traits'));
region.trt.net.basic_traits=rbind(region.trt.net.basic_traits,get(paste0(i,'net.basic_traits')))
}
rownames(region.trt.net.basic_traits)=region.trt

