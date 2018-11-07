#run this code if processing from raw copy-number data
#download copy-number segment table from https://www.synapse.org/#!Synapse:syn1710464 into data/ directory
tcga<-read.table("data/pancan12_absolute.segtab.txt",sep="\t",header=T,stringsAsFactors = F)

tcga_clin<-read.table("data/tcga_sample_info.tsv",stringsAsFactors = F,header=T,sep="\t")
tcga<-tcga[tcga$Sample%in%tcga_clin$bcr_aliquot_barcode,]

pcawg_clin<-read.table("data/pcawg_sample_info.tsv",sep="\t",header=T,stringsAsFactors = F)

#filter cases that appear in pcawg
tcga<-tcga[!substr(tcga$Sample,1,12)%in%pcawg_clin$submitted_donor_id,]

#download sample info file from https://www.synapse.org/#!Synapse:syn1710466 in data/ directory
tcga_info<-read.table("data/pancan12.sample_info.txt",sep="\t",header=T,stringsAsFactors = F)
rownames(tcga_info)<-tcga_info$tcga_id
tcga_info<-tcga_info[tcga_info$abs_call=="called",]
tcga_info<-tcga_info[tcga_info$tcga_id%in%tcga$Sample,]

tcga_segTabs<-list()
for(i in unique(tcga$Sample))
{
    tab<-tcga[tcga$Sample==i,c("Chromosome","Start","End","Expected_HSCN_a1","Expected_HSCN_a2")]
    tab$segVal<-tab$Expected_HSCN_a1+tab$Expected_HSCN_a2
    tab<-tab[,c(-4,-5)]
    colnames(tab)<-c("chromosome","start","end","segVal")
    tcga_segTabs[[i]]<-tab
}

tcga_CN_features<-extractCopynumberFeatures(tcga_segTabs)


tcga_CN_features<-readRDS("data/tcga_CN_features.rds")
tcga_sample_component_matrix<-generateSampleByComponentMatrix(tcga_CN_features,CN_components)
NMF::aheatmap(tcga_sample_component_matrix,Rowv=NULL, main="Sample x Component matrix")

tcga_ids<-rownames(tcga_sample_component_matrix)

tcga_sigs<-NMF::nmf(t(tcga_sample_component_matrix),nsig,seed=seed,nrun=1000,method=nmfalg,.opt = "p16")
coefmap(tcga_sigs,Colv="consensus",tracks=c("consensus:"), main="Sample x Signature matrix")#annCol=annCol,
basismap(tcga_sigs,Rowv=NA,main="Signature x Component matrix")


#--------- 从CN_feature data开始建模以及寻找最优signature数目
seed=77777
min_prior=0.001
model_selection="BIC"
nrep=1
niter=1000

dat<-as.numeric(CN_features[["segsize"]][,2])
segsize_mm<-fitComponent(dat,seed=seed,model_selection=model_selection,
                         min_prior=min_prior,niter=niter,nrep=nrep,min_comp=10,max_comp=10)

dat<-as.numeric(CN_features[["bp10MB"]][,2])
bp10MB_mm<-fitComponent(dat,dist="pois",seed=seed,model_selection=model_selection,
                        min_prior=min_prior,niter=niter,nrep=nrep,min_comp=3,max_comp=3)

dat<-as.numeric(CN_features[["osCN"]][,2])
osCN_mm<-fitComponent(dat,dist="pois",seed=seed,model_selection=model_selection,
                      min_prior=min_prior,niter=niter,nrep=nrep,min_comp=3,max_comp=3)

dat<-as.numeric(CN_features[["bpchrarm"]][,2])
bpchrarm_mm<-fitComponent(dat,dist="pois",seed=seed,model_selection=model_selection,
                          min_prior=min_prior,niter=niter,nrep=3,min_comp=2,max_comp=5)

dat<-as.numeric(CN_features[["changepoint"]][,2])
changepoint_mm<-fitComponent(dat,seed=seed,model_selection=model_selection,
                             min_prior=min_prior,niter=niter,nrep=nrep,min_comp=7,max_comp=7)

dat<-as.numeric(CN_features[["copynumber"]][,2])
copynumber_mm<-fitComponent(dat,seed=seed,model_selection=model_selection,
                            nrep=nrep,min_comp=2,max_comp=10,min_prior=0.005,niter=2000)

CN_components<-list(segsize=segsize_mm,bp10MB=bp10MB_mm,osCN=osCN_mm,changepoint=changepoint_mm,copynumber=copynumber_mm,bpchrarm=bpchrarm_mm)

britroc_sample_component_matrix<-generateSampleByComponentMatrix(CN_features,CN_components,cores=1,subcores=num_cores)
NMF::aheatmap(britroc_sample_component_matrix,fontsize = 7,Rowv=FALSE,Colv=FALSE,legend = T,breaks=c(seq(0,199,2),500),main="Component x Sample matrix")


nmfalg<-"brunet"
seed<-77777
estim.r <- NMF::nmfEstimateRank(t(britroc_sample_component_matrix), 3:12,seed = seed,nrun=1000,
                                verbose=F,method=nmfalg,.opt = paste0("p",num_cores))
V.random <- randomize(t(britroc_sample_component_matrix))
estim.r.random <- NMF::nmfEstimateRank(V.random, 3:12, seed =seed,nrun=1000,
                                       verbose=F,method=nmfalg,.opt = paste0("p",num_cores))
p<-plot(estim.r,estim.r.random,
        what = c("cophenetic", "dispersion","sparseness", "silhouette"),xname="Observed",yname="Randomised",main="")+
    theme(axis.text=element_text(size=5),axis.title=element_text(size=5),
          strip.text.x = element_text(size = 5),
          strip.text.y = element_text(size = 5),
          legend.text = element_text(size = 5),
          legend.title = element_text(size = 7))
g<-ggplotGrob(p)
g[["grobs"]][[2]]$children[[4]]$size[[1]]<-0.5
g[["grobs"]][[3]]$children[[4]]$size[[1]]<-0.5
g[["grobs"]][[4]]$children[[4]]$size[[1]]<-0.5
g[["grobs"]][[5]]$children[[4]]$size[[1]]<-0.5
grid::grid.newpage()
grid::grid.draw(g)
