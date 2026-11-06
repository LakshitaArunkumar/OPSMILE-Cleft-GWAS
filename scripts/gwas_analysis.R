#!/usr/bin/env Rscript
# Public reconstruction: prepare approved, de-identified phenotype tables and
# summarize PLINK 1.9 additive logistic association output on secure storage.
# Not a release of original participant-level data or clinical adjudication.
suppressPackageStartupMessages({library(data.table);library(ggplot2)})
args <- commandArgs(trailingOnly=TRUE)
if (!length(args) || args[1] %in% c('-h','--help')) {
 cat('Rscript scripts/gwas_analysis.R prepare --clinical FILE --fam FILE --pca FILE --outdir DIR\n',
     'Rscript scripts/gwas_analysis.R report --gwas-dir DIR --outdir DIR\n')
 quit(status=if(length(args))0 else 2)
}
mode <- args[1]; opts <- args[-1]
if (length(opts)%%2) stop('Use --key value pairs')
opts <- setNames(as.list(opts[seq(2,length(opts),2)]),sub('^--','',opts[seq(1,length(opts),2)]))
opt <- function(x) {if (is.null(opts[[x]])) stop('Missing --',x); opts[[x]]}
out <- opt('outdir');dir.create(out,recursive=TRUE,showWarnings=FALSE)
fwrite_tab <- function(x,p) fwrite(x,p,sep='\t',quote=FALSE,na='NA')
if (mode=='prepare') {
 clinical <- fread(opt('clinical'),colClasses='character')
 req <- c('FID','IID','case','subtype')
 if(!all(req %in% names(clinical)))stop('Clinical manifest requires FID,IID,case,subtype')
 if(anyDuplicated(clinical[,.(FID,IID)])) stop('Duplicate clinical sample keys')
 if(!all(na.omit(unique(clinical$case)) %in% c('0','1'))) stop('case must be 0 or 1')
 fam <- fread(opt('fam'),header=FALSE,colClasses='character')
 if(ncol(fam)<6)stop('Invalid FAM')
 setnames(fam,names(fam)[1:6],c('FID','IID','PID','MID','SEX','FAM_PHENO'))
 if(anyDuplicated(fam[,.(FID,IID)]))stop('Duplicate FAM sample keys')
 pca <- fread(opt('pca'),header=FALSE,colClasses='character')
 if(ncol(pca)<7)stop('Expected five PCs')
 setnames(pca,names(pca)[1:7],c('FID','IID',paste0('PC',1:5)))
 if(anyDuplicated(pca[,.(FID,IID)]))stop('Duplicate PCA sample keys')
 d <- merge(fam[,.(FID,IID,SEX)],clinical,by=c('FID','IID'),all.x=TRUE)
 d <- merge(d,pca[,c('FID','IID',paste0('PC',1:5)),with=FALSE],by=c('FID','IID'),all.x=TRUE)
 if(anyNA(d[,c('SEX',paste0('PC',1:5)),with=FALSE]))stop('Missing covariates: review before analysis')
 if(!all(d$SEX %in% c('1','2')))stop('Unknown/invalid sex needs adjudication')
 d[,subtype := trimws(subtype)]
 valid_case <- d$case=='1'
 valid_control <- d$case=='0' & (is.na(d$subtype)|d$subtype %in% c('','99'))
 valid_case[is.na(valid_case)]<-FALSE
 valid_control[is.na(valid_control)]<-FALSE
 # Documented phenotype definitions: unknown subtype (99/NA) in ALL only.
 groups <- list(ALL=as.character(c(1:7,99)),CL=c('1','2','6'),
                CLP=c('3','4','7'),CLPM=c('1','2','3','4','6','7'),CP='5')
 for (nm in names(groups)) {
  is_case <- valid_case & (d$subtype %in% groups[[nm]] |
                              (nm=='ALL' & (is.na(d$subtype)|d$subtype=='')))
  is_case[is.na(is_case)]<-FALSE
  keep <- is_case|valid_control
  sub <- d[keep,.(FID,IID)]
  sub[,PHENO:=ifelse(is_case[keep],2L,1L)]
  if(!nrow(sub)||!any(sub$PHENO==2)||!any(sub$PHENO==1)) stop('No usable cases/controls for ',nm)
  fwrite_tab(sub,file.path(out,paste0('pheno_',nm,'.tsv')))
 }
 for (n in c(3,5)) {
  names_needed <- c('FID','IID','SEX',paste0('PC',seq_len(n)))
  fwrite_tab(d[,..names_needed],file.path(out,paste0('covar_PC',n,'.tsv')))
 }
 message('Prepared PRIVATE per-subtype phenotypes and covariates')
} else if (mode=='report') {
 files <- list.files(opt('gwas-dir'),pattern='\\.assoc\\.logistic$',full.names=TRUE)
 if(!length(files))stop('No logistic results')
 summary <- list()
 for(f in files) {
  label <- sub('\\.assoc\\.logistic$','',basename(f));x <- fread(f)
  if(!all(c('CHR','BP','SNP','TEST','P') %in% names(x)))stop('Invalid PLINK output: ',f)
  x <- x[TEST=='ADD'];x[,P:=as.numeric(P)];x[,CHR:=as.integer(CHR)];x[,BP:=as.numeric(BP)]
  x <- x[is.finite(P)&P>0&P<=1&CHR %in% 1:22&is.finite(BP)]
  if(!nrow(x))next
  setorder(x,CHR,BP)
  lam <- median(qchisq(x$P,1,lower.tail=FALSE))/qchisq(.5,1,lower.tail=FALSE)
  summary[[label]] <- data.table(Model=label,Tested=nrow(x),Lambda_GC=lam,
    Genomewide=sum(x$P<=5e-8),Suggestive=sum(x$P>5e-8 & x$P<=1e-5),Minimum_P=min(x$P))
  chr <- x[,.(length=max(BP)),by=CHR][order(CHR)]
  chr[,offset:=shift(cumsum(length),fill=0)];chr[,center:=offset+length/2]
  plotdat <- merge(x,chr[,.(CHR,offset)],by='CHR')
  plotdat[,`:=`(position=BP+offset,neglogp=-log10(P),stripe=factor(CHR%%2))]
  p <- ggplot(plotdat,aes(position,neglogp,color=stripe))+geom_point(size=.35,alpha=.65)+
    scale_color_manual(values=c('grey40','grey70'),guide='none')+
    geom_hline(yintercept=-log10(5e-8),linetype='dashed')+
    geom_hline(yintercept=-log10(1e-5),linetype='dotted')+
    scale_x_continuous(breaks=chr$center,labels=chr$CHR)+
    labs(title=label,x='Chromosome',y=expression(-log[10](P)))+theme_bw()
  ggsave(file.path(out,paste0(label,'_manhattan.png')),p,width=12,height=5,dpi=180)
  n <- nrow(x);q <- data.table(expected=-log10(ppoints(n)),observed=-log10(sort(x$P)))
  z <- ggplot(q,aes(expected,observed))+geom_point(size=.45,alpha=.7)+
    geom_abline(slope=1,intercept=0,linetype='dashed')+theme_bw()+
    labs(title=paste(label,'QQ'),subtitle=sprintf('Lambda GC: %.3f',lam),
       x='Expected -log10(P)',y='Observed -log10(P)')
  ggsave(file.path(out,paste0(label,'_qq.png')),z,width=6,height=5,dpi=180)
  # PRIVATE outputs only: per-variant summary is NOT public-release approved.
  fwrite(x,file.path(out,paste0(label,'_association_results_PRIVATE.tsv')),sep='\t')
 }
 if(!length(summary))stop('No valid additive association results')
 fwrite(rbindlist(summary),file.path(out,'summary_PRIVATE.csv'))
 message('Saved outputs to PRIVATE directory; do not commit them')
} else stop('Unknown mode: ',mode)
