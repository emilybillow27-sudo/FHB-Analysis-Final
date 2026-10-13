# Historical training versus 2025 UIUC wheat lines. Base R only.
# Rscript scripts/population_comparison.R '/path/to/FHB Analysis Final'
# Optional second argument overrides the results root (not the tables folder).
# VCF diploid GT calls on wheat subgenome chromosomes are treated as disomic.
# Do not use this script for six-copy/polyploid GT calls.
args <- commandArgs(trailingOnly=TRUE)
project <- if (length(args)) args[1] else '.'
results_root <- if (length(args)>1) args[2] else file.path(project,'results')
table_dir <- file.path(results_root,'tables','population_comparison')
figure_dir <- file.path(results_root,'figures')
for (directory in c(table_dir,figure_dir)) dir.create(directory,recursive=TRUE,showWarnings=FALSE)
# Save only the four interpretation tables and the similarity figure.
# Historical is training-set membership, not necessarily an earlier year.
# SNP diversity is conditional on filtering; block CI is conditional on sampled lines.
# Physical blocks may not account fully for long-range LD or family dependence.
# Supplied VCF calls may have been imputed upstream; provenance must be checked.
# Method: https://doi.org/10.1111/j.1558-5646.1984.tb05657.x
max_sample_missing <- 0.20
max_marker_missing <- 0.20
min_maf <- 0.05
max_marker_heterozygosity <- 0.10 # Matches existing PCA marker preparation.
block_bp <- 5000000
bootstrap_replicates <- 1000
set.seed(20251009)
clean_id <- function(x) gsub('[^A-Z0-9]','',toupper(trimws(x)))
write_table <- function(x,name) write.csv(x,file.path(table_dir,paste0(name,'.csv')),row.names=FALSE)
read_data <- function(name) read.csv(file.path(project,'data',name),check.names=FALSE,stringsAsFactors=FALSE)
train <- read_data('FHB_Project_Training_Data.csv')
test <- read_data('FHB_Project_Testing_Data_2025.csv')
stopifnot('FullSampleName' %in% names(train),all(c('FullSampleName','EXPT','YR') %in% names(test)))
test <- test[!is.na(test$EXPT) & !is.na(test$YR) & toupper(trimws(test$EXPT))=='UIUC' & test$YR==2025,]
ids <- function(d) unique(clean_id(d$FullSampleName[!is.na(d$FullSampleName) & d$FullSampleName!='']))
training_ids <- ids(train); testing_ids <- ids(test)
vcf_file <- file.path(project,'data','fhb_analysis_2026_production_final.vcf.gz')
message('Reading original VCF calls...')
con <- gzfile(vcf_file,'rt'); lines <- readLines(con); close(con)
header <- lines[startsWith(lines,'#CHROM')]
if(length(header)!=1) stop('Expected one VCF header.')
samples <- strsplit(header,'\t',fixed=TRUE)[[1]][-(1:9)]
keys <- clean_id(samples)
if(anyDuplicated(keys)) stop('VCF IDs collide after normalization; resolve aliases manually.')
rows <- strsplit(lines[!startsWith(lines,'#')],'\t',fixed=TRUE)
if(any(lengths(rows)!=9+length(samples))) stop('Malformed VCF rows.')
# Restrict to biallelic single nucleotide variants on named wheat chromosomes.
valid <- vapply(rows,function(r) grepl('^[0-9]+[ABD]$',r[1]) && nchar(r[4])==1 && nchar(r[5])==1 &&
  r[4] %in% c('A','C','G','T') && r[5] %in% c('A','C','G','T') && r[7] %in% c('PASS','.'),logical(1))
rows <- rows[valid]
if(!length(rows)) stop('No supported wheat SNPs.')
selected <- which(keys %in% union(training_ids,testing_ids))
if(!length(selected)) stop('No VCF IDs match phenotype identifiers.')
meta <- data.frame(sample=samples[selected],normalized_id=keys[selected],
 historical=keys[selected] %in% training_ids,UIUC_2025=keys[selected] %in% testing_ids)
meta$membership <- ifelse(meta$historical & meta$UIUC_2025,'Shared',ifelse(meta$historical,'Historical_only','UIUC_2025_only'))
chr <- vapply(rows,`[`,character(1),1); pos <- as.numeric(vapply(rows,`[`,character(1),2))
marker <- vapply(rows,`[`,character(1),3)
G <- matrix(NA_real_,nrow=length(selected),ncol=length(rows))
for(j in seq_along(rows)) {
 r <- rows[[j]]; gt_index <- match('GT',strsplit(r[9],':',fixed=TRUE)[[1]])
 if(is.na(gt_index)) stop('VCF row has no GT field.')
 fields <- strsplit(r[9+selected],':',fixed=TRUE)
 gt <- vapply(fields,function(x) if(length(x)>=gt_index) x[gt_index] else '.',character(1))
 if(any(!grepl('^[01.][/|][01.]$|^\\.$',gt))) stop('Unsupported haploid, polyploid or multiallelic GT.')
 G[,j] <- ifelse(gt %in% c('0/0','0|0'),0,ifelse(gt %in% c('1/1','1|1'),2,
   ifelse(gt %in% c('0/1','1/0','0|1','1|0'),1,NA_real_)))
}
rm(lines,rows); invisible(gc(FALSE))
unmatched <- rbind(data.frame(population=rep('Historical',length(setdiff(training_ids,keys))),id=setdiff(training_ids,keys)),
 data.frame(population=rep('UIUC_2025',length(setdiff(testing_ids,keys))),id=setdiff(testing_ids,keys)))
if(nrow(unmatched)) warning('Unmatched phenotype IDs: ',paste(paste(unmatched$population,unmatched$id,sep=': '),collapse=', '))
meta$missing_fraction <- rowMeans(is.na(G))
meta$observed_heterozygosity_before_marker_QC <- rowMeans(G==1,na.rm=TRUE)
meta$pass_sample_QC <- meta$missing_fraction<=max_sample_missing
message('Shared lines before sample QC: ',sum(meta$membership=='Shared'),
 '; lines failing sample QC: ',sum(!meta$pass_sample_QC))
G <- G[meta$pass_sample_QC,,drop=FALSE]; meta <- meta[meta$pass_sample_QC,,drop=FALSE]
if(any(table(factor(meta$membership,levels=c('Historical_only','UIUC_2025_only')))<3)) stop('Too few unique lines after sample QC.')
p <- colMeans(G,na.rm=TRUE)/2
marker_qc <- data.frame(marker,chromosome=chr,position=pos,missing_fraction=colMeans(is.na(G)),
 MAF=pmin(p,1-p),heterozygosity=colMeans(G==1,na.rm=TRUE))
keep <- with(marker_qc,is.finite(MAF) & missing_fraction<=max_marker_missing & MAF>=min_maf & heterozygosity<=max_marker_heterozygosity)
message('SNPs failing marker QC: ',sum(!keep))
G <- G[,keep,drop=FALSE]; chr <- chr[keep]; pos <- pos[keep]
if(ncol(G)<10) stop('Too few markers after filtering.')
message('Retained ',nrow(G),' lines and ',ncol(G),' SNPs.')
# All diversity statistics refer to the SAME prefiltered SNP panel, not all bases.
groups <- list(Historical_all=which(meta$historical),UIUC_2025_all=which(meta$UIUC_2025),
 Historical_only=which(meta$membership=='Historical_only'),UIUC_2025_only=which(meta$membership=='UIUC_2025_only'))
diversity <- do.call(rbind,lapply(names(groups),function(name) {
 x <- G[groups[[name]],,drop=FALSE]; p <- colMeans(x,na.rm=TRUE)/2
 data.frame(population=name,n_lines=nrow(x),n_SNPs=ncol(x),
 mean_expected_heterozygosity=mean(2*p*(1-p),na.rm=TRUE),
 mean_observed_heterozygosity=mean(colMeans(x==1,na.rm=TRUE),na.rm=TRUE),
 fraction_polymorphic_among_callable_SNPs=mean(p>0 & p<1,na.rm=TRUE),missing_fraction=mean(is.na(x)))
}))
write_table(diversity,'SNP_panel_diversity')
# Weir & Cockerham (1984) variance components for two diploid populations.
# Aggregate as sum(a)/sum(a+b+c), never the mean of per-SNP ratios.
wc_components <- function(x,y) {
 n1 <- colSums(!is.na(x)); n2 <- colSums(!is.na(y)); nt <- n1+n2; nb <- nt/2
 p1 <- colMeans(x,na.rm=TRUE)/2; p2 <- colMeans(y,na.rm=TRUE)/2
 pb <- (n1*p1+n2*p2)/nt
 s2 <- (n1*(p1-pb)^2+n2*(p2-pb)^2)/nb
 hb <- (colSums(x==1,na.rm=TRUE)+colSums(y==1,na.rm=TRUE))/nt
 nc <- nt-(n1^2+n2^2)/nt
 a <- nb/nc*(s2-(pb*(1-pb)-s2/2-hb/4)/(nb-1))
 b <- nb/(nb-1)*(pb*(1-pb)-s2/2-(2*nb-1)*hb/(4*nb))
 c <- hb/2
 ok <- n1>=2 & n2>=2 & is.finite(a+b+c) & a+b+c>0
 data.frame(a=a[ok],denominator=(a+b+c)[ok],chr=chr[ok],pos=pos[ok])
}
z <- wc_components(G[groups$Historical_only,,drop=FALSE],G[groups$UIUC_2025_only,,drop=FALSE])
if(!nrow(z)) stop('No valid FST loci.')
fst <- sum(z$a)/sum(z$denominator)
# Resample physical blocks within chromosome. CI reflects marker/block uncertainty
# conditional on the sampled lines. It does not resample breeding families.
block <- paste(z$chr,floor((z$pos-1)/block_bp),sep=':')
blocks <- aggregate(z[c('a','denominator')],list(block=block,chromosome=z$chr),sum)
by_chr <- split(seq_len(nrow(blocks)),blocks$chromosome)
boot <- replicate(bootstrap_replicates,{
 k <- unlist(lapply(by_chr,function(ix) ix[sample.int(length(ix),length(ix),replace=TRUE)]),use.names=FALSE)
 sum(blocks$a[k])/sum(blocks$denominator[k])
})
ci <- quantile(boot,c(.025,.975),na.rm=TRUE)
write_table(data.frame(comparison='Historical_only_vs_UIUC_2025_only',FST=fst,CI_lower=ci[1],CI_upper=ci[2],
 n_historical=length(groups$Historical_only),n_UIUC_2025=length(groups$UIUC_2025_only),
 n_SNPs=nrow(z),n_blocks=nrow(blocks),block_bp,bootstrap_replicates),'FST_summary')
# Exploratory VanRaden-style genomic relationship matrix; mean impute missing
# dosages ONLY here. This is marker similarity, not pedigree kinship or IBD.
p <- colMeans(G,na.rm=TRUE)/2
X <- sweep(G,2,2*p,'-'); X[is.na(X)] <- 0
A <- tcrossprod(X)/sum(2*p*(1-p))
dimnames(A) <- list(meta$sample,meta$sample)
pair_rows <- function(i,j,label,same=FALSE) {
 M <- A[i,j,drop=FALSE]; vals <- if(same) M[upper.tri(M)] else as.vector(M)
 data.frame(pair_type=label,n_pairs=length(vals),mean=mean(vals),median=median(vals),
 q05=unname(quantile(vals,.05)),q95=unname(quantile(vals,.95)))
}
i <- groups$Historical_only; j <- groups$UIUC_2025_only
write_table(rbind(pair_rows(i,i,'Within historical-only',TRUE),pair_rows(j,j,'Within UIUC-only',TRUE),
 pair_rows(i,j,'Between unique sets')),'genomic_similarity_summary')
# Nearest historical line for every testing line; shared lines flagged explicitly.
h <- groups$Historical_all; t <- groups$UIUC_2025_all
nearest <- do.call(rbind,lapply(t,function(k) {
 candidates <- setdiff(h,k) # Exclude the same line from its own match.
 idx <- candidates[which.max(A[k,candidates])]
 data.frame(UIUC_line=meta$sample[k],also_in_historical=meta$historical[k],
 closest_other_historical_line=meta$sample[idx],genomic_relationship=A[k,idx])
}))
write_table(nearest,'UIUC_nearest_historical_lines')
pdf(file.path(figure_dir,'genomic_similarity.pdf'),width=8,height=5)
boxplot(list(Historical=A[i,i][upper.tri(A[i,i])],UIUC=A[j,j][upper.tri(A[j,j])],Between=as.vector(A[i,j])),
 ylab='Genomic relationship (relative to pooled allele frequencies)',main='Similarity among unique breeding lines',outline=FALSE,
 col=c('#214A68','#7699B2','#D4AD6A'))
dev.off()
message('Finished. FST = ',signif(fst,5),'; tables: ',normalizePath(table_dir),
 '; figure: ',file.path(normalizePath(figure_dir),'genomic_similarity.pdf'))
