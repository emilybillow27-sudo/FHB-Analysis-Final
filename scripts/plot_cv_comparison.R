# Plot existing CV results only; does not rerun genomic prediction.
library(dplyr)
library(ggplot2)
project_dir <- getOption('fhb.project_dir','/Users/emilybillow/Desktop/FHB Analysis Final')
figure_dir <- getOption('fhb.figure_dir',file.path(project_dir,'results','figures'))
dir.create(figure_dir,recursive=TRUE,showWarnings=FALSE)
h <- read.csv(file.path(project_dir,'results','tables','cross_validation','historical','historical_cv_results.csv'))
e <- read.csv(file.path(project_dir,'results','tables','cross_validation','expanded','cv_results_expanded.csv'))
d <- bind_rows(mutate(h,Population='Historical only'),mutate(e,Population='Expanded'))
stopifnot(all(c('TRAIT','REPETITION','N','PREDICTION_ABILITY') %in% names(d)))
d$TRAIT <- factor(d$TRAIT,levels=c('INC','SEV','DON','FDK'))
d$Population <- factor(d$Population,levels=c('Historical only','Expanded'))
s <- d |> group_by(TRAIT,Population) |> summarise(mean_r=mean(PREDICTION_ABILITY),n_min=min(N),n_max=max(N),reps=n(),.groups='drop')
stopifnot(all(s$reps==100),all(s$n_min==s$n_max))
s$label <- sprintf('%.2f',s$mean_r)
s$label_y <- vapply(seq_len(nrow(s)), function(i) {
  v <- d$PREDICTION_ABILITY[d$TRAIT == s$TRAIT[i] & d$Population == s$Population[i]]
  unname(boxplot.stats(v)$stats[5]) + 0.04
}, numeric(1))
p <- ggplot(d,aes(Population,PREDICTION_ABILITY,fill=Population)) +
 geom_hline(yintercept=0,color='#9CA7AE',linetype='dashed') +
 geom_boxplot(width=.55,outlier.size=1.3) +
 geom_text(data=s,mapping=aes(x=Population,y=label_y,label=label),
           inherit.aes=FALSE,size=5.5,fontface='bold',color='#214A68') +
 facet_wrap(~TRAIT,ncol=2) +
 scale_fill_manual(values=c('Historical only'='#214A68','Expanded'='#7699B2')) +
 scale_y_continuous(limits=c(-.2,.8),breaks=seq(-.2,.8,.2)) +
 labs(title='Historical and expanded population cross-validation',
 subtitle='100 repetitions of five-fold cross-validation',
 x=NULL,y='Predictive ability (Pearson r)',
 caption='Boxplots show variation across repetitions. Validation lines differ between populations.\nThese results do not evaluate prediction of the 2026 testing population.') +
 theme_minimal(base_size=19) +
 theme(legend.position='none',panel.grid.minor=element_blank(),panel.grid.major.x=element_blank(),
 axis.text.x=element_text(size=17),strip.text=element_text(face='bold',size=21),
 plot.title=element_text(face='bold',color='#214A68',size=21),plot.caption=element_text(hjust=0,size=12),panel.spacing=grid::unit(1.8,"lines"),
 panel.border=element_rect(color="#687985",fill=NA,linewidth=0.9),
 strip.background=element_rect(fill="#EDF1F4",color="#687985",linewidth=0.9))
ggsave(file.path(figure_dir,'CV_historical_expanded_comparison.png'),p,width=12,height=10,dpi=300,bg='white')
print(s)
