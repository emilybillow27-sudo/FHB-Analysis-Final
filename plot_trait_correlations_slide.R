# Standalone slide figures. Does not refit or modify your analysis.
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
project_dir <- getOption('fhb.project_dir', '/Users/emilybillow/Desktop/FHB Analysis Final')
out_dir <- getOption('fhb.figure_dir', file.path(project_dir, 'results', 'figures'))
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
traits <- c('INC','SEV','DON')
labels <- c(INC='Incidence', SEV='Severity', DON='DON')
bg <- '#F2F4F6'
prepare <- function(file, id) {
  d <- readRDS(file.path(project_dir,'results',file))
  stopifnot(all(c(id,'TRAIT','predicted.value') %in% names(d)))
  d <- d |> filter(TRAIT %in% traits)
  if(any(is.na(d[[id]]) | trimws(d[[id]])=='')) stop('Missing entry identifiers')
  # Native identifiers preserve entries that lack marker-matched FullSampleName.
  w <- d |> group_by(across(all_of(id)),TRAIT) |>
    summarise(value=if(all(!is.finite(predicted.value))) NA_real_ else
      mean(predicted.value[is.finite(predicted.value)]), .groups='drop') |>
    pivot_wider(names_from=TRAIT, values_from=value)
  m <- as.matrix(w[,traits])
  r <- cor(m, use='pairwise.complete.obs', method='pearson')
  n <- crossprod(1L*is.finite(m))
  list(wide=w, r=r, n=n)
}
train <- prepare('blues_me_train.rds','germplasmName')
test <- prepare('blues_se_test.rds','ID')
make_plot <- function(x,title,subtitle) {
  d <- as.data.frame(as.table(x$r)); names(d)<-c('x','y','r')
  d$x <- factor(d$x,levels=traits)
  d$y <- factor(d$y,levels=rev(traits))
  ggplot(d,aes(x,y,fill=r)) +
    geom_tile(color=bg,linewidth=2) +
    geom_text(aes(label=sprintf('%.2f',r),color=abs(r)>.65),size=7,fontface='bold',show.legend=FALSE) +
    scale_color_manual(values=c('FALSE'='#203542','TRUE'='white')) +
    scale_fill_gradient2(low='#326BA0',mid='white',high='#B43E42',midpoint=0,limits=c(-1,1),name='Pearson correlation',breaks=c(-1,0,1)) +
    scale_x_discrete(labels=c(INC='INC',SEV='SEV',DON='DON'),position='bottom',expand=c(0,0)) +
    scale_y_discrete(labels=labels,expand=c(0,0)) +
    coord_fixed() + labs(title=title,subtitle=subtitle,x=NULL,y=NULL) +
    theme_minimal(base_size=18) +
    theme(panel.grid=element_blank(),plot.background=element_rect(fill=bg,color=NA),
          panel.background=element_rect(fill=bg,color=NA),
          plot.title=element_text(face='bold',size=23,color='#214A68'),
          plot.subtitle=element_text(size=15,color='#536773',margin=margin(b=14)),
          axis.text=element_text(size=18,color='#203542'),
          axis.ticks=element_blank(),legend.position='bottom',
          legend.title=element_text(size=14),legend.text=element_text(size=13),
          plot.margin=margin(12,18,12,12))
}
p_train <- make_plot(train,'Historical training','Across-environment estimates')
p_test <- make_plot(test,'2025 UIUC','Adjusted estimates from\nreplicated nurseries')
combined <- (p_train | p_test) + plot_layout(guides='collect') + plot_annotation(theme=theme(plot.background=element_rect(fill=bg,color=NA))) & theme(legend.position='bottom')
ggsave(file.path(out_dir,'FHB_trait_correlations_side_by_side.png'),combined,width=13,height=5.5,dpi=300,bg=bg)
ggsave(file.path(out_dir,'FHB_trait_correlations_side_by_side.pdf'),combined,width=13,height=5.5,bg=bg)
for(pop in c('train','test')) {
  x <- get(pop)
  write.csv(x$wide,file.path(out_dir,paste0('FHB_correlations_',pop,'_input.csv')),row.names=FALSE)
  d <- as.data.frame(as.table(x$r));names(d)<-c('Trait_1','Trait_2','Pearson_r')
  d$Paired_entries <- as.vector(x$n)
  write.csv(d,file.path(out_dir,paste0('FHB_correlations_',pop,'_values.csv')),row.names=FALSE)
  print(d)
}
ggsave(file.path(out_dir,'FHB_correlations_historical.png'),p_train,width=6.5,height=5.5,dpi=300,bg=bg)
ggsave(file.path(out_dir,'FHB_correlations_UIUC_2025.png'),p_test,width=6.5,height=5.5,dpi=300,bg=bg)
