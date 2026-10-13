# University-city markers are approximate, not exact field-plot coordinates.
library(ggplot2)
library(maps)
project_dir <- getOption('fhb.project_dir','/Users/emilybillow/Desktop/FHB Analysis Final')
figure_dir <- getOption('fhb.figure_dir',file.path(project_dir,'results','figures'))
dir.create(figure_dir,recursive=TRUE,showWarnings=FALSE)
states <- map_data('state')
focus <- c('south dakota','nebraska','illinois')
sites <- data.frame(site=c('SDSU','UNL','UIUC'),
 city=c('Brookings, SD','Lincoln, NE','Urbana-Champaign, IL'),
 long=c(-96.80,-96.68,-88.23),lat=c(44.31,40.81,40.11),
 label_long=c(-99.0,-99.0,-86.7),label_lat=c(45.0,39.6,41.4))
state_labels <- data.frame(long=c(-101.3,-101.0,-89.5),
 lat=c(44.3,41.7,38.5),
 name=c('SOUTH DAKOTA','NEBRASKA','ILLINOIS'))
p <- ggplot() +
 geom_polygon(data=states,aes(long,lat,group=group),fill='#F1F3F5',color='white',linewidth=.6) +
 geom_polygon(data=subset(states,region %in% focus),aes(long,lat,group=group),fill='#D8E5EE',color='white',linewidth=.8) +
 geom_text(data=state_labels,aes(long,lat,label=name),color='#75828B',size=3.6) +
 geom_segment(data=sites,aes(x=long,y=lat,xend=label_long,yend=label_lat),color='#214A68',linewidth=.6) +
 geom_point(data=sites,aes(long,lat),shape=21,fill='#214A68',color='white',size=4.5,stroke=1) +
 geom_label(data=sites,aes(label_long,label_lat,label=paste(site,city,sep='\n')),
 color='#214A68',fill='white',fontface='bold',size=4.5,label.padding=grid::unit(.3,'lines'),linewidth=.3) +
 coord_quickmap(xlim=c(-104.5,-84),ylim=c(36.5,46.8),expand=FALSE) +
 labs(title='2026 testing locations',
 caption='Markers indicate university cities, not exact nursery coordinates. State boundaries: R maps package.') +
 theme_void(base_size=17) +
 theme(plot.title=element_text(face='bold',color='#214A68',size=24),
 plot.subtitle=element_text(color='#536773',size=16,margin=margin(b=12)),
 plot.caption=element_text(color='#687985',size=10,hjust=0,margin=margin(t=10)),
 plot.margin=margin(16,16,16,16),plot.background=element_rect(fill='white',color=NA))
ggsave(file.path(figure_dir,'2026_testing_locations_map.png'),p,width=10,height=8,dpi=300,bg='white')
