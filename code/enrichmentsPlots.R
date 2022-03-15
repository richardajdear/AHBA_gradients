suppressPackageStartupMessages(library(tidyverse))
library(pals)
library(ggridges)
library(ggrepel)
library(ggradar)



plot_enrichments <- function(enrichments, n=10, facet='v') {
    
    if (facet=='v') {
        p <- enrichments %>%
        ggplot(aes(enrichment, neglogFDR)) + 
        facet_wrap(~G, scales='free_x')
    } else {
        p <- enrichments %>%
        ggplot(aes(enrichment, neglogFDR)) + 
        facet_wrap(~G, scales='free')
    }
    
    p +
    # geom_point(aes(size=n_genes), alpha=.8, fill=brewer.blues(5)[2], shape=21) +
    geom_point(aes(size=n_genes, fill=direction), alpha=.8, shape=21) +
    geom_text_repel(aes(label=description), data = enrichments %>% filter(rank < n), force=100, force_pull=.1, nudge_y=-.1, size=9) +
    scale_y_continuous(name='FDR p-value', breaks=c(-log10(0.05),2,3,4,5), labels=function(x) 10^(-x)) +
    scale_size_continuous(name='# genes', range=c(1,20)) +
    scale_fill_manual(values=brewer.rdbu(5)[c(2,4)], guide='none') +
    xlab('Enrichment Ratio') +
    theme_minimal() + 
    theme(panel.grid = element_blank(),
          panel.border = element_rect(fill=NA),
          # legend.position=c(.9,.9),
          legend.position='bottom',
          strip.text.y = element_text(size=36, angle=0)
         )
    
}

plot_cell_dotplot <- function(null_p) {
    null_p %>% 
    ggplot(aes(x=celltype, y=z, group=rev(G), fill=G)) + 
    geom_col(aes(), width=.2, position=position_dodge(.5)) +
    geom_point(aes(alpha=!sig), shape=21, size=5, position=position_dodge(.5)) +
    geom_point(aes(alpha=sig), fill='white', shape=21, size=5, position=position_dodge(.5)) +
    scale_alpha_manual(values=c(1,0)) +
    geom_hline(yintercept=0) +
    scale_x_discrete(limits=rev, name='') +
    scale_y_continuous(limits=c(-20,20), name='Z-score of true mean weight vs permuted mean weights') +
    scale_fill_manual(values=brewer.rdylbu(5)[c(1,2,5)]) +
    coord_flip() +
    theme_minimal() +
    theme(panel.grid.minor=element_blank(),
          text=element_text(size=30)
         )
}

plot_cell_radar <- function(null_p) {
    null_p %>% 
    select(celltype, G, z) %>%
    spread(celltype, z) %>%
    ggradar(
            grid.min=-20, grid.max=20,
            values.radar = c("z = -20", "z = 0", "z = +20"),
            group.point.size = 3,
            base.size = 36,
            axis.label.size = 12,
            grid.label.size = 12,
            fill=T,
            fill.alpha=.05,
            group.colours=brewer.rdylbu(5)[c(1,2,5)],
            background.circle.colour = "white",
            gridline.mid.colour = "grey",
            legend.text.size = 36,
            legend.position = "left"
    )
}








plot_cell_enrichment <- function(true_scores, null_scores, null_p, how='mean', p_sig = .05/7) {

    mycolors = c(brewer.rdylbu(6)[1:3],brewer.rdylbu(5)[4:5])

    null_scores <- null_scores %>% gather(pc, score, -m)
    null_p <- null_p %>% rownames_to_column('m') %>% gather(pc, p, -m) %>% mutate(sig = ifelse(p < p_sig | p > 1-p_sig, T, F)) 
    null_scores <- null_p %>% left_join(null_scores, by=c('pc', 'm')) %>% mutate(m = factor(m, levels=rev(.$m %>% unique)))
    true_scores <- true_scores %>% rownames_to_column('m') %>% gather(pc, score, -m) %>% left_join(null_p, by=c('pc', 'm'))


    p <- null_scores %>% 
    ggplot() + 
    facet_grid(.~pc) +
    stat_density_ridges(geom = "density_ridges_gradient", calc_ecdf = TRUE, scale=.8, size=.2, color='grey',
                       aes(x=score, y=m, fill=stat(ecdf))) +
    scale_fill_gradientn(colors=rev(brewer.rdbu(100)[15:85]), name = 'Quantile') +
    geom_point(data=true_scores, aes(x=score, y=m, fill=p), shape=21, size=15) +
    geom_text(data=true_scores, aes(label=paste0('p=',round(p,2)), x=median(score), y=m), 
              size=12, hjust=-.2, nudge_y=.3) +
    guides(fill = guide_colorbar(barwidth=20, title.vjust=1))
    
    if (how == 'mean') {
        p <- p + xlab('Mean weight of genes associated with cell-type')
    } else if (how == 'median') {
        p <- p + xlab('Median rank of genes associated with cell-type')
    }
        
    p +
    theme_void() + 
    theme(
        axis.title.x = element_text(size=36, vjust=5), 
        axis.text.y = element_text(size=36), 
        strip.text.y.left = element_text(size=36, angle=0), 
        # strip.text.x = element_text(size=30, angle=0), 
        strip.text.x = element_blank(), 
        panel.spacing=unit(20, 'lines'), 
        text = element_text(size=30),
        legend.position='bottom'
    )
}



library(ggrepel)

plot_revigo_data <- function(revigo_R, text_threshold=.15, title='') {
    source(revigo_R)
    
    ex <- one.data [ one.data$dispensability < text_threshold, ];
    
    p1 <- ggplot( data = one.data ) + 
        # geom_vline(xintercept=0, colour = I (alpha ("black", 0.6) )) + 
        # geom_hline(yintercept=0, colour = I (alpha ("black", 0.6) )) +
        geom_point( aes(plot_X, plot_Y, size = log_size, fill=value), shape = 21, colour = I (alpha ("black", 0.6) ), alpha=.6) + 
        geom_label_repel(data = ex, aes(plot_X, plot_Y, label = description), colour = I(alpha("black", 0.85)), 
                         size = 8, box.padding = 1, force=5) +
        scale_size(range=c(5,30), name='GO terms') +
        # scale_fill_gradientn(colours = brewer.ylorrd(100)) +
        scale_fill_gradientn(colours = brewer.reds(100)) +
        guides(fill = guide_colorbar(barwidth=20, title='Log p-value', title.vjust=1), size='none') +
        # coord_cartesian(clip='off') +
        ggtitle(title) +
        theme_bw() +
        theme(panel.grid = element_blank(),
              panel.border = element_rect(colour = "grey", fill=NA, size=1),
              axis.text = element_blank(),
              axis.ticks = element_blank(),
              text = element_text(size=30),
              plot.title = element_text(hjust=.5, size=30),
              axis.title.x=element_text(vjust=5),
              legend.position='bottom'
              # aspect.ratio = 1
             ) + 
        # labs (y = "Semantic Space Y", x = "Semantic Space X") + 
        xlab('') + ylab('') +
        theme(legend.key = element_blank()) ;
    
    one.x_range = max(one.data$plot_X) - min(one.data$plot_X);
    one.y_range = max(one.data$plot_Y) - min(one.data$plot_Y);
    p1 <- p1 + xlim(min(one.data$plot_X)-one.x_range/10,max(one.data$plot_X)+one.x_range/10);
    p1 <- p1 + ylim(min(one.data$plot_Y)-one.y_range/10,max(one.data$plot_Y)+one.y_range/10);
    
    p1
}








########### OLD


plot_enrichment_nulls <- function(true_scores, null_scores, null_p, how='mean') {

    mycolors = c(brewer.rdylbu(6)[1:3],brewer.rdylbu(5)[4:5])

    true_scores <- true_scores %>% rownames_to_column('m') %>% gather(pc, score, -m)
    null_scores <- null_scores %>% gather(pc, score, -m)
    null_p <- null_p %>% rownames_to_column('m') %>% gather(pc, p, -m) %>% mutate(sig = p<.05)

    p <- null_scores %>% 
    left_join(null_p, by = c('m', 'pc')) %>% 
    ggplot() + 
    facet_grid(factor(m, levels=unique(null_p$m))~pc, switch='y') +
    geom_density(aes(score, alpha=sig), color=mycolors[5], fill=mycolors[4]) +
    scale_alpha_manual(values=c(.2,1), guide='none') +
    geom_point(data=true_scores, aes(x=score, y=0), color=mycolors[1], size=5)
    
    if (how == 'mean') {
        p <- p + 
        geom_text(data=true_scores, aes(label=paste0('s=',round(score,3)), x=0.005, y=200), size=8, hjust=0) +
        geom_text(data=null_p, aes(label=paste0('p=',round(p,3)), x=0.005, y=100), size=8, hjust=0)
    } else {
        p <- p + 
        geom_text(data=true_scores, aes(label=paste0('s=',round(score,2)), x=2000, y=0.0005), size=8, hjust=0) +
        geom_text(data=null_p, aes(label=paste0('p=',round(p,2)), x=2000, y=0.0012), size=8, hjust=0)
    }
        
    p +
    theme_void() + 
    theme(
        strip.text.y.left = element_text(size=30, angle=0), 
        strip.text.x = element_text(size=30, angle=0), 
        text = element_text(size=30))
}