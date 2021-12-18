
source("../code/brainPlots.R")

plot_ahba_bs_brains <- function(hcp_scores, corrs) {
    corrs <- as.character(round(corrs, 2))

    plot_hcp(hcp_scores, three=T) + 
    guides(fill=F) +
    # ggtitle("AHBA PCs with DS = 0.9 (1.6k genes) // Brainspan from rnaseq, with normalization") +
    # Hacky code to add correlations to plot
    facet_grid(component~version, 
               labeller=labeller(component=c(
                   `PC1`=corrs[1],`PC2`=corrs[2],`PC3`=corrs[3]))) & 
    theme(text=element_text(size=30))
    
    }

plot_ahba_bs_scatter <- function(cortex_scores, corrs) {
    ggplot(cortex_scores) + facet_grid(.~PC) +
    geom_abline(linetype=3) +
    geom_point(aes(AHBA_mean, Brainspan, color=AHBA_mean), size=5) +
    scale_color_gradientn(colors=rev(brewer.rdbu(100)), name='') +
    guides(color=guide_colorbar(barwidth=10)) +
    geom_text(aes(AHBA_mean, Brainspan, label=cortex), size=8, hjust=0, vjust=1.2) +
    geom_text(data=data.frame(x=-1.5, y=2, PC=c('PC1','PC2','PC3'), 
                              label=paste('corr =', round(corrs[1:3],2))), aes(x=x,y=y,label=label), size=10) +
    coord_fixed() +
    # ggtitle('Brainspan-AHBA cortex correlations, AHBA DS=0.5') +
    theme_minimal() + theme(text=element_text(size=30), legend.position='bottom')
    }

library(ggrepel)
library(lemon)

plot_bs_scatter_for_ohbm <- function(cortex_scores, corrs) {
    lim <- 2.8
    
    # ggplot(cortex_scores %>% filter(PC!='PC1')) + 
    ggplot(cortex_scores) + 
    facet_rep_grid(.~PC, repeat.tick.labels=T) +
    xlim(c(-lim,lim)) + ylim(c(-lim,lim)) +
    geom_abline(linetype=3) +
    geom_point(aes(AHBA_mean, Brainspan, color=cortex), size=5) +
    scale_color_manual(values=rev(brewer.rdylbu(13)[-c(6,7,8)]), name='') +
    # guides(color=guide_colorbar(barwidth=10)) +
    # geom_text(aes(AHBA_mean, Brainspan, label=cortex), 
    #                 size=7, hjust=0, vjust=1.2) +
    geom_text(
        data=data.frame(
            x=0, y=2.5, 
            PC=c('PC1','PC2','PC3'), 
            label=paste0('r = ', round(corrs[1:3],2))), 
        aes(x=x,y=y,label=label), size=10
    ) +
    coord_fixed() +
    guides(color=guide_legend(ncol=2)) +
    # ggtitle('Brainspan-AHBA cortex correlations, AHBA DS=0.5') +
    xlab('AHBA Region PC Score') + ylab('BrainSpan\nRegion PC Score') +
    theme_void() + 
    annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf) +
    annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf) +
    theme(axis.title = element_text(),
          panel.spacing = unit(3, 'lines'),
          axis.title.y = element_text(angle=90),
          aspect.ratio=1) +
    theme(text=element_text(size=30), legend.position=c(.5,-.3))
}