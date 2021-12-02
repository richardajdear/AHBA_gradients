
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
    geom_text(data=data.frame(x=-1.5, y=2, PC=c('PC1','PC2','PC3'), label=paste('corr =', round(corrs[1:3],2))), aes(x=x,y=y,label=label), size=10) +
    coord_fixed() +
    # ggtitle('Brainspan-AHBA cortex correlations, AHBA DS=0.5') +
    theme_minimal() + theme(text=element_text(size=30), legend.position='bottom')
    }