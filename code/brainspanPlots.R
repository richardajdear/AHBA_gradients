# Brainspan plots

library(pals)
library(shades)
# mycolors = brewer.rdylbu(5)[c(1,2,5)]
mycolors=brewer.set1(4)[c(1,4,2)]
# mycolors = c(brewer.rdylbu(6)[1:3],brewer.rdylbu(5)[4:5])
# mycolors = saturation(c(brewer.rdylbu(5)[1:2],brewer.rdylbu(6)[4:6]), delta(.2))
# mycolors2 = c(mycolors[1:2], saturation(mycolors[3], delta(.4)))
# mycolors3 = c(mycolors[1], 
#               saturation(mycolors[3], delta(.4)), 
#               saturation(mycolors[4], delta(.4))) 


plot_bs_hcp_mapping <- function(hcp_bs_mapping,title='',xlab='') {
    
    n_colors <- hcp_bs_mapping$structure_name %>% n_distinct
    cols = cols25(n_colors)

    hcp_bs_mapping <- hcp_bs_mapping %>% 
        mutate(structure_name = replace_na(structure_name, '__NA__')) %>%
        # mutate_at(vars(structure_name), ~ factor(., levels=unique(.))) %>% 
        rename(region=label)
        
    
    ggplot(hcp_bs_mapping) + 
    geom_brain(atlas=glasser, hemi='left', mapping=aes(fill=structure_name, geometry=geometry, hemi=hemi, side=side, type=type), position=position_brain(side+hemi~.)) +
    # guides(fill=guide_legend('')) +
    xlab(xlab) +
    scale_fill_manual(values=c('white', cols), guide='none') +
    theme_void() + 
    theme(text=element_text(size=20), legend.position='none',
          axis.title.x=element_text(size=20, color='grey7'),
         plot.title=element_text(size=20, vjust=-1, hjust=.5)) +
    ggtitle(title)
                         # legend.position=c(1.3,.5))
}

plot_bs_dk_mapping <- function(dk_bs_mapping,title='',xlab='') {
    
    n_colors <- dk_bs_mapping$structure_name %>% n_distinct
    cols = cols25(n_colors)

    dk_bs_mapping <- dk_bs_mapping %>% 
        # mutate_at(vars(cortex), ~ factor(., levels=unique(.))) %>%
        mutate(structure_name = replace_na(structure_name, '__NA__'))
    
    ggplot(dk_bs_mapping) + 
    geom_brain(atlas=dk, hemi='left', mapping=aes(fill=structure_name, geometry=geometry, hemi=hemi, side=side, type=type)) +
    # guides(fill=guide_legend('')) +
    xlab(xlab) +
    scale_fill_manual(values=c('white', cols), guide='none') +
    theme_void() + 
    theme(text=element_text(size=20), legend.position='none',
          axis.title.x=element_text(size=20, color='grey7'),
         plot.title=element_text(size=20, vjust=-1, hjust=.5)) +
    ggtitle(title)
                         # legend.position=c(1.3,.5))
}

plot_bs_scores_corr <- function(bs_scores_corr, title="", xint='Birth-3 yrs', rotate=F) {
    g <- bs_scores_corr %>% 
    mutate_at(vars(age), ~ factor(., levels=unique(.))) %>%
    ggplot() + 
    geom_hline(yintercept=0, color='grey') + 
    geom_vline(xintercept=xint, color='grey') +
    geom_line(aes(x=age, y=corr, color=G, group=G), size=1, alpha=1) + 
    geom_point(aes(x=age, y=corr, color=G), size=5) + 
    xlab("") + ylab("AHBA-BrainSpan correlation") +
    scale_color_manual(values=mycolors) +
    scale_y_continuous(limits=c(0,1), breaks=seq(0,1,.2)) +
    # scale_color_manual(values=brewer.rdylbu(4)) +
    ggtitle(title) +
    theme_minimal() + 
    theme(panel.grid.minor = element_blank(),
          #axis.title.y=element_text(angle=0, vjust=0.5, hjust=1),
          axis.text=element_text(size=16, color='grey7'),
          # legend.position = c(.95, .9),
          # legend.position=element_blank(),
          legend.title = element_blank()
         )
    
    if (rotate) {
        g + theme(text=element_text(size=8), axis.text.x=element_text(size=16, angle=30, hjust=.5))
    } else {
        g #+ theme(text=element_text(size=20), axis.text.x=element_text(size=16, hjust=.5))
    }
}




plot_ahba_bs_corr <- function(df) {
    p <- df %>% mutate(
        x=recode(x, `0`='PC1',`1`='PC2',`2`='PC3',`3`='PC4',`4`='PC5'),
        y=recode(y, `0`='PC1',`1`='PC2',`2`='PC3',`3`='PC4',`4`='PC5')
    ) %>% 
    ggplot() +
    geom_tile(aes(x,y, fill=corr)) +
    geom_text(aes(x,y, label=sprintf("%0.2f", round(corr, digits = 2))), size=8) +
    scale_fill_gradientn(colours=brewer.rdbu(100)[20:80], limits=c(-1,1), guide='colourbar') +
    guides(fill=guide_colourbar(title='Corr.', barheight=10)) +
    theme_minimal() + 
    coord_fixed() +
    xlab('AHBA') + ylab('Brainspan')
}


plot_ahba_bs_scatter <- function(both_scores, corrs, facet='h', size=4) {
    # lim <- 2.8
    
    g <- ggplot(both_scores, aes(AHBA, Brainspan)) + 
    # xlim(c(-lim,lim)) + ylim(c(-lim,lim)) +
    geom_point(aes(color=AHBA), size=size) +
    geom_smooth(method='lm', linetype=1, se=F, color='grey') +
    scale_color_gradientn(colors=rev(brewer.rdbu(100)), guide='none') +
    #scale_color_manual(values=cols25(), guide='none') +
    geom_text(data=data.frame(G=names(corrs), 
                              label=paste("r =", round(corrs, 2)))
              , aes(x=-.5,y=1, label=label), size=7
    ) +
    coord_fixed() +
    xlab('AHBA Score') + ylab('BrainSpan Score') +
    scale_y_continuous(breaks=0) +
    scale_x_continuous(breaks=0) +
    theme_minimal() + 
    theme(panel.grid.minor=element_blank(),
          axis.text = element_blank(),
          strip.text.x = element_blank(),
          strip.text.y = element_text(size=20, margin=margin(l=0, r=0)),          
          axis.title.x = element_text(size=20, angle=0),
          axis.title.y = element_text(size=20, angle=0, vjust=0.5),
          aspect.ratio=1
         )
    
    if (facet=='v') {
        g + facet_rep_grid(G~., repeat.tick.labels=T)
    } else {
        g + facet_rep_grid(.~G, repeat.tick.labels=T)
    }
}




plot_hcp_bs_mapping <- function(hcp_bs_mapping) {
    hcp_bs_mapping <- hcp_bs_mapping %>% 
        mutate_at(vars(structure_name), ~ factor(., levels=unique(.)))
    cols = as.character(cols25())
    
    
    g0 <- ggplot(hcp_bs_mapping) + 
    geom_brain(atlas=glasser, mapping=aes(fill=structure_name, geometry=geometry, hemi=hemi, side=side, type=type)) +
    guides(fill=guide_legend('')) +
    scale_fill_manual(values=cols) +
    theme_void() + ggtitle('HCP') + theme(text=element_text(size=20), plot.title=element_text(vjust=-1))

    g1 <- ggplot(hcp_bs_mapping) + 
    geom_brain(atlas=glasser, mapping=aes(fill=structure_name, geometry=geometry, hemi=hemi, side=side, type=type)) +
    guides(fill=guide_legend('')) +
    scale_fill_manual(values=cols) +
    theme_void() + ggtitle('Brainspan regions mapped to HCP') + theme(text=element_text(size=20), plot.title=element_text(vjust=-1))

    g0/g1
}




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

# plot_ahba_bs_scatter <- function(cortex_scores, corrs) {
#     ggplot(cortex_scores) + facet_grid(.~PC) +
#     geom_abline(linetype=3) +
#     geom_point(aes(AHBA_mean, Brainspan, color=AHBA_mean), size=5) +
#     scale_color_gradientn(colors=rev(brewer.rdbu(100)), name='') +
#     guides(color=guide_colorbar(barwidth=10)) +
#     geom_text(aes(AHBA_mean, Brainspan, label=cortex), size=8, hjust=0, vjust=1.2) +
#     geom_text(data=data.frame(x=-1.5, y=2, PC=c('PC1','PC2','PC3'), 
#                               label=paste('corr =', round(corrs[1:3],2))), aes(x=x,y=y,label=label), size=10) +
#     coord_fixed() +
#     # ggtitle('Brainspan-AHBA cortex correlations, AHBA DS=0.5') +
#     theme_minimal() + theme(text=element_text(size=30), legend.position='bottom')
#     }

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
    scale_color_manual(values=rev(brewer.rdylbu(14)[-c(6,7,8)]), name='') +
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