# R functions for other plots
suppressMessages(library(scales))
library(patchwork)
library(ggtext)
suppressMessages(library(lemon))
library(pals)
library(shades)
suppressMessages(library(tidyverse))

# mycolors = c(brewer.rdylbu(6)[1:3],brewer.rdylbu(5)[4:5])
# mycolorscale = brightness(rev(brewer.rdbu(100)[15:85]), delta(-.2))
# mycolors = brewer.rdylbu(5)[c(1,2,5)]
mycolors = brewer.set1(5)[c(1,4,2)]
mycolorscale = brightness(rev(brewer.rdbu(100)), delta(-.1))


plot_triplets_raster <- function(triplets_raster, n_components=3) {
    triplets_raster <- triplets_raster %>% 
    mutate(method=factor(method, ordered=T, levels=c('pca','dm'),
            labels=c('PCA', 'DME'))) %>%
    filter(component <= n_components) %>%
    mutate(component=factor(component, ordered=T, levels=seq(1,5),
            labels=c('PC1\n(G1)','PC2\n(G2)','PC3\n(G3)','PC4\n(G4)','PC5\n(G5)')))

    highlight = triplets_raster %>% 
    filter(
            ((method=='PCA') & (gene_filter==0.0) & (donors_filter==1)) |
            ((method=='PCA') & (gene_filter==0.7) & (donors_filter==3)) |
            ((method=='DME') & (gene_filter==0.5) & (donors_filter==3)) |
            ((method=='DME') & (gene_filter==0.8) & (donors_filter==1))
    ) %>% 
    arrange(component, method, gene_filter, donors_filter) %>%
    mutate(label=rep(c('(1)','(2)','(3)','(4)'), n_components))

    triplets_raster %>%
    ggplot(aes(x=gene_filter, y=donors_filter, fill=corr_abs)) + 
    geom_tile() +
    geom_text(aes(label=label), data=highlight, color='grey7', family='Calibri', size=8) +
    # geom_text(aes(label='★'), data=triplets_raster %>% filter(corr_abs>0.6), color='gray30', size=5) +
    # geom_tile(aes(x=gene_filter, y=donors_filter), fill=NA, color='white', size=1,
    #             data=triplets_raster %>% filter(corr_abs > 0.6)
    # ) +
    facet_grid(component~method, switch='y') +
    scale_fill_gradientn(
            colors=rev(brewer.spectral(100)),
            limits=c(0,1), breaks=seq(0,1,.2), 
            labels=c('0.0', '0.2', '0.4', '0.6', '0.8', '1.0'),
            name='Median \ntriplet \ncorrelation'
    ) +
    scale_x_continuous(
            breaks=c(0,.5,.9), minor_breaks=c(.3,.7),
            labels=c('All\ngenes','Top 50%\ngenes','Top 10%\ngenes'),
            name=''
    ) +
    scale_y_discrete(position='right',
            labels=c('All regions','2+ donor regions','3+ donor regions'),
            name=''
    ) +
    guides(fill=guide_colorbar(barheight=20)) +
    theme_minimal() +
    theme(
        text=element_text(size=22),
        axis.text = element_text(size=22, color='grey7', family='Calibri'),
        aspect.ratio=1/3.8,
        panel.grid=element_blank(),
        legend.title=element_text(vjust=3),
        # legend.position='bottom',
        strip.text.y.left=element_text(angle=0)
    )
}


plot_dist_donors_hcp <- function(df_donors) {
    df_donors %>%
    ggplot() + 
    geom_histogram(aes(x=count), alpha=1, binwidth=1, color='white', size=3, fill=brewer.blues(10)[4]) + 
    # scale_fill_manual(values=brewer.blues(10)[c(10,4)], name='', guide='none') +
    geom_vline(xintercept=2.5, linetype=2) +
    scale_x_continuous(breaks=seq(0,6,1)) +
    theme_minimal() + ylab('Regions') + xlab('Donors sampled') + 
    theme(legend.position=c(.2,.8), panel.grid=element_blank())    
}


plot_ds_dist_hcp <- function(df_stability) {
    ggplot(df_stability) + 
    geom_density(aes(x=ds), size=1, color=brewer.blues(10)[6], fill=brewer.blues(10)[4]) +
    geom_vline(xintercept=0.386, linetype=2) +
    annotate(x=.39,y=2.5,geom='text',label='Top 10%', hjust=-0.05, size=7) +
    ylab('Density') + xlab('Gene diff. stability') +
    theme_minimal() +
    theme(panel.grid=element_blank())
}


plot_weights_ds <- function(weights_ds, ncol=3) {
    df <- weights_ds %>% select(G1:G3, ds) %>% gather(G, weight, -ds)

    ggplot(df, aes(x=weight, y=ds)) + 
        facet_wrap(~G, ncol=ncol) +
        # geom_point(color=mycolors[5], alpha=.5, size=1) +
        geom_point(aes(color=weight), size=.01) +
        scale_color_gradientn(colors=mycolorscale, guide='none') +
    # scale_y_continuous(breaks=c(-.02,0,.02)) + 
    # scale_x_continuous(limits=c(0.1,.8), breaks=c(.1,.8)) +   
        xlab('Gene weight') + ylab("Gene DS") +
        theme_minimal() + 
        theme(
            panel.grid=element_blank(), 
            axis.text = element_blank(),
            strip.text = element_blank(),
            # axis.title.y = element_text(angle=0, vjust=0.5, hjust=1),
            aspect.ratio = 1
        )
}

plot_weights_dist <- function(weights_ds, ncol=1) {
    df <- weights_ds %>% select(G1:G3, ds) %>% gather(G, weight, -ds)

    ggplot(df) +
        facet_wrap(.~G, ncol=ncol) +
        geom_density(aes(weight), size=.2, alpha=.8, color=mycolors[3]) +
        # geom_density(size=1, alpha=.8, aes(color=stat(x))) +
        # stat_density_ridges(geom = "density_ridges_gradient", calc_ecdf = TRUE, scale=.8, size=.2, fill='transparent',
                       # aes(x=weight, y=1, color=stat(ecdf))) +
        # scale_color_gradientn(colors=rev(brewer.rdbu(100)[15:85]), guide='none') +
    coord_flip() +
        theme_void() +
        theme(strip.text=element_blank())
}


plot_triplets_v2 <- function(triplets_plot_v2, facet='h', ncol=4,
                             threshold=0.6, colors=mycolors, smooth=T,
                             xlab='% of genes retained by differential stability filter') {

    p <- triplets_plot_v2 %>%
    filter(component <= 3) %>%
    ggplot(aes(x=DS, y=corr_abs)) +
        geom_hline(yintercept=threshold, color=mycolors[1], linetype='dashed') +
        geom_jitter(aes(color=component), width=0.01) +
        # geom_ribbon(aes(group=component, fill=component, ymin=ymin, ymax=ymax), alpha=.3) +
        scale_color_manual(values=colors, name='',labels=c('G1','G2','G3')) +
        scale_fill_manual(values=colors, name='',labels=c('G1','G2','G3')) +
        scale_x_continuous(
            breaks=c(0,.5,.9), minor_breaks=c(.3,.7),
            labels=c('100%','50%','10%'),
            name=xlab
        ) +
        scale_y_continuous(breaks=seq(0,1,.2), 
                           # position='right',
                           name='Disjoint triplet correlation') +
        theme_minimal() + 
        theme(
            panel.grid.minor.y=element_blank(),
            axis.text=element_text(size=20, color='grey7', family='Calibri'),
              # legend.position=c(.9,.3),
              legend.position='bottom',
              # axis.text.y = element_text(margin=margin(r=-1)),
              panel.spacing=unit(2,'lines')
             )

    if (smooth) {
        p <- p + geom_smooth(aes(group=component, fill=component, color=component), method='loess', span=1)
    } else {
        p <- p + stat_summary(geom='line', fun.y=median, aes(group=component, fill=component, color=component))
    }
    
    if (facet=='h') {
        p + facet_grid(.~method)
    } else if (facet=='w') {
        p + facet_wrap(~method, ncol=ncol)
    }
    
}


plot_triplets_simple <- function(triplets_plot_v2, facet='h', colors=mycolors, size=1, alpha=.8) {

    p <- triplets_plot_v2 %>%
    filter(DS == 0) %>%
    filter(component <= 3) %>%
    ggplot(aes(x=component, y=corr_abs)) +
        geom_jitter(aes(color=component), width=0.2, size=size, alpha=alpha) +
        scale_color_manual(values=colors, name='',labels=c('G1','G2','G3')) +
        scale_fill_manual(values=colors, name='',labels=c('G1','G2','G3')) +
        scale_x_discrete(labels=c('PC1','PC2','PC3')) +
        scale_y_continuous(breaks=seq(0,1,.2), 
                           # position='right',
                           name='Disjoint triplet correlation') +
        theme_minimal() + 
        theme(panel.grid.minor=element_blank(),
              axis.text=element_text(size=20, color='grey7', family='Calibri'),
              # legend.position=c(.9,.3),
              legend.position='bottom',
              # axis.text.y = element_text(margin=margin(r=-1)),
              panel.spacing=unit(2,'lines')
             )
    
    if (facet=='h') {
        p + facet_grid(.~method)
    } else if (facet=='w') {
        p + facet_wrap(~method)
    }
}


plot_scatter_corrs <- function(versions_scatter, xlab='', ylab='', size=6) {
    corrs <- versions_scatter %>% group_by(G) %>% 
    summarize(cor(x, y, use='p')) %>%
    rename_with(function(x) c('G','r'), everything()) %>%
    mutate(label=paste('r =', round(r,3), '***'))

    versions_scatter %>%
    ggplot(aes(x=x, y=y)) + 
    facet_grid(G~.) +
    geom_point(alpha=.2) +
    geom_smooth(method='lm', aes(color=G), se=F) +
    geom_text(data=corrs, x=-2, y=2.3, aes(label=label), size=size, hjust=0) +
    xlab(xlab) + ylab(ylab) +
    scale_color_manual(values=mycolors) +
    guides(color='none') +
    coord_cartesian(clip='off') +
    theme_minimal() + 
    theme(
        # panel.border=element_rect(fill='transparent', color='grey'),
        panel.grid=element_blank(),
        strip.text.y=element_text(angle=0),
        axis.text=element_blank(),
        # axis.title.y=element_text(angle=0, vjust=.5),
        aspect.ratio=1
    )
}

plot_weight_corrs <- function(weight_corrs, spacing=8, ylab='') {
    weight_corrs %>% 
    mutate_at(vars(version), ~ factor(., levels=unique(.))) %>% 
    mutate_at(vars(y), ~ factor(., levels=rev(unique(.)))) %>% 
    ggplot() +
    facet_rep_grid(.~version, repeat.tick.labels=T) + #, switch='x') +
    geom_tile(aes(x,y, fill=corr)) +
    geom_text(aes(x,y, label=sprintf("%0.2f", round(corr, digits = 2))), size=10) +
    scale_fill_gradientn(colours=rev(brewer.rdbu(100)[20:80]), limits=c(-1,1), guide='colourbar') +
    guides(fill=guide_colourbar(title='Gene weight correlation', barwidth=20)) +
    scale_x_discrete(position = "top") +
    theme_minimal() + 
    theme(panel.spacing=unit(spacing,'lines'), 
          # axis.title.y=element_text(angle=0),
          strip.text.y=element_text(size=36),
          strip.text.x=element_text(size=36),
          strip.placement='outside',
          legend.position='bottom',
          legend.title = element_text(vjust=.5)
         ) +
    coord_fixed() + 
    xlab('') +
    ylab('')
}



plot_corrs <- function(df, facetting='h', xlab='', ylab='', size=8) {
    p <- df %>% 
    # filter(x<3, y<3) %>%
    mutate(
        x=recode(x, `0`='G1',`1`='G2',`2`='G3',`3`='G4',`4`='G5'),
        y=recode(y, `0`='G1',`1`='G2',`2`='G3',`3`='G4',`4`='G5')
    ) %>% 
    mutate_at(vars(version), ~ factor(., levels=unique(.))) %>% 
    ggplot() +
    geom_tile(aes(x,y, fill=corr)) +
    geom_text(aes(x,y, label=sprintf("%0.2f", round(corr, digits = 2))), size=size) +
#     scale_fill_gradient2(low=muted('red'),mid='white',high=muted('blue'), limits=c(-1,1)) +
    scale_fill_gradientn(colours=rev(brewer.rdbu(100)[20:80]), limits=c(-1,1), guide='colourbar') +
    scale_y_discrete(limits=rev) +
    guides(fill=guide_colourbar(title='Corr.', barwidth=20)) +
    theme_minimal() + 
    theme(panel.spacing=unit(4,'lines'), 
          text=element_text(size=30), 
          strip.text.y=element_text(size=30),
          strip.text.x=element_text(size=30),
          strip.placement='outside',
          legend.position='bottom',
          legend.title.align = 0.5
         ) +
    coord_fixed() +
    xlab(xlab) + ylab(ylab) + ggtitle('')
    
    if (facetting=='v') {
        p + facet_rep_grid(version~., repeat.tick.labels=T)
    } else if (facetting=='h') {
        p + facet_rep_grid(.~version, repeat.tick.labels=T)
    } else if (facetting=='w') {
        p + facet_rep_wrap(~version, repeat.tick.labels=T)
    } else {
        p + facet_rep_grid(how~version, repeat.tick.labels=T, switch='y')
    }        
}


plot_atlas_score_dots <- function(atlas_mean_scores_sig_colors, size=7) {
    atlas_mean_scores_sig_colors %>% 
    mutate(sig = case_when(q<0.001 ~ '***', q<0.01 ~ '**', q<0.05 ~ '*', TRUE ~ '')) %>%
    ggplot(aes(x=true_mean, y=G)) +
      facet_wrap(~atlas, ncol=3) +
      geom_point(aes(fill=I(color)), color='grey50', size=size, shape=21, alpha=.8) +
      geom_text(aes(label=sig), size=8, vjust=-.5) +
      scale_y_discrete(limits=rev, name='') +
      scale_x_continuous(breaks=c(-1,0,1), name='Mean axis z-score in atlas region') +
      coord_cartesian(clip='off') +
      theme_minimal() + 
      theme(panel.grid.minor=element_blank(),
            panel.spacing=unit(2,'lines'),
            axis.text=element_text(size=22, color='grey7', family='Calibri')
           )
}


plot_atlas_violins <- function(regions, sig, classes=Mesulam, 
                               classcolors=Mesulam_colors,
                               classlabels=Mesulam_names
                              ) {
    colors <- regions %>% select( {{classes}} ,  {{classcolors}} ) %>% unique() %>% 
        arrange( {{classes}} ) %>% drop_na() %>% pull( {{classcolors}} )
    
    labels <- regions %>% select( {{classes}} ,  {{classlabels}} ) %>% unique() %>% 
        arrange( {{classes}} ) %>% drop_na() %>% pull( {{classlabels}} )
    regions_labels <- data.frame(x=seq(1,4), label=labels)

    sig <- mesulam_scores_sig %>% mutate(m=as.character(label)) %>% 
            filter(m != 'NaN', m != 8, m != 0) %>%
            mutate(sig = case_when(q<0.001 ~ '***', q<0.01 ~ '**', q<0.05 ~ '*', TRUE ~ ''))

    
    regions %>% 
    mutate(m=factor( {{classes}} )) %>%
    filter(m != 'NaN', m != 8, m != 0) %>%
    select(G1:G3, m) %>%
    gather(G, score, -m) %>%
    group_by(G, m) %>%
    mutate(med = median(score, na.rm=T), mean = mean(score, na.rm=T)) %>%
    ggplot(aes(x=m, y=score)) + 
    facet_grid(.~G, switch='x') +
    # geom_hline(yintercept=0, color='grey50', size=.2) +
    geom_violin(aes(fill=m), width=.7, color='grey50') +
    # geom_point(aes(y=med), color='red') +
    # geom_point(aes(y=mean), color='blue') +
    geom_text(data=sig, aes(x=m, y=true_mean, label=sig), size=8, vjust=.75, hjust=.5) +
    # geom_boxplot(aes(m, score), width=.1) +
    scale_fill_manual(values=colors, guide='none') +
    scale_y_continuous(breaks=0) +
    # scale_fill_gradientn(colors=rev(brewer.rdbu(100)), limits=c(-2,2)) +
    coord_flip(ylim=c(-2.5,2.5), clip='off') +
    # geom_text(aes(x=x,y=-4,label=label),size=8,data=df_labels,angle=35,hjust=1) +
    theme_minimal() + 
    theme(panel.spacing=unit(2,'lines'),
          axis.text=element_blank(),
          axis.title=element_blank(),
          panel.grid.minor=element_blank(),
          panel.grid.major.x=element_line(color='gray50'),
          strip.placement='outside',
          plot.margin = unit(c(t=1,r=1,b=5,l=1), "lines"),
          text=element_text(size=30),
          aspect.ratio=1
         )
}


plot_violins_old <- function(df, classes=vonEconomo, 
                               classcolors=vonEconomo.colors,
                               classlabels=vonEconomo_labels
                              ) {
    colors <- df %>% select( {{classes}} ,  {{classcolors}} ) %>% unique() %>% 
        arrange( {{classes}} ) %>% drop_na() %>% pull( {{classcolors}} )
    
    labels <- df %>% select( {{classes}} ,  {{classlabels}} ) %>% unique() %>% 
        arrange( {{classes}} ) %>% drop_na() %>% pull( {{classlabels}} )
    df_labels <- data.frame(x=seq(1,7), label=labels)
    
    df %>% 
    mutate(m=factor( {{classes}} )) %>%
    filter(m != 'NaN', m != 8) %>%
    select(G1:G3,m) %>%
    gather(G, score, -m) %>%
    group_by(G, m) %>%
    mutate(Mean = mean(score)) %>%
    ggplot() + facet_grid(~G) +
    geom_hline(yintercept=0) +
    geom_violin(aes(m, score, fill=Mean)) +
    geom_boxplot(aes(m, score), width=.1) +
    scale_fill_gradientn(colors=rev(brewer.rdbu(100)), limits=c(-2,2)) +
    coord_cartesian(ylim=c(-2.5,2.5), clip='off') +
    geom_text(aes(x=x,y=-4,label=label),size=8,data=df_labels,angle=35,hjust=1) +
    theme_minimal() + 
    theme(panel.spacing=unit(2,'lines'),
          axis.text.x = element_text(color=colors, face='bold', size=30),
          plot.margin = unit(c(t=1,r=1,b=5,l=1), "lines"),
          text=element_text(size=30)
         ) +
    xlab('') + ylab('Region scores')
}


plot_xyz <- function(df, component=G1, title='G1', colors=brewer.set1(5)) {
    p <- df %>%
    mutate(colors=Lobe) %>%
#     select(colors, G1:G3, X=pos_all1, Y=pos_all2, Z=pos_all3) %>%
    dplyr::select(colors, {{ component }}, X=x, Y=y, Z=z) %>%
    gather(G,score,-(X:Z), -colors) %>%
    gather(dim,coord,-G,-score, -colors) %>%
    ggplot(aes(coord,score)) + 
    facet_grid(.~dim, scales='free_x') +
    geom_point(aes(color=colors), size=1, alpha=.5) +
    geom_smooth(method='lm', color='black', se=F) +
    scale_color_manual(values=colors, guide=guide_legend(ncol=1,title='Lobe')) +
#     coord_fixed() +
    theme_minimal() + 
    theme(panel.grid=element_blank(),
          plot.title=element_text(size=30,hjust=0.5,vjust=-4),
          axis.text=element_blank(),
          aspect.ratio=1) +
    xlab('XYZ coords') + 
    ggtitle(title)
    
    if (title=="G1") {
        p + ylab("Region scores")
    } else {
        p + ylab('')
    }
}



# plot_coefs_ds <- function(df_coefs_ds, facet='h') {
#     g <- coefs_ds %>% rownames_to_column %>% rename(gene=rowname) %>% 
#     gather(PC, coef, -gene, -DS) %>% 
#     ggplot(aes(coef, DS)) + 
#     geom_point(color=brewer.rdbu(100)[80], alpha=.3, size=.2) +
#     xlab('Gene Weight') + ylab("Differential Stability") +
#     theme_void() +
#     # annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf) +
#     # annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf) +
#     theme(axis.title.x = element_text(),
#           axis.title.y = element_text(angle=90),
#           strip.text.y = element_blank(),
#           panel.spacing = unit(3, 'lines'),
#           aspect.ratio=1)
    
#     if (facet=='v') {
#         g + facet_grid(PC~.)
#     } else {
#         g + facet_grid(.~PC)
#     }
# }








plot_var_exp <- function(df_var) {
    ggplot(df_var) + 
    geom_line(aes(PC, var, color=which, group=which),size=1) + 
    scale_color_manual(values=cols25(10), name='') +
    theme_minimal() +
    theme(panel.grid.minor = element_blank()) +
    theme(legend.position=c(.6,.8), legend.title=element_blank()) +
    ylab('Variance Explained') + xlab('') #+ ylim(c(0,.4))
}





plot_dist <- function(df_dist) {
    ggplot(df_dist) +
    geom_histogram(aes(x=count, fill=atlas), alpha=1, binwidth=1, position='identity') + 
    theme_minimal() + ylab('Regions') + xlab('Samples/Region') +
    theme(legend.position=c(.8,.8), panel.grid=element_blank()) +
#     scale_fill_manual(values=brewer.ylorrd(10)[c(9,4)], name='Parcellation')
    scale_fill_manual(values=brewer.blues(10)[c(10,4)], name='')
#     scale_fill_manual(values=brewer.rdbu(10)[c(3,8)], name='Parcellation')
}

plot_dist_donors <- function(df_donors) {
    ggplot(df_donors) + 
    geom_histogram(aes(x=count, fill=atlas), alpha=1, binwidth=1, position='dodge2', color='white', size=3) + 
    scale_fill_manual(values=brewer.blues(10)[c(10,4)], name='') +
#     scale_fill_manual(values=mycolors[c(1,3)], name='Parcellation') +
    scale_x_continuous(breaks=seq(0,6,1)) +
    theme_minimal() + ylab('Regions') + xlab('Donors/Region') + 
    theme(legend.position=c(.2,.8), panel.grid=element_blank())    
}




plot_exp_corr <- function(df_genes) {
    ggplot(df_genes) + 
    facet_grid(.~atlas) +
    geom_histogram(aes(corr, fill=atlas), binwidth=0.01) +
    scale_fill_manual(values=brewer.blues(10)[c(9,4)], guide='none') +
    scale_x_continuous(limits=c(0,1), breaks=c(0,0.5,1), labels=c(0,0.5,1)) + 
    scale_y_continuous(breaks=c(0,2000,4000), limits=c(0,4600)) + 
    theme_minimal() +
    theme(panel.grid.minor=element_blank(),
          panel.spacing=unit(1,'lines')
         ) +
    xlab('Correlation') + ylab('# genes')
}

plot_matching <- function(plot_df) {
    ggplot(plot_df %>% mutate_at(vars(version), ~ factor(., levels=unique(.)))) + 
     facet_grid(.~version) +
     geom_col(aes(x=reorder(labels, pct_count), y=pct_count, fill=version)) + 
#      scale_fill_manual(values=c(muted('red'),muted('blue')), guide='none') +
     scale_fill_manual(values=brewer.blues(10)[c(9,4)], guide='none') +
     scale_y_continuous(limits=c(-0.05,1), breaks=seq(0,1,.5)) +
     geom_hline(yintercept=1,size=.5) +
     geom_text(aes(x=reorder(labels, pct_count), y=pct_count, label=n), hjust=1.1, size=10) + 
     coord_flip() + 
     theme_minimal() + 
     theme(panel.grid=element_blank()) +
#      theme(text=element_text(size=20)) +
     xlab('') + ylab('% of samples')
}





plot_boot_test <- function(df_boot_test) {
    ggplot(df_boot_test %>% filter(type=='coefs')) +
    facet_grid(PC~atlas, scales='free_y') +
    geom_density(aes(x=corr,color=match),size=1) +
    scale_x_reverse() +
    scale_y_continuous(breaks=function(lims) {c(0, signif(lims[2]*3/5,1))}) +
    scale_color_manual(values=brewer.rdylbu(10)[c(3,8)]) +
    guides(color=guide_legend(title=element_blank())) +
    theme_minimal() + 
    theme(panel.grid.minor=element_blank()) +
    xlab('Gene weight correlation') + ylab('Density')
}

plot_missing <- function(df) {
    df %>% 
#     mutate(type = factor(type, levels=c('Gene corr.', 'Region corr.'))) %>%
#     gather(PC, value, -type, -atlas, -missing) %>%
    filter(atlas=='HCP') %>% select(-atlas) %>% 
    gather(PC, value, -type, -missing) %>%
    ggplot() +
    facet_grid(type~.) +
#     facet_grid(type~atlas) +
    geom_line(aes(missing, value, color=PC, group=PC), size=1) +
    theme_minimal() + theme(legend.title=element_blank()) +
    theme(panel.spacing=unit(2,'lines'),
          panel.grid.minor=element_blank(),
          strip.text.y = element_text(angle=0, size=30)
         ) +
    scale_y_continuous(limits=c(0,1),breaks=seq(0,1,.2)) +
    scale_x_continuous(breaks=seq(0,.5,.1)) +
    scale_color_manual(values=mycolors) +
#     scale_color_manual(values=brewer.rdylbu(6)) +
#     scale_color_manual(values=inferno(7)[-c(1,7)]) +
    ylab('') + xlab('% regions dropped')
}

plot_noise <- function(df) {
    df %>% 
#     mutate(type = factor(type, levels=c('Gene corr.', 'Region corr.'))) %>%
    filter(atlas=='HCP') %>% select(-atlas) %>%
    gather(PC, value, -type, -beta) %>%
#     gather(PC, value, -type, -atlas, -beta) %>%
    ggplot() +
#     facet_grid(type~atlas, scales='free') +
    facet_grid(type~., scales='free') +
    geom_line(aes(beta, value, color=PC, group=PC), size=1) +
    scale_y_continuous(limits=c(0,1), breaks=seq(0,1,.2)) +
    theme_minimal() + theme(legend.title=element_blank()) +
    theme(panel.spacing=unit(2,'lines'),
          panel.grid.minor=element_blank(),
          strip.text.y = element_text(angle=0, size=30)
         ) +
#     scale_x_continuous(breaks=seq(0,5,1)) +
#     scale_color_manual(values=brewer.rdylbu(6)) +
#     scale_color_manual(values=c(brewer.rdylbu(6)[1:3],brewer.rdylbu(5)[4:5])) +
    scale_color_manual(values=mycolors) +
    ylab('') + xlab('Added noise (# std deviations)')
}

plot_robustness <- function(df_versions) {
    df_versions %>%
    mutate(p=replace(p,p=='nan','no')) %>%
    gather(PC, corr, -type, -param, -p) %>%
    mutate(corr = abs(corr)) %>%
#     mutate(highlight = factor(ifelse(corr<.8,'<0.8','>0.8'))) %>%
    mutate_at(vars(param), ~ factor(., levels=unique(.))) %>% 
    ggplot() +
    geom_tile(aes(y=p, x=PC, fill=corr), color='grey', width=.7, height=.7, size=1) +
    facet_grid(param~type, scale='free', space='free', switch='y') +
    theme_minimal() +
    theme(strip.text.y.left = element_text(angle = 0, hjust=0, size=30), 
          strip.text.x = element_text(size=30),
          axis.title.y.left = element_text(vjust = -10),
          strip.placement = "outside",
          plot.title = element_text(hjust=.5),
          legend.position='bottom'
         ) +
    scale_fill_gradientn(colors=rev(brewer.ylorrd(100)), limits=c(0,1), breaks=seq(0,1,.2)) +
#     scale_fill_viridis(option='magma', limits=c(0,1), breaks=seq(0,1,.2)) +
    scale_color_manual(values=c(brewer.ylorrd(10)[9],'grey'), guide=guide_legend(title='')) +
    guides(fill=guide_colorbar(title="Correlation", barwidth=30, reverse=T, title.hjust=-.1, title.vjust=1)) +
    ylab('') + xlab('PC')
}



plot_triplets <- function(versions_df, title="", line=0.6) {
    versions_df %>%
    ggplot(aes(component, corr_abs)) + 
    facet_grid(how~version)  +
    geom_hline(yintercept=line,color=mycolors[1],size=.5,linetype=1) +
#     geom_violin(aes(fill=component), alpha=.2) + 
    geom_boxplot(aes(fill=component), alpha=.7) +
    geom_point() +
    # geom_jitter(size=.5) + 
    theme_minimal() +
    theme(panel.grid.minor=element_blank(), panel.spacing=unit(2,'lines'),
          strip.text.y = element_text(size=30), 
          strip.text.x = element_text(size=30, hjust=0)) +
    scale_fill_manual(values=mycolors,guide='none') +
    scale_y_continuous(limits=c(0,1), breaks=seq(0,1,.2)) +
#     scale_fill_manual(values=parula(5),guide='none') +
    ggtitle(title) + xlab('PC') + ylab('Matched correlation')
}




plot_regions_exist <- function(df_exist) {
    df_exist %>%
    mutate_at(vars(region), ~ factor(., levels=unique(.))) %>% 
    # mutate_at(vars(exist_pair), ~ factor(., levels=unique(.))) %>% 
    # mutate_at(vars(triplet), ~ factor(., levels=unique(.))) %>% 
    ggplot() +
    geom_raster(aes(x=region,y=triplet,fill=exist_pair),alpha=.7, color='white',size=.1) +
    scale_fill_manual(values=c('white',rev(brewer.blues(15))),guide='none') +
#     scale_fill_manual(values=c('white',brewer.blues(10)),guide='none') +
    theme_minimal() +
    theme(axis.text=element_blank(), legend.position='bottom', panel.grid=element_blank()) +
    xlab('Region') + ylab('Triplet')
}

plot_gene_exist <- function(df_exist) {
    df_exist %>%
    mutate_at(vars(gene_symbol), ~ factor(., levels=unique(.))) %>% 
    # mutate_at(vars(exist_pair), ~ factor(., levels=unique(.))) %>% 
    # mutate_at(vars(triplet), ~ factor(., levels=unique(.))) %>% 
    ggplot() +
    geom_raster(aes(x=gene_symbol,y=triplet,fill=exist_pair),alpha=.7, color='white',size=.1) +
    scale_fill_manual(values=c('white',rev(brewer.blues(15))),guide='none') +
#     scale_fill_manual(values=c('white',brewer.blues(10)),guide='none') +
    theme_minimal() +
    theme(axis.text=element_blank(), legend.position='bottom', panel.grid=element_blank()) +
    xlab('Gene number') + ylab('Triplet')
}


plot_ds_dist <- function(df_stability) {
    ggplot(df_stability) + 
    geom_density(aes(x=DS, group=triplet), size=.2, color=mycolors[5]) +
geom_vline(xintercept=0.31, linetype=2) +
annotate(x=.32,y=2.5,geom='text',label='0.8 quantile threshold', hjust=-0.05, size=8) +
theme_minimal() + ylab('Density') +
    theme_minimal() +
    theme(panel.grid.minor=element_blank())
}

plot_gene_intersection <- function(df_gene_intersection) {
    ggplot(df_gene_intersection) + 
    facet_grid(.~triplet) +
    geom_col(aes(x=gene, fill=fill, y=values), width=.8, position='identity') +
    scale_y_discrete(labels = NULL, breaks = NULL) +
#     scale_alpha_discrete(range=c(.7,.7,1)) +
#     scale_fill_manual(values=c(mycolors[c(1,5)],'white'),guide='none') +
 scale_fill_manual(values=c(brewer.blues(10)[c(10,4)],'white'),guide='none') +
#     scale_fill_manual(values=brewer.blues(10)[c(10,4)], name='Parcellation') +
    theme_minimal() + ylab('Gene number') + xlab('Triplets') +
    theme(panel.grid.major=element_blank(), axis.text=element_blank(),
         panel.spacing=unit(2,'lines'))
}


plot_gene_intersection_2 <- function(df_gene_intersection) {
    ggplot(df_intersection) + 
    facet_grid(.~which, scales='free_x', space='free_x') +
    geom_col(aes(x=DS, fill=fill, y=values), width=.08,alpha=1, position='identity') +
    scale_y_continuous(labels = NULL, breaks = NULL) +
    scale_x_reverse(breaks=rev(seq(0,.8,.1)), labels=percent(seq(.2,1,.1),1), 
                   expand=expansion(add = c(.025,.025))) +
    scale_fill_manual(values=c(brewer.blues(10)[c(10,4)],'white'),guide='none') +
    theme_minimal() + ylab('Gene number') + xlab('DS quantile filter for six brains') +
    coord_cartesian(clip='off') +
    theme(panel.grid = element_blank()) #, panel.spacing = unit(1,'lines'))
}

plot_gene_corrs <- function(df_gene_corrs) {
    ggplot(df_gene_corrs, aes(x=DS)) + 
    geom_ribbon(aes(ymax=p95,ymin=p05),alpha=.3,fill=brewer.blues(5)[2]) +
    geom_line(aes(y=mean), color=brewer.blues(5)[5],size=1) +
    scale_x_reverse(breaks=seq(.2,1,.1), labels=percent(rev(seq(.1,.9,.1)),1), limits=c(1,.2), 'DS quantile of genes') +
    scale_y_continuous(limits=c(0,1), breaks=seq(0,1,.2), name='Mean corr. with triplets') +
    theme_minimal() +
    theme(panel.grid.minor=element_blank())   
}




plot_triplet_ds_corrs <- function(df_ds_corrs) {
    ggplot(df_ds_corrs) + 
    facet_grid(.~how) +
    geom_vline(xintercept=0.5,size=40,alpha=0.3,color=brewer.blues(5)[2]) +
    geom_boxplot(aes(x=DS, y=Corr, fill=PC, group=interaction(DS,PC)),alpha=.5) +
    stat_summary(aes(x=DS, y=Corr, group=PC, color=PC), fun=mean, geom="line", size=1) +
    scale_x_reverse(breaks=rev(seq(0,0.8,0.1)), labels=percent(seq(.2,1,.1),1)) + 
    scale_y_continuous(breaks=seq(0,1,.2)) +
    scale_fill_manual(values=mycolors) +
    scale_color_manual(values=mycolors) +
#     scale_fill_manual(values=brewer.rdylbu(10)[c(1,4,8)]) +
#     scale_color_manual(values=brewer.rdylbu(10)[c(1,4,8)]) +
    theme_minimal() + 
    theme(panel.grid.minor=element_blank()) +
    xlab('DS quantile filter for six brains') + 
    ylab('Corr. with triplets')    
}



# library(gridExtra)
# plot_violins <- function(df, classes=vonEconomo, 
#                          classcolors=vonEconomo.colors,
#                          classlabels=vonEconomo_labels
#                         ) {
#     colors <- df %>% select( {{classes}} ,  {{classcolors}} ) %>% unique() %>% 
#         arrange( {{classes}} ) %>% drop_na() %>% pull( {{classcolors}} )
    
#     labels <- df %>% select( {{classes}} ,  {{classlabels}} ) %>% unique() %>% 
#         arrange( {{classes}} ) %>% drop_na() %>% pull( {{classlabels}} )
#     df_labels <- data.frame(x=seq(1,7), label=labels)
    
#     df %>% 
#     mutate(m=factor( {{classes}} )) %>%
#     filter(m != 'NaN', m != 8) %>%
#     select(PC1:PC3,m) %>%
#     gather(PC, score, -m) %>%
#     group_by(PC, m) %>%
#     mutate(Mean = mean(score)) %>%
#     ggplot() + facet_grid(~PC) +
#     geom_violin(aes(m, score, fill=Mean)) +
#     geom_boxplot(aes(m, score), width=.1) +
# #     geom_boxplot(aes(m, score, fill=Mean)) +
#     # geom_text(aes(x=m), y=-150) +
# #     scale_fill_gradient2(high='#3A3A98',mid='white',low='#832424', limits=c(-2,2)) +
#     scale_fill_gradientn(colors=rev(brewer.rdbu(100)), limits=c(-2,2)) +
#     coord_cartesian(ylim=c(-2.5,3.5), clip='off') +
#     theme_minimal() + 
#     theme(panel.spacing=unit(2,'lines')) +
#     theme(axis.text.x = element_text(color=colors, face='bold', size=30)) + 
#     geom_text(aes(x=x,y=-4,label=label),size=8,data=df_labels,angle=35,hjust=1) +
#     theme(plot.margin = unit(c(1,1,2,1), "lines")) + # This widens the right margin
# #     annotation_custom(geom = "text", x=2, y = -3, label = 'hello', size = 6) +
#     xlab('') + ylab('Region scores')
# }



# plot_xyz <- function(df, component=PC1, title='PC1', colors=brewer.set1(5)) {
#     p <- df %>%
#     mutate(colors=Lobe) %>%
# #     select(colors, PC1:PC3, X=pos_all1, Y=pos_all2, Z=pos_all3) %>%
#     select(colors, {{ component }}, X=x, Y=y, Z=z) %>%
#     gather(PC,score,-(X:Z), -colors) %>%
#     gather(dim,coord,-PC,-score, -colors) %>%
#     ggplot(aes(coord,score)) + 
#     facet_grid(.~dim, scales='free_x') +
#     geom_point(aes(color=colors), size=1, alpha=.5) +
#     geom_smooth(method='lm', color='black', se=F) +
#     scale_color_manual(values=colors, guide=guide_legend(ncol=1,title='Lobe')) +
# #     coord_fixed() +
#     theme_minimal() + 
#     theme(panel.grid=element_blank(),
#           plot.title=element_text(size=30,hjust=0.5,vjust=-4),
#           axis.text=element_blank(),
#           aspect.ratio=1) +
#     xlab('XYZ coords') + 
#     ggtitle(title)
    
#     if (title=="PC1") {
#         p + ylab("Region scores")
#     } else {
#         p + ylab('')
#     }
# }