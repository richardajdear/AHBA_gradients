library(ggseg)
# library(ggsegExtra)
library(ggsegGlasser)
# library(ggrepel)
library(ggh4x)

source("../code/brainPlots.R")


mycolors = c(brewer.rdylbu(6)[1:3],brewer.rdylbu(5)[4:5])
# mycolors = saturation(c(brewer.rdylbu(5)[1:2],brewer.rdylbu(6)[4:6]), delta(.2))
# mycolors2 = c(mycolors[1:2], saturation(mycolors[3], delta(.4)))
# mycolors3 = c(mycolors[1], 
#               saturation(mycolors[3], delta(.4)), 
#               saturation(mycolors[4], delta(.4))) 





plot_maps <- function(maps, title="", ncol=3, facet='w', spacing=0,
                      colors=rev(brewer.rdbu(100)), colorscale='symmetric', 
                      name='Z-score', labels='none') {
    df <- maps %>%
        rownames_to_column %>%
        mutate(region = recode(rowname,'7Pl'='7PL')) %>% 
        select(-rowname) %>%
        gather('map', 'value', -region) %>%
        mutate_at(vars(map), ~ factor(., levels=unique(.))) %>% 
        group_by(map)
    
    glasser$data <- glasser$data %>% filter(hemi=='left')
    
    # set scale limits at 99th percentile
    m_max <- pmax(
        df %>% .$value %>% quantile(.99, na.rm=T) %>% abs,
        df %>% .$value %>% quantile(.01, na.rm=T) %>% abs
    )
    if (colorscale=='symmetric') {
        m_min <- -m_max
    } else if (colorscale=='zero') {
        m_min <- 0
    } else if (colorscale=='absolute') {
        m_min <- df %>% .$value %>% quantile(.01)
    } else {
        m_min <- colorscale[1]
        m_max <- colorscale[2]
    }

    # set manual axis labels if desired
    if (labels=='none') {
        labels = c(round(m_min,2),round(m_max,2))
    } else if (labels=='centile') {
        labels = c(round(m_min+0.5,2),round(m_max+0.5,2))
    }

    p <- ggplot(df) + 
    geom_brain(
        atlas=glasser,
        mapping=aes(fill=value, geometry=geometry, hemi=hemi, side=side, type=type),
        colour='grey', size=.1,
        show.legend=T
        ) + 
    # facet_wrap(~map, ncol=ncol, dir="v") +
    scale_fill_gradientn(
        colors=colors, 
        limits=c(m_min,m_max), oob=squish, breaks=c(m_min,m_max), 
        labels=labels, 
        name=name
    ) +
    theme_void() + 
    theme(legend.position='bottom',
          legend.title=element_text(vjust=1),
          panel.spacing.x=unit(spacing,'lines'),
          panel.spacing.y=unit(spacing,'lines'),
          strip.text.x=element_text(vjust=1),
        #   strip.clip='off',
          plot.title=element_text(hjust=0.5)) +
ggtitle(title) + xlab("") + ylab("")
    
    if (facet=='h') {
        p + facet_grid(.~map)
    } else if (facet=='v') {
        p + facet_grid(map~.)
    } else if (facet=='w') {
        p + facet_wrap2(~map, ncol=ncol, dir="v",
                strip=strip_vanilla(clip='off')
        )
    }
}


plot_maps_dk <- function(maps, title="", ncol=3, facet='w', spacing=0,
                      colors=rev(brewer.rdbu(100)), colorscale='symmetric', 
                      name='Z-score', labels='none') {
    df <- maps %>%
        rownames_to_column %>%
        rename(label = rowname) %>%
        gather('map', 'value', -label) %>%
        mutate_at(vars(map), ~ factor(., levels=unique(.))) %>% 
        group_by(map)
    
    dk$data <- dk$data %>% filter(hemi=='left')
    
    # set scale limits at 99th percentile
    m_max <- pmax(
        df %>% .$value %>% quantile(.99) %>% abs,
        df %>% .$value %>% quantile(.01) %>% abs
    )
    if (colorscale=='symmetric') {
        m_min <- -m_max
    } else if (colorscale=='zero') {
        m_min <- 0
    } else if (colorscale=='absolute') {
        m_min <- df %>% .$value %>% quantile(.01)
    }

    # set manual axis labels if desired
    if (labels=='none') {
        labels = c(round(m_min,2),round(m_max,2))
    } else if (labels=='centile') {
        labels = c(round(m_min+0.5,2),round(m_max+0.5,2))
    }

    p <- ggplot(df) + 
    geom_brain(
        atlas=dk,
        mapping=aes(fill=value, geometry=geometry, hemi=hemi, side=side, type=type),
        colour='grey', size=.1,
        show.legend=T
        ) + 
    scale_fill_gradientn(
        colors=colors, 
        limits=c(m_min,m_max), oob=squish, breaks=c(m_min,m_max), 
        labels=labels, 
        name=name
    ) +
    guides(fill=guide_colorbar(title.vjust=1)) +
    theme_void() + 
    theme(legend.position='bottom',
          # text=element_text(color='grey30'),
          panel.spacing.x=unit(spacing,'lines'),
          panel.spacing.y=unit(spacing,'lines'),
          strip.text.x=element_text(vjust=1, size=20),
          strip.text.y.left=element_text(vjust=.5, size=10, angle=0),
          plot.title=element_text(hjust=0.5)) +
ggtitle(title) + xlab("") + ylab("")
    
    if (facet=='h') {
        p + facet_grid(.~map)
    } else if (facet=='v') {
        p + facet_grid(map~., switch='y')
    } else if (facet=='w') {
        p + facet_wrap(~map, ncol=ncol, dir="v")
    }
}



plot_map_corrs <- function(null_p, size=6) {
    lim <- max(abs(null_p$r))
    
    null_p %>%
    mutate(map = factor(map, ordered=T, levels=rev(unique(.$map)))) %>%
    mutate(p_level=ifelse(q<0.001, '***',
                   ifelse(q<0.01, '**',
                   ifelse(q<0.05, '*','')))) %>%
    ggplot(aes(x=G, y=map)) +
    geom_raster(aes(fill=r)) + 
    geom_text(aes(label=p_level), vjust=.5, hjust=.5, size=size, color='white') +
    scale_fill_gradientn(colors=rev(brewer.rdbu(100)), name='R', 
                        limits=c(-lim,lim), breaks=c(-floor(lim*10)/10,floor(lim*10)/10)) +
    scale_x_discrete(position='top') +
    guides(fill=guide_colorbar(barwidth=5)) +
    xlab('') + ylab('') +
    coord_fixed() +
    theme_minimal() +
    theme(panel.grid=element_blank(),
          axis.text = element_text(size=22, color='grey7', family='Calibri'),
          legend.title=element_text(vjust=1),
          legend.position='bottom'
         )
}

plot_pca_weights <- function(pca_weights, size=6) {
    lim <- max(abs(null_p$r))
    
    pca_weights %>%
    gather(PC, weight, -map_name) %>%
    mutate(map_name = factor(map_name, ordered=T, levels=rev(unique(.$map_name)))) %>%
    ggplot(aes(x=PC, y=map_name)) +
    geom_raster(aes(fill=weight)) + 
    scale_fill_gradientn(colors=rev(brewer.puor(100)[10:90]), name='Weight', 
                         limits=c(-lim,lim), breaks=round(c(-lim,lim),1)) +
    scale_x_discrete(position='top') +
    guides(fill=guide_colorbar(barwidth=5)) +
    xlab('') + ylab('') +
    coord_fixed() +
    theme_minimal() +
    theme(panel.grid=element_blank(),
          axis.text = element_text(size=20, color='black'),
          legend.title=element_text(vjust=1),
          legend.position='bottom'
         )
}


plot_maps_scatter <- function(maps_scatter, maps_scatter_corrs, facet='v', x=0, y=0,size=8,
                             xlab='', ylab='') {
    df <- maps_scatter
    # gather(G, G_score, -label, -map, -map_score)

    corrs <- maps_scatter_corrs %>%
    mutate(map = factor(map, ordered=T, levels=unique(.$map))) %>%
    mutate(sig_label=ifelse(q<0.001, '***',
                   ifelse(q<0.01, '**',
                   ifelse(q<0.05, '*','')))) %>%
    mutate(r_label=paste('R =', round(r,2), sig_label))
    # mutate(x=-1, y=2)

    # if (!is.null(which) ) {
    #     df <- df %>% filter(map %in% which)
    #     corrs <- corrs %>% filter(map %in% which)
    # }
    
    p <- df %>%
    mutate(map = factor(map, ordered=T, levels=unique(.$map))) %>%
    ggplot(aes(x=G_score, y=map_score)) +
    geom_point(alpha=.3) +
    geom_smooth(method='lm', color=brewer.puor(5)[5], size=2) +
    geom_text(data=corrs, aes(label=r_label, x=x, y=y), size=size, hjust=0.5, vjust=0.5,
             face='bold') +
    # geom_text(data=corrs, aes(label=p_label, x=x, y=y), size=6, hjust=0, vjust=1) +
    scale_x_continuous(breaks=0, position='top') + scale_y_continuous(breaks=0) +
    # xlim(c(-3,3)) + ylim(c(-3,3)) +
    xlab(xlab) + ylab(ylab) +
    coord_cartesian(clip='off') +
    theme_minimal() +
    theme(
    # strip.placement='outside',
          strip.text.x=element_text(),
          strip.text.y.left=element_text(angle=0),
          panel.grid.minor=element_blank(),
          axis.text = element_blank(),
          axis.title.y = element_text(angle=0, vjust=.5),
          plot.title = element_text(hjust=0.5, vjust=1),
          aspect.ratio=1
         )
    
    if (facet=='v') {
        p + facet_grid(map~G, switch='both')
    } else if (facet=='h') {
        p + facet_grid(G~map, switch='y')
    } else {
        p
    }
}




plot_maps_pca_scatter <- function(maps_pca_scatter, maps_pca_corrs, size=7, ncol=2, x=-1, y=2) {
    # corr <- both_scores_scatter %>% 
    #     select(-label) %>% 
    #     group_by(G) %>% 
    #     summarize(cor(`MRI PCA`, `AHBA DM`, use='p')) %>%
    #     rename_with(function(x) c('G','r'), everything()) %>%
    #     mutate(r = round(r, 2)) %>%
    #     mutate(x=-1, y=2)
    
    
    corrs <- maps_pca_corrs %>%
    mutate(map = factor(map, ordered=T, levels=unique(.$map))) %>%
    mutate(sig_label=ifelse(q<0.001, '***',
                   ifelse(q<0.01, '**',
                   ifelse(q<0.05, '*','')))) %>%
    mutate(r_label=paste('R =', round(r,2), sig_label))
    
    maps_pca_scatter %>%
    ggplot(aes(x=`AHBA DM`, y=`MRI PCA`)) +
    facet_wrap(~G, ncol=ncol, scales='free') +
    geom_point(alpha=.3) +
    geom_smooth(method='lm', color=brewer.set1(5)[4]) +
    geom_text(data=corrs, aes(label=r_label, x=x, y=y), size=size, hjust=0.5, vjust=0.5) +
    scale_x_continuous(breaks=0) + scale_y_continuous(breaks=0) +
          #strip.text.x=element_text(size=20),
    xlab('') + ylab('') +
    theme_minimal() +
    theme(strip.placement='outside',
          #strip.text.y.left=element_text(angle=0, size=20),
          axis.title.y = element_text(angle=0, vjust=0.5),
          panel.grid.minor=element_blank(),
          axis.text = element_blank(),
          aspect.ratio=1,
          #text=element_text(size=20)
         )
}



plot_ahba_mri_brains <- function(scores_df, colors=rev(brewer.rdbu(100)), ncol=2) {
    df <- scores_df %>% 
        mutate(region = recode(label,'7Pl'='7PL')) %>% select(-label) %>%
        gather('component', 'score', -region) %>%
        mutate(component = factor(component, ordered=T, levels=unique(.$component))) %>%
        group_by(component)

    m <- pmax(
        df %>% .$score %>% quantile(.99) %>% abs,
        df %>% .$score %>% quantile(.01) %>% abs
    )
    
    p <- ggplot(df) + 
    geom_brain(
        atlas=glasser,
        hemi='left',
        mapping=aes(fill=score, geometry=geometry, hemi=hemi, side=side, type=type),
        colour='grey', size=.1,
        show.legend=T
        ) + 
    theme_void() + 
    theme(legend.position='none',
          strip.text.x=element_text(vjust=1, size=20),
          strip.text.y.left = element_text(angle = 0, size=20),
          # panel.spacing.x = unit(spacing_x, 'lines'),
          plot.title=element_text(hjust=0.5, size=20)) +
    scale_fill_gradientn(colors=colors, 
                         limits=c(-m,m), oob=squish, breaks=c(-m,0,m), 
                         labels=c(round(-m,2),0,round(m,2)), name=''
                        ) +
    coord_sf(clip='off') +
    xlab("") + ylab("")
    
    p + facet_wrap(~component, ncol=ncol)
    
    # if (facet=='h') {
    #     p + facet_grid(component~version, switch=switch)
    # } else {
    #     p + facet_grid(version~component, switch=switch)   
    # }
}


                    
                    
##########################################
                    

# plot_map_scatters <- function(corrs, null_corrs, null_p) {
#     corrs <- corrs %>% rownames_to_column('map') %>% gather(G, corr, -map)
#     # null_corrs <- null_corrs %>% gather(G, corr, -map)
#     null_p <- null_p %>% rownames_to_column('map') %>% gather(G, p, -map) %>% mutate(sig = p<.05)

#     null_corrs %>% 
#     left_join(null_p, by = c('map', 'G')) %>% 
#     ggplot() + 
#     facet_grid(factor(map, levels=unique(null_p$map))~G, switch='y') +
#     # geom_histogram(aes(corr, alpha=sig)) +
#     geom_density(aes(corr, alpha=sig), color=mycolors[5], fill=mycolors[4]) +
#     scale_alpha_manual(values=c(.2,1), guide='none') +
#     ylim(c(-.6, NA)) +
#     ylab('') + 
#     geom_point(data=corrs, aes(x=corr, y=0), color=mycolors[1], size=5) +
#     # geom_vline(data=corrs, aes(xintercept=corr), color='red') +
#     geom_text(data=corrs, aes(label=paste0('r=',round(corr,2)), x=.4, y=1.4), size=8, hjust=0) +
#     geom_text(data=null_p, aes(label=paste0('p=',round(p,2)), x=.4, y=.8), size=8, hjust=0) +
#     # ggtitle('Distribution of map correlations over spin permutations') +
#     theme_void() + theme(strip.text = element_blank())
# }

# plot_maps_combined <- function(maps, corrs, null_corrs, null_p) {
#     g1 <- plot_null_corrs(corrs, null_corrs, null_p)
#     g2 <- plot_hcp_wide(scores_plot, spacing=0) + guides(fill='none')
#     g3 <- plot_maps(maps, colors=rev(brewer.puor(100)), ncol=1) + guides(fill='none')

#     (
#         ((plot_spacer() | g2) + plot_layout(widths=c(1,2))) / 
#         ((g3 | g1) + plot_layout(widths=c(1,2)))
#     ) + plot_layout(heights=c(1.2,12))       
# }



plot_corr_versions <- function(corr_versions, size=6, facet='h', nrow=1) {
    lim <- max(abs(corr_versions$r))
    
    p <- corr_versions %>%
    mutate(map = factor(map, ordered=T, levels=rev(unique(.$map)))) %>%
    mutate(version = factor(version, ordered=T, levels=unique(.$version))) %>%
    mutate(p_level=ifelse(p<0.001, '***',
                   ifelse(p<0.01, '**',
                   ifelse(p<0.05, '*','')))) %>%
    mutate(p_sig = paste(round(p,3), p_level)) %>%
    mutate(q_sig = ifelse(q<0.05, T, F)) %>%
    ggplot(aes(x=G, y=map)) +
    # facet_grid(.~version) +
    geom_raster(aes(fill=r)) + 
    geom_tile(aes(color=q_sig), fill='transparent', size=2) + 
    geom_text(aes(label=round(r,2)), vjust=-.5, size=size) +
    geom_text(aes(label=p_sig), vjust=1, size=size) +
    scale_fill_gradientn(colors=rev(brewer.rdbu(100)[10:90]), name='r', limits=c(-lim,lim), 
            breaks=c(-floor(lim*100)/100, floor(lim*100)/100)) +
    scale_color_manual(values=c('transparent', 'green'), labels=c('', 'FDR sig'), name='') +
    guides(fill=guide_colorbar(barwidth=10)) +
    xlab('') + ylab('') +
    coord_fixed() +
    theme_minimal() +
    theme(panel.grid=element_blank(),
          text = element_text(size=20),
          legend.title=element_text(vjust=1),
          legend.position='bottom'
         )
    
    if (facet=='h') {
        p + facet_grid(.~version)
    } else if (facet=='w') {
        p + facet_wrap(~version, nrow=nrow)
    } else {
        p
    }
}


plot_corr_versions_2 <- function(corr_versions, size=6) {
    lim <- max(abs(corr_versions$true_mean))
    
    corr_versions %>%
    mutate(map = factor(map, ordered=T, levels=rev(unique(.$map)))) %>%
    mutate(version = factor(version, ordered=T, levels=unique(.$version))) %>%
    mutate(sig_label = case_when(
        p < .001 ~ '***',p < .01 ~ '**',p < .05 ~ '*', TRUE ~ ''
        )) %>%
    # mutate(p_sig = paste(round(q,3), q_level)) %>%
    ggplot(aes(x=G, y=map)) +
    facet_grid(.~version) +
    geom_raster(aes(fill=true_mean)) + 
    
    geom_tile(aes(fill=true_mean, color=sig), size=1) +
    geom_text(aes(label=paste(round(true_mean, digits = 2), '\n', round(p, digits=3), sig_label)), size=8) +
    scale_color_manual(values=c('transparent','green'), name='FDR sig') +
    scale_fill_gradientn(colours=rev(brewer.rdbu(100)[20:80]), guide='colourbar') +
    guides(fill=guide_colourbar(title='Pearson r', barwidth=10)) +
    
#     geom_tile(aes(color=q_sig), fill='transparent', size=2) + 
#     geom_text(aes(label=round(true_mean,2)), vjust=-.5, size=size) +
#     geom_text(aes(label=p), vjust=1, size=size) +
#     scale_fill_gradientn(colors=rev(brewer.rdbu(100)[20:80]), name='Corr', limits=c(-lim,lim)) +
#     scale_color_manual(values=c('transparent', 'green'), labels=c('', 'FDR sig'), name='') +
#     guides(fill=guide_colorbar(barwidth=10)) +
    xlab('') + ylab('') +
    coord_fixed() +
    theme_minimal() +
    theme(panel.grid=element_blank(),
          text = element_text(size=20),
          legend.title=element_text(vjust=1),
          legend.position='bottom'
         )
}
                    
                    
                    
plot_null_corrs <- function(corrs, null_corrs, null_p) {

    null_corrs <- null_corrs %>% gather(G, corr, -map)
    corrs <- corrs %>% rownames_to_column('map') %>% gather(G, corr, -map) %>% left_join(null_p, by=c('G', 'map'))

    null_corrs %>% 
    left_join(null_p, by = c('map', 'G')) %>%
    ggplot() + 
    facet_grid(factor(map, levels=unique(null_p$map), ordered=T)~G, switch='y') +
    geom_density(aes(corr), color=mycolors[5], fill=mycolors[4]) +
    # facet_grid(.~G) +
    # stat_density_ridges(geom = "density_ridges_gradient", calc_ecdf = TRUE, scale=.8, size=.2, color='grey', alpha=.5,
    #                    aes(x=corr, y=factor(map, levels=unique(null_p$map), ordered=T), fill=stat(ecdf))) +
    # scale_fill_gradientn(colors=rev(brewer.rdbu(100)[15:85]), name = 'Quantile') +
    # geom_point(data=corrs, aes(x=corr, y=map, fill=p), shape=21, size=15) +
    
    ylim(c(-.6, NA)) +
    ylab('') + 
    geom_point(data=corrs, aes(x=corr, y=0, alpha=sig), color=brewer.puor(10)[9], size=10) +
    scale_alpha_manual(values=c(.4,1), guide='none') +

    geom_text(data=corrs, aes(label=paste0('r=',round(corr,2)), x=.4, y=1.4), size=8, hjust=0) +
    geom_text(data=null_p, aes(label=paste0('p=',round(p,2)), x=.4, y=.8), size=8, hjust=0) +
    theme_void() + theme(strip.text = element_blank())
}



# library(ggrepel)
# plot_map_weights <- function(map_weights) {
#     lim <- max(map_weights)+.1
#     # lim2 <- lim+0.1

#     df <- map_weights %>%
#     rownames_to_column('map')
    
#     df %>%
#     ggplot(aes(G1, G2)) +
#     geom_point() +
#     geom_text_repel(aes(label=map), size=5, data=df) + 
#     geom_vline(xintercept=0) + geom_hline(yintercept=0) +
#     lims(x=c(-lim,lim), y=c(-lim,lim)) +
#     coord_cartesian(clip='off') +
#     labs(x='',y='') +
#     annotate(geom='text', x=c(-lim, 0), y=c(0,-lim), label=c('PC1', 'PC2'), vjust=.5, hjust=.5, size=5, fontface='bold') +
#     theme_minimal() +
#     theme(aspect.ratio=1,
#           panel.grid=element_blank(),
#           axis.text=element_blank(),
#           text=element_text(size=20)
#          )
# }