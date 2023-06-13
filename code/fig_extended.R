theme_set(theme_classic())
theme_update(
    text = element_text(size=7, family = 'Calibri', color = 'grey7'),
    axis.text = element_text(size=7, family = 'Calibri', color = 'grey7'),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    legend.key.size = unit(1, "mm"),
    legend.position = 'bottom',
    strip.text = element_text(size=7, family = 'Calibri', color = 'grey7'),
    strip.background = element_blank(),
    strip.placement = 'outside',
    strip.clip = 'off',
    plot.title.position='panel',
    plot.title = element_text(size=7, family = 'Calibri', color = 'grey7', face='bold',
                    margin = margin(0,0,0,0,'mm')),
    plot.tag.position = c(0, 1),
    plot.tag = element_text(size=10, face='bold', family='Calibri', hjust=0)
)
colorbar_h <- guide_colorbar(barwidth=2, barheight=.5, ticks=FALSE, title.vjust=1)
colorbar_v <- guide_colorbar(barwidth=.5, barheight=2, ticks=FALSE)


## Raster plot

plot_triplets_raster <- function(triplets_raster, n_components=3, highlight = 'tick') {
    triplets_raster <- triplets_raster %>% 
    mutate(method=factor(method, ordered=T, levels=c('pca','dme'),
            labels=c('PCA', 'DME'))) %>%
    filter(component <= n_components) %>%
    mutate(component=factor(component, ordered=T, levels=seq(1,5),
            labels=c('C1','C2','C3','C4','C5')))

    # highlight_options = triplets_raster %>% 
    # filter(
    #         ((method=='PCA') & (gene_filter==0.0) & (donors_filter==1)) |
    #         ((method=='PCA') & (gene_filter==0.7) & (donors_filter==3)) |
    #         ((method=='DME') & (gene_filter==0.5) & (donors_filter==3)) |
    #         ((method=='DME') & (gene_filter==0.8) & (donors_filter==1))
    # ) %>% 
    # arrange(component, method, gene_filter, donors_filter) %>%
    # mutate(label=rep(c('1','2','3','4'), n_components))

    g <- triplets_raster %>%
    ggplot(aes(x=gene_filter, y=donors_filter, fill=corr_abs)) + 
    geom_tile() +
    facet_grid(component~method, switch='y') +
    scale_fill_gradientn(
            limits=c(0,1), breaks=seq(0,1,.2), 
            colors=rev(brewer.spectral(100)),
            labels=c('0.0', '0.2', '0.4', '0.6', '0.8', '1.0'),
            name='Generalisability\n(median \ntriplet \ncorrelation)'
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
    guides(fill=guide_colorbar(barwidth=.5, barheight=5, ticks=FALSE)) +
    # theme_minimal() +
    theme(
        # text=element_text(size=22),
        # axis.text = element_text(size=22, color='grey7', family='Calibri'),
        aspect.ratio=1/3.8,
        panel.grid=element_blank(),
        # legend.title=element_text(vjust=3, size=22, color='grey7', family='Calibri'),
        legend.position='right',
        strip.text.y.left=element_text(angle=0),
        # strip.text.x=element_text(angle=0,size=22, color='grey7', family='Calibri')
    )

    if (highlight=='options') {
        g + geom_text(aes(label=label), fontface='bold', data=highlight, color='grey7', family='Calibri', size=8)
        # geom_text(aes(label='★'), data=triplets_raster %>% filter(corr_abs>0.6), color='gray30', size=5) +
        # geom_tile(aes(x=gene_filter, y=donors_filter), fill=NA, color='white', size=1,
        #             data=triplets_raster %>% filter(corr_abs > 0.6)
        # ) +
    } else if (highlight=='tick') {
        g + geom_text(aes(label='✔'), data=triplets_raster %>% filter(corr_abs>0.6), color='grey7', size=2) +
            scale_fill_gradientn(
            colors=rev(brewer.spectral(100)),
            limits=c(0,1), breaks=seq(0,1,.2), 
            labels=c('0.0', '0.2', '0.4', '0.6 ✔', '0.8', '1.0'),
            name='Generalisability\n(median \ntriplet \ncorrelation)'
        ) 
    }
}



plot_brains <- function(maps, atlas='hcp',
                      title="", ncol=1, facet='w', spacing=0,
                      colors=rev(brewer.rdbu(100)), colorscale=c(-3,3),
                      name='', labels=c('-3σ','+3σ')) {
    
    if ("label" %in% colnames(maps)) {
        maps <- maps %>% remove_rownames %>% column_to_rownames('label')
    }

    df <- maps %>%
        rownames_to_column %>%
        mutate(region = recode(rowname,'7Pl'='7PL')) %>%
        select(-rowname) %>%
        gather('map', 'value', -region) %>%
        mutate_at(vars(map), ~ factor(., levels=unique(.))) %>%
        group_by(map)
    
    if (atlas %in% c('dk','dx','s400')) {
        df <- df %>% 
            rename(label = region)
    }

    # set scale limits at 99th percentile
    m_max <- pmax(
        df %>% .$value %>% quantile(.99, na.rm=T)
        # df %>% .$value %>% quantile(.01) %>% abs
    )

    if (length(colorscale > 1)) {
        m_min <- colorscale[1]
        m_max <- colorscale[2]
    } else if (colorscale=='symmetric') {
        m_min <- -m_max
    } else if (colorscale=='zero') {
        m_min <- 0
    } else if (colorscale=='absolute') {
        m_min <- df %>% .$value %>% quantile(.01)
    } else {
        print("Invalid colorscale")
    }

    # set manual axis labels if desired
    if (length(labels)>1) {
        labels = labels
    } else if (labels=='none') {
        labels = c(round(m_min,2),round(m_max,2))
    } else if (labels=='centile') {
        labels = c(round(m_min+0.5,2),round(m_max+0.5,2))
    }

    if (atlas=='dk') {
        p <- df %>% 
            ggseg(
                atlas=dk,
                hemi='left',
                mapping=aes(fill=value),
                colour='grey', size=.1,
                show.legend=T
                )
    } else if (atlas=='hcp') {
        p <- df %>% 
            ggseg(
                atlas=glasser,
                hemi='left',
                mapping=aes(fill=value),
                colour='grey', size=.1,
                show.legend=T
                )
    } else if (atlas=='dx') {
        p <- df %>% 
            ggseg(
                atlas=desterieux,
                hemi='left',
                mapping=aes(fill=value),
                colour='grey', size=.1,
                show.legend=T
                )
    } else if (atlas=='s400') {
        p <- df %>% 
            ggseg(
                atlas=schaefer17_400,
                hemi='left',
                mapping=aes(fill=value),
                colour='grey', size=.1,
                show.legend=T
                )
    }

    p <- p +
    scale_fill_gradientn(
        colors=colors, 
        limits=c(m_min,m_max), oob=squish, breaks=c(m_min,m_max), 
        labels=labels, 
        name=name
    ) +
    guides(fill = colorbar_v) +
    theme_void() + 
    theme(
        text = element_text(size=7, family='Calibri', color='grey7'),
        # strip.text = element_text(size=7, family='Calibri', color='grey7'),
        strip.text = element_blank(),
        strip.clip = 'off',
        legend.position = 'right',
        legend.text = element_text(size=6, family='Calibri', color='grey7'),
        plot.title = element_text(size=7, family = 'Calibri', color = 'grey7', face='bold',
                                  hjust=.5, margin = margin(0,0,1,0,'mm')),
        plot.tag.position = c(0, 0.95),
        plot.tag = element_text(size=10, face='bold', family='Calibri', hjust=0)
    )

    if (facet=='h') {
        p + facet_grid(.~map)
    } else if (facet=='v') {
        p + facet_grid(map~., switch='y')
    } else if (facet=='w') {
        p + facet_wrap(~map, ncol=ncol, dir="h")
    }
}


plot_donor_counts <- function(donor_counts) {

    donor_counts %>%
    mutate(region = recode(label,'7Pl'='7PL')) %>%
    select(-label) %>% 
    ggseg(
        atlas=glasser,
        hemi='left',
        mapping=aes(fill=count),
        colour='grey', size=.1,
        show.legend=T
        ) +
    scale_fill_gradientn(colors=rev(brewer.spectral(7)), name='Donors sampled') +
    guides(fill = colorbar_h) +
    # scale_fill_manual(values=brewer.pubugn(7), name='Donors sampled') +
    # guides(fill = guide_legend(nrow=1)) +
    theme_void() + 
    theme(
        text = element_text(size=7, family='Calibri', color='grey7'),
        strip.text = element_blank(),
        strip.clip = 'off',
        legend.position = 'bottom',
        legend.direction = 'horizontal',
        legend.margin = margin(-8,0,0,0,'mm'),
        # legend.key.size = unit(1,'mm'),
        legend.text = element_text(size=6, family='Calibri', color='grey7'),
        plot.title = element_text(size=7, family = 'Calibri', color = 'grey7', face='bold',
                                  hjust=.5, margin = margin(0,0,1,0,'mm')),
        # plot.tag.position = c(0, 0.95),
        plot.tag = element_text(size=10, face='bold', family='Calibri', hjust=0)
    )
}

plot_donor_counts_histogram <- function(donor_counts) {
    donor_counts %>%
    ggplot() + 
    geom_histogram(aes(x=count), fill=rev(brewer.spectral(7)), alpha=1, binwidth=1, color='white', size=.3) +#, fill=brewer.blues(10)[4]) + 
    stat_count(binwidth=1, geom="text", size=2.2,
        aes(x=count, label=..count..), vjust=-0.2) +
    # geom_vline(xintercept=2.5, linetype=2, size=.3) +
    scale_x_continuous(breaks=seq(0,6,1)) +
    ylim(0,55) +
    ylab('# regions') + xlab('Donors sampled') + 
    theme(
        legend.position=c(.2,.8), panel.grid=element_blank(),
        axis.text.y=element_blank(),
        axis.text.x = element_text(size=6),
        axis.title.y=element_text(margin=margin(0,-2,0,0, 'mm'))
    )    
}

plot_gene_stability_density <- function(df_stability) {
    ggplot(df_stability) + 
    geom_density(aes(x=ds), size=.3, color=brewer.blues(10)[6], fill=brewer.blues(10)[4]) +
    geom_vline(xintercept=0.386, linetype=2, size=.3) +
    annotate(x=.39, y=3.6, geom='text',label='10th', hjust=-0.05, size=2.6) +
    geom_vline(xintercept=0.120, linetype=2, size=.3) +
    annotate(x=.120, y=3.6, geom='text',label='50th', hjust=-0.05, size=2.6) +
    ylab('Density') + xlab('Gene differential stability') + 
    theme(
        axis.text.y = element_blank(),
        axis.text.x = element_text(size=6),
        axis.title.y=element_text(margin=margin(0,-2,0,0, 'mm'))
    )
}

plot_versions_scatter <- function(df, x=-2, y=2.5, size=2.6) {
    r = df %>% group_by(version_x, C_x, C_y) %>% summarize(r=cor(score_x, score_y)) %>% 
    mutate(r_label=paste0('r =', round(r,2))) %>% filter(C_x==C_y)

    df %>%
    filter(C_x==C_y) %>%
    ggplot(aes(x=score_x, y=score_y)) + 
    facet_grid(version_x~C_y, switch='y') + 
    geom_point(aes(fill=score_y), size=1, shape=21, stroke=.1) +
    geom_smooth(method='lm', color='darkgrey', size=.3, se=T) +
    geom_text(data=r, aes(label=r_label, x=x, y=y), size=size, hjust=0.5, vjust=0.5) +
    scale_fill_gradientn(colors=brewer.rdbu(100) %>% rev, name=NULL,
                         limits=c(-3,3), oob=squish, breaks=c(-3,3), labels=c('-3σ','+3σ')) +
    scale_y_continuous(position='right', limits=c(-3,3)) +
    scale_x_continuous(limits=c(-3,3)) +
    # guides(fill = colorbar_v) +
    guides(fill = 'none') +
    xlab('DME, top 50% genes, 3+ donor regions') +
    coord_cartesian(clip='off') +
    theme(
        axis.text = element_blank(), 
        axis.title.y = element_blank(),
        panel.border=element_rect(fill='transparent', color='grey', size=.3),
        legend.position='right',
        aspect.ratio=1
        # panel.margin=margin(l=2,r=2,b=2,t=2,'mm')
    )
}




plot_weight_heatmaps <- function(weight_corrs, spacing=1, ylab='') {
    weight_corrs %>% 
    mutate_at(vars(version), ~ factor(., levels=unique(.))) %>% 
    # mutate_at(vars(y), ~ factor(., levels=rev(unique(.)))) %>% 
    ggplot() +
    facet_grid(comparison~version, switch='both') +
    # facet_rep_grid(.~version, repeat.tick.labels=T) + #, switch='x') +
    geom_tile(aes(x,y, fill=corr)) +
    geom_text(aes(x,y, label=sprintf("%0.2f", round(corr, digits = 2))), size=2.6) +
    scale_fill_gradientn(colors=rev(brewer.rdbu(100)[10:90]), name='r', 
                        limits=c(-1,1), breaks=c(-1,1)
                        ) +
    guides(fill=colorbar_v) +
    # scale_x_discrete(position = "top") +
    theme(
        legend.position = 'right',
        legend.margin = margin(0,0,0,-3,'mm'),
        panel.spacing = unit(spacing, 'mm')
    ) +
    coord_fixed() + 
    xlab('') +
    ylab('')
}