theme_set(theme_classic())
theme_update(
    text = element_text(size=7, family = 'Calibri', color = 'grey7'),
    axis.text = element_text(size=7, family = 'Calibri', color = 'grey7'),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    strip.text = element_text(size=7, family = 'Calibri', color = 'grey7'),
    strip.background = element_blank(),
    strip.placement='outside',
    legend.key.size = unit(1, "mm"),
    plot.tag.position = c(0,1),
    plot.tag = element_text(size=10, face='bold', family='Calibri', hjust=0)
)
theme_colorbar <- guide_colorbar(barwidth=2, barheight=.5, ticks=FALSE)


plot_single_cell_posneg <- function(sc_projected_posneg_plot) {
    colors <- c(
        brewer.rdbu(21)[c(2,20)],
        brewer.piyg(11)[c(2,10)],
        brewer.brbg(11)[c(2,10)]
    ) %>% brightness(scalefac(1.5)) %>% c('grey30')

    p1 <- sc_projected_posneg_plot %>% 
    mutate(cell_type = factor(cell_type, ordered=T, levels=unique(.$cell_type))) %>% 
    ggplot(aes(positive, negative)) + 
    facet_wrap(~C, scales='free') + 
    geom_point(alpha=.1, size=.002, aes(color=cell_type)) +
    xlab("sc-RNAseq expression of C1/2/3 positive genes") +
    ylab("sc-RNAseq expression of\nC1/2/3 negative genes") +
    scale_color_manual(values=colors, name=NULL) +
    guides(colour = guide_legend(byrow=T, override.aes = list(size=1, alpha=.8))) +
    theme(
        strip.background = element_blank(),
        panel.grid = element_blank(),
        panel.spacing.x = unit(5,'mm'),
        legend.spacing.y = unit(1,'mm'),
        legend.text = element_text(size=7),
        axis.line = element_line(size=.2, color='darkgrey'),
        axis.text = element_blank(),
        axis.title.y = element_text(margin=margin(t=0,b=0,l=0,r=-2, unit='mm')),
        axis.ticks = element_blank()
    )

    single_inset <- function(C) {
        data_subset <- sc_projected_posneg_plot %>% 
        filter(subclass_label=='VIP', layer=='L2') %>% 
        filter(C==!!enquo(C))

        label <- paste0("r = ", round(cor(data_subset$positive, data_subset$negative),2))

        data_subset %>% 
        ggplot(aes(positive, negative)) +
        geom_smooth(method='lm', se=F, color='darkgrey', size=.3) + 
        geom_point(alpha=.1, size=.01, color=colors[2]) +
        annotate(geom='text', label=label, x=Inf, y=Inf, hjust=1.1, vjust=1.7, size=2,
                 family='Calibri', color='grey7') +
        theme_minimal() + 
        theme(
            aspect.ratio = 1,
            panel.border = element_rect(fill=NA, size=.2, color='darkgrey'),
            panel.grid = element_blank(),
            axis.text = element_blank(),
            axis.title = element_blank(),
            plot.title = element_text(hjust=.5, size=6, family='Calibri', color='gray7')
        ) +
        ggtitle("Layer 2 VIP\ninterneurons only")
    }

    inset_plots <- c('C1','C2','C3') %>% map(single_inset)
    inset_df <- tibble(
        x=c(1,1,1), 
        y=c(1,1,1),
        plot = inset_plots, C=c('C1','C2','C3'))

    p1 + geom_plot_npc(data=inset_df, 
        aes(npcx=x, npcy=y, label=plot), 
        vp.width=0.5, vp.height=0.5)
}


plot_brain_maps <- function(scores_df, 
                        facet = 'component~version', 
                        switch = 'y', spacing_x = 0,
                        colors = rev(brewer.rdbu(100))
                        ) {
    
    # Clean up labels and fix ordering
    df <- scores_df %>% 
        mutate_at(vars(version), ~ factor(., levels=unique(.))) %>% 
        mutate(region = recode(label,'7Pl'='7PL')) %>% select(-label) %>%
        filter(region != '10pp') %>% # region 10pp is not visualised in ggseg
        gather('component', 'score', -version, -region) %>%
        group_by(component, version) %>%
        mutate(component = factor(component, ordered=T, 
                                levels=unique(.$component), 
                                labels=unique(.$component)))

    # Set colorscale limits
    q99 <- pmax(
        df %>% .$score %>% quantile(.99, na.rm=T) %>% abs,
        df %>% .$score %>% quantile(.01, na.rm=T) %>% abs
    )
    breaks <- floor(q99)
    
    # Plot brains with ggseg
    p <- df %>% 
    ggseg(
        atlas = glasser,
        hemi = 'left',
        mapping = aes(fill = score),
        colour = 'grey', size = .05,
        show.legend = T
    ) + 
    facet_grid(facet, switch = switch) +
    scale_fill_gradientn(
        colors = colors, 
        limits = c(-q99,q99), oob = squish, 
        breaks = c(-breaks,breaks), 
        labels = paste0(c(-breaks,breaks),'σ'), 
        guide = theme_colorbar,
        name = NULL
    ) +
    theme_void() + 
    theme( 
        text = element_text(size=7, family='Calibri', color='grey7'),
        legend.text = element_text(size=6, family='Calibri', color='grey7'),
        strip.text = element_text(size=7, family='Calibri', color='grey7'),
        strip.text.x = element_text(vjust = 1, face='bold', margin=margin(t=0,b=2,l=0,r=0, unit='mm')),
        strip.text.y.left = element_text(angle = 0),
        legend.position = c(.5,-.1),
        legend.direction = 'horizontal',
        legend.title = element_text(vjust=1),
        plot.tag.position = c(0,1),
        plot.tag = element_text(size=10, face='bold', family='Calibri', hjust=0)
    )
}



plot_ahba_bs_scatter <- function(both_scores, corrs) {
    # lim <- 2.8
    
    g <- ggplot(both_scores, aes(AHBA, Brainspan)) + 
    facet_rep_grid(C~., repeat.tick.labels=T) +
    geom_point(aes(fill=AHBA), size=1.2, shape=21, stroke=.3, color='darkgrey') +
    geom_smooth(method='lm', linetype=1, se=FALSE, color='darkgrey', size=.3) +
    scale_fill_gradientn(colors=rev(brewer.rdbu(100)), guide='none') +
    geom_text(data=data.frame(C=names(corrs), 
                              label=paste("r =", round(corrs, 2)))
              , aes(x=-.5,y=1.2, label=label), size=2.6
    ) +
    coord_fixed() +
    xlab('AHBA Score') + ylab('BrainSpan Score') +
    scale_y_continuous(breaks=0) +
    scale_x_continuous(breaks=0) +
    theme(
        panel.grid=element_blank(),
        axis.text = element_blank(),
        strip.text = element_blank(),
        aspect.ratio=1
    )
}



plot_bs_scores_corr <- function(bs_scores_corr) {
    colors <- c(
        brewer.puor(10)[8], brewer.brbg(10)[8], brewer.rdbu(10)[2]
    )

    g <- bs_scores_corr %>% 
    mutate_at(vars(age), ~ factor(., levels=unique(.))) %>%
    ggplot() + 
    geom_line(aes(x=age, y=corr, color=C, group=C), size=.3, alpha=1) + 
    geom_point(aes(x=age, y=corr, color=C), size=1.2) + 
    xlab("") + ylab("AHBA-BrainSpan correlation") +
    scale_color_manual(values=colors, guide=guide_legend(byrow=T)) +
    scale_y_continuous(limits=c(0,1), breaks=seq(0,1,.2)) +
    theme(
        axis.line = element_blank(),
        axis.text.x = element_text(margin = margin(-1,0,0,0, 'mm')),
        axis.text.y = element_text(margin = margin(0,-4,0,0, 'mm')),
        panel.grid = element_blank(),
        legend.position = 'right',
        legend.title = element_blank(),
        legend.spacing.y = unit(1, 'mm'),
        legend.text = element_text(size=7)
    )
}




plot_quantile_curves <- function(quantile_curves, facet='~C', continuous=FALSE, ncol=3, which='pred') {
    n_quantiles <- quantile_curves$C_quantile %>% unique %>% length

    p <- quantile_curves %>% 
    mutate(curve = get(which)) %>% 
    arrange(desc(C_quantile)) %>% 
    mutate(C_quantile = factor(C_quantile, ordered=T, levels=unique(.$C_quantile))) %>% 
    ggplot(aes(x=age_log10, y=10**curve)) +
    facet_wrap(as.formula(facet), scales='fixed', ncol=ncol) +
    geom_line(aes(color=C_quantile, group=C_quantile, alpha=C_quantile), size=.3) +
    scale_x_continuous(
        breaks=log10(c(-0.5,0,1,5,14,40)*365+40*7),
        labels=function(x) round((10**x - 40*7)/365,2)) +
    scale_color_manual(values=brewer.rdbu(10), name='Decile', labels=seq(10,1)) +
    scale_alpha_manual(values=rep(1, n_quantiles), name='') +
    guides(color=guide_legend(override.aes = list(size=2)), alpha='none') +
    ylab('log10 RPKM') +
    xlab('Age') +
    coord_cartesian(clip='off') +
    theme(
        axis.text.y=element_blank(),
        axis.title.y = element_text(margin=margin(t=0,b=0,l=0,r=-5, unit='mm')),
        strip.text.x = element_text(margin=margin(0,0,-2,0,'mm')),
        legend.text = element_text(size=7),
        plot.title.position='plot'
    )
}
