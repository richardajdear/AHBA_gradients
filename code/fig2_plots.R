

plot_single_cell_posneg <- function(sc_axes_posneg_plot) {
    colors <- c(
        brewer.rdbu(21)[c(2,20)],
        brewer.piyg(11)[c(2,10)],
        brewer.brbg(11)[c(2,10)]
    ) %>% 
    brightness(scalefac(1.5)) %>% 
    # saturation(scalefac(1.2)) %>% 
    c('grey30')

    g1 <- sc_axes_posneg_plot %>% 
    mutate(cell_type = factor(cell_type, ordered=T, levels=unique(.$cell_type))) %>% 
    ggplot(aes(positive, negative)) + 
    facet_wrap(~G, scales='free') + 
    geom_point(alpha=.2, size=.2, aes(color=cell_type)) +
    xlab("sc-RNAseq expression of G1/2/3 positive genes") +
    ylab("sc-RNAseq expression of\nG1/2/3 negative genes") +
    scale_color_manual(values=colors, name=NULL) +
    # scale_color_manual(values=cols25(10), name=NULL) +
    # scale_color_manual(values=viridis(7), name=NULL) +
    guides(colour = guide_legend(override.aes = list(size=4, alpha=.8))) +
    theme_minimal() + 
    theme(
        # aspect.ratio=1,
        # panel.border=element_rect(fill=NA, size=1),
        strip.background = element_blank(),
        panel.grid = element_blank(),
        panel.spacing.x = unit(1,'cm'),
        text = element_text(size=20, family='Calibri', color='gray7'),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_text(size=20, family='Calibri', color='gray7'),
        axis.title.y = element_text(margin=margin(t=0,b=0,l=0,r=-1, unit='cm')),
        strip.text = element_text(size=20, family='Calibri', color='gray7')
    )

    single_inset <- function(G) {
        data_subset <- sc_axes_posneg_plot %>% 
        filter(subclass_label=='VIP', layer=='L2') %>% 
        filter(G==!!enquo(G))

        label <- paste0("R ", round(cor(data_subset$positive, data_subset$negative),2))

        data_subset %>% 
        ggplot(aes(positive, negative)) +
        geom_smooth(method='lm', se=F, color='black') + 
        geom_point(alpha=.2, size=.5, color=colors[2]) +
        annotate(geom='text', label=label, x=Inf, y=Inf, hjust=1.1, vjust=1.7, size=6,
                 family='Calibri', color='grey7') +
        theme_minimal() + 
        theme(
            aspect.ratio = 1,
            panel.border = element_rect(fill=NA, size=1, color='grey'),
            panel.grid = element_blank(),
            axis.text = element_blank(),
            axis.title = element_blank(),
            plot.title = element_text(hjust=.5, size=18, family='Calibri', color='gray7')
        ) +
        ggtitle("Layer 2 VIP\ninterneurons only")
    }

    inset_plots <- c('G1','G2','G3') %>% map(single_inset)
    inset_df <- tibble(
        x=c(1,1,1), 
        y=c(1,1,1),
        plot = inset_plots, G=c('G1','G2','G3'))


    g1 + geom_plot_npc(data=inset_df, 
        aes(npcx=x, npcy=y, label=plot), 
        vp.width=0.45, vp.height=0.45)
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
        colour = 'grey', size = .1,
        show.legend = T
    ) + 
    facet_grid(facet, switch = switch) +
    scale_fill_gradientn(
        colors = colors, 
        limits = c(-q99,q99), oob = squish, 
        breaks = c(-breaks,breaks), 
        labels = paste0(c(-breaks,breaks),'σ'), name = NULL
    ) +
    theme_void() + 
    theme( 
        text = element_text(size=20, family='Calibri', color='grey7'),
        legend.text = element_text(size=18, family='Calibri', color='grey7'),
        strip.text = element_text(size=20, family='Calibri', color='grey7'),
        strip.text.x = element_text(vjust = 1, face='bold', margin=margin(t=0,b=.5,l=0,r=0, unit='cm')),
        strip.text.y.left = element_text(angle = 0),
        legend.position = c(.5,-.1),
        legend.direction = 'horizontal',
        legend.title = element_text(vjust=1)
    )
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
    theme(panel.grid=element_blank(),
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



plot_bs_scores_corr <- function(bs_scores_corr, title="", xint='Birth-3 yrs', rotate=F) {
    colors <- c(
        brewer.puor(10)[8], brewer.brbg(10)[8], brewer.rdbu(10)[2]
    )

    g <- bs_scores_corr %>% 
    mutate_at(vars(age), ~ factor(., levels=unique(.))) %>%
    ggplot() + 
    # geom_hline(yintercept=0, color='grey') + 
    # geom_vline(xintercept=xint, color='grey') +
    geom_line(aes(x=age, y=corr, color=G, group=G), size=1, alpha=1) + 
    geom_point(aes(x=age, y=corr, color=G), size=5) + 
    xlab("") + ylab("AHBA-BrainSpan correlation") +
    scale_color_manual(values=colors) +
    scale_y_continuous(limits=c(0,1), breaks=seq(0,1,.2)) +
    # scale_color_manual(values=brewer.rdylbu(4)) +
    ggtitle(title) +
    theme_minimal() + 
    theme(panel.grid = element_blank(),
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




plot_quantile_curves <- function(quantile_curves, facet='~G', continuous=FALSE, ncol=3, which='pred') {
    n_labels <- quantile_curves$label %>% unique %>% length

    p <- quantile_curves %>% 
    mutate(curve = get(which)) %>% 
    arrange(desc(label)) %>% 
    mutate(label = factor(label, ordered=T, levels=unique(.$label))) %>% 
    ggplot(aes(x=age_log10, y=10**curve)) +
    geom_line(aes(color=label, group=label, alpha=label), size=1) +
    scale_x_continuous(
        breaks=log10(c(-0.5,0,1,5,14,40)*365+40*7),
        labels=function(x) round((10**x - 40*7)/365,2)) +
    guides(fill=guide_legend(reverse=T), alpha='none') +
    ylab('log10 RPKM') +
    xlab('Age') +
    theme_minimal() + 
    theme(
        text=element_text(size=20, color='grey7', family='Calibri'),
        axis.text.x=element_text(size=20, color='grey7', family='Calibri'), 
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.y = element_text(size=20, color='grey7', family='Calibri', margin=margin(t=0,b=0,l=0,r=-1, unit='cm')),
        panel.grid=element_blank(),
        strip.text.y.left = element_text(angle=0, size=20, color='grey7', family='Calibri'),
        strip.text.x = element_text(angle=0, size=20, color='grey7', family='Calibri'),
        strip.placement='outside',
        strip.background = element_blank(),
        plot.title.position='plot'
    )

    if (!continuous) {
        p <- p + 
            scale_color_manual(values=viridis(n_labels), name='') +
            scale_alpha_manual(values=rep(1,n_labels), name='')
    } else {
        p <- p + 
            scale_color_gradientn(colors=viridis(100), name='') + 
            scale_alpha_continuous(limits=c(1,1))
    }

    if (facet=='') {
        p
    } else (
        p + facet_wrap(as.formula(facet), scales='fixed', ncol=ncol)
    )
}
