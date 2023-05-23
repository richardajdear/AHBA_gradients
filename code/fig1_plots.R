

plot_generalisability <- function(triplets_df)
    triplets_df %>% 
        mutate(version=factor(version, ordered=T, levels=unique(.$version))) %>% 
        ggplot(aes(x=component, y=corr_abs, fill=version)) + 
        geom_point(aes(color=version), position=position_dodge(width = 0.5), size=3) +
        stat_summary(geom='line', size=1, fun=median, position=position_dodge(width = 0.5), 
                    aes(color=version, group=version)) +
        scale_color_manual(values=brewer.rdbu(10)[c(9,2)]) +
        scale_y_continuous(limits=c(0,1), breaks=seq(0,1,.2)) +
        geom_hline(yintercept=.6, linetype='dashed') +
        ylab('Generalisability') + xlab('Component') +
        theme_minimal() + 
        theme(
        # plot.tag.position = c(0,.95),
        # axis.text=element_text(size=22, color='grey7', family='Calibri'), 
        # axis.text.y=element_text(margin=margin(0,-20,0,0,'pt')),
        # axis.title.y=element_text(margin=margin(0,0,0,0,'pt')),
        text = element_text(size=20, family='Calibri', color='grey7'),
        axis.text = element_text(size=18, family='Calibri', color='grey7'),
        legend.text = element_text(size=20, family='Calibri', color='grey7'),
        panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.spacing.y = unit(3, 'pt'),
        legend.margin=margin(0,0,-3.5,0, 'cm'),
        legend.position='top'
        )


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


plot_coexp <- function(coexp_df) {
    # Take 99th percentile as limit for colorscale
    limit <- coexp_df %>% filter(x!=y) %>% .$r %>% abs %>% quantile(.99)

    coexp_df %>%
        ggplot(aes(x, rev(y))) +
        geom_raster(aes(fill=r)) +
        scale_fill_gradientn(colors=rev(brewer.rdbu(200)), name='R', 
                             limits=c(-limit,limit), 
                             breaks=c(-floor(limit*10)/10,floor(limit*10)/10),
                             labels=c(-floor(limit*10)/10,floor(limit*10)/10)
                             ) +
        coord_cartesian(clip='off') +
        theme_void() +
        ylab("G1-regressed\nregion-region\nco-expression") +
        guides(fill=guide_colorbar(barwidth=5, title.vjust=1)) +
        # ggtitle("G1-regressed\nregion-region\nco-expression") +
        theme(
            aspect.ratio=1,  
            text = element_text(size=20, family='Calibri', color='grey7'),
            axis.title.y = element_text(size=20, family='Calibri', color='grey7', angle=0, vjust=.8),
            # plot.title.position='panel',
            # plot.title = element_text(size=20, family='Calibri', color='grey7',
            #                           hjust=.5,
            #                           margin=margin(t=2,b=0,l=-5,r=0, unit='cm')),
            legend.text = element_text(size=18, family='Calibri', color='grey7'),
            # legend.position='right'
            legend.direction = 'horizontal',
            legend.position = c(-.3, .3)
        )
}


plot_go_enrichments <- function(enrichments) {
    df <- go_enrichments %>% 
    mutate(enrichment = ifelse(direction=='top', -enrichment, enrichment)) %>% 
    mutate(rank_directed = rank * ifelse(direction=='top', -1, 1)) %>% 
    arrange(-rank_directed)

    n_label <-  3
    label_size <- 7

    df_title <- go_enrichments %>% 
        group_by(G) %>% 
        summarize(n=n()) %>% 
        mutate(label = paste0(G, ": ", n, " enrichments"))

    df_label <- df %>% filter(FDR_rank <= n_label | 
            description %in% c(
                'Learning or memory'
                )
            )
    df_label_extra <- df %>% filter(
            description %in% c(
                # G3 top
                # 'Learning or memory',
                'Vocal learning',
                'Regulation of short-term plasticity',
                # G3 bottom
                'Adaptive immune response',
                'Immune response to tumor cell',
                # 'Regulation of vasculature development',
                'Regulation of macrophage migration',
                # G2 top
                'Regulation of chromatin silencing',
                # 'Regulation of mRNA splicing',
                # G2 bottom
                'Response to jasmonic acid',
                'Intestinal cholesterol absorption',
                'Regulation of mitochondrial fusion',
                # G1 top 
                '7-methylguanosine mRNA capping'
                # G1 bottom
                # 'Peptidyl-lysine modification'
            )
        )

    df %>% 
    ggplot(aes(x=neglogFDR, y=enrichment)) + 
    facet_grid(.~G) +
    geom_point(aes(fill=enrichment, size=n_genes), shape=21, color='grey7') + 
    # Top-n FDR labels
    geom_text_repel(aes(label=description), data=df_label %>% filter(direction=='top'), 
                    size=label_size, family='Calibri', segment.linetype=1, segment.color='grey',
                    force=10, force_pull=0, xlim=c(2,5), ylim=c(-2,-8), max.time=1, nudge_y=-2) + 

    geom_text_repel(aes(label=description), data=df_label %>% filter(direction=='bottom', G!='G1'), 
                    size=label_size, family='Calibri', segment.linetype=1, segment.color='grey',
                    force=10, force_pull=0, xlim=c(2,5), ylim=c(0,8), max.time=1, nudge_y=1.3) + 

    geom_text_repel(aes(label=description), data=df_label %>% filter(direction=='bottom', G=='G1'), 
                    size=label_size, family='Calibri', segment.linetype=1, segment.color='grey',
                    force=10, force_pull=0, xlim=c(2,5), ylim=c(0.5,8), max.time=1, nudge_y=-3) +
    # Top-n FDR labels
    geom_text_repel(aes(label=description), data=df_label_extra %>% filter(direction=='top', G!='G3'), 
                    size=label_size, family='Calibri', segment.linetype=1, segment.color='grey',
                    force=10, force_pull=0, xlim=c(2,5), ylim=c(-2,-9), max.time=1, nudge_y=-.5) +
    geom_text_repel(aes(label=description), data=df_label_extra %>% filter(direction=='bottom', G!='G3'), 
                    size=label_size, family='Calibri', segment.linetype=1, segment.color='grey',
                    force=10, force_pull=0, xlim=c(2,5), ylim=c(2,9), max.time=1, nudge_y=0) +
    geom_text_repel(aes(label=description), data=df_label_extra %>% filter(direction=='top', G=='G3'), 
                    size=label_size, family='Calibri', segment.linetype=1, segment.color='grey',
                    force=10, force_pull=0, xlim=c(1.9,5), ylim=c(-2,-9), max.time=1, nudge_y=0.1, nudge_x=1) +                
    geom_text_repel(aes(label=description), data=df_label_extra %>% filter(direction=='bottom', G=='G3'), 
                    size=label_size, family='Calibri', segment.linetype=1, segment.color='grey',
                    force=10, force_pull=0, xlim=c(2,5), ylim=c(2,9), max.time=1, nudge_y=1.5) +                
    # Titles in middle of plot
    geom_text(aes(label=label), data=df_title, x=-log10(0.001), y=0, size=7.5, 
                family='Calibri', color='grey7', fontface='bold') +
    scale_size_continuous(range=c(2,8)) + 
    scale_fill_gradientn(colors=rev(brewer.rdbu(100)), breaks=c(-7,7)) +
    scale_y_continuous(breaks=c(-4,0,4), name='Enrichment score') +
    scale_x_continuous(breaks = c(0,1.3,2,3,4), labels = c(0, 0.05, 0.01, 0.001, 0.0001), name='Enrichment FDR') +
    coord_cartesian(clip='off') +
    guides(
        fill=guide_colorbar(title='Enrichment score', title.vjust=1),
        size=guide_legend(title='# genes')
    ) +
    theme_minimal() + 
    theme(
        text = element_text(size=20, family='Calibri', color='grey7'),
        axis.text = element_text(size=18, family='Calibri', color='grey7'),
        axis.text.x = element_text(size=18, family='Calibri', color='grey7', margin=margin(t=-2.5,b=0,l=0,r=0,unit='cm')),
        axis.title.x = element_text(size=20, family='Calibri', color='grey7', margin=margin(t=-1.5,b=.2,l=0,r=0,unit='cm')),
        axis.title.y = element_text(size=20, family='Calibri', color='grey7', margin=margin(t=0,b=0,l=0,r=-1,unit='cm')),
        legend.text = element_text(size=18, family='Calibri', color='grey7'),
        legend.title.align=0.5,
        legend.position = 'bottom',
        legend.margin = margin(t=-1,b=.2,l=0,r=0,unit='cm'),
        strip.text=element_blank(),
        panel.grid = element_blank()
    )
}

plot_enrichment_bars_z <- function(stats, xlab='Enrichment z-score', facet='.~G') {
    
    lim <- max(abs(stats$z))*1.1
    
    stats %>% 
    mutate(G = factor(G, ordered = T, levels = unique(.$G))) %>%
    mutate(label = factor(label, ordered = T, levels = unique(.$label))) %>%
    mutate(sig  =  case_when(
        q < .001 ~ '***',q < .01 ~ '**',q < .05 ~ '*', TRUE ~ ''
        )) %>%
    mutate(hjust = ifelse(z < 0, 1, 0)) %>%
    ggplot(aes(x = z, y = label)) + 
    facet_grid(facet) +
    geom_col(aes(fill = z)) +
    geom_text(aes(label = sig, vjust = .7, hjust = hjust), size = 8) +
    scale_y_discrete(limits = rev, name = NULL) +
    scale_x_continuous(limits = c(-lim,lim), 
                       breaks = round(0.5*c(-lim,lim)), name = xlab) +
    scale_fill_gradientn(colors = rev(brewer.rdbu(100)), guide = 'none', limits = c(-lim,lim)) +
    coord_cartesian(clip = 'off') +
    theme_minimal() +
    theme(
        text = element_text(size = 20, family = 'Calibri', color = 'grey7'),
        axis.text.y = element_text(size = 20, family = 'Calibri', color = 'grey7'),
        strip.text.x = element_text(size = 20, family = 'Calibri', color = 'grey7'),
        axis.text.x = element_text(size = 18, family = 'Calibri', color = 'grey7'),
        plot.title = element_text(size = 20, family = 'Calibri', color = 'grey7', face='bold'),
        panel.grid = element_blank()
    )
}