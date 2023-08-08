library(ggseg)
library(ggtext)
library(ggsegGlasser)
library(ggrepel)
# library(ggh4x) # needed for facet_grid2 to not clip strip labels
library(ggpmisc)
library(eulerr)
suppressMessages(library(lemon))
suppressMessages(library(scales))
library(pals)
library(shades)
library(patchwork)
suppressMessages(library(tidyverse))


theme_set(theme_classic())
theme_update(
    text = element_text(size=7, family = 'Calibri', color = 'grey7'),
    axis.text = element_text(size=7, family = 'Calibri', color = 'grey7'),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    # legend.key.size = unit(1, "mm"),
    legend.position = 'bottom',
    legend.text = element_text(size=6, family = 'Calibri', color = 'grey7'),
    strip.text = element_text(size=7, family = 'Calibri', color = 'grey7'),
    strip.background = element_blank(),
    strip.placement = 'outside',
    strip.clip = 'off',
    plot.title.position='plot',
    plot.title = element_text(size=7, family = 'Calibri', color = 'grey7', face='bold',
                    hjust = 0.25, margin = margin(-1,0,0,0,'mm')),
    plot.tag.position = c(0, 1),
    plot.tag = element_text(size=10, face='bold', family='Calibri', hjust=0, color='grey7')
)
theme_colorbar <- guide_colorbar(barwidth=2, barheight=.5, ticks=FALSE)


plot_generalisability <- function(triplets_df)
    triplets_df %>% 
        mutate(version=factor(version, ordered=T, levels=unique(.$version))) %>% 
        ggplot(aes(x=component, y=corr_abs, fill=version)) + 
        geom_point(aes(color=version), position=position_dodge(width = 0.5), size=1, alpha=.5) +
        stat_summary(geom='line', size=.3, fun=median, position=position_dodge(width = 0.5), 
                    aes(color=version, group=version)) +
        scale_color_manual(values=brewer.rdbu(10)[c(9,2)]) +
        scale_y_continuous(limits=c(0,1), breaks=seq(0,1,.2)) +
        scale_x_discrete(labels=c('C1','C2','C3','C4','C5')) +
        geom_hline(size=.1, yintercept=.6, linetype='dashed') +
        ylab(expression(italic('g'))) +
        theme(
            axis.line = element_blank(),
            axis.title.y = element_text(angle=0, vjust=.5),
            axis.title.x = element_blank(),
            legend.text = element_text(size=7),
            legend.title = element_blank(),
            legend.spacing.y = unit(2, 'mm'),
            legend.margin = margin(0,0,-15,0, 'mm'),
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
        strip.text.x = element_text(vjust = 1, face='bold', margin=margin(t=0,b=3,l=0,r=0, unit='mm')),
        strip.text.y.left = element_text(angle = 0),
        legend.position = c(.5,-.1),
        legend.direction = 'horizontal',
        legend.title = element_text(vjust=1),
        plot.tag.position = c(0, 1),
        plot.tag = element_text(size=10, face='bold', family='Calibri', hjust=0, color='grey7')
    )
}



plot_go_enrichments <- function(enrichments) {
    df_title <- go_enrichments %>% 
        group_by(C) %>% 
        summarize(n=n()) %>% 
        mutate(label = paste0(C, ": ", n, " GO Enrichments"))

    df <- go_enrichments %>% 
    left_join(df_title, by='C') %>% 
    mutate(enrichment = ifelse(direction=='top', -enrichment, enrichment)) %>% 
    mutate(rank_directed = rank * ifelse(direction=='top', -1, 1)) %>% 
    arrange(-rank_directed) %>% 
    mutate(direction = ifelse(direction=='top', 'bottom', 'top')) # flip direction labels so top=positive

    n_label <-  1
    label_size <- 2.5

    df_label <- df %>% filter(FDR_rank <= n_label | 
            description %in% c(
                'Learning or memory',
                'Adaptive immune response',
                'Mitochondrion organization'
                # 'Steroid metabolic process'
                )
            )
    df_selected_labels <- df %>% filter(
            (C=='C3') & (description %in% c(
                # C3 positive
                # 'Learning or memory',
                'Vocal learning',
                # 'Learning or memory',
                # 'Regulation of synaptic plasticity',
                'Regulation of short-term plasticity',
                # C3 negative
                'Immune response to tumor cell',
                'Antigen processing'
                # 'Inflammatory response',
                # 'Regulation of vasculature development',
                # 'Regulation of macrophage migration',
            )) | (C=='C2') & (description %in% c(
                # C2 positive
                # 'Response to jasmonic acid',
                # 'Intestinal cholesterol absorption',
                'Regulation of mitochondrial fusion',
                # 'Generation of metabolites and energy',
                # C2 negative
                'Regulation of chromatin silencing'
                # 'Regulation of epigenetic expression',
                # 'Regulation of mRNA splicing',
            ))
                # C1 positive
                # 'Covalent chromatin modification'
                # C1 negative
                # '7-methylguanosine mRNA capping',
            )

    df %>% 
    ggplot(aes(x=neglogFDR, y=enrichment)) + 
    facet_grid(.~label) +
    geom_point(aes(fill=enrichment, size=n_genes), shape=21, stroke=.3, color='grey') + 
    # Top-n FDR labels
    ## Positive
    geom_text_repel(aes(label=description), data=df_label %>% filter(direction=='top', C=='C1'), 
                    size=label_size, family='Calibri', segment.linetype=1, segment.color='grey',
                    force=10, force_pull=0, xlim=c(2,5), ylim=c(0.5,8), max.time=1, nudge_y=-2.5) +
    geom_text_repel(aes(label=description), data=df_label %>% filter(direction=='top', C=='C2'), 
                    size=label_size, family='Calibri', segment.linetype=1, segment.color='grey',
                    force=10, force_pull=0, xlim=c(1.6,5), ylim=c(0,8), max.time=1, nudge_y=3) +
    geom_text_repel(aes(label=description), data=df_label %>% filter(direction=='top', C=='C3'), 
                    size=label_size, family='Calibri', segment.linetype=1, segment.color='grey',
                    force=10, force_pull=0, xlim=c(2,5), ylim=c(0,8), max.time=1, nudge_y=1.3) +
    ## Negative
    geom_text_repel(aes(label=description), data=df_label %>% filter(direction=='bottom', C=='C1'), 
                    size=label_size, family='Calibri', segment.linetype=1, segment.color='grey',
                    force=10, force_pull=0, xlim=c(2,5), ylim=c(0,-8), max.time=1, nudge_y=.5) + 
    geom_text_repel(aes(label=description), data=df_label %>% filter(direction=='bottom', C=='C2'), 
                    size=label_size, family='Calibri', segment.linetype=1, segment.color='grey',
                    force=10, force_pull=0, xlim=c(1.8,5), ylim=c(0,-8), max.time=1, nudge_y=-1) + 
    geom_text_repel(aes(label=description), data=df_label %>% filter(direction=='bottom', C=='C3'), 
                    size=label_size, family='Calibri', segment.linetype=1, segment.color='grey',
                    force=10, force_pull=0, xlim=c(2,5), ylim=c(0,-8), max.time=1, nudge_y=-1) +
    # Selected other labels
    ## Positive
    geom_text_repel(aes(label=description), data=df_selected_labels %>% filter(direction=='top', C=='C2'), 
                    size=label_size, family='Calibri', segment.linetype=1, segment.color='grey',
                    force=10, force_pull=0, xlim=c(2,5), ylim=c(0,9), max.time=1, nudge_y=1) +
    geom_text_repel(aes(label=description), data=df_selected_labels %>% filter(direction=='top', C=='C3'), 
                    size=label_size, family='Calibri', segment.linetype=1, segment.color='grey',
                    force=10, force_pull=0, xlim=c(2,5), ylim=c(2,9), max.time=1, nudge_y=.5) +
    ## Negative
    geom_text_repel(aes(label=description), data=df_selected_labels %>% filter(direction=='bottom', C=='C2'),
                    size=label_size, family='Calibri', segment.linetype=1, segment.color='grey',
                    force=10, force_pull=0, xlim=c(2,5), ylim=c(-2,-9), max.time=1, nudge_y=-.8) +                    
    geom_text_repel(aes(label=description), data=df_selected_labels %>% filter(direction=='bottom', C=='C3'), 
                    size=label_size, family='Calibri', segment.linetype=1, segment.color='grey',
                    force=10, force_pull=0, xlim=c(2,5), ylim=c(-2,-9), max.time=1, nudge_y=.7, nudge_x=1) +
    # Titles in middle of plot
    # geom_text(aes(label=label), data=df_title, x=-log10(0.001), y=0, size=7.5, 
    #             family='Calibri', color='grey7', fontface='bold') +
    scale_size_continuous(range=c(1,4)) + 
    scale_fill_gradientn(colors=rev(brewer.rdbu(100)), breaks=c(-7,7)) +
    scale_y_continuous(breaks=c(-3,0,3), name='Enrichment z-score') +
    scale_x_continuous(breaks = c(0,1.3,2,3,4), labels = c(0, 0.05, 0.01, 0.001, 0.0001), name='Enrichment FDR') +
    coord_cartesian(clip='off') +
    guides(
        fill=guide_colorbar(title='Enrichment\nz-score', title.vjust=1, order=1, barheight=2, barwidth=.5, ticks=FALSE,
                                                         title.theme = element_text(size=7, family='Calibri', color='grey7', 
                                                                                    margin=margin(0,0,4,0,'mm'), lineheight = 1)),
        size=guide_legend(title='# genes\nin GO term', title.vjust=1, order=2, byrow=T)
    ) +
    # geom_hline(yintercept=0, linetype='dashed', size=.1) + 
    # annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, size=.3, color='grey') +
    # annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, size=.3, color='grey') +
    theme(
        axis.text = element_text(size=6),
        strip.text = element_text(face='bold'),
        axis.title.y = element_text(margin=margin(t=0,b=0,l=0,r=-1, unit='mm')),
        legend.title = element_text(margin=margin(0,0,2,0,'mm')),
        legend.text = element_text(size=6, margin=margin(0,0,0,1,'mm')),
        legend.title.align = 0,
        legend.direction = 'vertical',
        legend.box = 'horizontal',
        legend.position = c(0.1, .85),
        legend.spacing.y = unit(-2.5, 'mm'),
        legend.spacing.x = unit(0, 'mm'),
        # legend.margin = margin(t=-1,b=.2,l=0,r=0,unit='cm'),
        # strip.text=element_blank(),
        panel.grid = element_blank(),
        # panel.grid.major = element_line(size=.1, color='grey'),
        # axis.line = element_line(color='grey', linewidth=.3)
    )
}


plot_enrichment_bars_z <- function(stats, xlab='Enrichment z-score', facet='.~C') {
    
    lim <- max(abs(stats$z))*1.1
    
    stats %>% 
    mutate(C = factor(C, ordered = T, levels = unique(.$C))) %>%
    mutate(label = factor(label, ordered = T, levels = unique(.$label))) %>%
    mutate(sig  =  case_when(
        q < .001 ~ '***',q < .01 ~ '**',q < .05 ~ '*', TRUE ~ ''
        )) %>%
    mutate(hjust = ifelse(z < 0, 1, 0)) %>%
    ggplot(aes(x = z, y = label)) + 
    facet_grid(facet) +
    geom_col(aes(fill = z)) +
    geom_text(aes(label = sig, vjust = .7, hjust = hjust), size = 2.5) +
    scale_y_discrete(limits = rev, name = NULL) +
    scale_x_continuous(limits = c(-lim,lim), 
                       breaks = round(0.5*c(-lim,lim)), name = xlab) +
    scale_fill_gradientn(colors = rev(brewer.rdbu(100)), guide = 'none', limits = c(-lim,lim)) +
    coord_cartesian(clip = 'off')
}

