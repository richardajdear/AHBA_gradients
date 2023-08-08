library(ggseg)
library(ggtext)
library(ggsegGlasser)
library(ggrepel)
# library(ggh4x) # needed for facet_grid2 to not clip strip labels
suppressMessages(library(ggpmisc))
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
    legend.key.size = unit(1, "mm"),
    legend.position = 'bottom',
    strip.text = element_text(size=7, family = 'Calibri', color = 'grey7'),
    strip.background = element_blank(),
    strip.placement='outside',
    plot.title.position='panel',
    plot.title = element_text(size=7, family = 'Calibri', color = 'grey7', face='bold',
                    margin = margin(0,0,0,0,'mm')),
    plot.tag.position = c(0, 1),
    plot.tag = element_text(size=10, face='bold', family='Calibri', hjust=0, color='grey7')
)
colorbar_h <- guide_colorbar(barwidth=2, barheight=.5, ticks=FALSE, title.vjust=1)
colorbar_v <- guide_colorbar(barwidth=.5, barheight=2, ticks=FALSE)



plot_brains <- function(maps, atlas='dk',
                      title="", ncol=3, facet='w', spacing=0,
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
    
    if (atlas=='dk') {
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
    } else if (atlas=='glasser') {
        p <- df %>% 
            ggseg(
                atlas=glasser,
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
        strip.text = element_text(size=7, family='Calibri', color='grey7'),
        strip.clip = 'off',
        legend.position = 'right',
        legend.text = element_text(size=6, family='Calibri', color='grey7'),
        plot.tag.position = c(0, 0.95),
        plot.tag = element_text(size=10, face='bold', family='Calibri', hjust=0, color='grey7')
    )

    if (facet=='h') {
        p + facet_grid(.~map)
    } else if (facet=='v') {
        p + facet_grid(map~., switch='y')
    } else if (facet=='w') {
        p + facet_wrap(~map, ncol=ncol, dir="h")
    }
}



plot_heatmap_maps <- function(null_p, size=4) {
    lim <- max(abs(null_p$r))
    
    null_p %>%
    mutate(map = factor(map, ordered=T, levels=rev(unique(.$map)))) %>%
    mutate(p_level=ifelse(q<0.001, '***',
                   ifelse(q<0.01, '**',
                   ifelse(q<0.05, '*','')))) %>%
    ggplot(aes(x=C, y=map)) +
    geom_raster(aes(fill=r)) + 
    geom_text(aes(label=p_level), vjust=.5, hjust=.5, size=size, color='white') +
    scale_fill_gradientn(colors=rev(brewer.rdbu(100)), name='r', 
                        limits=c(-lim,lim), breaks=c(-floor(lim*10)/10,floor(lim*10)/10)) +
    scale_x_discrete(position='top') +
    guides(fill = colorbar_h) +
    xlab('') + ylab('') +
    coord_fixed()
}


plot_heatmap_enrichments <- function(stats, size=4, limits=NULL) {
    if (is.null(limits)) {
        lim <- max(abs(stats$z))
    } else {
        lim <- limits
    }
    
    stats %>%
    mutate(label = factor(label, ordered=T, levels=rev(unique(.$label)))) %>%
    mutate(p_level=ifelse(q<0.001, '***',
                   ifelse(q<0.01, '**',
                   ifelse(q<0.05, '*','')))) %>%
    ggplot(aes(x=C, y=label)) +
    geom_raster(aes(fill=z)) + 
    geom_text(aes(label=p_level), vjust=.5, hjust=.5, size=size, color='white') +
    scale_fill_gradientn(colors=rev(brewer.rdbu(100)), name='Enrichment z-score', 
                        limits=c(-lim,lim), breaks=c(-floor(lim),floor(lim))
                        ) +
    scale_x_discrete(position='top') +
    guides(fill = colorbar_h) +
    xlab('') + ylab('') +
    coord_fixed() + 
    theme(
        plot.title.position = 'plot'
    )
}



plot_maps_scatter <- function(maps_scatter, maps_scatter_corrs, facet='v', switch='both',
                             x=0, y=0,
                             size=2.6, pointsize=1.2,
                             xlab='', ylab='') {
    corrs <- maps_scatter_corrs %>%
    mutate(map = factor(map, ordered=T, levels=unique(.$map))) %>%
    mutate(sig_label=ifelse(q<0.001, '***',
                   ifelse(q<0.01, '**',
                   ifelse(q<0.05, '*','')))) %>%
    mutate(r_label=paste('r =', round(r,2), sig_label))
    
    maps_scatter %>%
    mutate(map = factor(map, ordered=T, levels=unique(.$map))) %>%
    ggplot(aes(x=C_score, y=map_score)) +
    geom_point(aes(fill=C_score), shape=21, color='grey', stroke=.3, size=pointsize) +
    geom_smooth(method='lm', linetype=1, se=FALSE, color='darkgrey', size=.3) +
    geom_text(data=corrs, aes(label=r_label, x=x, y=y), size=size, hjust=0.5, vjust=0.5) +
    scale_fill_gradientn(colors=rev(brewer.rdbu(100)), guide='none') +
    scale_x_continuous(breaks=0) + 
    scale_y_continuous(breaks=0) +
    xlab(xlab) + ylab(ylab) +
    coord_cartesian(clip='off') +
    theme(
        #   strip.text.x=element_text(),
            strip.text = element_blank(),
            axis.line = element_line(size=.3, color='grey'),
            # panel.grid.major = element_line(),
            # strip.text.y.left=element_text(angle=0, hjust=1, margin(0,0,0.5,0,'mm')),
            axis.text = element_blank(),
            axis.title.y = element_text(angle=0, vjust=.5, margin=margin(0,0,0,0,'mm')),
            axis.title.x = element_text(angle=0, vjust=.5, margin=margin(-4,0,0,0,'mm')),
            # strip.text.x = element_text(margin=margin(-1,0,0,0,'mm')),
            strip.clip = 'off',
            plot.title = element_text(hjust=0.5, vjust=1)
            ) +
    facet_grid(C~map, switch=switch)

}



plot_venn <- function(overlap_deg_gwas, disorder, bottom_margin=1) {
    # title <- paste(disorder,'genes')
    title <- disorder
    disorder <- enquo(disorder)
    counts <- overlap_deg_gwas %>% 
    filter(disorder == !!disorder) %>% 
    group_by(gene) %>% 
    summarise(
        GWAS = sum(label=='GWAS'),
        DEG = sum(label=='DEG')
        ) %>% 
    group_by(GWAS, DEG) %>% 
    count() %>% ungroup() %>% 
    mutate(
        names = case_when(GWAS+DEG==2 ~ 'GWAS&DEG', GWAS==1 ~ 'GWAS', DEG==1 ~ 'DEG'),
        pct = n/sum(n)
    )
    x <- counts$n
    names(x) <- counts$names

    # colors <- c(
    #     brewer.puor(10)[8], brewer.brbg(10)[8], brewer.rdbu(10)[2]
    # )
    # colors <- brewer.puor(10)[c(3,8)]
    # colors <- brewer.piyg(11)[c(2,10)]
    colors <- brewer.brbg(11)[c(3,9)]

    plot <- plot(euler(x), 
                 edges = list(col='grey'),
                 fills = list(fill = colors, alpha = 0.5),
                 labels = list(col = "grey7", fontsize = 7, font='plain', cex=1, fontfamily='Calibri'),
                 quantities = list(col = "grey7", fontsize = 7, fontfamily='Calibri')
                ) %>% 
        wrap_elements +
        ggtitle(title) +
        coord_cartesian(clip='off') +
        theme( 
            # text=element_text(size=20),
            plot.margin = margin(t=0, r=0, b=bottom_margin, l=0, "mm"),
            # plot.title=element_text(hjust=.5, vjust=1, margin=margin(0,0,1,0,'mm'))
            plot.title=element_text(hjust=.5, vjust=.9, margin=margin(0,2,0,0,'mm'), face = 'plain')
        )
    return(plot)
}



plot_disorder_layer_enrichments <- function(
            layer_enrichments, 
            facet = '~disorder', 
            ylab = '% of C3+ GWAS/DEG genes linked to layer', 
            ncol = 1, 
            colors = NULL
            ) {
    # colors <- c(
    #     brewer.puor(10)[8], brewer.brbg(10)[8], brewer.rdbu(10)[2]
    # )
    # colors <- colors[c(1,3)]
    if (is.null(colors)) {
        # colors <- brewer.puor(10)[c(3,8)] %>% rev
        # colors <- brewer.piyg(11)[c(2,10)] %>% rev
        colors <- brewer.brbg(11)[c(3,9)] %>% rev
    }

    position = position_dodge2(reverse=TRUE, width=0.5)

    g <- layer_enrichments %>% 
        mutate(disorder = factor(disorder, ordered=T, levels=unique(.$disorder))) %>% 
        mutate(label = factor(label, ordered=T, levels=unique(.$label))) %>% 
        mutate(p_level=ifelse(q<0.001, '***',
                       ifelse(q<0.01, '**',
                       ifelse(q<0.05, '*','')))) %>%
        ggplot(aes(x=layer, y=pct, group=label)) + 
        geom_linerange(aes(color=label, x=layer, y=0, ymin=0, ymax=pct, alpha=sig), position=position, size=.3) +
        geom_point(aes(fill=label, size=n, alpha=sig), position=position, shape=21, color='grey', stroke=.3) +
        geom_text(aes(label = p_level), position=position, size=2.5, direction='mid', vjust=-.5) + 
        geom_hline(yintercept=0, color='grey', size=.1) +
        scale_alpha_manual(values=c(.5,1)) +
        # scale_alpha_continuous(range=c(.2,1)) +
        # guides(alpha=guide_legend()) +
        scale_fill_manual(values=colors, name=NULL) +
        scale_color_manual(values=colors, name=NULL) +
        scale_size_continuous(range=c(1,4), name='# genes') +
        scale_y_continuous(breaks=c(0,.1,.2), labels=percent) +
        ylab(ylab) +
        xlab('Cortical layer') +
        guides(fill=guide_legend(reverse=T, override.aes = list(size=2)), alpha='none', color='none') +
        coord_cartesian(clip='off') +
        theme(
            axis.text.y = element_text(size=6), 
            strip.text = element_text(margin=margin(5,0,-2,0,'mm')),
            # axis.title.y = element_text(angle=0, vjust=.5),
            legend.text = element_text(size=7),
            legend.position = 'right',
            plot.title = element_text(hjust=0.5, margin=margin(0,0,2,0,'mm')),
            plot.title.position='panel'
        )

    if (facet=='') {
        g
    } else (
        g + facet_wrap(as.formula(facet), ncol=ncol)
    )
}



# plot_disorder_curves <- function(disorder_curves, facet='~disorder', ncol=3) {
#     colors <- c(
#         brewer.puor(10)[8], brewer.brbg(10)[8], brewer.rdbu(10)[2]
#     )
#     colors <- colors[c(1,3)]

#     disorder_curves %>% 
#         ggplot(aes(x=age_log10, y=10**pred)) +
#         facet_wrap(as.formula(facet), scales='fixed', ncol=ncol) +
#         geom_line(aes(color=label, group=label), size=.3) +
#         scale_x_continuous(
#             breaks=log10(c(-0.5,0,1,5,14,40)*365+40*7),
#             labels=function(x) round((10**x - 40*7)/365,2)) +
#         scale_color_manual(values=colors) +
#         guides(color=guide_legend(reverse=T, override.aes = list(size=2)), alpha='none') +
#         ylab('log10 RPKM') +
#         xlab('Age') +
#         coord_cartesian(clip='off') +
#         theme(
#             axis.title.y = element_text(margin=margin(t=0,b=0,l=0,r=-5, unit='mm')),
#             strip.text.x = element_text(margin=margin(0,0,-2,0,'mm')),
#             legend.text = element_text(size=7),
#             legend.position = 'right',
#             plot.title.position='plot'
#         )
# }


