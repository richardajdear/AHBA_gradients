theme_set(theme_classic())
theme_update(
    text = element_text(size=7, family = 'Calibri', color = 'grey7'),
    axis.text = element_text(size=7, family = 'Calibri', color = 'grey7'),
    axis.line = element_line(size=.3, color='darkgrey'),
    axis.ticks = element_blank(),
    legend.key.size = unit(1, "mm"),
    legend.position = 'bottom',
    strip.text = element_text(size=7, family = 'Calibri', color = 'grey7'),
    strip.background = element_blank(),
    strip.placement='outside',
    plot.tag.position = c(0, 1),
    plot.tag = element_text(size=10, face='bold', family='Calibri', hjust=0, color='grey7')
)
theme_colorbar <- guide_colorbar(barwidth=2, barheight=.5, ticks=FALSE)



plot_brain_maps <- function(maps, ncol=1,
                      colors=rev(brewer.rdbu(100)), 
                      colorscale='fixed', limits=c(-3,3),
                      labels=c('-3σ','+3σ'), name=NULL) {
    
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
    
    # set scale limits at 99th percentile
    m_max <- pmax(
        df %>% .$value %>% quantile(.99, na.rm=T) %>% abs,
        df %>% .$value %>% quantile(.01, na.rm=T) %>% abs
    )

    if (colorscale == 'fixed') {
        m_min <- limits[1]
        m_max <- limits[2]
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

    df %>% ggseg(
        atlas=glasser,
        hemi='left',
        mapping=aes(fill=value),
        colour='grey', size=.05,
        show.legend=T
        ) + 
    facet_wrap(~map, ncol=ncol, dir="v") +
    scale_fill_gradientn(
        colors=colors, 
        limits=c(m_min,m_max), oob=squish, breaks=c(m_min,m_max), 
        labels=labels, 
        guide=theme_colorbar,
        name=name
    ) +
    theme_void() + 
    theme(
        text = element_text(size=7, family='Calibri', color='grey7'),
        strip.text = element_text(size=7, family='Calibri', color='grey7'),
        strip.clip = 'off',
        legend.position = 'bottom',
        legend.text = element_text(size=6, family='Calibri', color='grey7'),
        plot.tag.position = c(0, 1),
        plot.tag = element_text(size=10, face='bold', family='Calibri', hjust=0, color='grey7')
        # panel.spacing.x=unit(spacing_x,'lines'),
        # panel.spacing.y=unit(spacing_y,'lines'),
    )
}



plot_corrmat <- function(df, highlight_color='grey50') {
    ### Correlation matrix using continuous axes to support rectangle annotation of modules
    df_plot <- df %>%
    mutate(
        x = factor(x, ordered=T, levels=unique(.$x)),
        y = factor(y, ordered=T, levels=unique(.$y)),
        x1 = as.integer(x),
        y1 = as.integer(y)
    )
    # Set axes labels
    labels <- levels(df_plot$x)
    breaks <- seq_along(labels)

    # Get rectangle boundaries from module annotations of x label
    df_rectangles <- df_plot %>% 
    group_by(cluster) %>% 
    summarise(xmin=min(x1)-0.5, xmax=max(x1)+0.5, ymin=min(x1)-0.5, ymax=max(x1)+0.5)

    # # Get G1-3 for marker lines
    df_lines <- df_plot %>% 
    filter(x %in% c('C1','C2','C3')) %>% select(x, x1)

    df_plot %>% ggplot() + 
    geom_raster(aes(x1,y1, fill=r)) + 
    geom_rect(data=df_rectangles, aes(xmin=xmin, ymin=xmin, xmax=xmax, ymax=xmax),
              fill=NA, size=.3, color=highlight_color) +
    scale_x_continuous(breaks=breaks, labels=labels, limits=c(0.5,length(labels)+0.5)) +
    scale_y_continuous(breaks=breaks, labels=labels, limits=c(0.5,length(labels)+0.5)) +
    coord_cartesian(clip='off') +
    scale_fill_gradientn(colors=rev(brewer.rdbu(100)), limits=c(-1,1), breaks=c(-1,1), 
                         guide=theme_colorbar, name='r',) +
    xlab('') + ylab('') +
    theme(
        aspect.ratio = 1,
        axis.line = element_blank(),
        axis.text = element_text(face = ifelse(df_plot$x %in% c('C1','C2','C3'), 'bold', 'plain')),
        axis.text.x = element_text(angle=90, hjust=1, vjust=.5, margin=margin(t=-10, unit='pt')),
        axis.text.y = element_text(angle=0, margin=margin(r=-10, unit='pt')),
        legend.margin = margin(t=-20, unit='pt'),
        legend.title = element_text(vjust=1),
        legend.position = 'bottom'
    )
}

plot_brain_classes <- function(brain_classes, colors, names) {
    
    brain_classes <- brain_classes %>% 
    mutate(region = recode(region,'7Pl'='7PL')) %>% 
    select(
        colors = {{ colors }}, 
        names = {{ names }},
        region
    ) %>% 
    # mutate(colors = colors %>% saturation(scalefac(1.2))) %>% 
    filter(!is.na(names))

    brain_classes %>% 
    ggseg(
        atlas = glasser,
        hemi = 'left',
        mapping = aes(fill = colors),
        colour = 'grey', 
        alpha = 0.8,
        size = .05,
        show.legend = T
        ) + 
    scale_fill_identity(
        breaks = brain_classes$colors %>% unique,
        labels = brain_classes$names %>% unique,
        na.translate = F, na.value='grey',
        name=NULL,
        guide = 'none'
        # guide = guide_legend(ncol=2)
    ) +
    theme_void() + 
    theme(
        text = element_text(size=7, family='Calibri', color='grey7'),
        strip.text = element_text(size=7, family='Calibri', color='grey7'),
        strip.clip = 'off',
        legend.position = 'bottom',
        legend.text = element_text(size=7, family='Calibri', color='grey7'),
        legend.key.size = unit(2, "mm"),
        plot.title = element_text(size=7, family = 'Calibri', color = 'grey7', face='bold',
                                  hjust=0.5, margin = margin(-0.5,0,1,0,'mm')),
        plot.tag.position = c(0, 1),
        plot.tag = element_text(size=10, face='bold', family='Calibri', hjust=0, color='grey7')
    )
}

plot_class_violins <- function(regions, sig, 
                               classes = Mesulam,
                               colors = Mesulam_colors,
                               names = Mesulam_names,
                            #    legend_positions=c(-6.5, -6.8, -8.0, -7.1)
                            #    legend_positions=c(-5.2,-5.3,-6,-5.8)
                               legend_positions=c(-4.65,-4.8,-5.35,-5) * 1.3
                              ) {
    df_regions <- regions %>% 
    mutate(
        classes = factor( {{ classes }} ),
        colors = {{ colors }}, 
        names = {{ names }}
    ) %>% 
    filter(classes != 'NaN', classes != 8, classes != 0) %>%
    arrange(classes) %>% 
    mutate(names = factor(names, ordered=T, levels=unique(.$names))) %>% 
    select(C1:C3, classes, colors, names) %>% 
    gather(C, score, -classes, -colors, -names) %>% 
    group_by(C, classes, colors, names) %>%
    mutate(
        med = median(score, na.rm=T), 
        mean = mean(score, na.rm=T)
    ) 
    # %>% mutate(colors = colors %>% saturation(scalefac(1.2)) %>% as.character)

    regions_labels <- df_regions %>% 
        ungroup() %>% 
        select(names, colors) %>% 
        unique %>% 
        mutate(ypos=legend_positions %>% rev, C='C1')
    
    sig <- sig %>% 
        mutate(names = as.character(label)) %>% 
        mutate(sig = case_when(q<0.001 ~ '***', q<0.01 ~ '**', q<0.05 ~ '*', TRUE ~ '')) %>% 
        left_join(df_regions %>% select(C, names, med) %>% unique, by=c('C', 'names'))

    df_regions %>% 
    ggplot(aes(x=names, y=score)) + 
    facet_grid(.~C) +
    geom_hline(yintercept=0, color='darkgrey', size=.3) +
    geom_violin(aes(fill=colors), width=.7, color='grey', size=.3) +
    geom_text(data=sig, aes(x=names, y=med, label=sig), size=3, vjust=.75, hjust=.5) +
    geom_point(data=regions_labels, aes(y=ypos, fill=colors), size=4, shape=22, color='grey') +
    scale_fill_identity(guide='none') +
    scale_y_continuous(breaks=c(-2,2), labels=c('-2σ','+2σ'), name=NULL) +
    xlab('') +
    coord_flip(ylim=c(-2.5,2.5), clip='off') +
    theme(
        panel.spacing = unit(1,'mm'),
        axis.line = element_blank(),
        axis.text.x = element_text(size=6),
        plot.margin = margin(0,0,2,5,'mm')
         )
}




plot_scatter_with_colors <- function(data, corrs, 
                                     x_var, y_var, y_name, 
                                     color_var = NULL, 
                                     color_selection = NULL,
                                     with_densities = FALSE,
                                     ylab_adjust=0) {
    
    data <- data %>% 
        filter(!is.na(get(x_var)), !is.na(get(y_var))) %>% 
        arrange(get(color_var)) %>% 
        mutate(
            x=get(x_var), y=get(y_var),
            color = get(paste0(color_var,'_colors')),
            color_names = get(paste0(color_var,'_names'))
        ) %>% 
        mutate(
            color = color %>% saturation(scalefac(2))
        )

    if (!is.null(color_selection)) {
        data <- data %>% mutate(
            color = ifelse(color_names %in% color_selection, color, 'grey')
            )
    }

    corrs <- corrs %>%
        filter(C == !!enquo(x_var), map == !!enquo(y_var)) %>% 
        mutate(p_sig=ifelse(p<0.001, '***',
                    ifelse(p<0.01, '**',
                    ifelse(p<0.05, '*','')))) %>%
        mutate(q_sig=ifelse(q<0.001, '†',
                    ifelse(q<0.01, '†',
                    ifelse(q<0.05, '†','')))) %>%
        mutate(r_label=paste('r =', round(r,2), p_sig, q_sig))
    
    color_labels <- unique(data$color_names)

    scatter <- data %>%  
        ggplot() +
        geom_smooth(aes(x=x, y=y), method='lm', se=F, color='darkgrey', size=.3) +
        geom_point(aes(x=x, y=y, fill=color), size=1.2, stroke=.3, alpha=.8, color='grey50', shape=21) + 
        scale_fill_identity() +
        # geom_point(aes(x=x, y=y, color=color), size=1.2, stroke=1, alpha=.8, fill='darkgrey', shape=21) + 
        # scale_color_identity() +
        # guides(fill = guide_legend(name='')) +
        annotate(geom='text', label=corrs$r_label, 
                 size=2.6, color='grey7', family='Calibri',
                 x=-Inf, y=Inf, vjust=1.5, hjust=-0.1) +
        xlab(x_var) + ylab(y_name) +
        coord_cartesian(clip='off') +
        theme(
            legend.position='bottom',
            axis.text=element_blank(),
            axis.ticks=element_blank(),
            axis.title.y=element_text(angle=0, vjust=.5, margin=margin(t=0,r=ylab_adjust,b=0,l=0,'mm')),
        )

    if (!with_densities) {
        scatter
    } else {
        densities <- data %>% 
            filter(color!='grey') %>% 
            ggplot(aes(x=x, fill=color)) +
            geom_density(alpha=.3, color='grey', size=.2) +
            scale_fill_identity(
                breaks = data$color %>% unique,
                labels = data$names %>% unique,
                guide = 'none'
            ) +
            theme_void()

        densities / scatter + plot_layout(heights=c(1,3))
    }
}



plot_class_densities <- function(classes_with_scores, C, colors, names) {
    classes_with_scores <- classes_with_scores %>% 
    mutate(
        C = {{ C }},
        colors = {{ colors }}, 
        names = {{ names }}
    ) %>% 
    filter(!is.na(names))

    classes_with_scores %>% 
    ggplot(aes(x=C, fill=colors)) +
    geom_density(alpha=.3, color='grey', size=.2) +
    scale_fill_identity(
        breaks = classes_with_scores$colors %>% unique,
        labels = classes_with_scores$names %>% unique,
        na.translate = F, na.value='grey',
        guide = 'none'
    ) +
    theme_void() + 
    theme(
        legend.text=element_text(size=7, color='grey7', family='Calibri'),
        legend.key.height = unit(3, "mm"),
        legend.key.width = unit(3, "mm"),
    ) + plot_layout(tag_level='new')
}

# plot_scatter_with_colors <- function(data, corrs, x_name, x_var, y_var, color_var=NULL, left_margin=0) {
    
#     data <- data %>% 
#         filter(!is.na(get(x_var)), !is.na(get(y_var))) %>% 
#         arrange(get(color_var)) %>% 
#         mutate(
#             x=get(x_var), y=get(y_var),
#             color = get(paste0(color_var,'_colors')),
#             color_names = get(paste0(color_var,'_names'))
#         ) %>% 
#         mutate(
#             color = factor(color, ordered=T, levels=unique(.$color))
#         )

#     corrs <- corrs %>%
#         filter(map == !!enquo(x_var), C == !!enquo(y_var)) %>% 
#         mutate(p_sig=ifelse(p<0.001, '***',
#                     ifelse(p<0.01, '**',
#                     ifelse(p<0.05, '*','')))) %>%
#         mutate(q_sig=ifelse(q<0.001, '†',
#                     ifelse(q<0.01, '†',
#                     ifelse(q<0.05, '†','')))) %>%
#         mutate(r_label=paste('r =', round(r,2), p_sig, q_sig),
#                 # '\np = ', round(p,3), p_sig
#                 # '\nq = ', round(q,3), q_sig
#                 ) %>% 
#         mutate(
#             x = -Inf,
#             y = ifelse(r > 0, Inf, -Inf),
#             vjust = ifelse(r > 0, 1.5, -1)
#             )
    
#     color_labels <- unique(data$color_names)

#     scatter <- data %>%  
#         ggplot() +
#         geom_smooth(aes(x=x, y=y), method='lm', se=F, color='black', size=.3) +
#         geom_point(aes(x=x, y=y, fill=color), size=1.2, stroke=.3, alpha=.8, color='darkgrey', shape=21) + 
#         scale_fill_identity() +
#         guides(color=F) +
#         annotate(geom='text', label=corrs$r_label, size=2.6, color='grey7', family='Calibri',
#                  x=corrs$x, y=corrs$y, vjust=corrs$vjust, hjust=-.1) +
#         xlab(x_name) + ylab(y_var) +
#         theme(
#             axis.text=element_blank(),
#             axis.ticks=element_blank(),
#             axis.title.y=element_text(angle=0, vjust=.5, margin=margin(t=0,r=0,b=0,l=left_margin,'cm')),
#         )

#     densities <- data %>% 
#         ggplot(aes(x=y, fill=color)) +
#         geom_density(alpha=.3, color='grey', size=.2) +
#         scale_fill_identity(name=NULL, labels=color_labels, 
#             guide=guide_legend(reverse=T)) + 
#         coord_flip() + 
#         theme_void() + 
#         theme(
#             legend.text=element_text(size=7, color='grey7', family='Calibri'),
#             legend.key.height = unit(3, "mm"),
#             legend.key.width = unit(3, "mm"),
#         ) + plot_layout(tag_level='new')

#     scatter + densities + plot_layout(widths=c(3,1))
# }
