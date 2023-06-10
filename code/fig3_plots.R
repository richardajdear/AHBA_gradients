theme_set(theme_classic())
theme_update(
    text = element_text(size=7, family = 'Calibri', color = 'grey7'),
    axis.text = element_text(size=7, family = 'Calibri', color = 'grey7'),
    axis.line = element_line(size=.2),
    axis.ticks = element_blank(),
    legend.key.size = unit(1, "mm"),
    plot.tag.position = c(0,1),
    plot.tag = element_text(size=7, face='bold', family='Calibri', hjust=0)
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

    p <- df %>% ggseg(
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
        legend.text = element_text(size=6, family='Calibri', color='grey7')
        # panel.spacing.x=unit(spacing_x,'lines'),
        # panel.spacing.y=unit(spacing_y,'lines'),
    )
}


plot_scatter_with_colors <- function(data, corrs, x_name, x_var, y_var, color_var=NULL, left_margin=0) {
    
    data <- data %>% 
        filter(!is.na(get(x_var)), !is.na(get(y_var))) %>% 
        arrange(get(color_var)) %>% 
        mutate(
            x=get(x_var), y=get(y_var),
            color = get(paste0(color_var,'_colors')),
            color_names = get(paste0(color_var,'_names'))
        ) %>% 
        mutate(
            color = factor(color, ordered=T, levels=unique(.$color))
        )

    corrs <- corrs %>%
        filter(map == !!enquo(x_var), G == !!enquo(y_var)) %>% 
        mutate(p_sig=ifelse(p<0.001, '***',
                    ifelse(p<0.01, '**',
                    ifelse(p<0.05, '*','')))) %>%
        mutate(q_sig=ifelse(q<0.001, '†',
                    ifelse(q<0.01, '†',
                    ifelse(q<0.05, '†','')))) %>%
        mutate(r_label=paste('R =', round(r,2), p_sig, q_sig),
                # '\np = ', round(p,3), p_sig
                # '\nq = ', round(q,3), q_sig
                ) %>% 
        mutate(
            x = -Inf,
            y = ifelse(r > 0, Inf, -Inf),
            vjust = ifelse(r > 0, 1.5, -1)
            )
    
    color_labels <- unique(data$color_names)

    scatter <- data %>%  
        ggplot() +
        geom_smooth(aes(x=x, y=y), method='lm', se=F, color='black', size=.3) +
        geom_point(aes(x=x, y=y, fill=color), size=1.2, stroke=.3, alpha=.8, color='darkgrey', shape=21) + 
        scale_fill_identity() +
        guides(color=F) +
        annotate(geom='text', label=corrs$r_label, size=2.6, color='grey7', family='Calibri',
                 x=corrs$x, y=corrs$y, vjust=corrs$vjust, hjust=-.1) +
        xlab(x_name) + ylab(y_var) +
        theme(
            axis.text=element_blank(),
            axis.ticks=element_blank(),
            axis.title.y=element_text(angle=0, vjust=.5, margin=margin(t=0,r=0,b=0,l=left_margin,'cm')),
        )

    densities <- data %>% 
        ggplot(aes(x=y, fill=color)) +
        geom_density(alpha=.3, color='grey', size=.2) +
        scale_fill_identity(name=NULL, labels=color_labels, 
            guide=guide_legend(reverse=T)) + 
        coord_flip() + 
        theme_void() + 
        theme(
            legend.text=element_text(size=7, color='grey7', family='Calibri'),
            legend.key.height = unit(3, "mm"),
            legend.key.width = unit(3, "mm"),
        ) + plot_layout(tag_level='new')

    scatter + densities + plot_layout(widths=c(3,1))
}



plot_corrmat <- function(df, highlight_color='grey7') {
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
    filter(x %in% c('G1','G2','G3')) %>% select(x, x1)

    df_plot %>% ggplot() + 
    geom_raster(aes(x1,y1, fill=r)) + 
    geom_rect(data=df_rectangles, aes(xmin=xmin, ymin=xmin, xmax=xmax, ymax=xmax),
              fill=NA, size=.2, color=highlight_color) +
    scale_x_continuous(breaks=breaks, labels=labels, limits=c(0.5,length(labels)+0.5)) +
    scale_y_continuous(breaks=breaks, labels=labels, limits=c(0.5,length(labels)+0.5)) +
    coord_cartesian(clip='off') +
    scale_fill_gradientn(colors=rev(brewer.rdbu(100)), limits=c(-1,1), breaks=c(-1,1), 
                         guide=theme_colorbar, name='r',) +
    xlab('') + ylab('') +
    theme(
        aspect.ratio = 1,
        axis.line = element_blank(),
        axis.text = element_text(face = ifelse(df_plot$x %in% c('G1','G2','G3'), 'bold', 'plain')),
        axis.text.x = element_text(angle=90, hjust=1, vjust=.5, margin=margin(t=-10, unit='pt')),
        axis.text.y = element_text(angle=0, margin=margin(r=-10, unit='pt')),
        legend.margin = margin(t=-20, unit='pt'),
        legend.title = element_text(vjust=1),
        legend.position = 'bottom'
    )
}
