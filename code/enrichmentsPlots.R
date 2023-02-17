suppressPackageStartupMessages(library(tidyverse))
library(pals)
library(ggridges)
library(ggrepel)
library(ggradar)


plot_enrichment_bars <- function(null_p, xlab='Corr', lim='none') {
    
    if (lim=='none') {
        lim <- max(abs(null_p$r))
    }
    
    null_p %>% 
    mutate(sig = case_when(
        q < .001 ~ '***',q < .01 ~ '**',q < .05 ~ '*', TRUE ~ ''
        )) %>%
    mutate(hjust=ifelse(r < 0, 1, 0)) %>%
    ggplot(aes(x=r, y=label)) + 
    facet_grid(.~G) +
    geom_col(aes(fill=r)) +
    geom_text(aes(label=sig, vjust=.7, hjust=hjust), size=6) +
    # scale_alpha_manual(values=c(0.2,1)) +
    geom_hline(yintercept=0) +
    scale_y_discrete(limits=rev, name='') +
    scale_x_continuous(limits=1.2*c(-lim,lim), breaks=round(0.5*c(-lim,0,lim),1), name=xlab) +
    scale_fill_gradientn(colors=rev(brewer.rdbu(100)), guide='none', limits=c(-lim,lim)) +
    # ggtitle("Mean weight of cell genes vs permutations") +
    # coord_flip() +
    theme_minimal() +
    theme(panel.grid.minor=element_blank(),
          plot.title=element_text(size=20),
          axis.text=element_text(size=20),
          text=element_text(size=20)
         )
}


plot_enrichment_bars_z <- function(null_p, xlab='z-score') {
    
    lim <- max(abs(null_p$z))
    
    null_p %>% 
    mutate(G=factor(G, ordered=T, levels=unique(.$G))) %>%
    mutate(label=factor(label, ordered=T, levels=unique(.$label))) %>%
    mutate(sig = case_when(
        q < .001 ~ '***',q < .01 ~ '**',q < .05 ~ '*', TRUE ~ ''
        )) %>%
    mutate(hjust=ifelse(z < 0, 1, 0)) %>%
    ggplot(aes(x=z, y=label)) + 
    facet_grid(.~G) +
    # geom_vline(xintercept=0) +
    geom_col(aes(fill=z)) +
    geom_text(aes(label=sig, vjust=.7, hjust=hjust), size=8) +
    # scale_alpha_manual(values=c(0.2,1)) +
    scale_y_discrete(limits=rev, name='') +
    scale_x_continuous(limits=c(-lim,lim), breaks=round(0.5*c(-lim,lim)), name=xlab) +
    scale_fill_gradientn(colors=rev(brewer.rdbu(100)), guide='none', limits=c(-lim,lim)) +
    # ggtitle("Mean weight of cell genes vs permutations") +
    # coord_flip() +
    coord_cartesian(clip='off') +
    theme_minimal() +
    theme(panel.grid.minor=element_blank(),
          plot.title=element_text(size=20),
          axis.text=element_text(size=20, color='grey7'),
          strip.text=element_text(size=20, color='grey7'),
          text=element_text(size=20,color='grey7')
         )
}



plot_enrichment_heatmaps <- function(null_p_versions, ncol=3, scales=NULL) {
    lim <- max(abs(null_p_versions$z))
    
    null_p_versions %>%
    mutate(sig_label = case_when(
        p < .001 ~ '***',p < .01 ~ '**',p < .05 ~ '*', TRUE ~ ''
        )) %>%
    mutate_at(vars(label), ~ factor(., levels=rev(unique(.)))) %>% 
    mutate_at(vars(version), ~ factor(., levels=unique(.))) %>%         
    ggplot(aes(x=G, y=label)) + 
    facet_wrap(how~version, ncol=ncol, scales=scales) +
    geom_tile(aes(fill=z, color=sig), size=2) +
    geom_text(aes(label=paste(round(z, digits = 2), '\n', round(p, digits=3), sig_label)), size=8) +
    scale_color_manual(values=c('transparent','green'), name='FDR sig') +
    scale_fill_gradientn(colours=rev(brewer.rdbu(100)[20:80]), guide='colourbar', limits=c(-lim,lim)) +
    guides(fill=guide_colourbar(title='z-score of true mean gene weight vs random permutations', barheight=10)) +
    # scale_x_discrete(position = "top") +
    theme_minimal() + 
    theme(panel.spacing=unit(1,'lines'), 
          axis.title = element_blank(),
          strip.text.x=element_text(size=20),
          strip.placement='outside',
          aspect.ratio=1,
          text=element_text(size=20)
         )
    # coord_fixed()
}


plot_enrichment_heatmaps_2 <- function(null_p_versions, ncol=3, fill_name='pearson r') {
    null_p_versions %>%
    mutate(sig_label = case_when(
        p < .001 ~ '***',p < .01 ~ '**',p < .05 ~ '*', TRUE ~ ''
        )) %>%
    mutate_at(vars(label), ~ factor(., levels=rev(unique(.)))) %>% 
    mutate_at(vars(version), ~ factor(., levels=unique(.))) %>%     
    ggplot(aes(x=G, y=label)) + 
    facet_wrap(~version, ncol=ncol) +
    geom_tile(aes(fill=true_mean, color=sig), size=1) +
    geom_text(aes(label=paste(round(true_mean, digits = 2), '\n', round(p, digits=3), sig_label)), size=8) +
    scale_color_manual(values=c('transparent','green'), name='FDR sig') +
    scale_fill_gradientn(colours=rev(brewer.rdbu(100)[20:80]), guide='colourbar') +
    guides(fill=guide_colourbar(title=fill_name, barheight=10)) +
    scale_x_discrete(position = "top") +
    theme_minimal() + 
    theme(panel.spacing=unit(1,'lines'), 
          axis.title = element_blank(),
          strip.text.x=element_text(size=20),
          strip.placement='outside',
          text=element_text(size=20)
         ) +
    coord_fixed()
}




plot_weight_scatters_with_labels <- function(weights_labels, title='', colors=brewer.set1(10), drop='') {

    weights <- weights_labels %>% rename(G1=`0`, G2=`1`, G3=`2`) %>%
                select(label, contains('cluster'), contains('log2FC'), G1:G3)
    df <- rbind(
        weights %>% rename(x=G1, y=G2) %>% select(-G3) %>% mutate(which='G1 v G2'),
        weights %>% rename(x=G1, y=G3) %>% select(-G2) %>% mutate(which='G1 v G3'),
        weights %>% rename(x=G2, y=G3) %>% select(-G1) %>% mutate(which='G2 v G3')
    ) %>% 
    filter(!which %in% drop) %>%
    arrange(label, which) %>%
    mutate(which = factor(which, ordered=T, levels=unique(.$which)),
           label = factor(label, ordered=T, levels=unique(.$label)),
          )

    lim <- max(abs(weights %>% pivot_longer(G1:G3, names_to='G', values_to='weight') %>% .$weight))
    
    df_axs <- rbind(
        data.frame(x=c(-lim, 0), y=c(0, lim), ax=c('G1', 'G2'), which='G1 v G2', hjust=c(.5,1.2), vjust=c(-0.3,0)),
        data.frame(x=c(-lim, 0), y=c(0, lim), ax=c('G1', 'G3'), which='G1 v G3', hjust=c(.5,1.2), vjust=c(-0.3,0)),
        data.frame(x=c(-lim, 0), y=c(0, lim), ax=c('G2', 'G3'), which='G2 v G3', hjust=c(.5,1.2), vjust=c(-0.3,0))
    ) %>%
    filter(!which %in% drop) %>%    
    mutate(which = factor(which, ordered=T, levels=unique(.$which)),
           ax = factor(ax, ordered=T, levels=unique(.$ax)),
          )
    
    if ('cluster' %in% colnames(df)) {
        df <- df %>% mutate(colors=cluster)
    } else if (!'log2FC' %in% colnames(df)) {
        df <- df %>% mutate(colors=label)
    } else {
        df <- df %>% mutate(colors=log2FC)
        color_p95 <- quantile(abs(weights$log2FC), .95, na.rm=T)
    }
    
    p <- ggplot(df%>%filter(label!='none')) + 
    facet_grid2(which~label, strip=strip_vanilla(clip='off')) + 
    geom_hline(yintercept=0, color='grey7') + geom_vline(xintercept=0, color='grey7') +
    geom_text(data=df_axs, aes(x=x, y=y, label=ax, hjust=hjust, vjust=vjust), size=6.5, color='grey7') +
    geom_point(data=df%>%select(-label), aes(x,y), size=.1, color='lightgrey', alpha=.1) +
    geom_point(aes(x,y, color=colors), size=.5, alpha=.3) +
    scale_x_continuous(limits=c(-lim,lim)) + 
    scale_y_continuous(limits=c(-lim,lim)) +
    coord_cartesian(clip='off') +
    theme_void() + 
    ggtitle(title) +
    theme(aspect.ratio=1,
          text = element_text(size=20, color='grey7'),
          panel.spacing=unit(2,'lines'),
          strip.text.y = element_blank(),
          strip.text.x.top = element_text(size=20, vjust=3, color='grey7'),
          # strip.placement='outside', 
        #   strip.clip='off',
          # plot.margin = margin(t=10),
          # strip.placement='outside',
          legend.position=c(.5,-0.06),
          legend.text=element_text(size=20, color='grey7'),
          plot.title=element_text(size=24,hjust=0.5,vjust=3)
         )
    
    if ('cluster' %in% colnames(df)) {
        p + scale_color_manual(values=colors, name='')
    } else if ('log2FC' %in% colnames(df)) {
        p + scale_color_gradientn(
            colors=rev(brewer.rdbu(100)), 
            # colors=cubicl(100), 
            limits=c(-color_p95, color_p95), name='log2FC') + guides(color=guide_colorbar(barwidth=10)) + theme(legend.title=element_text(vjust=1))
    } else {
        p + scale_color_manual(values=colors)
    }
    
}







plot_go_enrichments <- function(enrichments, size=4, name='Enrichment FDR') {
    df <- enrichments %>%
    mutate(neglogFDR_fill = ifelse(direction=='top', -neglogFDR, neglogFDR)) %>%
    mutate(rank = ifelse(direction=='top', -rank, rank)) %>%
    # mutate(hjust = ifelse(direction=='top', 0, 1)) %>%
    arrange(neglogFDR)
    
    lim <- max(abs(df$neglogFDR))

    df %>%
    ggplot() + 
    facet_wrap(~G) +
    coord_flip(clip='off') +
    # geom_hline(yintercept=-log10(0.05), linetype=2, color=brewer.rdbu(5)[1]) +
    geom_col(aes(x=rank, y=neglogFDR, fill=neglogFDR_fill), alpha=.2) +
    geom_text(aes(x=rank, y=0, label=description, hjust=0), size=size, family='Calibri') +
    theme_minimal() +
    ylab(name) + xlab('') +
    # xlab('GO Biological\nProcess\nEnrichment') +
    scale_x_continuous(breaks=c(-5,5),labels=c('-','+')) +
    scale_y_continuous(breaks=c(0,1.3,3,4), labels=c(0, 0.05, 0.001, 0.0001)) +
    scale_fill_gradientn(colors=rev(brewer.rdbu(200)), guide='none',
                         limits=c(-lim,lim)) +
    theme(
          # axis.text.y=element_blank(),
          axis.text.y=element_text(size=20, color='grey7'),
          axis.text.x=element_text(size=20, color='grey7'),
          axis.line.y = element_line(arrow = grid::arrow(length = unit(0.3, "cm"), ends = "both")),
          axis.title.y = element_text(angle=0, vjust=0.5, color='grey7'),
          panel.grid.minor=element_blank(),
          panel.grid.major.y=element_blank(),
          panel.spacing=unit(10,'lines'),
          strip.placement='outside',
        #   strip.text=element_text(size=20, margin=margin(b=-10, unit='pt')),
          legend.position=c(0.1,0.8),
          text=element_text(size=20)
         )
}


plot_enrichments <- function(enrichments, size=4) {
    df <- enrichments %>%
    mutate(neglogFDR = ifelse(direction=='top', -neglogFDR, neglogFDR)) %>%
    mutate(rank = ifelse(direction=='top', -rank, rank)) %>%
    mutate(hjust = ifelse(direction=='top', 0, 1)) %>%
    arrange(neglogFDR)
    
    lim <- max(abs(df$neglogFDR))

    df %>%
    ggplot() + 
    facet_wrap(~G, switch='x') +
    coord_flip(clip='off') +
    geom_col(aes(x=rank, y=neglogFDR, fill=neglogFDR), alpha=.5) +
    geom_text(aes(x=rank, y=0, label=description, hjust=hjust), size=size) +
    theme_minimal() +
    ylab('-log(FDR)') + xlab('') +
    scale_y_continuous(breaks=c(-2,0,2,4), labels=c(0.01, 0, 0.01, 0.0001)) +
    scale_fill_gradientn(colors=rev(brewer.rdbu(200)), guide='none',
                         limits=c(-lim,lim)) +
    theme(axis.text.y=element_blank(),
          panel.grid=element_blank(),
          panel.spacing=unit(10,'lines'),
          strip.placement='outside',
          legend.position=c(0.1,0.8),
          text=element_text(size=20)
         )
}


plot_log2FC <- function(scatter, corr=NULL, x=0, y=1) {
    p <- scatter %>%
    mutate(label = factor(label, ordered=T, levels=rev(unique(.$label)))) %>%
    mutate(updown = ifelse(log2FC>0,'up','down')) %>%
    ggplot(aes(x=weight, y=log2FC)) + 
    facet_grid(label~G, switch='x') +
    geom_vline(xintercept=0, size=.2, linetype=1) +
    geom_point(alpha=.3, size=.5, aes(color=updown)) +
    scale_color_manual(values=c(brewer.rdbu(10)[c(9,2)])) +
    geom_smooth(method='lm', color='black') +
    theme_minimal() +
    theme(
        aspect.ratio=1,
        strip.placement='outside',
        panel.grid.minor=element_blank(),
        panel.border=element_rect(fill=NA),
        text=element_text(size=16),
        axis.title.y=element_text(size=20),
        axis.title.x=element_blank(),
        strip.text.y=element_text(angle=0),
        strip.text=element_text(size=20),
        legend.position = 'bottom',
        legend.text=element_text(size=20),
        legend.title=element_blank(),
        plot.title=element_text(size=22)
    )

    if (!is.null(corr)) {
        corr = corr %>%
        mutate(sig = case_when(q < .001 ~ '***',q < .01 ~ '**',q < .05 ~ '*', TRUE ~ '')) %>%
        mutate(text=paste('r =', round(true_mean,2), sig))

        p + geom_text(data=corr, aes(x=x, y=y, label=text), size=8)
    } else {
        p
    }
}

plot_gene_distributions <- function(dge, null_p=FALSE, split=FALSE, abs=FALSE) {
    if (null_p) {
        dge <- dge %>% 
            mutate(sig = ifelse(q<0.05,'FDR sig','not FDR sig')) %>% 
            mutate(p_sig = case_when(
                p < .0001 ~ '****',p < .001 ~ '***',p < .01 ~ '**',p < .05 ~ '*', TRUE ~ ''
            )) %>%
            mutate(q_sig = case_when(
                q < .0001 ~ '****',q < .001 ~ '***',q < .01 ~ '**',q < .05 ~ '*', TRUE ~ ''
            )) %>% 
            mutate(label = paste0(label, "\n(", n_matches,"/",n_genes, ")"))
    }

    if (split) {
        dge <- dge %>% mutate(updown = ifelse(log2FC>0,'up','down'))
    }

    if (abs) {
        dge <- dge %>% mutate(weight = abs(weight))
    }
    
    p <- dge %>% 
    mutate(label = factor(label, ordered=T, levels=rev(unique(.$label)))) %>%
    group_by(label, G) %>% 
    mutate(mean_weight = mean(weight)) %>% 
    ungroup() %>% 
    ggplot(aes(x=weight)) +
    facet_grid(label~G, switch='both') +
    # geom_vline(aes(xintercept=mean_weight), size=.5, linetype=2) +
    # geom_vline(aes(xintercept=0), size=.5, linetype=1) +
    theme_minimal() +
    theme(
        aspect.ratio=1,
        strip.placement='outside',
        panel.grid=element_blank(),
        panel.border=element_rect(fill=NA),
        text=element_text(size=16),
        # axis.title.y=element_text(size=20),
        axis.text=element_blank(),
        axis.ticks=element_blank(),
        axis.title=element_blank(),
        strip.text.y.left=element_text(angle=0),
        strip.text=element_text(size=20),
        # legend.position='none',
        legend.position = 'bottom',
        legend.text=element_text(size=20),
        legend.title=element_blank(),
        plot.title=element_text(size=22)
    )

    if (split) {
        p <- p + 
            geom_density(aes(y=after_stat(count), fill=updown), alpha=.3, size=.5) +
            geom_vline(aes(xintercept=mean_weight), linetype=2, size=1) +
            scale_fill_manual(values=c(brewer.rdbu(10)[c(9,2)]))
    } else if (!null_p) {
        p  <- p + 
            geom_density(aes(y=after_stat(count), fill=mean_weight), alpha=.3, size=.5) +
            scale_fill_gradientn(colors=rev(brewer.rdbu(100)), limits=c(-.005,.005))
    } else {
        data_for_annotate <- dge %>% 
            select(label,G,true_mean,null_mean,z,p,q,p_sig,q_sig) %>% 
            unique %>% 
            mutate(label = factor(label, ordered=T, levels=rev(unique(.$label))))

        p  <- p + 
            geom_density(aes(fill=z, color=factor(sig)), alpha=.3, size=1.5) +
            scale_color_manual(values=c('green','darkgrey')) +
            scale_fill_gradientn(colors=rev(brewer.rdbu(100)), limits=c(-10,10)) +
            geom_vline(data=data_for_annotate, aes(xintercept=null_mean), linetype=1, size=1) +
            geom_vline(data=data_for_annotate, aes(xintercept=true_mean), linetype=2, size=1) + 
            geom_text(data=data_for_annotate, x=0.005, y=30, hjust=0, size=6,
                    aes(label=paste0("z=",round(z,2),
                                    "\np=",round(p,2),p_sig,
                                    "\nq=",round(q,2),q_sig)
                    ))
    }

    p
}



# ########### OLD


# plot_cell_radar <- function(null_p, size=4) {
#     null_p %>% 
#     select(celltype, G, z) %>%
#     spread(celltype, z) %>%
#     ggradar(
#             grid.min=-20, grid.max=20,
#             values.radar = c("z = -20", "z = 0", "z = +20"),
#             group.point.size = 3,
#             base.size = 20,
#             axis.label.size = size,
#             grid.label.size = size,
#             fill=T,
#             fill.alpha=.05,
#             group.colours=brewer.rdylbu(5)[c(1,2,5)],
#             background.circle.colour = "white",
#             gridline.mid.colour = "grey",
#             legend.text.size = 20,
#             legend.position = "left"
#     )
# }


# plot_cell_gene_corrs <- function(corrs, sort=CT, ncol=3) {
#     corrs %>%
#     mutate(cell_type = factor(cell_type, ordered=T, levels=unique(.$cell_type))) %>%
#     arrange({{ sort }}) %>%
#     # mutate(gene = factor(gene, ordered=T, levels=unique(.$gene))) %>%
#     gather(map, r, -cell_type, -rank) %>%
#     mutate(map = factor(map, ordered=T, levels=unique(.$map))) %>%
#     ggplot() +
#     facet_wrap(~cell_type, ncol=ncol, scales='free') +
#     geom_raster(aes(map, rank, fill=r)) +
#     scale_fill_gradientn(colors=rev(brewer.rdbu(100))) +
#     ylab('ranked genes') +
#     theme_minimal() +
#     theme(text=element_text(size=15),
#           axis.text.y=element_blank(),
#           axis.text.x=element_text(size=10),
#           strip.text.x=element_text(size=20)
#          )
# }


# plot_cell_maps_scatters <- function(cell_maps_scatters, cell_maps_corrs) {
#     cell_maps_corrs <- cell_maps_corrs %>% 
#     rownames_to_column('cell_type') %>%
#     mutate(cell_type = factor(cell_type, ordered=T, levels=unique(.$cell_type))) %>%
#     gather(G, r, -cell_type) %>%
#     mutate(G = factor(G, ordered=T, levels=unique(.$G))) %>%
#     mutate(label = paste('r =', round(r, 3)))

#     cell_maps_scatters %>%
#     mutate(cell_type = factor(cell_type, ordered=T, levels=unique(.$cell_type))) %>%
#     mutate(G = factor(G, ordered=T, levels=unique(.$G))) %>%
#     ggplot(aes(x=g_score, y=cell_score)) + 
#     facet_grid(cell_type~G, switch='y') + 
#     geom_point(alpha=.5) +
#     geom_smooth(method='lm', color='skyblue') +
#     geom_text(data=cell_maps_corrs, x=1, y=.65, aes(label=label), size=6) +
#     xlab('') + ylab('') +
#     theme_minimal() +
#     theme(aspect.ratio=1, 
#           strip.placement='outside', 
#           strip.text.y.left=element_text(angle=0, size=20),
#           strip.text.x=element_text(size=20),
#           panel.grid=element_blank(),
#           panel.border=element_rect(color='grey', fill='transparent')
#          )
# }

# plot_cell_gene_scatters <- function(cell_gene_scatters, cell_gene_corrs) {
#     cell_gene_corrs <- cell_gene_corrs %>% 
#     mutate(cell_type = factor(cell_type, ordered=T, levels=unique(.$cell_type))) %>%
#     select(-gene) %>%
#     gather(G, r, -cell_type) %>%
#     mutate(label = paste('r =', round(r, 3)))

#     cell_gene_scatters %>%
#     mutate(cell_type = factor(cell_type, ordered=T, levels=unique(.$cell_type))) %>%
#     ggplot(aes(x=map_corr, y=weight)) + 
#     facet_grid(cell_type~G, switch='y') + 
#     geom_point(alpha=.5) +
#     geom_smooth(method='lm', color='purple') +
#     # geom_text(data=cell_gene_corrs, x=-0.5, y=.05, aes(label=label), size=6) +
#     xlab('') + ylab('') +
#     theme_minimal() +
#     theme(aspect.ratio=1, 
#           strip.placement='outside', 
#           strip.text.y.left=element_text(angle=0, size=20),
#           strip.text.x=element_text(size=20),
#           panel.grid=element_blank(),
#           panel.border=element_rect(color='grey', fill='transparent')
#          )
# }




# plot_cell_enrichment <- function(true_scores, null_scores, null_p, how='mean', p_sig = .05/7) {

#     mycolors = c(brewer.rdylbu(6)[1:3],brewer.rdylbu(5)[4:5])

#     null_scores <- null_scores %>% gather(pc, score, -m)
#     null_p <- null_p %>% rownames_to_column('m') %>% gather(pc, p, -m) %>% mutate(sig = ifelse(p < p_sig | p > 1-p_sig, T, F)) 
#     null_scores <- null_p %>% left_join(null_scores, by=c('pc', 'm')) %>% mutate(m = factor(m, levels=rev(.$m %>% unique)))
#     true_scores <- true_scores %>% rownames_to_column('m') %>% gather(pc, score, -m) %>% left_join(null_p, by=c('pc', 'm'))


#     p <- null_scores %>% 
#     ggplot() + 
#     facet_grid(.~pc) +
#     stat_density_ridges(geom = "density_ridges_gradient", calc_ecdf = TRUE, scale=.8, size=.2, color='grey',
#                        aes(x=score, y=m, fill=stat(ecdf))) +
#     scale_fill_gradientn(colors=rev(brewer.rdbu(100)[15:85]), name = 'Quantile') +
#     geom_point(data=true_scores, aes(x=score, y=m, fill=p), shape=21, size=15) +
#     geom_text(data=true_scores, aes(label=paste0('p=',round(p,2)), x=median(score), y=m), 
#               size=12, hjust=-.2, nudge_y=.3) +
#     guides(fill = guide_colorbar(barwidth=20, title.vjust=1))
    
#     if (how == 'mean') {
#         p <- p + xlab('Mean weight of genes associated with cell-type')
#     } else if (how == 'median') {
#         p <- p + xlab('Median rank of genes associated with cell-type')
#     }
        
#     p +
#     theme_void() + 
#     theme(
#         axis.title.x = element_text(size=36, vjust=5), 
#         axis.text.y = element_text(size=36), 
#         strip.text.y.left = element_text(size=36, angle=0), 
#         # strip.text.x = element_text(size=30, angle=0), 
#         strip.text.x = element_blank(), 
#         panel.spacing=unit(20, 'lines'), 
#         text = element_text(size=30),
#         legend.position='bottom'
#     )
# }



# library(ggrepel)

# plot_revigo_data <- function(revigo_R, text_threshold=.15, title='') {
#     source(revigo_R)
    
#     ex <- one.data [ one.data$dispensability < text_threshold, ];
    
#     p1 <- ggplot( data = one.data ) + 
#         # geom_vline(xintercept=0, colour = I (alpha ("black", 0.6) )) + 
#         # geom_hline(yintercept=0, colour = I (alpha ("black", 0.6) )) +
#         geom_point( aes(plot_X, plot_Y, size = log_size, fill=value), shape = 21, colour = I (alpha ("black", 0.6) ), alpha=.6) + 
#         geom_label_repel(data = ex, aes(plot_X, plot_Y, label = description), colour = I(alpha("black", 0.85)), 
#                          size = 8, box.padding = 1, force=5) +
#         scale_size(range=c(5,30), name='GO terms') +
#         # scale_fill_gradientn(colours = brewer.ylorrd(100)) +
#         scale_fill_gradientn(colours = brewer.reds(100)) +
#         guides(fill = guide_colorbar(barwidth=20, title='Log p-value', title.vjust=1), size='none') +
#         # coord_cartesian(clip='off') +
#         ggtitle(title) +
#         theme_bw() +
#         theme(panel.grid = element_blank(),
#               panel.border = element_rect(colour = "grey", fill=NA, size=1),
#               axis.text = element_blank(),
#               axis.ticks = element_blank(),
#               text = element_text(size=30),
#               plot.title = element_text(hjust=.5, size=30),
#               axis.title.x=element_text(vjust=5),
#               legend.position='bottom'
#               # aspect.ratio = 1
#              ) + 
#         # labs (y = "Semantic Space Y", x = "Semantic Space X") + 
#         xlab('') + ylab('') +
#         theme(legend.key = element_blank()) ;
    
#     one.x_range = max(one.data$plot_X) - min(one.data$plot_X);
#     one.y_range = max(one.data$plot_Y) - min(one.data$plot_Y);
#     p1 <- p1 + xlim(min(one.data$plot_X)-one.x_range/10,max(one.data$plot_X)+one.x_range/10);
#     p1 <- p1 + ylim(min(one.data$plot_Y)-one.y_range/10,max(one.data$plot_Y)+one.y_range/10);
    
#     p1
# }









# plot_enrichments_bubble <- function(enrichments, n=10, facet='v') {
    
#     if (facet=='v') {
#         p <- enrichments %>%
#         ggplot(aes(enrichment, neglogFDR)) + 
#         facet_wrap(~G, scales='free_x')
#     } else {
#         p <- enrichments %>%
#         ggplot(aes(enrichment, neglogFDR)) + 
#         facet_wrap(~G, scales='free')
#     }
    
#     p +
#     # geom_point(aes(size=n_genes), alpha=.8, fill=brewer.blues(5)[2], shape=21) +
#     geom_point(aes(size=n_genes, fill=direction), alpha=.8, shape=21) +
#     geom_text_repel(aes(label=description), data = enrichments %>% filter(rank < n), force=100, force_pull=.1, nudge_y=-.1, size=9) +
#     scale_y_continuous(name='FDR p-value', breaks=c(-log10(0.05),2,3,4,5), labels=function(x) 10^(-x)) +
#     scale_size_continuous(name='# genes', range=c(1,20)) +
#     scale_fill_manual(values=brewer.rdbu(5)[c(2,4)], guide='none') +
#     xlab('Enrichment Ratio') +
#     theme_minimal() + 
#     theme(panel.grid = element_blank(),
#           panel.border = element_rect(fill=NA),
#           # legend.position=c(.9,.9),
#           legend.position='bottom',
#           strip.text.y = element_text(size=36, angle=0)
#          )
    
# }
                       

# plot_enrichment_nulls <- function(true_scores, null_scores, null_p, how='mean') {

#     mycolors = c(brewer.rdylbu(6)[1:3],brewer.rdylbu(5)[4:5])

#     true_scores <- true_scores %>% rownames_to_column('m') %>% gather(pc, score, -m)
#     null_scores <- null_scores %>% gather(pc, score, -m)
#     null_p <- null_p %>% rownames_to_column('m') %>% gather(pc, p, -m) %>% mutate(sig = p<.05)

#     p <- null_scores %>% 
#     left_join(null_p, by = c('m', 'pc')) %>% 
#     ggplot() + 
#     facet_grid(factor(m, levels=unique(null_p$m))~pc, switch='y') +
#     geom_density(aes(score, alpha=sig), color=mycolors[5], fill=mycolors[4]) +
#     scale_alpha_manual(values=c(.2,1), guide='none') +
#     geom_point(data=true_scores, aes(x=score, y=0), color=mycolors[1], size=5)
    
#     if (how == 'mean') {
#         p <- p + 
#         geom_text(data=true_scores, aes(label=paste0('s=',round(score,3)), x=0.005, y=200), size=8, hjust=0) +
#         geom_text(data=null_p, aes(label=paste0('p=',round(p,3)), x=0.005, y=100), size=8, hjust=0)
#     } else {
#         p <- p + 
#         geom_text(data=true_scores, aes(label=paste0('s=',round(score,2)), x=2000, y=0.0005), size=8, hjust=0) +
#         geom_text(data=null_p, aes(label=paste0('p=',round(p,2)), x=2000, y=0.0012), size=8, hjust=0)
#     }
        
#     p +
#     theme_void() + 
#     theme(
#         strip.text.y.left = element_text(size=30, angle=0), 
#         strip.text.x = element_text(size=30, angle=0), 
#         text = element_text(size=30))
# }
                       
                       

# plot_enrichment_bars_2 <- function(null_p, xlab='z-score') {
    
#     lim <- max(abs(null_p$z))*1.2
    
#     null_p %>% 
#     mutate(sig_label = case_when(
#         q < .001 ~ '***',q < .01 ~ '**',q < .05 ~ '*', TRUE ~ ''
#         )) %>%
#     mutate(hjust=ifelse(z < 0, 1, 0)) %>%
#     ggplot(aes(x=z, y=label)) + 
#     facet_grid(.~G) +
#     geom_col(aes(fill=posneg, alpha=sig), width=.8, position=position_dodge(width=.8)) +
#     geom_text(aes(label=sig_label, vjust=.7, hjust=hjust, group=posneg), position=position_dodge(width=.8), size=6) +
#     scale_alpha_manual(values=c(0.3,1), guide='none') +
#     scale_y_discrete(limits=rev, name='') +
#     scale_x_continuous(limits=c(-lim,lim), breaks=round(0.5*c(-lim,0,lim)), name=xlab) +
#     scale_fill_manual(values=rev(brewer.rdbu(10))[c(3,8)], name='Gradient\ndirection') +
#     # scale_fill_gradientn(colors=rev(brewer.rdbu(100)), guide='none', limits=c(-lim,lim)) +
#     # ggtitle("Mean weight of cell genes vs permutations") +
#     # coord_flip() +
#     coord_cartesian(clip='off') +
#     theme_minimal() +
#     theme(panel.grid.minor=element_blank(),
#           plot.title=element_text(size=20),
#           axis.text=element_text(size=20),
#           text=element_text(size=20)
#          )
    
# }