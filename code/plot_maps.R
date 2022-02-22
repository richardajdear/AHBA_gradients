library(ggseg)
library(ggsegGlasser)
# library(ggrepel)

source("../code/brainPlots.R")


mycolors = c(brewer.rdylbu(6)[1:3],brewer.rdylbu(5)[4:5])
# mycolors = saturation(c(brewer.rdylbu(5)[1:2],brewer.rdylbu(6)[4:6]), delta(.2))
# mycolors2 = c(mycolors[1:2], saturation(mycolors[3], delta(.4)))
# mycolors3 = c(mycolors[1], 
#               saturation(mycolors[3], delta(.4)), 
#               saturation(mycolors[4], delta(.4))) 

plot_maps <- function(maps, title="", colors=rev(brewer.rdbu(100)), ncol=3) {
    df <- maps %>%
        rownames_to_column %>%
        mutate(region = recode(rowname,'7Pl'='7PL')) %>% 
        select(-rowname) %>%
        gather('map', 'value', -region) %>%
        mutate_at(vars(map), ~ factor(., levels=unique(.))) %>% 
        group_by(map)
    
    glasser$data <- glasser$data %>% filter(hemi=='left')
    
    # set scale limits at 99th percentile
    m <- pmax(
        df %>% .$value %>% quantile(.99) %>% abs,
        df %>% .$value %>% quantile(.01) %>% abs
    )
    
    ggplot(df) + 
    geom_brain(
        atlas=glasser,
        mapping=aes(fill=value, geometry=geometry, hemi=hemi, side=side, type=type),
        colour='grey', size=.1,
        show.legend=T
        ) + 
    facet_wrap(~map, ncol=ncol, dir="v") +
    scale_fill_gradientn(
        colors=colors, 
        limits=c(-m,m), oob=squish, breaks=c(-m,0,m), 
        labels=c(round(-m,2),0,round(m,2)), 
        name=''
    ) +
    theme_void() + 
    theme(legend.position='bottom',
          panel.spacing.x=unit(4,'lines'),
          strip.text.x=element_text(vjust=1, size=20),
          plot.title=element_text(hjust=0.5)) +
ggtitle(title) + xlab("") + ylab("")
}



plot_null_corrs <- function(corrs, null_corrs, null_p) {

    null_p <- null_p %>% rownames_to_column('map') %>% gather(pc, p, -map) %>% mutate(sig = p<.05)
    null_corrs <- null_corrs %>% gather(pc, corr, -map)
    corrs <- corrs %>% rownames_to_column('map') %>% gather(pc, corr, -map) %>% left_join(null_p, by=c('pc', 'map'))

    null_corrs %>% 
    left_join(null_p, by = c('map', 'pc')) %>%
    ggplot() + 
    facet_grid(factor(map, levels=unique(null_p$map), ordered=T)~pc, switch='y') +
    geom_density(aes(corr), color=mycolors[5], fill=mycolors[4]) +
    # facet_grid(.~pc) +
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



plot_map_scatters <- function(corrs, null_corrs, null_p) {
    corrs <- corrs %>% rownames_to_column('map') %>% gather(pc, corr, -map)
    # null_corrs <- null_corrs %>% gather(pc, corr, -map)
    null_p <- null_p %>% rownames_to_column('map') %>% gather(pc, p, -map) %>% mutate(sig = p<.05)

    null_corrs %>% 
    left_join(null_p, by = c('map', 'pc')) %>% 
    ggplot() + 
    facet_grid(factor(map, levels=unique(null_p$map))~pc, switch='y') +
    # geom_histogram(aes(corr, alpha=sig)) +
    geom_density(aes(corr, alpha=sig), color=mycolors[5], fill=mycolors[4]) +
    scale_alpha_manual(values=c(.2,1), guide='none') +
    ylim(c(-.6, NA)) +
    ylab('') + 
    geom_point(data=corrs, aes(x=corr, y=0), color=mycolors[1], size=5) +
    # geom_vline(data=corrs, aes(xintercept=corr), color='red') +
    geom_text(data=corrs, aes(label=paste0('r=',round(corr,2)), x=.4, y=1.4), size=8, hjust=0) +
    geom_text(data=null_p, aes(label=paste0('p=',round(p,2)), x=.4, y=.8), size=8, hjust=0) +
    # ggtitle('Distribution of map correlations over spin permutations') +
    theme_void() + theme(strip.text = element_blank())
}



plot_maps_combined <- function(maps, corrs, null_corrs, null_p) {
    g1 <- plot_null_corrs(corrs, null_corrs, null_p)
    g2 <- plot_hcp_wide(scores_plot, spacing=0) + guides(fill='none')
    g3 <- plot_maps(maps, colors=rev(brewer.puor(100)), ncol=1) + guides(fill='none')

    (
        ((plot_spacer() | g2) + plot_layout(widths=c(1,2))) / 
        ((g3 | g1) + plot_layout(widths=c(1,2)))
    ) + plot_layout(heights=c(1.2,12))       
}
