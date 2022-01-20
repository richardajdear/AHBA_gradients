library(ggseg)
library(ggsegGlasser)
library(ggrepel)

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


plot_spin_corrs <- function(corrs, spin_corrs, spin_p, spin_sig) {
    corrs <- corrs %>% rownames_to_column('map') %>% gather(pc, corr, -map)
    spin_corrs <- spin_corrs %>% gather(pc, corr, -map)
    spin_p <- spin_p %>% rownames_to_column('map') %>% gather(pc, p, -map)
    spin_sig <- spin_sig %>% rownames_to_column('map') %>% gather(pc, sig, -map)

    spin_corrs %>% 
    left_join(spin_sig, by = c('map', 'pc')) %>% 
    ggplot() + 
    facet_grid(factor(map, levels=unique(spin_sig$map))~pc, switch='y') +
    geom_histogram(aes(corr, alpha=sig)) +
    scale_alpha_manual(values=c(.4,1), guide='none') +
    ylab('') +
    geom_vline(data=corrs, aes(xintercept=corr), color='red') +
    geom_text(data=corrs, aes(label=paste0('r=',round(corr,3)), x=.5, y=120), size=8, hjust=0) +
    geom_text(data=spin_p, aes(label=paste0('p=',round(p,3)), x=.5, y=60), size=8, hjust=0) +
    # ggtitle('Distribution of map correlations over spin permutations') +
    theme_minimal() + 
    theme(
          strip.text.y.left = element_text(angle=0),
          panel.grid = element_blank(),
          axis.text.y = element_blank(),
          text = element_text(size=30)
         )
}