library(ggseg)
library(ggsegGlasser)
library(ggrepel)

plot_maps <- function(maps, title="", colors=rev(brewer.rdbu(100))) {
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
    facet_wrap(~map, ncol=3,dir="v") +
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