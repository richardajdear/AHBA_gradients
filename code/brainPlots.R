# R functions to plot brains using ggseg

library(ggseg)
library(ggsegGlasser)
suppressMessages(library(tidyverse))
suppressMessages(library(scales))
library(patchwork)
library(pals)

plot_dk <- function(scores_df, title="", three=F, switch=NULL) {
    df <- scores_df %>% 
        rename('PC1'='0', 'PC2'='1', 'PC3'='2', 'PC4'='3', 'PC5'='4') %>%
        mutate_at(vars(version), ~ factor(., levels=unique(.))) %>% 
        gather('component', 'score', -version, -label) %>%
        group_by(component, version)
    
    if (three) {
        df <- df %>% filter(component %in% c('PC1','PC2','PC3'))
    }
    
    m <- pmax(
#         df %>% filter(component=='PC1') %>% .$score %>% quantile(.95) %>% abs,
#         df %>% filter(component=='PC1') %>% .$score %>% quantile(.05) %>% abs
        df %>% .$score %>% quantile(.99) %>% abs,
        df %>% .$score %>% quantile(.01) %>% abs
    )

    dk$data <- dk$data %>% filter(hemi=='left')
    
    ggplot(df) + 
    geom_brain(
        atlas=dk,
        mapping=aes(fill=score, geometry=geometry, hemi=hemi, side=side, type=type),
        colour='grey', size=.1,
        show.legend=T
        ) + 
    theme_void() +
    facet_grid(component~version, switch=switch) +
    theme(legend.position='bottom', 
          strip.text.x=element_text(vjust=1),
          plot.title=element_text(hjust=0.5)) +
    #   scale_fill_cmocean(name='balance', limits=c(-m,m), oob=squish) +
#     scale_fill_gradient2(low=muted('red'), high=muted('blue'), 
    scale_fill_gradientn(colors=rev(brewer.rdbu(100)), 
                         limits=c(-m,m), oob=squish, breaks=c(-m,0,m), 
                         labels=c(round(-m,2),0,round(m,2)), name=''
                        ) +
ggtitle(title) + xlab("") + ylab("")
}



plot_hcp <- function(scores_df, title="", three=F, switch=NULL) {
    df <- scores_df %>% 
        rename('PC1'='0', 'PC2'='1', 'PC3'='2', 'PC4'='3', 'PC5'='4') %>%
        mutate_at(vars(version), ~ factor(., levels=unique(.))) %>% 
        mutate(region = recode(label,'7Pl'='7PL')) %>% select(-label) %>%
        gather('component', 'score', -version, -region) %>%
        group_by(component, version)

    if(three) {
        df <- df %>% filter(component %in% c('PC1','PC2','PC3'))
    }
    
    m <- pmax(
#         df %>% filter(component=='PC1') %>% .$score %>% quantile(.95) %>% abs,
#         df %>% filter(component=='PC1') %>% .$score %>% quantile(.05) %>% abs
        df %>% .$score %>% quantile(.99) %>% abs,
        df %>% .$score %>% quantile(.01) %>% abs
    )

    glasser$data <- glasser$data %>% filter(hemi=='left')
    
    ggplot(df) + 
    geom_brain(
        atlas=glasser,
        mapping=aes(fill=score, geometry=geometry, hemi=hemi, side=side, type=type),
        colour='grey', size=.1,
        show.legend=T
        ) + 
    theme_void() + 
    facet_grid(component~version, switch=switch) +
    theme(legend.position='bottom',
          strip.text.x=element_text(vjust=1),
          strip.text.y.left = element_text(angle = 0),
          plot.title=element_text(hjust=0.5)) +
    #   scale_fill_cmocean(name='balance', limits=c(-m,m), oob=squish) +
#     scale_fill_gradient2(low=muted('red'), high=muted('blue'), 
    scale_fill_gradientn(colors=rev(brewer.rdbu(100)), 
                         limits=c(-m,m), oob=squish, breaks=c(-m,0,m), 
                         labels=c(round(-m,2),0,round(m,2)), name=''
                        ) +
    coord_sf(clip='off') +
ggtitle(title) + xlab("") + ylab("")
}




plot_hcp_wide <- function(scores_df, title="", facet='h', spacing=4) {
    df <- scores_df %>% 
        rename('PC1'='0', 'PC2'='1', 'PC3'='2') %>%
        mutate(region = recode(label,'7Pl'='7PL')) %>% select(-label) %>%
        gather('component', 'score', -region) %>%
        group_by(component)

    m <- pmax(
#         df %>% filter(component=='PC1') %>% .$score %>% quantile(.95) %>% abs,
#         df %>% filter(component=='PC1') %>% .$score %>% quantile(.05) %>% abs
        df %>% .$score %>% quantile(.99) %>% abs,
        df %>% .$score %>% quantile(.01) %>% abs
    )

    glasser$data <- glasser$data %>% filter(hemi=='left')
    
    p <- ggplot(df) + 
    geom_brain(
        atlas=glasser,
        mapping=aes(fill=score, geometry=geometry, hemi=hemi, side=side, type=type),
        colour='grey', size=.1,
        show.legend=T
        ) + 
    theme_void() + 
    theme(strip.text=element_blank()) +
    theme(panel.spacing=unit(spacing,'lines')) +
    theme(legend.position='right', text=element_text(size=20),
          plot.title=element_text(hjust=0.5)) +
    #   scale_fill_cmocean(name='balance', limits=c(-m,m), oob=squish) +
    scale_fill_gradientn(colors=rev(brewer.rdbu(100)), 
                         limits=c(-m,m), oob=squish, breaks=c(-m,0,m), 
                         labels=c(round(-m,2),0,round(m,2)), name=''
                        ) +
ggtitle(title) + xlab("") + ylab("")
    
    if (facet=='h') {
        p + facet_grid(.~component)
    } else {
        p + facet_grid(component~.)   
    }
}




plot_dk_dist <- function(dist_dk, title="",
                         name='Samples', lim=c(0,100)) {
    dk$data <- dk$data %>% filter(hemi=='left')
    
    ggplot(dist_dk) + 
    geom_brain(
        atlas=dk,
        mapping=aes(fill=count, geometry=geometry),
        colour='grey', size=.1,
        show.legend=T
        ) + 
    theme_void() +
    theme(plot.title=element_text(hjust=0.5)) +
#     facet_grid(component~version) +
#     theme(legend.position='bottom', text=element_text(size=20), 
#           plot.title=element_text(hjust=0.5)) +
#     guides(fill=guide_colorbar(barwidth=30, title.hjust=-.2, title.vjust=1)) +
    #   scale_fill_cmocean(name='balance', limits=c(-m,m), oob=squish) +
    scale_fill_gradientn(colors=brewer.blues(100), 
                         limits=lim, oob=squish, name=name
                        ) +
    guides(fill=guide_colorbar(barheight=10)) +
#     scale_fill_viridis(option=palette, limits=lim, name=name, direction=-1) +
    ggtitle(title) + xlab("") + ylab("")
}

plot_hcp_dist <- function(dist_hcp, title="",
                          name='Samples', lim=c(0,100)) {
    glasser$data <- glasser$data %>% filter(hemi=='left')
    
    ggplot(dist_hcp %>% mutate(region = recode(label,'7Pl'='7PL')) %>% select(-label)) + 
    geom_brain(
        atlas=glasser,
        mapping=aes(fill=count, geometry=geometry),
        colour='grey', size=.1,
        show.legend=T
        ) + 
    theme_void() +
    theme(plot.title=element_text(hjust=0.5)) +
#     facet_grid(component~version) +
#     theme(legend.position='bottom', text=element_text(size=20), 
#           plot.title=element_text(hjust=0.5)) +
#     guides(fill=guide_colorbar(barwidth=30, title.hjust=-.2, title.vjust=1)) +
    #   scale_fill_cmocean(name='balance', limits=c(-m,m), oob=squish) +
    scale_fill_gradientn(colors=brewer.blues(100), 
                         limits=lim, oob=squish, name=name
                        ) +
    guides(fill=guide_colorbar(barheight=10)) +
#     scale_fill_viridis(option=palette, limits=lim, name=name, direction=-1) +
    ggtitle(title) + xlab("") + ylab("")
}



plot_hcp_classes <- function(df, classes=vonEconomo, classcolors=vonEconomo.colors) {
    colors <- df %>% select( {{classes}} ,  {{classcolors}} ) %>% unique() %>% 
        arrange( {{classes}} ) %>% drop_na() %>% pull( {{classcolors}} )
    names <- c('motor', 'assocation', 'association', 'sensory', 'sensory', 'limbic', 'insula')
    glasser$data <- glasser$data %>% filter(hemi=='left')

    df %>% 
    mutate(region = recode(label,'7Pl'='7PL')) %>% select(-label) %>%
    mutate( classes = factor( {{classes}} )) %>%
    filter( classes  != 'NaN') %>%
    ggplot() +
    geom_brain(
        atlas=glasser,
        mapping=aes(fill= classes , geometry=geometry, hemi=hemi, side=side, type=type),
        colour='grey', size=.1
    ) +
    theme_void() +
    scale_fill_manual(values=colors, guide="none")
}




plot_hcp_ranks <- function(scores_df, title="", three=F, switch=NULL) {
    df <- scores_df %>% 
        rename('PC1'='0', 'PC2'='1', 'PC3'='2', 'PC4'='3', 'PC5'='4') %>%
        mutate_at(vars(version), ~ factor(., levels=unique(.))) %>% 
        mutate(region = recode(label,'7Pl'='7PL')) %>% select(-label) %>%
        drop_na() %>%
        gather('component', 'rank', -version, -region) %>%
        group_by(component, version)

    if(three) {
        df <- df %>% filter(component %in% c('PC1','PC2','PC3'))
    }
    
#     m <- pmax(
# #         df %>% filter(component=='PC1') %>% .$score %>% quantile(.95) %>% abs,
# #         df %>% filter(component=='PC1') %>% .$score %>% quantile(.05) %>% abs
#         df %>% .$score %>% quantile(.99) %>% abs,
#         df %>% .$score %>% quantile(.01) %>% abs
#     )

    glasser$data <- glasser$data %>% filter(hemi=='left')
    
    ggplot(df) + 
    geom_brain(
        atlas=glasser,
        mapping=aes(fill=rank, geometry=geometry),
        colour='grey', size=.1,
        show.legend=T
        ) + 
    theme_void() + 
    facet_grid(component~version, switch=switch) +
    theme(legend.position='bottom',
          strip.text.x=element_text(vjust=1),
          plot.title=element_text(hjust=0.5)) +
ggtitle(title) + xlab("") + ylab("")
}

