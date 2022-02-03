
plot_hcp_bs_mapping <- function(hcp_bs_mapping) {
    hcp_bs_mapping <- hcp_bs_mapping %>% 
        mutate_at(vars(cortex), ~ factor(., levels=unique(.)))
    cols = as.character(cols25())
    
    
    g0 <- ggplot(hcp_bs_mapping) + 
    geom_brain(atlas=glasser, mapping=aes(fill=cortex, geometry=geometry, hemi=hemi, side=side, type=type)) +
    guides(fill=guide_legend('')) +
    scale_fill_manual(values=cols) +
    theme_void() + ggtitle('HCP cortex') + theme(text=element_text(size=20), plot.title=element_text(vjust=-1))

    g1 <- ggplot(hcp_bs_mapping) + 
    geom_brain(atlas=glasser, mapping=aes(fill=structure_name, geometry=geometry, hemi=hemi, side=side, type=type)) +
    guides(fill=guide_legend('')) +
    scale_fill_manual(values=cols) +
    theme_void() + ggtitle('Brainspan regions mapped to HCP cortex') + theme(text=element_text(size=20), plot.title=element_text(vjust=-1))

    g0/g1
}



plot_bs_pcs_corr <- function(bs_pcs_corr, title="", xint='4 mos') {
    bs_pcs_corr %>% 
    ggplot() + geom_hline(yintercept=0, color='grey') + geom_vline(xintercept=xint, color='grey') +
    geom_line(aes(x=age, y=corr, color=PC, group=PC)) + 
    ggtitle(title) +
    theme_minimal() + theme(text=element_text(size=20), axis.text.x=element_text(size=16, angle=30, hjust=1))
}