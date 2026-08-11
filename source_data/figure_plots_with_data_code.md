The misplaced mouse Pax6 interneuron subclass: A cross-species
transcriptomic reassignment
================
Jonathan Werner
2026-08-11

This markdown file contains the code for generating all plots for the
publication: Werner, Suresh, French, and Gillis. The misplaced mouse
Pax6 interneuron subclass: A cross-species transcriptomic reassignment.
2026.

All data for these plots is provided at:
<https://github.com/JonathanMWerner/mammalian_brain_consensus_taxonomy/source_data/data_for_plots>

``` r
library(ggplot2)
library(ggrepel)
library(dplyr)
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
library(MetaMarkers)
celltype_color_palette = c("Chandelier"='#F641A8', 
                           "Pvalb"='#D93137', 
                           "Sst"='#FF9900', 
                           "Sst Chodl"='#FFD700', 
                           "Lamp5"='#DA808C',
                           "Sncg"='#DF70FF',
                           "Vip"='#A45FBF', 
                           "Lamp5 Lhx6"='#935F50',
                           "Pax6"='#71238C',
                           'NA' = 'grey90')
```

## Figure 1C

``` r
load(file = 'data_for_plots/namesake_marker_de.Rdata')
load(file = 'data_for_plots/primate_subclass_meta_markers.Rdata')

primate_pax6_PAX6_auroc = primate_subclass_meta_markers %>% filter(cell_type == 'Pax6' & gene == 'PAX6') %>% pull(auroc)

ggplot(mouse_name_sake_genes, aes(x = primate_auroc, y = mouse_auroc, color = cell_type, label = plot_label)) +
  geom_point(size = 3, alpha = .75) + ylim(0,1) + xlim(0,1) +
  geom_label_repel(min.segment.length = 0, max.overlaps = Inf, hjust = 1, direction = 'y', nudge_x = -.4,box.padding = .1, show.legend = F) +
  geom_vline(xintercept = .5, color = 'black', linetype = 'dashed', alpha = .5) +
  geom_hline(yintercept = .5, color = 'black', linetype = 'dashed', alpha = .5) +
  scale_color_manual(values = celltype_color_palette, name = 'Subclass', 
                     breaks = c('Sst','Sst Chodl','Pvalb','Chandelier','Lamp5 Lhx6','Lamp5','Vip','Sncg')) +
  ylab('Mouse AUROC') + xlab('Primate AUROC') +
  theme_bw()+ ggtitle('Subclass namesake gene DE') +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  facet_wrap(~gene, nrow = 3) +
  geom_vline(data = subset(mouse_name_sake_genes, gene == "PAX6"), aes(xintercept = primate_pax6_PAX6_auroc), 
             color = 'red', linetype = 'longdash')
```

    ## Warning: Removed 63 rows containing missing values or values outside the scale range
    ## (`geom_label_repel()`).

![](figure_plots_with_data_code_files/figure-gfm/namesake_marker_de-1.png)<!-- -->

## Figure 1D

## Figure 1E

``` r
all_stats_df = readRDS('data_for_plots/all_stats_df.rds')

num_total_markers = c(1, 10, 25, 50, 100)
marker_colors = c('white','grey70','grey50','grey30','black')
names(marker_colors) = as.character(num_total_markers)

ggplot(all_stats_df, aes( x = avg_agg_exp, y =avg_enrich_exp , color =color_label, group = `Mouse Subclass`)) + 
  geom_point(shape = 21,size = 3, aes(fill = num_markers)) +
  geom_line(linewidth = 1) + geom_hline(yintercept = 1, color = 'black', linetype = 'dashed') +
  scale_color_manual(values = celltype_color_palette, name = 'Mouse Subclass') +
  scale_fill_manual(values = marker_colors, name = '# markers') +
  theme_bw() + ylab('Mouse expression enrichment') + xlab('Aggregate mouse expression') +
  facet_wrap(~`Primate Subclass Markers`, scales = 'free') + ggtitle('Primate marker expression in Mouse subclasses')
```

![](figure_plots_with_data_code_files/figure-gfm/mouse_agg_marker-1.png)<!-- -->

## Figure 1F

``` r
snapshot_df = all_stats_df %>% filter(`Primate Subclass Markers` %in% c('Primate Sncg','Primate Pax6') & `Mouse Subclass` == 'Sncg') %>%
  select(`Primate Subclass Markers`, `Mouse Subclass`, avg_enrich_exp, num_markers) %>% 
  tidyr::pivot_wider(names_from = `Primate Subclass Markers`, values_from = avg_enrich_exp)


ggplot(snapshot_df , aes(x = `Primate Sncg`, y = `Primate Pax6`, color = `Mouse Subclass`)) +
  geom_point(shape = 21,size = 5, aes(fill = num_markers)) +
  geom_abline(slope = 1, intercept = 0, color = 'red') +
  ylim(0, 2.8) + xlim(0, 2) +
  scale_color_manual(values = celltype_color_palette, name = 'Mouse Subclass') +
  scale_fill_manual(values = marker_colors, name = '# markers') +
  ggtitle('Mouse Sncg Subclass') + ylab('Primate Pax6 marker enrichment') + xlab('Primate Sncg marker enrichment') +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())
```

![](figure_plots_with_data_code_files/figure-gfm/mouse_primate_sncg_and_pax6_enrich-1.png)<!-- -->

## Figure 2H

``` r
load(file = 'data_for_plots/mouse_meta_markers.Rdata')

plot_pareto_summary(mouse_subclass_meta_markers %>% filter(cell_type %in% c('Pax6','Sncg')), min_recurrence=4) +
  theme_bw() + scale_color_manual(values = celltype_color_palette) + ggtitle('Homologous Mouse subclass markers') + ylim(.5, 1) + xlim(0,5)
```

![](figure_plots_with_data_code_files/figure-gfm/mouse_pareto-1.png)<!-- -->

## Figure 2I

``` r
load(file = 'data_for_plots/primate_subclass_meta_markers.Rdata')

human_de_counts = rowSums(primate_subclass_meta_markers  %>% select(human_10x, human_ss))
chimp_de_counts = rowSums(primate_subclass_meta_markers  %>% select(chimp_10x, chimp_ss))
gorilla_de_counts = rowSums(primate_subclass_meta_markers  %>% select(gorilla_10x, gorilla_ss))
macaca_de_counts = primate_subclass_meta_markers$macaca_10x
marmoset_de_counts = primate_subclass_meta_markers$marmoset_10x

#Getting genes DE in at least 1 dataset per species
cross_species_de_index = human_de_counts >= 1 & chimp_de_counts >= 1 & gorilla_de_counts >= 1 & macaca_de_counts >= 1 & marmoset_de_counts >= 1


plot_pareto_summary(primate_subclass_meta_markers[cross_species_de_index, ] %>% filter(cell_type %in% c('Pax6','Sncg')), min_recurrence=0) +
  theme_bw() + scale_color_manual(values = celltype_color_palette) + ggtitle('Homologous Primate subclass markers') + ylim(.5, 1) + xlim(0,4)
```

![](figure_plots_with_data_code_files/figure-gfm/primate_pareto-1.png)<!-- -->

## Figure 2J

``` r
load(file = 'data_for_plots/mammal_subclass_meta_markers.Rdata')

mouse_de_counts = rowSums(mammal_subclass_meta_markers %>% select(tasic18, zeng_10x_cells, zeng_10x_nuclei, zeng_smart_cells, zeng_smart_nuclei, zeng_10xv2_wba, zeng_10xv3_wba, macosko_wba))

human_de_counts = rowSums(mammal_subclass_meta_markers %>% select(human_10x, human_ss))
chimp_de_counts = rowSums(mammal_subclass_meta_markers %>% select(chimp_10x, chimp_ss))
gorilla_de_counts = rowSums(mammal_subclass_meta_markers %>% select(gorilla_10x, gorilla_ss))
macaca_de_counts = mammal_subclass_meta_markers$macaca_10x
marmoset_de_counts = mammal_subclass_meta_markers$marmoset_10x


#Getting genes DE in at least 1 dataset per species
cross_species_de_index = mouse_de_counts >= 1 & human_de_counts >= 1 & chimp_de_counts >= 1 & gorilla_de_counts >= 1 & macaca_de_counts >= 1 & marmoset_de_counts >= 1

plot_pareto_summary(mammal_subclass_meta_markers[cross_species_de_index, ] %>% filter(cell_type %in% c('Pax6','Sncg')), min_recurrence=6) +
  theme_bw() + scale_color_manual(values = celltype_color_palette) + ggtitle('Homologous Primate-Mouse subclass markers') + ylim(.5, 1) + xlim(0,4)
```

![](figure_plots_with_data_code_files/figure-gfm/mammal_pareto-1.png)<!-- -->

## Figure 6A mouse panels

``` r
load(file = 'data_for_plots/mouse_Vip_dot_plot_primate_enrich.Rdata')

order_vec = plotting_df %>% arrange(match(mapped_primate_subclass, c('Sst','Sst Chodl','Pvalb','Chandelier','Lamp5 Lhx6','Lamp5','Vip','Pax6','Sncg'))) %>%
  pull(study_primate_subclass)

plotting_df$study_primate_subclass = factor(plotting_df$study_primate_subclass, levels = order_vec)

mouse_vip_enrich_dot_plot = ggplot(plotting_df, aes(x = study_primate_subclass, y = mean_enrichment, color = mapped_primate_subclass)) +
  geom_point(size = 3, show.legend = F) + 
  geom_linerange(aes(ymin = mean_enrichment - sd_enrichment, ymax = mean_enrichment + sd_enrichment),
                 show.legend = F) +
  geom_hline(yintercept = 1, linetype = 'dashed', color = 'black') +
  ylab('Mean Primate Vip enrichment') + xlab('Mouse dataset Primate subclasses') + ylim(0,3) +
  scale_x_discrete(guide = guide_axis(n.dodge=2)) +
  scale_color_manual(values = celltype_color_palette) +
  theme_bw() + theme(panel.border = element_blank(), axis.line = element_line(),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),axis.text.x=element_blank())

mouse_vip_enrich_dot_plot
```

![](figure_plots_with_data_code_files/figure-gfm/mouse_primate_enrich_dot_plots-1.png)<!-- -->

``` r
load(file = 'data_for_plots/mouse_Pax6_dot_plot_primate_enrich.Rdata')
order_vec = plotting_df %>% arrange(match(mapped_primate_subclass, c('Sst','Sst Chodl','Pvalb','Chandelier','Lamp5 Lhx6','Lamp5','Vip','Pax6','Sncg'))) %>%
  pull(study_primate_subclass)

plotting_df$study_primate_subclass = factor(plotting_df$study_primate_subclass, levels = order_vec)

mouse_pax6_enrich_dot_plot = ggplot(plotting_df, aes(x = study_primate_subclass, y = mean_enrichment, color = mapped_primate_subclass)) +
  geom_point(size = 3, show.legend = F) + 
  geom_linerange(aes(ymin = mean_enrichment - sd_enrichment, ymax = mean_enrichment + sd_enrichment),
                 show.legend = F) +
  geom_hline(yintercept = 1, linetype = 'dashed', color = 'black') +
  ylab('Mean Primate Pax6 enrichment') + xlab('Mouse dataset Primate subclasses') + ylim(0,3) +
  scale_x_discrete(guide = guide_axis(n.dodge=2)) +
  scale_color_manual(values = celltype_color_palette) +
  theme_bw() + theme(panel.border = element_blank(), axis.line = element_line(),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),axis.text.x=element_blank())

mouse_pax6_enrich_dot_plot
```

![](figure_plots_with_data_code_files/figure-gfm/mouse_primate_enrich_dot_plots-2.png)<!-- -->

``` r
load(file = 'data_for_plots/mouse_Sncg_dot_plot_primate_enrich.Rdata')

order_vec = plotting_df %>% arrange(match(mapped_primate_subclass, c('Sst','Sst Chodl','Pvalb','Chandelier','Lamp5 Lhx6','Lamp5','Vip','Pax6','Sncg'))) %>%
  pull(study_primate_subclass)

plotting_df$study_primate_subclass = factor(plotting_df$study_primate_subclass, levels = order_vec)

mouse_sncg_enrich_dot_plot = ggplot(plotting_df, aes(x = study_primate_subclass, y = mean_enrichment, color = mapped_primate_subclass)) +
  geom_point(size = 3, show.legend = F) + 
  geom_linerange(aes(ymin = mean_enrichment - sd_enrichment, ymax = mean_enrichment + sd_enrichment),
                 show.legend = F) +
  geom_hline(yintercept = 1, linetype = 'dashed', color = 'black') +
  ylab('Mean Primate Sncg enrichment') + xlab('Mouse dataset Primate subclasses') + ylim(0,3) +
  scale_x_discrete(guide = guide_axis(n.dodge=2)) +
  scale_color_manual(values = celltype_color_palette) +
  theme_bw() + theme(panel.border = element_blank(), axis.line = element_line(),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),axis.text.x=element_blank())

mouse_sncg_enrich_dot_plot
```

![](figure_plots_with_data_code_files/figure-gfm/mouse_primate_enrich_dot_plots-3.png)<!-- -->

``` r
load(file = 'data_for_plots/mouse_Sncg_dot_plot_mouse_enrich.Rdata')

order_vec = plotting_df %>% arrange(match(mapped_primate_subclass, c('Sst','Sst Chodl','Pvalb','Chandelier','Lamp5 Lhx6','Lamp5','Vip','Pax6','Sncg'))) %>%
  pull(study_primate_subclass)

plotting_df$study_primate_subclass = factor(plotting_df$study_primate_subclass, levels = order_vec)

mouse_mouse_sncg_enrich_dot_plot = ggplot(plotting_df, aes(x = study_primate_subclass, y = mean_enrichment, color = mapped_primate_subclass)) +
  geom_point(size = 3, show.legend = F) + 
  geom_linerange(aes(ymin = mean_enrichment - sd_enrichment, ymax = mean_enrichment + sd_enrichment),
                 show.legend = F) +
  geom_hline(yintercept = 1, linetype = 'dashed', color = 'black') +
  ylab('Mean Mouse Sncg enrichment') + xlab('Mouse dataset Primate subclasses') + ylim(0,3.7) +
  scale_x_discrete(guide = guide_axis(n.dodge=2)) +
  scale_color_manual(values = celltype_color_palette) +
  theme_bw() + theme(panel.border = element_blank(), axis.line = element_line(),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),axis.text.x=element_blank())

mouse_mouse_sncg_enrich_dot_plot
```

![](figure_plots_with_data_code_files/figure-gfm/mouse_primate_enrich_dot_plots-4.png)<!-- -->
