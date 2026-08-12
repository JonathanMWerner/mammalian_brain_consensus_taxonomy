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

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

    ## Loading required package: grid

    ## ========================================
    ## ComplexHeatmap version 2.24.1
    ## Bioconductor page: http://bioconductor.org/packages/ComplexHeatmap/
    ## Github page: https://github.com/jokergoo/ComplexHeatmap
    ## Documentation: http://jokergoo.github.io/ComplexHeatmap-reference
    ## 
    ## If you use it in published research, please cite either one:
    ## - Gu, Z. Complex Heatmap Visualization. iMeta 2022.
    ## - Gu, Z. Complex heatmaps reveal patterns and correlations in multidimensional 
    ##     genomic data. Bioinformatics 2016.
    ## 
    ## 
    ## The new InteractiveComplexHeatmap package can directly export static 
    ## complex heatmaps into an interactive Shiny app with zero effort. Have a try!
    ## 
    ## This message can be suppressed by:
    ##   suppressPackageStartupMessages(library(ComplexHeatmap))
    ## ========================================

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

## Figure 2A

``` r
load(file = 'data_for_plots/avg_mouse_clust_primate_subclass_mappings_gaba.Rdata')
load(file = 'data_for_plots/cross_mammal_MN_results.Rdata')


keep_mclusters = max_avg_mapping %>% filter(orthology == 'high_orthology') %>% pull(mouse_mcluster)

r_names = rownames(mouse_clust_primate_sub_best_hits)
r_index = grepl('human', r_names) | grepl('chimp', r_names) | grepl('gorilla', r_names) | grepl('macaca', r_names) | grepl('marmoset', r_names) 
c_index = colnames(mouse_clust_primate_sub_best_hits) %in% keep_mclusters

test_metasubclass_hits = mouse_clust_primate_sub_best_hits[r_index, c_index]
m = test_metasubclass_hits 
#Cluster first, need to swap NAs to 0 for clustering
m[is.na(m)] <- 0
row_dendrogram <- stats::as.dendrogram(stats::hclust(dist(1-m), method = "average"))
col_dendrogram <- stats::as.dendrogram(stats::hclust(dist(1-t(m)), method = "average"))


#And the primate subclass colors and colors for each species
primate_subclass = rownames(test_metasubclass_hits)
primate_subclass = sapply(strsplit(primate_subclass, '|', fixed = T), '[[', 2)
primate_subclass = gsub('_',' ', primate_subclass)

primate_species = sapply(strsplit(rownames(test_metasubclass_hits), '|', fixed = T), '[[', 1)
primate_species = sapply(strsplit(primate_species, '_', fixed = T), '[[', 1)


row_ha = rowAnnotation(primate_subclass = primate_subclass,
                       primate_species = primate_species,
                       col = list(primate_subclass = celltype_color_palette,
                                  primate_species = c('human' = 'black', 'chimp' = 'grey25', 'gorilla' = 'grey40', 'macaca' = 'grey70', 'marmoset' = 'grey90')),
                       show_legend = F)

#Custom order for species
row_df = data.frame(species = primate_species, subclass = primate_subclass, row_index = 1:length(primate_subclass))
row_df = row_df %>% group_by(subclass) %>% arrange(match(species, c('human','chimp','gorilla','macaca','marmoset')), .by_group = T) %>%
  arrange(match(subclass, c("Pvalb", "Sst", 'Sst Chodl','Vip','Pax6', 'Sncg', 'Chandelier', 'Lamp5', 'Lamp5 Lhx6')))
row_custom_order = row_df$row_index



index = match(colnames(test_metasubclass_hits),  max_avg_mapping$mouse_mcluster )

col_ha = HeatmapAnnotation(orthologous_subclass = max_avg_mapping$primate_subclass[index],
                           col = list(orthologous_subclass = celltype_color_palette),
                           annotation_legend_param = list(orthologous_subclass = list(
                             at = c("Sst", "Sst Chodl", "Vip", 'Sncg', 'Pax6', 'Pvalb', 'Chandelier', 'Lamp5 Lhx6', 'Lamp5', 'none') )))

auroc_cols <- rev(grDevices::colorRampPalette(RColorBrewer::brewer.pal(11,"RdYlBu"))(100))

mouse_match_prim_sub = reshape2::melt(test_metasubclass_hits)
mouse_match_prim_sub$primate_subclass = sapply(strsplit(as.character(mouse_match_prim_sub$Var1), fixed = T, split = '|'), '[[', 2)

index = match(mouse_match_prim_sub$Var2, max_avg_mapping$mouse_mcluster)

mouse_match_prim_sub$orthologous_subclass = max_avg_mapping$primate_subclass[index]

mouse_match_prim_sub = mouse_match_prim_sub %>% group_by(Var2, Var1) %>% 
  summarise(mean_auroc = mean(value, na.rm = T),orthologous_subclass = orthologous_subclass) %>%
  filter(is.finite(mean_auroc)) %>%
  filter(mean_auroc == max(mean_auroc)) %>%
  filter(!duplicated(as.character(Var2))) %>%
  arrange(match(orthologous_subclass, 
                c("Pvalb", "Sst", 'Sst Chodl','Vip','Pax6', 'Sncg', 'Chandelier', 'Lamp5', 'Lamp5 Lhx6')), 
          desc(mean_auroc))
```

    ## `summarise()` has regrouped the output.
    ## ℹ Summaries were computed grouped by Var2 and Var1.
    ## ℹ Output is grouped by Var2.
    ## ℹ Use `summarise(.groups = "drop_last")` to silence this message.
    ## ℹ Use `summarise(.by = c(Var2, Var1))` for per-operation grouping
    ##   (`?dplyr::dplyr_by`) instead.

``` r
custom_order = match(as.character(mouse_match_prim_sub$Var2), colnames(test_metasubclass_hits))

ht_mc = ComplexHeatmap::Heatmap(test_metasubclass_hits, name = 'AUROC', col = auroc_cols,
                        row_order = row_custom_order, cluster_rows = F,
                        column_order = custom_order, cluster_columns = F,
                        column_names_gp = gpar(fontsize = 4),
                        row_names_gp = gpar(fontsize = 8),
                        left_annotation = row_ha, top_annotation = col_ha)

ht_mc = draw(ht_mc)
```

    ## Following `at` are removed: none, because no color was defined for
    ## them.
    ## Following `at` are removed: none, because no color was defined for
    ## them.
    ## Following `at` are removed: none, because no color was defined for
    ## them.

![](figure_plots_with_data_code_files/figure-gfm/mouse_primate_MN_heatmap-1.png)<!-- -->

## Figure 2 B, D, and E

The Seurat object is quite large to save as source data

\##Figure 2C

``` r
load(file = 'data_for_plots/subclass_sankey.Rdata')

ggplot(sankey_interneuron_df, aes(x = x, 
               next_x = next_x, 
               node = node, 
               next_node = next_node,
               fill = factor(node),
               label = node)) +
  geom_sankey(flow.alpha = 0.5, node.color = 1) +
  geom_sankey_label(size = 3, color = 1 , show.legend = F) +
  scale_fill_manual(values = celltype_color_palette) +
  theme_sankey(base_size = 16) +
  guides(fill = guide_legend(title = "Subclass"))
```

    ## Warning: The `size` argument of `element_rect()` is deprecated as of ggplot2 3.4.0.
    ## ℹ Please use the `linewidth` argument instead.
    ## ℹ The deprecated feature was likely used in the ggsankey package.
    ##   Please report the issue at <https://github.com/davidsjoberg/ggsankey/issues>.
    ## This warning is displayed once per session.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

![](figure_plots_with_data_code_files/figure-gfm/subclass_sankey-1.png)<!-- -->

## Figure 2F

``` r
load(file = 'data_for_plots/integrated_mouse_taxonomy.Rdata')

hclust_avg <- hclust(as.dist(1-centroid_corr), method = 'ward.D2')

clust_index = match(colnames(centroid_corr), dataset_subclasses$study_clust)
original_subclass = dataset_subclasses$cell_subclass[clust_index]
orthologous_subclass = dataset_subclasses$primate_subclass[clust_index]
study_label = dataset_subclasses$study[clust_index]

study_colors = MetBrewer::met.brewer("Veronese", n=8, type="continuous")
study_colors = study_colors[1:8] 
names(study_colors) = c('tasic18','zeng_10x_cells','zeng_10x_nuclei','zeng_smart_cells',
                        'zeng_smart_nuclei','zeng_10xv2_wba','zeng_10xv3_wba','macosko_wba')



column_ha = HeatmapAnnotation(original_subclass = original_subclass, 
                              orthologous_subclass = orthologous_subclass,
                              study = study_label,
                              col = list(original_subclass = celltype_color_palette, 
                                         orthologous_subclass = celltype_color_palette,
                                         study = study_colors ),
                              show_legend = c("orthologous_subclass" = FALSE))


corr_cols <- rev(grDevices::colorRampPalette(RColorBrewer::brewer.pal(11,"RdYlBu"))(100))
cent_hm = ComplexHeatmap::Heatmap(centroid_corr, col = corr_cols, name = 'spearman' ,
                                  cluster_columns = hclust_avg,
                                  cluster_rows = hclust_avg,
                                  top_annotation = column_ha, column_title = 'integrated_interneurons',
                                  show_row_dend = F, show_row_names = F, show_column_names = F)
cent_hm = ComplexHeatmap::draw(cent_hm)
```

![](figure_plots_with_data_code_files/figure-gfm/integrated_mouse_taxonomy-1.png)<!-- -->

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

## Figure 3A

Seurat object is quite large to save as source data

## Figure 3B

``` r
load( file = 'data_for_plots/main_ion_channel_heatmap.Rdata')

index = match(colnames(mean_exp_cpm_df) , clust_to_subclass_df$cell_cluster_study)

column_annot = HeatmapAnnotation(mouse_subclass = clust_to_subclass_df$cell_subclass[index],
                                 homologous_subclass = clust_to_subclass_df$primate_subclass[index],
                                 col = list(mouse_subclass = celltype_color_palette,
                                            homologous_subclass = celltype_color_palette))

set.seed(123)
h2 = Heatmap(mean_exp_cpm_df, name = 'z-scored Mean CPM',
        show_row_names = T, show_row_dend = F, show_column_names = F,
        clustering_method_rows = 'ward.D2', clustering_method_columns = 'ward.D2',
        top_annotation = column_annot, row_names_gp = grid::gpar(fontsize = 8),
        column_title = 'Ion channel GO term genes')
h2 = draw(h2)
```

![](figure_plots_with_data_code_files/figure-gfm/ion_channel_heatmap-1.png)<!-- -->

## Figure 3C and D

``` r
load(file = 'data_for_plots/omega_ephys_scatter.Rdata')

par(pty= 's')
plot(res_mouse$omega2, res_primate$omega2, main = 'Electrophysiology features', 
     ylab = 'Variance explained by homologous subclasses', 
     xlab = 'Variance explained by original subclasses', xlim = c(0,1), ylim = c(0,1) )
abline(a = 0, b = 1, col = 'red', lty = 'dashed')
```

![](figure_plots_with_data_code_files/figure-gfm/mouse_ephys_omega-1.png)<!-- -->

``` r
load(file = 'data_for_plots/omega_expression_scatter.Rdata')

par(pty= 's')
plot(res_exp_mouse$omega2, res_exp_primate$omega2, main = 'Transcriptomic features', 
     ylab = 'Variance explained by homologous subclasses', 
     xlab = 'Variance explained by original subclasses', xlim = c(0,1), ylim = c(0,1) )
abline(a = 0, b = 1, col = 'red', lty = 'dashed')
```

![](figure_plots_with_data_code_files/figure-gfm/mouse_ephys_omega-2.png)<!-- -->

## Figure 3E

``` r
load(file = 'data_for_plots/mouse_ephys_corr_heatmap.Rdata')

ephys_dist_mat = as.dist(1- ephys_correlation_mat)

index = match(colnames(ephys_correlation_mat), mouse_patch_data$cell_specimen_id)

corr_mouse_subclass = mouse_patch_data$mouse_subclass[index]
corr_primate_subclass = mouse_patch_data$primate_subclass[index]

top_annot = HeatmapAnnotation(mouse_subclass = corr_mouse_subclass, homologous_subclass = corr_primate_subclass,
                              col = list(mouse_subclass = celltype_color_palette, homologous_subclass = celltype_color_palette))



corr_cols <- rev(grDevices::colorRampPalette(RColorBrewer::brewer.pal(11,"RdYlBu"))(100))
Heatmap(ephys_correlation_mat , name = 'spearman', col = corr_cols,
        clustering_distance_rows = ephys_dist_mat, use_raster = T,
        clustering_distance_columns = ephys_dist_mat,
        clustering_method_rows = 'ward.D2', clustering_method_columns = 'ward.D2',
        show_row_names = F, show_column_names = F,
        top_annotation = top_annot, column_title = 'Cell similarity over electrophysiology features')
```

![](figure_plots_with_data_code_files/figure-gfm/ephy_corr_heatmap-1.png)<!-- -->

## Figure 3F

``` r
load(file = 'data_for_plots/exp_corr_heatmap.Rdata')

dist_exp_cor_mat = as.dist(1- exp_cor_mat)

index = match(colnames(exp_cor_mat), rownames(current_metadata))
exp_mouse_subclass = current_metadata$cell_subclass[index]
exp_primate_subclass = current_metadata$primate_subclass[index]


top_annot_exp = HeatmapAnnotation(mouse_subclass = exp_mouse_subclass, homologous_subclass = exp_primate_subclass,
                                   col = list(mouse_subclass = celltype_color_palette, homologous_subclass = celltype_color_palette))

corr_cols <- rev(grDevices::colorRampPalette(RColorBrewer::brewer.pal(11,"RdYlBu"))(100))
Heatmap(exp_cor_mat , name = 'spearman',use_raster = F, col = corr_cols,
        clustering_distance_rows = dist_exp_cor_mat,
        clustering_distance_columns = dist_exp_cor_mat,
        clustering_method_rows = 'ward.D2', clustering_method_columns = 'ward.D2',
        show_row_names = F, show_column_names = F,
        top_annotation = top_annot_exp, column_title = 'Cell similarity over transcriptomic features (top 25 hvgs)')
```

![](figure_plots_with_data_code_files/figure-gfm/exp_corr_heatmap-1.png)<!-- -->

## Figure 3 G and H

``` r
load( file = 'data_for_plots/cge_pax6_distinguishing_features.Rdata')

par(mar = c(4, 10, 2, 1))
barplot(-1*log10(plotting_feature_pval+1e-10), 
        names.arg = names(plotting_feature_pval),
        las = 2, horiz = T, cex.names = .75, main = 'All CGE vs Homologous Pax6')
abline(v = -1*log10(0.05), col = 'red')
```

![](figure_plots_with_data_code_files/figure-gfm/cge_pax6_ephy-1.png)<!-- -->

``` r
ggplot(sig_subset, aes(x = feature, y = value, color = label)) + geom_boxplot() + facet_wrap(~feature, scales = 'free') +
  scale_color_manual(values = c('Pax6' = celltype_color_palette[['Pax6']], 'non-Pax6 CGE' = 'grey20')) + theme_bw()
```

    ## Warning: Removed 5 rows containing non-finite outside the scale range
    ## (`stat_boxplot()`).

![](figure_plots_with_data_code_files/figure-gfm/cge_pax6_ephy-2.png)<!-- -->

## Figure 3I

``` r
load( file = 'data_for_plots/mge_pax6_distinguishing_features.Rdata')


ggplot(sig_subset, aes(x = feature, y = value, color = label)) + geom_boxplot() + facet_wrap(~feature, scales = 'free') +
  scale_color_manual(values = c('Pax6' = celltype_color_palette[['Pax6']], 'MGE' = 'grey20')) + theme_bw()
```

    ## Warning: Removed 15 rows containing non-finite outside the scale range
    ## (`stat_boxplot()`).

![](figure_plots_with_data_code_files/figure-gfm/mge_pax6_ephy-1.png)<!-- -->

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
