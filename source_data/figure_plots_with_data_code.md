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

## Figure 2 H

``` r
load(file = 'data_for_plots/mouse_meta_markers.Rdata')

plot_pareto_summary(mouse_subclass_meta_markers %>% filter(cell_type %in% c('Pax6','Sncg')), min_recurrence=4) +
  theme_bw() + scale_color_manual(values = celltype_color_palette) + ggtitle('Homologous Mouse subclass markers') + ylim(.5, 1) + xlim(0,5)
```

![](figure_plots_with_data_code_files/figure-gfm/mouse_pareto-1.png)<!-- -->
