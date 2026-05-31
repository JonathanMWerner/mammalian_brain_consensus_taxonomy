library(magrittr)
library(MetaMarkers)
library(stringr)
library(readxl)
library(readr)
library(dplyr)
library(tibble)
library(SingleCellExperiment)
library(zellkonverter)
library(Matrix)

# -----------------------------
# File paths
# -----------------------------

cell_metadata_file <- "/vault/lfrench/mouse_brain_consensus_taxonomy_spatial/data/Zeng/MERFISH-C57BL6J-638850/20230830/cell_metadata.csv"
ccf_coordinates_file <- "/vault/lfrench/mouse_brain_consensus_taxonomy_spatial/data/Zeng/from_API/ccf_coordinates_MERFISH-C57BL6J-638850.csv.gz"
ccf_expansions_file <- "/vault/lfrench/mouse_brain_consensus_taxonomy_spatial/data/Allen_CCFv3/1-s2.0-S0092867420304025-mmc2.xlsx"
cluster_annotations_file <- "/vault/lfrench/mouse_brain_consensus_taxonomy_spatial/data/Zeng/from_AWS/AIT21.0/AIT21_annotation_freeze_081523.tsv"
interneuron_metadata_file <- "/inkwell03/werner/mouse_brain_consensus_taxonomy/20_3_25_interneuron_mouse_metadata.Rdata"
h5ad_file <- "/vault/lfrench/whole_mouse_brain/zeng/MERFISH-C57BL6J-638850/20230630/C57BL6J-638850-raw.h5ad"


# -----------------------------
# Load cell metadata and CCF annotations
# -----------------------------

Zeng_calls_direct <- read_csv(cell_metadata_file, show_col_types = FALSE)

ccf_annotations <- read_csv(ccf_coordinates_file, show_col_types = FALSE)

Zeng_calls_direct <- Zeng_calls_direct %>%
  left_join(ccf_annotations, by = join_by(cell_label))


# -----------------------------
# Add CCF structure expansions and cortical layer labels
# -----------------------------

ccf_expansions <- read_xlsx(ccf_expansions_file, skip = 1) %>%
  select(
    parcellation_substructure = abbreviation,
    parcellation_substructure_expansion = `full structure name`
  ) %>%
  mutate(
    cortical_layer = parcellation_substructure_expansion %>%
      str_extract(regex(
        "\\blayer\\s*[0-6](?:\\s*/\\s*[0-6])?[ab]?\\b",
        ignore_case = TRUE
      )) %>%
      str_remove(regex("^layer\\s*", ignore_case = TRUE)) %>%
      str_replace_all("\\s+", "")
  )

Zeng_calls_direct <- Zeng_calls_direct %>%
  left_join(ccf_expansions, by = "parcellation_substructure") %>%
  mutate(
    cortical_layer = if_else(
      parcellation_division == "Isocortex",
      cortical_layer,
      NA_character_
    )
  )


# -----------------------------
# Add cluster annotations
# -----------------------------

cluster_annotations <- read_tsv(cluster_annotations_file, show_col_types = FALSE) %>%
  select(
    cluster_alias = cl,
    subclass_id,
    subclass_label,
    supertype_label,
    class_label
  )

Zeng_calls_direct <- Zeng_calls_direct %>%
  inner_join(cluster_annotations, by = "cluster_alias")


# -----------------------------
# Add primate subclass mappings
# -----------------------------

load(interneuron_metadata_file, verbose = TRUE)

primate_mappings <- interneuron_metadata %>%
  as_tibble() %>%
  filter(study %in% c("zeng_10xv2_wba", "zeng_10xv3_wba")) %>%
  select(
    supertype_label = cell_cluster,
    primate_subclass
  ) %>%
  distinct()

Zeng_calls_direct <- Zeng_calls_direct %>%
  left_join(primate_mappings, by = "supertype_label")


# -----------------------------
# Prepare annotations to add to SCE colData
# -----------------------------

ccf_ann_sub <- Zeng_calls_direct %>%
  select(
    cell_label,
    cluster_alias,
    primate_subclass,
    subclass_label,
    supertype_label,
    class_label,
    parcellation_division,
    parcellation_structure,
    parcellation_substructure,
    parcellation_substructure_expansion,
    cortical_layer
  )


# -----------------------------
# Load SCE and add annotations
# -----------------------------

sce <- readH5AD(h5ad_file)

cell_meta_merged <- as.data.frame(colData(sce)) %>%
  rownames_to_column("cell_label") %>%
  left_join(ccf_ann_sub, by = "cell_label") %>%
  column_to_rownames("cell_label")

colData(sce) <- S4Vectors::DataFrame(cell_meta_merged)

# -----------------------------
# Filter to cells with primate subclass annotations
# -----------------------------

sce <- sce[, !is.na(sce$primate_subclass)]

#mark layer one cells
sce$is_layer1 <- case_when(
  !is.na(sce$cortical_layer) & sce$cortical_layer == "1" ~ "layer 1",
  !is.na(sce$cortical_layer) ~ "other layers",
  TRUE ~ NA_character_
)


# -----------------------------
# Remove blank / control genes
# -----------------------------
blank_rows <- stringr::str_detect(rownames(sce), "^Blank-") |
  stringr::str_detect(rowData(sce)$gene_symbol, "^Blank-")

blank_rows[is.na(blank_rows)] <- FALSE

message("Removing ", sum(blank_rows), " blank genes")

sce <- sce[!blank_rows, ]

# -----------------------------
# CPM normalization using metamarkers
# -----------------------------
assay(sce, "cpm") = convert_to_cpm(assay(sce))

#dim: 550 119676 
sce



expr_mat <- assay(sce, "cpm")

primate_subclasses <- sce$primate_subclass %>%
  unique() %>%
  sort()

markers <- list()

for (subclass in primate_subclasses) {
  message("Computing layer 1 markers for: ", subclass)
  
  cells_keep <- sce$primate_subclass == subclass &
    !is.na(sce$is_layer1)
  
  labels <- sce$is_layer1[cells_keep]
  
  # Need both layer1 and not_layer1 cells for marker testing
  if (length(unique(labels)) < 2) {
    message("  Skipping: only one is_layer1 class present")
    next
  }
  
  subclass_key <- make.names(subclass)
  
  markers[[subclass_key]] <- compute_markers(
    expr_mat[, cells_keep, drop = FALSE],
    labels
  )
}

meta_markers <- make_meta_markers(
  markers,
  detailed_stats = TRUE
)

gene_lookup <- rowData(sce) %>%
  as.data.frame() %>%   rownames_to_column("gene") %>%
  tibble() 

meta_markers %<>% left_join(gene_lookup)
meta_markers %<>% dplyr::select(cell_type, gene_symbol, everything()) %>% dplyr::select(-transcript_identifier)

meta_markers

export_meta_markers(meta_markers, 
                    "/home/werner/projects/mouse_brain_consensus_taxonomy/data/L1_mouse_meta_markers.csv",
                    names(markers))
