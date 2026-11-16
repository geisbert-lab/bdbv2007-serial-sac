rm(list=ls(all.names=TRUE))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(ggrepel))
theme_set(ggpubr::theme_pubr() +
            theme(legend.position="right",
                  legend.background=element_blank()))

# helper variables
outcome.shapes <- c(Succumbed=21, Survived=22)

outcome.cols <- c("Survived"="black", "Succumbed"="#e41a1c")

group.cols <- c("Baseline"="white",
                "Pre-response"="#ffffcc",
                "Early"="#c7e9b4",
                "Middle"="#1d91c0",
                "Late"="#0c2c84",
                "Recovered"="grey70",
                "Terminal" = "black")

cols.reg <- c(Up="#e41a1c",
              Down="#377eb8")

sampling.days <- c(0, 3, 5, 7, 9, 12, 15, 21, 28)

#helper functions 
run.pca <- function(log.expr.mat) {
  
  pca_obj <- log.expr.mat %>%
    t() %>%
    prcomp()
  
  pcs <- summary(pca_obj)$importance["Proportion of Variance", 1:2]
  pcs <- round(100 * pcs)
  pcs <- paste0(names(pcs), " (", pcs, "%)")
  
  contribs <- pca_obj$rotation[, 1:2] %>%
    as.data.frame() %>%
    rownames_to_column("Analyte")
  
  list(
    scores = pca_obj$x,
    loadings = pca_obj$rotation,
    pca_obj = pca_obj,
    pcs = pcs
  )
}

get.results <- function(coef_name, diffexpr, alpha = 0.05) {
  topTable(diffexpr, coef = coef_name, number = Inf) %>%
    as.data.frame() %>%
    rownames_to_column("analyte") %>%
    mutate(
      coef = coef_name,
      Regulation = case_when(
        adj.P.Val < alpha & logFC > 0 ~ "Up",
        adj.P.Val < alpha & logFC < 0 ~ "Down",
        TRUE ~ "None"
      ),
      Regulation = factor(Regulation, levels = c("Up", "None", "Down")),
      Significant = adj.P.Val < alpha,
      lfc = logFC,
      padj = adj.P.Val
    )
}

# import data
data = read_excel("data.xlsx", sheet = "data")
meta = read_excel("data.xlsx", sheet = "meta")

# average replicates
data <- data %>%
  mutate(across(-Sample, ~ as.numeric(str_replace_all(., "[^0-9.]", "")))) %>%
  group_by(`Sample`) %>%
  summarise(across(everything(), ~ mean(., na.rm = TRUE))) %>%
  ungroup() %>%
  as.data.frame()

# disease staging
meta <- meta %>%
  mutate(Group=factor(Group, levels=names(group.cols)),
         NHP=factor(NHP, levels=paste0("CYNO2-", 1:12)))

# sampling plot 
sampling_plot <- meta %>%
  ggplot(aes(DPI, NHP)) +
  geom_line(aes(group=NHP)) +
  geom_point(aes(shape=Outcome, fill=Group), col="black", size=3) +
  scale_shape_manual(values=outcome.shapes) +
  scale_fill_manual(values=group.cols) +
  guides(shape=guide_legend(override.aes=list(fill="black")),
         fill=guide_legend("Stage", override.aes=list(pch=21))) +
  labs(x="Days postinfection") +
  theme(legend.box="horizontal")

ggsave("analysis/samples.png", units="in", width=7.5, height=3, dpi = 600)

# PCA
counts <- data %>%
  pivot_longer(
    cols = -Sample,
    names_to = "variable",
    values_to = "value"
  ) %>%
  pivot_wider(
    names_from = Sample,
    values_from = value
  ) %>%
  column_to_rownames(var = "variable")
rownames(meta) <- meta$Sample
x = intersect(colnames(counts),rownames(meta))
length(x)
counts <- counts[ , x]
counts_numeric <- as.data.frame(lapply(counts, as.numeric))
rownames(counts_numeric) <- rownames(counts)
counts <- as.matrix(counts_numeric)
cmat <- counts
cmat <- as.data.frame(cmat)
# meta$Outcome <- tools::toTitleCase(tolower(meta$Outcome))
cmat <- cmat %>%
  rownames_to_column("Analyte") %>%
  reshape2::melt(id.vars="Analyte",
                 variable.name="Sample",
                 value.name="Concentration") %>%
  left_join(meta, by="Sample") %>%
  group_by(Analyte, NHP, DPI,`Terminal Day`, 
           Outcome, Group) %>%
  summarise(Concentration=mean(Concentration),
            .groups="drop")
meta_pca <- cmat %>%
  select(-Analyte, -Concentration) %>%
  distinct() %>%
  rownames_to_column("ID") %>%
  mutate(ID=paste0("sample", ID)) %>%
  as.data.frame()
rownames(meta_pca) <- meta_pca$ID

#re-format concentration matrix using the new sample ID
cmat <- cmat %>%
  left_join(meta_pca) %>%
  reshape2::dcast(Analyte ~ ID, value.var="Concentration") %>%
  column_to_rownames("Analyte") %>%
  as.matrix()

# align rows and columns
x <- intersect(rownames(meta_pca), colnames(cmat))
meta_pca <- meta_pca[x, ]
cmat <- cmat[, x]
rm(x)

logcmat <- log2(cmat)
modmat <- model.matrix(~0 + Group, data=meta_pca)
pca <- logcmat %>% run.pca()

df_pca <- pca$scores %>%
  as.data.frame() %>%
  rownames_to_column("ID") %>%
  left_join(meta_pca, by = "ID")
df_pca$Group  <- factor(df_pca$Group)
df_pca$Outcome <- factor(df_pca$Outcome) 
df_pca <- df_pca %>%
  mutate(is_terminal = (DPI == `Terminal Day` & DPI != 28),
       fill_label = ifelse(is_terminal, "Terminal", as.character(Group)))

# plot pca 
pca_plot <- ggplot(df_pca, aes(PC1, PC2)) +
  geom_point(aes(shape = Outcome, fill = Group), size = 3) +
  scale_shape_manual(values=outcome.shapes) +
  scale_fill_manual(values=group.cols) +
  labs(x = pca$pcs[1], y = pca$pcs[2], title = "PCA clustering") +
  guides(shape=guide_legend(override.aes=list(fill="black")),
         fill=guide_legend("Stage", override.aes=list(pch=21))) +
  theme(axis.text=element_blank(),
        axis.ticks=element_blank(),
        legend.box="horizontal")

ggsave(filename = "analysis/pca.png", plot = pca_plot, dpi = 600,
       units = "in", width = 7.5, height = 3)

# plot spaghetti
analyte_names <- setdiff(names(data), "Sample")
output_dir <- "analysis/spaghetti_plots"
sampling_days <- c(0, 3, 5, 7, 9, 12, 15, 21, 28)
custom_labels <- ifelse(0:28 %in% sampling_days, 0:28, "")
for (analyte in analyte_names) {
  plot_data_long <- data %>%
    select(Sample, all_of(analyte)) %>% 
    left_join(meta, by = "Sample") 
  
  spaghetti_df <- plot_data_long %>%
    mutate(log2_val = log2(.data[[analyte]])) %>% 
    group_by(NHP) %>%
    mutate(
      baseline_val = log2_val[DPI == 0], 
      log2_fold_change = log2_val - baseline_val) %>%
    ungroup() %>%
    filter(!is.na(log2_fold_change))
  
  spaghetti_plot <- ggplot(spaghetti_df, aes(x = DPI, y = log2_fold_change)) +
    geom_line(aes(group = NHP, color = Outcome), linewidth = 0.6) +
    scale_color_manual(
      name = "Outcome", 
      values = outcome.cols
    ) +
    scale_x_continuous(
      breaks = sampling.days,
      labels = sampling.days, 
      limits = c(0, 28)
    ) +
    labs(
      x = "Days postinfection",
      y = "Fold Change (log2)",
      title = paste(analyte)
    ) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      legend.position = "none"
    )
  
  clean_name <- gsub("[^[:alnum:]_]", "_", analyte)
  file_path <- file.path(output_dir, paste0("spaghetti_", clean_name, ".png"))
  
  ggsave(filename = file_path, plot = spaghetti_plot, dpi = 600,
         units = "in", width = 7.5, height = 3)
}

# contrasts from succumbed and survived versus baseline
dpis.time <- c("0" = 0, "1" = 1, "3" = 3, "5" = 5, "7" = 7, 
               "9" = 9,"12" = 12, "15" = 15, "16" = 16, "17" = 17, "21" = 21, "28" = 28)

modmat <- interaction(meta_pca$Outcome, meta_pca$DPI) %>%
  droplevels()
modmat <- model.matrix(~0 + modmat)

meta_pca$Outcome <- factor(trimws(meta_pca$Outcome))
meta_pca$DPI <- as.integer(meta_pca$DPI)
meta_pca <- meta_pca[colnames(cmat), ]
grp <- interaction(meta_pca$Outcome, meta_pca$DPI, sep = ".", drop = TRUE)
modmat <- model.matrix(~0 + grp)
colnames(modmat) <- levels(grp)
colnames(modmat)
# format column names so they're human readable
colnames(modmat) <- colnames(modmat) %>%
  str_remove("^modmat") %>%
  str_remove("/.*")

# define all contrasts. Use all pre-infection (not bcx) baseline to increase N; combine all D0 samples (succumbed & survived) to get better distribution representation
ctrs <- makeContrasts(
  # 1. Survived over time
  Survived.1  - (Survived.0 + Succumbed.0)/2,
  Survived.3  - (Survived.0 + Succumbed.0)/2,
  Survived.5  - (Survived.0 + Succumbed.0)/2,
  Survived.7  - (Survived.0 + Succumbed.0)/2,
  Survived.9  - (Survived.0 + Succumbed.0)/2,
  Survived.12 - (Survived.0 + Succumbed.0)/2,
  Survived.15 - (Survived.0 + Succumbed.0)/2,
  Survived.21 - (Survived.0 + Succumbed.0)/2,
  Survived.28 - (Survived.0 + Succumbed.0)/2,
  
  # 2. Succumbed over time
  Succumbed.1  - (Survived.0 + Succumbed.0)/2,
  Succumbed.3  - (Survived.0 + Succumbed.0)/2,
  Succumbed.5  - (Survived.0 + Succumbed.0)/2,
  Succumbed.7  - (Survived.0 + Succumbed.0)/2,
  Succumbed.9  - (Survived.0 + Succumbed.0)/2,
  Succumbed.12 - (Survived.0 + Succumbed.0)/2,
  Succumbed.15 - (Survived.0 + Succumbed.0)/2,
  Succumbed.16 - (Survived.0 + Succumbed.0)/2,
  Succumbed.17 - (Survived.0 + Succumbed.0)/2,
  
  # 3. Survived vs Succumbed, controlling for BL
  Survived.0 - Succumbed.0,
  (Survived.1  - Survived.0) - (Succumbed.1  - Succumbed.0),
  (Survived.3  - Survived.0) - (Succumbed.3  - Succumbed.0),
  (Survived.5  - Survived.0) - (Succumbed.5  - Succumbed.0),
  (Survived.7  - Survived.0) - (Succumbed.7  - Succumbed.0),
  (Survived.9  - Survived.0) - (Succumbed.9  - Succumbed.0),
  (Survived.12 - Survived.0) - (Succumbed.12 - Succumbed.0),
  (Survived.15 - Survived.0) - (Succumbed.15 - Succumbed.0),
  
  levels = colnames(modmat)
)

# run limma pipeline
de <- cmat %>%
  voom(modmat) %>%
  lmFit(modmat) %>%
  contrasts.fit(ctrs) %>%
  eBayes(robust=TRUE)

# get DE tables for each comparison over time (comps. to 0 DPI)
res.time <- colnames(de)[1:18] %>%
  lapply(get.results, diffexpr = de) %>%
  do.call(rbind, .) %>%
  mutate(
    Outcome = str_extract(coef, "^[^\\.]+"),
    DPI = as.integer(str_extract(coef, "(?<=\\.)[^ ]+")),
    Daterange = factor(DPI, levels = dpis.time, labels = names(dpis.time))
  ) %>%
  select(analyte, Outcome, DPI, Regulation, lfc, padj, Significant)
write.csv(res.time, "results-timepoint.csv", row.names=FALSE)

# repeat for comparison between outcome at each time point
res.virs <- colnames(de)[19:26] %>%
  lapply(get.results, diffexpr = de) %>%
  do.call(rbind, .) %>%
  mutate(
    Regulation = factor(Regulation, levels = c("Up", "None", "Down"),
                        labels = c("Survived", "None", "Succumbed")),
    DPI = as.integer(str_extract(coef, "(?<=\\.)[^ ]+")),
    Daterange = factor(DPI, levels = dpis.time, labels = names(dpis.time))
  ) %>%
  select(analyte, DPI, Regulation, lfc, padj, Significant)
write.csv(res.virs, "results-outcome.csv", row.names=FALSE)

# clean up
rm(ctrs, de, modmat)

#volcano plots
res.time <- res.time %>%
  mutate(
    Phase = case_when(
      DPI %in% c(3, 5)          ~ "Early",
      DPI %in% c(7, 9) ~ "Middle",
      DPI %in% c(12, 15, 16, 17, 21)       ~ "Late",
      TRUE ~ NA_character_
    )
  )

res.time$Phase <- factor(res.time$Phase, levels = c("Early", "Middle", "Late"))

for (out in unique(res.time$Outcome)) {
  
  vdf <- res.time %>%
    filter(Outcome == out, !is.na(Phase)) %>%
    mutate(
      neg_log10_padj = -log10(padj),
      color_label = ifelse(padj < 0.05 & lfc > 0, "Up",
                           ifelse(padj < 0.05 & lfc < 0, "Down", "NS"))
    )
  
  if (nrow(vdf) == 0) next
  
  top5_up <- vdf %>%
    filter(padj < 0.05, lfc > 0) %>%
    group_by(Phase, analyte) %>%
    slice_min(order_by = padj, n = 1) %>%
    ungroup() %>%
    group_by(Phase) %>%
    slice_max(order_by = lfc, n = 5) %>%
    ungroup()
  
  top5_down <- vdf %>%
    filter(padj < 0.05, lfc < 0) %>%
    group_by(Phase, analyte) %>%
    slice_min(order_by = padj, n = 1) %>%
    ungroup() %>%
    group_by(Phase) %>%
    slice_min(order_by = lfc, n = 5) %>%  # most negative LFC
    ungroup()
  
  top5 <- bind_rows(top5_up, top5_down)
  
  volc <- ggplot(vdf, aes(x = lfc, y = neg_log10_padj, color = color_label)) +
    geom_point(aes(fill = color_label), shape = 21, size = 1.75, alpha = 0.75, color = "black", stroke = 0.3) +
    scale_fill_manual(values = c("Up" = "#e41a1c", "Down" = "#377eb8", "NS" = "grey50")) +
    geom_text_repel(data = top5, aes(label = analyte),
                    size = 2.5, max.overlaps = 10,
                    color = "black", fontface = "bold",
                    segment.color = "black",
                    segment.size = 0.4,
                    min.segment.length = 0,
                    box.padding = 0.5,
                    point.padding = 0.1,
                    force = 2,
                    force_pull = 0.1) +
    facet_wrap(~Phase, nrow = 1) +
    labs(
      title = paste0(out, " vs Baseline"),
      x = "Fold Change (log2)",
      y = "-log10(adj p-value)"
    ) +
    theme_classic(base_size = 11) +
    theme(
      strip.text = element_text(size = 14),
      axis.title = element_text(size = 11),
      axis.text  = element_text(size = 11),
      plot.title = element_text(size = 16),
      legend.position = ("none"),
      strip.background = element_rect(fill = "grey90", color = "black")
    )
  
  output_dir_volc <- "analysis/volcano_plots"
  dir.create(output_dir_volc, recursive = TRUE)
  fname <- file.path(output_dir_volc, paste0("volcano_", out, "_by_phase.jpg"))
  ggsave(fname, volc, dpi = 600, units = "in", width = 3.75, height = 2.5)
}

for (out in unique(res.time$Outcome)) {
  
  vdf <- res.time %>%
    filter(Outcome == out) %>%
    mutate(
      neg_log10_padj = -log10(padj),
      color_label = ifelse(padj < 0.05 & lfc > 0, "Up",
                           ifelse(padj < 0.05 & lfc < 0, "Down", "NS")),
      DPI = factor(DPI, levels = sort(unique(DPI)))
    )
  
  if (nrow(vdf) == 0) next
  
  top5_up <- vdf %>%
    filter(padj < 0.05, lfc > 0) %>%
    group_by(DPI) %>%
    slice_max(order_by = lfc, n = 5) %>%
    ungroup()
  
  top5_down <- vdf %>%
    filter(padj < 0.05, lfc < 0) %>%
    group_by(DPI) %>%
    slice_min(order_by = lfc, n = 5) %>%
    ungroup()
  
  top5 <- bind_rows(top5_up, top5_down)
  
  volc <- ggplot(vdf, aes(x = lfc, y = neg_log10_padj)) +
    geom_point(aes(fill = color_label), shape = 21, size = 1.5, alpha = 0.75, color = "black", stroke = 0.3) +
    scale_fill_manual(values = c("Up" = "#e41a1c", "Down" = "#377eb8", "NS" = "grey50")) +
    geom_text_repel(data = top5, aes(label = analyte),
                    size = 2.5, max.overlaps = 20,
                    color = "black", fontface = "bold",
                    segment.color = "black",
                    segment.size = 0.3,
                    min.segment.length = 0,
                    box.padding = 0.8,
                    point.padding = 0.5,
                    force = 2,
                    force_pull = 0.1) +
    facet_wrap(~DPI, nrow = 2, scales = "free_x") +
    labs(
      title = paste0(out, " vs Baseline — All Timepoints"),
      x = "log2 Fold Change",
      y = "-log10(adj p-value)"
    ) +
    theme_classic(base_size = 11) +
    theme(
      axis.title = element_text(size = 16),
      axis.text  = element_text(size = 14),
      strip.text = element_text(size = 14, face = "bold"),
      strip.background = element_rect(fill = "grey85", color = "black"),
      plot.title = element_text(size = 18, face = "bold"),
      legend.position = "none"
    )
  
  output_dir_volc <- "analysis/volcano_plots"
  dir.create(output_dir_volc, recursive = TRUE)
  fname <- file.path(output_dir_volc, paste0("volcano_", out, "All Timepoints.png"))
  ggsave(fname, volc, dpi = 600, units = "in", width = 7.5, height = 5)
  # print(volc)
}



