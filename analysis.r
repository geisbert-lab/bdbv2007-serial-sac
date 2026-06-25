#!/usr/bin/env Rscript

## setup -----------------------------------------------------------------------
rm(list=ls(all.names=TRUE))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(edgeR))
theme_set(ggpubr::theme_pubr() +
            theme(legend.position="right",
                  legend.background=element_blank()))

# helper variables
outcome.shapes <- c(Succumbed=21, Survived=22)
outcome.cols <- c(Succumbed="#e41a1c", Survived="black")
daterange.cols <- c("Baseline"="white", 
                    "1 DPI"="#ffffcc",
                    "3 DPI"="#c7e9b4",
                    "5 DPI"="#7fcdbb",
                    "7 DPI"="#41b6c4",
                    "9 DPI"="#1d91c0",
                    "12 DPI"="#225ea8",
                    "15 DPI"="#0c2c84",
                    "Terminal"="black",
                    "21 DPI"="grey40",
                    "28 DPI"="grey70")
stage.cols <- c("Baseline"="white",
                "Pre-response"="#ffffcc",
                "Early"="#c7e9b4",
                "Middle"="#1d91c0",
                "Late"="#0c2c84",
                "Recovered"="grey70")
regulation.cols <- c(Up="#e41a1c", Down="#377eb8")
sample.days <- c(0, 1, 3, 5, 7, 9, 12, 15, 21, 28)

# figure lists
fig.main <- list()
fig.sup1 <- list()
fig.sup2 <- list()

# helper functions
get.results <- function(coef, diffexpr, psig=0.05, fsig=1) {
  # get results and loosely format
  rmat <- topTable(de, coef=coef, number=Inf, sort.by="logFC") %>%
          rownames_to_column("Gene") %>%
          mutate(coef=coef,
                 psig=(adj.P.Val < psig),
                 fsig=(abs(logFC) > fsig),
                 Significant=(psig & fsig)) %>%
          rename(lfc=logFC, padj=adj.P.Val) %>%
          select(Gene, lfc, padj, Significant, coef)
  # make a "regulation" column
  rmat$Regulation <- "None"
  rmat$Regulation[rmat$Significant & rmat$lfc > 0] <- "Up"
  rmat$Regulation[rmat$Significant & rmat$lfc < 0] <- "Down"
  # return the formatted data.frame
  return(rmat)
}

## inputs ----------------------------------------------------------------------
# metadata
meta <- readxl::read_excel("data.xlsx", "animals") # individual NHP info
meta <- readxl::read_excel("data.xlsx", "nanostring.samples") %>%
        left_join(meta, by="NHP") %>%
        mutate(Daterange=factor(Daterange, levels=names(daterange.cols)),
               Stage=factor(Stage, levels=names(stage.cols)),
               NHP=factor(NHP, levels=paste0("CYNO2-", 1:12)))

# mRNA counts (Nanostring)
mrna <- readxl::read_excel("data.xlsx", "nanostring.thresholded") %>%
        as.data.frame() %>%
        # remove spike-in probes
        filter(!str_detect(Gene, "^POS|NEG")) %>%
        column_to_rownames("Gene") %>%
        as.matrix()

# protein concentration (LEGENDplex)
prot <- readxl::read_excel("data.xlsx", "legendplex.concentration") %>%
        reshape2::melt(id.vars="ID", 
                       variable.name="Analyte", 
                       value.name="Concentration") %>%
        group_by(Analyte, ID) %>%
        summarise(Concentration=mean(Concentration),
                  .groups="drop") 
# update protein IDs to align with mRNA IDs
prot <- readxl::read_excel("data.xlsx", "legendplex.samples") %>%
        select(ID, NHP, DPI) %>%
        left_join(rename(meta, NewID=ID), by=c("NHP", "DPI")) %>%
        select(ID, NewID) %>%
        # add in concentration data
        left_join(prot, by="ID") %>%
        # pivot to wide format and use new IDs for column names
        reshape2::dcast(Analyte ~ NewID, value.var="Concentration") %>%
        column_to_rownames("Analyte") %>%
        as.matrix()

# add LP batch info to metadata
meta <- readxl::read_excel("data.xlsx", "legendplex.samples") %>%
        select(NHP, DPI, Batch) %>%
        left_join(meta, by=c("NHP", "DPI")) %>%
        arrange(ID)

# align rows and columns
meta <- filter(meta, ID %in% colnames(mrna), ID %in% colnames(prot))
mrna <- mrna[, meta$ID]
prot <- prot[, meta$ID]

# format mRNA counts
mrna <- mrna %>%
        DGEList() %>%
        normLibSizes()

# save non-log CPM for CIBERSORTx
mrna %>%
  cpm(log=FALSE) %>%
  as.data.frame() %>%
  rownames_to_column("Gene") %>%
  write.table("analysis/cibersort-input.tsv", sep="\t", row.names=FALSE)

## clustering with PCA ---------------------------------------------------------
# build model for LP batch correction
pca <- model.matrix(~0 + Stage, data=meta)

# remove LP batch effect, then combine to logCPM, then run PCA
pca <- prot %>%
        log2() %>%
        removeBatchEffect(batch=meta$Batch, design=pca) %>%
        rbind(cpm(mrna, log=TRUE)) %>%
        t() %>%
        prcomp()
# get PC contributions
pcs <- summary(pca)$importance["Proportion of Variance", 1:2]
pcs <- paste0(names(pcs), " (", round(100*pcs), "%)")
# get contributions
contribs <- pca$rotation %>%
            as.data.frame() %>%
            rownames_to_column("Gene") %>%
            select(Gene, PC1, PC2)
# plot PCA by DPI
pca <- pca$x %>%
        as.data.frame() %>%
        rownames_to_column("ID") %>%
        select(ID, PC1, PC2) %>%
        left_join(meta, by="ID")
fig.sup1$a <- pca %>%
              ggplot(aes(PC1, PC2)) +
              geom_hline(yintercept=0, linetype=3, col="lightgrey") +
              geom_vline(xintercept=0, linetype=3, col="lightgrey") +
              geom_point(aes(shape=Outcome, fill=Daterange), size=3) +
              scale_shape_manual(values=outcome.shapes) +
              scale_fill_manual(values=daterange.cols) +
              labs(x=pcs[1], y=pcs[2], title="PCA clustering by DPI") +
              guides(shape=guide_legend(override.aes=list(fill="black")),
                     fill=guide_legend("DPI", override.aes=list(pch=21), 
                                       order=1)) +
              theme(legend.box="horizontal")
# plot PCA by stage
fig.sup1$b <- pca %>%
              ggplot(aes(PC1, PC2)) +
              geom_hline(yintercept=0, linetype=3, col="lightgrey") +
              geom_vline(xintercept=0, linetype=3, col="lightgrey") +
              geom_point(aes(shape=Outcome, fill=Stage), size=3) +
              scale_shape_manual(values=outcome.shapes) +
              scale_fill_manual(values=stage.cols) +
              labs(x=pcs[1], y=pcs[2], title="PCA clustering by disease stage") +
              guides(shape=guide_legend(override.aes=list(fill="black")),
                     fill=guide_legend("Stage", override.aes=list(pch=21), 
                                       order=1)) +
              theme(legend.box="horizontal")
# plot the top 30 contributors by calculating hypoteneuse/magnitude
contribs %>%
  mutate(Magnitude=sqrt(PC1^2 + PC2^2)) %>%
  slice_max(n=30, order_by=Magnitude, with_ties=FALSE) %>%
  ggplot() +
  geom_vline(xintercept=0, linetype=3, col="lightgrey") +
  geom_hline(yintercept=0, linetype=3, col="lightgrey") +
  geom_segment(aes(x=0, y=0, xend=PC1, yend=PC2), col="darkgrey",
               arrow=arrow(length=unit(0.05, "in"))) +
  ggrepel::geom_text_repel(aes(PC1, PC2, label=Gene), 
                           size=3, col="black") +
  labs(x=pcs[1], y=pcs[2], title="PCA contributors")

# samples per stage for comparison
meta %>%
  group_by(Outcome, Stage) %>%
  summarise(Samples=n(),
            .groups="drop")

# clean up
rm(contribs, pca, pcs)

## differential expression analyses --------------------------------------------
# running the same comparisons for both mRNA and protein
grps <- interaction(meta$Outcome, meta$Stage) %>%
        droplevels()
grps <- factor(str_remove(grps, "-"), 
        levels=str_remove(levels(grps), "-"))
# using all pre-infection baseline to increase N
ctrs <- c(#1. survived over time
          "Survived.Preresponse-(Survived.Baseline+Succumbed.Baseline)/2",
          "Survived.Early-(Survived.Baseline+Succumbed.Baseline)/2",
          "Survived.Middle-(Survived.Baseline+Succumbed.Baseline)/2",
          "Survived.Late-(Survived.Baseline+Succumbed.Baseline)/2",
          "Survived.Recovered-(Survived.Baseline+Succumbed.Baseline)/2",
          #2. succumbed over time
          "Succumbed.Preresponse-(Survived.Baseline+Succumbed.Baseline)/2",
          "Succumbed.Early-(Survived.Baseline+Succumbed.Baseline)/2",
          "Succumbed.Middle-(Survived.Baseline+Succumbed.Baseline)/2",
          "Succumbed.Late-(Survived.Baseline+Succumbed.Baseline)/2",
          #3. survived vs. succumbed
          "Survived.Baseline-Succumbed.Baseline",
          "Survived.Preresponse-Succumbed.Preresponse",
          "Survived.Early-Succumbed.Early",
          "Survived.Middle-Succumbed.Middle",
          "Survived.Late-Succumbed.Late")

# mRNA: build model matrix; no batch effects here
modmat <- model.matrix(~0 + grps)
# format column names so they're human readable
colnames(modmat) <- str_remove(colnames(modmat), "^grps")
# run limma pipeline
de <- mrna %>%
      voom(modmat) %>%
      lmFit(modmat) %>%
      contrasts.fit(makeContrasts(contrasts=ctrs, levels=levels(grps))) %>%
      eBayes(robust=TRUE)
# get DE tables for each comparison over time (comps. to 0 DPI)
ns.res <- list()
ns.res$time <- colnames(de)[1:9] %>%
               lapply(get.results, diffexpr=de) %>%
               do.call(rbind, .) %>%
                # format outcome and daterange as factors
                mutate(Outcome=str_extract(coef, "^[^\\.]+"),
                       Group=str_extract(coef, "(?<=\\.)[^- ]+"),
                       Group=if_else(Group=="Preresponse", 
                                     "Pre-response", 
                                     Group),
                       Group=factor(Group, levels=names(stage.cols))) %>%
                select(Gene, Outcome, Group, Regulation,
                       lfc, padj, Significant)
write.csv(ns.res$time, "analysis/nanostring-timepoint.csv", row.names=FALSE)
# repeat for comparison between outcome at each time point
ns.res$outcome <- colnames(de)[c(10:14)] %>%
                  lapply(get.results, diffexpr=de) %>%
                  do.call(rbind, .) %>%
                  # change regulation to "higher expr in surv/succ"
                  # and format daterange
                  mutate(Regulation=factor(Regulation, 
                                           levels=c("Up", "None", "Down"),
                                           labels=c("Survived", "None", "Succumbed")),
                         Group=str_extract(coef, "(?<=\\.)[^ -]+"),
                         Group=if_else(Group=="Preresponse", 
                                       "Pre-response", 
                                       Group),
                         Group=factor(Group, levels=names(stage.cols))) %>%
                  select(Gene, Group, Regulation, lfc, padj, Significant)
write.csv(ns.res$outcome, "analysis/nanostring-outcome.csv", row.names=FALSE)
rm(modmat, de)

# LEGENDplex: build model matrix accounting for batch effect
modmat <- model.matrix(~0 + grps + meta$Batch)
# format column names so they're human readable
colnames(modmat) <- colnames(modmat) %>%
                    str_remove("^grps") %>%
                    str_remove("^meta\\$")
# run limma pipeline
de <- prot %>%
      log2() %>%
      lmFit(modmat) %>%
      contrasts.fit(makeContrasts(contrasts=ctrs, levels=colnames(modmat))) %>%
      eBayes(robust=TRUE)
# get DE tables for each comparison over time (comps. to 0 DPI)
lp.res <- list()
lp.res$time <- colnames(de)[1:9] %>%
                lapply(get.results, diffexpr=de) %>%
                do.call(rbind, .) %>%
                # format outcome and daterange as factors
                mutate(Outcome=str_extract(coef, "^[^\\.]+"),
                       Group=str_extract(coef, "(?<=\\.)[^ -]+"),
                       Group=if_else(Group=="Preresponse", 
                                     "Pre-response", 
                                     Group),
                       Group=factor(Group, levels=names(stage.cols))) %>%
                select(Gene, Outcome, Group, Regulation,
                       lfc, padj, Significant)
write.csv(lp.res$time, "analysis/legendplex-timepoint.csv", row.names=FALSE)
# repeat for comparison between outcome at each time point
lp.res$outcome <- colnames(de)[c(10:14)] %>%
                  lapply(get.results, diffexpr=de) %>%
                  do.call(rbind, .) %>%
                  # change regulation to "higher expr in surv/succ"
                  # and format daterange
                  mutate(Regulation=factor(Regulation, 
                                           levels=c("Up", "None", "Down"),
                                           labels=c("Survived", "None", "Succumbed")),
                         Group=str_extract(coef, "(?<=\\.)[^ -]+"),
                         Group=if_else(Group=="Preresponse", 
                                       "Pre-response", 
                                       Group),
                         Group=factor(Group, levels=names(stage.cols))) %>%
                  select(Gene, Group, Regulation, lfc, padj, Significant)
write.csv(lp.res$outcome, "analysis/legendplex-outcome.csv", row.names=FALSE)
rm(modmat, de, grps, ctrs)

## plot DEGs -------------------------------------------------------------------
# total DEGs over time
# first, generate our groups
x <- expand.grid(Group=names(stage.cols),
                 Outcome=names(outcome.cols),
                 Regulation=names(regulation.cols)) %>%
      filter(!(Outcome=="Succumbed" & Group=="Recovered"))
# Nanostring
fig.main$a <- ns.res$time %>%
              filter(Significant) %>%
              group_by(Outcome, Group, Regulation) %>%
              summarise(Genes=n(),
                        .groups="drop") %>%
              full_join(x, by=c("Outcome", "Group", "Regulation")) %>%
              replace_na(list(Genes=0)) %>%
              ggplot(aes(Group, Genes)) +
              geom_line(aes(linetype=Outcome, 
                            group=interaction(Outcome, Regulation))) +
              geom_point(aes(shape=Outcome, fill=Regulation), size=2) +
              scale_shape_manual(NULL, values=outcome.shapes) +
              scale_fill_manual(NULL, values=regulation.cols) +
              scale_y_continuous("Total DE mRNAs", limits=c(NA, 150)) +
              guides(shape=guide_legend(override.aes=list(fill="black"),
                                        order=1),
                     linetype=guide_legend(title=NULL, order=1),
                     fill=guide_legend(override.aes=list(pch=21))) +
              labs(x="Disease stage", title="DE mRNAs over time") +
              theme(axis.text.x=element_text(angle=30, hjust=1),
                    legend.position=c(0.22, 0.57))
# LEGENDplex
fig.main$b <- lp.res$time %>%
              filter(Significant) %>%
              group_by(Outcome, Group, Regulation) %>%
              summarise(Genes=n(),
                        .groups="drop") %>%
              full_join(x, by=c("Outcome", "Group", "Regulation")) %>%
              replace_na(list(Genes=0)) %>%
              ggplot(aes(Group, Genes)) +
              geom_line(aes(linetype=Outcome, 
                            group=interaction(Outcome, Regulation))) +
              geom_point(aes(shape=Outcome, fill=Regulation), size=2) +
              scale_shape_manual(NULL, values=outcome.shapes) +
              scale_fill_manual(NULL, values=regulation.cols) +
              scale_y_continuous("Total DE proteins", limits=c(NA, 30)) +
              guides(shape=guide_legend(override.aes=list(fill="black"), 
                                        order=1),
                     linetype=guide_legend(title=NULL, order=1),
                     fill=guide_legend(override.aes=list(pch=21))) +
              labs(x="Disease stage", title="DE proteins over time") +
              theme(axis.text.x=element_text(angle=30, hjust=1),
                    legend.position=c(0.22, 0.57))
rm(x)

# volcano nanostring
topgenes <- ns.res$time %>%
            filter(Significant) %>%
            group_by(Regulation, Group, Outcome) %>%
            top_n(n=5, wt=abs(lfc)) %>%
            ungroup()
fig.sup2$a <- ns.res$time %>%
              ggplot(aes(lfc, -log10(padj))) +
              geom_point(aes(col=Regulation, alpha=Significant, size=Regulation)) +
              scale_size_manual(values=c(Up=1, Down=1), na.value=0.5) +
              scale_alpha_manual(values=c("TRUE"=0.8, "FALSE"=0.5)) +
              scale_color_manual(values=regulation.cols, na.value="grey80") + 
              ggrepel::geom_text_repel(data=topgenes, aes(label=Gene), size=2) +
              scale_x_continuous("Fold change (log2)", limits=c(-4, 7), 
                                 breaks=c(-3, 0, 3, 6)) +
              scale_y_continuous("Adjusted p-value (-log10)", limits=c(0, 30)) +
              facet_grid(Outcome~Group) +
              labs(title="DE mRNAs relative to baseline") +
              guides(alpha="none",
                     color=guide_legend(override.aes=list(size=3))) +
              theme(legend.position=c(0.9, 0.8))
# volcano legendplex
topgenes <- lp.res$time %>%
            filter(Significant) %>%
            group_by(Regulation, Group, Outcome) %>%
            top_n(n=5, wt=abs(lfc)) %>%
            ungroup()
fig.sup2$b <- lp.res$time %>%
              ggplot(aes(lfc, -log10(padj))) +
              geom_point(aes(col=Regulation, alpha=Significant, size=Regulation)) +
              scale_size_manual(values=c(Up=1, Down=1), na.value=0.5) +
              scale_alpha_manual(values=c("TRUE"=0.8, "FALSE"=0.5)) +
              scale_color_manual(values=regulation.cols, na.value="grey80") + 
              ggrepel::geom_text_repel(data=topgenes, aes(label=Gene), size=2) +
              scale_x_continuous("Fold change (log2)", limits=c(-4, 7), 
                                 breaks=c(-3, 0, 3, 6)) +
              scale_y_continuous("Adjusted p-value (-log10)", limits=c(0, 20)) +
              facet_grid(Outcome~Group) +
              labs(title="DE proteins relative to baseline") +
              guides(alpha="none",
                     color=guide_legend(override.aes=list(size=3))) +
              theme(legend.position=c(0.9, 0.8))

# outcome volcano: nanostring (late only)
topgenes <- ns.res$outcome %>%
            filter(Significant) %>%
            group_by(Group, Regulation) %>%
            top_n(n=8, wt=abs(lfc)) %>%
            ungroup()
fig.main$c <- ns.res$outcome %>%
              filter(Group %in% topgenes$Group) %>%
              ggplot(aes(lfc, -log10(padj))) +
              geom_point(aes(col=Regulation, alpha=Significant, size=Regulation)) +
              scale_size_manual("Higher expr. in", values=c(Survived=1, Succumbed=1), 
                                na.value=0.5) +
              scale_alpha_manual(values=c("TRUE"=0.8, "FALSE"=0.5)) +
              scale_color_manual("Higher expr. in", values=outcome.cols, 
                                 na.value="grey80") +
              ggrepel::geom_text_repel(data=topgenes, aes(label=Gene), size=3) +
              labs(x="Fold change (log2)",
                   y="Adjusted p-value (-log10)",
                   title="Late mRNA, survived vs. succumbed") +
              guides(alpha="none",
                     color=guide_legend(override.aes=list(size=3))) +
              theme(legend.position=c(0.4, 0.8))
# outcome volcano: legendplex (late only)
topgenes <- lp.res$outcome %>%
            filter(Significant) %>%
            group_by(Group, Regulation) %>%
            top_n(n=8, wt=abs(lfc)) %>%
            ungroup()
fig.main$d <- lp.res$outcome %>%
              filter(Group %in% topgenes$Group) %>%
              ggplot(aes(lfc, -log10(padj))) +
              geom_point(aes(col=Regulation, alpha=Significant, size=Regulation)) +
              scale_size_manual("Higher expr. in", values=c(Survived=1, Succumbed=1), 
                                na.value=0.5) +
              scale_alpha_manual(values=c("TRUE"=0.8, "FALSE"=0.5)) +
              scale_color_manual("Higher expr. in", values=outcome.cols, 
                                 na.value="grey80") +
              ggrepel::geom_text_repel(data=topgenes, aes(label=Gene), size=3) +
              labs(x="Fold change (log2)",
                   y="Adjusted p-value (-log10)",
                   title="Late protein, survived vs. succumbed") +
              guides(alpha="none",
                     color=guide_legend(override.aes=list(size=3))) +
              theme(legend.position=c(0.75, 0.7))
rm(topgenes)

## cibersort -------------------------------------------------------------------
cibx <- read.csv("analysis/cibersort-output.csv", check.names=FALSE) %>%
        rename(ID=Mixture) %>%
        filter(`P-value` < 0.05) %>%
        select(-`P-value`, -Correlation, -RMSE, -starts_with("Absolute score"),
               -`Naïve B`, -`Tfh`, -`Tregs`, -`Activated NK`) %>%
        reshape2::melt(id.vars="ID", 
                       variable.name="Celltype", 
                       value.name="Score") %>%
        left_join(meta, by="ID") %>%
        # remove outcome from naive samples
        mutate(Outcome=if_else(Stage=="Baseline", "Baseline", Outcome),
               Subgroup=interaction(Stage, Outcome))

# for each cell type, run KW test. Only run stats on cell types that pass
pval <- cibx %>%
        group_by(Celltype) %>%
        rstatix::kruskal_test(Score ~ Subgroup) %>%
        filter(p < 0.05) %>%
        select(Celltype) %>%
        # add back just the data from the filtered cell types
        left_join(cibx, by="Celltype") %>%
        group_by(Celltype) %>%
        # run Dunn's test for post-hoc
        rstatix::dunn_test(Score ~ Subgroup) %>%
        # keep only comparisons that are significant AND (1) comparison to 
        # baseline or (2) comparison of outcome at the same time
        mutate(Group.g1=str_extract(group1, "^[^\\.]+"),
               Group.g2=str_extract(group2, "^[^\\.]+"),
               Outcome.g1=str_extract(group1, "(?<=\\.).+$"),
               Outcome.g2=str_extract(group2, "(?<=\\.).+$")) %>%
        filter(p.adj < 0.05,
               (Group.g1=="Baseline" | Group.g1==Group.g2)) %>%
        # use group and outcome to define x position
        mutate(xmin=as.integer(factor(Group.g1, levels=names(stage.cols))),
               xmax=as.integer(factor(Group.g2, levels=names(stage.cols))),
               toggle.xmin=factor(Outcome.g1, 
                                  levels=c("Succumbed", "Survived", "Baseline"),
                                  labels=c(-0.2, 0.2, 0)),
               toggle.xmin=as.numeric(as.character(toggle.xmin)),
               toggle.xmax=factor(Outcome.g2, 
                                  levels=c("Succumbed", "Survived"),
                                  labels=c(-0.2, 0.2)),
               toggle.xmax=as.numeric(as.character(toggle.xmax)),
               xmin=xmin+toggle.xmin,
               xmax=xmax+toggle.xmax)

# plot each cell type with significant changes
cibx <- unique(pval$Celltype) %>%
        lapply(function(i) {
          # build the basic plot
          plt <- cibx %>%
            filter(Celltype==i) %>%
            ggplot(aes(Stage, Score)) +
            geom_boxplot(aes(group=Subgroup, col=Outcome, fill=Outcome), 
                         alpha=0.5, outlier.size=0.5, outlier.shape=21) 
          
          # subset the pvals and set our increment
          p <- filter(pval, Celltype==i)
          increment <- cibx[cibx$Celltype==i, "Score"] %>% 
            range() %>% 
            diff()
          increment <- 0.1 * increment
          
          # if there are outcome differences, add them
          if(nrow(filter(p, Group.g1 != "Baseline")) > 0) {
            comp <- p %>%
              filter(Group.g1 != "Baseline") %>%
              rstatix::add_y_position(scales="free_y", step.increase=0)
            plt <- plt + ggpubr::stat_pvalue_manual(data=comp, label="p.adj.signif", 
                                                    tip.length=0,
                                                    bracket.nudge.y=increment/2)
          }
          
          # if there are baseline differences, add them
          if("Baseline" %in% p$Group.g1) {
            comp <- filter(p, Group.g1=="Baseline")
            comp$y.position <- increment * 1:nrow(comp)
            comp$y.position <- comp$y.position + max(cibx[cibx$Celltype==i, "Score"])
            plt <- plt + ggpubr::stat_pvalue_manual(data=comp, label="p.adj.signif", 
                                                    tip.length=0, color="Outcome.g2")
          }
          
          # finish formatting the plot
          plt +
            scale_color_manual(values=outcome.cols, na.value="lightgrey") +
            scale_fill_manual(values=outcome.cols, na.value="lightgrey") +
            scale_y_continuous("Cell type score", expand=expansion(c(0.05, 0.1))) +
            labs(x=NULL, title=i) +
            theme(legend.position="none",
                  axis.ticks.y=element_blank(),
                  axis.text.y=element_blank(),
                  axis.text.x=element_text(angle=30, hjust=0.7, vjust=0.8))
        })
names(cibx) <- unique(pval$Celltype)

# plot all
cowplot::plot_grid(plotlist=cibx, labels="AUTO", ncol=3)

# clean up
rm(cibx, pval)

## selected expression ---------------------------------------------------------
# LEGENDplex: use log2-linearized, batch-corrected matrix
log.prot <- model.matrix(~0 + Stage, data=meta)
log.prot <- prot %>%
            log2() %>%
            removeBatchEffect(batch=meta$Batch, design=log.prot) %>%
            as.data.frame() %>%
            rownames_to_column("Analyte") %>%
            reshape2::melt(id.vars="Analyte",
                           variable.name="ID",
                           value.name="Concentration") %>%
            left_join(meta, by="ID")

# MPO
fig.main$e <- log.prot %>%
              filter(Analyte=="MPO") %>%
              ggplot(aes(Stage, Concentration)) +
              geom_boxplot(aes(group=interaction(Stage, Outcome), 
                               col=Outcome, fill=Outcome), 
                           alpha=0.5, outlier.size=1, outlier.shape=21) +
              scale_color_manual(NULL, values=outcome.cols) +
              scale_fill_manual(NULL, values=outcome.cols) +
              labs(x="Disease stage", y="Concentration (log2)", title="MPO") +
              theme(axis.text.x=element_text(angle=30, hjust=0.7, vjust=0.8),
                    legend.position=c(0.2, 0.8))

# MPO
fig.main$f <- log.prot %>%
              filter(Analyte=="IL-6") %>%
              ggplot(aes(Stage, Concentration)) +
              geom_boxplot(aes(group=interaction(Stage, Outcome), 
                               col=Outcome, fill=Outcome), 
                           alpha=0.5, outlier.size=1, outlier.shape=21) +
              scale_color_manual(NULL, values=outcome.cols) +
              scale_fill_manual(NULL, values=outcome.cols) +
              labs(x="Disease stage", y="Concentration (log2)", title="IL-6") +
              theme(axis.text.x=element_text(angle=30, hjust=0.7, vjust=0.8),
                    legend.position=c(0.2, 0.8))

# clean up
rm(log.prot)

## build and save figures ------------------------------------------------------
# main figure
cowplot::plot_grid(fig.main$a, fig.main$c, 
                   fig.main$b, fig.main$d, 
                   fig.main$e, fig.main$f, 
                   ncol=2, labels=c("a", "c", "b", "d", "e", "f"))
ggsave("analysis/figure-main.pdf", units="in", width=7.5, height=10)

# supplemental fig #1
cowplot::plot_grid(plotlist=fig.sup1, labels=names(fig.sup1), ncol=1)
ggsave("analysis/figure-sup1.pdf", units="in", width=7.5, height=6)

# supplemental fig #2
cowplot::plot_grid(plotlist=fig.sup2, labels=names(fig.sup2), ncol=1)
ggsave("analysis/figure-sup2.pdf", units="in", width=7.5, height=6)

## fin! ------------------------------------------------------------------------
sessionInfo()
