# LOAD LIBRARIES

library(tidyverse)
library(data.table); # data.table::setDTthreads(threads = 4)
library(RColorBrewer)
library(ggplot2)
library(edgeR)
library(patchwork)
library(scales)

options(scipen = 10)
options(rstudio.help.showDataPreview = FALSE)

# Working directory
setwd("~/HPC_data/ROHbin/")
# Read data & metadata
mtd = fread(file = "Metadata_shareJuan.csv")
mtd = mtd[order(mtd$Type),]

# Read HTZ data
reso = "250Kb" #"1Mb" # "100Kb" # "250Kb" # "1Mb" #  
reso_full = 250000

tt.raw = fread(file = paste0("Heterozygosity_window_",reso,".csv"),dec = ",")
tt0 = tt.raw[,-c(1:4)] # remove colums with chromosome info

# 1. FILTER DATA: ####

## 1.1 ALL CAPTIVE VS ALL WILD ####
ord = match(mtd$SampleID, colnames(tt0)) # Set the order of individuals
tt1 = tt0[,c(15,16,17,18,19,20,22,23,24,25,26,27,28,29,30,31,32,33,2,3 ,4 ,5 ,6 ,7 ,8 ,9 ,10,1,12,13,14,11,21,34,35,36,37)]
tt2 = cbind(tt.raw[,1:4],tt1)
# unique(colnames(tt1) == mtd$SampleID)

## 1.2 ALL WILD VS LATE CAPTIVE ####
# Filter for Late Captive vs. Wild
mtd = mtd[mtd$Category %in% c("Late_captive","Wild"),]
late_captive_wild = mtd$SampleID
tt1 = data.frame(tt0)[,c(late_captive_wild)]
tt2 = cbind(tt.raw[,1:4],tt1)
# colnames(tt1)
# unique(colnames(tt1) == mtd$SampleID)

# Chromosome list
chr_comp = unique(tt2$Scaffold)
# Entries per chromosome
sizes = tt2 %>% count(Scaffold)

# Functions: ####
# Plot distributions of EVs
distribEV1 <- function(tt = tt2){
  gg = reshape2::melt(tt[,-c(1,2,3,4)])
  gg$variable = factor(gg$variable)
  cnt_densities = ggplot(gg, aes(x = value, fill = variable)) +            # Draw two histograms in same plot
    geom_histogram(alpha = 0.5, position = "identity",bins=100) +
    geom_density(alpha = 0.3) +
    facet_wrap(~ variable) +
    theme_classic() + 
    guides(fill="none")
  return(cnt_densities)
}
# eBayes function
lmFit_ebayes_run <- function(tt.chr = tt.chr, design = NULL){
  fit0 = lmFit(tt.chr, design)
  fit = eBayes(fit0, trend = TRUE)
  # We can add the coefficient to just pick those contrasts affecting a certain sample. If we choose 1, equals the intercept.
  res.eb.full = topTable(fit, number = Inf, p.value = 1, adjust.method = "BH", confint = TRUE) # coef = 1
  # res.eb = res.eb.full[res.eb.full$adj.P.Val < pvaladj,]
  return(res.eb.full)
}
# Run it all the differential analysis for any condition, determined by the design matrix
differential_analysis_ebayes <- function(tt.chr = tt.chr, design = NULL, pval = pval){
  message("\nSelecting significant entries based on: ",colnames(design)[2]," | PVALUE < ",pval,"\n")
  # Run eBayes
  full.eb = lmFit_ebayes_run(tt.chr = tt.chr, design = design)
  # res.eb = full.eb[full.eb$adj.P.Val < pvaladj,]
  res.eb = full.eb[full.eb$P.Value < pval,]
  
  # How many variable bins between individuals?
  var.perc = round(dim(res.eb)[1]/dim(tt.chr)[1]*100, 2)
  message("\n",paste0("A total of ",var.perc,"% bins are significantly variable in terms of compartments in chr",chr.eb," WITH EFFECT SIZE > 0"))
  if (var.perc == 0){
    return(var.perc)
  } else {
    # Order by row, to see how continuous chunks are variable
    var.bins = sort(as.integer(rownames(res.eb)))
    
    # DBSCAN to get clusters of continuously variable compartments 
    bmag = dbscan_cluster(var.bins)
    return(list(bmag, var.bins, full.eb, var.perc))
  }
}
# Function for selecting the variable compartments in cluster.
select_variable_comps <- function(bmag = NULL, clust = NULL, var.bins = NULL,tt.chr = tt.chr){
  var.comps = which(bmag$cluster == clust)
  # Select components of that cluster.
  tr.comps.var = var.bins[var.comps]
  # Adding a few padding bins up and down, for context for variable compartments
  tr.comps = padding_compartments(tr.comps.var, 8, 0, max(as.numeric(rownames(tt.chr))))
  tt.reg = as.matrix(tt.chr[rownames(tt.chr) %in% tr.comps,])
  # colnames(tt.reg) = samps
  sel.comps = reshape2::melt(tt.reg)
  return(list(sel.comps,tr.comps.var,tr.comps))
}
# dbscan function
dbscan_cluster <- function(var.bins = var.bins, min_bins = 1){
  k = dist(var.bins, diag = FALSE, upper = FALSE);
  bmag = dbscan::dbscan(k, 1, min_bins); # 1 bin extra padding to get the continuous, and 3 minimum continuous bins to include
  nclust = max(bmag$cluster);
  return(bmag)
}
# Function for padding the variable bins
padding_compartments <- function(initial_series, consecutive, min_value, max_value) {
  extended_series <- seq(min(initial_series) - consecutive, max(initial_series) + consecutive)
  extended_series <- unique(pmax(extended_series, min_value))
  extended_series <- pmin(extended_series, max_value)
  return(extended_series)
}
# Plotting compartment variation, with the pvalue plot.
plot_view_variable <- function(sel.comps, tr.comps.var, tr.comps, eb.full = NULL, chr.eb = chr.eb, clust = clust, pval = pval, contrast, reso_full){
  zlim = max(abs(min(sel.comps$value)),abs(max(sel.comps$value)))
  sel.comps$xaxis = unique(tr.comps)[unique(tr.comps) != 0] * (reso_full/1000000)
  if(contrast == "LATE"){
    sepline = 6.5
    axis_y_name = "Late captive                                   Wild"
  } else {
    sepline = 18.5
    axis_y_name = "Captive                             Wild"
  }
  comp.pt = ggplot(sel.comps, aes(x = xaxis, y = Var2, fill = -log10(value))) +
    geom_tile() +
    scale_fill_gradient(low = "white", high = "firebrick") +
    labs(title = paste0("ROH change | In scaffold ", chr.eb, " | cluster: ", clust), 
         x = "Chromosome bin (Mb)", 
         y = axis_y_name) +
    labs(fill = "ROH\n-log10(Htz)") +
    geom_hline(yintercept = sepline, lty = "dashed", col="grey", lwd = 0.5) +
    geom_vline(xintercept = c((min(tr.comps.var)*(reso_full/1e6))-(reso_full/1000000)/2,max(tr.comps.var)*(reso_full/1e6)+(reso_full/1000000)/2), lty = 2, col="darkblue", lwd = 1) +
    theme_classic()
  # comp.pt
  # Adding pvalue line
  res.tmp.sorted = eb.full[order(as.numeric(row.names(eb.full))), ]
  res_ext = res.tmp.sorted[rownames(res.tmp.sorted) %in% tr.comps,] # $adj.P.Val
  pv.df = data.frame(as.numeric(rownames(res_ext)), res_ext$P.Value)
  colnames(pv.df) = c("bin","pv")
  pv.df$bin = pv.df$bin * (reso_full/1000000)
  pv.pt = ggplot(pv.df, aes(x = bin, y = -log10(pv))) + 
    geom_tile(aes(fill="black")) +
    scale_fill_manual(values="white") +
    geom_hline(yintercept=-log10(pval),col="red", lty=2)+
    geom_line(lty = 3) +
    geom_point() +
    guides(fill="none") +
    # geom_smooth() +
    theme_classic() 
  
  # Plot combined
  return(comp.pt / pv.pt)
}
# Generate plot automatically
generate_plot <- function(clust = cl, contrast = contrast){
  # Selecting the variable bins.
  sel.comps.list = select_variable_comps(bmag = bmag.chosen[[1]], clust = clust, var.bins = bmag.chosen[[2]], tt.chr) # Returns a list with 3 elements.
  sel.comps = sel.comps.list[[1]]
  tr.comps.var = sel.comps.list[[2]]
  tr.comps = sel.comps.list[[3]]
  # Order levels according to Day factor
  sel.comps$Var2 = factor(sel.comps$Var2, levels = mtd$SampleID)
  # Call plot
  roh_plot_pv = plot_view_variable(sel.comps, tr.comps.var, tr.comps, eb.full = bmag.chosen[[3]], chr.eb, clust, pval, contrast = contrast, reso_full)
  return(roh_plot_pv)
}

# 1. Distributions of values #####
htz_densities = distribEV1(tt = tt2)
print(htz_densities)

# 2. CONTRASTS #### 
# Design matrix type ALL captive - wild
abb.type = factor(mtd$Type, levels = unique(mtd$Type))
design.type = model.matrix(~abb.type)

# Design matrix type late captive - wild
abb.late = factor(mtd$Type, levels = unique(mtd$Type))
design.late = model.matrix(~abb.late)

# Design matrix ALL possible categories
abb.category = factor(mtd$Category, levels = unique(mtd$Category))
design.category = model.matrix(~abb.category)

# Design matrix combined
# sex35D = "Female"
# sex35DSpecific = factor(ifelse(mtd$Sex == sex35D & mtd$Day == "35",paste0(sex35D,"35"),"Other"))
# design.sex35D = model.matrix(~sex35DSpecific)

## 2.1. Choose contrast ####
contrast = "TYPE"; CHOSEN = design.type
contrast = "LATE"; CHOSEN = design.late
# contrast = "CATEGORY"; CHOSEN = design.category

# Order data frame for plotting
mtd.sort = mtd %>% arrange(desc(Type)) %>% dplyr::select(SampleID) %>% .$SampleID

# 3. RUN eBayes ####
pval = 0.05
# We can select the effect size change for our significant regions.
effect_size = 0 

# Run eBayes for all chromosomes, while also store the direction of change
ebayes_signif = list()
# To store the variable percentages
ebayes_var_perc = c()

for (chr.eb in chr_comp){
  # chr.eb = "NC_072852.1"
  tt = tt2[tt2$Scaffold == chr.eb,][,-c(1:4)]
  tt.chr = tt[complete.cases(tt),]
  # Running eBayes and selecting significant clusters
  bmag.chosen = differential_analysis_ebayes(tt.chr = tt.chr, design = CHOSEN, pval = pval); # bmag.chosen[[3]]  
  
  # Save variable compartments and the effect direction; otherwise save NA.
  if ( class(bmag.chosen) != "list" ) {
    
    # message("No significant bins computed for: chr",chr.eb)
    ebayes_signif[["posFC"]][[chr.eb]] = NA
    ebayes_signif[["negFC"]][[chr.eb]] = NA
    ret.perc = 0
  } else {
    
    # Select independently towards negative (B) and positive (A) switching compartments
    posFC.all = as.numeric(rownames(bmag.chosen[[3]][bmag.chosen[[3]]$logFC > effect_size,]))
    negFC.all = as.numeric(rownames(bmag.chosen[[3]][bmag.chosen[[3]]$logFC < -effect_size,]))
    posFC = bmag.chosen[[2]][bmag.chosen[[2]] %in% posFC.all]
    negFC = bmag.chosen[[2]][bmag.chosen[[2]] %in% negFC.all]
    ebayes_signif[["posFC"]][[chr.eb]] = posFC # More ROH in Captive
    ebayes_signif[["negFC"]][[chr.eb]] = negFC # More ROH in Wild
    ret.perc = bmag.chosen[[4]]
  }
  # Add the variable compartments % per chromosome
  ebayes_var_perc = c(ebayes_var_perc, ret.perc)
  
  if (class(bmag.chosen) != "list"){
    next
  } else {
    # Save all the plots to file
    for (cl in unique(bmag.chosen[[1]]$cluster)){
      # final_plot = generate_plot(clust = cl, contrast = contrast)
      # ggsave(filename = paste0("~/Projects/ROHbin/plots/","reso_",reso,"/",chr.eb,"_clust_",cl,"_contrast_",contrast,"_",reso,".pdf"),final_plot)
      }
  }
  }

# Summarize results
results = cbind(chr_comp,sizes,ebayes_var_perc)
results$bins_var = round((results$ebayes_var_perc * results$n)/100 )
fwrite(x = results, file = paste0("~/HPC_data/ROHbin/results_ROH_per_chrom", effect_size, ".pv_", scientific(pval), "_", contrast,"_",reso,".txt"),sep = "\t", row.names = F)
# results = fread(file = paste0("~/HPC_data/ROHbin/results_ROH_per_chrom", effect_size, ".pv_", scientific(pval), "_", contrast,"_",reso,".txt"))
# % of the genome is ROH:
roh_perc = sum(results$bins_var)/sum(results$n)*100
print(roh_perc)
# Mb totales in ROH
(sum(results$n)*reso_full/100) * roh_perc

# Save object in file
save(file = paste0("~/HPC_data/ROHbin/ebayes_signif.eff_", effect_size, ".pv_", scientific(pval), "_", contrast,"_",reso), ebayes_signif)

# ROH in captive
ROH_Captive = unname(unlist(ebayes_signif[1])) # positive LOG_FC
((length(ROH_Captive[!is.na(ROH_Captive)]) * reso_full) / 1112400596) * 100 # Het en Wild

# ROH in wild
ROH_Wild = unname(unlist(ebayes_signif[2])) # negative LOG_FC
((length(ROH_Wild[!is.na(ROH_Wild)]) * reso_full) / 1112400596) * 100 # Het en Captive

stop("THE END.")

# Check values by hand.
tt = tt2[tt2$Scaffold == chr.eb,]
tt[c(28:33),]
tt$EB8_S119


