library(rbioapi)
library(biomaRt)
library(GenomicRanges)
library(dplyr)
library(scales)

options(scipen = 10)
options(rstudio.help.showDataPreview = FALSE)

# Read and calculate enrichment
# Cross with genes in regions
## CUSTOM FUNCTION TO INTERSECT BEDTOOLS-LIKE REGIONS.
bedTools.2in = function(functionstring = "intersectBed",bed1,bed2,opt.string=""){
  #create temp files
  a.file=tempfile()
  b.file=tempfile()
  out   =tempfile()
  options(scipen =99) # not to use scientific notation when writing out
  
  #write bed formatted dataframes to tempfile
  write.table(bed1,file=a.file,quote=F,sep="\t",col.names=F,row.names=F)
  write.table(bed2,file=b.file,quote=F,sep="\t",col.names=F,row.names=F)
  
  # create the command string and call the command using system()
  command=paste(functionstring,"-a",a.file,"-b",b.file,opt.string,">",out,sep=" ")
  cat(command,"\n")
  try(system(command))
  
  res=read.table(out,header=F)
  unlink(a.file);unlink(b.file);unlink(out)
  return(res)
}
setwd("~/HPC_data/ROHbin/")
tss_genes = fread("genome_annotated.bed")
colnames(tss_genes) = c("chromosome_name","start_position","end_position","gene_name")
# head(tss_genes)

# 1. eBAYES SIGNIFICANT REGIONS FOR ENRICHMENT ####

## 1.1 Load data for the desired contrast. ####

contrast = "TYPE";
# contrast = "LATE"
reso = "250Kb" # " 1Mb" # "250Kb" 
reso_full = 250000 # 250000
overlap = 5000/reso_full
pvaladj = 0.05;
effect_size = 0;
# Load the object generated earlier with the significant bins at window size X.
load(file = paste0("~/HPC_data/ROHbin/ebayes_signif.eff_", effect_size, ".pv_", scientific(pvaladj), "_", contrast,"_",reso)) # ebayes_signif

## 1.2 Format coordinates ####
# This function adds the start and end coordinates
add_coords <- function(bins){
  colnames(bins) = c("chr","bin")
  bins$start = (bins$bin * reso_full) - reso_full
  bins$end = bins$start + reso_full + 1
  bins = bins[ , !names(bins) %in% c("bin")]
  return(bins)
}

# Adds coordinates to each of the compartments datasets
FC.bins = list()

# ROH IN CAPTIVE
posFC.bins = add_coords(bins = purrr::map_df(ebayes_signif[["posFC"]], ~as.data.frame(.x), .id="id"))
posFC.bins = posFC.bins[complete.cases(posFC.bins),]
FC.bins[["posFC"]] = posFC.bins

# ROH IN WILD
negFC.bins = add_coords(purrr::map_df(ebayes_signif[["negFC"]], ~as.data.frame(.x), .id="id"))
negFC.bins = negFC.bins[complete.cases(negFC.bins),]
FC.bins[["negFC"]] = negFC.bins

## 2.3 Overlap and enrichment genes ####
enriched.eb = list()

for (dd in names(FC.bins)){
  # dd = "posFC"
  FC = FC.bins[[dd]]
  # View(FC)
  eb.overlap = bedTools.2in("intersectBed", FC, tss_genes, opt.string = paste0("-wa -wb -f ",overlap," -F 1")) #-f 0.01 -F 0.1
  colnames(eb.overlap) = c("chrA","start_comp","end_comp","chrB","start_tss","end_tss","gene_name")
  # length(unique(eb.overlap$gene_name))
  eb.genes = unique(eb.overlap$gene_name[!is.na(eb.overlap$gene_name)])
  print(length(eb.genes))
  write.table(eb.overlap, file = paste0(dd, "_effect_size_",effect_size,"_contrast_",contrast,"_",reso,"_ebayes.txt"), quote = FALSE,sep = '\t', row.names = FALSE)
  
  reactome = rba_reactome_analysis(input = eb.genes, projection = TRUE, p_value = 0.05)
  enriched.eb[[dd]][["reactome"]] = reactome
  
  # enrichr = rba_enrichr(gene_list = eb.genes)
  # enriched.eb[[dd]][["enrichr"]] = enrichr
  
  # PANTHER DOES NOT WORK #
  # RUN MANUALLY ON THE WEB #
  
  # for (db in dbs.panther){
  #  db ="GO:0008150" # "ANNOT_TYPE_ID_PANTHER_PATHWAY" # "GO:0008150" # "ANNOT_TYPE_ID_REACTOME_PATHWAY"
  #  query.eb = rba_panther_enrich(genes = eb.genes, correction = "FDR", cutoff = 0.1, organism = 9031, annot_dataset = db) # 9606 (human) 9031 (chicken) # FDR 0.1 is the choice
  #  enriched.eb[[db]][[dd]] = query.eb$result
  # }
}


# 2.4 Enrichplot if something significant?:  ####
library(enrichplot)
# >> https://bioconductor.org/packages/release/bioc/html/enrichplot.html
# >> https://alexslemonade.github.io/refinebio-examples/03-rnaseq/pathway-analysis_rnaseq_01_ora.html#49_Run_ORA_using_the_enricher()_function

stop("END HERE")

# TRASH ####
# Genes associated with Rho GTPases:
# rhog = c("ANKLE2","DLC1","ABI2","PCDH7","ACTR3","SCAI","CPNE8","VAPB","PLEKHG1","GOLGA3","DOCK10","CENPP","DOCK1","PREX2","VAV2","BCR","PAK5","COPS2")
# eb.overlap[eb.overlap$chicken_gene %in% rhog,]
# 
# write.table(file = "~/Desktop/test111.txt",x = query.eb$input_list$mapped_id,quote = F,row.names = F)
# View(eb.overlap[order(eb.overlap$human_gene),])

# for (ch in chr_comp){
#   ch = "20"
#   x.coords = data.frame(unlist(lapply(ebayes_signif,`[`,ch)))
#   colnames(x.coords) = c("bin")
#   x.coords$chr = ch
#   x.coords$start = (x.coords$bin * reso) - reso
#   x.coords$end = x.coords$start + reso + 1
#   x.coords$effect = unlist(strsplit(rownames(x.coords),'\\.'))[c(TRUE,FALSE)]
# }