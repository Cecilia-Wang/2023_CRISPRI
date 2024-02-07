
# Prepare the environment

### load packages + install ones that aren't  already installed in the libraries

packs_required<-c("BiocManager","pacman")
if_insta_packs<-packs_required %in% rownames(installed.packages())
if (any(if_insta_packs == FALSE)) {
  install.packages(packs_required[!if_insta_packs])
}

# load CRAN and Bioconductor packages
pacman::p_load(tidyverse,data.table,ggrepel,ggpubr,DESeq2,viridis,RColorBrewer,ggupset,ggtext,ggExtra,rstatix,patchwork)

# load github packages where needed 
# pacman::p_load_gh("coolbutuseless/ggpattern", "davidsjoberg/ggsankey") # for example

### load exact test results and metadata

# Load exact results files without removing low abundant guides
files_path<-"~/Downloads/2023_CRISPRI-main/data/" # Modify this path to your path if necessary

# Optional, use this when loading exact results from multiple runs or projects
ex_res_all_multirun<-NULL


ex_file_all<-list.files(paste0(files_path,"exact_test"),pattern="tsv")
ex_res_all<-lapply(paste0(files_path,"exact_test/",ex_file_all), read.delim, header = TRUE,sep = "")
ex_res_all <- do.call("rbind", ex_res_all)

## add experiment id to differentiate results from different experiment or projects
ex_res_all$exp.id<-"ex03" # customise this to suit individual project
ex_res_all_multirun<-rbind(ex_res_all_multirun,ex_res_all)

# Load metadata
meta_group<-read.csv(paste0(files_path, "meta_group.csv"))

sum_guide<-read.csv(paste0(files_path,"summary_guide_data.csv"))

# Information regarding essential calls from Bosch 2021
rock_essM<-read.csv(paste0(files_path,"Mtb_H37rv_essentiality.csv"))

PATRIC_functions<-read.csv(paste0(files_path,"PATRIC_function_classification_customised.csv"))
