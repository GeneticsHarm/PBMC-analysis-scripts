###################################
##### Cell-Cycle List Builder #####
###################################
#
# Author - Martin Barron
#
# Date - 23/01/2016
#

###############################
### Part 1 - Load the Data ####
###############################

load("~/HScc_genes.rda")
load("~/MMcc_genes.rda")

#######################################
### Part 2 - Identifying Gene Names ###
#######################################

# Type 1 - Homo Sapiens 
# Type 2 - Mus Musculus
# name_type - 1 for ensembl, 2 for symbol, 3 for entrez, 4 for Unigene
# No unigene for mouse. 


gene_indexer <- function(gene_names, type, name_type){
  cc_gene_indices <- rep(NA, 1)
  if (type == 1){
    if (name_type == 1){
      cc_gene_indices <- na.omit(match(human_cell_cycle_genes[,1], gene_names))
    } else if (name_type == 2){
      cc_gene_indices <- na.omit(match(human_cell_cycle_genes[,2], gene_names))
    } else if (name_type == 3){
      cc_gene_indices <- na.omit(match(human_cell_cycle_genes[,3], gene_names))
    } else if (name_type == 4){
      cc_gene_indices <- na.omit(match(human_cell_cycle_genes[,4], gene_names))
    }
  } else if (type == 2) {
    if (name_type == 1){
      cc_gene_indices <- na.omit(match(mouse_cell_cycle_genes[,1], gene_names))
    } else if (name_type == 2){ 
      cc_gene_indices <- na.omit(match(mouse_cell_cycle_genes[,2], gene_names))
    } else if (name_type == 3) {
      cc_gene_indices <- na.omit(match(mouse_cell_cycle_genes[,3], gene_names))
    } else if (name_type == 4) {
      cc_gene_indices <- na.omit(match(mouse_cell_cycle_genes[,4], gene_names))
    }
  }
  return(cc_gene_indices)
}

##### Example #####

gene_names <- mouse_cell_cycle_genes[,4]
type <- 2
name_type <- 4

cc_gene_indices <- gene_indexer(gene_names, type, name_type)

summary(cc_gene_indices)
length(unique(cc_gene_indices))
length(unique(gene_names))
