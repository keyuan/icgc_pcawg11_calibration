morris_file = "data/morris/purity_ploidy.txt"
vanloo_wedge_file = "data/vanloo_wedge/1_purity_ploidy/purity_ploidy.txt"
peifer_file = "data/peifer/Purity_Ploidy/ICGC_Pilot63_Purity_Ploidy.txt"
sahinalp_file = "data/sahinalp/citup_pilot63/citup_purity_ploidy_v1.txt"
broad_file = "data/broad/titan_purity_ploidy.txt"

source("code/Parser.R")
source("code/Plotting.R")

morris = parse.purity.ploidy(morris_file)
vanloo_wedge = parse.purity.ploidy(vanloo_wedge_file)
peifer = parse.purity.ploidy.peifer(peifer_file, vanloo_wedge$sample)
sahinalp = parse.purity.ploidy(sahinalp_file)
# Replace 0 with NA
sahinalp$ploidy = rep(NA, nrow(sahinalp))
broad = parse.purity.ploidy(broad_file)
colnames(broad) = c("sample", "purity", "ploidy")
list_of_tables = list(morris, vanloo_wedge, peifer, sahinalp, broad)
vector_of_names = c("morris", "vanloo_wedge", "peifer", "sahinalp", "broad")

#' Convert a list of tables into a matrix with purity estimates only, placing NAs where a group has not reported anything
create.purity.table = function(list_of_tables, vector_of_names) {
  purity = list_of_tables[[1]][,c("sample", "purity")]
  for (i in 2:length(list_of_tables)) {
    # Add a new column for this purity table
    new.col = list_of_tables[[i]]$purity[match(purity$sample, list_of_tables[[i]]$sample)]
    purity = cbind(purity, new.col)
    
    # Create and add rows for samples not previously mentioned
    new.index = which(!(list_of_tables[[i]]$sample %in% purity$sample))
    new.rows = matrix(NA, ncol=ncol(purity), nrow=length(new.index))
    new.rows[,1] = list_of_tables[[i]]$sample[new.index]
    new.rows[,i+1] = list_of_tables[[i]]$purity[new.index]
    colnames(new.rows) = colnames(purity)
    purity = rbind(purity, new.rows)
  }
  colnames(purity) = c("sample", vector_of_names)
  purity = purity[order(purity$sample),]
  row.names(purity) = 1:nrow(purity)
  return(purity)
}

#' Create a ploidy overview table, inserting NAs where ploidy was not reported
create.ploidy.table = function(list_of_tables, vector_of_names) {
  ploidy = list_of_tables[[1]][,c("sample", "ploidy")]
  for (i in 2:length(list_of_tables)) {
    # Add a new column for this ploidy table
    new.col = list_of_tables[[i]]$ploidy[match(ploidy$sample, list_of_tables[[i]]$sample)]
    ploidy = cbind(ploidy, new.col)
    
    # Create and add rows for samples not previously mentioned
    new.index = which(!(list_of_tables[[i]]$sample %in% ploidy$sample))
    new.rows = matrix(NA, ncol=ncol(ploidy), nrow=length(new.index))
    new.rows[,1] = list_of_tables[[i]]$sample[new.index]
    new.rows[,i+1] = list_of_tables[[i]]$ploidy[new.index]
    colnames(new.rows) = colnames(ploidy)
    ploidy = rbind(ploidy, new.rows)
  }
  colnames(ploidy) = c("sample", vector_of_names)
  ploidy = ploidy[order(ploidy$sample),]
  row.names(ploidy) = 1:nrow(ploidy)
  return(ploidy)
}

#######################################################################
# Create tables
#######################################################################
purity = create.purity.table(list_of_tables, vector_of_names)
write.table(purity, "1_purity_ploidy/purity.tsv", sep="\t", quote=F, row.names=F)
ploidy = create.ploidy.table(list_of_tables, vector_of_names)
write.table(ploidy, "1_purity_ploidy/ploidy.tsv", sep="\t", quote=F, row.names=F)

#######################################################################
# Create heatmaps
#######################################################################
plot.purity.heatmap(purity, "1_purity_ploidy/purity.png", "Purity")
plot.purity.heatmap(ploidy, "1_purity_ploidy/ploidy.png", "Ploidy")

#######################################################################
# Create heatmap ploidy without outliers
#######################################################################
ploidy.no.outliers = ploidy[,2:ncol(ploidy)]
ploidy.no.outliers = apply(ploidy.no.outliers, 2, as.numeric)
ploidy.no.outliers[ploidy.no.outliers > 4.5] = NA
ploidy.no.outliers = data.frame(ploidy$sample, ploidy.no.outliers)
colnames(ploidy.no.outliers)[1] = "sample"
plot.purity.heatmap(ploidy.no.outliers, "1_purity_ploidy/ploidy_no_outliers.png", "Ploidy")
