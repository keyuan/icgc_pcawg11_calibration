####################################################################
# 1_purity_ploidy
####################################################################
#' Read in the purity/ploidy table
parse.purity.ploidy = function(infile) {
  return(read.table(infile, header=T, stringsAsFactors=F))
}

parse.purity.ploidy.peifer = function(infile, vector_of_samplenames) {
  d = read.table(infile, header=T, stringsAsFactors=F)
  for (i in 1:nrow(d)) {
    sampleid = unlist(strsplit(d$Sample[i], "_"))[2]
    
    selection = grepl(sampleid, vector_of_samplenames)
    if (sum(selection) != 1) {
      if (sampleid == "31f02f48") {
        d[i,1] = "31f02f48-44a4-445e-ac3d-e9bf3d8d25a2"
      } else if (sampleid == "f9c39eb7") {
        d[i,1] = "f9c39eb7-39a9-6626-e040-11ac0d4870c2"
      } else {
        print(paste("Found previously unknown sample that is not captured here:", sampleid))
      }
    } else {
      d[i,1] = purity$sample[grepl(sampleid, purity$sample)]
    }
  }
  colnames(d) = c("sample", "purity", "ploidy")
  return(d)
}

####################################################################
# 2_mutation_assignments
####################################################################
#' Read in a DP input file
parse.dp.input = function(infile) {
  d = read.table(infile, header=T)
  d = d[,c("chr", "end", "subclonal.fraction")]
  colnames(d) = c("Chromosome", "Position", "Subclonal.fraction")
  return(d)
}

#' Parse a mutation assignment file
parse.mut.assignments = function(infile) {
  return(read.table(infile, header=T, stringsAsFactors=F))
}

####################################################################
# 3_tree_structures
####################################################################
library(Matrix)
#' Parses the four tree mutation assignment files together with a dp_input file to construct four sparse matrices
parse.tree.structure = function(dp_input, mutation_names, identity, ancestor_child, child_ancestor, sibling) {
  #index = read.table(dp_input, header=T, stringsAsFactors=F)
  #no.muts = nrow(index)
  mutation_names = read.table(mutation_names, header=T, stringsAsFactors=F)
  no.muts = max(mutation_names$row)
  
  # TODO: sync the mutation names to the index? Or keep only the ones that are assigned?

  ident = sparse.input.to.Matrix(identity, no.muts)
  anc_child = sparse.input.to.Matrix(ancestor_child, no.muts)
  child_anc = sparse.input.to.Matrix(child_ancestor, no.muts)
  sibling = sparse.input.to.Matrix(sibling, no.muts)
  
  return(list(mut_names=mutation_names, ident=ident, anc_child=anc_child, child_anc=child_anc, sibling=sibling))
}

sparse.input.to.Matrix = function(infile, no.muts) {
  data = read.table(infile, header=T, stringsAsFactors=F)
  return(sparseMatrix(data$row, data$col, x=data$value, dims=c(no.muts, no.muts)))
}