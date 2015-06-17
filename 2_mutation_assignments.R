args = commandArgs(TRUE)

# samplename = "569393c8-e2fe-4580-a45b-81f1b1e01135"
# morris_file = paste("data/morris/mutation_assignment/mutation_assignment.all.", samplename, ".csv", sep="")
# vanloo_wedge_file = paste("data/vanloo_wedge/2_clustering/", samplename, "/", samplename, "_cluster_membership.txt", sep="")
# peifer_file_assignments = "data/peifer/Mutation_Clustering/ESAD_293a2f0a_cluster_assignments.txt"
# peifer_file_clusters = "data/peifer/Mutation_Clustering/ESAD_293a2f0a_mclusters.txt"
# sahinalp_file = paste("data/sahinalp/citup_pilot63/samples_v1/", samplename, "_cluster_membership.txt", sep="")

morris_file = toString(args[1])
vanloo_wedge_file = toString(args[2])
peifer_file_assignments = toString(args[3])
peifer_file_clusters = toString(args[4])
sahinalp_file = toString(args[5])

samplename = unlist(strsplit(vanloo_wedge_file, "/"))
samplename = unlist(strsplit(samplename[length(samplename)], "_"))[1]

source("code/Parser.R")
source("code/Plotting.R")
source("code/ArrayAnalysis.R")

if (file.exists(morris_file)) {
  morris = parse.mut.assignments(morris_file)
} else {
  morris = NULL
}

if (file.exists(vanloo_wedge_file)) {
  vanloo_wedge = parse.mut.assignments(vanloo_wedge_file)
} else {
  vanloo_wedge = NULL
}

if (file.exists(peifer_file_assignments)) {
  peifer = parse.mut.assignments.peifer(peifer_file_assignments, peifer_file_clusters)
} else {
  peifer = NULL
}

if (file.exists(sahinalp_file)) {
  sahinalp = parse.mut.assignments(sahinalp_file)
} else {
  sahinalp = NULL
}

vector_of_names = c("morris", "vanloo_wedge", "peifer", "sahinalp")
list_of_tables = list(morris, vanloo_wedge, peifer, sahinalp)

# vector_of_names = c("vanloo_wedge", "peifer")
# list_of_tables = list(vanloo_wedge, peifer)

#######################################################################
# Get which mutations were assigned by whome
#######################################################################
get.mut.assigned.by.who = function(master_table, list_of_tables, vector_of_names) {
  for (i in 1:length(list_of_tables)) {
    #is_assigned = array(F, nrow(master_table))
    sub = list_of_tables[[i]]
    if (!is.null(sub)) {
      num.cols = ncol(sub)
      # Subset such that all mentioned mutations are assigned
      is_assigned = sapply(1:nrow(sub), function(i, sub, num.cols) { !any(is.na(sub[i,3:num.cols])) }, sub=sub, num.cols=num.cols)
      sub = sub[is_assigned,]
      master_table = cbind(master_table, (master_table[,1] %in% sub[,1] & master_table[,2] %in% sub[,2]))
    } else {
      master_table = cbind(master_table, rep(NA, nrow(master_table)))
    }
  }
  colnames(master_table) = c("Chromosome", "Position", vector_of_names)
  return(master_table)
}

# Use this as the master table, the next step will only report mutations mentioned in this table
if(!is.null(vanloo_wedge)) {
  master_table = vanloo_wedge[,c(1,2)]
} else {
  master_table = peifer[,c(1,2)]
}
assignment_overview = get.mut.assigned.by.who(master_table, list_of_tables, vector_of_names)
write.table(assignment_overview, paste("2_mutation_assignments/tables/", samplename, "_mutation_assignment_inventory.tsv", sep=""), sep="\t", quote=F, row.names=F)

#######################################################################
# Synchronise mutation assignment tables
#######################################################################
#' Tests which mutations were assigned by all methods.
#' @return A vector with a boolean for each mutation.
mutation.assigned.by.all = function(assignment.tables, raw.data, vector_of_names) {
  sel = rep(T, nrow(raw.data))
  for (i in 1:length(vector_of_names)) {
    # Don't take into account empty assignment tables, which occurs when a method hasn't pruduced output for a sample
    if (!is.null(assignment.tables[[i]]) && nrow(assignment.tables[[i]]) > 0) {
	chrpos = paste(raw.data$Chromosome, raw.data$Position)
        chrpos_method = paste(assignment.tables[[i]][,1], assignment.tables[[i]][,2])
	sel = sel & chrpos %in% chrpos_method
      #sel = sel & (raw.data$Chromosome %in% assignment.tables[[i]][,1] & raw.data$Position %in% assignment.tables[[i]][,2])
    }
  }
  return(sel)
}

#' Synchronise the assignments based on whether a mutation is assigned by all methods. Also annotates Subclonal.fraction to the df.
sync.assignments = function(assignment.tables, raw.data, selection, vector_of_names) {
  dat.shared = list()
  for (i in 1:length(vector_of_names)) {
    dat = assignment.tables[[i]]
    if (is.null(dat) || nrow(dat) == 0) {
      # This happens when a method hasn't produced results for a sample. Reuse the empty data.frame
      dat.shared[[i]] = data.frame()
    } else {
	chrpos_method = paste(dat[,1], dat[,2])
    	chrpos = paste(raw.data$Chromosome[selection], raw.data$Position[selection])
	dat.sel = chrpos_method %in% chrpos
        dat.shared[[i]] = dat[dat.sel,]
        dat.shared[[i]]$Subclonal.fraction = raw.data$Subclonal.fraction[selection]
    }
  }
  return(dat.shared)
}

#' Calculate the identity matrix for all methods
calc.ident.matrices = function(dat.shared, vector_of_names, useprobs=T) {
  ident.matrices = list()
  for (i in 1:length(vector_of_names)) {
    if (is.null(dat.shared[[i]]) || nrow(dat.shared[[i]]) == 0) {
      # If there are no mutation assignments available, set the identity matrix to NA
      ident.matrices[[i]] = data.frame()
    } else {
      #ident.matrices[[i]] = GetIdentityArrayFromAssignments(dat.shared[[i]]$Cluster)
      if (useprobs) {
        # Select the columns from dat.shared (no chr/pos and CCF)
        ident.matrices[[i]] = GetIdentityArrayFromProbabilities(as.matrix(dat.shared[[i]][,3:(ncol(dat.shared[[i]])-1)]))
      } else {
	cols_select = 3:(ncol(dat.shared[[i]])-1)
	if (length(cols_select) == 1) {
		most.likely.node.assignments = rep(1, nrow(dat.shared[[i]]))
	} else {
        	most.likely.node.assignments = unlist(apply(dat.shared[[i]][,cols_select], 1, which.max))
	}
        ident.matrices[[i]] = GetIdentityArrayFromAssignments(most.likely.node.assignments)
      }
    }
  }
  return(ident.matrices)
}

d = parse.dp.input(paste("data/dirichlet_input/", samplename, "_allDirichletProcessInfo.txt", sep=""))

# Sync the mutation calls between tables
# Only take mutations that are used by all methods
sel = mutation.assigned.by.all(list_of_tables, d, vector_of_names)

# Subset the assignment tables using the above selection to get the shared mutations
dat.shared = sync.assignments(list_of_tables, d, sel, vector_of_names)

# Calculate the identity matrix
ident.matrices.probs = calc.ident.matrices(dat.shared, vector_of_names, useprobs=T)
ident.matrices.noprobs = calc.ident.matrices(dat.shared, vector_of_names, useprobs=F)

#######################################################################
# Calculate similarities between each pair of matrices
#######################################################################
# Probs
diff.mse = matrix(0, length(vector_of_names), length(vector_of_names))
for (i in 1:length(vector_of_names)) {
  for (j in i:length(vector_of_names)) {
    if (nrow(ident.matrices.probs[[i]]) > 0 & nrow(ident.matrices.probs[[j]]) > 0) {
      diff.mse[i, j] = mean((ident.matrices.probs[[i]] - ident.matrices.probs[[j]])^2)
    } else {
      diff.mse[i, j] = NA
    }
  }
}
diff.mse.d = data.frame(diff.mse)
colnames(diff.mse.d) = vector_of_names
write.table(diff.mse.d, file=paste("2_mutation_assignments/similarities/", samplename, "_probs.mse.txt", sep=""), sep="\t", row.names=F, quote=F)

# No probs
diff.mse = matrix(0, length(vector_of_names), length(vector_of_names))
for (i in 1:length(vector_of_names)) {
  for (j in i:length(vector_of_names)) {
    if (nrow(ident.matrices.noprobs[[i]]) > 0 & nrow(ident.matrices.noprobs[[j]]) > 0) {
      diff.mse[i, j] = mean((ident.matrices.noprobs[[i]] - ident.matrices.noprobs[[j]])^2)
    } else {
      diff.mse[i, j] = NA
    }
  }
}
diff.mse.d = data.frame(diff.mse)
colnames(diff.mse.d) = vector_of_names
write.table(diff.mse.d, file=paste("2_mutation_assignments/similarities/", samplename, "_noprobs.mse.txt", sep=""), sep="\t", row.names=F, quote=F)

#######################################################################
# Plot a heatmap with data ordered according to the first method, if no results use the second
#######################################################################
if (!is.null(dat.shared[[1]]) & nrow(dat.shared[[1]]) > 0) {
  cols_select = 3:(ncol(dat.shared[[1]])-1)
  if (length(cols_select) == 1) { 
	  most.likely.node.assignments = rep(1, nrow(dat.shared[[1]]))
  } else {
  	most.likely.node.assignments = apply(dat.shared[[1]][,3:(ncol(dat.shared[[1]])-1)], 1, which.max)
  }
} else if(!is.null(dat.shared[[2]]) & nrow(dat.shared[[2]]) > 0) {
  cols_select = 3:(ncol(dat.shared[[2]])-1)
  if (length(cols_select) == 1) {
	  most.likely.node.assignments = rep(1, nrow(dat.shared[[2]]))
  } else {
	most.likely.node.assignments = apply(dat.shared[[2]][,3:(ncol(dat.shared[[2]])-1)], 1, which.max)
  }
} else {
  cols_select = 3:(ncol(dat.shared[[3]])-1)
  if (length(cols_select) == 1) {
	  most.likely.node.assignments = rep(1, nrow(dat.shared[[3]]))
  } else {
	  most.likely.node.assignments = apply(dat.shared[[3]][,3:(ncol(dat.shared[[3]])-1)], 1, which.max)
  }
}

ord = order(most.likely.node.assignments)


# No probabilities
for (i in 1:length(vector_of_names)) {
  # Only plot heatmap when there is a matrix, i.e. when a method has produced mutation assignments for this sample
  if (!is.null(ident.matrices.noprobs[[i]]) & nrow(ident.matrices.noprobs[[i]]) > 0) {
    if (nrow(ident.matrices.noprobs[[i]]) < 10000) {
      plotHeatmapFull(ident.matrices.noprobs[[i]][ord,ord], paste("2_mutation_assignments/figures/", samplename, "_", vector_of_names[i], "_noprobs.png", sep=""))
    } else {
      plotHeatmapMedium(ident.matrices.noprobs[[i]][ord,ord], paste("2_mutation_assignments/figures/", samplename, "_", vector_of_names[i], "_noprobs.png", sep=""))
    }
  }
}

# Probabilities
for (i in 1:length(vector_of_names)) {
  # Only plot heatmap when there is a matrix, i.e. when a method has produced mutation assignments for this sample
  if (!is.null(ident.matrices.probs[[i]]) & nrow(ident.matrices.probs[[i]]) > 0) {
    if (nrow(ident.matrices.probs[[i]]) < 10000) {
      plotHeatmapFull(ident.matrices.probs[[i]][ord,ord], paste("2_mutation_assignments/figures/", samplename, "_", vector_of_names[i], "_probs.png", sep=""))
    } else {
      plotHeatmapMedium(ident.matrices.probs[[i]][ord,ord], paste("2_mutation_assignments/figures/", samplename, "_", vector_of_names[i], "_probs.png", sep=""))
    }
  }
}
#######################################################################
# Save a table with all most likely assignments
#######################################################################
# If first method hasn't produced assignments, use the second
if (!is.null(dat.shared[[1]]) & nrow(dat.shared[[1]]) > 0) {
	startpoint = 1
	most.likely.node.assignments = cbind(dat.shared[[startpoint]][,c(1,2)], most.likely.node.assignments)
} else if (!is.null(dat.shared[[2]]) & nrow(dat.shared[[2]]) > 0) {
	startpoint = 2
	most.likely.node.assignments = cbind(rep(NA, nrow(dat.shared[[startpoint]])), dat.shared[[startpoint]][,c(1,2)], most.likely.node.assignments)
} else {
	startpoint = 3
	most.likely.node.assignments = cbind(rep(NA, nrow(dat.shared[[startpoint]])), rep(NA, nrow(dat.shared[[startpoint]])), dat.shared[[startpoint]][,c(1,2)], most.likely.node.assignments)
}

#most.likely.node.assignments = cbind(dat.shared[[startpoint]][,c(1,2)], most.likely.node.assignments)
for (i in (startpoint+1):length(dat.shared)) {
  if (!is.null(dat.shared[[i]]) & nrow(dat.shared[[i]]) > 0) {
  	most.likely.node.assignments = cbind(most.likely.node.assignments, apply(dat.shared[[i]][,3:(ncol(dat.shared[[i]])-1)], 1, which.max))
  } else {
	  most.likely.node.assignments = cbind(most.likely.node.assignments, rep(NA, nrow(most.likely.node.assignments)))
  }
}
colnames(most.likely.node.assignments) = c("chr", "pos", vector_of_names)
write.table(most.likely.node.assignments, file=paste("2_mutation_assignments/tables/", samplename, "_most_likely_node_assignments.tsv", sep=""), sep="\t", col.names=T, quote=F, row.names=F)

#######################################################################
# Obtain cluster locations by taking the mean Subclonal fraction of its assigned mutations
#######################################################################
#' Get mean cluster locations based on the Subclonal fraction estimates of its assigned mutations
get.cluster.locations.shared.muts = function(dat.shared, vector_of_names) {
  mean.cluster.locs = list()
  for (i in 1:length(vector_of_names)) {
    dd = dat.shared[[i]]
    if (is.null(dd) || nrow(dd) == 0) {
      # No clusters available when a method hasn't produced mutation assignments.
      mean.cluster.locs[[i]] = data.frame()
    } else {
      cols_select = 3:(ncol(dd)-1)
      if (length(cols_select) == 1) {
	assignment = rep(1, nrow(dd))
      } else {
      	assignment = apply(dd[,cols_select], 1, which.max)
      }
      mean.subcl.frac = data.frame()
      for (cluster in unique(assignment)) {
        cluster.subcl.frac = dd$Subclonal.fraction[assignment==cluster]
        mean.subcl.frac = rbind(mean.subcl.frac, c(cluster, mean(cluster.subcl.frac), mean(cluster.subcl.frac)+sd(cluster.subcl.frac), sum(assignment==cluster)))
      }
      colnames(mean.subcl.frac) = c("Cluster", "Subclonal.fraction.mean", "Subclonal.fraction.sd", "Total.muts")
      mean.cluster.locs[[i]] = mean.subcl.frac
    }
  }
  return(mean.cluster.locs)
}

cluster.locations = get.cluster.locations.shared.muts(dat.shared, vector_of_names)
for (i in 1:length(vector_of_names)) {
  if (!is.null(cluster.locations[[i]]) & nrow(cluster.locations[[i]]) > 0) {
    write.table(file=paste("2_mutation_assignments/tables/", samplename, "_cluster_locations_shared_", vector_of_names[i], ".txt", sep=""), cluster.locations[[i]], row.names=F, quote=F, sep="\t")
  }
}

#' Obtain average CCF of clusters using all assigned mutations
get.cluster.locations.all.muts = function(list_of_tables, raw.data, vector_of_names) {
  mean.cluster.locs = list()
  for (i in 1:length(list_of_tables)) {
    dd = dat.shared[[i]]
    if (is.null(dd) || nrow(dd) == 0) {
	mean.cluster.locs[[i]] = data.frame()
    } else {
	chrpos = paste(raw.data[,1], raw.data[,2])    
    	chrpos_method = paste(dd[,1], dd[,2])
	selection = chrpos %in% chrpos_method
    	#selection = d$Chromosome %in% dd[,1] & d$Position %in% dd[,2]
    	dd = cbind(dd, raw.data$Subclonal.fraction[selection])
	cols_select = 3:(ncol(dd)-1)
	if (length(cols_select) == 1) {
		assignment = rep(1, nrow(dd))
	} else {
    		assignment = apply(dd[,cols_select], 1, which.max)
	}
    	
    	mean.subcl.frac = data.frame()
    	for (cluster in unique(assignment)) {
    	  cluster.subcl.frac = dd[assignment==cluster, ncol(dd)]
    	  mean.subcl.frac = rbind(mean.subcl.frac, c(cluster, mean(cluster.subcl.frac), mean(cluster.subcl.frac)+sd(cluster.subcl.frac), sum(assignment==cluster)))
    	}
    	colnames(mean.subcl.frac) = c("Cluster", "Subclonal.fraction.mean", "Subclonal.fraction.sd", "Total.muts")
    	mean.cluster.locs[[i]] = mean.subcl.frac
    }
  }
  return(mean.cluster.locs)
}

cluster.locations = get.cluster.locations.all.muts(list_of_tables, d, vector_of_names)
for (i in 1:length(vector_of_names)) {
  if (!is.null(cluster.locations[[i]]) & nrow(cluster.locations[[i]]) > 0) {
    write.table(file=paste("2_mutation_assignments/tables/", samplename, "_cluster_locations_all_", vector_of_names[i], ".txt", sep=""), cluster.locations[[i]], row.names=F, quote=F, sep="\t")
  }
}
