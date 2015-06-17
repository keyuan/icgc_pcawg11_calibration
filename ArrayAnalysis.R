#problems with array size if no.muts>30000, so we have to produce multiple arrays
GetIdentityArrayFromAssignments<-function(node.assignments,max.array.size = 30000){
	no.muts = length(node.assignments)
	if(no.muts>max.array.size){
		identity.strengths = list()
		no.arrays = ceiling(no.muts/max.array.size)
		print(paste("no.arrays=",no.arrays^2,sep=""))
		for(i in 1:(no.arrays^2)){
			print(paste("creating array ",i,sep=""))
			identity.strengths[[i]] = array(0,c(max.array.size,max.array.size))
		}
		#faster for large numbers of mutations (scales almost linearly rather than as no.muts^2)
		unique.clusters = unique(node.assignments)
		for(j in unique.clusters){
			indices = which(node.assignments==j)
			for(k in 1:no.arrays){
				for(l in 1:no.arrays){
					indices1 = indices[indices>(k-1)*max.array.size & indices<=k*max.array.size] - (k-1) * max.array.size
					indices2 = indices[indices>(l-1)*max.array.size & indices<=l*max.array.size] - (l-1) * max.array.size					
					identity.strengths[[(k-1)*no.arrays +l]][indices1,indices2] = 1
				}
			}			
		}
	}else{
		identity.strengths = array(0,c(no.muts,no.muts))
		if(no.muts>1000){
			#faster for large numbers of mutations (scales almost linearly rather than as no.muts^2)
			unique.clusters = unique(node.assignments)
			for(j in unique.clusters){
				indices = which(node.assignments==j)
				identity.strengths[indices,indices] = 1
			}
	
		}else{	
			for(m in 1:(no.muts-1)){
				identity.strengths[m,m] = 1
				for(n in (m+1):no.muts){
					identity.strengths[m,n] = identity.strengths[n,m] = sum(node.assignments[m] == node.assignments[n])
				}
			}
			identity.strengths[no.muts,no.muts] = 1
		}
	}
	return(identity.strengths)
}

#problems with array size if no.muts>30000, so we have to produce multiple arrays
GetIdentityArrayFromProbabilities<-function(node.probabilities,max.array.size = 30000){
	no.muts = nrow(node.probabilities)
	no.clusters = ncol(node.probabilities)	
	if(no.muts>max.array.size){
		identity.strengths = list()
		no.arrays = ceiling(no.muts/max.array.size)
		print(paste("no.arrays=",no.arrays^2,sep=""))
		for(i in 1:no.arrays){
			for(j in 1:no.arrays){
				#matrix cross product
				#identity.strengths[[(i-1)*no.arrays +j]] = node.probabilities[((i-1)*max.array.size+1):min(no.muts,i*max.array.size),] %*% t(node.probabilities[((j-1)*max.array.size+1):min(no.muts,j*max.array.size),])
				identity.strengths[[(i-1)*no.arrays +j]] = tcrossprod(node.probabilities[((i-1)*max.array.size+1):min(no.muts,i*max.array.size),], node.probabilities[((j-1)*max.array.size+1):min(no.muts,j*max.array.size),])
			}
		}			

	}else{
		
	  identity.strengths = tcrossprod(node.probabilities) # ke 
	  #matrix cross produc
		#identity.strengths = node.probabilities %*% t(node.probabilities)
	}
	return(identity.strengths)
}

