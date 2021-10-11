#This script is helpful for rasqualTools GC correction
#Sometimes there are issues with the large counts matrix and rasqualTools
#Sometimes rasqualCalculateSampleOffsets gives 'Error in Gs[bin + 1, ] : subscript out of bounds'

#Define Quantile and solveDuplication functions
Quantile <-
  function(x,k=20){x=rank(x,ties="random"); z=rep(0,length(x));for(i in 1:k){z=z+as.numeric(x<=quantile(x,i/k,na.rm=T))};k-z}
  
solveDuplication <- 
  function(x){
    x2=x
    for(i in as.numeric(names(table(x)[table(x)>1]))){
      ids=seq(length(x))[x==i];
      if(ids[1]>1){
        a = (x[ids[1]-1]+x[ids[1]])/2
      }else{
        a = x[1]
      }
      if(ids[length(ids)]<length(x)){
        b = (x[ids[length(ids)]]+x[ids[length(ids)]+1])/2
      }else{
        b = x[length(x)]
      }
      w = (b-a)/(length(ids)+1)
      x2[ids] = a+w*(1:length(ids))
    }
    x2
  }

#Calculate size factors
counts <- bulk_atac_counts2
library_size = colSums(counts)
size_factors = library_size/mean(library_size)

#Generate size matrix
size_matrix = matrix(rep(size_factors, nrow(counts)), nrow = nrow(counts), byrow = TRUE)
rownames(size_matrix) = rownames(counts)
colnames(size_matrix) = colnames(counts)

#Create gc vector
gene_data_offset <- as.data.frame(gene_data_offset)
rownames(gene_data_offset) = gene_data_offset$gene_id
gene_data_offset = gene_data_offset[rownames(counts), ]
gcvec = gene_data_offset$percentage_gc_content
bin = Quantile(gcvec, 200)

#Multiply size matrix by result matrix
x = sort(unlist(lapply(split(gcvec, bin), mean)))
x = solveDuplication(x)
S = apply(counts, 2, function(y) {
  unlist(lapply(split(y, bin), sum))[as.character(0:199)]
})
Fs = log(t(t(S)/apply(S, 2, sum))/apply(S, 1, sum) * sum(as.numeric(S)))
Gs = apply(Fs, 2, function(y) {
  smooth.spline(x, y, spar = 1)$y
})
result_matrix = exp(Gs[bin + 1, ])
rownames(result_matrix) = rownames(counts)

size_matrix = size_matrix * result_matrix

#Save size factor matrix
#saveRasqualMatrices(list(bulk = size_matrix), "/scratch/amt2ug/rasqualTools_bulk", file_suffix = "size_factors")
