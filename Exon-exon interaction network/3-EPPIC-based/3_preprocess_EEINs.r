##############################################################################################
# Purpose:  ## multiple instances of EEIs are present due to one uniprot mapping to more tan one PDB id
## This has no conseqeunce of EEI finding, but it can give different values for physiochemical
## properties of an EEI because of the differences in the numbers/types of amino acids resolved
## in a given PDB id. I take the EEI with the maximum number of resolved amino acids.
######################################################################################################

rm(list=ls())
library(Rcpp)

cppFunction("List filtereei(CharacterVector ex1, CharacterVector ex2, CharacterVector exon1, CharacterVector exon2, NumericVector allaa){

    int loop1 = ex1.size();
    int loop2 = exon1.size();
    NumericVector allpositions;

    for(int k=0; k<loop1; k++){

		NumericVector loopx;
		NumericVector tallaa;

        for(int j=0; j<loop2; j++){

            if((ex1[k] == exon1[j]) & (ex2[k] == exon2[j])){
                loopx.push_back(j);
                tallaa.push_back(allaa[j]);
            }

            if((ex1[k] == exon2[j]) & (ex2[k] == exon1[j])){
                loopx.push_back(j);
                tallaa.push_back(allaa[j]);
            }
        }

        int loop3 = tallaa.size();
        int MAX = 0;
        int position = 0;
        for(int j=0; j<loop3; j++){
        	if(tallaa[j] > MAX){
        		MAX = tallaa[j];
        		position = loopx[j];
        	}
            //Rcout << MAX << std::endl;
        }
        //Rcout << loop3 << std::endl;

        allpositions.push_back(position);
    }

    List L = List::create(allpositions);
	return L;
  
}")


temp <- data.table::fread('../data/EPPIC_EEIN.txt')
temp$allAA <- temp$AA1+temp$AA2

gnet <- igraph::as_data_frame(igraph::simplify(igraph::graph_from_data_frame(temp[,c(1,2)], directed=FALSE)))

tempf <- filtereei(gnet[[1]], gnet[[2]], temp$exon1, temp$exon2, temp$allAA)
ids <- tempf[[1]]+1
tempff <- temp[ids, ]
data.table::fwrite(tempff, paste0('../data/EPPIC_EEIN_filtered.txt'), row.names=FALSE, quote=FALSE, sep='\t')



