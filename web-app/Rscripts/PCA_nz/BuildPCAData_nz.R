###########################################################################
# Copyright 2008-2012 Janssen Research & Development, LLC.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
###########################################################################

###########################################################################
#BuildPCADataFile
###########################################################################

PCAData.build <- 
function
(
input.gexFile,
sample.subset1,
time.subset1,
tissues.subset1,
platform.subset1,
sample.subset2,
time.subset2,
tissues.subset2,
platform.subset2,
genes,
genes.aggregate,
output.dataFile="outputfile"
)
{
	##########################################
	print("-------------------")
	print("BuildPCAData.R")
	print("BUILDING PCA DATA")

	#Call a function to gather the gene expression data and filter the two subsets.
	geneExpressionMatrix <- pcaSubsetBuilder(
								GEXFile = input.gexFile,
								sample.subset1 = sample.subset1,
								time.subset1 = time.subset1,
								tissues.subset1 = tissues.subset1,
								platform.subset1 = platform.subset1,
								sample.subset2 = sample.subset2,
								time.subset2 = time.subset2,
								tissues.subset2 = tissues.subset2,
								platform.subset2 = platform.subset2,
								gene.list = genes,
								gene.aggregate = genes.aggregate							
								)
	
	#Pull out a subset of columns we are interested in.
	geneExpressionMatrix <- geneExpressionMatrix[c('PATIENT.ID','ZSCORE','PROBE.ID','GENE_SYMBOL','SUBSET','ASSAY.ID')]
	
	#Trim the probe.id field.
	geneExpressionMatrix$GENE_SYMBOL <- gsub("^\\s+|\\s+$", "",geneExpressionMatrix$GENE_SYMBOL)	
	
	#Pull the columns we are interested in out of the data.
	finalFrame <- geneExpressionMatrix[c('PATIENT.ID','ZSCORE','PROBE.ID','GENE_SYMBOL','SUBSET')]	

	colnames(finalFrame) <- c('PATIENT.ID','VALUE','PROBE.ID','GENE_SYMBOL','SUBSET')

	#We need MASS to dump the matrix to a file.
	require(MASS)	
	
	#Write the final data file.
	write.matrix(finalFrame,output.dataFile,sep = "\t")
	print("-------------------")
	##########################################
}


#**************************************************
pcaSubsetBuilder <-
function
(
	GEXFile = '',
	sample.subset1,
	time.subset1,
	tissues.subset1,
	platform.subset1,
	sample.subset2,
	time.subset2,
	tissues.subset2,
	platform.subset2,
	gene.list = '',
	gene.aggregate
)
{
  #allGEXData <- data.frame(read.delim(GEXFile));
      
	#We have two subsets that we have to filter on using the passed in criteria.
	#subset1 <- allGEXData[grep("^subset1",allGEXData$SUBSET),]
	#subset2 <- allGEXData[grep("^subset2",allGEXData$SUBSET),]
  
  subset1 <- getSubset(GEXFile, "SUBSET1")
  subset2 <- getSubset(GEXFile, "SUBSET2")
  
  
	subset1 <- gexPCABuilder(
				GEXData = subset1,
				sampleType = sample.subset1,
				timepointType = time.subset1,
				tissueType = tissues.subset1,
				platform.type = platform.subset1,
				gene.list = gene.list,
				gene.aggregate = gene.aggregate,
				probe.average = TRUE,
				subsetname.replace = TRUE,
				data.reduce = FALSE,
				subset.name = 'SUBSET1');

	if(nrow(subset2) > 0)
	{
		subset2 <- gexPCABuilder(
					GEXData = subset2,
					sampleType = sample.subset2,
					timepointType = time.subset2,
					tissueType = tissues.subset2,
					platform.type = platform.subset2,
					gene.list = gene.list,
					gene.aggregate = gene.aggregate,
					probe.average = TRUE,
					subsetname.replace = TRUE,
					data.reduce = FALSE,
					subset.name -> 'SUBSET2');

		geneExpressionMatrix <- rbind(subset1,subset2)
	}
	else
	{
		geneExpressionMatrix <- subset1
	}

	return(geneExpressionMatrix)
}
#########################


#**************************************************
#  extract subset data from Netezza directoly
#   added on 01-16-2014 by HX
#**************************************************

getSubset <- function(GEXFile, subset)
{  
  library(nzr);
  nzDisconnect();
  nzConnect('biomart_user', 'biomart_user', '172.20.5.8', 'tsmrt');
  
  mRNATableName <- extractTableName(GEXFile);
  sampleTableName <- gsub(' ', '', paste(mRNATableName, "_SAMPLE"));
  
  # extract subset
  mRNATableName1 <- gsub(' ', '', paste(mRNATableName, "_", toupper(subset)))
  sampleTableName1 <- gsub(' ', '', paste(sampleTableName, "_", toupper(subset)))

  subset <- data.frame()
  
  if(nzExistTable(mRNATableName1) && nzExistTable(sampleTableName1)){
      query <- c("select t1.sourcesystem_cd as patient_num,t2.sample_type,t2.TIMEPOINT,t2.TISSUE_TYPE,t1.gpl_id,t1.assay_id,t1.RAW_INTENSITY as value,t1.zscore,t1.LOG_INTENSITY as log2ed,t1.probe_id,t1.PROBESET_ID,t1.gene_id,t1.gene_symbol,t1.SEARCH_KEYWORD_ID as search_id, 'S1' as subset from ",mRNATableName1," t1, ", sampleTableName1, " t2 where t1.assay_id=t2.assay_id");
      subset <- nzQuery(query);
      #nzQuery(c("drop table ", sampleTableName1));
      #nzQuery(c("drop table ", mRNATableName1));
      nzDisconnect();
      
      colnames(subset)[1] <- c("PATIENT.ID");
      colnames(subset)[2] <- c("SAMPLE.TYPE");
      colnames(subset)[4] <- c("TISSUE.TYPE");
      colnames(subset)[5] <- c("GPL.ID");
      colnames(subset)[6] <- c("ASSAY.ID");
      colnames(subset)[10] <- c("PROBE.ID");
      colnames(subset)[11] <- c("PROBESET.ID");
 }

  return (subset)                         
}
