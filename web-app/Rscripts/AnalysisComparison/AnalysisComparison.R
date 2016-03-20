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
#BuildScatterDataFile
#Parse the i2b2 output file and create input files for a scatter plot.
###########################################################################

AnalysisComparison.Build <- 
function
(
input.dataFile,
concept.analysis1,
concept.analysis1.type = 'CLINICAL',
concept.analysis2,
concept.analysis2.type = 'CLINICAL'
)
{
	print("-------------------")
	print("AnalysisComparison.R")
	print("AnalysisComparison.Build")
	
	library(stringr)
	library(reshape2)
	library(plyr)
	
	#Read in the file with all the clinical data.
	dataFile <- data.frame(read.delim(input.dataFile))
	
	#Apply known column names so that we can access the columns easier.
	colnames(dataFile) <- defaultColumnList() 
	
	#Remove any of the security records.
	dataFile <- dataFile[dataFile$VALUE != "EXP:PUBLIC",]

	##############################################
	#VARIABLE PREPARATION
	
	#Split the others concept into a list on the "|" character.
	analysis1ConceptList <- strsplit(x=concept.analysis1,split="\\|")
	
	#Escape the special characters in analysis1ConceptList to use it in regexs.
	concept.analysis1.withescapes <- gsub("\\\\","\\\\\\\\",unlist(analysis1ConceptList))
	concept.analysis1.withescapes <- gsub("?","\\",unlist(concept.analysis1.withescapes), fixed=TRUE)
	concept.analysis1.withescapes <- gsub(".","\\.",unlist(concept.analysis1.withescapes), fixed=TRUE)
	concept.analysis1.withescapes <- gsub("(","\\(",unlist(concept.analysis1.withescapes), fixed=TRUE)
	concept.analysis1.withescapes <- gsub(")","\\)",unlist(concept.analysis1.withescapes), fixed=TRUE)
	concept.analysis1.withescapes <- gsub("-","\\-",unlist(concept.analysis1.withescapes), fixed=TRUE)
	concept.analysis1.withescapes <- gsub("+","\\+",unlist(concept.analysis1.withescapes), fixed=TRUE)
	concept.analysis1.withescapes <- gsub("*","\\*",unlist(concept.analysis1.withescapes), fixed=TRUE)
	concept.analysis1.withescapes <- gsub("/","\\/",unlist(concept.analysis1.withescapes), fixed=TRUE)
	
	#Split the others concept into a list on the "|" character.
	analysis2ConceptList <- strsplit(x=concept.analysis2,split="\\|")
	
	#Escape the special characters in analysis1ConceptList to use it in regexs.
	concept.analysis2.withescapes <- gsub("\\\\","\\\\\\\\",unlist(analysis2ConceptList))
	concept.analysis2.withescapes <- gsub("?","\\",unlist(concept.analysis2.withescapes), fixed=TRUE)
	concept.analysis2.withescapes <- gsub(".","\\.",unlist(concept.analysis2.withescapes), fixed=TRUE)
	concept.analysis2.withescapes <- gsub("(","\\(",unlist(concept.analysis2.withescapes), fixed=TRUE)
	concept.analysis2.withescapes <- gsub(")","\\)",unlist(concept.analysis2.withescapes), fixed=TRUE)
	concept.analysis2.withescapes <- gsub("-","\\-",unlist(concept.analysis2.withescapes), fixed=TRUE)
	concept.analysis2.withescapes <- gsub("+","\\+",unlist(concept.analysis2.withescapes), fixed=TRUE)
	concept.analysis2.withescapes <- gsub("*","\\*",unlist(concept.analysis2.withescapes), fixed=TRUE)
	concept.analysis2.withescapes <- gsub("/","\\/",unlist(concept.analysis2.withescapes), fixed=TRUE)

	#This is the baseline concept with the first two items trimmed from the path.
	concept.analysis1.trimmed <- sub(pattern="^\\\\(.*?\\\\){2}",replacement="\\\\",x=concept.analysis1,perl=TRUE)
	concept.analysis2.trimmed <- sub(pattern="^\\\\(.*?\\\\){2}",replacement="\\\\",x=concept.analysis2,perl=TRUE)
	##############################################

	##############################################
	#CONSOLIDATING DATA
	#For the time being the export includes all the leaf values, we need to filter out only the values under our input concepts.
	analysis1Data <- dataFile[grep(paste(unlist(concept.analysis1.withescapes),collapse="|"),dataFile$CONCEPT_PATH),]
	analysis2Data <- dataFile[grep(paste(unlist(concept.analysis2.withescapes),collapse="|"),dataFile$CONCEPT_PATH),]
	
	if(NROW(analysis1Data)==0) stop("||FRIENDLY||No Analysis 1 data found.")
	if(NROW(analysis2Data)==0) stop("||FRIENDLY||No Analysis 2 data found.")	
	
	#Reset our master data frame.
	dataFile <- rbind(analysis1Data,analysis2Data)
	##############################################
	
	##############################################
	#ATTRIBUTE FRAME
	#Create a frame of attributes (e.g. Visit, Sponsor..). 

	#Trim to the columns we are intested in.
	attributeFrame <- analysis2Data[,c('PATIENT_NUM','CONCEPT_PATH','VALUE')]
	
	#This is both the analysis concepts. We use this in a string replace statement.
	allConcepts <- paste(paste(concept.analysis1.withescapes,collapse="|"),"|",paste(concept.analysis2.withescapes,collapse="|"),sep="")

	#Two kinds of steps (Delimited by "\") exist under the passed in concept. Attributes and the actual value we are interested in.
	#In order to build the attribute frame we need to extract each step that is an attribute. 
	#This starts at the last step in the passed in concept and ends at the second to last concept in the CONCEPT_PATH_SHORT column.	
	numberOfStepsToTrim <- str_count(concept.analysis1,'\\\\') - 1

	#Take off the begining numberOfStepsToTrim levels to make the correct short path.
	attributeFrame$CONCEPT_PATH <- sub(pattern=paste("^\\\\(.*?\\\\){",numberOfStepsToTrim,"}",sep=""),replacement="\\\\",x=attributeFrame$CONCEPT_PATH,perl=TRUE)

	#If the value text appears in the concept_path, replace it.
	attributeFrame <- adply(attributeFrame, 1, function (x){x$CONCEPT_PATH <- gsub(pattern=paste(x$VALUE,"\\\\",sep=""),replacement="",x=x$CONCEPT_PATH,perl=TRUE)})
	
	attributeFrame$CONCEPT_PATH <- attributeFrame$V1
	attributeFrame$V1 <- NULL
	
	#Removing the trailing slash.
	attributeFrame$CONCEPT_PATH <- gsub("\\\\$","",attributeFrame$CONCEPT_PATH)
	
	#Take off the last item, it's the name of the value node.
	attributeFrame$CONCEPT_PATH <- sub(pattern="[^\\\\]*$",replacement="",x=attributeFrame$CONCEPT_PATH,perl=TRUE)

	#Remove the first slash.
	attributeFrame$CONCEPT_PATH <- gsub("^\\\\","",attributeFrame$CONCEPT_PATH)

	#We want to remove records that don't have a concept path at this point.
	attributeFrame <- attributeFrame[attributeFrame$CONCEPT_PATH != "",]
	
	#Convert our concept path into seperate column names.
	attributeFrame <- cbind(attributeFrame, do.call("rbind", strsplit(attributeFrame[, 2], "\\\\")))

	#Remove the concept path column, it is no longer needed.
	attributeFrame$CONCEPT_PATH <- NULL
	attributeFrame$VALUE <- NULL

	attributeFrame <- unique(attributeFrame)
	##############################################
		
	##############################################
	#VALUE CONCEPTS
	#Get a distinct list of the actual data values (e.g. Frequency, mutation..) we are dealing with.
	#This is achieved by taking the last part of each of the short concept paths and unique'ing it.
	listOfValueConcepts <- as.character(dataFile$CONCEPT_PATH_SHORT)
	listOfValueConcepts <- unique(gsub(".*\\\\","",listOfValueConcepts))
	##############################################

	#We initialize this so we can add records to it later. This will be the final data frame we write to a file.
	finalData <- data.frame()

	#For each of these we need to add 4 coumns to finalData. The name, The analysis1 value, the analysis2 value, and a difference (If the type is numeric)
	for(valueConcept in listOfValueConcepts)
	{
	  for (i in 1:nrow(attributeFrame))
	  {
		currentDifference <- 0
		currentDifferencePercent <- 100
		
		#Initialize this row of the data frame.
		tempFrameRow <- data.frame( PATIENT_NUM = attributeFrame[i, "PATIENT_NUM"] )
		
		#This list will be a list of all the attribute values.
		attributeValueList <- character(0)
		
		#For each attribute (Column in the attribute frame not including the patient ID column).		
		for(columnName in colnames(attributeFrame))
		{
			#Make a list of all the attributes as we loop through.
			attributeValueList <- c(attributeValueList,as.character(attributeFrame[i, columnName]))
			
			if(columnName != "PATIENT_NUM")
			{
				tempFrameRow[[columnName]] <- attributeFrame[i, columnName]
			}
		}
		
		tempFrameRow$DATA_LABEL=valueConcept
		tempFrameRow$ANALYSIS1VALUE=''
		tempFrameRow$ANALYSIS2VALUE=''
		tempFrameRow$DIFFERENCE=0
		tempFrameRow$DIFFERENCEPERCENT=0

		#Reconstruct the name of the concept.
		if(length(attributeValueList) > 2)
		{
			currentFullAnalysis1Attribute <- paste(concept.analysis1.trimmed,paste(attributeValueList[2:length(attributeValueList)],collapse="\\"),"\\",valueConcept,sep="")
		}
		else
		{
			currentFullAnalysis1Attribute <- paste(concept.analysis1.trimmed,valueConcept,sep="")
		}
		
		currentFullAnalysis2Attribute <- paste(concept.analysis2.trimmed,paste(attributeValueList[2:length(attributeValueList)],collapse="\\"),"\\",valueConcept,sep="")

		currentAnalysis1 <- as.character(dataFile[dataFile$PATIENT_NUM == tempFrameRow$PATIENT_NUM & dataFile$CONCEPT_PATH_SHORT == currentFullAnalysis1Attribute,c('VALUE')])
		currentAnalysis2 <- as.character(dataFile[dataFile$PATIENT_NUM == tempFrameRow$PATIENT_NUM & dataFile$CONCEPT_PATH_SHORT == currentFullAnalysis2Attribute,c('VALUE')])	

		#As long as currentAnalysis1 is numeric, trim it to 2 decimal places.
		if(!identical(currentAnalysis1, character(0)))
		{
			if(suppressWarnings(!is.na(as.numeric(currentAnalysis1))) & length(currentAnalysis1) != 0)
			{
				currentAnalysis1 <- format(round(as.numeric(currentAnalysis1), 3), nsmall = 3)
			}
		}
		
		#As long as currentAnalysis2 is numeric, trim it to 2 decimal places.
		if(!identical(currentAnalysis2, character(0)))
		{
			if(suppressWarnings(!is.na(as.numeric(currentAnalysis2))) & length(currentAnalysis2) != 0)
			{
				currentAnalysis2 <- format(round(as.numeric(currentAnalysis2), 3), nsmall = 3)
			}
		}		

		#Calculate the percent difference field.
		#If the currentAnalysis1 is blank, we calculate the difference using a zero.
		if(identical(currentAnalysis1, character(0)))
		{
			if(length(currentAnalysis2) != 0)
			{
				if(suppressWarnings(!is.na(as.numeric(currentAnalysis2))))
				{
					currentDifference <- as.numeric(currentAnalysis2)
					currentDifferencePercent <- "Infinity"
				}
			}
		}
		else if(!identical(currentAnalysis1, character(0)) & !identical(currentAnalysis2, character(0)))
		{
			if(length(currentAnalysis1) != 0 & length(currentAnalysis2) != 0)
			{
				if(suppressWarnings(!is.na(as.numeric(currentAnalysis1))) & suppressWarnings(!is.na(as.numeric(currentAnalysis2))))
				{
					currentDifference <- as.numeric(currentAnalysis2) - as.numeric(currentAnalysis1)
					currentDifferencePercent <- ((as.numeric(currentAnalysis2) - as.numeric(currentAnalysis1)) /  as.numeric(currentAnalysis1)) * 100

					#As long as currentAnalysis2 and currentAnalysis1 are numeric, trim the difference to 2 decimal places.
					currentDifference <- format(round(currentDifference, 3), nsmall = 3)
					currentDifferencePercent <- format(round(currentDifferencePercent, 3), nsmall = 3)
				}
				else
				{
					if(paste(as.character(currentAnalysis1),collapse=",") != paste(as.character(currentAnalysis2),collapse=","))
					{
						currentDifference <- 1
						currentDifferencePercent <- 100
					}
				}
			}
		}
		
		if(length(currentAnalysis1) != 0 | length(currentAnalysis2) != 0)
		{
			tempFrameRow$ANALYSIS1VALUE 	<- paste(currentAnalysis1,collapse=",")
			tempFrameRow$ANALYSIS2VALUE 	<- paste(currentAnalysis2,collapse=",")
			tempFrameRow$DIFFERENCE 		<- currentDifference
			tempFrameRow$DIFFERENCEPERCENT 	<- currentDifferencePercent
			
			finalData <- rbind(finalData,tempFrameRow)
		}

	  }
	}
	
	finalColumnNames <- c('Patient Number')
	
	columnAttributeCounter <- 1

	#For each attribute (Column in the attribute frame not including the patient ID column).		
	for(columnName in colnames(attributeFrame))
	{	
		if(columnName != "PATIENT_NUM")
		{
			finalColumnNames <- c(finalColumnNames,paste('Attribute',columnAttributeCounter,sep=""))
			columnAttributeCounter <- columnAttributeCounter + 1
		}
	}
	
	finalColumnNames <- c(finalColumnNames,'Data Label','Analysis 1 Value','Analysis 2 Value','Difference (Analysis 2 Value - Analysis 1 Value)','Difference Percent')
	
	colnames(finalData) <- finalColumnNames
	
	#If after all that we don't have any values in the "New value" column, throw an error.
	if(length(unique(finalData[['Analysis 2 Value']])) == 1 && unique(finalData[['Analysis 2 Value']]) == "") stop("||FRIENDLY||There were no matching new values found for the specified comparison.")
	
	#We need MASS to dump the matrix to a file.
	require(MASS)

	#Write the final data file.
	write.table(finalData,'analysisComparison.txt',sep = "\t", row.names = FALSE, quote=FALSE)
	##########################################	
	print("-------------------")
	##########################################
}