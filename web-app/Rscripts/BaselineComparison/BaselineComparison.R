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

BaselineComparison.Build <- 
function
(
input.dataFile,
concept.baseline,
concept.baseline.type = 'CLINICAL',
concept.others,
concept.others.type = 'CLINICAL'
)
{
	print("-------------------")
	print("BaselineComparison.R")
	print("BaselineComparison.Build")
	
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
	#This is the concept with slashes doubled up so the string can be used in REGEXP.
	concept.baseline.withescapes <- gsub("\\\\","\\\\\\\\",concept.baseline)
	
	#This is the baseline concept with the first two items trimmed from the path.
	concept.baseline.trimmed <- sub(pattern="^\\\\(.*?\\\\){2}",replacement="\\\\",x=concept.baseline,perl=TRUE)
	
	#Remove the first attribute from the end of the string.
	concept.baseline.trimmed.noattr <- sub(pattern="[^\\\\]*\\\\$",replacement="",x=concept.baseline.trimmed,perl=TRUE)	
	
	#Split the others concept into a list on the "|" character.
	otherConceptList <- strsplit(x=concept.others,split="\\|")
	
	#Escape the special characters in otherConceptList to use it in regexs.
	concept.others.withescapes <- gsub("\\\\","\\\\\\\\",unlist(otherConceptList))
	concept.others.withescapes <- gsub("?","\\",unlist(concept.others.withescapes), fixed=TRUE)
	concept.others.withescapes <- gsub(".","\\.",unlist(concept.others.withescapes), fixed=TRUE)
	concept.others.withescapes <- gsub("(","\\(",unlist(concept.others.withescapes), fixed=TRUE)
	concept.others.withescapes <- gsub(")","\\)",unlist(concept.others.withescapes), fixed=TRUE)
	concept.others.withescapes <- gsub("-","\\-",unlist(concept.others.withescapes), fixed=TRUE)
	concept.others.withescapes <- gsub("+","\\+",unlist(concept.others.withescapes), fixed=TRUE)
	concept.others.withescapes <- gsub("*","\\*",unlist(concept.others.withescapes), fixed=TRUE)
	concept.others.withescapes <- gsub("/","\\/",unlist(concept.others.withescapes), fixed=TRUE)
	##############################################

	##############################################
	#CONSOLIDATING DATA
	#For the time being the export includes all the leaf values, we need to filter out only the values under our input concepts.
	baselineData <- dataFile[grep(concept.baseline,dataFile$CONCEPT_PATH, fixed = TRUE),]
	otherData <- dataFile[grep(paste(unlist(concept.others.withescapes),collapse="|"),dataFile$CONCEPT_PATH),]
	
	if(NROW(baselineData)==0) stop("||FRIENDLY||No baseline data found.")
	if(NROW(otherData)==0) stop("||FRIENDLY||No data found for the 'Other' Variables.")	
	
	#Reset our master data frame.
	dataFile <- rbind(baselineData,otherData)
	##############################################
	
	##############################################
	#DETERMINE MISSING BASELINE
	#Create a list of patients who are missing baseline measurements. We handle them differently later.
	baselinePatients <- baselineData$PATIENT_NUM
	otherDataPatients <- otherData$PATIENT_NUM
	
	baselinePatients <- unique(baselinePatients)
	otherDataPatients <- unique(otherDataPatients)
	
	patientsWithoutBaseline <- setdiff(otherDataPatients,baselinePatients)
	##############################################		
	
	##############################################
	#ATTRIBUTE FRAME
	#Create a frame of attributes (e.g. Visit, Sponsor..). 

	#Trim to the columns we are intested in.
	attributeFrame <- otherData[,c('PATIENT_NUM','CONCEPT_PATH','VALUE')]
	
	#This is both the baseline and other concepts. We use this in a string replace statement.
	allConcepts <- paste(concept.baseline.withescapes,"|",paste(concept.others.withescapes,collapse="|"),sep="")
	
	#Two kinds of steps (Delimited by "\") exist under the passed in concept. Attributes and the actual value we are interested in.
	#In order to build the attribute frame we need to extract each step that is an attribute. 
	#This starts at the last step in the passed in concept and ends at the second to last concept in the CONCEPT_PATH_SHORT column.	
	numberOfStepsToTrim <- str_count(concept.baseline,'\\\\') - 2
	
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

	#For each of these we need to add 4 coumns to finalData. The name, The baseline value, the new value, and a difference (If the type is numeric)
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
		tempFrameRow$BASELINEVALUE=''
		tempFrameRow$NEWVALUE=''
		tempFrameRow$DIFFERENCE=0
		tempFrameRow$DIFFERENCEPERCENT=0
	
		#Reconstruct the name of the concept. The baseline lacks the first attribute.
		if(length(attributeValueList) > 2)
		{
			currentFullBaselineAttribute <- paste(concept.baseline.trimmed,paste(attributeValueList[3:length(attributeValueList)],collapse="\\"),"\\",valueConcept,sep="")
		}
		else
		{
			currentFullBaselineAttribute <- paste(concept.baseline.trimmed,valueConcept,sep="")
		}
		
		currentFullOtherAttribute <- paste(concept.baseline.trimmed.noattr,paste(attributeValueList[2:length(attributeValueList)],collapse="\\"),"\\",valueConcept,sep="")

		#Extract baseline and current value as numerics.
		#If the patient was one of the ones missing baseline measurements, indicate that here.
		if(tempFrameRow$PATIENT_NUM %in% patientsWithoutBaseline)
		{
			currentBaseline <- 'Baseline not measured'
		}
		else
		{
			currentBaseline <- as.character(dataFile[dataFile$PATIENT_NUM == tempFrameRow$PATIENT_NUM & dataFile$CONCEPT_PATH_SHORT == currentFullBaselineAttribute,c('VALUE')])
		}
		currentNewvalue <- as.character(dataFile[dataFile$PATIENT_NUM == tempFrameRow$PATIENT_NUM & dataFile$CONCEPT_PATH_SHORT == currentFullOtherAttribute,c('VALUE')])	

		#As long as currentBaseline is numeric, trim it to 2 decimal places.
		if(!identical(currentBaseline, character(0)))
		{
			if(suppressWarnings(!is.na(as.numeric(currentBaseline))) & length(currentBaseline) != 0)
			{
				currentBaseline <- format(round(as.numeric(currentBaseline), 3), nsmall = 3)
			}
		}
		
		#As long as currentNewvalue is numeric, trim it to 2 decimal places.
		if(!identical(currentNewvalue, character(0)))
		{
			if(suppressWarnings(!is.na(as.numeric(currentNewvalue))) & length(currentNewvalue) != 0)
			{
				currentNewvalue <- format(round(as.numeric(currentNewvalue), 3), nsmall = 3)
			}
		}		

		#Calculate the percent difference field.
		#If the currentBaseline is blank, we calculate the difference using a zero.
		if(identical(currentBaseline, character(0)))
		{
			if(length(currentNewvalue) != 0)
			{
				if(suppressWarnings(!is.na(as.numeric(currentNewvalue))))
				{
					currentDifference <- as.numeric(currentNewvalue)
					currentDifferencePercent <- "Infinity"
				}
			}
		}
		else if(!identical(currentBaseline, character(0)) & !identical(currentNewvalue, character(0)))
		{
			if(length(currentBaseline) != 0 & length(currentNewvalue) != 0)
			{
				if(suppressWarnings(!is.na(as.numeric(currentBaseline))) & suppressWarnings(!is.na(as.numeric(currentNewvalue))))
				{
					currentDifference <- as.numeric(currentNewvalue) - as.numeric(currentBaseline)
					currentDifferencePercent <- ((as.numeric(currentNewvalue) - as.numeric(currentBaseline)) /  as.numeric(currentBaseline)) * 100

					#As long as currentNewvalue and currentBaseline are numeric, trim the difference to 2 decimal places.
					currentDifference <- format(round(currentDifference, 3), nsmall = 3)
					currentDifferencePercent <- format(round(currentDifferencePercent, 3), nsmall = 3)
				}
				else
				{
					if(paste(as.character(currentBaseline),collapse=",") != paste(as.character(currentNewvalue),collapse=","))
					{
						currentDifference <- 1
						currentDifferencePercent <- 100
					}
				}
			}
		}
		
		if(length(currentBaseline) != 0 | length(currentNewvalue) != 0)
		{
			tempFrameRow$BASELINEVALUE 		<- paste(currentBaseline,collapse=",")
			tempFrameRow$NEWVALUE 			<- paste(currentNewvalue,collapse=",")
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
	
	finalColumnNames <- c(finalColumnNames,'Data Label','Baseline Value','New Value','Difference (New Value - Baseline)','Difference Percent')
	
	colnames(finalData) <- finalColumnNames
	
	#If after all that we don't have any values in the "New value" column, throw an error.
	if(length(unique(finalData[['New Value']])) == 1 && unique(finalData[['New Value']]) == "") stop("||FRIENDLY||There were no matching new values found for the specified baseline.")
	
	#We need MASS to dump the matrix to a file.
	require(MASS)

	#Write the final data file.
	write.table(finalData,'baselineComparison.txt',sep = "\t", row.names = FALSE, quote=FALSE)
	##########################################	
	print("-------------------")
	##########################################
}