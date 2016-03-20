/*************************************************************************   
* Copyright 2008-2012 Janssen Research & Development, LLC.
*
* Licensed under the Apache License, Version 2.0 (the "License");
* you may not use this file except in compliance with the License.
* You may obtain a copy of the License at
*
*     http://www.apache.org/licenses/LICENSE-2.0
*
* Unless required by applicable law or agreed to in writing, software
* distributed under the License is distributed on an "AS IS" BASIS,
* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
* See the License for the specific language governing permissions and
* limitations under the License.
******************************************************************/

package com.recomdata.transmart.data.association

import org.apache.commons.io.FileUtils;
import org.codehaus.groovy.grails.commons.ConfigurationHolder;

import groovy.util.ConfigObject;

class RModulesOutputRenderService {

    static transactional = true
	static scope= "request"

	def grailsApplication
	def zipService
	
	def config = ConfigurationHolder.config
	def tempFolderDirectory = config.RModules.tempFolderDirectory
	
	//This is a boolean indicating if we need to move the file before serving it to the user.
	boolean transferImageFile = config.RModules.transferImageFile
	
	//This is the URL we use to serve the user the image.
	def imageURL = config.RModules.imageURL
	
	def tempDirectory = ""
	def jobName = ""
	def jobTypeName = ""
	
	def zipLink = ""
	
    def initializeAttributes(jobName,jobTypeName,linksArray) {
		
		this.jobName = jobName
		this.jobTypeName = jobTypeName
		
		//Create the string that represents the directory to the temporary files.
		this.tempDirectory = "${tempFolderDirectory}${jobName}" + File.separator + "workingDirectory" + File.separator
	
		def zipLocation = ""
	
		//If we need to use a different location so that the image is under a web path, use the config here.
		if(transferImageFile)
		{
			String tempImageFolder = config.RModules.temporaryImageFolder
			String tempImageJobFolder = "${tempImageFolder}" + File.separator + "${this.jobName}" + File.separator
			
			File createDirectory = new File(tempImageJobFolder)
			
			//Determine if the folder for this job exists in the temp image directory.
			if(!createDirectory.exists())
			{
				createDirectory.mkdir()
			}
			
			def tempDirectoryFile = new File(this.tempDirectory)
	
			tempDirectoryFile.traverse(nameFilter:~/.*${jobTypeName}.*\.png/) 
			{
				currentImageFile -> 
			
				//For each of the image files we find, move them to the new directory.
				String tempImageLocation = "${tempImageJobFolder}" + File.separator + currentImageFile.name
	
				//Move the image to a location where we can actually render it.
				File oldImage = new File(currentImageFile.path);
				File newImage = new File(tempImageLocation);
				//TODO move FileUtils to Core
				FileUtils.copyFile(oldImage,newImage)
				
				String currentLink = "${imageURL}${jobName}/${currentImageFile.name}"
				linksArray.add(currentLink)
			};
			
			zipLocation = "${tempImageJobFolder}" + File.separator + "zippedData.zip"
			this.zipLink = "${imageURL}${jobName}/zippedData.zip"
			
			zipService.zipFolder(tempDirectory,zipLocation)
			
			zipLink = "${imageURL}${jobName}/zippedData.zip"
		}
    }
	
	def String fileParseLoop(tempDirectoryFile,fileNamePattern,fileNameExtractionPattern,fileParseFunction)
	{
		//This is the string we return.
		String parseValueString = ""
		
		//Reinitialize the text files array list.
		def ArrayList<String> txtFiles = new ArrayList<String>()
		
		//Loop through the directory create an array of txt files to be parsed.
		tempDirectoryFile.traverse(nameFilter:~fileNamePattern)
		{
			currentTextFile ->
			
			txtFiles.add(currentTextFile.path)
		}

		//Loop through the file path array and parse each of the files. We do this to make different tables if there are multiple files.
		txtFiles.each
		{
			//Parse out the name of the group from the name of the text file.
			def matcher = (it =~ fileNameExtractionPattern)
			
			if (matcher.matches() && txtFiles.size > 1)
			{
				//Add the HTML that will separate the different files.
				parseValueString += "<br /><br /><span class='AnalysisHeader'>${matcher[0][1]}</span><hr />"
			}
			
			//Create objects for the output file.
			File parsableFile = new File(it);
			
			//Parse the output files.
			parseValueString += fileParseFunction.call(parsableFile.getText())
		}
		
		return parseValueString
	}
	
	
	def parseVersionFile()
	{
		def tempDirectoryFile = new File(tempDirectory)
		String versionData = fileParseLoop(tempDirectoryFile,/.*sessionInfo.*\.txt/,/.*sessionInfo(.*)\.txt/,parseVersionFileClosure)
		
		return versionData
	}
	
	def parseLegendTable()
	{
		def tempDirectoryFile = new File(tempDirectory)
		String legendData = fileParseLoop(tempDirectoryFile,/.*legend.*\.txt/,/.*legend(.*)\.txt/,parseLegendTableClosure)
		
		return legendData
	}
	
	def parseVersionFileClosure = {
		statsInStr ->
		
		//Buffer that will hold the HTML we output.
		StringBuffer buf = new StringBuffer();
		
		buf.append("<br /><a href='#' onclick='\$(\"versionInfoDiv\").toggle()'><span class='AnalysisHeader'>R Version Information</span></a><br /><br />")
		
		buf.append("<div id='versionInfoDiv' style='display: none;'>")
		
		//This will tell us if we are printing the contents of the package or the session info. We will print the package contents in a table.
		Boolean packageCommand = false
		Boolean firstPackageLine = true
		
		statsInStr.eachLine
		{
			
			if(it.contains("||PACKAGEINFO||"))
			{
				packageCommand = true
				return;
			}
			
			if(!packageCommand)
			{
				buf.append(it)
				buf.append("<br />")
			}
			else
			{
				def currentLine = it.split("\t")
				
				if(firstPackageLine)
				{
					buf.append("<br /><br /><table class='AnalysisResults'>")
					buf.append("<tr>")
					currentLine.each()
					{
						currentSegment ->
						
						buf.append("<th>${currentSegment}</th>")
						
					}
					buf.append("</tr>")
					
					firstPackageLine = false
				}
				else
				{
					buf.append("<tr>")
					currentLine.each()
					{
						currentSegment ->
						
						buf.append("<td>${currentSegment}</td>")
						
					}
					buf.append("</tr>")
				}
			}
		}
		
		buf.append("</table>")
		buf.append("</div>")
		
		buf.toString();
	}
	

	
	def parseLegendTableClosure =
	{
		legendInStr ->
		
		//Buffer that will hold the HTML we output.
		StringBuffer buf = new StringBuffer();
		
		buf.append("<span class='AnalysisHeader'>Legend</span><br /><br />")
		buf.append("<table class='AnalysisResults'>")
		
		legendInStr.eachLine
		{

			//Start a new row.
			buf.append("<tr>")
			
			//Split each line.
			String[] strArray = it.split("\t");
			
			Integer cellCounter = 0;
			
			strArray.each
			{
				tableValue ->
				
				tableValue = tableValue.replace("|","<br />")
				
				buf.append("<th>${tableValue}</th>")
			}
			
			//End this row.
			buf.append("</tr>")
		}
		
		buf.append("</table><br />")
		//################################
		
		buf.toString();
	}

    /*
    * This moves an image file to the temporary directory so it can be rendered to the user.
    */
    def moveImageFile(currentFileLocation, newImageFileName, moveToDirectoryName)
    {
        String tempImageFolder = config.RModules.temporaryImageFolder
        String tempImageJobFolder = "${tempImageFolder}" + File.separator + "${moveToDirectoryName}" + File.separator

        //For each of the image files we find, move them to the new directory.
        String tempImageLocation = "${tempImageJobFolder}" + File.separator + newImageFileName

        //Move the image to a location where we can actually render it.
        File oldImage = new File(currentFileLocation);
        File newImage = new File(tempImageLocation);
        //TODO move FileUtils to Core
        FileUtils.copyFile(oldImage,newImage)

        String currentLink = "${imageURL}${moveToDirectoryName}/${newImageFileName}"

        //Delete the old directory.

        return currentLink
    }
	
}
