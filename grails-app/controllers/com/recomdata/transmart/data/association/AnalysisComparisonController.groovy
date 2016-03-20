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
import java.util.ArrayList;

class AnalysisComparisonController {

	def RModulesOutputRenderService
	
	def analysisComparisonOut = 
	{
		//This will be the array of image links.
		def ArrayList<String> imageLinks = new ArrayList<String>()
		
		//This will be the array of text file locations.
		def ArrayList<String> txtFiles = new ArrayList<String>()
		
		//Grab the job ID from the query string.
		String jobName = params.jobName
		
		//Gather the image links.
		RModulesOutputRenderService.initializeAttributes(jobName,"AnalysisComparison",imageLinks)
		
		String tempDirectory = RModulesOutputRenderService.tempDirectory
		
		//Traverse the temporary directory for the LinearRegression files.
		def tempDirectoryFile = new File(tempDirectory)
		
		//This string will be the HTML that represents our Linear Regression data.
		String analysisComparisonGrid = ""
		String rVersionInfo = ""
		
		analysisComparisonGrid = RModulesOutputRenderService.fileParseLoop(tempDirectoryFile,/.*analysisComparison.*\.txt/,/.*analysisComparison(.*)\.txt/,parseAnalysisComparisonStr)
		rVersionInfo = RModulesOutputRenderService.parseVersionFile()
		
		render(template: "/plugin/analysisComparison_out", model:[analysisComparisonGrid:analysisComparisonGrid,zipLink:RModulesOutputRenderService.zipLink,rVersionInfo:rVersionInfo], , contextPath:pluginContextPath)

	}
	
	def parseAnalysisComparisonStr =
	{
		legendInStr ->
		
		//Buffer that will hold the HTML we output.
		StringBuffer buf = new StringBuffer();
		
		buf.append("<table class='AnalysisResults' id='analysisComparisonTable'>")
		
		boolean firstLine = true;
		
		legendInStr.eachLine
		{
			if(firstLine)
			{
				//Start a new row.
				buf.append("<thead>")
				buf.append("<tr>")
				
				//Split each line.
				String[] strArray = it.split("\t");
				
				strArray.each
				{
					tableValue ->
					
					buf.append("<th>${tableValue}</th>")
				}
				
				//End this row.
				buf.append("</tr>")
				buf.append("</thead>")
				

				
				buf.append("<tbody>")
			}
			else
			{
				//Start a new row.
				buf.append("<tr>")
				
				//Split each line.
				String[] strArray = it.split("\t");
				
				strArray.each
				{
					tableValue ->
					if(tableValue == "") tableValue = "&nbsp;"
					buf.append("<td>${tableValue}</td>")
				}
				
				//End this row.
				buf.append("</tr>")
			}
			
			firstLine = false
		}
		
		buf.append("</tbody></table><br />")
		//################################
		
		buf.toString();
	}
	
}
