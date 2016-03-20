<!--
 Copyright 2008-2012 Janssen Research & Development, LLC.

 Licensed under the Apache License, Version 2.0 (the "License");
 you may not use this file except in compliance with the License.
 You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

 Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the License for the specific language governing permissions and
 limitations under the License.
-->

<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN">
<html>
<head>
<title>subsetPanel.html</title>

<meta http-equiv="keywords" content="keyword1,keyword2,keyword3">
<meta http-equiv="description" content="this is my page">
<meta http-equiv="content-type" content="text/html; charset=ISO-8859-1">
<link rel="stylesheet" type="text/css" href="${resource(dir:'css', file:'datasetExplorer.css')}">

</head>

<body>
	<form>	
		<span class='AnalysisHeader'>
			Analysis Comparison
			
			<a href='JavaScript:D2H_ShowHelp(1512,helpURL,"wndExternal",CTXT_DISPLAY_FULLHELP )'>
				<img src="${resource(dir:'images',file:'help/helpicon_white.jpg')}" alt="Help" border=0 width=18pt style="margin-top:1pt;margin-bottom:1pt;margin-right:18pt;"/>
			</a>				
		</span><br />
		
		<br />
		
			<table class="AnalysisResults" style="font: 12px verdana,arial,helvetica,sans-serif;">
				<tr>
					<th>
						Difference Filter
					</th>
				</tr>
				<tr>
					<td>
						Only display entries where the difference is greater than : <input id="txtDiffGreaterThanAnalysis" type="text" size="4"></input> <br />
					</td>
				</tr>
				<tr>
					<td>
						<input type="radio" name="grpDifferenceTypeAnalysis" id="grpDifferenceType1Analysis" value="DIFF" checked>Difference
						<input type="radio" name="grpDifferenceTypeAnalysis" id="grpDifferenceType2Analysis" value="DIFF_PERCENT">Difference Percent<br />
					</td>
				</tr>
				<tr>
					<td>
						Use absolute value of difference : <input id="chkAbsoluteValueAnalysis" type = "checkbox"></input><br />				
					</td>
				</tr>
			</table>

			<br />
			
			<table class="AnalysisResults" style="font: 12px verdana,arial,helvetica,sans-serif;">
				<tr>
					<th>
						Other Filters
					</th>
				</tr>
				<tr>
					<td>
						Filter out categorical variables that remain constant across analyses : <input id="chkRemoveConstantAnalysis" type = "checkbox"></input><br />			
					</td>
				</tr>
			</table>		
		
		<br />
		<br />
		
		
		${analysisComparisonGrid}
		
		<br />
		<br />
		<br />
		
		<a class='AnalysisLink' href="${zipLink}">Download raw R data</a>
		
		<br />
		${rVersionInfo}
		
	</form>
</body>

</html>