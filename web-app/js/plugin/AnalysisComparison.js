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

function submitAnalysisComparison(form){
	
	var analysis1VariableEle = Ext.get("divAnalysis1Variable");
	var analysis2VariableEle = Ext.get("divAnalysis2Variable");
	
	var analysis1VariableCode = "";
	var analysis2VariableCode = "";	
	
	var analysis1VariableConceptCode = "";
	var analysis2VariableConceptCode = "";
	
	var tempReturnArray;
	
	tempReturnArray = readConceptVariables2("divAnalysis1Variable");
		
	analysis1VariableConceptCode = tempReturnArray[0];
	analysis1VariableCode = tempReturnArray[1];
	
	tempReturnArray = readConceptVariables2("divAnalysis2Variable");
	
	analysis2VariableConceptCode = tempReturnArray[0];
	analysis2VariableCode = tempReturnArray[1];	

	
	/////////////////////////////////////////
	//Combine the different arrays so we can make sure the type matches across all input boxes.
	var finalNodeType = []

	//This is what we use to determine if we are running a modifier_cd analysis or a concept_cd analysis.
	var codeType
	
	var analysis1NodeType = createNodeTypeArrayFromDiv(analysis1VariableEle,"concepttablename")
	var analysis2NodeType = createNodeTypeArrayFromDiv(analysis2VariableEle,"concepttablename")
	
	if(analysis1NodeType[0] && analysis1NodeType[0] != "null") finalNodeType.push(analysis1NodeType[0]) 
	if(analysis2NodeType[0] && analysis2NodeType[0] != "null") finalNodeType.push(analysis2NodeType[0])
	
	//Distinct this final list.
	finalNodeType = finalNodeType.unique()
	
	if(finalNodeType.length > 1)
	{
		Ext.Msg.alert('Wrong input', 'You have selected inputs from different ontology trees, please only select nodes from the \'Navigate By Study\' or \'Across Trial\' tree.');
		return;			
	}
	
	if(finalNodeType[0] == "CONCEPT_DIMENSION")
	{
		codeType = "Concept"
	}
	
	if(finalNodeType[0] == "MODIFIER_DIMENSION")
	{
		codeType = "Modifier"
	}
	
	var variablesConceptCode = ""
	
	if(codeType == "Concept")
	{
		//Create a string of all the concept paths that we need to convert to codes.
		variablesConceptCode = analysis1VariableConceptCode+"|"+analysis2VariableConceptCode;	
		
		//Sloppy, but for now reassign the codes if we want the correct concept path.
		analysis1VariableCode = analysis1VariableConceptCode
		analysis2VariableCode = analysis2VariableConceptCode
		
		GLOBAL.codeType = "Concept"
	}
	else if(codeType == "Modifier")
	{
		//Create a string of all the modifier codes so we can put them in the clinical data query.
		variablesConceptCode = analysis1VariableCode+"|"+analysis2VariableCode;	
	
		analysis1VariableCode = analysis1VariableConceptCode
		analysis2VariableCode = analysis2VariableConceptCode	
		
		GLOBAL.codeType = "Modifier"
	}	
	/////////////////////////////////////////
	
	var formParams = {
			analysis1Variable:			analysis1VariableConceptCode,
			analysis2Variable:			analysis2VariableConceptCode,
			jobType:					'AnalysisComparison',
			codeType : 					codeType,
			analysis1VariableCode:		analysis1VariableCode,
			analysis2VariableCode:		analysis2VariableCode					
		};
	
	submitJob(formParams);
}

/**
 * Register drag and drop.
 * Clear out all gobal variables and reset them to blank.
 */
function loadAnalysisComparisonView(){
	registerAnalysisComparisonDragAndDrop();
}

function loadAnalysisComparisonOutput()
{
	//Create a jQuery table from the html one.
	var analysisTable = $j("#analysisComparisonTable").dataTable({
		"aoColumnDefs": [
		                 { "sType": "numeric", "aTargets": [ -1 ] }
		               ]
	});
	
	//Add an event handler to the text box for the difference filter.
	$j('#txtDiffGreaterThanAnalysis').keyup( function() { analysisTable.fnDraw(); } );
	
	//Add an event handler to the column type radio buttons for the difference filter.
	$j("input[name='grpDifferenceTypeAnalysis']").change( function() { analysisTable.fnDraw();} );
	
	//Add an event handler to the absolute value checkbox for the difference filter.
	$j("#chkAbsoluteValueAnalysis").change( function() { analysisTable.fnDraw();} );
	
	//Add an event handler to the checkbox for filtering the constant analysis1/analysis2 value.
	$j("#chkRemoveConstantAnalysis").change( function() { analysisTable.fnDraw();} );
	
}

/**
 * Clear the variable selection box
 * Clear all selection stored in global variables
 * Clear the selection display
 * @param divName
 */
function clearGroupBaseline(divName)
{
	//Clear the drag and drop div.
	var qc = Ext.get(divName);
	
	for(var i=qc.dom.childNodes.length-1;i>=0;i--)
	{
		var child=qc.dom.childNodes[i];
		qc.dom.removeChild(child);
	}	
	clearHighDimDataSelections(divName);
	clearSummaryDisplay(divName);
}

function registerAnalysisComparisonDragAndDrop()
{
	
	var analysis1NodeDiv = Ext.get("divAnalysis1Variable");
	var analysis2NodeDiv = Ext.get("divAnalysis2Variable");
	
	//Add the drop targets and handler function.
	dtgD = new Ext.dd.DropTarget(analysis1NodeDiv,{ddGroup : 'makeQuery'});
	dtgD.notifyDrop =  dropOntoVariableSelection;
	
	dtgD = new Ext.dd.DropTarget(analysis2NodeDiv,{ddGroup : 'makeQuery'});
	dtgD.notifyDrop =  dropOntoVariableSelection;
	
}

/* Custom filtering function which will filter data in the difference or difference percent column based on input parameters */
$j.fn.dataTableExt.afnFiltering.push(
	function( oSettings, aData, iDataIndex ) {
		var keepRow = true;
		
		//Only apply this row filtering if this is the baseline comparison grid.
		if(oSettings.sTableId == "analysisComparisonTable"){
			keepRow = false;
			
			/////////////////////////////////////////////
			//Difference Filtering
			/////////////////////////////////////////////
			//This is the minimum value from the text box.
			var iMin = document.getElementById('txtDiffGreaterThanAnalysis').value * 1;
			
			//This radio button tells us which column to filter on.
			var columnFilter = $j('input[@name="grpDifferenceTypeAnalysis"]:checked').val();
			
			//Difference = aData.length - 2
			//Difference Percent = aData.length - 1
			var filterColumn = aData.length - 2
			
			if(columnFilter == "DIFF_PERCENT")
			{
				filterColumn = aData.length - 1
			}
		
			//This is the data from the current row.
			var iDataValue = aData[filterColumn] == "-" ? 0 : aData[filterColumn]*1;
			
			//If the absolute value checkbox is checked, use the absolute value for the comparison.
			if(document.getElementById('chkAbsoluteValueAnalysis').checked == true)
			{
				iDataValue = Math.abs(iDataValue)
			}
			
			//Return a boolean if the value in the row is greater than the text entered into the text input box.
			if ( iMin == "" )
			{
				keepRow = true;
			}
			else if ( iMin < iDataValue )
			{
				keepRow = true;
			}		
			/////////////////////////////////////////////
			
			/////////////////////////////////////////////
			//Other Filtering
			/////////////////////////////////////////////
			//Continue as long as we are already keeping this row and the remove constant checkbox is checked.
			if(keepRow && document.getElementById('chkRemoveConstantAnalysis').checked == true)
			{
				var baselineColumn = aData.length - 4
				var newValueColumn = aData.length - 3
				
				//If Baseline and New Value are not parsable as numbers and they are the same, remove this record.
				if(isNaN(aData[baselineColumn]) && isNaN(aData[newValueColumn]) && (aData[baselineColumn] == aData[newValueColumn]))
				{
					keepRow = false;
				}
			}
			/////////////////////////////////////////////
	
		}
		return keepRow;
	}
);
