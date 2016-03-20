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

function submitBaselineComparison(form){
	
	var baselineVariableEle = Ext.get("divBaselineVariable");
	var otherVariableEle = Ext.get("divOtherVariable");
	
	var baselineVariableCode = "";
	var otherVariableCode = "";	
	
	var baselineVariableConceptCode = "";
	var otherVariableConceptCode = "";
	
	var tempReturnArray;
	
	tempReturnArray = readConceptVariables2("divBaselineVariable");
		
	baselineVariableConceptCode = tempReturnArray[0];
	baselineVariableCode = tempReturnArray[1];
	
	tempReturnArray = readConceptVariables2("divOtherVariable");
	
	otherVariableConceptCode = tempReturnArray[0];
	otherVariableCode = tempReturnArray[1];	

	/////////////////////////////////////////
	//Validation
	//Make sure the user entered some items into the variable selection boxes.
	if(baselineVariableConceptCode == '' && otherVariableConceptCode == '')
	{
		Ext.Msg.alert('Missing input', 'Please drag at least one concept into the Baseline Variable and Other Variable boxes.');
		return;
	}
	if(baselineVariableConceptCode == '')
	{
		Ext.Msg.alert('Missing input', 'Please drag at least one concept into the Baseline Variable box.');
		return;
	}
	if(otherVariableConceptCode == '')
	{
		Ext.Msg.alert('Missing input', 'Please drag at least one concept into the Other Variable box.');
		return;
	}	
	
	//Baseline variable can only have one input.
	if(baselineVariableConceptCode.indexOf("|") != -1)
	{
		Ext.Msg.alert('Incorect input', 'Please drag only 1 concept into the Baseline Variable box.');
		return;		
	}
	
	//Split the other box.
	var otherValues = otherVariableConceptCode.split("\|");
	
	//Count the number of slashes in the first item.
	var firstSlashCount = otherValues[0].split("\\").length - 1;
	
	//If there are other items, loop through them.
	for (var currentOther = 0; currentOther < otherValues.length; currentOther++) 
	{
		//If the item length is different from the first item, alert the user.
		if( (otherValues[currentOther].split("\\").length - 1) != firstSlashCount)
		{
			Ext.Msg.alert('Incorect input', 'All concepts in the Other Variable box must come from the same level in the Dataset Explorer Tree.');
			return;					
		}  
	}	
	
	//Make sure each of the concepts that were dragged in has the same path length by counting the separators.
	if((baselineVariableConceptCode.split("\\").length - 1) != firstSlashCount)
	{
		Ext.Msg.alert('Incorect input', 'All concepts in the Baseline Variable and Other Variable box must come from the same level in the Dataset Explorer Tree.');
		return;				
	}
	
	/////////////////////////////////////////
	
	/////////////////////////////////////////
	//Combine the different arrays so we can make sure the type matches across all input boxes.
	var finalNodeType = []

	//This is what we use to determine if we are running a modifier_cd analysis or a concept_cd analysis.
	var codeType
	
	var baselineNodeType = createNodeTypeArrayFromDiv(baselineVariableEle,"concepttablename")
	var otherNodeType = createNodeTypeArrayFromDiv(otherVariableEle,"concepttablename")
	
	if(baselineNodeType[0] && baselineNodeType[0] != "null") finalNodeType.push(baselineNodeType[0]) 
	if(otherNodeType[0] && otherNodeType[0] != "null") finalNodeType.push(otherNodeType[0])
	
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
		variablesConceptCode = baselineVariableConceptCode+"|"+otherVariableConceptCode;	
		
		//Sloppy, but for now reassign the codes if we want the correct concept path.
		baselineVariableCode = baselineVariableConceptCode
		otherVariableCode = otherVariableConceptCode
		
		GLOBAL.codeType = "Concept"
	}
	else if(codeType == "Modifier")
	{
		//Create a string of all the modifier codes so we can put them in the clinical data query.
		variablesConceptCode = baselineVariableCode+"|"+otherVariableCode;	
	
		baselineVariableCode = baselineVariableConceptCode
		otherVariableCode = otherVariableConceptCode	
		
		GLOBAL.codeType = "Modifier"
	}	
	/////////////////////////////////////////
	
	var formParams = {
			baselineVariable:			baselineVariableConceptCode,
			otherVariable:				otherVariableConceptCode,
			jobType:					'BaselineComparison',
			codeType : 					codeType,
			baselineVariableCode:		baselineVariableConceptCode,
			otherVariableCode:			otherVariableConceptCode					
		};
	
	submitJob(formParams);
}

/**
 * Register drag and drop.
 * Clear out all gobal variables and reset them to blank.
 */
function loadBaselineComparisonView(){
	registerBaselineComparisonDragAndDrop();
}

function loadBaselineComparisonOutput()
{
	//Create a jQuery table from the html one.
	var baselineTable = $j("#timepointComparisonTable").dataTable({
		"aoColumnDefs": [
		                 { "sType": "numeric", "aTargets": [ -1 ] }
		               ]
	});
	
	//Add an event handler to the text box for the difference filter.
	$j('#txtDiffGreaterThan').keyup( function() { baselineTable.fnDraw(); } );
	
	//Add an event handler to the column type radio buttons for the difference filter.
	$j("input[name='grpDifferenceType']").change( function() { baselineTable.fnDraw();} );
	
	//Add an event handler to the absolute value checkbox for the difference filter.
	$j("#chkAbsoluteValue").change( function() { baselineTable.fnDraw();} );
	
	//Add an event handler to the checkbox for filtering the constant baseline/new value.
	$j("#chkRemoveConstant").change( function() { baselineTable.fnDraw();} );
	
	//Add an event handler to the checkbox for filtering the missing baseline value.
	$j("#chkRemoveBaseline").change( function() { baselineTable.fnDraw();} );
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

function registerBaselineComparisonDragAndDrop()
{
	
	var baselineNodeDiv = Ext.get("divBaselineVariable");
	var otherNodeDiv = Ext.get("divOtherVariable");
	
	//Add the drop targets and handler function.
	dtgD = new Ext.dd.DropTarget(baselineNodeDiv,{ddGroup : 'makeQuery'});
	dtgD.notifyDrop =  dropOntoVariableSelection;
	
	dtgD = new Ext.dd.DropTarget(otherNodeDiv,{ddGroup : 'makeQuery'});
	dtgD.notifyDrop =  dropOntoVariableSelection;
	
}

/* Custom filtering function which will filter data in the difference or difference percent column based on input parameters */
$j.fn.dataTableExt.afnFiltering.push(
	function( oSettings, aData, iDataIndex ) {
		var keepRow = true;
		
		//Only apply this row filtering if this is the baseline comparison grid.
		if(oSettings.sTableId == "timepointComparisonTable"){
			keepRow = false;
			
			/////////////////////////////////////////////
			//Difference Filtering
			/////////////////////////////////////////////
			//This is the minimum value from the text box.
			var iMin = document.getElementById('txtDiffGreaterThan').value * 1;
			
			//This radio button tells us which column to filter on.
			var columnFilter = $j('input[@name="grpDifferenceType"]:checked').val();
			
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
			if(document.getElementById('chkAbsoluteValue').checked == true)
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
			if(keepRow && document.getElementById('chkRemoveConstant').checked == true)
			{
				var baselineColumn = aData.length - 4
				var newValueColumn = aData.length - 3
				
				//If Baseline and New Value are not parsable as numbers and they are the same, remove this record.
				if(isNaN(aData[baselineColumn]) && isNaN(aData[newValueColumn]) && (aData[baselineColumn] == aData[newValueColumn]))
				{
					keepRow = false;
				}
			}
			
			if(keepRow && document.getElementById('chkRemoveBaseline').checked == true)
			{
				var baselineColumn = aData.length - 4
				
				if(aData[baselineColumn] == "Baseline not measured")
				{
					keepRow = false;
				}
				
			}
			/////////////////////////////////////////////
	
		}
		return keepRow;
	}
);



