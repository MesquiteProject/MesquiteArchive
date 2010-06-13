/* Mesquite (package mesquite.io).  Copyright 2000-2009 D. Maddison and W. Maddison. Version 2.72, December 2009.Disclaimer:  The Mesquite source code is lengthy and we are few.  There are no doubt inefficiencies and goofs in this code. The commenting leaves much to be desired. Please approach this source code with the spirit of helping out.Perhaps with your help we can be more than a few, and make Mesquite better.Mesquite is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY.Mesquite's web site is http://mesquiteproject.orgThis source code and its compiled class files are free and modifiable under the terms of GNU Lesser General Public License.  (http://www.gnu.org/copyleft/lesser.html)*/package mesquite.io.InterpretFastaProtein;/*~~  */import java.util.*;import java.awt.*;import mesquite.lib.*;import mesquite.lib.characters.*;import mesquite.lib.duties.*;import mesquite.categ.lib.*;import mesquite.io.lib.*;/* ============  a file interpreter for Protein  Fasta files ============*/public class InterpretFastaProtein extends InterpretFasta {/*.................................................................................................................*/	public boolean canExportEver() {  		 return true;  //	}/*.................................................................................................................*/	public boolean canExportProject(MesquiteProject project) {  		 return project.getNumberCharMatrices(ProteinState.class) > 0;  //	}/*.................................................................................................................*/	public boolean canExportData(Class dataClass) {  		return (dataClass==ProteinState.class);	}	/*.................................................................................................................*/	public boolean canImport(Class dataClass){		return (dataClass==ProteinState.class);	}/*.................................................................................................................*/	public CharacterData createData(CharactersManager charTask, Taxa taxa) {  		 return charTask.newCharacterData(taxa, 0, ProteinData.DATATYPENAME);  //	}/*.................................................................................................................*/	public CharacterData findDataToExport(MesquiteFile file, String arguments) { 		return getProject().chooseData(containerOfModule(), file, null, ProteinState.class, "Select data to export");	}/*.................................................................................................................*/	public void setFastaState(CharacterData data, int ic, int it, char c) { 		if (!(data instanceof ProteinData))			return;		((ProteinData)data).setState(ic, it, c);    // setting state to that specified by character c	}	/*.................................................................................................................*/	public  String getUnassignedSymbol(){		return "X";	}/*.................................................................................................................*/    	 public String getName() {		return "FASTA (protein)";   	 }/*.................................................................................................................*/ 	/** returns an explanation of what the module does.*/ 	public String getExplanation() { 		return "Imports and exports FASTA files that consist of amino acid sequence data." ;   	 }   	 }	