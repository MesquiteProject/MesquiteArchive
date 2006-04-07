/* Mesquite source code.  Copyright 1997-2005 W. Maddison and D. Maddison. Version 1.06, September 2005.Disclaimer:  The Mesquite source code is lengthy and we are few.  There are no doubt inefficiencies and goofs in this code. The commenting leaves much to be desired. Please approach this source code with the spirit of helping out.Perhaps with your help we can be more than a few, and make Mesquite better.Mesquite is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY.Mesquite's web site is http://mesquiteproject.orgThis source code and its compiled class files are free and modifiable under the terms of GNU Lesser General Public License.  (http://www.gnu.org/copyleft/lesser.html)*/package mesquite.align2.lib; import java.util.*;import java.awt.*;import java.io.*;import java.awt.image.*;import mesquite.lib.*;import mesquite.lib.characters.*;import mesquite.lib.duties.*;import mesquite.lib.table.*;import mesquite.align2.lib.*;import mesquite.categ.lib.*;/* ======================================================================== */public  class RestEnzyme implements Listable {	boolean hasBases = false;	String name;	String bases;	String sources;	int length;	int cutPosition;	long[] stateArray =null;	long[] reverseComplementArray = null;				/*.................................................................................................................*/ 	public RestEnzyme(String name, String bases, String cut, String sources) {		super();		this.name = name;		this.bases = bases;		length = bases.length();		if (cut.equalsIgnoreCase("?"))			cutPosition=-1;		else			cutPosition = MesquiteInteger.fromString(cut);		this.sources = sources;		createArrays();	}   	/*.................................................................................................................*/ 	public void createArrays() {		hasBases = false;		if (!StringUtil.blank(bases)) {			stateArray = new long[length];			for (int i = 0; i<length; i++)				stateArray[i] = DNAState.fromCharStatic(bases.charAt(i));			reverseComplementArray = new long[length];			for (int i = 0; i<length; i++)				reverseComplementArray[i] = DNAData.complement(stateArray[length-i-1]);			hasBases = true;		}	}  	/*.................................................................................................................*/  	/** Returns true iff dataArray is cut by this enzyme.  If the position requiredPosition is non-negative,   	then the position requiredPosition must be included within recognition sequence */ 	public boolean cutsAtPosition(long[] dataArray, int requiredPosition, boolean allowReverseComplement) { 		if (dataArray==null || stateArray == null) { 			return false; 		}// Debugg.println(name + " " + bases); 		if (dataArray.length<length) {// Debugg.println(name + ", length of data shorter than enzyme recog pattern:  "+ bases); 			return false; 		} 		for (int start = 0; start<dataArray.length-length; start++) { 			boolean match = true; 			for (int i =0; i<length; i++) {// Debugg.println("     " + i); 				if (!DNAState.isSubset(dataArray[start+i],stateArray[i])) { 					match = false; 					break; 				} 			} 			if (match) {// Debugg.println("     MATCHED"); 				if (requiredPosition<0) 					return true; 				else if (requiredPosition>=start && requiredPosition<=start+length-1) 					return true; 			} 		}		return false;	}   	/*.................................................................................................................*/ 	public String getSources() {		return sources;	}   	/*.................................................................................................................*/ 	public String getName() {		return name;	}   	/*.................................................................................................................*/ 	public String getBasesAsString() {		return bases;	}	 }