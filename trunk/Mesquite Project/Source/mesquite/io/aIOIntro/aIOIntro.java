/* Mesquite (package mesquite.io).  Copyright 2000-2006 D. Maddison and W. Maddison. Version 1.11, June 2006.Disclaimer:  The Mesquite source code is lengthy and we are few.  There are no doubt inefficiencies and goofs in this code. The commenting leaves much to be desired. Please approach this source code with the spirit of helping out.Perhaps with your help we can be more than a few, and make Mesquite better.Mesquite is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY.Mesquite's web site is http://mesquiteproject.orgThis source code and its compiled class files are free and modifiable under the terms of GNU Lesser General Public License.  (http://www.gnu.org/copyleft/lesser.html)*/package mesquite.io.aIOIntro;/*~~  */import java.util.*;import java.awt.*;import mesquite.lib.*;import mesquite.lib.characters.*;import mesquite.lib.duties.*;/* ======================================================================== */public class aIOIntro extends PackageIntro {	/*.................................................................................................................*/	public boolean startJob(String arguments, Object condition, boolean hiredByName) { 		return true;  	 }  	 public Class getDutyClass(){  	 	return aIOIntro.class;  	 }	/*.................................................................................................................*/    	 public String getExplanation() {		return "Serves as an introduction to the import/export package for Mesquite.";   	 }   	/*.................................................................................................................*/    	 public String getName() {		return "Import and Export Package Introduction";   	 }	/*.................................................................................................................*/	/** Returns the name of the package of modules (e.g., "Basic Mesquite Package", "Rhetenor")*/ 	public String getPackageName(){ 		return "Import and Export Package"; 	}	/*.................................................................................................................*/	/** Returns citation for a package of modules*/ 	public String getPackageCitation(){ 		return "Maddison, D.R. & W.P. Maddison. 2006.  Import-export package for Mesquite, version 1.1."; 	}	/*.................................................................................................................*/  	 public String getPackageVersion() {		return "1.11";   	 }	/*.................................................................................................................*/  	 public String getPackageAuthors() {		return "David R. Maddison and Wayne P. Maddison";   	 }	/*.................................................................................................................*/	/** Returns whether there is a splash banner*/	public boolean hasSplash(){ 		return false; 	}}