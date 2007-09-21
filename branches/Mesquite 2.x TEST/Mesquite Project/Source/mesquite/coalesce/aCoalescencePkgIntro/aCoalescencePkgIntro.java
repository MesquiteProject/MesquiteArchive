/* Mesquite source code.  Copyright 1997-2007 W. Maddison. Version 2.0, September 2007.Disclaimer:  The Mesquite source code is lengthy and we are few.  There are no doubt inefficiencies and goofs in this code. The commenting leaves much to be desired. Please approach this source code with the spirit of helping out.Perhaps with your help we can be more than a few, and make Mesquite better.Mesquite is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY.Mesquite's web site is http://mesquiteproject.orgThis source code and its compiled class files are free and modifiable under the terms of GNU Lesser General Public License.  (http://www.gnu.org/copyleft/lesser.html)*/package mesquite.coalesce.aCoalescencePkgIntro;/*~~  */import java.util.*;import java.awt.*;import mesquite.lib.*;import mesquite.lib.duties.*;import mesquite.coalesce.lib.*;/* ======================================================================== */public class aCoalescencePkgIntro extends PackageIntro {	/*.................................................................................................................*/	public boolean startJob(String arguments, Object condition, boolean hiredByName) { 		return true;  	 }  	 public Class getDutyClass(){  	 	return aCoalescencePkgIntro.class;  	 }	/*.................................................................................................................*/	public boolean isSubstantive(){		return false;	}	/*.................................................................................................................*/    	 public String getExplanation() {		return "Serves as an introduction to the coalescence package.";   	 }   	/*.................................................................................................................*/    	 public String getName() {		return "Coalescence Package Introduction";   	 }	/*.................................................................................................................*/	/** Returns citation for a package of modules*/ 	public String getPackageCitation(){ 		return "Maddison, W.P. 2006.  Coalescence Package for Mesquite.  Version 1.1. http://mesquiteproject.org"; 	}	/*.................................................................................................................*/  	 public String getPackageVersion() {		return "1.11";   	 }	/*.................................................................................................................*/  	 public String getPackageAuthors() {		return "W. Maddison";   	 }	/*.................................................................................................................*/  	 public boolean hasSplash() {		return true;   	 }	/*.................................................................................................................*/	/** Returns the name of the package of modules (e.g., "Basic Mesquite Package", "Rhetenor")*/ 	public String getPackageName(){ 		return "Coalescence Package"; 	}}