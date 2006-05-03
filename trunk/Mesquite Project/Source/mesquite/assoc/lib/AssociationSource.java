/* Mesquite source code.  Copyright 1997-2005 W. Maddison and D. Maddison. Version 1.06, August 2005.Disclaimer:  The Mesquite source code is lengthy and we are few.  There are no doubt inefficiencies and goofs in this code. The commenting leaves much to be desired. Please approach this source code with the spirit of helping out.Perhaps with your help we can be more than a few, and make Mesquite better.Mesquite is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY.Mesquite's web site is http://mesquiteproject.orgThis source code and its compiled class files are free and modifiable under the terms of GNU Lesser General Public License.  (http://www.gnu.org/copyleft/lesser.html)*/package mesquite.assoc.lib;import java.awt.*;import java.util.*;import mesquite.lib.*;import mesquite.lib.duties.*;/* ======================================================================== *//**Supplies TaxaAssociation's for instance from a file or simulated.*/public abstract class AssociationSource extends MesquiteModule  {   	 public Class getDutyClass() {   	 	return AssociationSource.class;   	 } 	public String getDutyName() { 		return "Taxa Association Source";   	 }   	    	 public String[] getDefaultModule() {   	 	return new String[] {"#StoredAssociations"};   	 }   	 /** THE following assume a single taxa block is used as anchor point, and looks for others.     	  * These routines are independent of the two-taxa block routines, and snapshotting will only be done on the single taxa block routines */   	 /**Returns number of TaxaAssociation available.  If TaxaAssociation can be supplied indefinitely, returns MesquiteInteger.infinite*/    	public abstract int getNumberOfAssociations(Taxa taxa, CommandRecord commandRec);  	 /**Returns indexth TaxaAssociation*/   	public abstract TaxaAssociation getAssociation(Taxa taxa, int index, CommandRecord commandRec);   	   	 /**Returns current TaxaAssociation.*/   	public abstract TaxaAssociation getCurrentAssociation(Taxa taxa, CommandRecord commandRec);   	/**Returns String naming list number index*/   	public String getAssociationNameString(Taxa taxa, int index, CommandRecord commandRec) {   		return "";   	} 	 /** THE following assume two taxa blocks are specified.    	  * These routines are independent of the one-taxa block routines, and snapshotting will only be done on the single taxa block routines */  	 	 /**Returns number of TaxaAssociation available.  If TaxaAssociation can be supplied indefinitely, returns MesquiteInteger.infinite*/   	public abstract int getNumberOfAssociations(Taxa taxa1, Taxa taxa2, CommandRecord commandRec);  	 /**Returns indexth TaxaAssociation*/   	public abstract TaxaAssociation getAssociation(Taxa taxa1, Taxa taxa2, int index, CommandRecord commandRec);    	}