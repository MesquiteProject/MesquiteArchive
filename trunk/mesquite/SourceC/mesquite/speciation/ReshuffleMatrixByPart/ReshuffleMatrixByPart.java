/* Mesquite source code.  Copyright 1997-2004 W. Maddison and D. Maddison. Version 1.01, January 2004.Disclaimer:  The Mesquite source code is lengthy and we are few.  There are no doubt inefficiencies and goofs in this code. The commenting leaves much to be desired. Please approach this source code with the spirit of helping out.Perhaps with your help we can be more than a few, and make Mesquite better.Mesquite is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY.Mesquite's web site is http://mesquiteproject.orgThis source code and its compiled class files are free and modifiable under the terms of GNU Lesser General Public License.  (http://www.gnu.org/copyleft/lesser.html)*/package mesquite.speciation.ReshuffleMatrixByPart;/*~~  */import java.util.*;import java.awt.*;import mesquite.lib.*;import mesquite.lib.characters.*;import mesquite.lib.duties.*;import mesquite.charMatrices.lib.RandomMatrixModifier;/* ======================================================================== */public class ReshuffleMatrixByPart extends RandomMatrixModifier {	/*.................................................................................................................*/	public boolean startJob(String arguments, Object condition, CommandRecord commandRec, boolean hiredByName) {  	 	return true;   	 }	/*.................................................................................................................*/  	public void modifyMatrix(MCharactersDistribution matrix, MAdjustableDistribution modified, RandomBetween rng, CommandRecord commandRec){		if (matrix==null || modified == null)			return;					Taxa taxa = matrix.getTaxa();		TaxaPartition partition = (TaxaPartition)taxa.getCurrentSpecsSet(TaxaPartition.class);		TaxaGroup[] groups = null;		if (partition != null)			groups = partition.getGroups();					int numTaxa = matrix.getNumTaxa();		int numChars = matrix.getNumChars();		if (modified.getNumTaxa()!=numTaxa || modified.getNumChars()!=numChars)			modified.setSize(numChars, numTaxa);		if (groups == null){			for (int ic = 0; ic<numChars; ic++) {		   		modified.transferFrom(ic, matrix.getCharacterDistribution(ic));		   		for (int i=0; i < (numTaxa-1); i++) {		   			int sh = rng.randomIntBetween(i, numTaxa-1);		   			if (sh != i)		   				modified.tradeStatesBetweenTaxa(ic, i, sh);		   		}	 	 	} 	 	} 	 	else {			for (int ic = 0; ic<numChars; ic++) {		   		modified.transferFrom(ic, matrix.getCharacterDistribution(ic));		   		for (int g = 0; g< groups.length; g++){ //examining group		   			int numTaxaInGroup = partition.getNumberInGroup(groups[g]);		   						   		for (int i=0; i < (numTaxaInGroup-1); i++) {			   			int sh = rng.randomIntBetween(i, numTaxaInGroup-1);			   			if (sh != i)			   				modified.tradeStatesBetweenTaxa(ic, getTaxonNumber(i, groups[g], taxa, partition), getTaxonNumber(sh, groups[g], taxa, partition));			   		}		   		}	 	 	} 	 	}   	}   	int getTaxonNumber(int numInGroup, TaxaGroup group, Taxa taxa, TaxaPartition partition){   		int count = 0;   		for (int it = 0; it<taxa.getNumTaxa(); it++){   			if (partition.getTaxaGroup(it) == group){   				if (count == numInGroup)   					return it;   				count++;   			}   		}   		return -1;   	}   	   		/*.................................................................................................................*/    	 public String getName() {   		return "Reshuffle Matrix Within Taxa Partitions";   	 }	/*.................................................................................................................*/  	 public boolean isPrerelease() {		return true;   	 }	/*.................................................................................................................*/  	 public boolean showCitation() {		return true;   	 }	/*.................................................................................................................*/ 	/** returns an explanation of what the module does.*/ 	public String getExplanation() { 		return "Shuffles (permutes) character states among taxa within each taxa partition, for each character." ;   	 }   	 }