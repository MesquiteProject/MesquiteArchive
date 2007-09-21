/* Mesquite source code.  Copyright 1997-2007 W. Maddison and D. Maddison.Version 2.0, September 2007.Disclaimer:  The Mesquite source code is lengthy and we are few.  There are no doubt inefficiencies and goofs in this code. The commenting leaves much to be desired. Please approach this source code with the spirit of helping out.Perhaps with your help we can be more than a few, and make Mesquite better.Mesquite is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY.Mesquite's web site is http://mesquiteproject.orgThis source code and its compiled class files are free and modifiable under the terms of GNU Lesser General Public License.  (http://www.gnu.org/copyleft/lesser.html)*/package mesquite.lib.characters; import java.awt.*;import mesquite.lib.duties.*;import mesquite.lib.*;/*======================================================================== *//** a designation of codon positions for characters. */public class CodonPositionsSet  extends CharNumSet {		public CodonPositionsSet (String name, int numChars, CharacterData data) {		super(name, numChars, new MesquiteNumber(0), data);	}	public SpecsSet cloneSpecsSet(){		CodonPositionsSet ms = new CodonPositionsSet(new String(name), getNumberOfParts(),  data);		MesquiteNumber position = new MesquiteNumber();		for (int i=0; i< getNumberOfParts(); i++) {			placeValue(i, position);			ms.setValue(i, position);		}		return ms;	}	public SpecsSet makeSpecsSet(AssociableWithSpecs parent, int numParts){		if (!(parent instanceof CharacterData))			return null;		return new CodonPositionsSet("Untitled Codon Positions", numParts, (CharacterData)parent);	} 	/*.................................................................................................................*/	/** Add num parts just after "starting" (filling with default values)  */  	public boolean addParts(int starting, int num){  		boolean success = super.addParts(starting, num);		MesquiteNumber position = new MesquiteNumber(0);		if (getInt(starting)!=MesquiteInteger.unassigned)					for (int i=0; i< num; i++) 			setValue(i+starting+1, position);				return success;	}}