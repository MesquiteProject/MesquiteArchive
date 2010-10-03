/* Mesquite source code.  Copyright 1997-2010 W. Maddison and D. Maddison.Version 2.74, October 2010.Disclaimer:  The Mesquite source code is lengthy and we are few.  There are no doubt inefficiencies and goofs in this code. The commenting leaves much to be desired. Please approach this source code with the spirit of helping out.Perhaps with your help we can be more than a few, and make Mesquite better.Mesquite is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY.Mesquite's web site is http://mesquiteproject.orgThis source code and its compiled class files are free and modifiable under the terms of GNU Lesser General Public License.  (http://www.gnu.org/copyleft/lesser.html)*/package mesquite.categ.lib;import java.awt.*;import java.util.*;import mesquite.lib.*;import mesquite.lib.characters.*;import mesquite.lib.duties.*;/* ======================================================================== *//**A class for an array of  categorical character states for many characters, at each of the taxa  or nodes.*/public class MCategoricalEmbedded extends MCategoricalDistribution {	public MCategoricalEmbedded (CharacterData data) {		super(data.getTaxa());		this.data = data;	}	/*..........................................  MCategoricalEmbedded  ..................................................*/	/**returns state set of character ic in taxon */	public long getState (int ic, int it){		return CategoricalState.dataBitsMask &((CategoricalData)data).getState(ic, it);	}	/**returns raw state set of character ic in taxon */	public long getStateRaw (int ic, int it){		return ((CategoricalData)data).getStateRaw(ic, it);	}	/**return CharacterDistribution object for character ic */	public CharacterDistribution getCharacterDistribution (int ic){		return data.getCharacterDistribution(ic);	}	public int getNumTaxa(){		return data.getNumTaxa();	}	public int getNumChars(){		return data.getNumChars();	}	public String getName(){		return data.getName();	}}