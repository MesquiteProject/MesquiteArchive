/* Mesquite source code.  Copyright 1997-2010 W. Maddison and D. Maddison. 
Version 2.74, October 2010.
Disclaimer:  The Mesquite source code is lengthy and we are few.  There are no doubt inefficiencies and goofs in this code. 
The commenting leaves much to be desired. Please approach this source code with the spirit of helping out.
Perhaps with your help we can be more than a few, and make Mesquite better.

Mesquite is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY.
Mesquite's web site is http://mesquiteproject.org

This source code and its compiled class files are free and modifiable under the terms of 
GNU Lesser General Public License.  (http://www.gnu.org/copyleft/lesser.html)
 */
package mesquite.molec.lib;

import mesquite.categ.lib.*;
import mesquite.molec.lib.*;

public class GenCodeEuplotid extends GeneticCode {

	public void setCode() {
		setStandardCode();
		setCode(U,G,A, ProteinData.CYS);
	}

	public String getName (){
		return "The Euplotid Nuclear Code";
	}

	public static String getShortName (){
		return "Euplotid Nuclear";
	}

	public String getNEXUSName(){
		return "nuc.euplotid";
	}

	public int getNCBITranslationTableNumber(){
		return 10;
	}

}
