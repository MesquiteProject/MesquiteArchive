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
package mesquite.categ.lib;

import java.awt.*;
import java.util.*;
import mesquite.lib.*;
import mesquite.lib.characters.*;
import mesquite.lib.characters.CharacterData;
import mesquite.lib.duties.*;

/* ======================================================================== */
/** An object a module can create and pass back to store in module info.  Tests whether module will be compatible with
passed object.  Classes of modules will have known ways of responding to particular classes of objects, e.g. character sources
should test whether they can handle given CharacterState types.*/
public class RequiresAnyMolecularData extends CompatibilityTest {
	public boolean isCompatible(Object obj, MesquiteProject project, EmployerEmployee prospectiveEmployer){
		if (obj == null)
			return true;
		if (obj instanceof RequiresAnyMolecularData)
			return true;
		if (obj instanceof CompatibilityTest)
			return false;
		if (!(obj instanceof Class))
			return true;
		if (!CharacterState.class.isAssignableFrom((Class)obj) &&  !CharacterData.class.isAssignableFrom((Class)obj))
			return true;
		//To pass the test, the datatype MUST be a molecular or a subclass of molecular
		return  (MolecularState.class.isAssignableFrom((Class)obj) || MolecularData.class.isAssignableFrom((Class)obj));
	}
	public Class getAcceptedClass(){
		return MolecularState.class;
	}
}


