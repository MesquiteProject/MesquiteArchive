/* Mesquite (package mesquite.io).  Copyright 2000-2009 D. Maddison and W. Maddison. 
Version 2.6, January 2009.
Disclaimer:  The Mesquite source code is lengthy and we are few.  There are no doubt inefficiencies and goofs in this code. 
The commenting leaves much to be desired. Please approach this source code with the spirit of helping out.
Perhaps with your help we can be more than a few, and make Mesquite better.

Mesquite is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY.
Mesquite's web site is http://mesquiteproject.org

This source code and its compiled class files are free and modifiable under the terms of 
GNU Lesser General Public License.  (http://www.gnu.org/copyleft/lesser.html)
*/

package mesquite.io.lib;

import java.awt.*;

import mesquite.lib.*;
	
public class PhylipExporterDialog extends ExporterDialog {
	int tnl = 10;
	IntegerField f;
	public PhylipExporterDialog (InterpretPhylip module, MesquiteWindow parent, String title, MesquiteInteger buttonPressed) {
		super(module, parent, title, buttonPressed);
		this.tnl = module.taxonNameLength;
	}
	/*.................................................................................................................*/
	public void completeAndShowDialog (boolean dataSelected, boolean taxaSelected) {
		 f = addIntegerField ("Maximum length of taxon names", tnl, 4, 1, 40);
		super.completeAndShowDialog(dataSelected, taxaSelected);
	}
	public int getTaxonNamesLength(){
		return f.getValue();
	}
}