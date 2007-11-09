/* Mesquite source code.  Copyright 1997-2007 W. Maddison and D. Maddison.
Version 2.0, September 2007.
Disclaimer:  The Mesquite source code is lengthy and we are few.  There are no doubt inefficiencies and goofs in this code. 
The commenting leaves much to be desired. Please approach this source code with the spirit of helping out.
Perhaps with your help we can be more than a few, and make Mesquite better.

Mesquite is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY.
Mesquite's web site is http://mesquiteproject.org

This source code and its compiled class files are free and modifiable under the terms of 
GNU Lesser General Public License.  (http://www.gnu.org/copyleft/lesser.html)
*/
package mesquite.trees.RandomlyRotate;

import mesquite.lib.*;
import mesquite.treefarm.lib.RndTreeModifier;

public class RandomlyRotate extends RndTreeModifier {

	/*.................................................................................................................*/
	public boolean startJob(String arguments, Object condition, boolean hiredByName) {
		return true;
	}
	
   	public boolean isPrerelease(){
   		return true;
   	}
	/*.................................................................................................................*/
   	 public void modifyTree(Tree tree, MesquiteTree modified, RandomBetween rng){
   		if (modified == null)
   			return;
   		modified.randomlyRotateDescendants(modified.getRoot(), rng, true);
   	}
	/*.................................................................................................................*/
    	 public String getName() {
		return "Randomly Rotate Nodes";
   	 }
	/*.................................................................................................................*/
  	 public String getExplanation() {
		return "Randomly rearranges the immediate descendants of each node; thus, it does not change the topology, just its appearance.";
   	 }
   	 
}
