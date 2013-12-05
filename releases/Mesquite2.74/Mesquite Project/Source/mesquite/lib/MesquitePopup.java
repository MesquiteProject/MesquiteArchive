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
package mesquite.lib;

import java.awt.*;
import java.awt.event.*;
import mesquite.lib.duties.*;

public class MesquitePopup extends PopupMenu {
	Container c;
	public MesquitePopup(Container c){
		super();
		this.c = c;
	}
	Container getComponent(){
		return c;
	}
	
	public static Polygon getDropDownTriangle(){
		Polygon dropDownTriangle=new Polygon();
		dropDownTriangle.xpoints = new int[4];
		dropDownTriangle.ypoints = new int[4];
		dropDownTriangle.npoints=0;
		dropDownTriangle.addPoint(0, 0);
		dropDownTriangle.addPoint(6, 0);
		dropDownTriangle.addPoint(3,3);
		dropDownTriangle.addPoint(0, 0);
		dropDownTriangle.npoints=4;
		return dropDownTriangle;
	}
	public void addItem(String label, MesquiteModule module, MesquiteCommand respondCommand, String argument){
		MesquiteMenuItem m = new MesquiteMenuItem(label, module, respondCommand, argument);
		add(m);
	}
	public void addItem(String label, MesquiteCommand respondCommand, String argument){
		addItem(label, null, respondCommand, argument);
	}
	public void addItem(String label, MesquiteModule module, MesquiteCommand respondCommand){
		addItem(label, module, respondCommand, null);
	}


	public void addItems(String[] labels, MesquiteModule[] modules, MesquiteCommand[] respondCommands){
		if (labels!=null){
			for (int i=0; i<labels.length; i++) 
				addItem(labels[i], modules[i], respondCommands[i], null);
		}
	}


	public void showPopup(int x, int y){
		
		c.add(this);
		if (c instanceof MousePanel)
			((MousePanel)c).suppressEvents = true;
		show(c, x,y);
		if (c instanceof MousePanel)
			((MousePanel)c).suppressEvents = false;
		//c.remove(this);
	}
	

}



