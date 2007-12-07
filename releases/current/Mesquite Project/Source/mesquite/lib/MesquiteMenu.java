/* Mesquite source code.  Copyright 1997-2007 W. Maddison and D. Maddison. Version 2.01, December 2007. Disclaimer:  The Mesquite source code is lengthy and we are few.  There are no doubt inefficiencies and goofs in this code.  The commenting leaves much to be desired. Please approach this source code with the spirit of helping out. Perhaps with your help we can be more than a few, and make Mesquite better.  Mesquite is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY. Mesquite's web site is http://mesquiteproject.org  This source code and its compiled class files are free and modifiable under the terms of  GNU Lesser General Public License.  (http://www.gnu.org/copyleft/lesser.html) */package mesquite.lib;import java.awt.*;import java.awt.event.*;import mesquite.lib.duties.*;/* ======================================================================== *//** A menu.*/public class MesquiteMenu extends Menu implements Commandable, Listable{	public static int totalSubmenus = 0;	public boolean recycle = false;	MesquiteMenuSpec spec;	//static MesquiteSubmenu[] menus = new MesquiteSubmenu[255];	long id = 0;	static long numInstances = 0;		public MesquiteMenu(MesquiteMenuSpec spec) {		super(spec.getLabel(), true);  // true to designate as tearoff; doesn't seem to work on macos		id = numInstances++;		if (spec.getLabel() == null) {			MesquiteMessage.println("menu with no name: ");			setEnabled(false);		}		if (!spec.isEnabled())			setEnabled(false);		this.spec = spec;	}	public MesquiteMenu(String label) {		super(label, true);  // true to designate as tearoff; doesn't seem to work on macos		id = numInstances++;		if (label == null) {			MesquiteMessage.println("menu with no name: ");			setEnabled(false);		}		this.spec = null;	}	public MesquiteMenuSpec getSpecification(){		return spec;	}	public static MesquiteMenu getMenu(MesquiteMenuSpec spec) {				return new MesquiteMenu(spec);	}	public long getID(){		return id;	}	public String getName(){		if (this instanceof MesquiteSubmenu)			return "Submenu: " + getLabel();		else			return "Menu: " + getLabel();	}	public void listItems(){		String output = "";		if (this instanceof MesquiteSubmenu)			output += ("  Submenu: " + this.getLabel() +"\n");		else			output += ("  Menu: " + this.getLabel() +"\n");		for (int j = 0; j<this.getItemCount(); j++) {			if (!("-".equals(this.getItem(j).getLabel()))){				output += ("      " + (j+1) + "  --  " + this.getItem(j).getLabel());				if(this.getItem(j) instanceof Menu)					output += " >";				output += "\n";			}		}		System.out.println(output);		System.out.println("Enter number to select item");	}	/*.................................................................................................................*/	/** A request for the object to perform a command.  It is passed two strings, the name of the command and the arguments.*/	public Object doCommand(String commandName, String arguments, CommandChecker checker) { 		if (checker.compare(getClass(), null, null, commandName, "show")) {			listItems();		}		else {			MesquiteInteger pos = new MesquiteInteger();			int im = MesquiteInteger.fromFirstToken(commandName, pos);			if (MesquiteInteger.isCombinable(im)){				im--;				if (im >=0 && im<getItemCount()){					MenuItem item = getItem(im);					if (item instanceof MesquiteMenu){						System.out.println(((MesquiteMenu)item).getName() + " selected");						((MesquiteMenu)item).listItems();						ConsoleThread.setConsoleObjectCommanded(item, true, false);					}					else if (item instanceof MesquiteMenuItem){						if (this instanceof MesquiteSubmenu && ((MesquiteSubmenu)this).getCommand()!= null){							System.out.println("Item " + (im+1) + " selected");							((MesquiteSubmenu)this).chooseItem(im);						}						else  {							System.out.println(((MesquiteMenuItem)item).getLabel() + " selected");							((MesquiteMenuItem)item).chooseItem(arguments);  						}					}				}			}		}		return null;			}	boolean itemWithSameLabelExists(String label){		for (int i=0; i<getItemCount(); i++)			if (getItem(i).getLabel().equals(label))				return true;		return false;	}	static boolean itemWithSameLabelExists(Menu menu, String label){		if (menu == null)			return false;		for (int i=0; i<menu.getItemCount(); i++)			if (menu.getItem(i).getLabel().equals(label))				return true;		return false;	}		public MenuItem add(MenuItem mmi) {		if (mmi==null)			return null;		if (mmi instanceof Menu)			totalSubmenus++;		return super.add(mmi);	}	public static void add(Menu menu, MenuItem mmi) {		if (mmi==null  || menu == null)			return;				menu.add(mmi);	}}