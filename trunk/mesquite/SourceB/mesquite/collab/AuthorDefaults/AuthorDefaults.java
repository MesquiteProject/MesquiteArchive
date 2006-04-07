/* Mesquite source code.  Copyright 1997-2005 W. Maddison and D. Maddison. 
Version 1.06, September 2005.
Disclaimer:  The Mesquite source code is lengthy and we are few.  There are no doubt inefficiencies and goofs in this code. 
The commenting leaves much to be desired. Please approach this source code with the spirit of helping out.
Perhaps with your help we can be more than a few, and make Mesquite better.

Mesquite is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY.
Mesquite's web site is http://mesquiteproject.org

This source code and its compiled class files are free and modifiable under the terms of 
GNU Lesser General Public License.  (http://www.gnu.org/copyleft/lesser.html)
*/
package mesquite.collab.AuthorDefaults;

import java.util.*;
import java.awt.*;
import mesquite.lib.*;
import mesquite.lib.duties.*;

/* ======================================================================== */
public class AuthorDefaults extends DefaultsAssistant {
	boolean nowarn = false;
	/*.................................................................................................................*/
	public boolean startJob(String arguments, Object condition, CommandRecord commandRec, boolean hiredByName) {
		loadPreferences();
		if (!nowarn && MesquiteModule.author.hasDefaultSettings()){
			if (!AlertDialog.query(containerOfModule(), "Author", "The Author for this account and machine has not yet been set.  Please go to the Set Author... menu item in the Defaults menu of the Log window to set an author name.  For the code, please indicate a short code unique in your collaborative group.  If you wish, turn off this warning.", "OK", "Turn off Warning")){
				
				nowarn = true;
				storePreferences();
			}
		}
		//following for collab only
		addMenuItem( "Set Author...", makeCommand("setAuthor",  this));
		return true;
  	 }
	public void processPreferencesFromFile (String[] prefs) {
		if (prefs!=null && prefs.length>0) {
				if ("-*".equals(prefs[0]))
					nowarn = true;
				else {
					MesquiteModule.author.setName(prefs[0]);
					if (prefs.length>1)
						MesquiteModule.author.setCode(prefs[1]);
    	 			}
		}
	}
	/*.................................................................................................................*/
	public String[] preparePreferencesForFile () {
		if (MesquiteModule.author.hasDefaultSettings() && nowarn){
			return (new String[] {"-*"});
		}
		else if (!StringUtil.blank(MesquiteModule.author.getName())) {
			if (MesquiteModule.author.getCode() != null) 
				return (new String[] {MesquiteModule.author.getName(), MesquiteModule.author.getCode()});
			else
				return (new String[] {MesquiteModule.author.getName()});
		}
		return null;
	}
	MesquiteInteger pos = new MesquiteInteger();
	/*.................................................................................................................*/
    	 public Object doCommand(String commandName, String arguments, CommandRecord commandRec, CommandChecker checker) {
    	 	if (checker.compare(MesquiteWindow.class, "Sets the author for this account and machine", null, commandName, "setAuthor")) {
    	 		MesquiteBoolean answer = new MesquiteBoolean(false);
    	 		MesquiteString resp1 = new MesquiteString(MesquiteModule.author.getName());
    	 		MesquiteString resp2 = new MesquiteString(MesquiteModule.author.getCode());
			MesquiteString.queryTwoStrings(containerOfModule(), "Set Author", "Author", "Author code (do not use a number!)", answer, resp1, resp2, false);
			if (answer.getValue()){
				MesquiteModule.author.setName(resp1.getValue());
				MesquiteModule.author.setCode(resp2.getValue());
			}
			storePreferences();
			setCurrentAllProjects();
			return null;

    	 	}
    	 	else
    	 		return super.doCommand(commandName, arguments, commandRec, checker);
   	 }
  	 
  	 private void setCurrentAllProjects(){
  	 //go through all projects changing current author or adding current author
  	 	Projects p = MesquiteTrunk.mesquiteTrunk.getProjectList();
  	 	for (int i=0; i<p.getNumProjects(); i++){
  	 		MesquiteProject proj = p.getProject(i);
  	 		ListableVector v = proj.getAuthors();
  	 		boolean found = false;
  	 		for (int ia = 0; ia< v.size(); ia++){
  	 			Author au = (Author)v.elementAt(ia);
 	 			if (au.isCurrent()) {
  	 				au.setName(MesquiteModule.author.getName());
  	 				au.setCode(MesquiteModule.author.getCode());
  	 				found = true;
  	 			}
  	 		}
  	 		if (!found){
				Author a = new Author();
				a.setName(MesquiteModule.author.getName());
				a.setCode(MesquiteModule.author.getCode());
				a.setCurrent(true);
				v.addElement(a, true);
  	 		}
  	 	}
  	 }
	/*.................................................................................................................*/
    	 public String getName() {
		return "Set Author";
   	 }
	/*.................................................................................................................*/
  	 public String getExplanation() {
		return "Sets the author for this machine and account.";
   	 }
}