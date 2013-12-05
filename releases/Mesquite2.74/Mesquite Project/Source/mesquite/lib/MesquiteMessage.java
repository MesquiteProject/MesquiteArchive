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

 import java.util.*;
 
/*=======================*/
/** A utility class with method to duplicate System.out.println, or to generate user alerts.
These methods currently all call System.out.println, but will be differentiable in the future (e.g., some warnings will
go only to the console, while others will bring up dialog boxes).*/
public class MesquiteMessage extends Debugg {
	
	public static void warnProgrammer(String s) {
		MesquiteModule.showLogWindow(true);
		if (MesquiteTrunk.mesquiteTrunk !=null)
			MesquiteTrunk.mesquiteTrunk.logln(s);
		else
			System.out.println(s);
	}
	public static void notifyProgrammer(String s) {
		MesquiteModule.showLogWindow(true);
		if (MesquiteTrunk.mesquiteTrunk !=null)
			MesquiteTrunk.mesquiteTrunk.logln(s);
		else
			System.out.println(s);
	}
	public static void warnUser(String s) {
		if (MesquiteTrunk.mesquiteTrunk !=null)
			MesquiteTrunk.mesquiteTrunk.logln(s);
		else
			System.out.println(s);
	}
	public static void discreetNotifyUser( String s) {
		if (MesquiteTrunk.mesquiteTrunk !=null)
			MesquiteTrunk.mesquiteTrunk.discreetAlert( s);
	}
	public static void notifyUser(String s) {
		if (MesquiteTrunk.mesquiteTrunk !=null)
			AlertDialog.notice(MesquiteTrunk.mesquiteTrunk.containerOfModule(), "Notice", s);
	}
	public static void beep() {
           try {
        	   if (MesquiteWindow.GUIavailable)
        	   java.awt.Toolkit.getDefaultToolkit().beep();
           }
        	   catch(Throwable e){
        	   }
	}
	
	
	public static void logCurrentTime(String message) {
		String s = getCurrentTime();
		if (MesquiteTrunk.mesquiteTrunk !=null)
			MesquiteTrunk.mesquiteTrunk.logln(message + " " + s);
		else
			System.out.println(message + " " + s);
		
	}
	public static String getCurrentTime() {
		Date dnow = new Date(System.currentTimeMillis());
		return dnow.toString();
	}
	
	
	
}

