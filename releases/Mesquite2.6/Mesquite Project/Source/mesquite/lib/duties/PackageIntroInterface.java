/* Mesquite source code.  Copyright 1997-2009 W. Maddison and D. Maddison.Version 2.6, January 2009.Disclaimer:  The Mesquite source code is lengthy and we are few.  There are no doubt inefficiencies and goofs in this code. The commenting leaves much to be desired. Please approach this source code with the spirit of helping out.Perhaps with your help we can be more than a few, and make Mesquite better.Mesquite is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY.Mesquite's web site is http://mesquiteproject.orgThis source code and its compiled class files are free and modifiable under the terms of GNU Lesser General Public License.  (http://www.gnu.org/copyleft/lesser.html) */package mesquite.lib.duties;import java.awt.*;import mesquite.lib.*;/* ======================================================================== *//**Serves to introduce a package*/public interface PackageIntroInterface  {	/*.................................................................................................................*/	/** Returns whether there is a splash screen banner and possibly also URL*/	public boolean hasSplash();	/*.................................................................................................................*/	/** Returns the URL of document shown when splash screen icon touched. By default, returns path to module's manual*/	public String getSplashURL();	/*.................................................................................................................*/	/** Returns citation for a package of modules*/	public String getPackageCitation();	/*.................................................................................................................*/	/** Returns version for a package of modules*/	public String getPackageVersion();	/*.................................................................................................................*/	/** Returns authors for a package of modules*/	public String getPackageAuthors();	/*.................................................................................................................*/	/** Returns date released for a package of modules*/	public String getPackageDateReleased();	/*.................................................................................................................*/	/** Returns the name of the package of modules*/	public String getPackageName();	/*.................................................................................................................*/	/** Returns version for a package of modules as an integer*/	public int getPackageVersionInt();	/*.................................................................................................................*/	/** Returns the URL for the web page about the package*/	public String getPackageURL(); 	/*.................................................................................................................*/	/** Returns the integer version of Mesquite for which this package was first released*/	public int getVersionOfFirstRelease();		public boolean getHideable();}