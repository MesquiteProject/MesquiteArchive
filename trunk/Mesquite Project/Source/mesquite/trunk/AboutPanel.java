/* Mesquite source code.  Copyright 1997-2006 W. Maddison and D. Maddison. Version 1.11, June 2006.Disclaimer:  The Mesquite source code is lengthy and we are few.  There are no doubt inefficiencies and goofs in this code. The commenting leaves much to be desired. Please approach this source code with the spirit of helping out.Perhaps with your help we can be more than a few, and make Mesquite better.Mesquite is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY.Mesquite's web site is http://mesquiteproject.orgThis source code and its compiled class files are free and modifiable under the terms of GNU Lesser General Public License.  (http://www.gnu.org/copyleft/lesser.html)*/package mesquite.trunk;import java.awt.*;import java.awt.event.*;import java.awt.image.*;import java.util.*;import mesquite.lib.*;/* ======================================================================== *//** The Panel containing the Mesquite logo on the startup window */public class AboutPanel extends Panel {	Image logo;	HPanel superimposed = null;	public AboutPanel (Image logo) {		this.logo = logo;		setBackground(Color.white);	}	void superimposePanel(HPanel p){		superimposed = p;		add(p);		p.setSize(getBounds().width, AboutWindow.aboutHeight - MesquiteModule.textEdgeCompensationHeight - 16);	}	/*.................................................................................................................*/	public void paint(Graphics g) {	   	if (superimposed!=null || MesquiteWindow.checkDoomed(this))	   		return;		g.drawImage(logo,0,0,(ImageObserver)this);				g.setColor(ColorDistribution.lightYellow);		g.drawString("Version " + MesquiteModule.getMesquiteVersion() + MesquiteModule.getBuildVersion() , 8,95); //was 15		if (MesquiteTrunk.mesquiteTrunk.isPrerelease())			g.drawString(MesquiteModule.getBuildDate() , 8,110); //was 15		if (MesquiteTrunk.substantivePrereleasesExist) {			g.drawString( "Some installed modules", 21,12);			g.drawString("are pre-release versions.", 21,27);			//StringUtil.highlightString(g, "Touch on red alert symbol in windows for information.", 25,45);			g.drawImage(InfoBar.prereleaseImage,3,3,(ImageObserver)this);		}		g.drawString("http://mesquiteproject.org", 5,210);		g.drawString("Copyright (c) 1997-2006 W. & D. Maddison.", 5,225);				/*StringUtil.highlightString(g, "Version " + MesquiteModule.getMesquiteVersion() + MesquiteModule.getBuildVersion() , 8,95, Color.blue, Color.yellow); //was 15		if (MesquiteTrunk.mesquiteTrunk.isPrerelease())			StringUtil.highlightString(g, MesquiteModule.getBuildDate() , 8,110, Color.blue, Color.yellow); //was 15		if (MesquiteTrunk.substantivePrereleasesExist) {			StringUtil.highlightString(g, "Some installed modules", 21,30, Color.blue, Color.yellow);			StringUtil.highlightString(g, "are pre-release versions.", 21,45, Color.blue, Color.yellow);			//StringUtil.highlightString(g, "Touch on red alert symbol in windows for information.", 25,45, Color.blue, Color.yellow);			g.drawImage(InfoBar.prereleaseImage,3,21,(ImageObserver)this);		}		StringUtil.highlightString(g, "http://mesquiteproject.org", 4,10, Color.blue, Color.yellow);		StringUtil.highlightString(g, "Copyright (c) 1997-2006 W. & D. Maddison.", 5,195, Color.blue, Color.yellow);		*/		MesquiteWindow.uncheckDoomed(this);	}}