/* Mesquite source code.  Copyright 1997-2006 W. Maddison and D. Maddison.Version 1.11, June 2006.Disclaimer:  The Mesquite source code is lengthy and we are few.  There are no doubt inefficiencies and goofs in this code. The commenting leaves much to be desired. Please approach this source code with the spirit of helping out.Perhaps with your help we can be more than a few, and make Mesquite better.Mesquite is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY.Mesquite's web site is http://mesquiteproject.orgThis source code and its compiled class files are free and modifiable under the terms of GNU Lesser General Public License.  (http://www.gnu.org/copyleft/lesser.html)Modified May 02 especially for annotations*/package mesquite.lib;import java.awt.*;import java.util.*;import java.awt.print.*;import java.awt.geom.*;import java.io.*;/* ======================================================================== */public class MesquitePrintJob {	public static final int AUTOFIT = 3;	public static final int FIT_LANDSCAPE = 1;	public static final int FIT_PORTRAIT = 0;	public static final int NATURAL = -1;		int fitToPage = AUTOFIT; ;	Component component;	Font font;	PrintJob job1;	PrinterJob job2;	String name;	MesquiteFrame frame;	Dimension dimension;	P2 p2;	PageFormat pf = null;			protected MesquitePrintJob(MesquiteFrame frame, String name){		this.frame = frame;		if (name == null)			name = "Print";		this.name = name;		if (MesquiteWindow.Java2Davailable){			try {				p2 = new P2();			}			catch (Throwable e){				MesquiteFile.throwableToLog(this, e);				MesquiteTrunk.mesquiteTrunk.alert("Exception or Error in making P2; details in Mesquite log file");			}		}	}	public static MesquitePrintJob getPrintJob(MesquiteFrame frame, String name, int fitToPage){ 		MesquiteTrunk.mesquiteTrunk.logln("Getting print job");		MesquitePrintJob job = new MesquitePrintJob(frame, name);		if (job == null)			return null;		if (job.preparePrint(fitToPage))			return job;		else {			return null;		}	}	public boolean preparePrint(int fitToPage) {		this.fitToPage = fitToPage;				if (!MesquiteWindow.Java2Davailable){			job1 = Toolkit.getDefaultToolkit().getPrintJob(frame, name, (Properties)null);			if (job1 == null)				MesquiteMessage.warnProgrammer("Error: Java 1.1 print job not obtained.");			return job1 !=null;		}	        job2 = PrinterJob.getPrinterJob();        	if (job2 ==null) {        		MesquiteTrunk.mesquiteTrunk.alert("Error: no printer job returned in preparePrint");        		return false;	        }	        job2.setPrintable(p2);	        		if (fitToPage == AUTOFIT){  //fits landscape or portrait depending on what requires least shrinkage	        	pf = job2.defaultPage();	        	if (pf ==null) {        			MesquiteTrunk.mesquiteTrunk.alert("Error: no defaultPage returned in preparePrint (autofit)");	        		return false;	        	}		}		else if (fitToPage == FIT_LANDSCAPE){ //fits into landscape	        	pf = job2.defaultPage();	        	if (pf ==null) {        			MesquiteTrunk.mesquiteTrunk.alert("Error: no defaultPage returned in preparePrint (landscape)");	        		return false;	        	}	        	pf.setOrientation(PageFormat.LANDSCAPE);	        }		else if (fitToPage == FIT_PORTRAIT){ //fits into portrait			pf = job2.defaultPage();	        	if (pf ==null) {        			MesquiteTrunk.mesquiteTrunk.alert("Error: no defaultPage returned in preparePrint (portrait)");	        		return false;	        	}	        	pf.setOrientation(PageFormat.PORTRAIT);	        }	        else { //prints at current size using page setup (thus, possibly over multiple pages)			pf = job2.pageDialog(job2.defaultPage());			if (pf == null) {        			MesquiteTrunk.mesquiteTrunk.alert("Error: no page dialog returned in preparePrint");				return false;			}		}		job2.validatePage(pf);		job2.setPrintable(p2, pf);  		if (MesquiteThread.isScripting())			return true;		else			return job2.printDialog();	}	public void printComponent(Component component, Dimension dim, Font font) {		this.font = font;		this.component = component;				if (dim == null)			dimension = component.getSize();		else			dimension = dim;		if (!MesquiteWindow.Java2Davailable){			if (job1 == null)				return;			//dimension = component.getSize();			Object fitReturned = null;			if (fitToPage >= 0 && frame instanceof Fittable){  //fits landscape or portrait depending on what requires least shrinkage				if (dimension.width<=0 || dimension.height <=0)					return;				Dimension page = job1.getPageDimension();				double shrinkWidth = page.width*1.0/dimension.width;				double shrinkHeight = page.height*1.0/dimension.height;				double shrink;				if (shrinkWidth> shrinkHeight)					shrink = shrinkHeight;				else					shrink = shrinkWidth;				fitReturned = ((Fittable)frame).fit(new Dimension((int)(dimension.width*shrink), (int)(dimension.height*shrink)));				Graphics pg = job1.getGraphics();	 			if (pg!=null) {	 				if (font !=null)	 					pg.setFont(font);	 				component.printAll(pg);	 				pg.dispose();	 			}				((Fittable)frame).unfit(fitReturned);			}			else {				Dimension page = job1.getPageDimension();				for (int vertical = 0; vertical<dimension.height; vertical+= page.height) {					for (int horizontal = 0; horizontal<dimension.width; horizontal+= page.width) {						Graphics pg = job1.getGraphics();			 			if (pg!=null) {  			 				if (font !=null)			 					pg.setFont(font);		    	 				pg.translate(-horizontal, -vertical);			    	 			component.printAll(pg);			 		    	 	pg.dispose();		   	 			}		   	 		}    	 			}			}			return;		}	        		if (job2==null)			return;		if (fitToPage == AUTOFIT){  //fits landscape or portrait depending on what requires least shrinkage			if (dimension.width<=0 || dimension.height <=0)				return;	        	PageFormat pf = job2.defaultPage();	        	pf.setOrientation(PageFormat.LANDSCAPE);			double shrinkWidth = pf.getImageableWidth()*1.0/dimension.width;			double shrinkHeight = pf.getImageableHeight()*1.0/dimension.height;			double shrinkRatioLANDSCAPE;			if (shrinkWidth< shrinkHeight)				shrinkRatioLANDSCAPE = shrinkHeight/shrinkWidth;			else				shrinkRatioLANDSCAPE = shrinkWidth/shrinkHeight;	        	pf.setOrientation(PageFormat.PORTRAIT);			shrinkWidth = pf.getImageableWidth()*1.0/dimension.width;			shrinkHeight = pf.getImageableHeight()*1.0/dimension.height;			double shrinkRatioPORTRAIT;			if (shrinkWidth< shrinkHeight)				shrinkRatioPORTRAIT = shrinkHeight/shrinkWidth;			else				shrinkRatioPORTRAIT = shrinkWidth/shrinkHeight;							if (shrinkRatioPORTRAIT>shrinkRatioLANDSCAPE)	        		pf.setOrientation(PageFormat.LANDSCAPE);			job2.validatePage(pf);			job2.setPrintable(p2, pf);  		}				try {			job2.print();  		} 		catch (Exception ex) {	              //  ex.printStackTrace();		}			}		public void printText(String s, Font font) {		if (s == null || font==null)			return;		if (!MesquiteWindow.Java2Davailable){			if (job1 == null)				return;			Dimension page = job1.getPageDimension();			StringBuffer sB = new StringBuffer(s);			StringInABox sBox = new StringInABox(sB, font, page.width);  			int tot =  sBox.getHeight();			int lastY = 0;			boolean done = false; 			while (lastY<tot && !done){				int r = sBox.getRemainingHeight(lastY);				if (r> page.height)					r=page.height;				Graphics pg = job1.getGraphics();	 			if (pg!=null) {	 				lastY = sBox.draw(pg, 0, 0, lastY, r);  	 				pg.dispose();				}				else done = true;			}			return;		}	        		if (job2==null)			return;        	if (pf == null){			pf = job2.defaultPage(); //currently disallow page setup because at least on OS X cropping prevents the 			//pf = job2.pageDialog(job2.defaultPage());        	}		job2.validatePage(pf);		job2.setPrintable(new P2Text(s, font), pf);  				try {			job2.print();  		} 		catch (Exception ex) {	              //  ex.printStackTrace();		}			}	public void end(){		if (!MesquiteWindow.Java2Davailable){			if (job1 == null)				return;			job1.end();		}	}	public class P2Text implements Printable {		StringBuffer sB;		Font font;		StringInABox sBox;		public P2Text(String s, Font font){			sB = new StringBuffer(s);			this.font = font;			sBox = new StringInABox(sB, font, 1000);  		}		public int print(Graphics g, PageFormat pf, int pi) throws PrinterException {			Graphics2D g2 = (Graphics2D)g;			if (font !=null)				g.setFont(font);			AffineTransform at = g2.getTransform();			double scale = at.getScaleX();			int effectivePageWidth = (int)(pf.getImageableWidth()/scale);			sBox.setWidth(effectivePageWidth);				double dPagesHigh =  (sBox.getHeight() * scale)/pf.getImageableHeight(); //need to add pixels for number of pages in case lines cut			int pagesHigh = (int)dPagesHigh;			if ((double)pagesHigh != dPagesHigh)				pagesHigh++;						if (pi >= pagesHigh) {			    return Printable.NO_SUCH_PAGE;			}				double dEffectivePageHeight = pf.getImageableHeight()/scale;						g2.scale(1/scale, 1/at.getScaleY());  //why this is & following line are needed isn't clear, but without it OS X 10.2 crops instead of scales			g2.scale(scale, at.getScaleY());//why this is & preceding line are needed isn't clear, but without it OS X 10.2 crops instead of scales			g2.translate((int)(pf.getImageableX() + 0.5),(int)(-dEffectivePageHeight*pi + pf.getImageableY() + 0.5));	 		sBox.draw(g2, 0, 0, 0, 99999999); 		        return Printable.PAGE_EXISTS;		}	}	public class P2 implements Printable {		public int print(Graphics g, PageFormat pf, int pi) throws PrinterException {			Graphics2D g2 = (Graphics2D)g;			if (font !=null)				g.setFont(font);			AffineTransform at = g2.getTransform();			double scale = at.getScaleX();			if (fitToPage >= 0){ // not using natural size; fit into page				if (pi >= 1 || dimension == null || dimension.width <= 0 || dimension.height <= 0) {				    return Printable.NO_SUCH_PAGE;				}				double shrinkWidth = pf.getImageableWidth()*1.0/dimension.width;				double shrinkHeight = pf.getImageableHeight()*1.0/dimension.height;				double shrink;				if (shrinkWidth< shrinkHeight)					shrink = shrinkWidth;				else					shrink = shrinkHeight;				g2.translate((int)(pf.getImageableX() + 0.5),(int)(pf.getImageableY() + 0.5));				g2.scale(shrink, shrink);			}			else {  //print at current size using page setup information				if (dimension == null) {					dimension = component.getSize();				}				double dPagesWide =  (dimension.width * scale)/pf.getImageableWidth();				double dPagesHigh =  (dimension.height * scale)/pf.getImageableHeight();				int pagesWide = (int)dPagesWide;				int pagesHigh = (int)dPagesHigh;								if ((double)pagesWide != dPagesWide)					pagesWide++;				if ((double)pagesHigh != dPagesHigh)					pagesHigh++;								if (pi >= pagesWide*pagesHigh) {				    return Printable.NO_SUCH_PAGE;				}								double dEffectivePageWidth = pf.getImageableWidth()/scale;				double dEffectivePageHeight = pf.getImageableHeight()/scale;				int piW = pi % pagesWide;				int piH = pi / pagesWide;								g2.scale(1/scale, 1/at.getScaleY());  //why this is & following line are needed isn't clear, but without it OS X 10.2 crops instead of scales				g2.scale(scale, at.getScaleY());//why this is & preceding line are needed isn't clear, but without it OS X 10.2 crops instead of scales				g2.translate((int)(-dEffectivePageWidth*piW + pf.getImageableX() + 0.5),(int)(-dEffectivePageHeight*piH + pf.getImageableY() + 0.5));			}						component.printAll((Graphics2D) g);		        return Printable.PAGE_EXISTS;		    }	    }}