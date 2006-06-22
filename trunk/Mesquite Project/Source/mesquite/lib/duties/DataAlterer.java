/* Mesquite source code.  Copyright 1997-2006 W. Maddison and D. Maddison.Version 1.11, June 2006.Disclaimer:  The Mesquite source code is lengthy and we are few.  There are no doubt inefficiencies and goofs in this code. The commenting leaves much to be desired. Please approach this source code with the spirit of helping out.Perhaps with your help we can be more than a few, and make Mesquite better.Mesquite is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY.Mesquite's web site is http://mesquiteproject.orgThis source code and its compiled class files are free and modifiable under the terms of GNU Lesser General Public License.  (http://www.gnu.org/copyleft/lesser.html)*/package mesquite.lib.duties;import java.awt.*;import mesquite.lib.*;import mesquite.lib.table.*;import mesquite.lib.characters.*;/* ======================================================================== *//**This is superclass of modules to alter a data matrix.*/public abstract class DataAlterer extends MesquiteModule  {   	 public Class getDutyClass() {   	 	return DataAlterer.class;   	 } 	public String getDutyName() { 		return "Data Alterer";   	}   	   	/** if returns true, then requests to remain on even after operateData is called.  Default is false*/   	public boolean pleaseLeaveMeOn(){   		return false;   	}	/*.................................................................................................................*/   	/** Called to alter data in those cells selected in table.  Returns true if data altered*/   	public abstract boolean alterData(mesquite.lib.characters.CharacterData data, MesquiteTable table, CommandRecord commandRec);   		/*.................................................................................................................*/   	/** Called to alter the data in a single cell.  If you use the alterContentOfCells method of this class,    	then you must supply a real method for this, not just this stub. */   	public void alterCell(mesquite.lib.characters.CharacterData data, int ic, int it, CommandRecord commandRec){   	}   		/*.................................................................................................................*/   	/** Called to alter data in cells in table. This is used if the altering procedure can be done on one cell   	at a time, independent of all other cells.  If the altering procedure involves dependencies between cells,   	then a different method must be built.  */   	public boolean alterContentOfCells(mesquite.lib.characters.CharacterData data, MesquiteTable table, CommandRecord commandRec){		boolean did=false; 		if (table==null && data!=null){    // alter entire matrix			for (int i=0; i<data.getNumChars(); i++)				for (int j=0; j<data.getNumTaxa(); j++) {						alterCell(data,i,j, commandRec);				}			return true; 		} 		else if (table!=null && data !=null){   	 		boolean[][] done = new boolean[table.getNumColumns()][table.getNumRows()];   			if (table.anyCellSelected()) {				for (int i=0; i<table.getNumColumns(); i++)					for (int j=0; j<table.getNumRows(); j++)						if (!done[i][j] && table.isCellSelected(i,j)) {							alterCell(data,i,j,  commandRec);							did = true;							done[i][j]=true;						}			}			if (table.anyRowSelected()) {				for (int j=0; j<table.getNumRows(); j++) {					if (table.isRowSelected(j))						for (int i=0; i<table.getNumColumns(); i++) {							if (!done[i][j]) {								alterCell(data,i,j, commandRec);								done[i][j]=true;							}							did = true;						}				}			}			if (table.anyColumnSelected()) {				for (int i=0; i<table.getNumColumns(); i++){					if (table.isColumnSelected(i))						for (int j=0; j<table.getNumRows(); j++) {							if (!done[i][j]) {								alterCell(data,i,j, commandRec);								done[i][j]=true;							}							did=true;						}				}			}			if (!table.anythingSelected()) {				for (int i=0; i<data.getNumChars(); i++)					for (int j=0; j<data.getNumTaxa(); j++) {							if (!done[i][j])								alterCell(data,i,j, commandRec);					}				return true;			}		}		return did;   	}	/*.................................................................................................................*/	/** Returns CompatibilityTest so other modules know if this is compatible with some object. */	public CompatibilityTest getCompatibilityTest(){		return new CharacterStateTest();	}}