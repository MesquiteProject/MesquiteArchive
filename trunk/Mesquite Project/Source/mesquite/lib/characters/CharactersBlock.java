/* Mesquite source code.  Copyright 1997-2007 W. Maddison and D. Maddison.Version 2.01, December 2007.Disclaimer:  The Mesquite source code is lengthy and we are few.  There are no doubt inefficiencies and goofs in this code. The commenting leaves much to be desired. Please approach this source code with the spirit of helping out.Perhaps with your help we can be more than a few, and make Mesquite better.Mesquite is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY.Mesquite's web site is http://mesquiteproject.orgThis source code and its compiled class files are free and modifiable under the terms of GNU Lesser General Public License.  (http://www.gnu.org/copyleft/lesser.html)*/package mesquite.lib.characters; import java.awt.*;import mesquite.lib.duties.*;import mesquite.lib.*;/*Last documented:  April 2003 *//*===============================================*//** A block in a NEXUS file representing a CharacterData */public class CharactersBlock extends NexusBlock {	CharacterData data = null;	public CharactersBlock(MesquiteFile f, MesquiteModule mb){		super(f, mb);	}	public boolean getWritable(){		if (data == null)			return false;		return data.getWritable();	}	public boolean contains(FileElement e) {		return e != null && data == e;	}	public boolean mustBeAfter(NexusBlock block){ 		if (block==null)			return false;		if (data!=null && block instanceof TaxaBlock) {			return data.getTaxa() == ((TaxaBlock)block).getTaxa();		}		return (block.getBlockName().equalsIgnoreCase("TAXA"));			}	public String getBlockName(){		return "CHARACTERS";	}	public void setData(CharacterData data) {		this.data = data;	}	public CharacterData getData() {		return data;	}	public void written() {		data.setDirty(false);	}	public String getName(){		if (data==null)			return "empty characters block";		else			return "Characters block: " + data.getName();	}	public void writeNEXUSBlock(MesquiteFile file, ProgressIndicator progIndicator){		if (data==null)			return;		if (data.getMatrixManager()!=null) {			data.getMatrixManager().writeCharactersBlock(data, this, file, progIndicator);			data.resetChangedSinceSave();		}	}}	