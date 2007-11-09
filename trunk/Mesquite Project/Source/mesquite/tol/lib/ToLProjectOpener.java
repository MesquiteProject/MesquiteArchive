/* Mesquite source code.  Copyright 1997-2003 W. Maddison and D. Maddison. Version 0.996+. August 2003.Disclaimer:  The Mesquite source code is lengthy and we are few.  There are no doubt inefficiencies and goofs in this code. The commenting leaves much to be desired. Please approach this source code with the spirit of helping out.Perhaps with your help we can be more than a few, and make Mesquite better.Mesquite is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY.Mesquite's web site is http://mesquiteproject.orgThis source code and its compiled class files are free and modifiable under the terms of GNU Lesser General Public License.  (http://www.gnu.org/copyleft/lesser.html) */package mesquite.tol.lib;/*~~  */import java.net.*;import java.util.*;import java.io.*;import mesquite.lib.*;import mesquite.lib.duties.*;import org.jdom.*;import org.jdom.input.SAXBuilder;//import org.tolweb.treegrow.tree.*;import org.tolweb.treegrow.main.*;//import org.jaxen.jdom.JDOMXPath;//import org.jaxen.XPath;public class ToLProjectOpener  {	public static final String GROUP_SEARCH_SRV_URL = "http://tolweb.org/onlinecontributors/app?service=external&page=xml/GroupSearchService&group=";	public static final String TREE_STRUCTURE_SRV_URL = "http://newsystem.tolweb.org/onlinecontributors/app?service=external&page=xml/TreeStructureService";		/*.................................................................................................................*/	public MesquiteProject establishProject(MesquiteModule ownerModule, String arguments, int pageDepth){		FileCoordinator fileCoord = ownerModule.getFileCoordinator();		MesquiteFile thisFile = new MesquiteFile();						// the tol web services requires us to look up the nodeId via the Group Search before getting the tree structure		int nodeId = retrieveNodeIdFromGroupSearchResults(arguments, ownerModule);				// if a nodeId was not found, nodeId < 0 is true		// TODO handle the case when we fail to retrieve the nodeId from the group search service 				String treeServiceURL = TREE_STRUCTURE_SRV_URL;		if (MesquiteInteger.isCombinable(pageDepth)&& pageDepth>0)			treeServiceURL += "&page_depth=" + pageDepth;		treeServiceURL += "&node_id=" +  nodeId;		MesquiteMessage.println("Request to the Tree of Life Web Project for the following URL:\n"+ treeServiceURL + "\n");				// verifyURLValid wraps the MesquiteFile call that was here previously		thisFile = verifyURLValid(treeServiceURL, arguments, ownerModule);				if (thisFile == null) {			// if verifyURLValid failed, the Exception will already be reported through the owner-module and the debugg output			MesquiteMessage.println("URL for request for tree from ToL not valid");			return null;		}					org.jdom.Element root = getDocumentRootFromRemoteURL(treeServiceURL, ownerModule);		// if call fails, the exception will likely be reported twice 		if (root == null) {			ownerModule.discreetAlert( "Sorry, no tree was obtained from the database");			return null;		}				int numTaxa = ToLUtil.countTerminals(root, "  ");		if (numTaxa == 0) {			ownerModule.discreetAlert( "Sorry, no tree was obtained from the database");			return null;		}			//looks as if tree was recovered properly; prepare project		MesquiteProject p = fileCoord.initiateProject(thisFile.getFileName(), thisFile);		MesquiteFile sf = CommandRecord.getScriptingFileS();		if (MesquiteThread.isScripting())			CommandRecord.setScriptingFileS(thisFile);		//getting taxon names & building Taxa block		String[] names= new String[numTaxa];		boolean[] leaves = new boolean[numTaxa];		boolean[] hasChildren = new boolean[numTaxa];		ToLUtil.getTerminals(root, names, leaves, hasChildren, new MesquiteString(), new MesquiteInteger(0));		TaxaManager taxaTask = (TaxaManager)ownerModule.findElementManager(Taxa.class);		Taxa taxa = taxaTask.makeNewTaxa("Taxa from ToL", numTaxa, false);		NameReference notesRef = NameReference.getNameReference("notes");//		NameReference leavesRef = NameReference.getNameReference("ToLLeaves");//		NameReference hasChildrenRef = NameReference.getNameReference("ToLHadChildren");		for (int i = 0; i<numTaxa; i++){			Taxon t = taxa.getTaxon(i);			t.setName(names[i]);			//taxa.setAnnotation(i, names[i]);							//AttachedNotesVector attachedNotes = new AttachedNotesVector(taxa);			//AttachedNote newNote = new AttachedNote();			//newNote.setAuthor("tolweb.org");			//newNote.setComment(names[i], true);			//attachedNotes.addNote(newNote, true);			taxa.setAssociatedObject(NameReference.getNameReference("ToLLeaves"), i, new MesquiteBoolean(leaves[i]));			taxa.setAssociatedObject(NameReference.getNameReference("ToLHasChildren"), i, new MesquiteBoolean(hasChildren[i]));		}		taxa.addToFile(thisFile, p, taxaTask);		//getting tree structure		MesquiteTree tree = new MesquiteTree(taxa);		ToLUtil.buildTree(true,root, tree, tree.getRoot(), names, new MesquiteInteger(0));		tree.setName("Tree for " + arguments);		TreeVector trees = new TreeVector(taxa);		trees.addElement(tree, false);		trees.addToFile(thisFile,p,ownerModule.findElementManager(TreeVector.class));			trees.setName("Trees for " + arguments);		//cleaning up and scripting the windows to show the tree		CommandRecord.setScriptingFileS(sf);		MesquiteModule treeWindowCoord = ownerModule.getFileCoordinator().findEmployeeWithName("#BasicTreeWindowCoord");		if (treeWindowCoord!=null){			String commands = "makeTreeWindow " + p.getTaxaReference(taxa) + "  #BasicTreeWindowMaker; tell It; ";			commands += "getEmployee #mesquite.tol.SearchToLTaxon.SearchToLTaxon; tell It; enableTools; endTell;";			commands += "setTreeSource  #StoredTrees; tell It; setTaxa " + p.getTaxaReference(taxa) + " ;  setTreeBlock 1; endTell; ";			commands += "getTreeDrawCoordinator #mesquite.trees.BasicTreeDrawCoordinator.BasicTreeDrawCoordinator;";			commands += "tell It; suppress; setTreeDrawer  #mesquite.trees.SquareTree.SquareTree; tell It; orientRight; endTell; desuppress; endTell;";			commands += "getWindow; tell It; setActive; setSize 600 600; getToolPalette; tell It; setTool mesquite.tol.SearchToLTaxon.SearchToLTaxonToolExtra.goToToLTaxon; endTell; endTell;";			commands += "  showWindow; endTell; ";			MesquiteInteger pos = new MesquiteInteger(0);			Puppeteer pup = new Puppeteer(ownerModule);			CommandRecord oldCR = MesquiteThread.getCurrentCommandRecord();			MesquiteThread.setCurrentCommandRecord(new CommandRecord(true));			pup.execute(treeWindowCoord, commands, pos, null, false);			MesquiteThread.setCurrentCommandRecord(oldCR);		}		return p;	}	/*--------------------------*/	private int retrieveNodeIdFromGroupSearchResults(String groupName, MesquiteModule ownerModule) {		String serviceURL = GROUP_SEARCH_SRV_URL + groupName;		MesquiteFile theFile = verifyURLValid(serviceURL, groupName, ownerModule);		if (theFile == null) {			return Integer.MIN_VALUE;		}		org.jdom.Element root = getDocumentRootFromRemoteURL(serviceURL, ownerModule);		try {			int count = Integer.parseInt(root.getAttributeValue("COUNT"));			if (count == 1) {				return Integer.parseInt(((Element)root.getChildren().get(0)).getAttributeValue("ID"));			} else {				return Integer.MIN_VALUE;			}		} catch (Exception e) {			ownerModule.discreetAlert("Sorry, the Group ID Service appears to have failed");			Debugg.println("Exception " + e);			return Integer.MIN_VALUE;		}	}	/*--------------------------*/	private MesquiteFile verifyURLValid(String url, String groupName, MesquiteModule ownerModule) {		MesquiteFile theFile = new MesquiteFile();			//the following shouldn't be needed but reflects inertia in the MesquiteFile class		try {			theFile.setLocs(false, new URL(url), groupName, null);			return theFile;		}		catch (MalformedURLException e){			ownerModule.discreetAlert( "Sorry, the URL appears malformed");			return null;		}			}	/*--------------------------*/	private org.jdom.Element getDocumentRootFromRemoteURL(String url, MesquiteModule ownerModule) {		//preparing XML parsing		SAXBuilder saxBuilder = new SAXBuilder();		Document jdomDocument;		org.jdom.Element root = null;		try {			jdomDocument = saxBuilder.build(url);			return jdomDocument.getRootElement();		}		catch (IOException e){			ownerModule.discreetAlert( "Sorry, the database was inaccessible");			Debugg.println("IOException " + e);			return null;		}		catch (JDOMException e){			ownerModule.discreetAlert( "Sorry, there has been a JDOMException");			Debugg.println("JDOMException " + e);			return null;		}			}	/*--------------------------*/}