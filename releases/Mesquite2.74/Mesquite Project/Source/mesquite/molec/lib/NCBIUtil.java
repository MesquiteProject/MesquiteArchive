/* Mesquite source code.  Copyright 1997-2010 W. Maddison and D. Maddison.
Disclaimer:  The Mesquite source code is lengthy and we are few.  There are no doubt inefficiencies and goofs in this code. 
The commenting leaves much to be desired. Please approach this source code with the spirit of helping out.
Perhaps with your help we can be more than a few, and make Mesquite better.

Mesquite is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY.
Mesquite's web site is http://mesquiteproject.org

This source code and its compiled class files are free and modifiable under the terms of 
GNU Lesser General Public License.  (http://www.gnu.org/copyleft/lesser.html)
 */
package mesquite.molec.lib; 


import java.net.*;
import java.io.*;
import mesquite.lib.*;
import mesquite.lib.characters.*;
import mesquite.categ.lib.*;
import mesquite.io.InterpretFastaDNA.InterpretFastaDNA;   //is this guaranteed to be an installed package?
import mesquite.io.InterpretFastaProtein.InterpretFastaProtein;   //is this guaranteed to be an installed package?


/* ======================================================================== */
/** A set of static methods for NCBI*/
public class NCBIUtil {
	/*.................................................................................................................*/
	public static boolean responseSaysBLASTIsReady(String response) {
		if (StringUtil.blank(response))
			return false;
		String  s = StringUtil.getAllAfterSubString(response, "QBlastInfoBegin");
		if (!StringUtil.blank(s)) {
			Parser parser = new Parser(s);
			String token = parser.getNextToken();
			while (!StringUtil.blank(token)) {
				if ("Status".equalsIgnoreCase(token)) {
					token = parser.getNextToken();  // =
					token = parser.getNextToken();
					return ("READY".equalsIgnoreCase(token));
				}
			}
		}
		return true;
	}
	/*.................................................................................................................*/
	protected static String getTaxID(String response){
		Parser parser = new Parser();
		//		String[] stringList=null;
		parser.setString(response);
		if (!parser.isXMLDocument(false))   // check if XML
			return null;
		MesquiteString nextTag = new MesquiteString();
		MesquiteString attributes = new MesquiteString();
		String tagContent = parser.getNextXMLTaggedContent(nextTag);
		while (!StringUtil.blank(nextTag.getValue())) {
			if ("eSummaryResult".equalsIgnoreCase(nextTag.getValue())) {  //make sure it has the right root tag
				parser.setString(tagContent);
				tagContent = parser.getNextXMLTaggedContent(nextTag);
				while (!StringUtil.blank(nextTag.getValue())) {
					if ("DocSum".equalsIgnoreCase(nextTag.getValue())) {
						parser.setString(tagContent);
						tagContent = parser.getNextXMLTaggedContent(nextTag, attributes);
						while (!StringUtil.blank(nextTag.getValue())) {
							if ("Item".equalsIgnoreCase(nextTag.getValue()) && Parser.attributesContains(attributes.getValue(), "Name", "TaxID")) {
								return tagContent;
							}
							tagContent = parser.getNextXMLTaggedContent(nextTag, attributes); 
						}
					}
					tagContent = parser.getNextXMLTaggedContent(nextTag); 
				}
			} 
			tagContent = parser.getNextXMLTaggedContent(nextTag);
		}
		return "";
	}
	/*.................................................................................................................*/
	public static String fetchGenBankTaxID(String id, CharacterData data,  MesquiteModule mod, boolean writeLog, StringBuffer report){ 
		try {
			URL queryURL = getFetchTaxIDAddress(id, data);
			URLConnection connection = queryURL.openConnection();
			InputStream in = connection.getInputStream();

			StringBuffer fetchBuffer = new StringBuffer();
			int c;
			int count = 0;
			while ((c = in.read()) != -1) {
				fetchBuffer.append((char) c);
				count++;
				if (count % 200==0 && writeLog && mod!=null)
					mod.log(".");
			}
			in.close();
			String taxID = getTaxID(fetchBuffer.toString());
			return taxID;

		} catch ( Exception e ){
			// give warning
			return "";
		}
	}
	/*.................................................................................................................*/
	protected static String getLineage(String response){
		Parser parser = new Parser();
		parser.setString(response);
		if (!parser.isXMLDocument(false))   // check if XML
			return null;
		MesquiteString nextTag = new MesquiteString();
		String tagContent = parser.getNextXMLTaggedContent(nextTag);
		while (!StringUtil.blank(nextTag.getValue())) {
			if ("TaxaSet".equalsIgnoreCase(nextTag.getValue())) {  //make sure it has the right root tag
				parser.setString(tagContent);
				tagContent = parser.getNextXMLTaggedContent(nextTag);
				while (!StringUtil.blank(nextTag.getValue())) {
					if ("Taxon".equalsIgnoreCase(nextTag.getValue())) {
						parser.setString(tagContent);
						tagContent = parser.getNextXMLTaggedContent(nextTag);
						while (!StringUtil.blank(nextTag.getValue())) {
							if ("Lineage".equalsIgnoreCase(nextTag.getValue())) {
								return tagContent;
							}
							tagContent = parser.getNextXMLTaggedContent(nextTag); 
						}
					}
					tagContent = parser.getNextXMLTaggedContent(nextTag); 
				}
			} 
			tagContent = parser.getNextXMLTaggedContent(nextTag);
		}
		return "";
	}
	/*.................................................................................................................*/
	public static String fetchGenBankTaxonomy(String id, CharacterData data,  MesquiteModule mod, boolean writeLog, StringBuffer report){ 
		try {
			URL queryURL = getFetchTaxonomyAddress(id, data);
			URLConnection connection = queryURL.openConnection();
			InputStream in = connection.getInputStream();

			StringBuffer fetchBuffer = new StringBuffer();
			int c;
			int count = 0;
			while ((c = in.read()) != -1) {
				fetchBuffer.append((char) c);
				count++;
				if (count % 200==0 && writeLog && mod!=null)
					mod.log(".");
			}
			in.close();
			return getLineage(fetchBuffer.toString());
		} catch ( Exception e ){
			// give warning
			return "";
		}
	}
	/*.................................................................................................................*/
	public static String fetchTaxonomyList(String accession, CharacterData data,   MesquiteModule mod, boolean writeLog, StringBuffer report){ 
		String id = getGenBankID(accession, data instanceof DNAData,  null, false);
		if (StringUtil.blank(id))
			return "";
		String taxID = fetchGenBankTaxID(id,data,  mod, writeLog, report);
		if (StringUtil.blank(taxID))
			return "";
		return fetchGenBankTaxonomy(taxID,data,  mod, writeLog, report);
	}
	/*.................................................................................................................*/
	public static String getPutQueryURL(CharacterData data, int it, int icStart, int icEnd, int maxHits,  StringBuffer report){
		if (data==null)
			return null;
		StringBuffer searchBuffer = new StringBuffer(data.getNumChars());
		String s = data.getTaxa().getTaxonName(it);
		if (!StringUtil.blank(s))
			searchBuffer.append("%3E" + StringUtil.encodeForURL(s) + "%0D%0A");  //to make it a FASTA format
		for (int ic = icStart; ic<=icEnd; ic++) {
			data.statesIntoStringBuffer(ic, it, searchBuffer, false, false, false);
		}
		String seq = searchBuffer.toString();
		if (!StringUtil.blank(seq) && (seq.length()>49)) {
			String url = "http://www.ncbi.nlm.nih.gov/blast/Blast.cgi?"+ getMesquiteGenBankURLMarker();
			url += "&DATABASE=nr&FORMAT_TYPE=HTML";
			if (data instanceof ProteinData) 
				url += "&PROGRAM=blastp";
			else
				url += "&PROGRAM=blastn";
			url += "&CLIENT=web&SERVICE=plain&PAGE=";
			if (isNucleotides(data))
				url+="Nucleotides";
			else
				url+="Protein";
			url += "&HITLIST_SIZE="+ maxHits + "&CMD=Put&QUERY=";
			return url+seq;
		}
		else {
			if (report!=null)
				report.append("Sorry, to use the BLAST search you need to have one or more regions of 50 or more nucleotides or amino acids selected.\n");
			return null;
		}
	}
	/*.................................................................................................................*/
	public static boolean isNucleotides(CharacterData data){
		return (data instanceof DNAData);
	}
	/*.................................................................................................................*/
	public static String getMesquiteGenBankURLMarker(){
		return "tool=mesquite" + MesquiteTrunk.mesquiteTrunk.getVersionInt() + "&email=info@mesquiteproject.org";
	}
	/*.................................................................................................................*/
	public static URL getFetchTaxonomyAddress(String taxid, CharacterData data)
	throws MalformedURLException {
		String query = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?"+getMesquiteGenBankURLMarker() +"&db=taxonomy&id="+taxid+"&retmode=xml";
		return new URL(query);
	}
	/*.................................................................................................................*/
	public static URL getFetchTaxIDAddress(String uid, CharacterData data)
	throws MalformedURLException {
		String query = getMesquiteGenBankURLMarker() + "&db=" ;
		if (isNucleotides(data))
			query += "nucleotide";
		else
			query += "protein";
		query += "&id="+uid+"&retmode=xml";

		return new URL("http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?" + query);
	}
	/*.................................................................................................................*/
	public static URL getFetchSequenceAddress(String uid, CharacterData data)
	throws MalformedURLException {
		String query = getMesquiteGenBankURLMarker() + "&db=" ;
		if (isNucleotides(data))
			query += "nucleotide";
		else
			query += "protein";
		query += "&id="+uid+"&rettype=fasta";

		return new URL("http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?" + query);
	}
	/*.................................................................................................................*/
	public static URL getESearchAddress(String accessionNumber, boolean nucleotides)
	throws MalformedURLException {
		String query = getMesquiteGenBankURLMarker() + "&db=" ;
		if (nucleotides)
			query += "nucleotide";
		else
			query += "protein";
		query+= "&term="+accessionNumber;

		return new URL("http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?" + query);
	}
	/*.................................................................................................................*/
	public static String[] getGIs(String eSearchResponse){
		Parser parser = new Parser();
		String[] stringList=null;
		parser.setString(eSearchResponse);
		if (!parser.isXMLDocument(false))   // check if XML
			return null;
		MesquiteString nextTag = new MesquiteString();
		String tagContent = parser.getNextXMLTaggedContent(nextTag);
		while (!StringUtil.blank(nextTag.getValue())) {
			if ("eSearchResult".equalsIgnoreCase(nextTag.getValue())) {  //make sure it has the right root tag
				parser.setString(tagContent);
				tagContent = parser.getNextXMLTaggedContent(nextTag);
				while (!StringUtil.blank(tagContent)) {
					if ("idList".equalsIgnoreCase(nextTag.getValue())) {
						int count=0;
						String idListContent=tagContent;
						parser.setString(idListContent);
						tagContent = parser.getNextXMLTaggedContent(nextTag);
						while (!StringUtil.blank(tagContent)) {
							if ("id".equalsIgnoreCase(nextTag.getValue())) {
								count++;
							}
							tagContent = parser.getNextXMLTaggedContent(nextTag);
						}
						if (count>0)
							stringList=new String[count];
						if (stringList!=null) {
							count=0;
							parser.setString(idListContent);
							tagContent = parser.getNextXMLTaggedContent(nextTag);
							while (!StringUtil.blank(tagContent)) {
								if ("id".equalsIgnoreCase(nextTag.getValue())) {
									stringList[count]=tagContent;
									count++;
								}
								tagContent = parser.getNextXMLTaggedContent(nextTag);
							}
						}
					}
					tagContent = parser.getNextXMLTaggedContent(nextTag); 
				}
				return stringList;
			} 
			tagContent = parser.getNextXMLTaggedContent(nextTag);
		}
		return null;
	}
	/*.................................................................................................................*/
	public static String[] getGenBankIDs(String[] accessionNumbers, boolean nucleotides,  MesquiteModule mod, boolean writeLog){ 
		try {
			String searchString="";
			for (int i=0; i<accessionNumbers.length;i++) {
				searchString+=accessionNumbers[i];
				if (i<accessionNumbers.length-1)
					searchString+="+OR+";
			}
			if (writeLog && mod!=null)
				mod.log(".");

			StringBuffer sb = new StringBuffer();
			URL queryURL = getESearchAddress(searchString, nucleotides);
			if (writeLog && mod!=null)
				mod.log(".");
			URLConnection connection = queryURL.openConnection();
			if (writeLog && mod!=null)
				mod.log(".");
			InputStream in = connection.getInputStream();

			int c;
			int count=0;
			while ((c = in.read()) != -1) {
				sb.append((char) c);
				count++;
				if (count % 100==0 && writeLog && mod!=null)
					mod.log(".");
			}
			in.close();

			String[] idList = getGIs(sb.toString());

			return idList;
		} catch ( Exception e ){
			// give warning
			return null;
		}

	}
	/*.................................................................................................................*/
	public static String getGenBankID(String accessionNumber, boolean nucleotides,  MesquiteModule mod, boolean writeLog){ 
		String[] accessionNumbers = new String[1];
		accessionNumbers[0]=accessionNumber;
		String[] idList = getGenBankIDs(accessionNumbers,nucleotides,mod, writeLog);
		if (idList!=null && idList.length>0)
			return idList[0];
		else
			return null;
	}
	/*.................................................................................................................*/
	public static void importFASTASequence(CharacterData data, String sequence, MesquiteModule mod,StringBuffer report){
		data.setCharNumChanging(true);
		if (data instanceof ProteinData) {
			InterpretFastaProtein importer = new InterpretFastaProtein();
			importer.readString(data,sequence);
		} else {
			InterpretFastaDNA importer = new InterpretFastaDNA();
			importer.readString(data,sequence);
		}
		data.setCharNumChanging(false);
		
		Taxa taxa = data.getTaxa();
		String s = taxa.getTaxon(taxa.getNumTaxa()-1).getName();
		if (report!=null)
			report.append("Acquired: "+s+"\n\n");
		taxa.notifyListeners(mod, new Notification(MesquiteListener.PARTS_ADDED));
		data.notifyListeners(mod, new Notification(MesquiteListener.PARTS_ADDED));
	}
	/*.................................................................................................................*/
	public static boolean fetchGenBankSequence(String id, CharacterData data,  MesquiteModule mod, boolean writeLog, StringBuffer report){ 
		try {
			URL queryURL = getFetchSequenceAddress(id, data);
			URLConnection connection = queryURL.openConnection();
			InputStream in = connection.getInputStream();

			StringBuffer fetchBuffer = new StringBuffer();
			int c;
			int count = 0;
			while ((c = in.read()) != -1) {
				fetchBuffer.append((char) c);
				count++;
				if (count % 200==0 && writeLog && mod!=null)
					mod.log(".");
			}
			in.close();
			importFASTASequence(data,fetchBuffer.toString(),  mod, report);
			return true;
		} catch ( Exception e ){
			// give warning
			return false;
		}

	}
	/*.................................................................................................................*/
	public static boolean fetchGenBankSequences(String[] idList, CharacterData data,  MesquiteModule mod, boolean writeLog, StringBuffer report){ 
		for (int i=0; i<idList.length; i++) {
			if (!StringUtil.blank(idList[i])) {

				if (writeLog && mod!=null){
					mod.log("Fetching " + idList[i]);
					mod.log(".");
				}
				if (!fetchGenBankSequence(idList[i], data,  mod, writeLog, report))
					return false;
				mod.logln("");
			}
		}
		return true;

	}
	/*.................................................................................................................*/
	public static boolean fetchGenBankSequencesFromAccessions(String[] accessionNumbers, CharacterData data,  MesquiteModule mod, boolean writeLog, StringBuffer report){ 
		String[] idList = getGenBankIDs(accessionNumbers, data instanceof DNAData,  mod, writeLog);
		if (idList==null)
			return false;
	/*	if (writeLog && mod!=null)
			for (int i=0; i<idList.length; i++) 
				if (!StringUtil.blank(idList[i])) 				
					mod.logln ("To Fetch " + idList[i]);
*/
		boolean acquired = fetchGenBankSequences(idList,data, mod, writeLog, report);
		if (mod!=null && writeLog)
			mod.logln("");
		return acquired;

	}
}





