package mesquite.tol.iSpeciesURLServer;

import mesquite.lib.*;
import mesquite.tol.lib.TaxonOnWebServer;

public class iSpeciesURLServer extends TaxonOnWebServer {

	public String getURL(String taxonName) {
		String openName = StringUtil.encodeForURL(StringUtil.replace(taxonName,' ', '+'));
		return "http://darwin.zoology.gla.ac.uk/~rpage/ispecies/?q=" + openName + "&submit=Go";
	
	}

	public String getName() {
		return "iSpecies";
	}
	
	public String getExplanation() {
		return "Provides a URL to a taxon in iSpecies";
	}


}
