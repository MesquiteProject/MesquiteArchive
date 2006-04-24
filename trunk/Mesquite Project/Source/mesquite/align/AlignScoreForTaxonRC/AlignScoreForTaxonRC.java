package mesquite.align.AlignScoreForTaxonRC;

import mesquite.align.lib.*;
import mesquite.categ.lib.*;
import mesquite.lib.*;
import mesquite.lib.characters.*;
import mesquite.lib.duties.*;

/* TODO: 
 *   have it so the reference sequence can be chosen
 *   report pairwise alignment parameters
 *   snapshot, etc.
 */

public class AlignScoreForTaxonRC extends AlignScoreForTaxonGen {
	/*.................................................................................................................*/
	protected  void getAlignmentScore(DNAData data, MCategoricalDistribution observedStates, int it1, int it2, MesquiteNumber score, CommandRecord commandRec) {
		if (aligner==null)
			return;
		int firstSite = 0;
		int lastSite = data.getNumChars()-1;
		int numChars = lastSite - firstSite+1;
		
		long[] extracted1 = new long[numChars];
		long[] extracted2 = new long[numChars];
		
		for (int ic = firstSite; ic<=lastSite; ic++){
			extracted1[ic] = data.getState(ic, it1);
			extracted2[ic] = data.getState(ic, it2);
		}
		MesquiteNumber alignScore = new MesquiteNumber();
  		CategoricalState state = data.getNewState();
		aligner.alignSequences(extracted1, extracted2, false, alignScore);

		for (int ic = firstSite; ic<=lastSite; ic++){
			extracted1[lastSite-ic] = DNAData.complement(data.getState(ic, it1));
		}
		MesquiteNumber alignRCScore = new MesquiteNumber();
		aligner.alignSequences(extracted1, extracted2, false, alignRCScore);
		alignScore.divideBy(alignRCScore);
		score.setValue(alignScore);

   	 }
	/*.................................................................................................................*/
	 public String getName() {
	return "Align Score/RC Align Score To First";  
	 }
	/*.................................................................................................................*/
	 public String getScoreName() {
	return "Alignment Score / Reverse Complement Alignment Score";  
	 }
  	 
}



