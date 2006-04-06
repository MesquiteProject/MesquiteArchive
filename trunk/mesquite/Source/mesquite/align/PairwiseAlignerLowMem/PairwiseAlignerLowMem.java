package mesquite.align.PairwiseAlignerLowMem;

import mesquite.align.lib.AlignmentHelperQuadraticSpace;
import mesquite.align.lib.PairwiseAlignerBasicHelper;
import mesquite.align.lib.AlignmentHelperLinearSpace;
import mesquite.align.lib.PairwiseAlignerLowMemHelper;
import mesquite.align.lib.TwoSequenceAligner;
import mesquite.categ.lib.CategoricalState;
import mesquite.lib.CommandRecord;
import mesquite.lib.Debugg;
import mesquite.lib.MesquiteInteger;
import mesquite.lib.MesquiteNumber;


public class PairwiseAlignerLowMem extends TwoSequenceAligner {
	
	AlignmentHelperLinearSpace helper;
	PairwiseAlignerBasicHelper pa;
	boolean isMinimize = true;	
    //	internal representation of gap costs is that first gap char costs gapOpen + gapExtend, and each additional character costs gapExtend 
	
	
	public boolean startJob(String arguments, Object condition, CommandRecord commandRec, boolean hiredByName) {
		//TODO gather preferences
		return true;
	}

	public String getName() {
		return "Low Memory Pairwise Sequence Aligner";
	}

	/** returns whether this module is requesting to appear as a primary choice */
   	public boolean requestPrimaryChoice(){
   		return true;  
   	}

	
	public long[][] alignSequences(long[] sequence1, long[] sequence2, boolean returnAlignment, MesquiteNumber score, CommandRecord commandRec) {

		pa = new PairwiseAlignerBasicHelper(sequence1,sequence2,false);
		helper = pa.getLinearSpaceAlignmentHelper(false, !returnAlignment);
		
		if (!returnAlignment) { // no need to  recover the alignment
			helper.fillForward(0,pa.lengthB,0,pa.lengthA,helper.noGap);			
			int myScore = Math.min(helper.fH[pa.lengthA], Math.min (helper.fD[pa.lengthA], helper.fV[pa.lengthA])) ;
			if (score != null)
				score.setValue( myScore );
			return null;  // no alignment
		} else {
			PairwiseAlignerLowMemHelper lowmemHelper =new PairwiseAlignerLowMemHelper(pa,helper); 
			int myScore = lowmemHelper.recursivelyFillArray(0, pa.lengthA, 0, pa.lengthB, helper.noGap, helper.noGap);			
			if (score != null)
				score.setValue(myScore);
			return lowmemHelper.recoverAlignment();
		}
	}		

}
