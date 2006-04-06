package mesquite.align.lib;

import mesquite.categ.lib.CategoricalState;
import mesquite.lib.CommandRecord;
import mesquite.lib.MesquiteModule;
import mesquite.lib.MesquiteNumber;

public class PairwiseAlignerBasicHelper  {

	boolean preferencesSet = false;
	boolean isMinimize = true;	

	
    //	internal representation of gap costs is that first gap char costs gapOpen + gapExtend, and each additional character costs gapExtend 
	public int gapOpen;
	public int gapExtend;
	boolean keepGaps;
	public int totalGapChars = 0;
	public int subs[][];
	
	public int A[];
	public int B[];
	public int followsGapSize[] = null; // unused if "keepGaps" is false
	public int lengthA;
	public int lengthB;
	
	public boolean seqsWereExchanged=false;
	
	public PairwiseAlignerBasicHelper (long[] A_withGaps, long[] B_withGaps, boolean keepGaps) {
		this.keepGaps = keepGaps;
		getCosts();
		totalGapChars = stripGaps(A_withGaps, B_withGaps, keepGaps);	
	}

	public long[][] alignSequences( boolean returnAlignment, MesquiteNumber score) {

		if ( returnAlignment) { // fast (but quadratic space) alignment 
			AlignmentHelperQuadraticSpace helper = new AlignmentHelperQuadraticSpace(A, B, lengthA, lengthB, subs, gapOpen, gapExtend);
			return helper.doAlignment(returnAlignment,score,keepGaps, followsGapSize, totalGapChars);
		} else { //linear space	

			AlignmentHelperLinearSpace helper = getLinearSpaceAlignmentHelper(keepGaps, true);
			helper.fillForward(0,lengthB,0,lengthA,helper.noGap);			
			int myScore = Math.min(helper.fH[lengthA], Math.min (helper.fD[lengthA], helper.fV[lengthA])) ;
			
			if (score != null)
				score.setValue( myScore );

			return null;  // no alignment
		}
	}
	

	public int stripGaps (long[] A_withGaps, long[] B_withGaps, boolean keepGaps) {
		int i;
		int totalGapChars = 0;
		
		A = new int[A_withGaps.length];
		B = new int[B_withGaps.length];
		
		if (keepGaps) { // only do this in case where gaps in A are tracked
			followsGapSize = new int[A_withGaps.length];
			//translate sequences to ints, and remove gaps
			for (i=0; i<A_withGaps.length; i++) 
				followsGapSize[i] = 0;
		}
		
		//translate sequences to ints, and remove gaps
		for (i=0; i<A_withGaps.length; i++) {
			if (!CategoricalState.isInapplicable(A_withGaps[i])) {
				int state = (int)CategoricalState.getOnlyElement(A_withGaps[i]);
				if (state >=0)  //Travis: added this to protect against ambiguity codes; ideally do this better
					A[lengthA]= state;
				lengthA++;
			}	else if (keepGaps) {
				followsGapSize[lengthA]++;
				totalGapChars++;	
			}
		}

		for (i=0; i<B_withGaps.length; i++) {
			if (!CategoricalState.isInapplicable(B_withGaps[i])) {
				int state = (int)CategoricalState.getOnlyElement(B_withGaps[i]);
				if (state>=0)  //Travis: added this to protect against ambiguity codes; ideally do this better
						B[lengthB]= state;
				lengthB++;
			}
		}	
		
		
		//Just to decrease the amount of memory used to O(min of A_length and B_length)
		//note: the "withGaps" version won't do this (at least at first), since then the "A" that needs to have gaps retained would be wrong. 
		if ( !keepGaps && lengthB < lengthA) {
			int [] tmp2 = A;
			A = B;
			B = tmp2;

			int tmp = lengthA;
			lengthA = lengthB;
			lengthB = tmp;
			
			seqsWereExchanged = true;
		}		
		
		return totalGapChars;
	}
	
	private boolean getCosts(){
	    //	internal representation of gap costs is that first gap char costs gapOpen+gapExtend, and each additional character costs gapExtend
		gapOpen = 2;
		gapExtend = 2;
		
		return getCostMatrix();
		
	}

	private boolean getCostMatrix(){
		int alphabetLength = 4;
		subs = subs =  new int[alphabetLength][alphabetLength];
		subs[0][0] = subs[1][1] = subs[2][2] = subs[3][3] = 0;
		subs[0][1] = subs[1][0] = subs[2][3] = subs[3][2] = 1;
		subs[0][2] = subs[2][0] = subs[1][2] = subs[2][1] = 2;
		subs[0][3] = subs[3][0] = subs[1][3] = subs[3][1] = 2;
		return true;
	}
	
	public AlignmentHelperLinearSpace getLinearSpaceAlignmentHelper (boolean keepGaps, boolean justForwardArray) {
		return new AlignmentHelperLinearSpace(A, B, lengthA, lengthB, subs, gapOpen, gapExtend, justForwardArray);
	}
}
