package mesquite.align.lib;

import mesquite.categ.lib.CategoricalState;
import mesquite.lib.CommandRecord;
import mesquite.lib.MesquiteLong;
import mesquite.lib.MesquiteModule;
import mesquite.lib.MesquiteNumber;

public class PairwiseAlignerLowMemHelper  {

	boolean preferencesSet = false;
	boolean isMinimize = true;	

	private int lastAWhenBAligned[];
	private int shapeWhenBAligned[];

    //	internal representation of gap costs is that first gap char costs gapOpen + gapExtend, and each additional character costs gapExtend 
	PairwiseAlignerBasicHelper pa;
	AlignmentHelperLinearSpace helper;
	boolean keepGaps = false;
	boolean gapInsertionArray[];
	
	public PairwiseAlignerLowMemHelper (PairwiseAlignerBasicHelper pa, AlignmentHelperLinearSpace helper, boolean keepGaps) {
		initialize(pa,helper,keepGaps);
	}

	public PairwiseAlignerLowMemHelper (PairwiseAlignerBasicHelper pa, AlignmentHelperLinearSpace helper) {
		initialize(pa,helper,false);
	}	

	private void initialize(PairwiseAlignerBasicHelper pa, AlignmentHelperLinearSpace helper, boolean keepGaps) {
		this.pa = pa;
		this.helper = helper;
		this.keepGaps = keepGaps;
		lastAWhenBAligned = new int[pa.lengthB +1];
		shapeWhenBAligned  = new int[pa.lengthB +1];
	}
	
	public int recursivelyFillArray(int firstColumn, int lastColumn, int firstRow, int lastRow, int precedingShape, int succeedingShape) {
		
		helper.fillArrays(firstRow, lastRow, firstColumn, lastColumn, precedingShape, succeedingShape );
		
		int midRow = (firstRow + lastRow)/2;
		
		int i, verticalColScore, diagonalColScore; 
		int bestColScore = helper.bigNumber , bestCol = -1, bestColShape = helper.noGap;
		for (i=firstColumn; i<=lastColumn; i++) {			
			// best of the ways of leaving the i,j cell of the full DP table with a vertical edge 
			verticalColScore = Math.min(helper.fH[i] + helper.rH[i] + pa.gapExtend + pa.gapOpen,
							Math.min(helper.fH[i] + helper.rD[i] + pa.gapExtend + pa.gapOpen,
							Math.min(helper.fH[i] + helper.rV[i] + pa.gapExtend,
							Math.min(helper.fD[i] + helper.rH[i] + pa.gapExtend + pa.gapOpen,
							Math.min(helper.fD[i] + helper.rD[i] + pa.gapExtend + pa.gapOpen,
							Math.min(helper.fD[i] + helper.rV[i] + pa.gapExtend,			
							Math.min(helper.fV[i] + helper.rH[i] + pa.gapExtend,
							Math.min(helper.fV[i] + helper.rD[i] + pa.gapExtend,
										helper.fV[i] + helper.rV[i] + pa.gapExtend - pa.gapOpen))))))));			

			// best of the ways of leaving the i,j cell of the full DP table with a diagonal edge
			if (i == lastColumn) { 
				diagonalColScore = verticalColScore + 1;
			} else {
				int a = pa.A[i];
				int b = pa.B[midRow];
				int s = pa.subs[a][b];
				diagonalColScore = Math.min(helper.fH[i], Math.min (helper.fD[i], helper.fV[i])) +
										Math.min(helper.rH[i+1], Math.min (helper.rD[i+1], helper.rV[i+1])) +
										s;
			}
			// no need to track horizontal edges - those are already stored in the forward and reverse arrays
			
			if (verticalColScore < diagonalColScore) {
				if ( verticalColScore < bestColScore) {
					bestColScore = verticalColScore;
					bestCol = i;
					bestColShape = helper.gapInA;
				}
			} else {
				if ( diagonalColScore < bestColScore) {
					bestColScore = diagonalColScore;
					bestCol = i;
					bestColShape = helper.noGap;
				}				
			}
		}
		
		lastAWhenBAligned[midRow] = bestCol;
		shapeWhenBAligned[midRow] = bestColShape;
				
		//	Recurse to find the full list of cells through which the alignment passes.
		if ( firstRow != midRow){
			recursivelyFillArray(firstColumn, bestCol, firstRow, midRow, precedingShape, bestColShape);
		}
		if (midRow+1 != lastRow){
			if (bestColShape == helper.noGap){
				recursivelyFillArray(bestCol+1, lastColumn, midRow+1, lastRow, bestColShape, succeedingShape); 
			} else {//gapInA
				recursivelyFillArray(bestCol, lastColumn, midRow+1, lastRow, bestColShape, succeedingShape);
			}		
		}
		
		return bestColScore; // useful for the first level of recursion - returns the total alignment cost
	}
	
	
	public long[][] recoverAlignment () {
		int i=0,j, k=0;
		
		long[][] alignment = new long[pa.lengthA + pa.lengthB][2];
		for (j = 0; j<pa.lengthB; j++) {

			//gap in B
			while (i < lastAWhenBAligned[j]) {
				alignment[k][0] = CategoricalState.makeSet(pa.A[i]);
				alignment[k][1] = CategoricalState.inapplicable;
				i++;
				k++;					
			}

			//now we're ready to burn off a letter from B, and possibly a letter from A if diagonal.
			if (shapeWhenBAligned[j] == helper.noGap) {
				alignment[k][0] = CategoricalState.makeSet(pa.A[i]);
				alignment[k][1] = CategoricalState.makeSet(pa.B[j]);
				i++;
			} else {		// gapInA								
				alignment[k][0] = CategoricalState.inapplicable;
				alignment[k][1] = CategoricalState.makeSet(pa.B[j]);
			}
			k++;
		}

		while (i < pa.lengthA) {
			alignment[k][0] = CategoricalState.makeSet(pa.A[i]);
			alignment[k][1] = CategoricalState.inapplicable;
			i++;
			k++;					
		}		
		
		//trim off all the empty space at the end
		long seq2return[][] = new long[k][2];
		int ii=0; 
		
		for (i=0; i<k; i++) {
			if (pa.seqsWereExchanged) {//exhange the sequences
				seq2return[i][0] = alignment[i][1];
				seq2return[i][1] = alignment[i][0];				
			} else {
				seq2return[i][0] = alignment[i][0];
				seq2return[i][1] = alignment[i][1];
			}
		}
				
		

		if (keepGaps) {
			long finalSeq2return[][] = new long[k+pa.totalGapChars][2];
			gapInsertionArray = new boolean[k+pa.totalGapChars];
			for(i=0; i<k+pa.totalGapChars; i++) {
				gapInsertionArray[i] =false;
			}
			
			int usedGaps=0;
			int recentGapRunLength=0;
			j=0; // counts the number of letters in A seen so far
			for (i=0; i<k; i++) {
				if(seq2return[i][0] == CategoricalState.inapplicable) {
					recentGapRunLength++;
					gapInsertionArray[i+usedGaps]=true;
				} else {
					for (int m=0 ; m < pa.followsGapSize[j]-recentGapRunLength; m++){
						finalSeq2return[i+usedGaps][0] =  CategoricalState.inapplicable; 
						finalSeq2return[i+usedGaps][1] =  CategoricalState.inapplicable; 
						usedGaps++;
					}
					j++;
					recentGapRunLength=0;
				}				
				finalSeq2return[i+usedGaps][0] = seq2return[i][0] ;
				finalSeq2return[i+usedGaps][1] = seq2return[i][1] ;									
			}		
			return finalSeq2return;
		} 
		
		return seq2return;
	}


	//for now, there's not a good way to get to this array; I'll add it when this interface gets worked out
	public boolean[] getGapInsertionArray () {
		return gapInsertionArray;
	}	
}
