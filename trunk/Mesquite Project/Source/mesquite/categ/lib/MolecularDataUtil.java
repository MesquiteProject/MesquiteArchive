package mesquite.categ.lib;

import java.awt.Rectangle;

import mesquite.basic.ManageSetsBlock.ManageSetsBlock;
import mesquite.lib.*;
import mesquite.lib.characters.CodonPositionsSet;
import mesquite.lib.duties.ElementManager;
import mesquite.lists.lib.ListModule;
import mesquite.align.lib.*;

/** This is a utility class that provides static methods to do various jobs with molecular data  */

public class MolecularDataUtil {

	/*.................................................................................................................*
	private boolean alignTouchedToDropped(int rowToAlign, int recipientRow){
		MesquiteNumber score = new MesquiteNumber();
		if (aligner==null) {
			aligner = new PairwiseAligner(true,allowNewGaps.getValue(), subs,gapOpen.getValue(), gapExtend.getValue(), gapOpenTerminal.getValue(), gapExtendTerminal.getValue(), alphabetLength);
			//aligner.setUseLowMem(true);
		}
		if (aligner!=null){
			//aligner.setUseLowMem(data.getNumChars()>aligner.getCharThresholdForLowMemory());
			originalCheckSum = ((CategoricalData)data).storeCheckSum(0, data.getNumChars()-1,rowToAlign, rowToAlign);
			aligner.setAllowNewInternalGaps(allowNewGaps.getValue());
			long[][] aligned = aligner.alignSequences((MCategoricalDistribution)data.getMCharactersDistribution(), recipientRow, rowToAlign,MesquiteInteger.unassigned,MesquiteInteger.unassigned,true,score);
			if (aligned==null) {
				logln("Alignment failed!");
				return false;
			}
			logln("Align " + (rowToAlign+1) + " onto " + (recipientRow+1));
			long[] newAlignment = Long2DArray.extractRow(aligned,1);

			int[] newGaps = aligner.getGapInsertionArray();
			if (newGaps!=null)
				alignUtil.insertNewGaps((MolecularData)data, newGaps);
			Rectangle problem = alignUtil.forceAlignment((MolecularData)data, 0, data.getNumChars()-1, rowToAlign, rowToAlign, 1, aligned);

			((CategoricalData)data).examineCheckSum(0, data.getNumChars()-1,rowToAlign, rowToAlign, "Bad checksum; alignment has inapproppriately altered data!", warnCheckSum, originalCheckSum);
			return true;
		}
		return false;
	}
	/*.................................................................................................................*/
	public static void pairwiseAlignMatrix(MesquiteModule module, MolecularData data, int referenceTaxon, int itStart, int itEnd, boolean allowNewGaps) {
		MesquiteNumber score = new MesquiteNumber();
		AlignUtil alignUtil = new AlignUtil();
		PairwiseAligner aligner = PairwiseAligner.getDefaultAligner(true,data);
		aligner.setAllowNewInternalGaps(allowNewGaps);
		for (int it=itStart; it<=itEnd; it++) {
			if (it!=referenceTaxon) {
				long[][] aligned = aligner.alignSequences((MCategoricalDistribution)data.getMCharactersDistribution(), referenceTaxon, it,MesquiteInteger.unassigned,MesquiteInteger.unassigned,true,score);
				int[] newGaps = aligner.getGapInsertionArray();
				if (newGaps!=null)
					alignUtil.insertNewGaps((MolecularData)data, newGaps);
				Rectangle problem = alignUtil.forceAlignment((MolecularData)data, 0, data.getNumChars()-1, it, it, 1, aligned);
			}
		}

	}
	/*.................................................................................................................*/
	public static void pairwiseAlignMatrix(MesquiteModule module, MolecularData data, int referenceTaxon, boolean allowNewGaps) {
		pairwiseAlignMatrix(module, data,referenceTaxon, 0, data.getNumTaxa()-1, allowNewGaps);
	}
	/*.................................................................................................................*/
	public static void shiftAlignTaxon(MolecularData data, int referenceTaxon, int taxonToAlign) {
		pairwiseAlignMatrix(null, data, referenceTaxon, taxonToAlign, taxonToAlign, false);
	}

	/*.................................................................................................................*/
	public static void setCodonPositions(DNAData data, CodonPositionsSet modelSet,  int position,  boolean calc, boolean notify){
		if (modelSet==null)
			return;
		MesquiteNumber num = new MesquiteNumber();
		num.setValue(position);
		if (modelSet != null) {
			for (int i=0; i<data.getNumChars(); i++) {
				modelSet.setValue(i, num);
				if (calc) {
					num.setValue(num.getIntValue()+1);
					if (num.getIntValue()>3)
						num.setValue(1);
				}
			}


		}
	}
	/*.................................................................................................................*/
	public static int getMinimumStops(DNAData data, int it, CodonPositionsSet modelSet){
		int minStops = -1;
		for (int i = 1; i<=3; i++) {
			setCodonPositions(data,modelSet, i,true,false);  //set them temporarily
			int totNumStops = ((DNAData)data).getAminoAcidNumbers(it,ProteinData.TER);					 
			if (minStops<0 || totNumStops<minStops) {
				minStops = totNumStops;
			}
		}
		return minStops;
	}
	/*.................................................................................................................*
	public static int getShiftForMinimumStops(DNAData data, int it, CodonPositionsSet modelSet){
		int minStops = -1;
		for (int i = 1; i<=3; i++) {
			setCodonPositions(data,modelSet, i,true,false);  //set them temporarily
			int totNumStops = ((DNAData)data).getAminoAcidNumbers(it,ProteinData.TER);					 
			if (minStops<0 || totNumStops<minStops) {
				minStops = totNumStops;
			}
		}
		return minStops;
	}

	/*.................................................................................................................*/
	protected static double alignmentScoreRatioToRCScore(DNAData data, int it1, int it2) {
		PairwiseAligner aligner = PairwiseAligner.getDefaultAligner(false,data);
		if (aligner==null)
			return 0.0;
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
		aligner.alignSequences(extracted1, extracted2, false, alignScore);

		for (int ic = firstSite; ic<=lastSite; ic++){
			extracted2[lastSite-ic] = DNAState.complement(data.getState(ic, it2));
		}
		MesquiteNumber alignRCScore = new MesquiteNumber();
		aligner.alignSequences(extracted1, extracted2, false, alignRCScore);
		alignScore.divideBy(alignRCScore);
		return alignScore.getDoubleValue();

   	 }

	/*.................................................................................................................*/
	public static void reverseComplementSequencesIfNecessary(DNAData data, MesquiteModule module, Taxa taxa, int itStart, int itEnd, boolean baseOnStopCodons) {

		if (baseOnStopCodons) {
			CodonPositionsSet modelSet = (CodonPositionsSet) data.getCurrentSpecsSet(CodonPositionsSet.class);
			if (modelSet == null) {
				modelSet= new CodonPositionsSet("Codon Positions", data.getNumChars(), data);
				modelSet.addToFile(data.getFile(), module.getProject(), module.findElementManager(CodonPositionsSet.class)); //THIS
				data.setCurrentSpecsSet(modelSet, CodonPositionsSet.class);
			}
			for (int it = itStart; it<taxa.getNumTaxa() && it<=itEnd; it++) {
				int stops = getMinimumStops(data, it, modelSet);
				data.reverseComplement(0, data.getNumChars()-1, it, false, false);
				int stopsRC = getMinimumStops(data, it, modelSet);
				if (stops<=stopsRC) {
					data.reverseComplement(0, data.getNumChars(), it, false, false);  // then we need to reverse them back.
				} else
					module.logln("Reverse complemented sequence " + (it+1));
			}
		} else {
			for (int it = itStart; it<taxa.getNumTaxa() && it<=itEnd; it++) {
				double score = alignmentScoreRatioToRCScore((DNAData)data, 0, it);
				if (score>1.0){
					data.reverseComplement(0, data.getNumChars(), it, false, false);  // then we need to reverse them back.
					module.logln("Reverse complemented sequence " + (it+1));
				}
			}
		
		}
	}
	
	/*.................................................................................................................*/
	public static int numApplicable(long[] sequence) {
		int count = 0;
		for (int i=0; i<sequence.length; i++) {
			if (!DNAState.isInapplicable(sequence[i]))
				count++;
		}
		return count;
	}
	/*.................................................................................................................*/
	public static double numNucleotide(long[] sequence, int nucleotideState) {
		double count = 0.0;
		for (int i=0; i<sequence.length; i++) {
			if (DNAState.isElement(sequence[i], nucleotideState)) {
				count+= 1.0/DNAState.cardinality(sequence[i]);
			}
		}
		return count;
	}
	/*.................................................................................................................*/
	public static double numA(long[] sequence) {
		return numNucleotide(sequence,0);
	}
	/*.................................................................................................................*/
	public static double numC(long[] sequence) {
		return numNucleotide(sequence,1);
	}
	/*.................................................................................................................*/
	public static double numG(long[] sequence) {
		return numNucleotide(sequence,2);
	}
	/*.................................................................................................................*/
	public static double numT(long[] sequence) {
		return numNucleotide(sequence,3);
	}
	/** calculates melting temperature of sequence based upon simple formula presented on
	 * http://www.promega.com/techserv/tools/biomath/calc11.htm#disc
	 * which is cited as being from
	 * Rychlik, W. and Rhoads, R.E. (1989) Nucl. Acids Res. 17, 8543.
	/*.................................................................................................................*/
	public static double getMeltingTemperature(long[] sequence) {  
		double tm=0.0;
		int numSites = numApplicable(sequence);
		if (numSites==0) return 0.0;
		if (numSites<14) {
			tm =  4 * (numC(sequence) + numG(sequence)) + 2 * (numA(sequence) + numT(sequence)) ;
		} else {
			tm =  64.9 + 41 * (numC(sequence) + numG(sequence) - 16.4)/numSites;
		}
		return tm;
	}

	
	/*.................................................................................................................*/
	public static void setCodonPositions(DNAData data, MesquiteModule module, Taxa taxa, int itStart, int itEnd, int startPos) {

		CodonPositionsSet modelSet = (CodonPositionsSet) data.getCurrentSpecsSet(CodonPositionsSet.class);
		if (modelSet == null) {
			modelSet= new CodonPositionsSet("Codon Positions", data.getNumChars(), data);
			modelSet.addToFile(data.getFile(), module.getProject(), module.findElementManager(CodonPositionsSet.class)); //THIS
			data.setCurrentSpecsSet(modelSet, CodonPositionsSet.class);
		} 
		setCodonPositions(data,modelSet, startPos,true,false);  
	}
	/*.................................................................................................................*/
	public static void setCodonPositionsToMinimizeStops(DNAData data, MesquiteModule module, Taxa taxa, int itStart, int itEnd) {

		CodonPositionsSet modelSet = (CodonPositionsSet) data.getCurrentSpecsSet(CodonPositionsSet.class);
		if (modelSet == null) {
			modelSet= new CodonPositionsSet("Codon Positions", data.getNumChars(), data);
			modelSet.addToFile(data.getFile(), module.getProject(), module.findElementManager(CodonPositionsSet.class)); //THIS
			data.setCurrentSpecsSet(modelSet, CodonPositionsSet.class);
		} 
		int minStops = -1;
		int minPosStart = -1;
		for (int i = 1; i<=3; i++) {
			setCodonPositions(data,modelSet, i,true,false);  //set them temporarily
			int totNumStops = 0;
			for (int it=itStart; it<=itEnd; it++)
				totNumStops+=((DNAData)data).getAminoAcidNumbers(it,ProteinData.TER);					 
			if (minStops<0 || totNumStops<minStops) {
				minStops = totNumStops;
				minPosStart=i;
			}
		}
		if (minPosStart>=0)
			setCodonPositions(data,modelSet, minPosStart,true,false);  
	}
	/*.................................................................................................................*/
	public static void shiftToMinimizeStops(DNAData data, MesquiteModule module, Taxa taxa, int itStart, int itEnd) {

		CodonPositionsSet modelSet = (CodonPositionsSet) data.getCurrentSpecsSet(CodonPositionsSet.class);
		if (modelSet == null) {
			modelSet= new CodonPositionsSet("Codon Positions", data.getNumChars(), data);
			modelSet.addToFile(data.getFile(), module.getProject(), module.findElementManager(CodonPositionsSet.class)); //THIS
			data.setCurrentSpecsSet(modelSet, CodonPositionsSet.class);
		}

		for (int it = itStart; it<taxa.getNumTaxa() && it<itEnd; it++) {
			int stops = getMinimumStops(data, it, modelSet);
		}
	}

}
