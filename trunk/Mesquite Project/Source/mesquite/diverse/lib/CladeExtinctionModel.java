package mesquite.diverse.lib;

import mesquite.lib.MesquiteDouble;


public class CladeExtinctionModel implements DESpeciationSystem {
	
	
	private double e0;   //extinction rate in state 0
	private double s0;   //speciation rate in state 0
	private double e1;   //extinction rate in state 1
	private double s1;   //speciation rate in state 1
	private double t01;  //transition rate from state 0 to state 1
	private double t10;  //transition rate from state 1 to state 0
	

	public CladeExtinctionModel(double e0, double s0, double e1, double s1, double t01, double t10){
		this.e0 = e0;
		this.s0 = s0;
		this.e1 = e1;
		this.s1 = s1;
		this.t01 = t01;
		this.t10 = t10;
	}

	
	/*
	 *  The probability values are passed the probs array in the following order
	 *  probs[0] = E0 = extinction probability in lineage starting at t in the past at state 0
	 *  probs[1] = E1 = extinction probability in lineage starting at t in the past at state 1
	 *  probs[2] = P0 = probability of explaining the data, given system is in state 0 at time t
	 *  probs[3] = P1 = probability of explaining the data, given system in in state 1 at time t
	 * @see mesquite.correl.lib.DESystem#calculateDerivative(double, double[])
	 */
	public String toString(){
		return "Model s0=" + MesquiteDouble.toString(s0, 4) +" s1=" + MesquiteDouble.toString(s1, 4) +" e0=" + MesquiteDouble.toString(e0, 4) +" e1=" + MesquiteDouble.toString(e1, 4) +" t01=" + MesquiteDouble.toString(t01, 4) +" t10=" + MesquiteDouble.toString(t10, 4);
	}
	public double[]calculateDerivative(double t,double probs[],double[] result){
		// for clarity
		double extProb0 = probs[0];
		double extProb1 = probs[1];
		double dataProb0 = probs[2];
		double dataProb1 = probs[3];
        result[0] = -(e0+t01+s0)*extProb0 + s0*extProb0*extProb0 + e0 + t01*extProb1; 
        result[1] = -(e1+t10+s1)*extProb1 + s1*extProb1*extProb1 + e1 + t10*extProb0;            
		result[2] = -(e0+t01+s0)*dataProb0 + 2*s0*extProb0*dataProb0 + t01*dataProb1;
		result[3] = -(e1+t10+s1)*dataProb1 + 2*s1*extProb1*dataProb1 + t10*dataProb0;
		return result;
	}
    
    public void setE0(double e0){
        this.e0 = e0;
    }
    
    public void setS0(double s0){
        this.s0 = s0;
    }
    
    public void setE1(double e1){
        this.e1 = e1;
    }
    
    public void setS1(double s1){
        this.s1 = s1;
    }
    
    public void setT01(double t01){
        this.t01 = t01;
    }
    
    public void setT10(double t10){
        this.t10 = t10;
    }

	public void resetParameters(double e0, double s0, double e1, double s1, double t01, double t10){
		this.e0 = e0;
		this.s0 = s0;
		this.e1 = e1;
		this.s1 = s1;
		this.t01 = t01;
		this.t10 = t10;
	}


    public double getSRate(int state) {
        if (state == 0)
            return s0;
        else if (state == 1)
            return s1;
        else return MesquiteDouble.unassigned;
    }
	
    public double getERate(int state) {
        if (state == 0)
            return e0;
        else if (state == 1)
            return e1;
        else return MesquiteDouble.unassigned;
    }
	
}