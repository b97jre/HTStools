package general;

import java.util.ArrayList;
import java.util.Hashtable;

import flanagan.analysis.Regression;

public class Statistics {

	public static double[] calculateAUC(double[][] ROCs){
		double[] AUC = null;

		for(int i = 0; i < ROCs.length; i++){
			if(ROCs[i] != null){
				//System.out.println("nr of genomes : "+  (i+1));
				AUC = Functions.addDouble(AUC, calculateAUC(ROCs[i]));
			}
		}
		return AUC;
	}
	
	public static double calculateAUC(double[] ROCs){
		double AUC = 0;
		ROCs = Functions.sortScores(ROCs);
		//System.out.println("nr of interactions : "+  ROCs.length);
		for(int i = 0; i < ROCs.length; i++){
			AUC += ((double)(1)/(double)(ROCs.length))*ROCs[i];
		}
		return AUC;
	}
	
	public static double[][] getScores(ArrayList<double[][]> targetScores){
		
		int maxNrOfGenomes = getMaximuNumberOfGenomes(targetScores);
		double [][] sRNAscores = new double[maxNrOfGenomes][targetScores.size()];
		for(int j = 0; j < targetScores.size();j++){
			for(int k = 0; k < maxNrOfGenomes; k++){
				double[][] interactionScores = targetScores.get(j);
				if(interactionScores.length > k){
					sRNAscores[k][j] = interactionScores[k][0];
				}
			}
		}
		return sRNAscores;
	}
	
	public static int getMaximuNumberOfGenomes(ArrayList<double[][]> targetScores){
		int maximumNumberOfGenomes = 0;
		for(int j = 0; j < targetScores.size();j++){
				double[][] interactionScores = targetScores.get(j);
				if(interactionScores.length > maximumNumberOfGenomes){
					maximumNumberOfGenomes = interactionScores.length;
				}
			}
		return maximumNumberOfGenomes;
	}
		
/*
	public static void setExtremeValueDistributionValues(InteractionList interactions){
		if(interactions == null)
			return ;

		ArrayList<String> sRNANames = new ArrayList<String>();
		ArrayList<ArrayList<String>> mRNANames = new ArrayList<ArrayList<String>>();
		ArrayList<ArrayList<double[][]>> scores = new ArrayList<ArrayList<double[][]>>();
		interactions.getScores(mRNANames,scores,sRNANames);
		
		
		for(int i = 0; i < sRNANames.size(); i++){

			String sRNAName = sRNANames.get(i);
			if(!interactions.getsRNAs().EVDisSet(sRNAName)){
				double [][] sRNAscores = Statistics.getScores(scores.get(i));
				int maxNrOfGenomes = Statistics.getMaximuNumberOfGenomes(scores.get(i));
				for(int k = 0; k< maxNrOfGenomes; k++){
					double[] EVD1 = Statistics.getExtremeValueDistribution1(sRNAscores[k], sRNAscores[k].length, sRNAName, k);
					double[] EVD2 = Statistics.getExtremeValueDistribution2(sRNAscores[k], sRNAscores[k].length, sRNAName, k);
					if(EVD1 != null && EVD2 != null){
						if(EVD1[2] > EVD2[2]){
							interactions.getsRNAs().setXi(sRNAName, EVD1[0], k);
							interactions.getsRNAs().setTheta(sRNAName, EVD1[1], k);
						}
						else{
							interactions.getsRNAs().setXi(sRNAName, EVD2[0], k);
							interactions.getsRNAs().setTheta(sRNAName, EVD2[1], k);
						}
					}
					else
					{
						interactions.getsRNAs().setXi(sRNAName, 13-k, k);
						interactions.getsRNAs().setTheta(sRNAName, 1, k);
					}
					
				}
			}
		}
	}
	
	public static void setNormalDistributionValues(InteractionList interactions){
		if(interactions == null)
			return ;

		ArrayList<String> sRNANames = new ArrayList<String>();
		ArrayList<ArrayList<String>> mRNANames = new ArrayList<ArrayList<String>>();
		ArrayList<ArrayList<double[][]>> scores = new ArrayList<ArrayList<double[][]>>();
		interactions.getScores(mRNANames,scores,sRNANames);

		for(int i = 0; i < sRNANames.size(); i++){
			String sRNAName = sRNANames.get(i);
			double [][] sRNAscores = Statistics.getScores(scores.get(i));
			int maxNrOfGenomes = Statistics.getMaximuNumberOfGenomes(scores.get(i));
			for(int k = 0; k< maxNrOfGenomes; k++){
				double[] Norm = Statistics.getNormalDistribution(sRNAscores[k], sRNAscores[k].length, sRNAName, k);
				interactions.getsRNAs().setMean(sRNAName, Norm[0], k);
				interactions.getsRNAs().setSd(sRNAName, Norm[1], k);
			}
		}
	}
	*/
	
	public static double[] getExtremeValueDistribution1(double[] scores,int nrOfTargets,String sRNAName, int nrOfGenomes){
		
		
		double[] HitScores = Functions.sortScores(scores);
		int position = Functions.findLessThan(HitScores,0);
		HitScores = Functions.getSubarray(HitScores,position);
		HitScores = Functions.changePre(HitScores);
		//System.out.println(HitScores.length);

		if(HitScores == null )
			return null;

		Regression WD = new Regression(HitScores,0.2,0);
		//Regression.fitOneOrSeveralDistributions(HitScores);
		WD.setTitle(sRNAName+" "+nrOfGenomes);
		WD.gumbelMax();
		
		double[] gumbelValues = WD.getCoeff();
		double statValuesSD = WD.getSampleR();;
//		System.out.println("GumbellStats");
//		System.out.println("Xi    : "+gumbelValues[0]);
//		System.out.println("Theta : "+gumbelValues[1]);
//		System.out.println("Linear corr coef: " +statValuesSD);
		
		double[] StatValues = new double[3];
		StatValues[0] = gumbelValues[0];
		StatValues[1] = gumbelValues[1];
		StatValues[2] = statValuesSD;
		return StatValues;
		
		


	}
	
	
	public static double[] getExtremeValueDistribution2(double[] scores,int nrOfTargets,String sRNAName, int nrOfGenomes){
		scores = Functions.changePre(scores);
		double[] HitScores = Functions.sortScores(scores);
		

		double RValue = 0;
		double slope = 0;
		double intercept = 0;
		double threshhold = 1;
		int position = Functions.findLessThan(HitScores,threshhold);
		HitScores = Functions.getSubarray(HitScores,0,  position);
		double[] dens = new double[HitScores.length];
		
		double[] LLD = new double[HitScores.length];
		for(int i = 0; i < HitScores.length; i++){
			double density = (double)(i) / (double)(HitScores.length);
			LLD[HitScores.length-i-1] = java.lang.Math.log(-(java.lang.Math.log(density)));
			dens[HitScores.length-i-1] = density;
		}
		
//		HitScores = Functions.getSubarray(HitScores,HitScores.length*1/100, position);
//		LLD = Functions.getSubarray(LLD,LLD.length*1/100,position);
		
		
		
	 while(RValue < 0.90 && HitScores.length > 200){
			SimpleLinearRegressionModel slr = new SimpleLinearRegressionModel(HitScores,LLD);
			slope = slr.getSlope();
			intercept = slr.getIntercept();
			RValue = slr.getRSquared();
//			System.out.println("Model is y = " + intercept + 
//				" + " + slope 
//				+ " * x with R-Squared = " + RValue+" with "+HitScores.length  );
			HitScores = Functions.getSubarray(HitScores,HitScores.length*10/100,HitScores.length*99/100);
			LLD = Functions.getSubarray(LLD,LLD.length*10/100,LLD.length*99/100);
		}
	 	double Theta = -(1/ slope);
		double Xi = Theta *intercept;
//		System.out.println("Theta  = " + Theta + 
//						   "  Xi = " + Xi+
//						   "  with R-Squared2 = " + RValue);
		
		
		
		if(RValue < 0.90)
			return null;
		//System.out.println(" Mean : "+Functions.getMean(HitScores));
		//System.out.println(" Std : "+Functions.getStd(HitScores));

		double[] StatValues = new double[3];
		StatValues[0] = Xi;
		StatValues[1] = Theta;
		StatValues[2] = RValue;
		return StatValues;
	}	

	
	
	public static double[] getNormalDistribution(double[] scores,int nrOfTargets,String sRNAName, int nrOfGenomes){
		double[] HitScores = Functions.sortScores(scores);
		int position = Functions.findLessThan(HitScores,0);
		HitScores = Functions.getSubarray(HitScores,position);

		double mean = Functions.getMean(HitScores);
		double std = Functions.getStd(HitScores);
		System.out.println(" Mean : "+mean);
		System.out.println(" Std : "+std);

		double[] StatValues = new double[2];
		StatValues[0] = mean;
		StatValues[1] = std;
		return StatValues;
	}	

	
}
