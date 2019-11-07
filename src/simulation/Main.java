package simulation;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;

import org.jfree.ui.RefineryUtilities;

public class Main {
	public static void main(String[] args) {
		long startTIme = System.currentTimeMillis();
		// Read land parcel file to generate agents, land parcels are agents
		String DataDir = System.getProperty("user.dir") + "\\Data";
		String LandParcelFile = DataDir + "\\SamplePrivateOwnedParcelsPt2.txt";
		String SurveyFile = DataDir + "\\Survey.txt";
		String results1 = DataDir + "\\FinalAdopters.txt";
		String results2 = DataDir + "\\AdoptersDynamism.txt";
		String results3 = DataDir + "\\OpTh.txt";
		String results4 = DataDir + "\\opinions.txt";
		String results5 = DataDir + "\\result.txt";
		Sim sim_case = new Sim();
		ArrayList<ArrayList<String>> finalAdoptersOutput = new ArrayList<ArrayList<String>>();
		ArrayList<ArrayList<String>> AdoptersOutput = new ArrayList<ArrayList<String>>();
		/* ArrayList<double[][]> OpinionOutput = new ArrayList<double[][]>(); */
		/* Monte Carlo Simulation - Model Parameters */
		int MCSTimes = 50;
		double IncreaseStep = 0.015;
		double AvgAvailabledissenterPercentage = 0.37;
		int size = 0;
		int steps = 1000;
		double meanCC = 0.2;
		double stdDevCC = 0.1;
		double meanOT = 2.5;
		double stdDevOT = 0.5;
		DecimalFormat df2 = new DecimalFormat(".000");
		int[][] OpChangeCount = new int[MCSTimes][];
		double[][] finalOpDevs = new double[MCSTimes][];
		ArrayList<Double> AvgOpChangeCount = new ArrayList<Double>();
		ArrayList<Double> AvgFinalOpDevs = new ArrayList<Double>();
		/* Start Monte Carlo */
		for (int i = 0; i < MCSTimes; i++) {
			System.out.println("Start Monte Carlo Simulation #-" + i);
			List<Agent> agents = sim_case.InitializeAgents(LandParcelFile, SurveyFile, meanCC, stdDevCC, meanOT, stdDevOT);
			size = agents.size();
			int[] opChanges = new int[(int) (((int) (size * AvgAvailabledissenterPercentage) * 1.0)
					/ ((int) (size * IncreaseStep) * 1.0)) + 1];
			double[] finalOpDev = new double[(int) (((int) (size * AvgAvailabledissenterPercentage) * 1.0)
					/ ((int) (size * IncreaseStep) * 1.0)) + 1];

			List<Agent> dissenters = sim_case.GetInitialDissenters(agents);
			int[] rank = sim_case.sortOpTh(dissenters);
			try {
				BufferedWriter bw3 = new BufferedWriter(new FileWriter(new File(results3)));
				for (int o = 0; o < agents.size(); o++)
					bw3.write(agents.get(o).getOpinion() + "\t" + agents.get(o).getOpinionThreshold() + "\t"
							+ agents.get(o).getDegree() + "\n");
				bw3.close();
			} catch (IOException e) {
				e.printStackTrace();
			}
			// sim.outputDegreeInConsole(agents);
			ArrayList<String> tempFinalAdopter = new ArrayList<String>();
			tempFinalAdopter.add(Integer.toString(i));
			/* One simulation in Monte Carlo Simulation */
			/* for(double k = 0.0; k < 1.0; k += IncreasedAdoptersPer){ */
			for (int k = 0; k < (int) (size * AvgAvailabledissenterPercentage); k += (int) (size * IncreaseStep)) {
				System.out.println("Start Simulation " + df2.format(k) + " in Monte Carlo #-" + i);
				List<Agent> mirrorAgents = Ultility.cloneAgents(agents);
				List<Agent> mirrordissenters = sim_case.GetInitialDissenters(mirrorAgents);
				/*
				 * sim.IncreaseInitialAdoptersByThreshold(mirrordissenters,
				 * (int)(mirrordissenters.size() * k), rank);
				 */
				sim_case.IncreaseInitialAdoptersByThreshold(mirrordissenters, k, rank);
				double initialAdoptersPer = sim_case.GetAdoptersPercentage(mirrorAgents);
				// Set two dimensional list and one dimensional list to track
				// elements in recursive simulation
				List<List<Agent>> tempList = new ArrayList<List<Agent>>();
				ArrayList<String> tempAdoptersPer = new ArrayList<String>();
				tempAdoptersPer.add(i + "_" + df2.format(k));
				tempAdoptersPer.add(df2.format(initialAdoptersPer));

				// Get initial opinions
				double[][] opinions = new double[mirrorAgents.size()][steps + 1];
				for (int m = 0; m < mirrorAgents.size(); m++) {
					opinions[m][0] = mirrorAgents.get(m).getOpinion();
				}

				mirrorAgents = sim_case.Start(mirrorAgents, 1, steps, tempList, tempAdoptersPer, opinions);
				/*
				 * mirrorAgents = sim.simulation(mirrorAgents, 1, steps, tempList,
				 * tempAdoptersPer);
				 */

				/* Count opinions changes and final opinion standard deviation */
				opChanges[k / (int) (size * IncreaseStep)] = sim_case.CountOpinionChangeTimes(opinions);
				finalOpDev[k / (int) (size * IncreaseStep)] = sim_case.FinalOpinionDev(opinions);

				int finalAdopters = 0;
				for (Agent ag : mirrorAgents) {
					if (ag.getOpinion() > 3.0)
						finalAdopters++;
				}
				double finalAdoptersPer = finalAdopters * 1.0 / (agents.size() * 1.0);
				System.out.println("end: " + df2.format(finalAdoptersPer));
				tempFinalAdopter.add(df2.format(finalAdoptersPer));
				AdoptersOutput.add(tempAdoptersPer);
				/* OpinionOutput.add(opinions); */
				Visual vs = new Visual("Opinion Dynamism " + i + "_" + df2.format(k),
				opinions); vs.pack(); RefineryUtilities.centerFrameOnScreen(vs);
				vs.setVisible(true);
				 
			}
			OpChangeCount[i] = opChanges;
			finalOpDevs[i] = finalOpDev;
			finalAdoptersOutput.add(tempFinalAdopter);
		}
		for (int i = 0; i < OpChangeCount[0].length; i++) {
			double sumOpChange = 0.0;
			double sumFinalOpDevs = 0.0;
			for (int k = 0; k < OpChangeCount.length; k++) {
				sumOpChange += OpChangeCount[k][i] * 1.0;
				sumFinalOpDevs += finalOpDevs[k][i];
			}
			AvgOpChangeCount.add((sumOpChange / (OpChangeCount.length * 1.0)));
			AvgFinalOpDevs.add((sumFinalOpDevs / (finalOpDevs.length * 1.0)));
		}

		/* Write outputs */
		try {
			File f1 = new File(results1);
			File f2 = new File(results2);
			File f4 = new File(results4);
			File f5 = new File(results5);
			if (f1.exists())
				f1.delete();
			if (f2.exists())
				f2.delete();
			if (f4.exists())
				f4.delete();
			if (f5.exists())
				f5.delete();
			// Write final adopters
			BufferedWriter bw = new BufferedWriter(new FileWriter(f1));
			bw.write("MeanConvergenceCoefficient," + meanCC + ",StdDevConvergenceCoefficient," + stdDevCC
					+ ",MeanOpinionThreshold," + meanOT + ",StdDevOpinionThreshold," + stdDevOT + "\n");
			bw.write(",");
			// for(double k = 0; k < 1; k += IncreasedAdoptersPer){
			for (int k = 0; k < (int) (size * AvgAvailabledissenterPercentage); k += (int) (size * IncreaseStep)) {
				bw.write(df2.format(k) + ",");
			}
			bw.write("\n");
			for (ArrayList<String> arr : finalAdoptersOutput) {
				for (String s : arr) {
					bw.write(s + ", ");
				}
				bw.write("\n");
			}
			bw.close();
			// count succussful case for each condition
			BufferedWriter bw1 = new BufferedWriter(new FileWriter(f5));
			bw1.write("MeanConvergenceCoefficient," + meanCC + ",StdDevConvergenceCoefficient," + stdDevCC
					+ ",MeanOpinionThreshold," + meanOT + ",StdDevOpinionThreshold," + stdDevOT + "\n");
			int[] count = new int[(int) (((int) (size * AvgAvailabledissenterPercentage) * 1.0)
					/ ((int) (size * IncreaseStep) * 1.0)) + 1];
			for (int p = 0; p < count.length; p++)
				count[p] = 0;
			for (int k = 0; k < (int) (size * AvgAvailabledissenterPercentage); k += (int) (size * IncreaseStep)) {
				bw1.write(df2.format(k) + ",");
			}
			bw1.write("\n");
			for (ArrayList<String> arr : finalAdoptersOutput) {
				for (int m = 1; m < arr.size(); m++) {
					if (Double.parseDouble(arr.get(m)) > 0.5)
						count[m - 1]++;
				}
			}
			for (int c : count) {
				bw1.write(c + ", ");
			}
			bw1.write("\n");
			for (double d : AvgOpChangeCount) {
				bw1.write(df2.format(d) + ", ");
			}
			bw1.write("\n");
			for (double d : AvgFinalOpDevs) {
				bw1.write(df2.format(d) + ", ");
			}
			bw1.write("\n");
			bw1.close();
			// Writer adopters
			BufferedWriter bw2 = new BufferedWriter(new FileWriter(f2));
			bw2.write("MeanConvergenceCoefficient," + meanCC + ",StdDevConvergenceCoefficient," + stdDevCC
					+ ",MeanOpinionThreshold," + meanOT + ",StdDevOpinionThreshold," + stdDevOT + "\n");
			for (int k = 0; k <= steps; k++) {
				bw2.write(k + ",");
			}
			bw2.write("\n");
			for (ArrayList<String> arr : AdoptersOutput) {
				for (String s : arr) {
					bw2.write(s + ", ");
				}
				bw2.write("\n");
			}
			bw2.close();
			/*
			 * //Write opinoin distribution BufferedWriter bw4 = new BufferedWriter(new
			 * FileWriter(f4)); bw4.write("MeanConvergenceCoefficient," + meanCC +
			 * ",StdDevConvergenceCoefficient," + stdDevCC + ",MeanOpinionThreshold," +
			 * meanOT + ",StdDevOpinionThreshold," + stdDevOT + "\n"); for(int k = 0; k <
			 * MCSTimes; k++){ //for(int m = 0; m < 1; m += IncreasedAdoptersPer){ for(int m
			 * = 0; m < 69; m+=2){ bw4.write(k + "_" + m + ": " + "\n"); for(int l = 0; l <=
			 * steps; l++){ bw4.write(l + ", "); } bw4.write("\n"); //double[][] temp =
			 * OpinionOutput.get(k * (int)(1.0 / IncreasedAdoptersPer) + m); double[][] temp
			 * = OpinionOutput.get(k * 69/2+ m/2); for(int l = 0; l < temp.length; l++){
			 * for(int p = 0; p < temp[l].length; p++){ bw4.write(temp[l][p] + ", "); }
			 * bw4.write("\n"); } bw4.write("\n"); } bw4.write("\n\n\n\n"); } bw4.close();
			 */
		} catch (IOException e) {
			e.printStackTrace();
		}
		long endTime = System.currentTimeMillis();
		System.out.println("Running time: " + new DecimalFormat("0.").format(((endTime - startTIme) * 0.001)) + "s\n");
	}
}
