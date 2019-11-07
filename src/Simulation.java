package Simulation;

import java.text.DecimalFormat;
import java.util.*;
import java.io.*;

import org.jfree.ui.RefineryUtilities;

/*This is main class*/
public class Simulation {
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
		Simulation sim = new Simulation();
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
			List<Agent> agents = sim.InitializeAgents(LandParcelFile,
					SurveyFile, meanCC, stdDevCC, meanOT, stdDevOT);
			size = agents.size();
			int[] opChanges = new int[(int) (((int) (size * AvgAvailabledissenterPercentage) * 1.0) / ((int) (size * IncreaseStep) * 1.0)) + 1];
			double[] finalOpDev = new double[(int) (((int) (size * AvgAvailabledissenterPercentage) * 1.0) / ((int) (size * IncreaseStep) * 1.0)) + 1];

			List<Agent> dissenters = sim.GetInitialDissenters(agents);
			int[] rank = sim.sortOpTh(dissenters);
			try {
				BufferedWriter bw3 = new BufferedWriter(new FileWriter(
						new File(results3)));
				for (int o = 0; o < agents.size(); o++)
					bw3.write(agents.get(o).getOpinion() + "\t"
							+ agents.get(o).getOpinionThreshold() + "\t"
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
				System.out.println("Start Simulation " + df2.format(k)
						+ " in Monte Carlo #-" + i);
				List<Agent> mirrorAgents = Ultility.cloneAgents(agents);
				List<Agent> mirrordissenters = sim
						.GetInitialDissenters(mirrorAgents);
				/*
				 * sim.IncreaseInitialAdoptersByThreshold(mirrordissenters,
				 * (int)(mirrordissenters.size() * k), rank);
				 */
				sim.IncreaseInitialAdoptersByThreshold(mirrordissenters, k,
						rank);
				double initialAdoptersPer = sim
						.GetAdoptersPercentage(mirrorAgents);
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

				mirrorAgents = sim.simulation(mirrorAgents, 1, steps, tempList,
						tempAdoptersPer, opinions);
				/*
				 * mirrorAgents = sim.simulation(mirrorAgents, 1, steps,
				 * tempList, tempAdoptersPer);
				 */

				/* Count opinions changes and final opinion standard deviation */
				opChanges[k / (int) (size * IncreaseStep)] = sim
						.CountOpinionChangeTimes(opinions);
				finalOpDev[k / (int) (size * IncreaseStep)] = sim
						.FinalOpinionDev(opinions);

				int finalAdopters = 0;
				for (Agent ag : mirrorAgents) {
					if (ag.getOpinion() > 3.0)
						finalAdopters++;
				}
				double finalAdoptersPer = finalAdopters * 1.0
						/ (agents.size() * 1.0);
				System.out.println("end: " + df2.format(finalAdoptersPer));
				tempFinalAdopter.add(df2.format(finalAdoptersPer));
				AdoptersOutput.add(tempAdoptersPer);
				/* OpinionOutput.add(opinions); */
				/*
				 * Visual vs = new Visual("Opinion Dynamism " + i + "_" +
				 * df2.format(k), opinions); vs.pack();
				 * RefineryUtilities.centerFrameOnScreen(vs);
				 * vs.setVisible(true);
				 */
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
			bw.write("MeanConvergenceCoefficient," + meanCC
					+ ",StdDevConvergenceCoefficient," + stdDevCC
					+ ",MeanOpinionThreshold," + meanOT
					+ ",StdDevOpinionThreshold," + stdDevOT + "\n");
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
			bw1.write("MeanConvergenceCoefficient," + meanCC
					+ ",StdDevConvergenceCoefficient," + stdDevCC
					+ ",MeanOpinionThreshold," + meanOT
					+ ",StdDevOpinionThreshold," + stdDevOT + "\n");
			int[] count = new int[(int) (((int) (size * AvgAvailabledissenterPercentage) * 1.0) / ((int) (size * IncreaseStep) * 1.0)) + 1];
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
			bw2.write("MeanConvergenceCoefficient," + meanCC
					+ ",StdDevConvergenceCoefficient," + stdDevCC
					+ ",MeanOpinionThreshold," + meanOT
					+ ",StdDevOpinionThreshold," + stdDevOT + "\n");
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
			 * //Write opinoin distribution BufferedWriter bw4 = new
			 * BufferedWriter(new FileWriter(f4));
			 * bw4.write("MeanConvergenceCoefficient," + meanCC +
			 * ",StdDevConvergenceCoefficient," + stdDevCC +
			 * ",MeanOpinionThreshold," + meanOT + ",StdDevOpinionThreshold," +
			 * stdDevOT + "\n"); for(int k = 0; k < MCSTimes; k++){ //for(int m
			 * = 0; m < 1; m += IncreasedAdoptersPer){ for(int m = 0; m < 69;
			 * m+=2){ bw4.write(k + "_" + m + ": " + "\n"); for(int l = 0; l <=
			 * steps; l++){ bw4.write(l + ", "); } bw4.write("\n"); //double[][]
			 * temp = OpinionOutput.get(k * (int)(1.0 / IncreasedAdoptersPer) +
			 * m); double[][] temp = OpinionOutput.get(k * 69/2+ m/2); for(int l
			 * = 0; l < temp.length; l++){ for(int p = 0; p < temp[l].length;
			 * p++){ bw4.write(temp[l][p] + ", "); } bw4.write("\n"); }
			 * bw4.write("\n"); } bw4.write("\n\n\n\n"); } bw4.close();
			 */
		} catch (IOException e) {
			e.printStackTrace();
		}
		long endTime = System.currentTimeMillis();
		System.out.println("Running time: "
				+ new DecimalFormat("0.")
						.format(((endTime - startTIme) * 0.001)) + "s\n");
	}

	/**
	 * Simulate the model in recursive way until reach the required amount of
	 * simulation time
	 * 
	 * @param agents
	 *            a list of agents
	 * @param steps
	 *            number of steps for this simulation
	 */
	public List<Agent> simulation(List<Agent> agents, int startStep, int steps,
			List<List<Agent>> temp, ArrayList<String> tempAdoptersPer,
			double[][] opinions) {
		/*
		 * public List<Agent> simulation(List<Agent> agents, int startStep, int
		 * steps, List<List<Agent>> temp, ArrayList<String> tempAdoptersPer){
		 */
		/*
		 * Opinion exchange model - concept of bounded confidence; asymmetric
		 * updating
		 */
		if (steps > 0) {
			List<Agent> UpdatedList = new ArrayList<Agent>();
			/*
			 * In each setup, each agent is randomly picked and assumed to
			 * randomly communicate with one of its neighbor
			 */
			Random rd = new Random();
			while (agents.size() > 0) {
				if (agents.size() > 1) {
					int index = rd.nextInt(agents.size());
					Agent ag = agents.get(index);
					if (ag.getNeighbors().size() > 0) {
						Agent neighbor = ag.getNeighbors().get(
								rd.nextInt(ag.getNeighbors().size()));
						ag.communicate(neighbor);
						opinions[ag.getParcelID() - 1][startStep] = ag
								.getOpinion();
						opinions[neighbor.getParcelID() - 1][startStep] = neighbor
								.getOpinion();
					} else {
						opinions[ag.getParcelID() - 1][startStep] = ag
								.getOpinion();
					}
					agents.remove(index);
					UpdatedList.add(ag);
				} else {
					Agent ag = agents.get(0);
					if (ag.getNeighbors().size() > 0) {
						Agent neighbor = ag.getNeighbors().get(
								rd.nextInt(ag.getNeighbors().size()));
						ag.communicate(neighbor);
						opinions[ag.getParcelID() - 1][startStep] = ag
								.getOpinion();
						opinions[neighbor.getParcelID() - 1][startStep] = neighbor
								.getOpinion();
					} else {
						opinions[ag.getParcelID() - 1][startStep] = ag
								.getOpinion();
					}
					agents.remove(0);
					UpdatedList.add(ag);
				}
			}
			int adopters = 0;
			for (Agent ag : UpdatedList) {
				if (ag.getOpinion() > 3.0)
					adopters++;
			}
			steps--;
			startStep++;
			temp.add(UpdatedList);
			tempAdoptersPer.add(new DecimalFormat(".000")
					.format(GetAdoptersPercentage(UpdatedList)));
			simulation(UpdatedList, startStep, steps, temp, tempAdoptersPer,
					opinions);
			/* simulation(UpdatedList, startStep, steps, temp, tempAdoptersPer); */
		}
		return temp.get(temp.size() - 1);
	}

	/**
	 * Setup Initial condition
	 * 
	 * @param ParcelFile
	 *            Path to land parcel file
	 * @param SurveyFile
	 *            Path to survey file
	 * @param outFIle
	 *            output path of file for initial land area and opinion setup
	 * @param meanConvergenceCoefficient
	 *            mean of the Normal distribution of opinion convergence
	 *            coefficient
	 * @param stdDevConvergenceCoefficient
	 *            standard deviation of the Normal distribution of opinion
	 *            convergence coefficient
	 * @param meanOpinionThreshold
	 *            mean of the Normal distribution of opinion threshold
	 * @param stdDevOpinionThreshold
	 *            standard deviation of the Normal distribution of opinion
	 *            threshold
	 * @return a list of agents with initial condition set
	 */
	public List<Agent> InitializeAgents(String ParcelFile, String SurveyFile,
			double meanConvergenceCoefficient,
			double stdDevConvergenceCoefficient, double meanOpinionThreshold,
			double stdDevOpinionThreshold) {
		/*
		 * Read survey data to setup correlation between opinion and land area
		 * owned
		 */
		String[][] OpLand = Ultility.ReadSurveyData(SurveyFile);
		int[] OpDis = new int[OpLand.length]; // Opinion value: 1-Strongly
												// disagree; 2-Disagree;
												// 3-Neutral; 4-Agree;
												// 5-Strongly Agree
		int[] LandDis = new int[OpLand.length]; // Land area owned
		for (int i = 0; i < OpLand.length; i++) {
			OpDis[i] = Integer.parseInt(OpLand[i][0]);
			LandDis[i] = Integer.parseInt(OpLand[i][1]);
		}
		double[][] OpDisLand = OpDisLandCalc(OpDis, LandDis); // Opinion
																// probability
																// in each land
																// area category
		/* Read parcel data to setup agents */
		List<Agent> agents = new ArrayList<Agent>(); // Agents - LandParcels
		String[][] parcelsInfo = Ultility.ReadParcelData(ParcelFile);
		for (String[] parcelInfo : parcelsInfo) {
			int id = Integer.parseInt(parcelInfo[0]);
			int areaCategory = Integer.parseInt(parcelInfo[1]);
			double x = Double.parseDouble(parcelInfo[2]);
			double y = Double.parseDouble(parcelInfo[3]);
			int opinion = randomOpCalc(OpDisLand, areaCategory);
			double rdz = 9999.0;
			while (Math.abs(rdz) > meanConvergenceCoefficient
					/ stdDevConvergenceCoefficient)
				rdz = new Random().nextGaussian();
			double covergenceCoefficient = rdz * stdDevConvergenceCoefficient
					+ meanConvergenceCoefficient;
			rdz = 9999.0;
			while (Math.abs(rdz) > meanOpinionThreshold
					/ stdDevOpinionThreshold)
				rdz = new Random().nextGaussian();
			double opinionThreshold = rdz * stdDevOpinionThreshold
					+ meanOpinionThreshold;
			/*
			 * double covergenceCoefficient = 0.2; double opinionThreshold =
			 * 1.5;
			 */
			Agent lp = new Agent(id, areaCategory, x, y, opinion,
					covergenceCoefficient, opinionThreshold);
			agents.add(lp);
		}
		/* Set neighborhoods of agents */
		BuildSocialNetwork(agents);
		return agents;
	}

	public int CountOpinionChangeTimes(double[][] opinion) {
		int count = 0;
		for (int i = 0; i < opinion.length; i++) {
			for (int k = 0; k < opinion[i].length - 1; k++) {
				if ((opinion[i][k] > 3.0) && (opinion[i][k + 1] < 3.0))
					count++;
				else if ((opinion[i][k] < 3.0) && (opinion[i][k + 1] > 3.0))
					count++;
				else {
				}
			}
		}
		return count;
	}

	public double FinalOpinionDev(double[][] opinion) {
		double stddev = 0.0;
		double mean = 0.0;
		double sum = 0.0;
		double numerator = 0.0;
		double n = opinion.length;
		for (int i = 0; i < opinion.length; i++) {
			sum += opinion[i][opinion[i].length - 1];
		}
		mean = sum / n;
		for (int i = 0; i < opinion.length; i++) {
			numerator += Math.pow(opinion[i][opinion[i].length - 1] - mean, 2);
		}
		stddev = Math.sqrt(numerator / (n - 1));
		return stddev;
	}

	/**
	 * Setup the amount of persuaded initial adopters
	 * 
	 * @param agents
	 *            a list of policy dissenters
	 * @param amount
	 *            the number of persuaded initial adopters
	 */
	public void IncreaseInitialAdopters(List<Agent> disenters, int amount) {
		while (amount > 0) {
			int index = new Random().nextInt(disenters.size());
			disenters.get(index).setOpinion(4.0);
			disenters.remove(index);
			amount--;
		}
	}

	public void IncreaseInitialAdoptersByThreshold(List<Agent> dissenters,
			int amount, int[] rank) {
		for (int i = 0; i < amount; i++) {
			dissenters.get(rank[i]).setOpinion(4.0);
		}
	}

	public int[] sortOpTh(List<Agent> dissenters) {
		int[] rank = new int[dissenters.size()];
		for (int i = 0; i < rank.length; i++) {
			rank[i] = i;
		}
		for (int i = 0; i < rank.length - 1; i++) {
			for (int k = i + 1; k < rank.length; k++) {
				if (dissenters.get(rank[k]).getOpinionThreshold() < dissenters
						.get(rank[i]).getOpinionThreshold()) {
					int temp = rank[i];
					rank[i] = rank[k];
					rank[k] = temp;
				}
			}
		}
		return rank;
	}

	/* return total number of initial agents */
	public double GetAdoptersPercentage(List<Agent> agents) {
		int count = 0;
		for (Agent ag : agents) {
			if (ag.getOpinion() > 3.0) {
				count++;
			}
		}
		return (count * 1.0) / (agents.size() * 1.0);
	}

	/**
	 * Get initial policy dissenters
	 * 
	 * @param agents
	 *            a list of agents
	 * @return a list of policy dissenters
	 */
	public List<Agent> GetInitialDissenters(List<Agent> agents) {
		List<Agent> disenters = new ArrayList<Agent>();
		for (Agent ag : agents) {
			if (ag.getOpinion() < 3.0) {
				disenters.add(ag);
			}
		}
		return disenters;
	}

	/**
	 * Opinion distribution in each land area category
	 * 
	 * @param OpDist
	 *            Opinion distribution
	 * @param LandDist
	 *            Land area owned distribution
	 * @return Opinion distribution in each land area category: row - land area
	 *         owned; col - opinion
	 */
	private double[][] OpDisLandCalc(int[] OpDist, int[] LandDist) {
		double[][] OpDisLand = new double[5][5];
		// Calculate opinion distribution in each land area category
		for (int i = 0; i < OpDisLand.length; i++) {
			int count_i = 0;
			ArrayList<Integer> tempOp = new ArrayList<Integer>();
			for (int k = 0; k < LandDist.length; k++) {
				if (LandDist[k] == (i + 1)) {
					count_i++;
					tempOp.add(OpDist[k]);
				}
			}
			for (int o = 0; o < OpDisLand[i].length; o++) {
				int count_o = 0;
				for (int p : tempOp) {
					if (p == (o + 1))
						count_o++;
				}
				OpDisLand[i][o] = count_o * 1.0 / count_i;
			}
		}
		// Calculate cumulative opinion distribution in each land area category
		for (int i = 0; i < OpDisLand.length; i++) {
			for (int k = 1; k < OpDisLand[i].length; k++) {
				OpDisLand[i][k] = OpDisLand[i][k] + OpDisLand[i][k - 1];
				if (OpDisLand[i][k] > 0.999999) {
					OpDisLand[i][k] = 1.0;
				}
			}
			// for(int k = 0; k < OpDisLand[i].length; k++) System.out.println(i
			// + ": " + OpDisLand[i][k]);
		}
		return OpDisLand;
	}

	/**
	 * Randomly generate Opinion based on land area agent owned
	 * 
	 * @param OpDisLand
	 *            Opinion distribution in each land area category
	 * @param areaCategory
	 *            Land area agent owned, categorized into integer level
	 * @return opinion value of agent: 1-Strongly disagree; 2-Disagree;
	 *         3-Neutral; 4-Agree; 5-Strongly Agree
	 */
	private int randomOpCalc(double[][] OpDisLand, int areaCategory) {
		int opinion = 0;
		double rd = new Random().nextDouble();
		for (int k = 0; k < OpDisLand[areaCategory - 1].length; k++) {
			if (rd <= OpDisLand[areaCategory - 1][k]) {
				opinion = k + 1;
				break;
			}
		}
		return opinion;
	}

	/**
	 * Build social network by assign neighbors to each agent in the input list
	 * of agents
	 * 
	 * @param a
	 *            list of agents
	 */
	private void BuildSocialNetwork(List<Agent> agents) {
		double[] x = new double[agents.size()];
		double[] y = new double[agents.size()];
		for (int i = 0; i < agents.size(); i++) {
			x[i] = agents.get(i).getX();
			y[i] = agents.get(i).getY();
		}
		double[][] ProbConn = new double[x.length][x.length];
		for (int i = 0; i < ProbConn.length; i++) {
			for (int k = i + 1; k < ProbConn[i].length; k++) {
				// Initial prob = 1 / distance
				ProbConn[i][k] = 1 / Math.sqrt(Math.pow(x[i] - x[k], 2)
						+ Math.pow(y[i] - y[k], 2));
			}
			// Get max and min of each agent's ProbConn
			double max = 0.0;
			double min = 999999999999999999.0;
			for (int k = i + 1; k < ProbConn[i].length; k++) {
				if (max < ProbConn[i][k])
					max = ProbConn[i][k];
				if ((min > ProbConn[i][k]))
					min = ProbConn[i][k];
			}
			// Normalize prob by min-max scale times 95%
			for (int k = i + 1; k < ProbConn[i].length; k++) {
				ProbConn[i][k] = (ProbConn[i][k] - min) / (max - min) * 0.90;
			}
		}
		// build social network
		for (int i = 0; i < agents.size(); i++) {
			for (int k = i + 1; k < ProbConn[i].length; k++) {
				double rd = new Random().nextDouble();
				if (rd < ProbConn[i][k]) {
					agents.get(i).addNeighbors(agents.get(k));
					agents.get(k).addNeighbors(agents.get(i));
				}
			}
		}
	}

	/* Output degree setup */
	public void outputDegreeInConsole(List<Agent> agents) {
		int sum_degree = 0;
		double avg_degree = 0;
		double std_degree = 0;
		int min_degree = 999999;
		int max_degree = 0;
		for (int i = 0; i < agents.size(); i++) {
			int temp = agents.get(i).getDegree();
			sum_degree += temp;
			if (min_degree > temp)
				min_degree = temp;
			if (max_degree < temp)
				max_degree = temp;
		}
		avg_degree = sum_degree * 1.0 / (agents.size() * 1.0);
		double temp_sum = 0;
		for (int i = 0; i < agents.size(); i++) {
			temp_sum += Math.pow(agents.get(i).getDegree() - avg_degree, 2);
		}
		std_degree = Math.sqrt(temp_sum / ((agents.size() - 1) * 1.0));
		System.out.println("avg_degree: " + avg_degree + "\nstd_degree: "
				+ std_degree + "\nmin: " + min_degree + "\nmax: " + max_degree);
	}

	/* output agent Opinion & Land attributes */
	public void outputAgentsOpLandText(String outFile, List<Agent> agents) {
		try {
			File outF = new File(outFile);
			if (outF.exists())
				outF.delete();
			BufferedWriter bw = new BufferedWriter(new FileWriter(new File(
					outFile)));
			for (Agent a : agents) {
				bw.write(a.getOpinion() + "\t" + a.getAreaCategory() + "\n");
			}
			bw.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

}
