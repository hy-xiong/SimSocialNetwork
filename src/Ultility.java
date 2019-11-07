package Simulation;

import java.io.*;
import java.util.*;

public class Ultility {

	public static String[][] ReadParcelData(String inFile) {
		ArrayList<String[]> content = new ArrayList<String[]>();
		String[][] output;
		// Read LandParcel file
		try {
			BufferedReader br = new BufferedReader(new FileReader(new File(
					inFile)));
			String line = br.readLine();
			while ((line = br.readLine()) != null) {
				if (line.startsWith("\"")) {
					String[] temp1 = line.split("\",");
					String[] temp2 = temp1[1].split(",\"");
					String[] temp = new String[] {
							temp1[0].replace("\"", "").replace(",", ""),
							temp2[0], temp2[1].replace(",", ""),
							temp1[2].replace("\"", "").replace(",", "") };
					content.add(temp);
				} else {
					String[] temp1 = line.split(",", 3);
					String[] temp2 = temp1[2].split("\",\"");
					String[] temp = new String[] { temp1[0], temp1[1],
							temp2[0].replace("\"", "").replace(",", ""),
							temp2[1].replace("\"", "").replace(",", "") };
					content.add(temp);
				}
			}
			br.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		// Convert to 2d String Array to output
		output = new String[content.size()][];
		for (int i = 0; i < output.length; i++) {
			output[i] = content.get(i);
		}
		return output;
	}

	public static String[][] ReadSurveyData(String inFile) {
		String[][] output;
		ArrayList<String[]> content = new ArrayList<String[]>();
		try {
			BufferedReader br = new BufferedReader(new FileReader(new File(
					inFile)));
			String line = br.readLine();
			while ((line = br.readLine()) != null) {
				String[] temp = line.split("\t");
				content.add(temp);
			}
			br.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		output = new String[content.size()][];
		for (int i = 0; i < output.length; i++) {
			output[i] = content.get(i);
		}
		return output;
	}

	public static List<Agent> cloneAgents(List<Agent> agents) {
		List<Agent> listCopy = new ArrayList<Agent>();
		// Clone agents
		for (int i = 0; i < agents.size(); i++) {
			listCopy.add(agents.get(i).cloneExceptNeighbors());
		}
		// Clone social network
		for (int i = 0; i < agents.size(); i++) {
			for (int k = 0; k < agents.get(i).getNeighbors().size(); k++) {
				listCopy.get(i).addNeighbors(
						listCopy.get(agents.indexOf(agents.get(i)
								.getNeighbors().get(k))));
			}
		}
		return listCopy;
	}
}
