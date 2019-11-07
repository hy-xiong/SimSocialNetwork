package Simulation;

import java.awt.Color;

import javax.swing.JPanel;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.axis.ValueAxis;
import org.jfree.chart.axis.Axis.*;
import org.jfree.chart.plot.FastScatterPlot;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.data.xy.XYDataset;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.ui.ApplicationFrame;
import org.jfree.ui.RectangleInsets;
import org.jfree.ui.RefineryUtilities;

public class Visual extends ApplicationFrame {

	public Visual(String title, double[][] xy) {
		super(title);

		/*
		 * float[][] data = new float[xy.length][xy[0].length]; for(int i = 0; i
		 * < data.length; i++){ for(int k = 0; k < data[i].length; k++)
		 * data[i][k] = (float) xy[i][k]; } ValueAxis domainAxis = new
		 * NumberAxis("Steps"); ValueAxis rangeAxis = new NumberAxis("Opinion");
		 * FastScatterPlot fsp = new FastScatterPlot(data, domainAxis,
		 * rangeAxis); JFreeChart chart = new JFreeChart("Fast Scatter Plot",
		 * fsp);
		 */
		XYDataset dataset = createDataset(xy);
		JFreeChart chart = createChart(dataset);
		ChartPanel chartpanel = new ChartPanel(chart);
		chartpanel.setPreferredSize(new java.awt.Dimension(1000, 500));
		setContentPane(chartpanel);
	}

	private XYDataset createDataset(double[][] xy) {
		XYSeriesCollection dataset = new XYSeriesCollection();
		for (int i = 0; i < xy.length; i++) {
			XYSeries serie = new XYSeries(Integer.toString(i), false);
			for (int k = 0; k < xy[i].length; k++) {
				serie.add((double) k, xy[i][k], false);
			}
			dataset.addSeries(serie);
		}
		return dataset;
	}

	private JFreeChart createChart(XYDataset dataset) {
		JFreeChart chart = ChartFactory.createScatterPlot("Opinion Dynamism",
				"Steps", "Opinion", dataset, PlotOrientation.VERTICAL, false,
				true, false);
		/*
		 * chart.setBackgroundPaint(Color.white); XYPlot plot = (XYPlot)
		 * chart.getPlot(); plot.setBackgroundPaint(Color.lightGray);
		 * plot.setAxisOffset(new RectangleInsets(5.0, 5.0, 5.0, 5.0));
		 * plot.setDomainGridlinePaint(Color.white);
		 * plot.setRangeGridlinePaint(Color.white); XYLineAndShapeRenderer
		 * renderer = (XYLineAndShapeRenderer) plot.getRenderer();
		 * renderer.setShapesVisible(true); renderer.setShapesFilled(true);
		 * NumberAxis rangeAxis = (NumberAxis) plot.getRangeAxis();
		 * rangeAxis.setStandardTickUnits(NumberAxis.createIntegerTickUnits());
		 */
		return chart;
	}
	/*
	 * public static JPanel createDemoPanel(double[][] xy) { JFreeChart chart =
	 * createChart(createDataset(xy)); return new ChartPanel(chart); }
	 */
}
