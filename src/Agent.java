package Simulation;

import java.util.ArrayList;
import java.util.List;

/*This is a class for agents*/
public class Agent {
	private int _ParcelID;
	private int _AreaCategory;
	private double _ParcelX;
	private double _ParcelY;
	private double _Opinion;
	private List<Agent> _NeighborParcels;
	private int _Degree;
	private double _ConvergenceCoefficient;
	private double _OpinionThreshold;

	public Agent(int id, int areaCategory, double x, double y, double opinion,
			double converageCoefficient, double opinionThreshold) {
		_ParcelID = id;
		_AreaCategory = areaCategory;
		_ParcelX = x;
		_ParcelY = y;
		_Opinion = opinion;
		_NeighborParcels = new ArrayList<Agent>();
		_Degree = 0;
		_ConvergenceCoefficient = converageCoefficient;
		_OpinionThreshold = opinionThreshold;
	}

	public Agent cloneExceptNeighbors() {
		return new Agent(this._ParcelID, this._AreaCategory, this._ParcelX,
				this._ParcelY, this._Opinion, this._ConvergenceCoefficient,
				this._OpinionThreshold);
	}

	public int getParcelID() {
		return this._ParcelID;
	}

	public int getAreaCategory() {
		return this._AreaCategory;
	}

	public double getX() {
		return this._ParcelX;
	}

	public double getY() {
		return this._ParcelY;
	}

	public double getOpinion() {
		return this._Opinion;
	}

	public void setOpinion(double Op) {
		this._Opinion = Op;
	}

	public void addNeighbors(Agent lp) {
		_NeighborParcels.add(lp);
		_Degree++;
	}

	public int getDegree() {
		return _Degree;
	}

	public List<Agent> getNeighbors() {
		return _NeighborParcels;
	}

	public double getCovergenceCoefficient() {
		return _ConvergenceCoefficient;
	}

	public double getOpinionThreshold() {
		return _OpinionThreshold;
	}

	public void communicate(Agent ag) {
		double currentOwnOpinion = this._Opinion;
		double currentTargetOpinion = ag.getOpinion();
		if (Math.abs(currentTargetOpinion - currentOwnOpinion) <= this._OpinionThreshold)
			this._Opinion = currentOwnOpinion + this._ConvergenceCoefficient
					* (currentTargetOpinion - currentOwnOpinion);
		if (Math.abs(currentTargetOpinion - currentOwnOpinion) <= ag
				.getOpinionThreshold())
			ag.setOpinion(currentTargetOpinion + ag.getCovergenceCoefficient()
					* (currentOwnOpinion - currentTargetOpinion));
	}
}
