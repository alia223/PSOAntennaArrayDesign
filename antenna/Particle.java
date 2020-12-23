package antenna;

public class Particle {
	//current position of particle
	private double[] pos;
	//personal best position of particle
	private double[] pBestPos;
	//personal best SSL of particle
	private double pBestCost;
	//current velocity of particle
	private double[] velocity;
	
	/**
	 * Constructs AntennaAray design problem
	 * @param a AntennaArray instance
	 * @param position Position of particle
	 */
	public Particle(AntennaArray a, double[] position) {
		this.pos = position;
		this.pBestPos = position;
		this.pBestCost = a.evaluate(this.pBestPos);
		this.velocity = new double[a.getN_antennae()];
		double[] randomPosition = a.randomPosition(a);
		//initial velocity set to half distance between initial random position and a random second position
		for(int i = 0;i < a.getN_antennae();i++) {
			this.velocity[i] = (double)((randomPosition[i] - position[i])/2);
		}
	}
	
	/**
	 * Gets position of particle
	 * @return pos Current position of particle
	 */
	public double[] getPos() {
		return this.pos;
	}
	
	/**
	 * Sets position of particle
	 * @param pos New position to set position of particle to 
	 */
	public void setPos(double[] pos) {
		this.pos = pos;
	}
	
	/*
	 * Gets personal best position of particle
	 * @return pBestPos Personal best position of particle
	 */
	public double[] getPBestPos() {
		return this.pBestPos;
	}
	
	/**
	 * Sets personal best position of particle
	 * @param pBestPos New personal best position to set personal best position of particle to
	 */
	public void setPBestPos(double[] pBestPos) {
		this.pBestPos = pBestPos;
	}
	
	/**
	 * Gets velocity of particle
	 * @return velocity Velocity of Particle
	 */
	public double[] getVelocity() {
		return this.velocity;
	}
	
	/**
	 * Sets velocity of particle
	 * @param velocity New velocity to set velocity of particle to
	 */
	public void setVelocity(double[] velocity) {
		this.velocity = velocity;
	}
	
	/*
	 * Gets personal best cost of particle
	 * @return pBestPos Personal best cost of particle
	 */
	public double getPBestCost() {
		return this.pBestCost;
	}
	
	/**
	 * Sets personal best cost of particle
	 * @param pBestCost New personal best cost to set personal best cost of particle to
	 */
	public void setPBestCost(double pBestCost) {
		this.pBestCost = pBestCost;
	}
	
}
