package antenna;
import java.util.*;
import java.util.concurrent.ThreadLocalRandom;

/** Antenna array design problem */
public class AntennaArray {
    /** Minimum spacing permitted between antennae. */
    public static final double MIN_SPACING = 0.25;

    /**
     * Construct an antenna design problem.
     * @param n_ant Number of antennae in our array.
     * @param steering_ang Desired direction of the main beam in degrees.
     */
    public AntennaArray(int n_ant,double steering_ang) {
        n_antennae = n_ant;
        steering_angle = steering_ang;
    }
    
	public int getN_antennae() {
		return n_antennae;
	}
	
  /**
   * Rectangular bounds on the search space.
   * @return Vector b such that b[i][0] is the minimum permissible value of the
   * ith solution component and b[i][1] is the maximum.
   */
    public double[][] bounds() {
        double[][] bnds = new double[getN_antennae()][2];
        double[] dim_bnd = {0.0,((double)getN_antennae())/2.0};
        for(int i = 0;i<getN_antennae();++i)
            bnds[i] = dim_bnd;
        return bnds;
    }
    /**
     * Check whether an antenna design lies within the problem's feasible
     * region.
     * A design is a vector of n_antennae anntena placements.
     * A placement is a distance from the left hand side of the antenna array.
     * A valid placement is one in which
     *   1) all antennae are separated by at least MIN_SPACING
     *   2) the aperture size (the maximum element of the array) is exactly
     *      n_antennae/2.
     */
    public boolean is_valid(double[] design) {
        if(design.length != getN_antennae()) return false;
        double[] des = new double[design.length];
        System.arraycopy(design,0,des,0,design.length);
        Arrays.sort(des);

        //Aperture size is exactly n_antennae/2
        if(Math.abs(des[des.length - 1] - ((double)getN_antennae()) / 2.0)>1e-10)
            return false;
        //All antennae lie within the problem bounds
        for(int i = 0;i<des.length-1;++i)
            if(des[i] < bounds()[i][0] || des[i] > bounds()[i][1] )
                return false;
        //All antennae are separated by at least MIN_SPACING
        for(int i = 0;i<des.length-1;++i)
            if(des[i+1] - des[i] < MIN_SPACING || des[i+1] - des[i] < 0)
                return false;
        return true;
    }
    /**
     * Evaluate an antenna design returning peak SSL.
     * Designs which violate problem constraints will be penalised with extremely
     * high costs.
     * @param design A valid antenna array design.
     */
    public double evaluate(double[] design) {
        if(design.length != getN_antennae())
            throw new RuntimeException(
                    "AntennaArray::evaluate called on design of the wrong size. Expected: " + getN_antennae() +
                    ". Actual: " +
                    design.length
            );
        if(!is_valid(design)) return Double.MAX_VALUE;

        class PowerPeak {
            public double elevation;
            public double power;

            public PowerPeak(double e,double p){
                elevation = e;
                power = p;
            }
        }

        //Find all the peaks in power
        List<PowerPeak> peaks = new ArrayList<PowerPeak>();
        PowerPeak prev = new PowerPeak(0.0,Double.MIN_VALUE);
        PowerPeak current = new PowerPeak(0.0,array_factor(design,0.0));
        for(double elevation = 0.01; elevation <= 180.0; elevation += 0.01){
            PowerPeak next = new PowerPeak(elevation,array_factor(design,elevation));
            if(current.power >= prev.power && current.power >= next.power)
                peaks.add(current);
            prev = current;
            current = next;
        }
        peaks.add(new PowerPeak(180.0,array_factor(design,180.0)));

        Collections.sort(peaks,(PowerPeak l,PowerPeak r) -> l.power > r.power ? -1 : 1);

        //No side-lobes case
        if(peaks.size()<2) return Double.MIN_VALUE;
        //Filter out main lobe and then return highest lobe level
        final double distance_from_steering = Math.abs(peaks.get(0).elevation - steering_angle);
        for(int i=1;i<peaks.size();++i)
            if(Math.abs(peaks.get(i).elevation - steering_angle) < distance_from_steering)
                return peaks.get(0).power;
        return peaks.get(1).power;
    }

    private int n_antennae;
    private double steering_angle;

    private double array_factor(double[] design,double elevation) {
        double steering = 2.0*Math.PI*steering_angle/360.0;
        elevation = 2.0*Math.PI*elevation/360.0;
        double sum = 0.0;
        for(double x : design){
            sum += Math.cos(2 * Math.PI * x * (Math.cos(elevation) - Math.cos(steering)));
        }
        return 20.0*Math.log(Math.abs(sum));
    }
    
    public double[] randomPosition(AntennaArray a) {
    	double[] position = new double[a.getN_antennae()];
    	int count = 0;
    	while(count == 0) {
	    	for(int i = position.length-1;i >= 0;i--) {
	    		if(i == position.length-1) {
	    			position[i] = (double)a.getN_antennae()/2;
	    		}
	    		else {
	    			if(position[i+1] > 0.25) {
	    				position[i] = ThreadLocalRandom.current().nextDouble(0, position[i+1] - 0.25);
	    			}
	    		}
	    	}
	    	if(is_valid(position)) {
	    		count++;
	    	}
    	}
    	return position;
    }
    
    /**
     * Run AntennaArray design problem
     * @param a AntennaArray instance
     * @param numberOfParticles Number of particles in initial population
     * @param numberOfIterations Number of iteration to be done
     */
    public void run(AntennaArray a, int numberOfParticles, int terminationCondition) {
    	//Set up initial population with random particles in random position
		double[] gBestPos = new double[a.getN_antennae()];
		for(int i = 0;i < gBestPos.length;i++) {
			gBestPos[i] = Integer.MAX_VALUE;
		}
		int populationCount = 0;
		int particleCount = 1;
		Particle[] population = new Particle[numberOfParticles];
		while(populationCount < numberOfParticles) {
			population[populationCount] = new Particle(a, randomPosition(a));
			populationCount++;
		}
		//Print out all particles in initial population
		System.out.println("INITIAL POPULATION:");
		for(Particle p: population) {
			System.out.print("Particle " + particleCount + ": ");
			particleCount++;
			for(int i = 0;i<p.getPos().length;i++) {
				System.out.print(p.getPos()[i] + " ");
			}
			System.out.println(" = " + evaluate(p.getPos()));
		}
		System.out.println("-----------------------------------------------------------------------------");
		int numberOfIterations = 0;
		while(numberOfIterations < terminationCondition) {
			System.out.println("ITERATION " + numberOfIterations);
			System.out.println("GLOBAL BEST SSL: " + evaluate(gBestPos) + "\n");
			//iterate through all particles in population
			for(Particle currentParticle: population) {
				//Current position of particle
				double[] currentPos = currentParticle.getPos();
				//move particle to new position based on velocity
				double[] nextPos = new double[a.getN_antennae()];
				double[] nextVelocity = new double[a.getN_antennae()];
				for(int j = 0;j < currentPos.length;j++) {
					nextPos[j] = currentPos[j] + currentParticle.getVelocity()[j];
				}
				currentParticle.setPos(nextPos);
				//caluclate new velocity based on previous position so new position for particle can be calculated
				for(int j = 0;j < nextVelocity.length;j++) {
					nextVelocity[j] =  (double)((0.721 * currentParticle.getVelocity()[j]) + (1.1193 * ThreadLocalRandom.current().nextDouble(0, 1) * (currentParticle.getPBestPos()[j] - currentParticle.getPos()[j])) + (1.1193 * ThreadLocalRandom.current().nextDouble(0, 1) * (gBestPos[j] - currentParticle.getPos()[j])));
				} 
				currentParticle.setVelocity(nextVelocity);

				//Invisible wall: if solution isn't valid, don't evaluate position of particle
				if(is_valid(currentParticle.getPos())) {
					double sslOfCurrentPos = evaluate(currentParticle.getPos());
					if (sslOfCurrentPos < currentParticle.getPBestCost()) {
						currentParticle.setPBestPos(currentParticle.getPos());
						currentParticle.setPBestCost(sslOfCurrentPos);
						if (sslOfCurrentPos < evaluate(gBestPos)) {
							gBestPos = currentParticle.getPos();
						}
					}
				}
				System.out.println("Personal best cost: " + currentParticle.getPBestCost());
				System.out.print("Current Position: ");
				for(int cp = 0;cp < currentParticle.getPos().length;cp++) {
					System.out.print(currentParticle.getPos()[cp] + " ");
				}
				System.out.print("\nCurrent position velocity: ");
				for(int cp = 0;cp < currentParticle.getVelocity().length;cp++) {
					System.out.print(currentParticle.getVelocity()[cp] + " ");
				}
				System.out.print("\nCurrent position SSL: " + evaluate(currentParticle.getPos()));
				System.out.println("\n-----------------------------------------------------------------------------");
			}
			numberOfIterations++;
		}
		//Print global best position and SSL
		System.out.print("[ ");
		for(int i = 0;i < gBestPos.length;i++) {
			System.out.print(gBestPos[i] + " ");
		}
		System.out.print("]" + " = " + a.evaluate(gBestPos));
    }
}
