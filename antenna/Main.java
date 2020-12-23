package antenna;

import java.util.concurrent.ThreadLocalRandom;

public class Main {

	public static void main(String[] args) {
		//parameters:number of antennae, steering angle
		AntennaArray a = new AntennaArray(4, 100.0);
		//parameters: AntennaArray instance, number of particles in initial population, termination condition(number of iterations)
		a.run(a, 10, 20);
	}
}
