package poc.statistics_in_java;

import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.distribution.TDistribution;
import org.apache.commons.math3.stat.StatUtils;
import org.apache.commons.math3.stat.inference.TTest;
import org.apache.commons.math3.stat.inference.TestUtils;
import org.springframework.boot.SpringApplication;
import org.springframework.boot.autoconfigure.SpringBootApplication;

import java.io.IOException;
import java.net.InetAddress;

@SpringBootApplication
public class StatisticsInJavaApplication {

	/**
	 * Keep in mind that the Java stats library does not have a direct way to test if the observed values are
	 * greater or less than a given mean, which is easily possible in R
	 * @param args, command line args if needed
	 */
	public static void main(String[] args) {
		SpringApplication.run(StatisticsInJavaApplication.class, args);
		System.out.println("\n***Starting Statistics test***");

		try {
			String target = "127.0.0.1";
			int runs = 1;
			int measurements = 1000;

			// gather data
			double[] latencies = new double[runs*measurements];
			double[] riggedLatencies = new double[runs*measurements];
			double riggedMean = 15000000.0;
			double riggedStandardDeviation = 10000.0;
			double riggedModifier = 1.01;  // scales the sampled values by a little bit for sanity testing

			NormalDistribution normalDistribution = new NormalDistribution(riggedMean, riggedStandardDeviation);

			// outer loop which collects the measurements per run (resp. per experiment)
			for (int run = 0; run < runs; run++) {
				Thread.sleep(25);

				// inner loop which collects the measurements for statistical data analysis due to the stochastic nature
				for (int measurement = 0; measurement < measurements; measurement++) {
					// Sleep for a while before sending the next packet (optional)
					Thread.sleep(25);

					long startTime = System.nanoTime();
					boolean isReachable = InetAddress.getByName(target).isReachable(1000); // Timeout in milliseconds
					long endTime = System.nanoTime();

					if (isReachable) {
						latencies[measurement + measurements*run] = Math.log(endTime - startTime);
						riggedLatencies[measurement + measurements*run] = Math.log(normalDistribution.sample()*riggedModifier);
					}

				}
			}

			// perform t-Test
			double alpha = 0.01;
			double populationMean = Math.log(riggedMean);
			double degreesOfFreedom = runs*measurements - 1; // df = N - #groups, #groups = 1 (one-sample t-test)
			TTest tTest = new TTest();
			TDistribution tDistribution = new TDistribution(degreesOfFreedom);

			var isSignificant = tTest.tTest(populationMean, latencies, alpha*2);
			// Calculate the critical value for the upper tail
			double criticalValue = tDistribution.inverseCumulativeProbability(1 - (alpha/2));  // alpha/2 due to being one-sided

			System.out.println("\n***Perform test with observed values***");
			// false => sample mean is drawn from the same population as populationMean, which is good
			// thus, the sample mean is not significantly different from populationMean under a given alpha
			// (and confidence interval with confidence level = 1-alpha)
			System.out.println("Result from t-Test: " + isSignificant);
			System.out.println("Mean of observed values: " + StatUtils.mean(latencies));
			System.out.println("SD of observed values: " + Math.sqrt(StatUtils.variance(latencies)));
			System.out.println("is significant (H_0 is rejected): " + TestUtils.tTest(populationMean, latencies, alpha*2));  // *2 because it is one-sided
			System.out.println("p-value: " + TestUtils.tTest(populationMean, latencies));
			System.out.println("t-statistic: " + TestUtils.t(populationMean, latencies));
			// calculate the critical values to better check if the observed values are above or below
			System.out.println("Critical Value: " + criticalValue);
			System.out.println("Comparing the means: computedMean/populationMean = " + StatUtils.mean(latencies)/populationMean);
			System.out.println("t-Statistic > Critical Value: " + (TestUtils.t(populationMean, latencies) > criticalValue));


			System.out.println("\n***Perform test with RIGGED values***");
			System.out.println("Mean of observed values: " + StatUtils.mean(riggedLatencies));
			System.out.println("is significant (H_0 is rejected): " + TestUtils.tTest(populationMean, riggedLatencies, alpha*2));  // *2 because it is one-sided
			System.out.println("p-value: " + TestUtils.tTest(populationMean, riggedLatencies));  // *2 because it is one-sided
			System.out.println("t-statistic: " + TestUtils.t(populationMean, riggedLatencies));
			// calculate the critical values to better check if the observed values are above or below
			System.out.println("Critical Value: " + criticalValue);
			System.out.println("Comparing the means: computedMean/populationMean = " + StatUtils.mean(riggedLatencies)/populationMean);
			System.out.println("t-Statistic > Critical Value: " + (TestUtils.t(populationMean, riggedLatencies) > criticalValue));



		} catch (InterruptedException | IOException e) {
            throw new RuntimeException(e);
        }

    }

}
