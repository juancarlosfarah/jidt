package infodynamics.measures.discrete;

import infodynamics.utils.RandomGenerator;
import junit.framework.TestCase;
import java.util.Arrays;

/**
 * @author Juan Carlos Farah (<a href="farah.juancarlos at gmail.com">email</a>,
 * <a href="http://juancarlosfarah.com/">www</a>)
 */
public class IntegratedInformationTester extends TestCase {

    public void testAddObservations() {

        int[] var0 = { 0, 0, 0, 0, 0, 1, 1, 0, 0, 1 };
        int[] var1 = { 1, 1, 1, 0, 0, 1, 1, 0, 0, 1 };
        int[] var2 = { 1, 1, 1, 0, 0, 1, 1, 0, 0, 1 };
        int[] var3 = { 0, 0, 1, 0, 1, 1, 1, 0, 0, 1 };
        int[] var4 = { 1, 1, 1, 0, 1, 0, 1, 1, 0, 1 };
        int[] var5 = { 0, 0, 1, 0, 1, 1, 1, 0, 0, 1 };
        int[][] states0 = {var0, var1, var2, var3, var4, var5};
        int tau = 1;
        IntegratedInformationCalculatorDiscrete eicd;
        eicd = new IntegratedInformationCalculatorDiscrete(2, tau);
        eicd.addObservations(states0);
        int[][] states1 = eicd.getData();
        assertTrue(Arrays.equals(states0, states1));
    }

    public void testGenerativeModel() {
        // Example adapted from Wikipedia article on IIT.
        // https://en.wikipedia.org/wiki/Integrated_information_theory

        // Generate model.
        RandomGenerator rg = new RandomGenerator();
        int duration = 1000000;

        int[] var0a = rg.generateRandomInts(duration, 2);
        var0a[0] = 0;
        int[] var1a = new int[duration];
        var1a[0] = 0;
        System.arraycopy(var0a, 0, var1a, 1, duration - 1);
        int[][] states = {var0a, var1a};

        int tau = 1;
        IntegratedInformationCalculatorDiscrete iicd;
        iicd = new IntegratedInformationCalculatorDiscrete(2, tau);
        iicd.addObservations(states);
        iicd.computePossiblePartitions();

        // The integrated information for this model should be very close
        // to 1, but there is a margin of error that needs to be considered.
        double marginOfError = 0.001;
        double diff = Math.abs(iicd.compute() - 1);
        assertTrue(diff < marginOfError);
    }

}
