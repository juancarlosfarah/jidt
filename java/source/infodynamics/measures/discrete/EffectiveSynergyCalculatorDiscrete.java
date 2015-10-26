package infodynamics.measures.discrete;

import infodynamics.utils.MatrixUtils;
import infodynamics.utils.Input;
import java.util.Arrays;

/**
 * Created by juancarlosfarah on 22/05/15.
 *
 */
public class EffectiveSynergyCalculatorDiscrete
       extends EffectiveMeasureCalculatorDiscrete {

    public EffectiveSynergyCalculatorDiscrete(int base, int tau) {
        super(base, tau);
    }

    // TODO: Refactor this out to EffectiveMeasureCalculatorDiscrete.
    public double computeForSystem() {

        try {
            MutualInformationCalculatorDiscrete micd;

            // Calculate MI for whole system.
            Input sys = new Input(data, base);
            int sysBase = sys.getReducedBase();
            int[][] sysPaired = sys.pair(tau);
            micd = new MutualInformationCalculatorDiscrete(sysBase, 0);
            micd.initialise();
            micd.addObservations(sysPaired[0], sysPaired[1]);
            systemInformation = micd.computeAverageLocalOfObservations();

        } catch (Exception e) {
            e.printStackTrace();
        }

        return systemInformation;
    }

    // TODO: Refactor this out to EffectiveMeasureCalculatorDiscrete.
    public double computeForPartition(int[] p1) {
        double rvalue = 0.0;

        try {

            double sum = 0;

            // data is in data[var][t] indexing. Transpose it to put it in
            // [t][var] indexing
            int[][] dataT = MatrixUtils.transpose(data);
            int observations = dataT.length;
            int numVars = dataT[0].length;

            int[] reducedSystem = MatrixUtils.computeCombinedValues(dataT, 2);
            int[] future = Arrays.copyOfRange(reducedSystem, tau, observations);
            int rBase = (int)(Math.pow(base, numVars));

            int[][] sysP1 = MatrixUtils.selectColumns(dataT, p1);
            int[] reducedP1 = MatrixUtils.computeCombinedValues(sysP1, 2);
            int[] past1 = Arrays.copyOfRange(reducedP1, 0, observations - tau);

            int[][] sysP2 = MatrixUtils.selectColumns(dataT, allExcept(p1, numVars));
            int[] reducedP2 = MatrixUtils.computeCombinedValues(sysP2, 2);
            int[] past2 = Arrays.copyOfRange(reducedP2, 0, observations - tau);

            // Declare mutual info calculator
            MutualInformationCalculatorDiscrete micd;
            micd = new MutualInformationCalculatorDiscrete(rBase);

            // Calculate MI for the first partition.
            micd.initialise();
            micd.addObservations(future, past1);
            sum += micd.computeAverageLocalOfObservations();

            // Calculate MI for the second partition.
            micd.initialise();
            micd.addObservations(future, past2);
            sum += micd.computeAverageLocalOfObservations();

            // Calculate redundancy between both partitions
            RedundantInformationCalculatorDiscrete riCalc = new
              RedundantInformationCalculatorDiscrete(rBase, 2);
            riCalc.initialise();
            int[][] bothSources = new int[observations - tau][2];
            MatrixUtils.copyIntoColumn(bothSources, 0, past1);
            MatrixUtils.copyIntoColumn(bothSources, 1, past2);
            riCalc.addObservations(future, bothSources);
            sum -= riCalc.compute();

            // Subtract sum of MI of partitions from the MI of system.
            rvalue = systemInformation - sum;

        } catch (Exception e) {
            e.printStackTrace();
        }

        return rvalue;
    }

    private static boolean[] allExcept(int[] idx, int N) {
        boolean[] v = new boolean[N];
        Arrays.fill(v, true);
        for (int i = 0; i < idx.length; i++) {
            v[idx[i]] = false;
        }
        return v;
    }

}
