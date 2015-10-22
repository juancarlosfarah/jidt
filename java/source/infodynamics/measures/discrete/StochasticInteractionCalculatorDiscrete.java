package infodynamics.measures.discrete;

import infodynamics.utils.MatrixUtils;
import infodynamics.utils.Input;

/**
 * Created by juancarlosfarah on 22/05/15.
 *
 */
public class StochasticInteractionCalculatorDiscrete
       extends EffectiveMeasureCalculatorDiscrete {

    public StochasticInteractionCalculatorDiscrete(int base, int tau) {
        super(base, tau);
    }

    public void addObservations(int[][] states) {
        data = states;
    }

    public double computeForSystem() {

        try {
            ConditionalEntropyCalculatorDiscrete cecd;

            // Calculate Conditional Entropy for whole system.
            Input sys = new Input(data, base);
            int sysBase = sys.getReducedBase();
            int[][] sysPaired = sys.pair(tau);
            cecd = new ConditionalEntropyCalculatorDiscrete(sysBase);
            cecd.addObservations(sysPaired[0], sysPaired[1]);
            systemInformation = cecd.compute();

        } catch (Exception e) {
            e.printStackTrace();
        }

        return systemInformation;
    }

    public double computeForPartition(int[] p1) {
        double rvalue = 0.0;

        try {

            ConditionalEntropyCalculatorDiscrete cecd;

            double sum = 0;

            // Calculate Conditional Entropy for the first partition.
            int[][] part1 = MatrixUtils.selectRows(data, p1);
            Input partition1 = new Input(part1, base);
            int p1Base = partition1.getReducedBase();
            int[][] p1Paired = partition1.pair(tau);
            cecd = new ConditionalEntropyCalculatorDiscrete(p1Base);
            cecd.addObservations(p1Paired[0], p1Paired[1]);
            double p1Ei = cecd.compute();
            sum += p1Ei;

            // Calculate Conditional Entropy for the second partition.
            int[][] part2 = MatrixUtils.selectAllRowsExcept(data, p1);
            Input partition2 = new Input(part2, base);
            int p2Base = partition2.getReducedBase();
            int[][] p2Paired = partition2.pair(tau);
            cecd = new ConditionalEntropyCalculatorDiscrete(p2Base);
            cecd.addObservations(p2Paired[0], p2Paired[1]);
            double p2Ei = cecd.compute();
            sum += p2Ei;

            // Subtract sum of MI of partitions from the MI of system.
            rvalue = sum - systemInformation;

        } catch (Exception e) {
            e.printStackTrace();
        }

        return rvalue;
    }


}
