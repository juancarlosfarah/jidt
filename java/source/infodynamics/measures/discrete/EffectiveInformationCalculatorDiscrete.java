package infodynamics.measures.discrete;

import infodynamics.utils.MatrixUtils;
import infodynamics.utils.Input;

/**
 * Created by juancarlosfarah on 22/05/15.
 *
 */
public class EffectiveInformationCalculatorDiscrete
       extends EffectiveMeasureCalculatorDiscrete {

    public EffectiveInformationCalculatorDiscrete(int base, int tau) {
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

            MutualInformationCalculatorDiscrete micd;

            double sum = 0;

            // Calculate MI for the first partition.
            int[][] part1 = MatrixUtils.selectRows(data, p1);
            Input partition1 = new Input(part1, base);
            int p1Base = partition1.getReducedBase();
            int[][] p1Paired = partition1.pair(tau);
            micd = new MutualInformationCalculatorDiscrete(p1Base, 0);
            micd.initialise();
            micd.addObservations(p1Paired[0], p1Paired[1]);
            double p1Ei = micd.computeAverageLocalOfObservations();
            sum += p1Ei;

            // Calculate MI for the second partition.
            int[][] part2 = MatrixUtils.selectAllRowsExcept(data, p1);
            Input partition2 = new Input(part2, base);
            int p2Base = partition2.getReducedBase();
            int[][] p2Paired = partition2.pair(tau);
            micd = new MutualInformationCalculatorDiscrete(p2Base, 0);
            micd.initialise();
            micd.addObservations(p2Paired[0], p2Paired[1]);
            double p2Ei = micd.computeAverageLocalOfObservations();
            sum += p2Ei;

            // Subtract sum of MI of partitions from the MI of system.
            rvalue = systemInformation - sum;

        } catch (Exception e) {
            e.printStackTrace();
        }

        return rvalue;
    }
}
