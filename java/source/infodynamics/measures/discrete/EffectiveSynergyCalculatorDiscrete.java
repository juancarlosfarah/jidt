package infodynamics.measures.discrete;

import infodynamics.utils.MatrixUtils;
import infodynamics.utils.Input;

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

        System.out.println("Method not implemented");

        /*
        try {

            MutualInformationCalculatorDiscrete micd;

            double sum = 0;

            // Calculate MI for the first partition.
            micd.initialise();
            micd.addObservations();
            sum += micd.computeAverageLocalOfObservations();

            // Calculate MI for the second partition.
            micd.initialise();
            micd.addObservations();
            sum += micd.computeAverageLocalOfObservations();

            // Calculate redundancy between both partitions
            RedundantInformationCalculatorDiscrete riCalc = new
              RedundantInformationCalculatorDiscrete(someBase);
            ri.initialise();
            ri.addObservations(future, [past1, past2]);
            sum -= ri.compute();

            // Subtract sum of MI of partitions from the MI of system.
            rvalue = systemInformation - sum;

        } catch (Exception e) {
            e.printStackTrace();
        }
        */

        return rvalue;
    }
}
