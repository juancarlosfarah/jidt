package infodynamics.measures.discrete;

import infodynamics.utils.Input;
import infodynamics.utils.MathsUtils;
import infodynamics.utils.MatrixUtils;

import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;

/**
 * Created by juancarlosfarah on 22/05/15.
 *
 */
public class IntegratedInformationEmpiricalCalculatorDiscrete {

    private int base;
    private int tau;
    private int[][] data;
    private Set<int[]> partitions;
    private int[] minimumInformationPartition;
    private int minimumInformationPartitionSize;
    private double minimumInformationPartitionScore;
    private double mutualInformation;


    public IntegratedInformationEmpiricalCalculatorDiscrete(int base, int tau) {
        this.base = base;
        this.tau = tau;
        partitions = new HashSet<int[]>();
        minimumInformationPartitionScore = Double.POSITIVE_INFINITY;
    }

    public void addObservations(int[][] data) {
        this.data = data;
    }

    public void computePossiblePartitions() {
        try {

            for (int i = 1; i < data.length; i++) {
                int[][] sets = MathsUtils.generateAllSets(data.length, i);
                partitions.addAll(Arrays.asList(sets));
            }

        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public double compute() {

        double integratedInformation = 0.0;
        EffectiveInformationCalculatorDiscrete eicd;
        eicd = new EffectiveInformationCalculatorDiscrete(base, tau);
        eicd.addObservations(data);
        mutualInformation = eicd.computeMutualInformationForSystem();

        for (int[] partition : partitions) {

            double k = computeNormalizationFactor(partition);
            double ei = eicd.computeForBipartition(partition);

            // If k = 0, it means that one of the partitions has an entropy
            // of 0, which means that it doesn't tell us anything about the
            // rest of the system. Return 0 otherwise return normalised EI.
            double mipScore = (k == 0) ? 0 : ei / k;

            if (mipScore < minimumInformationPartitionScore) {
                minimumInformationPartition = partition;
                minimumInformationPartitionSize = partition.length;
                minimumInformationPartitionScore = mipScore;
                integratedInformation = ei;
            }

        }
        return integratedInformation;
    }

    public double computeNormalizationFactor(int[] partition) {

        int[][] part1 = MatrixUtils.selectRows(this.data, partition);
        int[][] part2 = MatrixUtils.selectAllRowsExcept(this.data, partition);

        return computeNormalizationFactor(part1, part2);
    }


    public double computeNormalizationFactor(int[][] part1, int[][] part2) {

        // TODO: Fix EntropyCalculatorDiscrete for 2D inputs.

        // Prepare input for entropy calculator.
        Input input1 = new Input(part1, base);
        int[] p1 = input1.getReducedArray();
        int rBase1 = input1.getReducedBase();

        // Calculate entropy.
        EntropyCalculatorDiscrete ecd1 = new EntropyCalculatorDiscrete(rBase1);
        ecd1.initialise();
        ecd1.addObservations(p1);
        double entropy1 = ecd1.computeAverageLocalOfObservations();

        // Prepare input for entropy calculator.
        Input input2 = new Input(part2, base);
        int[] p2 = input2.getReducedArray();
        int rBase2 = input2.getReducedBase();

        // Calculate entropy.
        EntropyCalculatorDiscrete ecd2 = new EntropyCalculatorDiscrete(rBase2);
        ecd2.initialise();
        ecd2.addObservations(p2);
        double entropy2 = ecd2.computeAverageLocalOfObservations();

        return Math.min(entropy1, entropy2);

    }

    public double getMutualInformation() {
        return mutualInformation;
    }

    public int getMinimumInformationPartitionSize() {
        return minimumInformationPartitionSize;
    }

    public int[] getMinimumInformationPartition() {
        return minimumInformationPartition;
    }

}
