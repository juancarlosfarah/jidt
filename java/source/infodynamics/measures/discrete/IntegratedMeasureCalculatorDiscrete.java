package infodynamics.measures.discrete;

import infodynamics.utils.Input;
import infodynamics.utils.MathsUtils;
import infodynamics.utils.MatrixUtils;

import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;

/**
 * <p>
 * Base class for our integrated information measure calculators on discrete
 * (int[]) data, providing common functionality for user-level measure classes.
 * </p>
 *
 * <p>
 * Usage of the child classes extending this class is intended to follow this
 * paradigm:
 * </p>
 * <ol>
 * 		<li>Construct the calculator;</li>
 * </ol>
 *
 * @author Juan Carlos Farah (<a href="farah.juancarlos at gmail.com">email</a>,
 * <a href="http://juancarlosfarah.com/">www</a>)
 */
public abstract class IntegratedMeasureCalculatorDiscrete {

    /**
     * Number of available quantised states for each variable
     * (ie binary is base-2).
     */
    protected int base;
    /**
     * Time step.
     */
    protected int tau;
    /**
     * Data.
     */
    protected int[][] data;
    /**
     * Set of all partitions (in our case bipartitions
     * for efficiency's sake) that the system can have.
     */
    protected Set<int[]> partitions;
    /**
     * Stores the indexes of the variables making
     * up the minimum information partition (MIP).
     */
    protected int[] minimumInformationPartition;
    /**
     * Stores the size of the MIP.
     */
    protected int minimumInformationPartitionSize;
    /**
     * Initialise the MIP score to positive infinity so that you return
     * the minimum score as scores for partitions start coming in.
     */
    protected double minimumInformationPartitionScore = Double.POSITIVE_INFINITY;
    /**
     * Stores the value of basic information measure (mutual information
     * or conditional entropy) for the system so that it can be reused.
     */
    protected double systemInformation;
    /**
     * Calculator to base the integrated measure on.
     */
    protected EffectiveMeasureCalculatorDiscrete baseCalculator;

    /**
     * Constructor.
     * @param base
     * @param tau
     */
    protected IntegratedMeasureCalculatorDiscrete(int base, int tau) {
        this.base = base;
        this.tau = tau;
        partitions = new HashSet<int[]>();
    }

    public void addObservations(int[][] data) {
        this.data = data;
    }

    public void computePossiblePartitions() {
        try {
            partitions.clear();
            for (int i = 1; i <= Math.floor(data.length / 2); i++) {
                int[][] sets = MathsUtils.generateAllSets(data.length, i);
                partitions.addAll(Arrays.asList(sets));
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public double compute() {

        double integratedInformation = 0.0;
        minimumInformationPartitionScore = Double.POSITIVE_INFINITY;
        baseCalculator.addObservations(data);
        systemInformation = baseCalculator.computeForSystem();

        for (int[] partition : partitions) {

            double k = computeNormalizationFactor(partition);
            double ei = baseCalculator.computeForPartition(partition);

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

    public double getSystemInformation() {
        return systemInformation;
    }

    public int getMinimumInformationPartitionSize() {
        return minimumInformationPartitionSize;
    }

    public int[] getMinimumInformationPartition() {
        return minimumInformationPartition;
    }

    /**
     * Returns observations added to the calculator.
     * @return data
     */
    public int[][] getData() {
        return data;
    }

}
