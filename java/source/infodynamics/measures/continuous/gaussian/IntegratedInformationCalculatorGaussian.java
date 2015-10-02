package infodynamics.measures.continuous.gaussian;
import infodynamics.utils.MathsUtils;
import infodynamics.utils.MatrixUtils;
import infodynamics.measures.continuous.gaussian.MutualInfoCalculatorMultiVariateGaussian;
import infodynamics.measures.continuous.gaussian.EntropyCalculatorMultiVariateGaussian;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;
/**
* Created by pmediano on 26/09/15.
*
*/
public class IntegratedInformationCalculatorGaussian {
  private int tau;
  private double[][] data;
  private int totalObservations;
  private int dimensions;
  private Set<int[]> partitions;
  private int[] minimumInformationPartition;
  private int minimumInformationPartitionSize;
  private double minimumInformationPartitionScore;
  private double mutualInformation;

  public IntegratedInformationCalculatorGaussian(int tau) {
    this.tau = tau;
    partitions = new HashSet<int[]>();
    minimumInformationPartitionScore = Double.POSITIVE_INFINITY;
  }

  public void addObservations(double[][] data) {
    this.data = MatrixUtils.transpose(data);
    this.dimensions = this.data[0].length;
    this.totalObservations = this.data.length;

    if (dimensions > totalObservations) {
      System.out.printf("The number of dimensions is smaller than the number"+
          " of observations. You're probably doing something wrong.");
    }
  }

  public void computePossiblePartitions() {
    try {
      for (int i = 1; i < Math.floor(dimensions/2); i++) {
        int[][] sets = MathsUtils.generateAllSets(dimensions, i);
        partitions.addAll(Arrays.asList(sets));
      }
    } catch (Exception e) {
      e.printStackTrace();
    }
  }

  public double compute() {
    double integratedInformation = 0.0;
    EffectiveInformationCalculatorGaussian eicg;
    eicg = new EffectiveInformationCalculatorGaussian(tau);
    eicg.addObservations(data);
    mutualInformation = eicg.computeMutualInformationForSystem();
    for (int[] partition : partitions) {
      double k = computeNormalizationFactor(partition);
      double ei = eicg.computeForBipartition(partition);

      // If k = 0, it means that one of the partitions has an entropy
      // of 0, which means that it doesn't tell us anything about the
      // rest of the system. Return 0 otherwise return normalised EI.
      double mipScore = (k == 0 || Double.isNaN(ei)) ? 0 : ei / k;

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

    double[][] part1 = MatrixUtils.selectColumns(this.data, partition);

    int[] rest = new int[dimensions - partition.length];
    int index = 0;
    for (int i = 0; i < dimensions; i++) {
      if (!MatrixUtils.contains(partition, i)) {
        rest[index] = i;
        index++;
      }
    }
    double[][] part2 = MatrixUtils.selectColumns(this.data, rest);
    return computeNormalizationFactor(part1, part2);
  }

  public double computeNormalizationFactor(double[][] part1, double[][] part2) {

    double entropy1 = 0, entropy2 = 0;

    try {
      EntropyCalculatorMultiVariateGaussian ecg =
              new EntropyCalculatorMultiVariateGaussian();
      int dimensionsPart1 = part1[0].length;
      ecg.initialise(dimensionsPart1);
      ecg.setObservations(part1);
      entropy1 = ecg.computeAverageLocalOfObservations();

      int dimensionsPart2 = part2[0].length;
      ecg.initialise(dimensionsPart2);
      ecg.setObservations(part2);
      entropy2 = ecg.computeAverageLocalOfObservations();

    } catch(Exception e) {
      e.printStackTrace();
    }

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
