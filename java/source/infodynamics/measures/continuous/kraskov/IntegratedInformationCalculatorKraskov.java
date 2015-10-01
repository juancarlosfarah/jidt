package infodynamics.measures.continuous.kraskov;
import infodynamics.utils.MathsUtils;
import infodynamics.utils.MatrixUtils;
import infodynamics.measures.continuous.kozachenko.EntropyCalculatorMultiVariateKozachenko;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;
/**
* Created by pmediano on 26/09/15.
*
*/
public class IntegratedInformationCalculatorKraskov {
  private int tau;
  private double[][] data;
  private int totalObservations;
  private int dimensions;
  public Set<int[]> partitions;
  public int[][] minimumInformationPartition;
  private double minimumInformationPartitionValue;
  private double mutualInformation;

  public IntegratedInformationCalculatorKraskov(int tau) {
    this.tau = tau;
    partitions = new HashSet<int[]>();
    minimumInformationPartition = new int[2][];
    minimumInformationPartitionValue = Double.POSITIVE_INFINITY;
  }

  public void addObservations(double[][] data) {
    this.data = MatrixUtils.transpose(data);
    this.totalObservations = data.length;
    this.dimensions = data[0].length;

    if (dimensions > totalObservations) {
      System.out.printf("The number of dimensions is smaller than the number"+
          " of observations. You're probably doing something wrong.");
    }
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
    EffectiveInformationCalculatorKraskov eick;
    eick = new EffectiveInformationCalculatorKraskov(tau);
    eick.addObservations(data);
    mutualInformation = eick.computeMutualInformationForSystem();
    for (int[] partition : partitions) {
      double k = computeNormalizationFactor(partition);
      double ei = eick.computeForBipartition(partition);
      // If k = 0, it means that one of the partitions has an entropy
      // of 0, which means that it doesn't tell us anything about the
      // rest of the system. Return 0 otherwise return normalised EI.
      double mipScore = (k == 0) ? 0 : ei / k;
      if (mipScore < minimumInformationPartitionValue) {
        minimumInformationPartition[0] = partition;
        int[] partition2 = new int[data.length - partition.length];
        int index = 0;
        for (int i = 0; i < data.length; i++) {
          if (!MatrixUtils.contains(partition, i)) {
            partition2[index] = i;
            index++;
          }
        }
        minimumInformationPartition[1] = partition2;
        minimumInformationPartitionValue = mipScore;
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
    EntropyCalculatorMultiVariateKozachenko eck = new EntropyCalculatorMultiVariateKozachenko();
    int dimensionsPart1 = part1.length;
    eck.initialise(dimensionsPart1);
    eck.setObservations(part1);
    double entropy1 = eck.computeAverageLocalOfObservations();
    int dimensionsPart2 = part2.length;
    eck.initialise(dimensionsPart2);
    eck.setObservations(part2);
    double entropy2 = eck.computeAverageLocalOfObservations();
    return Math.min(entropy1, entropy2);
  }

  public double getMutualInformation() {
    return mutualInformation;
  }

}
