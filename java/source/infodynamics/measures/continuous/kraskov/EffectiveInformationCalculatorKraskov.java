package infodynamics.measures.continuous.kraskov;
import infodynamics.measures.continuous.kraskov.MutualInfoCalculatorMultiVariateKraskov1;
import infodynamics.utils.MatrixUtils;

/**
* Created by pmediano on 26/09/15.
*
*/
public class EffectiveInformationCalculatorKraskov {
  private double[][] data;
  private int tau;
  private int totalObservations;
  private int dimensions;
  private double systemMutualInformation;

  public EffectiveInformationCalculatorKraskov(int tau) {
    this.tau = tau;
  }

  public void addObservations(double[][] states) {
    data = states;
    totalObservations = data.length;
    dimensions = data[0].length;
    if (dimensions > totalObservations) {
      System.out.printf("The number of dimensions is smaller than the number"+
          " of observations. You're probably doing something wrong.");
    }
  }

  public double computeMutualInformationForSystem() {
    try {
      MutualInfoCalculatorMultiVariateKraskov1 mick;
      // Calculate MI for whole system.
      double[][] paired0 = MatrixUtils.selectRows(data, 0, totalObservations - tau);
      double[][] paired1 = MatrixUtils.selectRows(data, tau, totalObservations - tau);
      mick = new MutualInfoCalculatorMultiVariateKraskov1();
      mick.initialise(dimensions, dimensions);
      mick.addObservations(paired0, paired1);
      systemMutualInformation = mick.computeAverageLocalOfObservations();

    } catch (Exception e) {
      e.printStackTrace();
    }

    return systemMutualInformation;
  }

  public double computeForBipartition(int[] p1) {
    double rvalue = 0.0;
    try {
      MutualInfoCalculatorMultiVariateKraskov1 mick;
      double sum = 0;

      // Calculate MI for the first partition.
      double[][] part1 = MatrixUtils.selectColumns(data, p1);
      double[][] p1Paired0 = MatrixUtils.selectRows(part1, 0, totalObservations - tau);
      double[][] p1Paired1 = MatrixUtils.selectRows(part1, tau, totalObservations - tau);
      mick = new MutualInfoCalculatorMultiVariateKraskov1();
      mick.initialise(p1.length, p1.length);
      mick.addObservations(p1Paired0, p1Paired1);
      double p1Ei = mick.computeAverageLocalOfObservations();
      sum += p1Ei;

      // Calculate MI for the second partition.
      int[] p2 = new int[dimensions - p1.length];
      int index = 0;
      for (int i = 0; i < dimensions; i++) {
        if (!MatrixUtils.contains(p1, i)) {
          p2[index] = i;
          index++;
        }
      }

      double[][] part2 = MatrixUtils.selectColumns(data, p2);
      double[][] p2Paired0 = MatrixUtils.selectRows(part2, 0, totalObservations - tau);
      double[][] p2Paired1 = MatrixUtils.selectRows(part2, tau, totalObservations - tau);
      mick = new MutualInfoCalculatorMultiVariateKraskov1();
      mick.initialise(p2.length, p2.length);
      mick.addObservations(p2Paired0, p2Paired1);
      double p2Ei = mick.computeAverageLocalOfObservations();
      sum += p2Ei;

      // Subtract sum of MI of partitions from the MI of system.
      rvalue = systemMutualInformation - sum;
    } catch (Exception e) {
      e.printStackTrace();
    }
    return rvalue;
  }



}
