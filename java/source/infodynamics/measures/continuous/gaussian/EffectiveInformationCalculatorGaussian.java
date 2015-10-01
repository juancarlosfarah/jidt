package infodynamics.measures.continuous.gaussian;
import infodynamics.measures.continuous.gaussian.MutualInfoCalculatorMultiVariateGaussian;
import infodynamics.utils.MatrixUtils;

/**
* Created by pmediano on 26/09/15.
*
*/
public class EffectiveInformationCalculatorGaussian {
  private double[][] data;
  private int tau;
  private int totalObservations;
  private int dimensions;
  private double systemMutualInformation;

  public EffectiveInformationCalculatorGaussian(int tau) {
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
      MutualInfoCalculatorMultiVariateGaussian micg;
      // Calculate MI for whole system.
      double[][] paired0 = MatrixUtils.selectRows(data, 0, totalObservations - tau);
      double[][] paired1 = MatrixUtils.selectRows(data, tau, totalObservations - tau);
      micg = new MutualInfoCalculatorMultiVariateGaussian();
      micg.initialise(dimensions, dimensions);
      micg.setObservations(paired0, paired1);
      systemMutualInformation = micg.computeAverageLocalOfObservations();

    } catch (Exception e) {
      e.printStackTrace();
    }

    return systemMutualInformation;
  }

  public double computeForBipartition(int[] p1) {
    double rvalue = 0.0;
    try {
      MutualInfoCalculatorMultiVariateGaussian micg;
      double sum = 0;

      // Calculate MI for the first partition.
      double[][] part1 = MatrixUtils.selectColumns(data, p1);
      double[][] p1Paired0 = MatrixUtils.selectRows(part1, 0, totalObservations - tau);
      double[][] p1Paired1 = MatrixUtils.selectRows(part1, tau, totalObservations - tau);
      micg = new MutualInfoCalculatorMultiVariateGaussian();
      micg.initialise(p1.length, p1.length);
      micg.setObservations(p1Paired0, p1Paired1);
      double p1Ei = micg.computeAverageLocalOfObservations();
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
      micg = new MutualInfoCalculatorMultiVariateGaussian();
      micg.initialise(p2.length, p2.length);
      micg.setObservations(p2Paired0, p2Paired1);
      double p2Ei = micg.computeAverageLocalOfObservations();
      sum += p2Ei;

      // Subtract sum of MI of partitions from the MI of system.
      rvalue = systemMutualInformation - sum;
    } catch (Exception e) {
      e.printStackTrace();
    }
    return rvalue;
  }



}
