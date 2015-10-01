package infodynamics.measures.discrete;

import infodynamics.measures.discrete.MutualInformationCalculatorDiscrete;
import infodynamics.measures.discrete.MultiInformationCalculatorDiscrete;
import infodynamics.utils.MathsUtils;
import infodynamics.utils.MatrixUtils;
import infodynamics.utils.Input;

import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;

/**
 * Created by pmediano on 18/09/15.
 *
 */
public class InteractionComplexityCalculatorDiscrete {

    private int base;
    private int[][] data;
    private int numVars;
    private double mutualInf;
    private double multiInf;
    private double IC;


    public InteractionComplexityCalculatorDiscrete(int base) {
        this.base = base;
    }

    public void addObservations(int[][] data) {
        this.data = MatrixUtils.transpose(data);
        this.numVars = this.data[0].length;
    }


    public double compute() {

        // Necessary to run MultiInfoCalculator
        int[] groupOffsets = new int[numVars];
        Arrays.fill(groupOffsets, 1);

        // 1. Calculate MultiInformation
        MultiInformationCalculatorDiscrete multiInfCalc = new MultiInformationCalculatorDiscrete(base, numVars);
        multiInfCalc.initialise();
        multiInfCalc.addObservations(data, groupOffsets);
        multiInf = multiInfCalc.computeAverageLocalOfObservations();


        // 2. Calculate MI for all bipartitions of size (1,N-1) 
        double mutualInf = 0;
        try {
          MutualInformationCalculatorDiscrete mutualInfCalc = new MutualInformationCalculatorDiscrete((int) Math.pow(base, numVars-1), 0);
          boolean[] idx = new boolean[numVars];
          
          for (int i = 0; i < numVars; i++) {
              mutualInfCalc.initialise();

              // Now idx is an array with all 1's and a 0 in position i
              Arrays.fill(idx, true);
              idx[i] = false;

              // Separate the data and reduce the N-1 remaining variables to a single array
              int[][] part = MatrixUtils.selectColumns(data, idx); // All variables except i
              Input partition = new Input(part, base);
              int[] var1 = partition.getReducedArray();
              int[] var2 = MatrixUtils.selectColumn(data, i); // Only variable i

              mutualInfCalc.addObservations(var1, var2);
              mutualInf += mutualInfCalc.computeAverageLocalOfObservations();

          }
        } catch (Exception e) {
          e.printStackTrace();
        }

        // 3. Calculate InteractionComplexity
        IC = mutualInf/numVars - multiInf;

        return IC;

    }

    public double getMultiInformation() {
      return multiInf;
    }

    public double getAverageMutualInformation() {
        return (mutualInf/numVars);
    }

}
