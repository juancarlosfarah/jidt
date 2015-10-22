package infodynamics.measures.discrete;

import infodynamics.measures.discrete.MutualInformationCalculatorDiscrete;
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
public class NeuralComplexityCalculatorDiscrete {

    private int base;
    private int[][] data;
    private int numVars;
    private double NC;


    public NeuralComplexityCalculatorDiscrete(int base) {
        this.base = base;
    }

    public void addObservations(int[][] data) {
        this.data = MatrixUtils.transpose(data);
        this.numVars = this.data[0].length;
    }

    private static boolean[] allExcept(int[] idx, int N) {
        boolean[] v = new boolean[N];
        Arrays.fill(v, true);
        for (int i = 0; i < idx.length; i++) {
            v[idx[i]] = false;
        }
        return v;
    }


    public double compute() {

        NC = 0;

        int maxSubsetSize = (int) Math.floor(numVars/2);
        Set<int[]> partitions = new HashSet<int[]>();

        // Loop for all possible subset sizes
        try {
          for (int k = 1; k < maxSubsetSize; k++) {

              // Set up calculator for MI with the right base for this bipartition size
              int pBase = (int) Math.pow(base, numVars - k);
              MutualInformationCalculatorDiscrete micd = new MutualInformationCalculatorDiscrete(pBase, 0);

              // Get all bipartitions of this size
              int[][] sets = MathsUtils.generateAllSets(numVars, k);
              partitions.clear();
              partitions.addAll(Arrays.asList(sets));

              double tmpMI = 0;

              // Loop over all bipartitions of this size
              for (int[] idx : partitions) {
                micd.initialise();

                int[][] part1 = MatrixUtils.selectColumns(data, idx);
                Input partition1 = new Input(MatrixUtils.transpose(part1), base);
                int[] var1 = partition1.getReducedArray();

                int[][] part2 = MatrixUtils.selectColumns(data, allExcept(idx, numVars));
                Input partition2 = new Input(MatrixUtils.transpose(part2), base);
                int[] var2 = partition2.getReducedArray();

                micd.addObservations(var1, var2);
                tmpMI += micd.computeAverageLocalOfObservations();

              }

              // To calculate NC, add average of MI for each subset size
              NC += tmpMI/partitions.size();

          }
        } catch (Exception e) {
          e.printStackTrace();
        }

        return NC;


    }

}

