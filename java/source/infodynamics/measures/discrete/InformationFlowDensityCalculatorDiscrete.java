package infodynamics.measures.discrete;

import infodynamics.measures.discrete.ActiveInformationCalculatorDiscrete;
import infodynamics.measures.discrete.TransferEntropyCalculatorDiscrete;
import infodynamics.utils.MathsUtils;
import infodynamics.utils.MatrixUtils;
import infodynamics.utils.Input;

import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;

/**
 * Created by pmediano on 18/09/15.
 *
 * Names:
 *  - AIF: apparent information flow
 *  - CIF: complete information flow
 *  - AIS: active information storage
 *  - CER: community entropy rate
 *  - GER: global entropy rate
 *
 */
public class InformationFlowDensityCalculatorDiscrete {


    private int base;
    private int tau;
    private int[][] data;
    private int dimensions;
    private double AIF;
    private double CIF;
    private double AIS;
    private double CER;
    private double GER;
    private boolean isAIFcomputed;
    private boolean isCIFcomputed;
    private boolean isAIScomputed;
    private boolean isCERcomputed;
    private boolean isGERcomputed;


    public InformationFlowDensityCalculatorDiscrete(int base, int tau) {
      this.base = base;
      this.tau = tau;
    }

    public void addObservations(int[][] data) {
      this.data = MatrixUtils.transpose(data);
      this.dimensions = this.data[0].length;
      isAIFcomputed = false;
      isCIFcomputed = false;
      isAIScomputed = false;
      isCERcomputed = false;
      isGERcomputed = false;
    }

    public void computeBothFlows() {

      TransferEntropyCalculatorDiscrete tecd = new TransferEntropyCalculatorDiscrete(base, tau, tau);
      ConditionalTransferEntropyCalculatorDiscrete ctecd = new ConditionalTransferEntropyCalculatorDiscrete(base, tau, dimensions - 2);
      AIF = 0.0;
      CIF = 0.0;

      int[] others;

      for (int i = 0; i < dimensions; i++) {
        for (int j = 0; j < dimensions; j++) {
          if (i==j) continue;
          AIF += tecd.computeAverageLocal(data, i, j);

          others = allExcept(i, j, dimensions);
          CIF += ctecd.computeAverageLocal(data, i, j, others);
        }
      }

      AIF /= dimensions;
      CIF /= dimensions;

      isAIFcomputed = true;
      isCIFcomputed = true;
    }

    public double computeApparentFlow() {
      if (!isAIFcomputed) {
        computeBothFlows();
      }
      return AIF;
    }

    public double computeCompleteFlow() {
      if (!isCIFcomputed) {
        computeBothFlows();
      }
      return CIF;
    }



    public void computeStorageAndRates() {

      ActiveInformationCalculatorDiscrete aicd = new ActiveInformationCalculatorDiscrete(base, tau);
      AIS = 0.0;
      CER = 0.0;

      for (int i = 0; i < dimensions; i++) {
        aicd.addObservations(data, i);
        AIS += aicd.computeAverageLocalOfObservations();
        CER += aicd.computeAverageLocalEntropyRateOfObservations();
      }

      AIS /= dimensions;
      CER /= dimensions;

      isAIScomputed = true;
      isCERcomputed = true;

    }

    public double computeActiveStorage() {
      if (!isAIScomputed) {
        computeStorageAndRates();
      }
      return AIS;
    }

    public double computeCommunityEntropyRate() {
      if (!isCERcomputed) {
        computeStorageAndRates();
      }
      return CER;
    }

    public double computeGlobalEntropyRate() {
      if (!isGERcomputed) {

        Input input1 = new Input(MatrixUtils.transpose(data), base);
        int rBase = input1.getReducedBase();
        int[] rArray = input1.getReducedArray();
        EntropyRateCalculatorDiscrete ercd = new EntropyRateCalculatorDiscrete(rBase, tau);
        ercd.addObservations(rArray);
        GER = ercd.computeAverageLocalOfObservations();
        isGERcomputed = true;

      }

      return GER;
    }

    private static int[] allExcept(int i, int j, int N) {
      int[] res = new int[N - 2];
      int idx = 0;
      for (int k = 0; k < N; k++) {
        if ((k != i) && (k != j)) {
          res[idx] = k;
          idx++;
        }
      }
      return res;
    }
}

