/*
 *  Java Information Dynamics Toolkit (JIDT)
 *  Copyright (C) 2012, Joseph T. Lizier
 *  
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *  
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *  
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

package infodynamics.measures.discrete;

import infodynamics.utils.MatrixUtils;
import infodynamics.measures.discrete.MutualInformationCalculatorDiscrete;

/**
 * <p>Redundant information (I_min) calculator for univariate discrete (int[]) data.</p>
 * 
 * <p><b>References:</b><br/>
 * <ul>
 * 	<li>P. Williams, 'Information Dynamics: Its Theory and Application to
 * 	Embodied Cognitive Systems' (2011).</li>
 * </ul>
 * 
 * @author Pedro Mediano
 */
public class RedundantInformationCalculatorDiscrete {

  private int numSources = 0;
  private int observations = 0;
  private int base = 0;
  private int[] targetCount;

  private MutualInformationCalculatorDiscrete[] miCalc;
	
	
	/**
	 * Create a new redundant information calculator
	 * 
	 * @param base number of symbols for each variable.
	 *        E.g. binary variables are in base-2.
	 * @param numSources number of sources
	 */
	public RedundantInformationCalculatorDiscrete(int base, int numSources) {
		this.base = base;
		this.numSources = numSources;
    targetCount = new int[base];

    try {
      miCalc = new MutualInformationCalculatorDiscrete[numSources];
      for (int i = 0; i < numSources; i++) {
        miCalc[i] = new MutualInformationCalculatorDiscrete(base);
      }
    } catch (Throwable e) {
      e.printStackTrace();
    }

	}

	public void initialise() {
    observations = 0;
		MatrixUtils.fill(targetCount, 0);

    for (int i = 0; i < numSources; i++) {
      miCalc[i].initialise();
    }

	}
	
	/**
	 * Add observations for target and all sources simultaneously.
	 * 
   * @param target 1D array of states of the target variable.
   * @param sources 2D array with all source variables, indexed
   *        sources[time][var]. Rows are samples, columns are variables.
	 *
	 */
	public void addObservations(int[] target, int[][] sources) {

		observations += target.length;
		
    if (sources.length != target.length) {
      throw new RuntimeException("Target and all sources must have same number of samples");
    }
    if (sources[0].length != numSources) {
      throw new RuntimeException("Number of sources is different from the previously specified");
    }

    for (int t = 0; t < target.length; t++) {
      targetCount[target[t]]++;
    }

    for (int i = 0; i < numSources; i++) {
      miCalc[i].addObservations(target, MatrixUtils.selectColumn(sources, i));
    }

	}

	
  /**
   * Compute the redundant information I_min for the given
   * target and set of source variables.
   *
   */
	public double compute() {
		double ri = 0.0;
		double[] riTmp = new double[numSources];

    for (int s = 0; s < base; s++) {

      if (targetCount[s] == 0) {
        continue;
      }

      for (int i = 0; i < numSources; i++) {
        riTmp[i] = miCalc[i].computeSpecificSurprise(s);
      }

      ri +=  ( (double) targetCount[s] / (double) observations) * MatrixUtils.min(riTmp);

    }

    return ri;

	}


}

