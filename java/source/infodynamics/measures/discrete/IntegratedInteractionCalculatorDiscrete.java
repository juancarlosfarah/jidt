package infodynamics.measures.discrete;

/**
 * Created by juancarlosfarah on 22/05/15.
 *
 */
public class IntegratedInteractionCalculatorDiscrete
       extends IntegratedMeasureCalculatorDiscrete {

    /**
     * Constructor.
     * @param base
     * @param tau
     */
    public IntegratedInteractionCalculatorDiscrete(int base, int tau) {
        super(base, tau);
        baseCalculator = new StochasticInteractionCalculatorDiscrete(base, tau);
        shuffledCalculator = new StochasticInteractionCalculatorDiscrete(base,
                                                                         tau);
    }
}
