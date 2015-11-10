package infodynamics.measures.discrete;

/**
 * Created by juancarlosfarah on 22/05/15.
 *
 */
public class IntegratedInformationCalculatorDiscrete
       extends IntegratedMeasureCalculatorDiscrete {

    /**
     * Constructor.
     * @param base
     * @param tau
     */
    public IntegratedInformationCalculatorDiscrete(int base, int tau) {
        super(base, tau);
        baseCalculator = new EffectiveInformationCalculatorDiscrete(base, tau);
        shuffledCalculator = new EffectiveInformationCalculatorDiscrete(base,
                                                                        tau);
    }

}
