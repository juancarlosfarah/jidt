package infodynamics.measures.discrete;

/**
 * Created by juancarlosfarah on 22/05/15.
 *
 */
public class IntegratedSynergyCalculatorDiscrete
       extends IntegratedMeasureCalculatorDiscrete {

    /**
     * Constructor.
     * @param base
     * @param tau
     */
    public IntegratedSynergyCalculatorDiscrete(int base, int tau) {
        super(base, tau);
        baseCalculator = new EffectiveSynergyCalculatorDiscrete(base, tau);
    }

}

