package infodynamics.measures.discrete;

/**
 * Created by juancarlosfarah on 22/10/15.
 *
 */
public abstract class EffectiveMeasureCalculatorDiscrete {
    protected int[][] data;
    protected int base;
    protected int tau;
    protected double systemInformation;

    /**
     * Constructor.
     * @param base
     * @param tau
     */
    public EffectiveMeasureCalculatorDiscrete(int base, int tau) {
        this.base = base;
        this.tau = tau;
    }

    /**
     * Adds observations.
     * @param states
     */
    public void addObservations(int[][] states) {
        data = states;
    }

    // TODO: Implement here.
    public abstract double computeForSystem();

    // TODO: Implement here.
    public abstract double computeForPartition(int[] p1);
}
