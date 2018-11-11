package integrator;

import org.apache.commons.math3.ode.FirstOrderIntegrator;
import org.apache.commons.math3.ode.nonstiff.AdamsMoultonIntegrator;

public class AdamsMoultonFactory implements IntegratorFactory {

    double hmin;
    double hmax;
    
    
    /**
     * @param hmin the min step size for the method
     * @param hmax the max step size for the method
     */
    public AdamsMoultonFactory(double hmin, double hmax) {
        super();
        this.hmin = hmin;
        this.hmax = hmax;
    }

    @Override
    public FirstOrderIntegrator newIntegrator() {
        return new AdamsMoultonIntegrator(2, hmin, hmax, 1.0e-5, 1.0e-5);
    }

}
