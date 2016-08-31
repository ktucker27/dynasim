package integrator;

import org.apache.commons.math3.ode.FirstOrderIntegrator;
import org.apache.commons.math3.ode.nonstiff.AdamsMoultonIntegrator;

public class AdamsMoultonFactory implements IntegratorFactory {

    double h;
    
    
    /**
     * @param h the min and max step size for the method
     */
    public AdamsMoultonFactory(double h) {
        super();
        this.h = h;
    }

    @Override
    public FirstOrderIntegrator newIntegrator() {
        return new AdamsMoultonIntegrator(2, 1.0e-18, h, 1.0e-3, 1.0e-2);
    }

}
