package integrator;

import org.apache.commons.math3.ode.FirstOrderIntegrator;
import org.apache.commons.math3.ode.nonstiff.EulerIntegrator;

public class EulerIntegratorFactory implements IntegratorFactory {

    double h;
    
    
    /**
     * @param h the step size for Euler's method
     */
    public EulerIntegratorFactory(double h) {
        super();
        this.h = h;
    }

    @Override
    public FirstOrderIntegrator newIntegrator() {
        return new EulerIntegrator(h);
    }

}
