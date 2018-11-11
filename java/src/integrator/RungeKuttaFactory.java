package integrator;

import org.apache.commons.math3.ode.FirstOrderIntegrator;
import org.apache.commons.math3.ode.nonstiff.ClassicalRungeKuttaIntegrator;

public class RungeKuttaFactory implements IntegratorFactory {
    
    double h;
    
    public RungeKuttaFactory(double h) {
        this.h = h;
    }

    @Override
    public FirstOrderIntegrator newIntegrator() {
        return new ClassicalRungeKuttaIntegrator(h);
    }

}
