package integrator;

import org.apache.commons.math3.ode.FirstOrderIntegrator;

public interface IntegratorFactory {
    FirstOrderIntegrator newIntegrator();
}
