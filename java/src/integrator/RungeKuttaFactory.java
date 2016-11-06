package integrator;

import org.apache.commons.math3.ode.FirstOrderIntegrator;
import org.apache.commons.math3.ode.nonstiff.ClassicalRungeKuttaIntegrator;

public class RungeKuttaFactory implements IntegratorFactory {
    
    double[][] a;
    double[] b;
    double[] c;
    double h;
    
    public RungeKuttaFactory(double h) {
        this.h = h;
        
        c = new double[3];
        c[0] = 0.5; 
        c[1] = 0.5; 
        c[2] = 1.0;
        
        b = new double[4];
        b[0] = 1.0/6.0;
        b[1] = 1.0/3.0;
        b[2] = 1.0/3.0;
        b[3] = 1.0/6.0;
        
        a = new double[3][3];
        a[0][0] = 0.5; 
        a[1][0] = 0.0; 
        a[1][1] = 0.5; 
        a[2][0] = 0.0; 
        a[2][1] = 0.0; 
        a[2][2] = 1.0; 
    }

    @Override
    public FirstOrderIntegrator newIntegrator() {
        return new ClassicalRungeKuttaIntegrator(h);
    }

}
