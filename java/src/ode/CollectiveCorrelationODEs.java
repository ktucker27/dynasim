package ode;

import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.exception.MaxCountExceededException;

import utils.DynaComplex;

public class CollectiveCorrelationODEs implements DynaComplexODEs {

    private int n;
    private double gamma;
    private double w;
    private DynaComplex alpha;
    
    double sz, szz, pm;
    
    DynaComplex t1, t2, t3;
    
    public CollectiveCorrelationODEs(SystemParams params, double sz, double szz, double pm) {
        n = params.getN();
        gamma = params.getGamma();
        w = params.getW();
        alpha = new DynaComplex(params.getAlpha());
        
        this.sz = sz;
        this.szz = szz;
        this.pm = pm;
        
        t1 = new DynaComplex();
        t2 = new DynaComplex();
        t3 = new DynaComplex();
    }
    
    @Override
    public void computeDerivatives(double t, DynaComplex[] z, DynaComplex[] zDot)
            throws MaxCountExceededException, DimensionMismatchException {
        // <sigma_a^+(t + tau) sigma_b^-(t)>
        t1.set(z[2]).multiply(n - 2).add(z[3]).multiply(t2.set(alpha).conjugate().multiply(0.5*gamma));
        zDot[0].set(z[0]).multiply(-0.5*(gamma + w)).add(t1);
        
        // <sigma_a^+(t + tau) sigma_a^-(t)>
        t1.set(z[4]).multiply(t2.set(alpha).conjugate().multiply(0.5*gamma*(n-1)));
        zDot[1].set(z[1]).multiply(-0.5*(gamma + w)).add(t1);
        
        // <sigma_a^z(t + tau) sigma_b^+(t + tau) sigma_c^-(t)>
        t1.set(z[0]).multiply(szz - 2.0*sz*sz).add(t3.set(z[2]).multiply(2*sz)).multiply(n - 3);
        t2.set(z[1]).multiply(szz - 2.0*sz*sz).add(t3.set(z[3]).multiply(2*sz));
        zDot[2].set(t1).add(t2).multiply(t3.set(alpha).conjugate()).multiply(0.5*gamma);
        zDot[2].subtract(t1.set(z[2]).multiply(1.5*(gamma + w)));
        zDot[2].subtract(t1.set(z[0]).multiply(gamma - w));
        zDot[2].subtract(t1.set(z[0]).multiply(alpha).multiply(0.5*gamma));
        zDot[2].subtract(t1.set(z[2]).multiply(gamma*alpha.getReal()));
        
        zDot[2].subtract(t1.set(alpha).multiply(gamma*(n-2)).multiply(t2.set(z[0]).multiply(2*pm)));
        zDot[2].subtract(t1.set(alpha).conjugate().multiply(gamma).multiply(t2.set(z[1]).add(z[0]).multiply(pm)));
        zDot[2].subtract(t1.set(alpha).conjugate().multiply(gamma*(n-3)).multiply(t2.set(z[0]).multiply(2*pm)));
        
        // <sigma_a^z(t + tau) sigma_b^+(t + tau) sigma_b^-(t)>
        t1.set(z[0]).multiply(szz - 2.0*sz*sz).add(t2.set(z[2]).multiply(sz)).add(t2.set(z[4]).multiply(sz));
        zDot[3].set(t1).multiply(t1.set(alpha).conjugate()).multiply(0.5*gamma*(n-2));
        zDot[3].subtract(t1.set(z[3]).multiply(1.5*(gamma + w)));
        zDot[3].subtract(t1.set(z[1]).multiply(gamma - w));
        zDot[3].subtract(t1.set(z[0]).multiply(alpha).multiply(0.5*gamma));
        zDot[3].subtract(t1.set(z[4]).multiply(gamma*alpha.getReal()));
        
        zDot[3].subtract(t1.set(alpha).multiply(gamma*(n-2)).multiply(t2.set(z[0]).add(z[1]).multiply(pm)));
        zDot[3].subtract(t1.set(alpha).conjugate().multiply(gamma*(n-2)).multiply(t2.set(z[0]).add(z[1]).multiply(pm)));

        // <sigma_a^z(t + tau) sigma_b^+(t + tau) sigma_a^-(t)>
        t1.set(z[0]).multiply(szz - 2.0*sz*sz).add(t2.set(z[2]).multiply(sz)).add(t2.set(z[4]).multiply(sz));
        zDot[4].set(t1).multiply(t1.set(alpha).conjugate()).multiply(0.5*gamma*(n-2));
        zDot[4].subtract(t1.set(z[4]).multiply(1.5*(gamma + w)));
        zDot[4].subtract(t1.set(z[0]).multiply(gamma - w));
        zDot[4].subtract(t1.set(z[1]).multiply(alpha).multiply(0.5*gamma));
        zDot[4].subtract(t1.set(z[3]).multiply(gamma*alpha.getReal()));

        zDot[4].subtract(t1.set(alpha).multiply(gamma*(n-2)).multiply(t2.set(z[0]).add(z[1]).multiply(pm)));
        zDot[4].subtract(t1.set(alpha).conjugate().multiply(gamma*(n-2)).multiply(t2.set(z[0]).multiply(2*pm)));
    }

    @Override
    public int getDimension() {
        return 5;
    }

}
