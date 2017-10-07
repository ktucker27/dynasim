package exe;

import org.apache.commons.math3.complex.Complex;

import ode.CumulantAllToAllODEs;
import ode.CumulantParams;
import ode.SynchMeanFieldODEs;
import utils.DynaComplex;
import utils.SynchUtils;
import eval.CumulantEval;
import eval.MeanFieldEval;

public class MeanFieldTest {

    /**
     * @param args
     */
    public static void main(String[] args) {
        int n = 2;
        double gamma = 1.0;
        double f = 1.0;
        double g = 30.0;
        double delta = 0.0;
        double w = 3.0;
        
        double[] d = new double[n];
        SynchUtils.detuneGauss(delta, d);
        
        CumulantParams dparams = new CumulantParams(n, gamma, w, delta, new DynaComplex(f,g), d);
        
        MeanFieldEval meval = new MeanFieldEval(n);
        double[] y = new double[meval.getRealDimension()];
        double[] yDot = new double[meval.getRealDimension()];
        
        CumulantEval ceval = new CumulantEval(n);
        DynaComplex[] z = new DynaComplex[ceval.getDimension()];
        DynaComplex[] zDot = new DynaComplex[ceval.getDimension()];
        for(int i = 0; i < ceval.getDimension(); ++i) {
            z[i] = new DynaComplex();
            zDot[i] = new DynaComplex();
        }

        double dtheta = 0.1;
        double dphi = 0.1;
        double maxVal = 0.0;
        for(double theta = 0.1; theta < Math.PI; theta += dtheta) {
            for(double phi = 0; phi <= 2*Math.PI; phi += dphi) {
                meval.initialize(y, theta, phi, MeanFieldEval.InitPhaseType.CONST);
                SynchMeanFieldODEs modes = new SynchMeanFieldODEs(dparams);
                modes.computeDerivatives(0, y, yDot);

                ceval.initialize(z, theta, phi, false);
                CumulantAllToAllODEs codes = new CumulantAllToAllODEs(dparams);
                codes.computeDerivatives(0, z, zDot);

                for(int i = 0; i < n; ++i) {
//                    System.out.println(zDot[ceval.getStartIdx()[2] + i].getReal() - yDot[i]);
                    maxVal = Math.max(maxVal,  Math.abs(zDot[ceval.getStartIdx()[2] + i].getReal() - yDot[i]));
                    Complex val = new Complex(0.5*yDot[n+i], 0.5*y[n+i]*yDot[2*n+i]);
                    Complex val2 = new Complex(0, y[2*n+i]);
                    val2 = val2.exp();
                    val = val.multiply(val2);
//                    System.out.println(val.getReal() - zDot[i].getReal());
//                    System.out.println(val.getImaginary() - zDot[i].getImaginary());
//                    System.out.println(y[2*n+i]);
                    maxVal = Math.max(maxVal, Math.abs(val.getReal() - zDot[i].getReal()));
                    maxVal = Math.max(maxVal, Math.abs(val.getImaginary() - zDot[i].getImaginary()));
                }
            }
        }
        
        System.out.println("maxVal: " + maxVal);
    }

}
