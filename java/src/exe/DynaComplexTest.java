package exe;

import org.apache.commons.math3.complex.Complex;

import utils.DynaComplex;

public class DynaComplexTest {

    /**
     * @param args
     */
    public static void main(String[] args) {
        Complex z = new Complex(1, 0);
        DynaComplex dz = new DynaComplex(1, 0);
        double s1 = 0.0;
        double s2 = 0.0;
        
        Complex a = new Complex(1.0/Math.sqrt(2.0), 1.0/Math.sqrt(2.0));
        DynaComplex da = new DynaComplex(1.0/Math.sqrt(2.0), 1.0/Math.sqrt(2.0));
        double sa = 1.0/Math.sqrt(2.0);
        
        int N = 100000000;
        
        long startTime = System.nanoTime();
        for(int i = 0; i < N; ++i) {
            //z = z.add(a);
            //dz.add(da);
            //s1 += sa;
            //s2 += sa;
            
            //dz.set(da);

            //z = z.multiply(a);
            //dz.multiply(da);
            //s1 += sa;
            //s2 += sa;
        }
        long endTime = System.nanoTime();
        
        System.out.println("Run time: " + (endTime - startTime)/1.0e9 + " seconds | (" + z.getReal() + ", " + z.getImaginary() + ")");
        System.out.println("Run time: " + (endTime - startTime)/1.0e9 + " seconds | (" + dz.getReal() + ", " + dz.getImaginary() + ")");
        System.out.println(s1 + ", " + s2);
    }

}
