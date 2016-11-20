package exe;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.Scanner;

import ode.DynaCumulantODEs;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.complex.ComplexFormat;

import utils.DynaComplex;
import coupling.DynaConstCoupling;

public class DynaCumulantODEsTest {

    /**
     * @param args
     * @throws FileNotFoundException 
     */
    public static void main(String[] args) throws FileNotFoundException {
        int n = 3;
        int dim = n*(n-1) + 3*n*(n-1)/2 + 2*n;
        double w = 2.5;
        double gamma = 1.0;
        
        // Get natural frequencies from file
        DynaComplex[] z = new DynaComplex[dim];
        ComplexFormat format = new ComplexFormat();
        Scanner inputStream = new Scanner(new File("/Users/kristophertucker/output/z_in.txt"));
        int idx = 0;
        String line;
        Complex temp;
        while(inputStream.hasNext()) {
            line = inputStream.nextLine();
            temp = format.parse(line);
            z[idx] = new DynaComplex(temp.getReal(), temp.getImaginary());
            ++idx;
        }
        
        double[] d = new double[n];
        d[0] = -4.8676032;
        d[1] = 8.8587453;
        d[2] = -5.6883997;
        
        DynaConstCoupling coupling = new DynaConstCoupling(1.0, 0.2);
        
        DynaCumulantODEs codes = new DynaCumulantODEs(n, gamma, w, coupling, d);

        // Calculate complex derivative values
        DynaComplex[] zDot = new DynaComplex[dim];
        for(int i = 0; i < dim; ++i) {
            zDot[i] = new DynaComplex(0, 0);
        }
        
        codes.computeDerivatives(0, z, zDot);
        
        // DEBUGG
        for(int i = 0; i < z.length; ++i) {
            System.out.println(z[i].getReal() + " + j*" + z[i].getImaginary());
        }
        
        System.out.println("");

        for(int i = 0; i < z.length; ++i) {
            System.out.println(zDot[i].getReal() + " + j*" + zDot[i].getImaginary());
        }
    }

}
