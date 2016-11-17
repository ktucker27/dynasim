package exe;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.Scanner;

import ode.SynchCumulantlODEs;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.complex.ComplexFormat;

import coupling.ConstCoupling;

public class CumulantODEsTest {

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
        Complex[] z = new Complex[dim];
        ComplexFormat format = new ComplexFormat();
        Scanner inputStream = new Scanner(new File("/Users/kristophertucker/output/z_in.txt"));
        int idx = 0;
        String line;
        while(inputStream.hasNext()) {
            line = inputStream.nextLine();
            z[idx] = format.parse(line);
            ++idx;
        }
        
        double[] d = new double[n];
        d[0] = -4.8676032;
        d[1] = 8.8587453;
        d[2] = -5.6883997;
        
        ConstCoupling coupling = new ConstCoupling(1.0, 0.2);
        
        SynchCumulantlODEs codes = new SynchCumulantlODEs(n, gamma, w, coupling, d);

        // Calculate complex derivative values
        Complex[] zDot = new Complex[dim];
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
