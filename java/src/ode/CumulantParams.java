package ode;

import java.io.File;

import utils.DynaComplex;

public class CumulantParams {
    private int n;
    private double gamma;
    private double w;
    private double delta;
    private DynaComplex alpha;
    private double[] d;
    File resDir;
    
    public CumulantParams(int n, double gamma, double w, double delta, DynaComplex alpha, double[] d) {
        this.n = n;
        this.gamma = gamma;
        this.w = w;
        this.delta = delta;
        this.alpha = new DynaComplex(alpha);
        this.d = d;
        
        resDir = new File(String.format("/N%d/D%.1f/g%.1f", n, delta, alpha.getImaginary()).replace('.', 'p'));
    }

    public int getN() {
        return n;
    }

    public double getGamma() {
        return gamma;
    }

    public double getW() {
        return w;
    }
    
    public double getDelta() {
        return delta;
    }

    public DynaComplex getAlpha() {
        return alpha;
    }

    public double[] getD() {
        return d;
    }
    
    public File getResultsDir() {
        return resDir;
    }
    
    public String toString() {
        String out = "N = " + n + "\n";
        out += "Delta = " + delta + "\n";
        out += "W = " + w + "\n";
        out += "g = " + alpha.getImaginary() + "\n";
        out += "f = " + alpha.getReal() + "\n";
        out += "Gamma = " + gamma + "\n";
        return out;
    }
}
