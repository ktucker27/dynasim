package ode;

import java.io.File;

import utils.DynaComplex;

public class CumulantParams implements Comparable<CumulantParams> {
    private int n;
    private double gamma;
    private double w;
    private double gel;
    private double delta;
    private DynaComplex alpha;
    private double[] d;
    
    public CumulantParams(int n, double gamma, double w, double delta, DynaComplex alpha, double[] d) {
        this.n = n;
        this.gamma = gamma;
        this.w = w;
        this.gel = 0.0;
        this.delta = delta;
        this.alpha = new DynaComplex(alpha);
        this.d = new double[d.length];
        for(int i = 0; i < d.length; ++i) {
            this.d[i] = d[i];
        }
    }
    
    public CumulantParams(CumulantParams rhs) {
        this.n = rhs.n;
        this.gamma = rhs.gamma;
        this.w = rhs.w;
        this.gel = rhs.gel;
        this.delta = rhs.delta;
        this.alpha = new DynaComplex(rhs.alpha);
        this.d = new double[rhs.d.length];
        for(int i = 0; i < d.length; ++i) {
            this.d[i] = rhs.d[i];
        }
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
    
    public void setGel(double gel) {
        this.gel = gel;
    }
    
    public double getGel() {
        return gel;
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
        String resDirStr = String.format("/N%d/D%.1f/g%.1f", n, delta, alpha.getImaginary()).replace('.', 'p');
        if(gel != 0) {
            resDirStr += String.format("/gel%.2f", gel).replace('.','p');
        }
        return new File(resDirStr);
    }
    
    public String getFilename() {
        String filename = String.format("N%d_D%.1f_g%.1f_w%.2f_f%.1f", n, delta, alpha.getImaginary(), w, alpha.getReal()).replace('.', 'p');
        if(gel != 0) {
            filename += String.format("_gel%.2f", gel).replace('.','p');
        }
        filename += ".txt";
        return filename;
    }
    
    public String toString() {
        String out = "N = " + n + "\n";
        out += "Delta = " + delta + "\n";
        out += "W = " + w + "\n";
        out += "g = " + alpha.getImaginary() + "\n";
        out += "f = " + alpha.getReal() + "\n";
        out += "Gamma = " + gamma + "\n";
        out += "gel = " + gel + "\n";
        return out;
    }

    @Override
    public int compareTo(CumulantParams rhs) {
        if(this.equals(rhs)) return 0;
        
        int sgn;
        if((sgn = (int)Math.signum(this.n - rhs.n)) != 0) return sgn;
        if((sgn = (int)Math.signum(this.gamma - rhs.gamma)) != 0) return sgn;
        if((sgn = (int)Math.signum(this.w - rhs.w)) != 0) return sgn;
        if((sgn = (int)Math.signum(this.gel - rhs.gel)) != 0) return sgn;
        if((sgn = (int)Math.signum(this.delta - rhs.delta)) != 0) return sgn;
        if((sgn = (int)Math.signum(this.alpha.getReal() - rhs.alpha.getReal())) != 0) return sgn;
        if((sgn = (int)Math.signum(this.alpha.getImaginary() - rhs.alpha.getImaginary())) != 0) return sgn;
        
        for(int i = 0; i < this.n; ++i) {
            if((sgn = (int)Math.signum(this.d[i] - rhs.d[i])) != 0) return sgn;
        }
        
        // This should never happen
        return 0;
    }
    
    @Override
    public boolean equals(Object arhs) {
        if(arhs.getClass() != CumulantParams.class) return false;
        
        CumulantParams rhs = (CumulantParams)arhs;
        if (this.n != rhs.n ||
            this.gamma != rhs.gamma ||
            this.w != rhs.w ||
            this.gel != rhs.gel ||
            this.delta != rhs.delta ||
            this.alpha.getReal() != rhs.alpha.getReal() ||
            this.alpha.getImaginary() != rhs.alpha.getImaginary()) {
            return false;
        }
        
        for(int i = 0; i < this.n; ++i) {
            if(this.d[i] != rhs.d[i]) {
                return false;
            }
        }
        
        return true;
    }
}
