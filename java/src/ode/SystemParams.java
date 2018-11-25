package ode;

import java.io.File;

import utils.DynaComplex;

public class SystemParams implements Comparable<SystemParams> {
    private int n;
    private double gamma;
    private double w;
    private double omega;
    private double gel;
    private double delta;
    private double faa;
    private double fab;
    private double gaa;
    private double gab;
    private double[] d;
    private DynaComplex alpha;
    
    public SystemParams(int n, double gamma, double w, double omega, double gel, double delta,
                        double faa, double fab, double gaa, double gab, double[] d) {
        this.n = n;
        this.gamma = gamma;
        this.w = w;
        this.omega = omega;
        this.gel = gel;
        this.delta = delta;
        this.faa = faa;
        this.fab = fab;
        this.gaa = gaa;
        this.gab = gab;
        this.d = new double[d.length];
        for(int i = 0; i < d.length; ++i) {
            this.d[i] = d[i];
        }
        this.alpha = new DynaComplex(this.fab, this.gab);
    }
    
    public SystemParams(int n, double gamma, double w, double delta, DynaComplex alpha, double[] d) {
        this.n = n;
        this.gamma = gamma;
        this.w = w;
        this.omega = 0.0;
        this.gel = 0.0;
        this.delta = delta;
        this.faa = alpha.getReal();
        this.fab = alpha.getReal();
        this.gaa = 0.0;
        this.gab = alpha.getImaginary();
        this.d = new double[d.length];
        for(int i = 0; i < d.length; ++i) {
            this.d[i] = d[i];
        }
        this.alpha = new DynaComplex(alpha);
    }
    
    public SystemParams(SystemParams rhs) {
        this.n = rhs.n;
        this.gamma = rhs.gamma;
        this.w = rhs.w;
        this.omega = rhs.omega;
        this.gel = rhs.gel;
        this.delta = rhs.delta;
        this.faa = rhs.faa;
        this.fab = rhs.fab;
        this.gaa = rhs.gaa;
        this.gab = rhs.gab;
        this.d = new double[rhs.d.length];
        for(int i = 0; i < d.length; ++i) {
            this.d[i] = rhs.d[i];
        }
        this.alpha = new DynaComplex(rhs.alpha);
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
    
    public double getOmega() {
        return omega;
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
    
    public double getFaa() {
        return faa;
    }

    public double getFab() {
        return fab;
    }
    
    public double getGaa() {
        return gaa;
    }
    
    public double getGab() {
        return gab;
    }
    
    public DynaComplex getAlpha() {
        return alpha;
    }

    public double[] getD() {
        return d;
    }
    
    public File getResultsDir() {
        String resDirStr = String.format("/N%d/D%.1f/g%.1f", n, delta, gab).replace('.', 'p');
        if(gel != 0) {
            resDirStr += String.format("/gel%.2f", gel).replace('.','p');
        }
        return new File(resDirStr);
    }
    
    public String getFilename() {
        String pumpStr = String.format("o%.2f", omega).replace('.','p');;
        if(w != 0.0) {
            pumpStr = String.format("w%.2f", w).replace('.','p');            
        }
        
        String filename = String.format("N%d_D%.1f_g%.1f_%s_f%.1f", n, delta, gab, pumpStr, fab).replace('.', 'p');
        
        if(gel != 0) {
            filename += String.format("_gel%.3f", gel).replace('.','p');
        }
        
        if(faa != fab) {
            filename += String.format("_faa%.2f", faa).replace('.','p');
        }
        
        if(gaa != 0) {
            filename += String.format("_gaa%.2f", gaa).replace('.','p');
        }
        
        filename += ".txt";
        return filename;
    }
    
    public String getMcwfFilename() {
        String pumpStr = String.format("o%.2f", omega).replace('.','p');
        
        String filename = String.format("N%d_chi%.1f_%s_f%.1f", n, gab/2, pumpStr, fab).replace('.', 'p');
        
        if(gel != 0) {
            filename += String.format("_gel%.3f", gel).replace('.','p');
        }
        
        if(faa != fab) {
            filename += String.format("_gs%.2f", gamma*(faa - fab)).replace('.','p');
        }
        
        return filename;
    }
    
    public String toString() {
        String out = "N = " + n + "\n";
        out += "W = " + w + "\n";
        out += "Omega = " + omega + "\n";
        out += "gel = " + gel + "\n";
        out += "Delta = " + delta + "\n";
        out += "faa = " + faa + "\n";
        out += "fab = " + fab + "\n";
        out += "gaa = " + gaa + "\n";
        out += "gab = " + gab + "\n";
        out += "Gamma = " + gamma + "\n";
//        for(int j = 0; j < d.length; ++j) {
//            out += "d[" + j + "]: " + d[j] + "\n";
//        }

        return out;
    }
    
    public static String getHeader() {
        return "n, w, omega, gel, delta, faa, fab, gaa, gab, gamma";
    }
    
    public String getLine() {
        String out = n + ", ";
        out += w + ", ";
        out += omega + ", ";
        out += gel + ", ";
        out += delta + ", ";
        out += faa + ", ";
        out += fab + ", ";
        out += gaa + ", ";
        out += gab + ", ";
        out += gamma;
        
        return out;
    }

    @Override
    public int compareTo(SystemParams rhs) {
        if(this.equals(rhs)) return 0;
        
        int sgn;
        if((sgn = (int)Math.signum(this.n - rhs.n)) != 0) return sgn;
        if((sgn = (int)Math.signum(this.gamma - rhs.gamma)) != 0) return sgn;
        if((sgn = (int)Math.signum(this.w - rhs.w)) != 0) return sgn;
        if((sgn = (int)Math.signum(this.omega - rhs.omega)) != 0) return sgn;
        if((sgn = (int)Math.signum(this.gel - rhs.gel)) != 0) return sgn;
        if((sgn = (int)Math.signum(this.delta - rhs.delta)) != 0) return sgn;
        if((sgn = (int)Math.signum(this.faa - rhs.faa)) != 0) return sgn;
        if((sgn = (int)Math.signum(this.fab - rhs.fab)) != 0) return sgn;
        if((sgn = (int)Math.signum(this.gaa - rhs.gaa)) != 0) return sgn;
        if((sgn = (int)Math.signum(this.gab - rhs.gab)) != 0) return sgn;
        
        for(int i = 0; i < this.n; ++i) {
            if((sgn = (int)Math.signum(this.d[i] - rhs.d[i])) != 0) return sgn;
        }
        
        // This should never happen
        return 0;
    }
    
    @Override
    public boolean equals(Object arhs) {
        if(arhs.getClass() != SystemParams.class) return false;
        
        SystemParams rhs = (SystemParams)arhs;
        if (this.n != rhs.n ||
            this.gamma != rhs.gamma ||
            this.w != rhs.w ||
            this.omega != rhs.omega ||
            this.gel != rhs.gel ||
            this.delta != rhs.delta ||
            this.faa != rhs.faa ||
            this.fab != rhs.fab ||
            this.gaa != rhs.gaa ||
            this.gab != rhs.gab) {
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
