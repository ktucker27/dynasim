package ode;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Scanner;

import coupling.DynaComplexCoupling;
import coupling.DynaConstCoupling;
import eval.CumulantEval;
import utils.DynaComplex;

public class FOCorrelationODEsCorrected implements DynaComplexODEs {
    
    /**
     * Fourth order cumulant version of the two-time correlation solver
     */

    private int n;
    private double gamma;
    private double w;
    private DynaComplexCoupling coupling;
    private double[] d;
    
    CumulantEval eval;
    DynaComplex[] steady;
    DynaComplex[] szs;
    DynaComplex[][] zzs;
    DynaComplex[][] pms;
    
    DynaComplex[] c1;
    DynaComplex[] c2;
    
    DynaComplex t1, t2;
    DynaComplex sum, sum2, sum3;
    
    ArrayList<DynaComplex[]> corrections;
    
    /**
     * @param n number of atoms
     * @param gamma spontaneous decay rate
     * @param w pumping rate
     * @param coupling coupling coefficients
     * @param d size n vector of natural frequencies
     * @param szs size n vector of the steady state values of sigma_z
     * @throws FileNotFoundException 
     */
    public FOCorrelationODEsCorrected(SystemParams params, DynaComplex[] z) throws FileNotFoundException {
//        public FOCorrelationODEsCorrected(int n, double gamma, double w, DynaComplexCoupling coupling, double[] d, DynaComplex[] z) throws FileNotFoundException {
        super();
        this.n = params.getN();
        this.gamma = params.getGamma();
        this.w = params.getW();
        this.coupling = new DynaConstCoupling(params.getAlpha().getReal(), params.getAlpha().getImaginary());
        this.d = params.getD();
        
        eval = new CumulantEval(n);
        this.steady = z;
        
        szs = new DynaComplex[n];
        zzs = new DynaComplex[n][n];
        pms = new DynaComplex[n][n];
        
        t1 = new DynaComplex(0, 0);
        t2 = new DynaComplex(0, 0);
        
        sum = new DynaComplex(0, 0);
        sum2 = new DynaComplex(0, 0);
        sum3 = new DynaComplex(0, 0);
        
        initConstants();
    }
    
    private void initConstants() throws FileNotFoundException {
        c1 = new DynaComplex[n];
        c2 = new DynaComplex[n];
        for(int i = 0; i < n; ++i) {
            c1[i] = new DynaComplex(-0.5*(gamma+w), -d[i]);
            c2[i] = new DynaComplex(-1.5*(gamma+w), -d[i]);
            
//            szs[i] = new DynaComplex(eval.getSingle(2, i, steady, t1));
//            szs[i] = new DynaComplex(0.487907742998352, 0);
            szs[i] = new DynaComplex(0.452811301424729, 0);
            
            for(int j = 0; j < n; ++j) {
                if(i == j) continue;
                
//                zzs[i][j] = new DynaComplex(eval.getDouble(2, 2, i, j, steady, t1));
//                pms[i][j] = new DynaComplex(eval.getDouble(0, 1, i, j, steady, t1));
//                zzs[i][j] = new DynaComplex(0.301087314662273, 0);
//                pms[i][j] = new DynaComplex(0.093410214168040, 0);
                zzs[i][j] = new DynaComplex(0.268825428205640, 0);
                pms[i][j] = new DynaComplex(0.091992936609545, 0);
            }
        }
        
        corrections = null;
        corrections = new ArrayList<DynaComplex[]>();
        Scanner fileReader = new Scanner(new File("/Users/tuckerkj/output/foc/foc.csv"));
        while(fileReader.hasNext()) {
            String line = fileReader.nextLine();
            String[] tokens = line.split(",");
            DynaComplex[] cline = new DynaComplex[(tokens.length-1)/2];
            int idx = 0;
            for(int i = 1; i < tokens.length; i += 2) {
                double rp = Double.parseDouble(tokens[i]);
                double ip = Double.parseDouble(tokens[i+1]);
                cline[idx] = new DynaComplex(rp, ip);
//                cline[idx] = new DynaComplex(0, 0);
                ++idx;
            }
            corrections.add(cline);
        }
        fileReader.close();
        System.out.println(corrections.get(0).length + ", " + corrections.size());
    }
    
    private int getTripleIndex(int a, int b, int c) {
        int idx = n*n + c*n*(n-1) + a*(n-1) + b;
        
        if(b > a) --idx;
        
        return idx;
    }
    
    @Override
    public void computeDerivatives(double t, DynaComplex[] z, DynaComplex[] zDot) {
        int idx = 0;
        int tidx = (int)Math.round(t/0.0001);
        
        // <sigma_a^+(t + tau) sigma_b^-(t)>
        for(int a = 0; a < n; ++a) {
            for(int b = 0; b < n; ++b) {
                
                sum.set(0,0);
                for(int j = 0; j < n; ++j) {
                    if(j == a) continue;
                    
                    t1.set(z[getTripleIndex(a,j,b)]).multiply(t2.set(coupling.getAlpha(a, j)).conjugate());
                    sum.add(t1);
                }
                sum.multiply(gamma*0.5);
                
                zDot[idx].set(c1[a]).multiply(z[idx]);
                zDot[idx].add(sum);
                
                ++idx;
            }
        }
        
        // <sigma_a^z(t + tau) sigma_b^+(t + tau) sigma_c^-(t)>
        DynaComplex t3 = new DynaComplex();
        for(int c = 0; c < n; ++c) {
            for(int a = 0; a < n; ++a) {
                for(int b = 0; b < n; ++b) {
                    if(a == b) continue;
                    
                    sum.set(0,0);
                    for(int j = 0; j < n; ++j) {
                        if(j == a || j == b) continue;
                        
                        t1.set(z[getTripleIndex(b, j, c)]).multiply(szs[a]);
                        t1.add(t2.set(z[getTripleIndex(a, j, c)]).multiply(szs[b]));
                        t1.add(t2.set(z[j*n+c]).multiply(zzs[a][b]));
                        t1.subtract(t2.set(z[j*n+c]).multiply(szs[a]).multiply(szs[b]).multiply(2.0));
                        
                        // Correction
                        if(corrections != null) {
                            if(c == a) {
                                t3.set(corrections.get(tidx)[0]);
                            } else if(c == b) {
                                t3.set(corrections.get(tidx)[1]);
                            } else if(c == j) {
                                t3.set(corrections.get(tidx)[2]);
                            } else {
                                t3.set(corrections.get(tidx)[3]);
                            }
                            t1.add(t3);
                        }
                        
                        t1.multiply(t2.set(coupling.getAlpha(b, j)).conjugate());
                        sum.add(t1);
                    }
                    sum.multiply(0.5*gamma);
                    
                    sum2.set(0,0);
                    for(int j = 0; j < n; ++j) {
                        if(j == a || j == b) continue;
                        
                        t1.set(z[a*n+c]).multiply(pms[b][j]);
                        t1.add(t2.set(z[b*n+c]).multiply(pms[a][j]));
                        
                        // Correction
                        if(corrections != null) {
                            if(c == a) {
                                t3.set(corrections.get(tidx)[4]);
                            } else if(c == b) {
                                t3.set(corrections.get(tidx)[4]);
                            } else if(c == j) {
                                t3.set(corrections.get(tidx)[5]);
                            } else {
                                t3.set(corrections.get(tidx)[6]);
                            }
                            t1.add(t3);
                        }
                        
                        t1.multiply(coupling.getAlpha(a, j));
                        sum2.add(t1);
                    }
                    sum2.multiply(-1.0*gamma);
                    
                    sum3.set(0,0);
                    for(int j = 0; j < n; ++j) {
                        if(j == a || j == b) continue;
                        
                        t1.set(z[j*n+c]).multiply(pms[b][a]);
                        t1.add(t2.set(z[b*n+c]).multiply(pms[j][a]));
                        
                        // Correction
                        if(corrections != null) {
                            if(c == a) {
                                t3.set(corrections.get(tidx)[5]);
                            } else if(c == b) {
                                t3.set(corrections.get(tidx)[4]);
                            } else if(c == j) {
                                t3.set(corrections.get(tidx)[4]);
                            } else {
                                t3.set(corrections.get(tidx)[6]);
                            }
                            t1.add(t3);
                        }
                        
                        t1.multiply(t2.set(coupling.getAlpha(a, j)).conjugate());
                        sum3.add(t1);
                    }
                    sum3.multiply(-1.0*gamma);
                    
                    zDot[idx].set(c2[b]).multiply(z[getTripleIndex(a, b, c)]);
                    zDot[idx].subtract(t1.set(z[b*n+c]).multiply(gamma - w));
                    zDot[idx].subtract(t1.set(z[a*n+c]).multiply(coupling.getAlpha(a, b)).multiply(0.5*gamma));
                    zDot[idx].subtract(t1.set(z[getTripleIndex(b, a, c)]).multiply(gamma*coupling.getAlpha(a, b).getReal()));
                    zDot[idx].add(sum).add(sum2).add(sum3);
                    
                    ++idx;
                }
            }
        }
    }

    @Override
    public int getDimension() {
        return n*n + n*n*(n-1);
    }

}
