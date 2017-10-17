package utils;

/**
 * 
 * @author kristophertucker
 *
 * Class representing a tensor product space operator.  Offers fast routines
 * for applying strings of Pauli operators.
 * 
 */
public class TPSOperator {

    int n;
    int dim;
    DynaComplex[] vals;
    DynaComplex t1;
    
    public enum PauliOp {
        IDENTITY,
        PLUS,
        MINUS,
        X,
        Y,
        Z
    }
    
    public TPSOperator(int n) {
        this.n = n;
        dim = SynchUtils.pow(2, n);
        vals = new DynaComplex[dim*dim];
        t1 = new DynaComplex();
        
        for(int i = 0; i < dim; ++i) {
            for(int j = 0; j < dim; ++j) {
                vals[i*dim + j] = new DynaComplex(0,0);
            }
        }
    }
    
    public TPSOperator(DynaComplex[] z) {
        this.dim = (int)Math.sqrt(z.length);
        this.n = (int)(Math.log(dim)/Math.log(2));
        vals = z;
        t1 = new DynaComplex();
    }
    
    public TPSOperator(TPSOperator rhs) {
        this.n = rhs.n;
        this.dim = rhs.dim;
        this.vals = new DynaComplex[dim*dim];
        set(rhs);
        t1 = new DynaComplex();
    }
    
    public void set(DynaComplex z) {
        for(int i = 0; i < dim; ++i) {
            for(int j = 0; j < dim; ++j) {
                vals[i*dim + j].set(z);
            }
        }
    }
    
    public void set(double x) {
        for(int i = 0; i < dim; ++i) {
            for(int j = 0; j < dim; ++j) {
                vals[i*dim + j].set(x,0);
            }
        }
    }
    
    public void set(TPSOperator rhs) {
        for(int i = 0; i < dim; ++i) {
            for(int j = 0; j < dim; ++j) {
                vals[i*dim + j].set(rhs.getVal(i,j));
            }
        }
    }
    
    public int getN() {
        return n;
    }
    
    public DynaComplex getVal(int i, int j) {
        return vals[i*dim + j];
    }
    
    public DynaComplex[] getVals() {
        return vals;
    }
    
    public void multiply(DynaComplex z) {
        for(int i = 0; i < dim; ++i) {
            for(int j = 0; j < dim; ++j) {
                this.vals[i*dim + j].multiply(z);
            }
        }
    }
    
    public void multiply(double x) {
        for(int i = 0; i < dim; ++i) {
            for(int j = 0; j < dim; ++j) {
                this.vals[i*dim + j].multiply(x);
            }
        }
    }
    
    public void add(TPSOperator rhs) {
        for(int i = 0; i < dim; ++i) {
            for(int j = 0; j < dim; ++j) {
                this.vals[i*dim + j].add(rhs.getVal(i,j));
            }
        }
    }
    
    public void subtract(TPSOperator rhs) {
        for(int i = 0; i < dim; ++i) {
            for(int j = 0; j < dim; ++j) {
                this.vals[i*dim + j].subtract(rhs.getVal(i,j));
            }
        }
    }
    
    public void trace(DynaComplex sum) {
        sum.set(0,0);
        for(int i = 0; i < dim; ++i) {
            sum.add(vals[i*dim + i]);
        }
    }
    
//    public void pauliLeft(PauliOp[] op) {
//        for(int i = 0; i < n; ++i) {
//            if(op[i] == PauliOp.IDENTITY) continue;
//            pauliLeft(op[i], i);
//        }
//    }
    
//    public void pauliLeft(PauliOp[] op, TPSOperator rho) {
//        for(int i = 0; i < n; ++i) {
//            if(op[i] == PauliOp.IDENTITY) continue;
//            pauliLeft(op[i], i, rho);
//        }
//    }
    
//    public void pauliRight(PauliOp[] op) {
//        for(int i = 0; i < n; ++i) {
//            if(op[i] == PauliOp.IDENTITY) continue;
//            pauliRight(op[i], i);
//        }
//    }
    
//    public void pauliRight(PauliOp[] op, TPSOperator rho) {
//        for(int i = 0; i < n; ++i) {
//            if(op[i] == PauliOp.IDENTITY) continue;
//            pauliRight(op[i], i, rho);
//        }
//    }
    
    public void pauliLeft(PauliOp op, int idx) {
        int mask = SynchUtils.pow(2,idx);
        for(int i = 0; i < dim; ++i) {
            for(int j = 0; j < dim; ++j) {
                switch(op) {
                case PLUS:
                    if((i & mask) != 0) {
                       vals[(i - mask)*dim + j].set(vals[i*dim + j]);
                       vals[i*dim + j].set(0,0);
                    }
                    break;
                case MINUS:
                    if((i & mask) == 0) {
                        vals[(i + mask)*dim + j].set(vals[i*dim + j]);
                        vals[i*dim + j].set(0,0);
                     }
                    break;
                case X:
                    if((i & mask) == 0) {
                        t1.set(vals[i*dim + j]);
                        vals[i*dim + j].set(vals[(i + mask)*dim + j]);
                        vals[(i + mask)*dim + j].set(t1);
                    }
                    break;
                case Y:
                    if((i & mask) == 0) {
                        t1.set(vals[i*dim + j]);
                        vals[i*dim + j].set(vals[(i + mask)*dim + j]);
                        vals[(i + mask)*dim + j].set(t1);
                        
                        t1.set(0,1);
                        vals[i*dim + j].multiply(t1);
                        vals[(i + mask)*dim + j].multiply(t1.multiply(-1));
                    }
                    break;
                case Z:
                    if((i & mask) != 0) {
                        vals[i*dim + j].multiply(-1.0);
                    }
                    break;
                case IDENTITY:
                    break;
                default:
                    throw new UnsupportedOperationException("TPSOperator does not yet support left multiplication by " + op.name());
                }
            }
        }
    }
    
    public void pauliLeft(PauliOp op, int idx, DynaComplex z, final TPSOperator rho) {
        int mask = SynchUtils.pow(2,idx);
        for(int i = 0; i < dim; ++i) {
            for(int j = 0; j < dim; ++j) {
                switch(op) {
                case PLUS:
                    if((i & mask) != 0) {
                       vals[(i - mask)*dim + j].add(t1.set(rho.vals[i*dim + j]).multiply(z));
                    }
                    break;
                case MINUS:
                    if((i & mask) == 0) {
                        vals[(i + mask)*dim + j].add(t1.set(rho.vals[i*dim + j]).multiply(z));
                     }
                    break;
                case Z:
                    if((i & mask) != 0) {
                        vals[i*dim + j].subtract(t1.set(rho.vals[i*dim + j]).multiply(z));
                    } else {
                        vals[i*dim + j].add(t1.set(rho.vals[i*dim + j]).multiply(z));
                    }
                    break;
                case IDENTITY:
                    vals[i*dim + j].add(t1.set(rho.vals[i*dim + j]).multiply(z));
                    break;
                default:
                    throw new UnsupportedOperationException("TPSOperator does not yet support left multiplication by " + op.name());
                }
            }
        }
    }
    
    public void pauliRight(PauliOp op, int idx) {
        int mask = SynchUtils.pow(2,idx);
        for(int i = 0; i < dim; ++i) {
            for(int j = 0; j < dim; ++j) {
                switch(op) {
                case PLUS:
                    if((j & mask) == 0) {
                        vals[i*dim + j + mask].set(vals[i*dim + j]);
                        vals[i*dim + j].set(0,0);
                    }
                    break;
                case MINUS:
                    if((j & mask) != 0) {
                        vals[i*dim + j - mask].set(vals[i*dim + j]);
                        vals[i*dim + j].set(0,0);
                    }
                    break;
                case Z:
                    if((j & mask) != 0) {
                        vals[i*dim + j].multiply(-1.0);
                    }
                    break;
                case IDENTITY:
                    break;
                default:
                    throw new UnsupportedOperationException("TPSOperator does not yet support right multiplication by " + op.name());
                }
            }
        }
    }
    
    public void pauliRight(PauliOp op, int idx, DynaComplex z, final TPSOperator rho) {
        int mask = SynchUtils.pow(2,idx);
        for(int i = 0; i < dim; ++i) {
            for(int j = 0; j < dim; ++j) {
                switch(op) {
                case PLUS:
                    if((j & mask) == 0) {
                        vals[i*dim + j + mask].add(t1.set(rho.vals[i*dim + j]).multiply(z));
                    }
                    break;
                case MINUS:
                    if((j & mask) != 0) {
                        vals[i*dim + j - mask].add(t1.set(rho.vals[i*dim + j]).multiply(z));
                    }
                    break;
                case Z:
                    if((j & mask) != 0) {
                        vals[i*dim + j].subtract(t1.set(rho.vals[i*dim + j]).multiply(z));
                    } else {
                        vals[i*dim + j].add(t1.set(rho.vals[i*dim + j]).multiply(z));
                    }
                    break;
                case IDENTITY:
                    vals[i*dim + j].add(t1.set(rho.vals[i*dim + j]).multiply(z));
                    break;
                default:
                    throw new UnsupportedOperationException("TPSOperator does not yet support right multiplication by " + op.name());
                }
            }
        }
    }
    
    public String toString() {
        String ans = "";
        for(int i = 0; i < dim; ++i) {
            for(int j = 0; j < dim; ++j) {
                ans += getVal(i,j).toString();
                if(j < dim - 1) ans += ",  ";
            }
            ans += "\n";
        }
        
        return ans;
    }
}
