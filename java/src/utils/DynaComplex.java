package utils;

public class DynaComplex {

    double x;
    double y;
    double temp;
    
    public DynaComplex() {
        x = 0.0;
        y = 0.0;
    }
    
    public DynaComplex(DynaComplex z) {
        set(z);
    }
    
    public DynaComplex(double x) {
        this.x = x;
        this.y = 0.0;
    }
    
    public DynaComplex(double x, double y) {
        this.x = x;
        this.y = y;
    }
    
    public DynaComplex set(double x, double y) {
        this.x = x;
        this.y = y;
        return this;
    }
    
    public DynaComplex set(DynaComplex z) {
        x = z.x;
        y = z.y;
        return this;
    }
    
    public double getReal() {return x;}
    public double getImaginary() {return y;}
    
    public void setReal(double x) {
        this.x = x;
    }
    
    public void setImaginary(double y) {
        this.y = y;
    }
    
    public DynaComplex conjugate() {
        y *= -1;
        return this;
    }
    
    public DynaComplex multiply(DynaComplex z) {
        temp = x*z.y + y*z.x;
        x = x*z.x - y*z.y;
        y = temp;
        return this;
    }
    
    public DynaComplex multiply(double c) {
        x *= c;
        y *= c;
        return this;
    }
    
    public DynaComplex add(DynaComplex z) {
        x += z.x;
        y += z.y;
        return this;
    }
    
    public DynaComplex add(double c) {
        x += c;
        return this;
    }
    
    public DynaComplex subtract(DynaComplex z) {
        x -= z.x;
        y -= z.y;
        return this;
    }
    
    public double mod() {
        return Math.sqrt(x*x + y*y);
    }
    
    public String toString() {
        return x + " + 1i*" + y;
    }
}
