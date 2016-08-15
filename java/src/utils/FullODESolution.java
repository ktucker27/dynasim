package utils;

public class FullODESolution implements ODESolution {

    double[] y;
    
    public FullODESolution() {
        super();
        y = null;
    }
    
    public double[] getSolution() {
        return y;
    }
    
    @Override
    public void setSolution(double[] y) {
        this.y = y;
    }

}
