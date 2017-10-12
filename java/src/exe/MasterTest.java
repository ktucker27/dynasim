package exe;

import utils.DynaComplex;
import utils.TPSOperator;
import utils.TPSOperator.PauliOp;

public class MasterTest {

    /**
     * @param args
     * @throws InterruptedException 
     */
    public static void main(String[] args) throws InterruptedException {
        int n = 2;
        DynaComplex t1 = new DynaComplex(0,0);
        
        TPSOperator rho = new TPSOperator(n);
        rho.set(1.0/Math.pow(2,n));
        
//        System.out.println(rho);
        
        rho.trace(t1);
        System.out.println("trace: " + t1);
        
        rho.pauliLeft(PauliOp.Z, 0);
        System.out.println(rho);
        
//        for(int i = 0; i < 100; ++i) {
//            Thread.sleep(1000);
//        }
        
//        TPSOperator sum = new TPSOperator(n);
//        sum.set(0);
//        t1.set(1,0);
//        sum.pauliLeft(PauliOp.Z, 1, t1, rho);
//        System.out.println(sum);
    }

}
