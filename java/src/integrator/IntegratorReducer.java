package integrator;

import utils.ODESolution;

public interface IntegratorReducer {
    public void reduce(ODESolution soln);
}
