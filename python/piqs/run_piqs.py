from qutip import *
from piqs import *
import matplotlib.pyplot as plt
import argparse
import time
import math

class RunParams:
    def __init__(self, args):
        self.n     = int(args.n)
        self.o     = float(args.o)
        self.chi   = float(args.chi)
        self.gamma = float(args.gamma)
        self.gs    = float(args.gs)
        self.gel   = float(args.gel)
        self.tmax  = float(args.tmax)
        self.it    = float(args.it)
        self.ip    = float(args.ip)

    def to_string(self):
        return 'n     = ' + str(self.n) + '\n' + \
               'o     = ' + str(self.o) + '\n' + \
               'chi   = ' + str(self.chi) + '\n' + \
               'gamma = ' + str(self.gamma) + '\n' + \
               'gs    = ' + str(self.gs) + '\n' + \
               'gel   = ' + str(self.gel) + '\n' + \
               'tmax  = ' + str(self.tmax) + '\n' + \
               'it    = ' + str(self.it) + '\n' + \
               'ip    = ' + str(self.ip)

def parse_args():
    parser = argparse.ArgumentParser(description="Simulate the collective driven dissapative spin system using PIQS (QuTIP)")
    
    parser.add_argument("outfile", help="Full path of the output expected value file")

    parser.add_argument("-n", help="Number of atoms", action="store", required=True)
    parser.add_argument("-o", help="Coherent drive", action="store", required=True)
    parser.add_argument("-chi", help="Elastic interaction rate", action="store", required=True)
    parser.add_argument("-gamma", help="Collective decay rate", action="store", required=True)
    parser.add_argument("-gs", help="Individual atom decay rate (default 0)", action="store", required=False, default=0)
    parser.add_argument("-gel", help="Dephasing rate (default 0)", action="store", required=False, default=0)
    parser.add_argument("-tmax", help="Simulation time", action="store", required=True)
    parser.add_argument("-it", help="Initial CSS theta (default = pi/2)", action="store", required=False, default=math.pi*0.5)
    parser.add_argument("-ip", help="Initial CSS phi (default = pi)", action="store", required=False, default=math.pi)

    args = parser.parse_args()

    return [RunParams(args), args.outfile]

def write_results_csv(filepath, t, results):
    '''Assumes collective, real valued expectations occupy results.expect[0:3]'''
    all_data = np.vstack((t, results.expect[0], results.expect[1], results.expect[2]))
    for ii in range(3, len(results.expect)):
        all_data = np.vstack((all_data, results.expect[ii].real, results.expect[ii].imag))
    
    np.savetxt(filepath, np.transpose(all_data), delimiter=',')

def run_sim(outfile, run_params):

    # Define collective spin operators
    jx, jy, jz = jspin(run_params.n)
    jp = jspin(run_params.n, "+")
    jm = jspin(run_params.n, "-")

    # Define Hamiltonian
    ham = run_params.chi*jp*jm + run_params.o*jx

    # Get the Liouvillian
    ensemble = Dicke(run_params.n, hamiltonian=ham, collective_emission=run_params.gamma, emission=run_params.gs, dephasing=run_params.gel)
    liouv = ensemble.liouvillian()

    # Set the initial conditions
    #rho0 = ground(run_params.n)
    rho0 = css(run_params.n, run_params.it, run_params.ip, coordinates='polar')

    # Set the simulation parameters
    # TODO - Make nt configurable
    nt = 10001
    tmax = run_params.tmax
    t = np.linspace(0, tmax, nt)

    # Run the solver
    tic = time.time()
    result = mesolve(liouv, rho0, t, [], e_ops = [jx, jy, jz, jx**2, jx*jy, jx*jz, jy*jx, jy**2, jy*jz, jz*jx, jz*jy, jz**2],
                     options = Options(store_states=True))
    toc = time.time()
    print('Simulation time:', toc - tic, 'seconds')

    # Write results to the output file
    write_results_csv(outfile, t, result)

def main():
    [run_params, outfile] = parse_args()

    print('Running with params:')
    print(run_params.to_string())
    print('Writing results to:', outfile)

    run_sim(outfile, run_params)

if __name__ == "__main__":
    main()