package exe;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.UnsupportedEncodingException;
import java.math.BigDecimal;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.Scanner;
import java.util.UUID;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.MissingOptionException;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

import mcwf.MCWFThreadPoolIntegrator;
import mcwf.MCWFWriter;
import mcwf.QuantumTrajectory;
import ode.SystemParams;
import utils.DynaComplex;

public class MCWF {

    private static final int MAX_VALIDATION_MSGS = 10;
    
    private static void printUsage(Options options) {
        System.out.println("java -jar MCWF.jar outdir [options]\n");

        System.out.println("Options (all are required unless otherwise noted):");
        
        Iterator<Option> iter = options.getOptions().iterator();
        while(iter.hasNext()) {
            Option opt = iter.next();
            System.out.printf("%5s: %s\n", opt.getOpt(), opt.getDescription());
        }
    }
    
    private static Options setupOptions() {
        Options options = new Options();
        options.addRequiredOption("n", "n", true, "Number of particles");
        options.addRequiredOption("o", "o", true, "Coherent pumping");
        options.addRequiredOption("f", "f", true, "Inelastic interaction term");
        options.addRequiredOption("faa", "faa", true, "Onsite inelastic interaction term (ind. decay term is (faa - f))");
        options.addRequiredOption("chi", "chi", true, "Elastic interaction term");
        //options.addRequiredOption("gel", "gel", true, "Single particle dephasing term");
        options.addRequiredOption("dt", "dt", true, "Time spacing");
        options.addRequiredOption("tmax", "tmax", true, "Maximum simulation time");
        options.addRequiredOption("nt", "nt", true, "Number of threads to use for integration");
        options.addRequiredOption("traj", "traj", true, "Number of trajectories");
        options.addOption("evt", true, "Time between expected value evaluations (optional, defaults to dt)");
        options.addOption("debug", false, "Enable debug features");
        options.addOption("upx", false, "Use spin-up in the X direction as the initial condition");
        options.addOption("dnx", false, "Use spin-down in the X direction as the initial condition");
        options.addOption("wf", true, "Write final trajectory states to a file. 'def' puts it in the default output directory location with prefix 'final'");
        options.addOption("icf", true, "Read initial trajectory states from a file");
        
        return options;
    }
    
    public static void main(String[] args) throws FileNotFoundException, UnsupportedEncodingException, ParseException {
        // Parse command line options
        Options options = setupOptions();

        if(args.length == 0) {
            printUsage(options);
            return;
        }
        
        CommandLineParser parser = new DefaultParser();
        CommandLine cmd;
        try {
            cmd = parser.parse(options, args);
        } catch(MissingOptionException ex) {
            System.out.println("Please specify the following command line options:");
            Iterator<?> iter = ex.getMissingOptions().iterator();
            while(iter.hasNext()) {
                System.out.println(iter.next());
            }
            return;
        }
        
        // Get the output directory
        String outdir = "";
        if(cmd.getArgList().size() > 0) {
            outdir = cmd.getArgList().get(0);
        } else {
            System.err.println("Must provide an output directory as an argument");
            return;
        }
        
        File fdir = new File(outdir);
        fdir.mkdirs();
        
        // Set simulation parameters
        int numTrajectories = Integer.parseInt(cmd.getOptionValue("traj"));
        int numThreads = Integer.parseInt(cmd.getOptionValue("nt"));
        
        int n = Integer.parseInt(cmd.getOptionValue("n"));
        double gaa = 2.0*Double.parseDouble(cmd.getOptionValue("chi"));
        double gab = gaa;
        double o = Double.parseDouble(cmd.getOptionValue("o"));
        double w = 0.0;
        double faa = Double.parseDouble(cmd.getOptionValue("faa"));
        double fab = Double.parseDouble(cmd.getOptionValue("f"));
        //double gel = Double.parseDouble(cmd.getOptionValue("gel"));
        double gel = 0.0;
        double gamma = 1.0;
        
        double dt = Double.parseDouble(cmd.getOptionValue("dt"));
        double tf = Double.parseDouble(cmd.getOptionValue("tmax"));
        int numTimes = (int)Math.ceil(tf/dt + 1);
        
        double evDelta = dt;
        if(cmd.hasOption("evt")) {
            evDelta = Double.parseDouble(cmd.getOptionValue("evt"));
        }
        
        boolean debug = cmd.hasOption("debug");
        
        long timeout = 7*24*3600; // One week
        
        double[] d = new double[n];
        for(int i = 0; i < n; ++i) {
            d[i] = 0.0;
        }
        
        SystemParams params = new SystemParams(n, gamma, w, o, gel, 0, faa, fab, gaa, gab, d);
        
        // Define the initial condition
        // TODO - Enable arbitrary Bloch sphere initial conditions
        DynaComplex[] initialState = new DynaComplex[n+1];
        if(cmd.hasOption("upx")) {
            initialState[0] = new DynaComplex(Math.sqrt(1.0/Math.pow(2.0, n)), 0);
            initialState[n] = new DynaComplex(Math.sqrt(1.0/Math.pow(2.0, n)), 0);
            initialState[1] = new DynaComplex((1.0/getAjmp(n/2, n/2 - 1))*((n/2)*initialState[0].getReal()), 0);
            int m = n/2 - 2;
            for(int i = 2; i < n; ++i) {
                initialState[i] = new DynaComplex((1.0/getAjmp(n/2,m))*((n/2)*initialState[i-1].getReal() - getAjmp(n/2, m+1)*initialState[i-2].getReal()), 0);
                --m;
            }
            
//            for(int i = 0; i < n+1; ++i) {
//                System.out.println(initialState[i]);
//            }
        } else if(cmd.hasOption("dnx")) {
            BigDecimal[] ic = new BigDecimal[n+1];
            BigDecimal bd2 = new BigDecimal(Math.sqrt(2.0));
            ic[0] = new BigDecimal("-1.0");
            for(int i = 0; i < n; ++i) {
                ic[0] = ic[0].divide(bd2, 2000, BigDecimal.ROUND_HALF_UP);
            }
            ic[n] = ic[0];
            
            ic[1] = new BigDecimal((1.0/getAjmp(n/2, n/2 - 1))*(-n/2));
            ic[1] = ic[1].multiply(ic[0]);
            ic[n-1] = ic[1];
            
            int m = n/2 - 2;
            for(int i = 2; i <= n/2; ++i) {
                BigDecimal bdt1 = new BigDecimal(-getAjmp(n/2, m+1));
                bdt1 = bdt1.multiply(ic[i-2]);
                
                BigDecimal bdt2 = new BigDecimal(-n/2);
                bdt2 = bdt2.multiply(ic[i-1]);
                
                bdt2 = bdt1.add(bdt2);
                
                ic[i] = new BigDecimal(1.0/getAjmp(n/2,m));
                ic[i] = ic[i].multiply(bdt2);
                
                if(i < n/2) {
                    ic[n-i] = ic[i];
                }
                
                --m;
            }
            
            for(int i = 0; i <= n; ++i) {
                initialState[i] = new DynaComplex(ic[i].doubleValue(), 0);
            }
            
//            initialState[0] = new DynaComplex(-Math.sqrt(1.0/Math.pow(2.0, n)), 0);
//            initialState[n] = new DynaComplex(-Math.sqrt(1.0/Math.pow(2.0, n)), 0);
//            initialState[1] = new DynaComplex((1.0/getAjmp(n/2, n/2 - 1))*((-n/2)*initialState[0].getReal()), 0);
//            initialState[n-1] = new DynaComplex((1.0/getAjmp(n/2, n/2 - 1))*((-n/2)*initialState[0].getReal()), 0);
//            int m = n/2 - 2;
//            for(int i = 2; i <= n/2; ++i) {
//                initialState[i] = new DynaComplex((1.0/getAjmp(n/2,m))*((-n/2)*initialState[i-1].getReal() - getAjmp(n/2, m+1)*initialState[i-2].getReal()), 0);
//                if(i < n/2) {
//                    initialState[n-i] = new DynaComplex((1.0/getAjmp(n/2,m))*((-n/2)*initialState[i-1].getReal() - getAjmp(n/2, m+1)*initialState[i-2].getReal()), 0);
//                }
//                --m;
//            }
            
//            for(int i = 0; i < n+1; ++i) {
//                System.out.println(initialState[i].getReal());
//            }
//            return;
        } else {
            for(int i = 0; i < n+1; ++i) {
                initialState[i] = new DynaComplex(0,0);
            }
            initialState[n].set(1,0);
        }
        
        // Create the integrator
        MCWFThreadPoolIntegrator integrator = new MCWFThreadPoolIntegrator(numTrajectories, numTimes, dt, evDelta, params, initialState, numThreads, debug);
        
        if(cmd.hasOption("icf")) {
            // Override initial conditions from a file
            setICsFromFile(integrator.getTrajectories(), cmd.getOptionValue("icf"));
        }
        
        // Integrate
        long startTime = System.nanoTime();
        integrator.start();
        boolean success = integrator.waitForFinished(timeout);

        long endTime = System.nanoTime();
        System.out.println("Run time: " + (endTime - startTime)/1.0e9 + " seconds");
        
        if(!success) {
            System.err.println("Received failure from integrator");
            return;
        }
        
        // Perform validation
        if(debug) {
            QuantumTrajectory[] trajectories = integrator.getTrajectories();
            ArrayList<String> msgs = new ArrayList<String>();
            int numFailed = 0;
            for(int i = 0; i < trajectories.length; ++i) {
                msgs.clear();
                if(!trajectories[i].getStats().validate(msgs)) {
                    if(numFailed < MAX_VALIDATION_MSGS) {
                        System.err.println("Trajectory " + i + " failed validation. Messages:");
                        for(int j = 0; j < msgs.size(); ++j) {
                            System.err.println(msgs.get(j));
                        }
                    }
                    ++numFailed;

                    if(numFailed == MAX_VALIDATION_MSGS) {
                        System.err.println("Supressing further validation error output");
                    }
                }
            }

            if(numFailed == 0) {
                System.out.println("Validation OK");
            } else {
                System.err.println(numFailed + " trajectories failed validation");
            }
        }
        
        // Write results
        String uuidStr = UUID.randomUUID().toString();
        String outfile = outdir + "/mcwf" + String.format("_traj1-%d_", numTrajectories) + params.getMcwfFilename() + "_" + uuidStr + ".txt";
        MCWFWriter writer = new MCWFWriter();
        writer.write(integrator.getAggregator(), outfile);
        
        // Write debug info
        if(debug) {
            String dbgdir = outdir + "/debug";
            File dfdir = new File(dbgdir);
            dfdir.mkdirs();
            
            String jumpfile = dbgdir + "/jumps_" + uuidStr + ".txt";
            writer.writeJumps(integrator.getAggregator(), jumpfile);
            
            String statefile = dbgdir + "/state_" + uuidStr + ".txt";
            writer.writeState(integrator.getAggregator(), statefile);
            
            writer.writeHusimi(integrator.getAggregator(), dbgdir, uuidStr);
        }
        
        // Write final states if requested
        if(cmd.hasOption("wf")) {
            String finalfile;
            if(cmd.getOptionValue("wf").equals("def")) {
                finalfile = outdir + "/final_" + uuidStr + ".txt";
            } else {
                finalfile = cmd.getOptionValue("wf");
            }
            
            writer.writeAllStates(integrator.getTrajectories(), finalfile);
        }
    }

    private static double getAjmp(int j, int m) {
        return 0.5*Math.sqrt((double)(j-m)*(double)(j+m+1));
    }
    
    private static void setICsFromFile(QuantumTrajectory[] trajectories, String filename) throws FileNotFoundException {
        DynaComplex[] newState = new DynaComplex[trajectories[0].getState().length];
        for(int i = 0; i < newState.length; ++i) {
            newState[i] = new DynaComplex(0,0);
        }
        
        int trajIdx = 0;
        Scanner fileReader = new Scanner(new File(filename));
        while(fileReader.hasNextLine()) {
            QuantumTrajectory traj = trajectories[trajIdx];
            
            String line = fileReader.nextLine();
            String[] lineVec = line.split(",");
            if(lineVec.length != traj.getState().length*2 + 1) {
                fileReader.close();
                throw new UnsupportedOperationException("Initial state " + trajIdx + " contained incorrect number of elements " + lineVec.length + ". Expected: " + traj.getState().length*2 + 1);
            }

            int jval = Integer.parseInt(lineVec[0]);
            for(int i = 0; i < newState.length; ++i) {
                newState[i].set(Double.parseDouble(lineVec[2*i + 1]), Double.parseDouble(lineVec[2*i + 2]));
            }
            traj.setState(newState, jval);
            
            ++trajIdx;
        }
        fileReader.close();
        
        if(trajIdx != trajectories.length) {
            throw new UnsupportedOperationException("Initial condition file had incorrect number of trajectories");
        }
    }
}
