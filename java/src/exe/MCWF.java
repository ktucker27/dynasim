package exe;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.UnsupportedEncodingException;
import java.util.Iterator;
import java.util.UUID;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.MissingOptionException;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

import mcwf.MCWFIntegrator;
import mcwf.MCWFWriter;
import ode.SystemParams;
import utils.DynaComplex;

public class MCWF {

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
        
        long timeout = 7*24*3600; // One week
        
        double[] d = new double[n];
        for(int i = 0; i < n; ++i) {
            d[i] = 0.0;
        }
        
        SystemParams params = new SystemParams(n, gamma, w, o, gel, 0, faa, fab, gaa, gab, d);
        
        // Define the initial condition
        DynaComplex[] initialState = new DynaComplex[n+1];
        for(int i = 0; i < n+1; ++i) {
            initialState[i] = new DynaComplex(0,0);
        }
        initialState[n].set(1,0);
        
        // Create the integrator
        MCWFIntegrator integrator = new MCWFIntegrator(numTrajectories, numTimes, dt, params, initialState, numThreads);
        
        // Integrate
        long startTime = System.nanoTime();
        integrator.start();
        integrator.waitForFinished(timeout);

        long endTime = System.nanoTime();
        System.out.println("Run time: " + (endTime - startTime)/1.0e9 + " seconds");
        
        // Write results
        String outfile = outdir + "/mcwf" + String.format("_traj1-%d_", numTrajectories) + params.getMcwfFilename() + "_" + UUID.randomUUID().toString() + ".txt";
        MCWFWriter writer = new MCWFWriter(outfile);
        writer.write(integrator.getAggregator());
    }

}
