package exe;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.Scanner;
import java.util.TreeMap;
import java.util.TreeSet;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;
import org.apache.commons.math3.ode.nonstiff.AdamsMoultonIntegrator;
import org.apache.commons.math3.ode.sampling.StepHandler;

import eval.CumulantEval;
import eval.MasterEval;
import eval.MeanFieldEval;
import eval.RPAEval;
import eval.SystemEval;
import eval.SystemEval.InitAngleType;
import handlers.CumulantSteadyStateTerminator;
import handlers.DataRecorder;
import handlers.SummaryWriter;
import handlers.WriteHandlerCorr;
import handlers.WriteHandlerMaster;
import handlers.WriteHandlerMeanField;
import handlers.WriteHandlerRPA;
import integrator.AdamsMoultonFactory;
import integrator.IntegratorRequest;
import integrator.ThreadPoolIntegrator;
import ode.CumulantAllToAllODEs;
import ode.CumulantParams;
import ode.DynaComplexODEAdapter;
import ode.MasterAllToAllODEs;
import ode.RPAAllToAllODEs;
import ode.SynchMeanFieldODEs;
import utils.DynaComplex;
import utils.SynchReducer;
import utils.SynchSolution;
import utils.SynchUtils;

public class RunSim {
    
    public enum Simulator {
        MEAN_FIELD,
        CUMULANT,
        MASTER,
        RPA
    }
    
    public static class ICAngleParams {
        public double initAngle;
        public InitAngleType type;

        public ICAngleParams() {
            initAngle = 0.0;
            type = InitAngleType.CONST;
        }
    }

    private static Options setupOptions(TreeSet<String> ignore) {
        Options options = new Options();
        options.addOption("n", true, "Number of particles");
        options.addOption("w", true, "Incoherent pumping");
        options.addOption("d", true, "Disorder sigma");
        options.addOption("f", true, "Elastic interaction term");
        options.addOption("g", true, "Inelastic interaction term");
        options.addOption("t", false, "Output results versus time");
        options.addOption("zd", false, "Output collective z distribution when running MASTER");
        options.addOption("df", true, "File of detunings to use");
        options.addOption("l", false, "Use Lorentzian detunings (default is Gaussian)");
        options.addOption("tmin", true, "Minimum simulation time");
        options.addOption("tmax", true, "Maximum simulation time");
        options.addOption("iz", true, "Initial zenith angle (CONST value in degrees | EQUAL_SPACING | RANDOM)");
        options.addOption("ip", true, "Initial phase angle (CONST value in degrees | EQUAL_SPACING | RANDOM)");
        options.addOption("wic", true, "Write initial condition Bloch vectors to the given file");
        options.addOption("c", false, "Carry over initial conditions from the previous solution");
        options.addOption("tt", false, "Perform two-time correlation simulation");
        options.addOption("nt", true, "Number of threads to use for integration");
        options.addOption("h", false, "Print this help message");
        
        // Options to ignore when setting up run parameters
        ignore.add("t");
        ignore.add("zd");
        ignore.add("df");
        ignore.add("l");
        ignore.add("tmin");
        ignore.add("tmax");
        ignore.add("iz");
        ignore.add("ip");
        ignore.add("wic");
        ignore.add("c");
        ignore.add("tt");
        ignore.add("nt");
        
        return options;
    }
    
    private static void printUsage(Options options) {
        System.out.println("RunSim [options] outdir SIMULATOR\n");
        System.out.println("Simulator types:");
        for(int i = 0; i < Simulator.values().length; ++i) {
            System.out.println(Simulator.values()[i]);
        }
        
        System.out.println("\nAll numeric options can be single values or ranges specified as min:inc:max");
        System.out.println("\nOptions:");
        
        Iterator<Option> iter = options.getOptions().iterator();
        while(iter.hasNext()) {
            Option opt = iter.next();
            System.out.printf("%5s: %s\n", opt.getOpt(), opt.getDescription());
        }
    }
    
    public static void main(String[] args) throws ParseException, FileNotFoundException, UnsupportedEncodingException {
        // Defaults
        Simulator sim = Simulator.MEAN_FIELD;
        String outdir = "/Users/kristophertucker/output/temp/";
        int n = 2;
        double h = 0.001;
        double gamma = 1.0;
        double f = 0.0;
        double g = 0.0;
        double delta = 10.0;
        double w = 0.0;
        
        double tmin = 2;
        double tmax = 20;
        
        long timeout = 7*24*3600; // One week
        
        // Initial condition defaults
        ICAngleParams izParams = new ICAngleParams();
        izParams.initAngle = Math.PI/2.0;
        izParams.type = InitAngleType.CONST;
        
        ICAngleParams ipParams = new ICAngleParams();
        ipParams.initAngle = 0.0;
        ipParams.type = InitAngleType.EQUAL_SPACING;
        
        // Parse the command line
        TreeSet<String> optIgnore = new TreeSet<String>();
        Options options = setupOptions(optIgnore);
        
        CommandLineParser parser = new DefaultParser();
        CommandLine cmd = parser.parse(options, args);
        
        if(cmd.hasOption("h")) {
            printUsage(options);
            return;
        }
        
        boolean outputVsTime = cmd.hasOption("t");
        boolean outputZdist = cmd.hasOption("zd");
        boolean useLorDetunings = cmd.hasOption("l");
        boolean carryOverIC = cmd.hasOption("c");
        boolean correlate = cmd.hasOption("tt");
        
        int numThreads = 1;
        if(cmd.hasOption("nt")) {
            numThreads = Integer.parseInt(cmd.getOptionValue("nt"));
        }
        
        // Get initial conditions output file path if provided
        boolean outputIC = false;
        String icOutFile = "";
        if(cmd.hasOption("wic")) {
            outputIC = true;
            icOutFile = cmd.getOptionValue("wic");
        }
        
        // Process initial condition overrides
        if(cmd.hasOption("iz")) {
            if(!parseInitAngle(cmd.getOptionValue("iz"), "iz", izParams)) {
                return;
            }
        }
        
        if(cmd.hasOption("ip")) {
            if(!parseInitAngle(cmd.getOptionValue("ip"), "ip", ipParams)) {
                return;
            }
        }
        
        System.out.println("IC - zenith: " + izParams.type.name() + ", " + izParams.initAngle + 
                           "; phase: " + ipParams.type.name() + ", " + ipParams.initAngle);

        // Override tmin/tmax
        if(cmd.hasOption("tmin")) {
            try {
                tmin = Double.parseDouble(cmd.getOptionValue("tmin"));
            } catch(NumberFormatException ex) {
                System.err.println("Failed to parse tmin value " + cmd.getOptionValue("tmin"));
                return;
            }
        }
        
        if(cmd.hasOption("tmax")) {
            try {
                tmax = Double.parseDouble(cmd.getOptionValue("tmax"));
            } catch(NumberFormatException ex) {
                System.err.println("Failed to parse tmax value " + cmd.getOptionValue("tmax"));
                return;
            }
        }
        
        // Read detunings from a file if requested
        boolean detuningsFromFile = cmd.hasOption("df");
        ArrayList<Double> fileDetunings = null;
        if(detuningsFromFile) {
            fileDetunings = new ArrayList<Double>();
            Scanner inputStream = new Scanner(new File(cmd.getOptionValue("df")));
            while(inputStream.hasNext()) {
                fileDetunings.add(Double.parseDouble(inputStream.next()));
            }
            inputStream.close();
        }

        // Default detunings
        double[] d = new double[n];
        generateDetunings(delta, d, useLorDetunings);
        
        CumulantParams dparams = new CumulantParams(n, gamma, w, delta, new DynaComplex(f,g), d);
        
        TreeMap<String, double[]> optMap = new TreeMap<String, double[]>();
        if(!parseOptions(options, optIgnore, args, cmd, optMap)) {
            printUsage(options);
            return;
        }
        
        // Get the output directory
        if(cmd.getArgList().size() > 0) {
            outdir = cmd.getArgList().get(0);
        }
        File fdir = new File(outdir);
        fdir.mkdirs();
        
        // Get the simulator type
        if(cmd.getArgList().size() > 1) {
            try {
                sim = Simulator.valueOf(cmd.getArgList().get(1));
            } catch(IllegalArgumentException ex) {
                System.err.println("Invalid simulator type.  Please choose from the following:");
                for(int i = 0; i < Simulator.values().length; ++i) {
                    System.err.println(Simulator.values()[i].name());
                }
                return;
            }
        }
        
        // Create params list
        ArrayList<CumulantParams> params = new ArrayList<CumulantParams>();
        createParams(optMap, dparams, params, useLorDetunings);
        
        // Initialize the summary writer
        SummaryWriter summWriter = null;
        if(params.size() > 1 || !outputVsTime) {
            summWriter = new SummaryWriter(outdir + "/" + sim.name().toLowerCase() + "_summary.txt");
        }
        
        // Prepare the thread pool integrator
        ThreadPoolIntegrator tpIntegrator = null;
        ArrayList<SynchSolution> solns = new ArrayList<SynchSolution>();
        if(numThreads > 1) {
            SynchReducer reducer = null;
            if(summWriter != null) {
                reducer = new SynchReducer(summWriter);
                // TODO - Set live update to false on summWriter
            }
            
            tpIntegrator = new ThreadPoolIntegrator(numThreads, new AdamsMoultonFactory(h*1.0e-6, h), reducer);
        }
        
        long startTime = System.nanoTime();
        
        int lastn = -1;
        double[] y0 = null;
        for(int i = 0; i < params.size(); ++i) {
            n = params.get(i).getN();
            
            // Override the detunings if provided in a file
            if(detuningsFromFile) {
                if(n != fileDetunings.size()) {
                    System.err.println("Number of detunings in provided file does not match number of particles");
                    break;
                }
                
                for(int j = 0; j < fileDetunings.size(); ++j) {
                    params.get(i).getD()[j] = fileDetunings.get(j);
                }
                
//                String filename = cmd.getOptionValue("df") + "N" + n + "_D" + (int)(params.get(i).getDelta()) + ".txt";
//                SynchUtils.detuneFile(filename, params.get(i).getD());
            }
            
            if(numThreads == 1) {
                System.out.println(params.get(i).toString());
            }

            String timefile = outdir + "/" + sim.name().toLowerCase() + "_" + params.get(i).getFilename();
            
            // Initialize the ODE object
            StepHandler writeHandler = null;
            CumulantSteadyStateTerminator term = null;
            SystemEval eval;
            FirstOrderDifferentialEquations odes;
            switch(sim) {
            case MEAN_FIELD:
                if(outputVsTime) {
                    writeHandler = new WriteHandlerMeanField(timefile, n);
                }
                eval = new MeanFieldEval(n);
                odes = new SynchMeanFieldODEs(params.get(i));
                break;
            case CUMULANT:
                if(outputVsTime) {
                    writeHandler = new WriteHandlerCorr(timefile, n);
                }
                term = new CumulantSteadyStateTerminator(tmin, 0.015, 50, 1000000, 0.0025, n);
                term.setQuietMode(numThreads > 1);
                eval = new CumulantEval(n);
                CumulantAllToAllODEs codes = new CumulantAllToAllODEs(params.get(i));
                odes = new DynaComplexODEAdapter(codes);
                break;
            case MASTER:
                if(outputVsTime) {
                    writeHandler = new WriteHandlerMaster(timefile, n);
                }
                eval = new MasterEval(n);
                MasterAllToAllODEs modes = new MasterAllToAllODEs(params.get(i));
                odes = new DynaComplexODEAdapter(modes);
                break;
            case RPA:
                if(outputVsTime) {
                    writeHandler = new WriteHandlerRPA(timefile, n);
                }
                eval = new RPAEval(n);
                odes = new RPAAllToAllODEs(params.get(i));
                break;
            default:
                System.err.println("Simulator type not yet supported");
                return;
            }
            
            // Setup the data recorder
            DataRecorder recorder = new DataRecorder(eval, 0.01, 1.0);
            
            // Setup initial conditions
            if(n != lastn) {
                lastn = n;
                y0 = new double[eval.getRealDimension()];
                if(sim == Simulator.MASTER) {
                    // TODO - Remove this once MasterEval has implemented the initialize method
                    eval.initSpinUpX(y0);
                } else {
                    eval.initialize(y0, izParams.initAngle, ipParams.initAngle, izParams.type, ipParams.type);
                }
                
                if(outputIC) {
                    if(i > 0) {
                        System.err.println("ERROR: Multiple initial conditions used with wic flag");
                        break;
                    }
                    
                    outputInitialBVs(eval, y0, icOutFile);
                }
            }

            // Setup integrator
            double[] y = null;
            if(numThreads == 1) {
                // We are running single threaded, so simply integrate now
                AdamsMoultonIntegrator integrator = new AdamsMoultonIntegrator(2, h*1.0e-6, h, 1.0e-3, 1.0e-2);
                integrator.addStepHandler(recorder);

                if(writeHandler != null) {
                    integrator.addStepHandler(writeHandler);
                }

                if(term != null) {
                    integrator.addStepHandler(term.getDetector());
                    integrator.addEventHandler(term, Double.POSITIVE_INFINITY, 1.0e-12, 100);
                }

                // Perform the integration
                y = new double[eval.getRealDimension()];
                integrator.integrate(odes, 0, y0, tmax, y);
                
                SynchSolution soln = new SynchSolution(params.get(i), recorder, eval);
                soln.setSolution(y);
                solns.add(soln);
            } else {
                // We are running multi-threaded, so queue up the request
                SynchSolution soln = new SynchSolution(params.get(i), recorder, eval);
                solns.add(soln);
                
                IntegratorRequest request = new IntegratorRequest(odes, 0, y0, tmax, soln);
                
                request.addStepHandler(recorder);

                if(writeHandler != null) {
                    request.addStepHandler(writeHandler);
                }

                if(term != null) {
                    request.addStepHandler(term.getDetector());
                    request.addEventHandler(term);
                }
                
                tpIntegrator.addIvp(request);
            }
            
            // Carry over the initial conditions if requested
            if(carryOverIC) {
                if(numThreads != 1) {
                    throw new UnsupportedOperationException("Initial condition carry over not supported in multi-threaded mode");
                }
                
                for(int j = 0; j < y.length; ++j) {
                    y0[j] = y[j];
                    if(Math.abs(y0[j]) < 1.0e-16) {
                        y0[j] = 0.0;
                    }
                    y0[j] += 1.0e-10;
                }
            }
            
            // Add the results to the summary writer
            if(numThreads == 1 && summWriter != null) {
                summWriter.addVals(params.get(i), recorder);
            }
        }
        
        if(numThreads > 1) {
            tpIntegrator.waitForFinished(timeout);
        }
        
        long simEndTime = System.nanoTime();
        
        System.out.println("Simulation time: " + (simEndTime - startTime)/1.0e9 + " seconds");
        
        // Write the summary
        // TODO - Only do this for multiple threads if single threaded mode is updating live
        if(summWriter != null) {
            summWriter.writeToFile();
        }
        
        // Post processing
        for(int i = 0; i < solns.size(); ++i) {
            SynchSolution soln = solns.get(i);
            
            if(outputZdist && sim == Simulator.MASTER) {
                ((MasterEval)soln.getEval()).writeZDist(soln.getSolution(), outdir + "zdist_" + params.get(i).getFilename());
            }
            
//            if(sim == Simulator.MASTER) {
//                ((MasterEval)eval).writeTriples(y);
//            }
            
            // Compute the correlation function if requested
            if(correlate) {
                if(sim == Simulator.CUMULANT) {
                    SynchUtils.compCorr(params.get(i), soln.getSolution(), outdir + "time_corr_" + soln.getParams().getFilename());
                } else {
                    throw new UnsupportedOperationException("Two-time correlation currently only supported for CUMULANT simulator");
                }
            }
        }
        
        long endTime = System.nanoTime();
        
        System.out.println("Run time: " + (endTime - startTime)/1.0e9 + " seconds");
    }
    
    private static boolean parseOptions(Options options, TreeSet<String> ignore, String[] args, CommandLine cmd, TreeMap<String, double[]> optMap) {
        Iterator<Option> iter = options.getOptions().iterator();
        while(iter.hasNext()) {
            Option opt = iter.next();
            
            if(ignore.contains(opt.getOpt())) continue;

            if(cmd.hasOption(opt.getOpt())) {
                double[] dvals = parseOption(cmd.getOptionValue(opt.getOpt()), opt.getOpt());
                if(dvals == null) {
                    return false;
                }

                optMap.put(opt.getOpt(), dvals);
                
//                System.out.println(opt.getOpt() + " vals:");
//                for(int i = 0; i < dvals.length; ++i) {
//                    System.out.println(dvals[i]);
//                }
            }
        }
        
        return true;
    }
    
    private static void createParams(TreeMap<String, double[]> optMap, CumulantParams dparams, ArrayList<CumulantParams> params, boolean useLorDetunings) {
        double[] nvals = getValList(optMap, "n", dparams.getN());
        double[] wvals = getValList(optMap, "w", dparams.getW());
        double[] dvals = getValList(optMap, "d", dparams.getDelta());
        double[] fvals = getValList(optMap, "f", dparams.getAlpha().getReal());
        double[] gvals = getValList(optMap, "g", dparams.getAlpha().getImaginary());
        
        for(int nidx = 0; nidx < nvals.length; ++nidx) {
            int n = (int)nvals[nidx];
            for(int widx = 0; widx < wvals.length; ++widx) {
                for(int didx = 0; didx < dvals.length; ++didx) {
                    double[] d = new double[n];
                    generateDetunings(dvals[didx], d, useLorDetunings);
                    for(int fidx = 0; fidx < fvals.length; ++fidx) {
                        for(int gidx = 0; gidx < gvals.length; ++gidx) {
                            params.add(new CumulantParams(n, dparams.getGamma(),
                                                          wvals[widx], dvals[didx],
                                                          new DynaComplex(fvals[fidx], gvals[gidx]), d));
                        }
                    }
                }
            }
        }
    }
    
    private static double[] getValList(TreeMap<String, double[]> optMap, String key, double dval) {
        if(optMap.containsKey(key)) {
            return optMap.get(key);
        }
        
        double[] vals = new double[1];
        vals[0] = dval;
        
        return vals;
    }
    
    private static void generateDetunings(double delta, double[] d, boolean useLorDetunings) {
        if(useLorDetunings) {
            SynchUtils.detuneLor(delta, d);
        } else {
            SynchUtils.detuneGauss(delta, d);
        }
    }
    
    private static double[] parseOption(String opt, String name) {
        String[] optvec = opt.split(":");
        
        double[] vals = null;
        try {
            if(optvec.length == 1) {
                vals = new double[1];
                vals[0] = Double.parseDouble(optvec[0]);
            } else if(optvec.length == 3) {
                double[] incvals = new double[3];
                for(int i = 0; i < incvals.length; ++i) {
                    incvals[i] = Double.parseDouble(optvec[i]);
                }
                
                int num = (int)((incvals[2] - incvals[0])/incvals[1]) + 1;
                vals = new double[num];
                for(int i = 0; i < num; ++i) {
                    vals[i] = incvals[0] + i*incvals[1];
                }
            } else {
                System.err.println("Option " + name + " incorrectly formatted");
            }
        } catch(NumberFormatException ex) {
            System.err.println("Failed to parse option " + name);
            return null;
        }
        
        return vals;
    }
    
    private static boolean parseInitAngle(String optVal, String optName, ICAngleParams params) {
        try {
            params.initAngle = Double.parseDouble(optVal);
            params.initAngle *= Math.PI/180.0;
            params.type = InitAngleType.CONST;
        } catch(NumberFormatException ex) {
            params.initAngle = 0.0;
            try {
                params.type = InitAngleType.valueOf(optVal);
                if(params.type == InitAngleType.CONST) {
                    System.err.println("Must specify a const value in degrees for option " + optName);
                    return false;
                }
            } catch(IllegalArgumentException ex2) {
                System.err.println("Invalid value of option " + optName + ".  Choices are:");
                for(int i = 0; i < InitAngleType.values().length; ++i) {
                    if(InitAngleType.values()[i].equals(InitAngleType.valueOf("CONST"))) {
                        System.err.println("Initial const value (in degrees)");
                    } else {
                        System.err.println(InitAngleType.values()[i]);
                    }
                }
                
                return false;
            }
        }
        
        return true;
    }
    
    private static void outputInitialBVs(SystemEval eval, double[] y0, String filepath) throws FileNotFoundException {
        int n = eval.getN();
        double[] xs = new double[n];
        double[] ys = new double[n];
        double[] zs = new double[n];
        
        eval.getBlochVectors(y0, xs, ys, zs);
        
        PrintWriter writer = new PrintWriter(filepath);
        for(int i = 0; i < n; ++i) {
            writer.write(xs[i] + ", " + ys[i] + ", " + zs[i] + "\n");
        }
        writer.close();
    }
}
