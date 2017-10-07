package exe;

import handlers.SummaryWriter;
import handlers.WriteHandlerCorr;
import handlers.WriteHandlerMaster;
import handlers.WriteHandlerMeanField;

import java.io.FileNotFoundException;
import java.io.UnsupportedEncodingException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.TreeMap;
import java.util.TreeSet;

import ode.CumulantAllToAllODEs;
import ode.CumulantParams;
import ode.DynaComplexODEAdapter;
import ode.MasterAllToAllODEs;
import ode.SynchMeanFieldODEs;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;
import org.apache.commons.math3.ode.nonstiff.AdamsMoultonIntegrator;
import org.apache.commons.math3.ode.sampling.StepHandler;

import utils.DynaComplex;
import utils.SynchUtils;
import eval.CumulantEval;
import eval.MasterEval;
import eval.MeanFieldEval;
import eval.SystemEval;

public class RunSim {
    
    public enum Simulator {
        MEAN_FIELD,
        CUMULANT,
        MASTER
    }

    private static Options setupOptions() {
        Options options = new Options();
        options.addOption("n", true, "Number of particles");
        options.addOption("w", true, "Incoherent pumping");
        options.addOption("d", true, "Disorder sigma");
        options.addOption("f", true, "Elastic interaction term");
        options.addOption("g", true, "Inelastic interaction term");
        options.addOption("t", false, "Output results versus time");
        options.addOption("h", false, "Print this help message");
        
        return options;
    }
    
    private static void printUsage(Options options) {
        System.out.println("RunSim [options] SIMULATOR outdir\n");
        System.out.println("Simulator types:");
        for(int i = 0; i < Simulator.values().length; ++i) {
            System.out.println(Simulator.values()[i]);
        }
        
        System.out.println("\nAll numeric options can be single values or ranges specified as min:inc:max");
        System.out.println("\nOptions:");
        
        Iterator<Option> iter = options.getOptions().iterator();
        while(iter.hasNext()) {
            Option opt = iter.next();
            System.out.println("    " + opt.getOpt() + ": " + opt.getDescription());
        }
    }
    
    public static void main(String[] args) throws ParseException, FileNotFoundException, UnsupportedEncodingException {
        // Defaults
        Simulator sim = Simulator.MEAN_FIELD;
        String outdir = "/Users/kristophertucker/output/mf/";
        int n = 2;
        double h = 0.001;
        double gamma = 1.0;
        double f = 1.0;
        double g = 30.0;
        double delta = 0.0;
        double w = 3.0;
        
        // TODO - Steady state detection and more general tmax
        double tmax = 5;
        
        // TODO - Figure out a more general way to handle detunings
        double[] d = new double[n];
        generateDetunings(delta, d);
        
        CumulantParams dparams = new CumulantParams(n, gamma, w, delta, new DynaComplex(f,g), d);
        
        // Parse the command line
        Options options = setupOptions();
        
        CommandLineParser parser = new DefaultParser();
        CommandLine cmd = parser.parse(options, args);
        
        if(cmd.hasOption("h")) {
            printUsage(options);
            return;
        }
        
        boolean outputVsTime = cmd.hasOption("t");
        
        TreeMap<String, double[]> optMap = new TreeMap<String, double[]>();
        if(!parseOptions(options, args, cmd, optMap)) {
            printUsage(options);
            return;
        }
        
        // Get the simulator type
        if(cmd.getArgList().size() > 0) {
            try {
                sim = Simulator.valueOf(cmd.getArgList().get(0));
            } catch(IllegalArgumentException ex) {
                System.err.println("Invalid simulator type.  Please choose from the following:");
                for(int i = 0; i < Simulator.values().length; ++i) {
                    System.err.println(Simulator.values()[i].name());
                }
                return;
            }
        }
        
        // Get the output directory
        if(cmd.getArgList().size() > 1) {
            outdir = cmd.getArgList().get(1);
        }
        
        // Create params list
        ArrayList<CumulantParams> params = new ArrayList<CumulantParams>();
        createParams(optMap, dparams, params);
        
        // Initialize the summary writer
        SummaryWriter summWriter = new SummaryWriter(outdir + "/summary.txt");
        
        long startTime = System.nanoTime();
        
        for(int i = 0; i < params.size(); ++i) {
            System.out.println(params.get(i).toString());

            n = params.get(i).getN();
            
            // Initialize the ODE object
            StepHandler writeHandler = null;
            SystemEval eval;
            FirstOrderDifferentialEquations odes;
            switch(sim) {
            case MEAN_FIELD:
                if(outputVsTime) {
                    writeHandler = new WriteHandlerMeanField(outdir + "/mf_" + params.get(i).getFilename(), n);
                }
                eval = new MeanFieldEval(n);
                odes = new SynchMeanFieldODEs(params.get(i));
                break;
            case CUMULANT:
                if(outputVsTime) {
                    writeHandler = new WriteHandlerCorr(outdir + "/cumulant_" + params.get(i).getFilename(), n);
                }
                eval = new CumulantEval(n);
                CumulantAllToAllODEs codes = new CumulantAllToAllODEs(params.get(i));
                odes = new DynaComplexODEAdapter(codes);
                break;
            case MASTER:
                if(outputVsTime) {
                    writeHandler = new WriteHandlerMaster(outdir + "/master_" + params.get(i).getFilename(), n);
                }
                eval = new MasterEval(n);
                MasterAllToAllODEs modes = new MasterAllToAllODEs(params.get(i));
                odes = new DynaComplexODEAdapter(modes);
                break;
            default:
                System.err.println("Simulator type not yet supported");
                return;
            }
            
            // Setup initial conditions
            double[] y0 = new double[eval.getRealDimension()];
            eval.initSpinUpX(y0);

            // Perform the integration
            double[] y = new double[eval.getRealDimension()];
            AdamsMoultonIntegrator integrator = new AdamsMoultonIntegrator(2, h*1.0e-4, h, 1.0e-3, 1.0e-2);
            
            if(writeHandler != null) {
                integrator.addStepHandler(writeHandler);
            }
            
            integrator.integrate(odes, 0, y0, tmax, y);
            
            // Add the results to the summary writer
            summWriter.addVals(eval, params.get(i), y);
        }
        
        long endTime = System.nanoTime();
        
        System.out.println("Run time: " + (endTime - startTime)/1.0e9 + " seconds");
        
        // Write the summary
        summWriter.writeToFile();
    }
    
    private static boolean parseOptions(Options options, String[] args, CommandLine cmd, TreeMap<String, double[]> optMap) {
        TreeSet<String> ignore = new TreeSet<String>();
        ignore.add("t");
        
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
    
    private static void createParams(TreeMap<String, double[]> optMap, CumulantParams dparams, ArrayList<CumulantParams> params) {
        double[] nvals = getValList(optMap, "n", dparams.getN());
        double[] wvals = getValList(optMap, "w", dparams.getW());
        double[] dvals = getValList(optMap, "d", dparams.getDelta());
        double[] fvals = getValList(optMap, "f", dparams.getAlpha().getReal());
        double[] gvals = getValList(optMap, "g", dparams.getAlpha().getImaginary());
        
        for(int nidx = 0; nidx < nvals.length; ++nidx) {
            for(int widx = 0; widx < wvals.length; ++widx) {
                for(int didx = 0; didx < dvals.length; ++didx) {
                    for(int fidx = 0; fidx < fvals.length; ++fidx) {
                        for(int gidx = 0; gidx < gvals.length; ++gidx) {
                            int n = (int)nvals[nidx];
                            double[] d = new double[n];
                            generateDetunings(dvals[didx], d);
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
    
    private static void generateDetunings(double delta, double[] d) {
        SynchUtils.detuneGauss(delta, d);
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
}
