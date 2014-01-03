////////////////////////////////////////////////////////////////////////////////////////////////////
// main.cc
// Program to use a spring system to converge on a fit for an idealized 2d model
// Note that this code, and all code included in this directory, was written with the Schira model
// (from Schira, Tyler, Spehar, and Breakspear, 2010, DOI:10.1371/journal.pcbi.1000651) in mind.
// by Noah Benson

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <sys/stat.h>

#include "springs.hh"
#include "schira.hh"
#include "params.hh"
#include "vdw.hh"
#include "angle_springs.hh"
#include "coord_spring.hh"

using namespace std;
using namespace springs;



////////////////////////////////////////////////////////////////////////////////////////////////////
// Before main, we want to declare all the parameter types that can be declared

class param_positive_real: public param<double>
{public:
   param_positive_real(const string& name, char abbrev, const string& doc, const string& descr, 
                       const double& default_value):
      param<double>(name, abbrev, doc, descr, default_value)
   {
   }
 protected:
   void process()
   {
      if (m_value <= 0)
         throw (string("Parameter") + name() + " must be a positive real number");
   }
};
class param_nonnegative_real: public param<double>
{public:
   param_nonnegative_real(const string& name, char abbrev, const string& doc, const string& descr, 
                          const double& default_value):
      param<double>(name, abbrev, doc, descr, default_value)
   {
   }
 protected:
   void process()
   {
      if (m_value < 0)
         throw (string("Parameter") + name() + " must be a non-negative real number");
   }
};
class param_probability: public param<double>
{public:
   param_probability(const string& name, char abbrev, const string& doc, const string& descr, 
                       const double& default_value):
      param<double>(name, abbrev, doc, descr, default_value)
   {
   }
 protected:
   void process()
   {
      if (m_value < 0 || m_value > 1)
         throw (string("Parameter") + name() + " must be a real number between 0 and 1");
   }
};
class param_positive_integer: public param<int>
{public:
   param_positive_integer(const string& name, char abbrev, const string& doc, const string& descr, 
                          const int& default_value):
      param<int>(name, abbrev, doc, descr, default_value)
   {
   }
 protected:
   void process()
   {
      if (m_value <= 0)
         throw (string("Parameter") + name() + " must be a positive integer");
   }
};
class param_nonnegative_integer: public param<int>
{public:
   param_nonnegative_integer(const string& name, char abbrev, const string& doc, const string& descr, 
                             const int& default_value):
      param<int>(name, abbrev, doc, descr, default_value)
   {
   }
 protected:
   void process()
   {
      if (m_value < 0)
         throw (string("Parameter") + name() + " must be a non-negative integer");
   }
};
class param_journal: public param<string>
{public:
   param_journal(const string& name, char abbrev, const string& doc, const string& descr, 
                          const string& default_value):
      param<string>(name, abbrev, doc, descr, default_value)
   {
   }
 protected:
   void post_process()
   {
      struct stat st;
      int r;
      if (!provided()) return;
      r = stat(m_value.data(), &st);
      if (r != 0) {
         // we need to create it...
         if (mkdir(m_value.data(), 0755) != 0) {
            ostringstream ss;
            ss << "Could not make directory " << m_value
               << " for journaling!";
            throw ss.str();
         }
         r = stat(m_value.data(), &st);
         if (r != 0) {
            ostringstream ss;
            ss << "Created journal directory " << m_value 
               << " but could not stat it!";
            throw ss.str();
         }
      }
      // okay, we should have stat structure now; make sure we can read/write it
      if (!S_ISDIR(st.st_mode)) {
         ostringstream ss;
         ss << "Journal directory " << m_value
            << " is not a directory!";
         throw ss.str();
      }
   }
};



////////////////////////////////////////////////////////////////////////////////////////////////////
// The Main Function!

int main(int argc, char** argv)
{
   int i, j;

   /////////////////////////////////////////////////////////////////////////////////////////////////
   // Argument parsing...

   param_positive_real 
      timestep("timestep", 't',
               "The duration of each iteration of the numerical integration loop",
               "The simulation's timestep indicates how long each loop of the simulation ostensibly"
               " lasts.  A large timestep will cause the simulation to run more quickly (i.e.,"
               " fewer loops will be required for the same amount of time to pass), but a timestep"
               " that is too large will cause numerical drift in the energies.",
               0.002);
   param_positive_integer
      steps("steps", 's',
            "The number of steps to run during each simulated annealing iteration",
            "This parameter sets the number of steps to be run during each annealing iteration."
            " The amount of time passed during each annealing iteration is timesteps * steps.", 
            5000);
   param_probability
      dampen("dampen", 'd',
             "The dampening coefficient for the simulation",
             "Dampening is the extend to which the simulation is cooled every integration step.  A"
             " dampening coefficient of 1 indicates that no dampening occurs (v(t) = d*v(t-1)); a"
             " dampening coefficient near 0 will cause the simulation to act as a minimizer rather"
             " than a proper simulation.",
             1.0);
   param_positive_real
      spring_str("spring-strength", 'k',
                 "The strength constant for anatomical springs",
                 "The spring-strength sets the const k in the spring equation f = k * d; setting"
                 " this parameter to a high value will cause the anatomical springs to be very"
                 " strong while setting it to a low value will cause them to be weak.",
                 1.0);
   param_positive_real
      schira_str("schira-strength", 'b',
                 "The strength constant for the Schira model springs",
                 "The schira-strength sets the const k in the spring equation f = k * d for the"
                 " Schira model potential function.  Setting this parameter to a high value will"
                 " cause the Schira model to be very strong while setting it to a low value will"
                 " cause it to be weak.  Generally speaking, it is a good idea to keep this value"
                 " lower than the spring-strength.",
                 0.5);
   param_positive_real
      schira_cut("schira-cutoff", 'x',
                 "The distance cutoff for the Schira model force field",
                 "Because the Schira model's attraction grows stronger the farther a vertex is from"
                 " its Schira predicted origin, it is a good idea to give the field a maximum"
                 " distance over which it can operate.  This prevents noisy vertices far from the"
                 " known extent of the Schira model from having enormous energy during the"
                 " simulation.",
                 0.2);
   param_nonnegative_real
      energy("initial-energy", 'H',
             "The amount of energy to add to the system before beginning each annealing step",
             "Each step of the simulated annealing algorithm begins by heating the system by adding"
             " some amount of energy to it.  This parameter determines how much energy this is."
             " Energy is added by giving each vertex a small randomly oriented vector of a length"
             " equal to H/k where H is the initial-energy value and k is the number of vertices."
             " The net vector (sum of all random vectors) is then divided by k and subtracted from"
             " each vertex so that there is no drift in the system.",
             20.0);
   param_nonnegative_integer
      nthreads("threads", 'n',
               "The number of threads to use when running the simulation",
               "If this number is 1, the simulation is run single-threadedlyl if the number is 0,"
               " the simulation is run on k-1 threads where k is the number of cores on the host."
               " Otherwise, the simulation is run on the number of threads provided.",
               0);
   param_nonnegative_integer
      fixfc("fix-foveal-confluence", 'f',
            "Fixes the given vertex to the foveal confluence in the Schira model",
            "If this value is 0, then no foveal confluence is fixed; otherwise, fixes the foveal"
            " confluence vertex (numbered by the argument given to this option) to the point in"
            " the Schira model that should be the FC.",
            0);
   param_nonnegative_integer
      fixpd("fix-peripheral-diffluence", 'F',
            "Fixes the given vertex to the peripheral diffluence in the Schira model",
            "If this value is 0, then no peripheral diffluence is fixed; otherwise, fixes the"
            " perpheral vertex (numbered by the argument given to this option) to the point in"
            " the Schira model that should be the PD.",
            0);
   param<int>
      random_seed("random-seed", 'R',
                  "Sets the random seed for the simulation",
                  "If the argument to this option is < 0, then time(0) is used as the random seed;"
                  " otherwise, the provided value is used.",
                  -1);
   param_positive_real
      spring_cutoff("spring-cutoff", 'c',
                    "Sets the cutoff distance for creating a spring between vertices",
                    "This distance determines how close two vertices must be in order for there to"
                    " be an anatomical spring connecting them.  All vertices that are at least as"
                    " close as this cutoff value will be connected by springs.",
                    0.0125);
   param_journal
      journal_dir("journal", 'j',
                  "Sets the directory in which to journal the simulation",
                  "If provided, all relevant data about the simulation will be written to the"
                  " given directory.  This includes arguments, logs, RMSD, initial potential,"
                  " simulated annealing minima, and periodic position values.",
                  string());
   param_nonnegative_real
      schira_min_stepsz("schira-minimize-stepsize", 'm',
                        "Sets the stepsize when minimizing the Schira model to the data",
                        "If the provided value is 0, then the Schira model is not minimized to the"
                        " vertex data; otherwise, it is minimized using the given stepsize.",
                        0);
   param_nonnegative_real
      data_min_stepsz("data-minimize-stepsize", 'z',
                      "Sets the stepsize when minimizing the data to the Schira model",
                       "If the provided value is 0, then the data is not minimized to the"
                      " Schira model before and after simulated annealing iterations.  Otherwise,"
                      " it is minimized using the given stepsize.",
                      0.01);
   param_nonnegative_integer
      annealing_steps("annealing-iterations", 'a',
                      "Sets the number of simulated annealing iterations to perform",
                      "If this number is 0 or 1, do not perform simulated annealing; instead run"
                      " the simulation only once. Otherwise, run it k times (where k is the given"
                      " argument), and use the one with the lowest PE as the result.",
                      5);
   param_nonnegative_real
      max_energy_drift("max-energy-drift", 'J',
                       "The maximum energy drift allowed before energies are rescaled",
                       "Numerical drift forces the simulation to rescale energies periodically;"
                       " otherwise kinetic energy could explode.  This sets the maximum difference"
                       " between the initial energy and the energy at a given step that the"
                       " simulation may acheive before energies are rescaled.",
                       2.0);
   param<string> 
      schira_params_file("schira-param-file", 'p',
                         "The Schira model parameters filename",
                         "If provided, uses the given file as a Schira model parameters file.  This"
                         " file must list each parameter with its value in a tab separated table.",
                         string());
   param<bool> no_vdw("no-vdw", 'l',
                      "Simulate without van der Waals forces",
                      "This option speeds up the simulation considerably, but runs it without"
                      " van der Waals forces, which can result in crossing over and inversion of"
                      " the cortical surface.",
                      false);
   param<bool> verbose_arg("verbose", 'v',
                           "Use verbose output (implied by journaling)",
                           "Print verbose output to stdout (if not journaling) or to the log (if"
                           " journaling).",
                       false);
   param<bool> debug_arg("debug", 'D',
                         "Use verbose output (implied by journaling) to stderr",
                         "Print verbose output to stderr (if not journaling) or to the log (if"
                         " journaling).",
                         false);

   // add all parameters to the parameter space...
   param_space* PS = 0;
   try {
      PS = new param_space(&timestep,
                           &steps,
                           &dampen,
                           &spring_str,
                           &schira_str,
                           &schira_cut,
                           &energy,
                           &nthreads,
                           &fixfc, 
                           &fixpd,
                           &random_seed,
                           &spring_cutoff,
                           &journal_dir,
                           &schira_min_stepsz,
                           &schira_params_file,
                           &data_min_stepsz,
                           &annealing_steps,
                           &max_energy_drift,
                           &no_vdw,
                           &verbose_arg,
                           &debug_arg,
                           (void*)0);
      // parse...
      ostringstream ss;
      if (!PS->parse(argc, argv, ss)) {
         cerr << "Syntax: schira [options] <input file> <output file>" << endl
              << ss.str() << endl;
         return 0;
      }
   } catch (string s) {
      cerr << s << endl;
      return 1;
   }
   // if there are too many input args, take care of that now...
   if (argc < 2) {
      cerr << "Syntax: schira [options] <input file> <output file>" << endl;
      return 1;      
   } else if (argc > 2) {
      cerr << "Unrecognized options or too many inputs:";
      while (--argc > 0) {
         cerr << " \"" << *(argv++) << '"';
      }
      cerr << endl;
      return 1;
   }

   // At this point arguments have been parsed; let's process them

   // First, process the journal arg
   bool journal;
   if (journal_dir.provided())
      journal = true;
   else
      journal = false;
   bool verbose = (journal || *debug_arg || *verbose_arg);

   // Next, log, input, and output args...
   ofstream ofstream_err(
     journal
     ? (*journal_dir + "/log.txt").data()
     : "/dev/null");
   if (journal && !ofstream_err.is_open()) {
      cerr << "Could not open error log!" << endl;
      return 3;
   }
   ostream& err = (journal
                   ? ofstream_err
                   : (*debug_arg
                      ? cerr 
                      : (*verbose_arg
                         ? cout
                         : ofstream_err)));
   // start by stamping the log
   if (journal) {
      time_t t;
      time(&t);
      char buf[26];
      ctime_r(&t, buf);
      buf[24] = 0;
      err << "================================================================================" 
          << endl
          << "Beginning simulation on " << buf << endl;
      // also write the command line file
      ofstream f((*journal_dir + "/options.txt").data());
      if (!f.is_open()) {
         cerr << "Unable to open options file in journal!" << endl;
         return 4;
      }
      f << *PS << endl;
      f.close();
   }
   ifstream ifstream_in(*argv);
   if (!ifstream_in.is_open()) {
      cerr << "Error: could not open file: " << *argv << endl;
      return 2;
   }
   istream& in = ifstream_in;
   if (verbose) err << "Input file:  " << *argv << endl;
   argv++;
   ofstream ofstream_out(*argv);
   if (!ofstream_out.is_open()) {
      cerr << "Error: could not open file: " << *argv << endl;
      return 2;
   }
   ofstream& out = ofstream_out;
   if (verbose) err << "Output file: " << *argv << endl;

   // now, save the timestep...
   double dt = *timestep; // the duration of one step

   // and random seed...
   if (random_seed.provided() && *random_seed >= 0)
      srand(*random_seed);
   else
      srand(time(0));


   /////////////////////////////////////////////////////////////////////////////////////////////////
   // Read in the input and setup the springs
   
   // first, we get the number of points and the number of springs
   int n;
   in.read((char*)&n, sizeof(int));
   if (n < 1) {
      cerr << "Error: Input file specifies fewer than 1 atom" << endl;
      return 3;
   } else if (!in) {
      cerr << "Error reading file before loading points" << endl;
      return 3;
   }

   // next, we get all the points
   struct atom_data {double xy[2], th, r, conf;};
   atom_data* pts = new atom_data[n];
   in.read((char*)pts, sizeof(atom_data)*n);
   if (!in) {
      cerr << "Error reading file while reading point data" << endl;
      return 3;
   }

   // at this point, we're done with the input, so go ahead and clean it up of need be
   if (dynamic_cast<ofstream*>(&in))
      dynamic_cast<ofstream*>(&in)->close();

   // now, put all the data into the system...
   springs::system<2> U(*nthreads);
   for (i = 0; i < n; ++i) {
      // make sure each point is between pi/2 and -pi/2
      U.add(new springs::system<2>::atom(pts[i].xy));
   }

   // add springs for everything within cutoff distance
   i = U.make_springs(*spring_cutoff, *spring_str);
   if (verbose) err << "Added " << i << " springs to " << n << " atoms ("
                    << (double)i/(double)n << " springs per atom)." << endl;

   // add the Schira field
   vector<double> tmp_r(n);
   vector<double> tmp_th(n);
   vector<double> tmp_conf(n);
   for (i = 0; i < n; ++i) {
      tmp_r[i] = pts[i].r;
      tmp_th[i] = pts[i].th;
      tmp_conf[i] = pts[i].conf;
   }

   // add the schira model field
   schira* model = ((*schira_params_file).empty()? 
                    new schira(tmp_th, tmp_r, tmp_conf, *schira_str, 0, *schira_cut)
                    : new schira(tmp_th, tmp_r, tmp_conf, *schira_str, *schira_params_file,
                                 *schira_cut));

   // Before the other fields, we want to add the fixing field if requested
   if (*fixfc > 0)
      U.add(new coord_spring<2>((unsigned)*fixfc - 1, model->parameters().fc()));
   if (*fixpd > 0)
      U.add(new coord_spring<2>((unsigned)*fixpd - 1, model->parameters().pd()));

   U.add(model);

   // add the angle springs!
   //if (verbose) err << "Adding angle springs..." << endl;
   //angle_field* angles = new angle_field();
   //U.add(angles);
   // also, add the van-der-waals force...
   simple_vdw<2>* vdw = 0;
   if (!*no_vdw) {
      vdw = new simple_vdw<2>(-1.0);
      U.add(vdw); // the negative will force it to auto-choose a cutoff
   }


   /////////////////////////////////////////////////////////////////////////////////////////////////
   // Run the spring system

   double ke, pe, pe0, ke0, e0, *rmsd;
   if (journal) rmsd = new double[*steps];
   else rmsd = 0;
   // calculate the potential energy
   if (verbose) err << "Starting simulation..." << endl;
   pe0 = U.start(verbose, err);
   // output the potentials if desired
   if (journal) {
      double tmpx[2];
      for (i = 0; i < n; ++i) {
         if (pts[i].r >= 0 && pts[i].conf > 0) {
            tmpx[0] = pts[i].xy[0];
            tmpx[1] = pts[i].xy[1];
            pts[i].conf = model->predict(tmpx, M_PI_2 - pts[i].th*M_PI/180.0, pts[i].r, pts[i].xy);
         } else
            pts[i].conf = 0;
      }
      ofstream fl((*journal_dir + "/initial_schira_pe.bin").data());
      if (!fl.is_open()) {
         cerr << "Could not open potential file!" << endl;
         return 5;
      }
      fl.write((char*)&n, sizeof(int));
      fl.write((char*)pts, sizeof(atom_data) * n);
      fl.close();
      // return the pts[i].conf values...
      for (i = 0; i < n; ++i)
         pts[i].conf = tmp_conf[i];
   }
   if (verbose) err << "Initial potential energy: PE0 = " << pe0 << endl;
   if (verbose) model->report(err);
   // write out the initial schira model if were journaling
   if (journal && !model->write_params((*journal_dir + "/initial-schira-params.txt").data()))
      err << "Warning: Failed to output initial Schira params to file!" << endl;
   if (*schira_min_stepsz > 0) {
      if (verbose) err << "Minimizing Schira model to data..." << endl;
      model->explore_params(*nthreads, *schira_min_stepsz);
      if (verbose) err << model->parameters() << endl;
      // if we are saving the schira model params, we do that now...
      if (journal && !model->write_params((*journal_dir + "/minimized-schira-params.txt").data()))
         err << "Warning: Failed to output Schira params to file!" << endl;
      pe0 = U.PE();
      if (verbose) err << "Final potential energy: PE0 = " << pe0 << endl;
   }

   // before we start, we want to minimize...
   if (*data_min_stepsz > 0) {
      err << "Minimizing data to Schira model..." << endl;
      pe = U.minimize(*data_min_stepsz, 500);
      pe0 = U.PE();
      err << "\tNew PE0 = " << pe0 << " (" << pe << ")" << endl;
   }

   // In case we're journaling, we need to keep track of a file...
   ofstream* journal_file = 0;
   double* journal_dat = (journal? new double[2*n] : 0);
   // These variables keep track of the best simulated annealing cycle
   unsigned anneal = *annealing_steps;
   unsigned startAnneal = anneal, bestAnneal = 0;
   double bestPE = 10.0 * fabs(pe0);
   // we store the best xy positions in pts, which is initialized to the starting arrangement,
   // which also has potential energy pe0

   // The main loop!
   do {
      // first off, wipe all velocities so that all atoms are still
      U.wipe(0);
      // now, bring the system up to temperature... 
      if (*energy > 0) {
         if (verbose) err << "Heating system to goal KE = " << *energy << endl;
         U.heat(*energy);
      }
      // store the current energies...
      pe0 = U.PE();
      ke0 = U.KE();
      e0 = pe0 + ke0;
      // report, if needed
      if (verbose) {
         err << "Running simulated annealing step " << (startAnneal - anneal + 1) << ":" << endl
             << "\tKE0 = " << ke0 << ";\tPE0 = " << pe0 << ";\tE0 = " << e0 << endl;
      }
      // now, simulate it for the required number of steps
      for (i = 0; i < *steps; ++i) {
         // If we're journaling, we want to deal with coordinate dumping first.
         if (journal) {
            if (i % 1000 == 0) {
               ostringstream ss;
               unsigned jwidth = (unsigned)std::floor(std::log10((double)*steps / 1000.0)) + 1;
               ss << *journal_dir << "/t=" << std::setfill('0') << std::setw(2) << (startAnneal - anneal) 
                  << '.' << std::setfill('0') << std::setw(jwidth) << (i + 1)/1000 << ".bin";
               // we want to start a new journal file
               if (journal_file) {
                  journal_file->close();
                  delete journal_file;
                  if (verbose) {
                     if (i+1 == *steps) err << "Ending journal." << endl;
                     else err << "Transitioning to journal file: " << ss.str().data() << endl;
                  }
               } else if (verbose)
                  err << "Opening journal file: " << ss.str().data() << endl;
               journal_file = new ofstream(ss.str().data());
               if (!journal_file->is_open()) {
                  err << "Error: could not open journal file (" << ss.str().data() << ") at step " << (i+1) << endl;
                  journal_file = 0;
               } else
                  journal_file->write((char*)&n, sizeof(int));
            }
            if (journal_file && i % 10 == 0) {
               for (j = 0; j < n; ++j) {
                  journal_dat[2*j + 0] = U[j]->x()[0];
                  journal_dat[2*j + 1] = U[j]->x()[1];
               }
               journal_file->write((char*)journal_dat, sizeof(double)*2*n);
            }
         }
         // And! Step by dt!
         ke = U.step(dt, (rmsd? rmsd + i : 0));
         // dampen, if we're doing that
         if (*dampen < 1.0) U.wipe(*dampen);
         // and check for energy growth
         if (i % 10 == 9) {
            pe = U.PE();
            // if the energy has grown, rescale the kinetic energy
            if (ke + pe - e0 > *max_energy_drift) {
               err << "Rescaling energies at step " << (i + 1) << " (T = " << U.t() 
                   << "): KE = "  << ke << "; PE = " << pe << endl;
               // if the potential energy now is greater than e0, we need to minimize
               if (pe > e0) {
                  pe = U.minimize(*data_min_stepsz <= 0? 0.2 : *data_min_stepsz, 500);
                  if (verbose) err << "\tMinimized PE by " << pe << ";\t new PE = ";
                  pe = U.PE();
                  if (verbose) err << pe << endl;
               }
               // cool the velocity if we must
               if (pe + ke > e0) {
                  double tmp = ke;
                  ke = U.cool(e0 - pe);
                  if (verbose) err << "\tCooled KE by " << (ke - tmp) 
                                   << "; new KE = " << ke << endl;
               }
               // if e is still too high, we need to break; we can't do anything about this
               if (pe + ke - e0 > *max_energy_drift) {
                  err << "Energy growth uncontrolable; failing simulation." << endl;
                  break;
               }
            } else if (verbose && i % 1000 == 999) {
               err << "Simulation stable at step " << (i+1) << " (T = " << U.t() 
                   << "):\tKE = " << ke << ";\tPE = " << pe
                   << "\t[PEf = " << U.fieldPE()
                   << "; PEs = " << U.springPE() << "]" << endl;
            }
         }
      }
      // if there's a journal file open, finish it off...
      if (journal_file) {
         // write the last set of data
         for (j = 0; j < n; ++j) {
            journal_dat[2*j + 0] = U[j]->x()[0];
            journal_dat[2*j + 1] = U[j]->x()[1];
         }
         journal_file->write((char*)journal_dat, sizeof(double)*2*n);
         // and close it up!
         if (verbose) err << "Closing journal for annealing step." << endl;
         journal_file->close();
         delete journal_file;
         journal_file = 0;
      }
      // we minimize...
      if (*data_min_stepsz > 0) {
         pe = U.minimize(*data_min_stepsz, 500);
         if (verbose) err << "\tMinimized PE by " << pe << ";\t new PE = ";
         pe = U.PE();
         if (verbose) err << pe << endl;
      } else
         pe = U.PE();
      // if we're journaling, stamp this out as a normal potential file
      if (journal) {
         ostringstream ss;
         ss << *journal_dir << "/min" << (startAnneal - anneal + 1) << ".bin";
         ofstream f(ss.str().data());
         if (!f.is_open())
            err << "Unable to open journal file: " << ss.str().data() << endl;
         else {
            if (verbose) err << "Writing journal anneal minimum: " << ss.str().data() << endl;
            atom_data* ad = new atom_data[n];
            for (i = 0; i < n; ++i) {
               std::memcpy(ad + i, pts + i, sizeof(atom_data));
               ad[i].xy[0] = U[i]->x()[0];
               ad[i].xy[1] = U[i]->x()[1];
            }
            f.write((char*)&n, sizeof(int));
            f.write((char*)ad, sizeof(atom_data)*n);
            f.close();
            delete[] ad;
         }
         // we also write out an rmsd file
         ss.str("");
         ss << *journal_dir << "/rmsd" << (startAnneal - anneal + 1) << ".bin";
         f.open(ss.str().data());
         if (!f.is_open())
            err << "Unable to open RMSD file: " << ss.str().data() << endl;
         else {
            if (verbose) err << "Writing journal RMSD: " << ss.str().data() << endl;
            f.write((char*)rmsd, sizeof(double)*(*steps));
            f.close();
         }
      }
      // see how our state is... if our PE is low enough, save this state
      if (pe < bestPE) {
         bestPE = pe;
         for (i = 0; i < n; ++i) {
            pts[i].xy[0] = U[i]->x()[0];
            pts[i].xy[1] = U[i]->x()[1];
         }
         bestAnneal = startAnneal - anneal;
      }
   } while (--anneal > 0);

   // report, if verbose
   if (verbose)
      err << startAnneal << " annealing iterations finished with " << *steps << " steps each."
          << endl;

   /////////////////////////////////////////////////////////////////////////////////////////////////
   // Write output

   // we just dump the positions of the nodes in order; very simply
   //for (i = 0; i < n; ++i) {
   //   pts[i].xy[0] = U[i]->x()[0];
   //   pts[i].xy[1] = U[i]->x()[1];
   //}
   out.write((char*)&n, sizeof(int));
   out.write((char*)pts, sizeof(atom_data)*n);

   /////////////////////////////////////////////////////////////////////////////////////////////////
   // Cleanup!
   if (rmsd)
      delete[] rmsd;
   delete[] pts;
   //delete[] ss;

   return 0;
}
