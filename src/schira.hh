////////////////////////////////////////////////////////////////////////////////////////////////////
// schira.hh
// Specification of the spring system field of the Schira model's potential function
// by Noah Benson

#ifndef _____nben_SCHIRA_HH
#define _____nben_SCHIRA_HH

#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>
#include <pthread.h>
#include "springs.hh"

class schira:
   public springs::system<2>::field
{
 public:
   // the class that holds parameter sets for the schira model...
   // a parameter set class
   class param_set
   {
    private:
      // our parameters for the model
      double m_a; // the a parameter of the transform
      double m_b; // the b parameter of the transform
      double m_logab; // log(a/b)
      double m_psi; // the angle of rotation
      double m_rot[2][2]; // the rotation matrix (from psi)
      double m_irot[2][2]; // the inverse rotation matrix (from psi)
      double m_fc[2]; // the center (foveal confluence)
      double m_pd[2]; // the opposite end from the fc
      double m_scale[2]; // the x and y stretch
      double m_shear[2][2]; // the shear matrix [1 x; y 1]
      double m_ishear[2][2]; // the inverted shear matrix
      double m_sizes[4]; // sizes of various regions in the intermediate representation
      double m_lambda; // the foveal shift
      // pointers to the actual parameters (used in minimization)
      double** m_param;
      unsigned m_params;
      double* m_grad;
      double* m_steps;
      // min and max values for the parameters
      double* m_min;
      double* m_max;

      // once params are set/read, this inits the min/max/etc variables
      void init()
      {
         // fill in the param pointers
         m_params = 13;
         m_param = new double*[m_params];
         m_grad = new double[m_params];
         m_steps = new double[m_params];
         m_min = new double[m_params];
         m_max = new double[m_params];
         // param a
         m_param[0] = &m_a;
         m_steps[0] = 0.1;
         m_min[0] = 0;  m_max[0] = 50.0;
         // param b
         m_param[1] = &m_b;
         m_steps[1] = 1.0;
         m_min[1] = 10.0;
         m_max[1] = 120.0;
         // param psi (rotation)
         m_param[2] = &m_psi;
         m_steps[2] = 0.005;
         m_min[2] = -M_PI_2;
         m_max[2] = M_PI_2;
         // foveal confluence params
         m_param[3] = &m_fc[0];
         m_steps[3] = 0.01;
         m_min[3] = -1.0;
         m_max[3] = 0.0;
         m_param[4] = &m_fc[1];
         m_steps[4] = 0.01;
         m_min[4] = -0.25;
         m_max[4] = 0.25;
         m_pd[0] = -m_fc[0];
         m_pd[1] = -m_fc[1];
         // scale factors (x and y)
         m_param[5] = &m_scale[0];
         m_steps[5] = 0.005;
         m_min[5] = -2.0;
         m_max[5] = 2.0;
         m_param[6] = &m_scale[1];
         m_steps[6] = 0.005;
         m_min[6] = -2.0;
         m_max[6] = 2.0;
         // shear parameters...
         m_param[7] = &m_shear[0][1];
         m_steps[7] = 0.002;
         m_min[7] = -0.8;
         m_max[7] = 0.9;
         m_param[8] = &m_shear[1][0];
         m_steps[8] = 0.002;
         m_min[8] = -0.8;
         m_max[8] = 0.9;
         // sizes: V1, V2, V3, hV4/V3a
         m_param[9] = &m_sizes[0];
         m_steps[9] = 0.01;
         m_min[9] = M_PI_2 / 4.0;
         m_max[9] = M_PI;
         m_param[10] = &m_sizes[1];
         m_steps[10] = 0.01;
         m_min[10] = M_PI_2 / 4.0;
         m_max[10] = M_PI;
         m_param[11] = &m_sizes[2];
         m_steps[11] = 0.01;
         m_min[11] = M_PI_2 / 4.0;
         m_max[11] = M_PI;
         // Finally, lambda
         m_param[12] = &m_lambda;
         m_steps[12] = 0.01;
         m_min[12] = 0;
         m_max[12] = 2.0;
      }

      // recalculate dependency parameters given new parameters
      void recalculate()
      {
         // we want to sanity check as we go
         if (m_a <= 0) m_a = 0.001;
         if (m_b <= m_a) m_b = m_a + 0.001;
         m_logab = std::log(m_a/m_b);
         if (m_psi > 0.05) m_psi = 0.05;
         else if (m_psi < -0.05) m_psi = -0.05;
         m_rot[0][0] = std::cos(m_psi);
         m_rot[0][1] = -std::sin(m_psi);
         m_rot[1][0] = -m_rot[0][1];
         m_rot[1][1] = m_rot[0][0];
         m_irot[0][0] = m_rot[0][0];
         m_irot[0][1] = -m_rot[0][1];
         m_irot[1][0] = -m_rot[1][0];
         m_irot[1][1] = m_rot[1][1];
         if (m_shear[0][1] > 0.5) m_shear[0][1] = 0.5;
         else if (m_shear[0][1] < -0.5) m_shear[0][1] = -0.5;
         if (m_shear[1][0] > 0.5) m_shear[1][0] = 0.5;
         else if (m_shear[1][0] < -0.5) m_shear[1][0] = -0.5;
         m_ishear[0][0] = 1.0/(1.0 - m_shear[1][0]*m_shear[0][1]);
         m_ishear[0][1] = -m_shear[0][1] / (1.0 - m_shear[1][0]*m_shear[0][1]);
         m_ishear[1][0] = -m_shear[1][0] / (1.0 - m_shear[1][0]*m_shear[0][1]);
         m_ishear[1][1] = 1.0 / (1.0 - m_shear[1][0]*m_shear[0][1]);
         if (m_sizes[0] > M_PI_2) m_sizes[0] = M_PI_2;
         else if (m_sizes[0] < 0.1) m_sizes[0] = 0.1;
         if (m_sizes[1] > M_PI / 3.0) m_sizes[1] = M_PI / 3.0;
         else if (m_sizes[1] < 0.1) m_sizes[1] = 0.1;
         if (m_sizes[2] > M_PI / 4.0) m_sizes[2] = M_PI / 4.0;
         else if (m_sizes[2] < 0.1) m_sizes[2] = 0.1;
         if (m_sizes[3] > M_PI / 4.0) m_sizes[3] = M_PI / 4.0;
         else if (m_sizes[3] < 0.05) m_sizes[3] = 0.05;
         if (m_sizes[1] > m_sizes[0]) m_sizes[1] = m_sizes[0] - 0.001;
         if (m_sizes[2] > m_sizes[1]) m_sizes[2] = m_sizes[1] - 0.001;
         if (m_sizes[3] > m_sizes[2]) m_sizes[3] = m_sizes[2] - 0.001;
         if (m_lambda < 0.1) m_lambda = 0.1;
         else if (m_lambda > 4) m_lambda = 4;
         m_pd[0] = -m_fc[0];
         m_pd[1] = -m_fc[1];
      }

      void cleanup()
      {
         delete[] m_param;
         delete[] m_grad;
         delete[] m_steps;
         delete[] m_min;
         delete[] m_max;
      }

      // we're friends with the schira class
      friend class schira;

      // here we have the minimize and related functions for minimizing the parameters
      void minimize(const schira& sch, double stepsz = 0.001)
      {
         unsigned i;
         double tmp, d=0, pe;
         // basic idea: use finite differences to figure out which way to move the model
         // First off, we'll want to find the gradient.
         // for each parameter, change it a little and see what happens...
         // we treat psi separately because rot and psi are dependent
         for (i = 0; i < m_params; ++i) {
            tmp = *m_param[i];
            *m_param[i] += m_steps[i];
            recalculate();
            m_grad[i] = sch.potential(this);
            *m_param[i] = tmp - m_steps[i];
            recalculate();
            m_grad[i] -= sch.potential(this);
            m_grad[i] /= (2.0 * m_steps[i]);
            *m_param[i] = tmp;
            // remember the length of the gradient vector
            d += m_grad[i]*m_grad[i];
            // adjust the step
            //tmp = std::fabs(m_grad[i]);
            //tmp = tmp > 2.0? 2.0 : (tmp < 0.5? 0.5 : tmp);
            //m_steps[i] /= tmp;
            //if (m_steps[i] < 0.001) m_steps[i] = 0.001;
            //else if (m_steps[i] > 1.0) m_steps[i] = 1.0;
         }
         recalculate();
         pe = sch.potential(this);
         d = std::sqrt(d);
         // if the gradient is small enough, we just ignore it and don't minimize
         if (d > 0.01) {
            // we want to step half-way to the zero intersection
            tmp = pe / d / 2.0;
            for (i = 0; i < m_params; ++i)
               *m_param[i] -= (m_grad[i] / d) * tmp * stepsz;
         }
      }


    public:
      // init/create/destroy
      param_set(const char* init_file_name = 0)
      {
         // default values
         double 
            a = 2.5,
            b = 40.0,
            psi = 0, 
            fcx = -0.59731, fcy = 0,
            scalex = 0.423, scaley = 0.250,
            shearx = 0, sheary = 0,
            v1sz = 1.5709, v2sz = 0.7854, v3sz = 0.4, hv4sz = 0.8000,
            lambda = 0.8;
         if (init_file_name) {
            std::ifstream fl(init_file_name);
            if (!fl.is_open()) {
               std::cerr << "Could not open init file: " << init_file_name << endl;
               std::exit(1);
            }
            std::string name;
            double val;
            int lineno = 0;
            while (fl >> name >> val) {
               lineno++;
               if (name == "a") a = val;
               else if (name == "b") b = val;
               else if (name == "psi") psi = val;
               else if (name == "FCx") fcx = val;
               else if (name == "FCy") fcy = val;
               else if (name == "Scalex") scalex = val;
               else if (name == "Scaley") scaley = val;
               else if (name == "Shearx") shearx = val;
               else if (name == "Sheary") sheary = val;
               else if (name == "V1Size") v1sz = val;
               else if (name == "V2Size") v2sz = val;
               else if (name == "V3Size") v3sz = val;
               else if (name == "HV4Size") hv4sz = val;
               else if (name == "lambda") lambda = val;
               else if (name == "PE" || name == "potential") ; // do nothing (but not an error)
               else {
                  std::cerr << "Error reading " << init_file_name << " (line " << lineno << ")" << endl;
                  exit(1);
               }
            }
            if (!fl.eof())
               std::cerr << "Warning: Did not finish reading file " << init_file_name << "!" << endl;
            fl.close();
         }
         m_a = a;
         m_b = b;
         m_psi = psi;
         m_fc[0] = fcx;
         m_fc[1] = fcy;
         m_scale[0] = scalex;
         m_scale[1] = scaley;
         m_shear[0][0] = 1;
         m_shear[0][1] = sheary;
         m_shear[1][0] = shearx;
         m_shear[1][1] = 1;
         m_sizes[0] = v1sz;
         m_sizes[1] = v2sz;
         m_sizes[2] = v3sz;
         m_sizes[3] = hv4sz;
         m_lambda = lambda;
         // init the parameter variables...
         init();
         // make sure meta-params are correct
         recalculate();
      }
      param_set(const param_set& p)
      {
         m_a = p.m_a;
         m_b = p.m_b;
         m_psi = p.m_psi;
         m_fc[0] = p.m_fc[0];
         m_fc[1] = p.m_fc[1];
         m_scale[0] = p.m_scale[0];
         m_scale[1] = p.m_scale[1];
         m_shear[0][0] = 1;
         m_shear[0][1] = p.m_shear[0][1];
         m_shear[1][0] = p.m_shear[1][0];
         m_shear[1][1] = 1;
         m_sizes[0] = p.m_sizes[0];
         m_sizes[1] = p.m_sizes[1];
         m_sizes[2] = p.m_sizes[2];
         m_sizes[3] = p.m_sizes[3];
         m_lambda = p.m_lambda;
         // init the parameter variables...
         init();
         // make sure meta-params are correct
         recalculate();
      }
      ~param_set()
      {
         cleanup();
      }

      // copy operator...
      param_set& operator = (const param_set& p)
      {
         cleanup();
         m_a = p.m_a;
         m_b = p.m_b;
         m_psi = p.m_psi;
         m_fc[0] = p.m_fc[0];
         m_fc[1] = p.m_fc[1];
         m_scale[0] = p.m_scale[0];
         m_scale[1] = p.m_scale[1];
         m_shear[0][0] = 1;
         m_shear[0][1] = p.m_shear[0][1];
         m_shear[1][0] = p.m_shear[1][0];
         m_shear[1][1] = 1;
         m_sizes[0] = p.m_sizes[0];
         m_sizes[1] = p.m_sizes[1];
         m_sizes[2] = p.m_sizes[2];
         m_sizes[3] = p.m_sizes[3];
         m_lambda = p.m_lambda;
         init();
         recalculate();
         return *this;
      }

      // accessors
      double a() const {return m_a;}
      double b() const {return m_b;}
      double psi() const {return m_psi;}
      double xshear() const {return m_shear[0][1];}
      double yshear() const {return m_shear[1][0];}
      const double* sizes() const {return m_sizes;}
      double lambda() const {return m_lambda;}
      const double* fc() const {return m_fc;}
      const double* pd() const {return m_pd;}
      const double* scale() const {return m_scale;}
      unsigned count() const {return m_params;}
      double get(unsigned param_no) const {return *m_param[param_no];}

      // print out the parameter set...
      friend std::ostream& operator << (std::ostream& o, const param_set& p);
      // or save it in a schira-params file format
      std::ostream& save(std::ostream& o)
      {
         return o << "a\t" << m_a << std::endl
                  << "b\t" << m_b << std::endl
                  << "psi\t" << m_psi << std::endl
                  << "FCx\t" << m_fc[0] << std::endl
                  << "FCy\t" << m_fc[1] << std::endl
                  << "Scalex\t" << m_scale[0] << std::endl
                  << "Scaley\t" << m_scale[1] << std::endl
                  << "Shearx\t" << m_shear[1][0] << std::endl
                  << "Sheary\t" << m_shear[0][1] << std::endl
                  << "V1Size\t" << m_sizes[0] << std::endl
                  << "V2Size\t" << m_sizes[1] << std::endl
                  << "V3Size\t" << m_sizes[2] << std::endl
                  << "HV4Size\t" << m_sizes[3] << std::endl
                  << "lambda\t" << m_lambda << std::endl;
      }
   };

 private:
   // the atoms that matter charge-wise
   std::vector<springs::system<2>::atom*> m_atoms;
   std::vector<unsigned> m_atoms_idx;
   // and their charges
   std::vector<double> m_theta;
   std::vector<double> m_r;
   // and the confidence of each
   std::vector<double> m_str;

   // Finally, the strength of the field's general spring; this is the max force on an obj
   double m_strength;
   // If we have a distance cutoff this will be > 0
   double m_dist_cutoff;

   // the parameters!
   param_set m_params;

 private:

   // a way to approximate the sech function
   inline double sech(double x) const {
      x = std::exp(x);
      return 2.0 / (x + 1.0 / x);
   }
   inline double sech_pow(double x, double y) const {
      x = std::exp(x);
      return std::pow(2.0 / (x + 1.0 / x), y);
   }

   // calculate the x/y values for a give r/theta, ignoring v1/v2/v3 magnification.
   // the result is stored in the xy vector
   void forward(double th, double r, double* xy, const param_set* params = 0) const
   {
      const param_set& p = (params == 0? m_params : *params);
      double x[2];
      double tmp[2];
      x[0] = r * std::cos(th);
      x[1] = r * std::sin(th);
      // first we do banding...
      if (x[0] >= 0) x[0] += p.m_lambda;
      else x[0] += 2.0 * p.m_lambda * (1.0 - (th > 0? th : -th) / M_PI);
      th = std::atan2(x[1], x[0]);
      r = std::sqrt(x[0]*x[0] + x[1]*x[1]);
      // then, we do the two double-sech transforms
      //   store the abs for now in r and the exponent of the transform in x[0] and x[1]
      x[0] = 0.1821 * sech(0.76 * std::log(r/p.m_a));
      x[1] = 0.1821 * sech(0.76 * std::log(r/p.m_b));
      //   convert theta into the corrected value for theta
      x[0] = th * sech_pow(th, x[0]);
      x[1] = th * sech_pow(th, x[1]);
      // now we do the log-polar transform
      // we can be tricky with this keeping in mind a few things:
      //  1) x[0] is really the theta for z_a
      //  2) x[1] is really the theta for z_b
      //  3) r is the abs for both z_a and z_b
      //  4) log(r*exp(i*th)) == log(r) + i*th
      tmp[0] = r * std::cos(x[0]) + p.m_a;
      tmp[1] = r * std::cos(x[1]) + p.m_b;
      x[0] = r * std::sin(x[0]);
      x[1] = r * std::sin(x[1]);
      xy[0] = 0.5 * std::log((tmp[0]*tmp[0] + x[0]*x[0]) / (tmp[1]*tmp[1] + x[1]*x[1]));
      xy[1] = std::atan2(x[0], tmp[0]) - std::atan2(x[1], tmp[1]);
      // at the end, we perform the simpler transformations...
      x[0] = (xy[0] - p.m_logab) * p.m_scale[0];
      x[1] = xy[1] * p.m_scale[1];
      xy[0] = p.m_shear[0][0] * x[0] + p.m_shear[0][1] * x[1];
      xy[1] = p.m_shear[1][0] * x[0] + p.m_shear[1][1] * x[1];
      x[0] = p.m_rot[0][0] * xy[0] + p.m_rot[0][1] * xy[1];
      x[1] = p.m_rot[1][0] * xy[0] + p.m_rot[1][1] * xy[1];
      xy[0] = x[0] + p.m_fc[0];
      xy[1] = x[1] + p.m_fc[1];
   }

   // these functions are the same as forward but assume that th is in v1 v2 or v3
   inline void forward1(double th, double r, double* xy, const param_set* params = 0) const
   {
      const param_set& p = (params == 0? m_params : *params);
      forward(th / M_PI_2 * p.m_sizes[0], r, xy, params);
   }
   inline void forward2(double th, double r, double* xy, const param_set* params = 0) const
   {
      const param_set& p = (params == 0? m_params : *params);
      forward((th > 0? M_PI_2 - th : -M_PI_2 - th) / M_PI_2 * p.m_sizes[1] +
                 (th > 0? p.m_sizes[0] : -p.m_sizes[0]),
              r, xy, params);
   }
   inline void forward3(double th, double r, double* xy, const param_set* params = 0) const
   {
      const param_set& p = (params == 0? m_params : *params);
      forward(th / M_PI_2 * p.m_sizes[2] + 
                 (th > 0? p.m_sizes[0] + p.m_sizes[1] : -(p.m_sizes[0] + p.m_sizes[1])),
              r, xy, params);
   }
   inline void forward4(double th, double r, double* xy, const param_set* params = 0) const
   {
      const param_set& p = (params == 0? m_params : *params);
      forward((th > 0? M_PI_2 - th : -M_PI_2 - th) / M_PI_2 * p.m_sizes[3] +
                 (th > 0? 1 : -1) * (p.m_sizes[0] + p.m_sizes[1] + p.m_sizes[2]),
              r, xy, params);
      //forward(-th / M_PI * p.m_sizes[3] - 
      //           (p.m_sizes[3]/2.0 + p.m_sizes[0] + p.m_sizes[1] + p.m_sizes[2]), 
      //        r, xy, params);
   }

 public:
   // given an atom, this finds the location of the best x/y value for it's r/theta
   // (ie, the closest x/y value) and returns that value by reference in xy; the
   // potential between them is returned
   double potential(const double* x0, double th, double r, double s, 
                    double* xy, double* vec = 0, const param_set* params = 0) const
   {
      double bestd;
      double x[2], tmp;
      forward1(th, r, xy, params);
      bestd = (xy[0] - x0[0])*(xy[0] - x0[0]) + (xy[1] - x0[1])*(xy[1] - x0[1]);
      forward2(th, r, x, params);
      tmp = (x[0] - x0[0])*(x[0] - x0[0]) + (x[1] - x0[1])*(x[1] - x0[1]);
      if (tmp < bestd) {
         bestd = tmp;
         xy[0] = x[0];
         xy[1] = x[1];
      }
      forward3(th, r, x, params);
      tmp = (x[0] - x0[0])*(x[0] - x0[0]) + (x[1] - x0[1])*(x[1] - x0[1]);
      if (tmp < bestd) {
         bestd = tmp;
         xy[0] = x[0];
         xy[1] = x[1];
      }
      forward4(th, r, x, params);
      tmp = (x[0] - x0[0])*(x[0] - x0[0]) + (x[1] - x0[1])*(x[1] - x0[1]);
      if (tmp < bestd) {
         bestd = tmp;
         xy[0] = x[0];
         xy[1] = x[1];
      }
      bestd = std::sqrt(bestd);
      // we implement a distance cutoff
      //if (vec) {
      //   // distance factor is included in the (xy - x0) part
      //   vec[0] = (xy[0] - x0[0]) * s * m_strength;
      //   vec[1] = (xy[1] - x0[1]) * s * m_strength;
      //}
      //return s * m_strength * bestd*bestd;
      if (vec) {
         vec[0] = 2.0 * (xy[0] - x0[0]) * exp(-64.0 * (xy[0] - x0[0])*(xy[0] - x0[0])) * s * m_strength;
         vec[1] = 2.0 * (xy[1] - x0[1]) * exp(-64.0 * (xy[1] - x0[1])*(xy[1] - x0[1])) * s * m_strength;
      }
      return s * m_strength / 81.0 * (1.0 - exp(-64.0 * bestd*bestd));
   }
 private:
   inline double potential(unsigned atom_idx, double* xy, double* grad = 0,
                           const param_set* params = 0) const
   {
      const springs::system<2>::atom* a = m_atoms[atom_idx];
      double r = m_r[atom_idx];
      double th = m_theta[atom_idx];
      double s = m_str[atom_idx];
      return potential(a->x(), th, r, s, xy, grad, params);
   }
   inline double potential(unsigned atom_idx, const param_set* params = 0) const
   {
      const springs::system<2>::atom* a = m_atoms[atom_idx];
      double r = m_r[atom_idx];
      double th = m_theta[atom_idx];
      double s = m_str[atom_idx];
      double xy[2];
      return potential(a->x(), th, r, s, xy, 0, params);
   }

 public:
   schira(const std::vector<double>& th, const std::vector<double>& r, const std::vector<double>& conf,
          double strength, const char* init_file_name = 0, double dist_cutoff = 0):
      m_theta(th),
      m_r(r),
      m_str(conf),
      m_strength(strength),
      m_dist_cutoff(dist_cutoff),
      m_params(init_file_name)
   {
   }
   schira(const std::vector<double>& th, const std::vector<double>& r, const std::vector<double>& conf,
          double strength, const std::string& init_file_name, double dist_cutoff = 0):
      m_theta(th),
      m_r(r),
      m_str(conf),
      m_strength(strength),
      m_dist_cutoff(dist_cutoff),
      m_params(init_file_name.data())
   {
   }
   ~schira()
   {
   }

   // get the parameters...
   const param_set& parameters() const {return m_params;}
   
   // these are handy utility functions for getting the ends of the v1 region
   double* getFC(double* x) const
   {
      if (x) {
         // this is in the parameters...
         x[0] = m_params.m_fc[0];
         x[1] = m_params.m_fc[1];
      }
      return x;
   }

   // return the distance cutoff
   double distance_cutoff() const {return m_dist_cutoff;}

   // write out a parameters file
   bool write_params(const char* flnm)
   {
      std::ofstream f(flnm);
      if (!f.is_open()) return false;
      m_params.save(f);
      if (f.fail() || f.bad()) {
         f.close();
         return false;
      }
      f.close();
      return true;
   }

   // initialize by finding all atoms with r/theta values
   void init(const std::vector<springs::system<2>::atom*>& all, 
             const std::vector<springs::system<2>::spring*>& all_springs)
   {
      double xy[2], d;
      unsigned k = 0;
      for (unsigned i = 0; i < all.size(); ++i) {
         if (m_str[i] > 0) {
            m_theta[i] = M_PI_2 - m_theta[i] * M_PI / 180.0;
            predict(all[i]->x(), m_r[i], m_theta[i], m_str[i], xy);
            d = std::sqrt((xy[0] - all[i]->x()[0])*(xy[0] - all[i]->x()[0]) +
                          (xy[1] - all[i]->x()[1])*(xy[1] - all[i]->x()[1]));
            if (m_dist_cutoff <= 0 || d <= m_dist_cutoff) {
               m_r[k] = m_r[i];
               m_theta[k] = m_theta[i];
               m_str[k] = m_str[i];
               m_atoms.push_back((springs::system<2>::atom*)all[i]);
               m_atoms_idx.push_back(k);
               ++k;
               continue;
            }
         }
         m_atoms_idx.push_back(all.size() + 1);
      }
      m_r.resize(k);
      m_theta.resize(k);
      m_str.resize(k);
   }

   void update(unsigned idx)
   {
      double x[2], grad[2]; // r, th;
      unsigned i = m_atoms_idx[idx];
      if (i > m_atoms_idx.size()) return;
      springs::system<2>::atom* a = m_atoms[i];
      //r = m_r[i];
      //th = m_theta[i];
      // First, calculate the best schira prediction for the x/y value of this atom, based
      // on this atom's r/theta values
      potential(i, x, grad);
      // that gives the potential and the gradient (as grad)
      a->a()[0] += grad[0];
      a->a()[1] += grad[1];
   }

   double potential(const param_set* params) const
   {
      double tot = 0;
      for (unsigned i = 0; i < m_atoms.size(); ++i)
         tot += potential(i, params);
      return tot;
   }
   double potential() const
   {
      return potential(&m_params);
   }

   // minimize attempts to change the schira parameters to make the model fit better...
   void minimize(param_set* p, double stepsz = 0.001)
   {
      double old_dist_cutoff = m_dist_cutoff;
      m_dist_cutoff = 0;
      p->minimize(*this, stepsz);
      m_dist_cutoff = old_dist_cutoff;
   }
   void minimize(double stepsz = 0.001)
   {
      double old_dist_cutoff = m_dist_cutoff;
      m_dist_cutoff = 0;
      m_params.minimize(*this, stepsz);
      m_dist_cutoff = old_dist_cutoff;
   }

   // minimize until the variance of the last 5 RSS measurements is near 0
   double full_minimize(param_set* params, double stepsz = 0.001)
   {
      double last5steps[5];
      double tmp1=0, tmp2=0, pe=0, pe0;
      unsigned ii, jj;

      pe0 = potential(params);
      for (ii = 0, jj = 0; ii < 200; ++ii, jj = (jj + 1) % 5) {
         minimize(params, stepsz);
         pe = potential(params);
         if (ii >= 5) {
            tmp1 -= last5steps[jj];
            tmp2 -= last5steps[jj]*last5steps[jj];
         }
         last5steps[jj] = pe / pe0;
         tmp1 += last5steps[jj];
         tmp2 += last5steps[jj]*last5steps[jj];
         // see if the standard deviation is below 0.01
         if (!std::isfinite(tmp2) || (ii >= 5 && tmp2 / 5.0 - tmp1*tmp1 / 25.0 < 0.0001))
            break;
      }
      return potential(params);
   }
   double full_minimize(double stepsz = 0.001)
   {
      double last5steps[5];
      double tmp1, tmp2, pe, pe0;
      unsigned ii, jj;

      pe0 = potential();
      for (ii = 0, jj = 0; ii < 200; ++ii, jj = (jj + 1) % 5) {
         minimize(stepsz);
         pe = potential();
         if (ii >= 5) {
            tmp1 -= last5steps[jj];
            tmp2 -= last5steps[jj]*last5steps[jj];
         }
         last5steps[jj] = pe / pe0;
         tmp1 += last5steps[jj];
         tmp2 += last5steps[jj]*last5steps[jj];
         // see if the standard deviation is below 0.01
         if (!std::isfinite(tmp2) || (ii >= 5 && tmp2 / 5.0 - tmp1*tmp1 / 25.0 < 0.0001))
            break;
      }
      return potential();
   }

  private:
   // these data structures and functions are used when  exploring parameters...
   struct thread_job
   {
      schira::param_set params;
      double stepsz;
      thread_job* next;
      thread_job(const schira::param_set& p, double ss, thread_job* n = 0):
         params(p), stepsz(ss), next(n)
      {
      }
   };
   struct thread_info 
   {
      // basic data for passing along
      schira* caller;
      pthread_cond_t* cond;
      pthread_mutex_t* mux;
      pthread_cond_t* scond;
      int* status;

      // the linked list of jobs
      thread_job** head;

      // the best score and best params
      double* best_score;
      schira::param_set* best_params;

      thread_info(schira* s, pthread_cond_t* c, pthread_mutex_t* m, pthread_cond_t* sc, 
                  int* stat, thread_job** h, double* bestsc, schira::param_set* bestp):
         caller(s), cond(c), mux(m), scond(sc), status(stat), head(h), 
         best_score(bestsc), best_params(bestp)
      {
      }
   };
   void explore_params_worker(pthread_cond_t* cond, pthread_mutex_t* mux, pthread_cond_t* scond,
                              int& status, thread_job*& head,
                              double& best_score, schira::param_set& best_params)
   {
      double sc;
      thread_job* j;
      pthread_mutex_lock(mux);
      while (status) {
         // see if there's a job available
         if (head) {
            // there is a job... go ahead and grab it...
            j = head;
            head = head->next;
            // we can release the mutex now and run the job...
            pthread_mutex_unlock(mux);

            // here, we actually run the job we've been given...
            sc = full_minimize(&j->params, j->stepsz);

            // done--reacquire the mutex to see if there is more to do
            pthread_mutex_lock(mux);
            // we have the mutex... if the best score is an improvement, we edit this...
            if (sc < best_score) {
               best_score = sc;
               best_params = j->params;
            }
            // we can clean up that job...
            delete j;
         } else {
            // we're about to stop, so decrement status
            status--;
            // if we just finished the last job, we signal the other cond
            if (status == 1) {
               pthread_cond_signal(scond);
            }
            // then, wait for a signal that tells us a job is available
            if (pthread_cond_wait(cond, mux) != 0)
               break;
            // and we woke up, so increment status
            if (status) status++;
         }         
      }
      pthread_mutex_unlock(mux);
   }
   static void* explore_thread(thread_info* ti)
   {
      ti->caller->explore_params_worker(ti->cond, ti->mux, ti->scond, *ti->status,
                                        *ti->head, *ti->best_score, *ti->best_params);
      delete ti;
      return 0;
   }

  public:
   // this function allows us to explore a wide range of possible parameters using 
   // the minimize function
   void explore_params(unsigned threads = 0, double stepsz = 0.001)
   {
      param_set best = m_params;
      double pe0 = potential();
      double sc = pe0;
      double tmp1;
      unsigned ii, jj, i;
      int thread_status = 1 + threads;

      // before we get started, we want to allocate threads (we may have none)
      thread_job* jobs = 0, *j;
      pthread_t* thread = 0;
      thread_info** thread_data = 0;
      pthread_cond_t* cond = 0;
      pthread_cond_t* scond = 0;
      pthread_mutex_t* mux = 0;
      if (threads > 0) {
         thread = new pthread_t[threads];
         thread_data = new thread_info*[threads];
         cond = new pthread_cond_t;
         mux = new pthread_mutex_t;
         scond = new pthread_cond_t;
         pthread_cond_init(cond, 0);
         pthread_mutex_init(mux, 0);
         pthread_cond_init(scond, 0);
         pthread_mutex_lock(mux);
         for (i = 0; i < threads; ++i) {
            thread_data[i] = new thread_info(this, cond, mux, scond, &thread_status, 
                                             &jobs, &sc, &best);
            pthread_create(&thread[i], 0, (void*(*)(void*))&explore_thread, (void*)thread_data[i]);
         }
         pthread_mutex_unlock(mux);
      }

      // we start off with the best as the current configuration (first 3 lines of this fun)

      // okay, for some number of iterations...
      unsigned num_iterations = 12;
      // we want to create a bunch of jobs starting at the current best position...
      unsigned num_directions = 28;
      // first, grab the mutex...
      if (threads > 0) pthread_mutex_lock(mux);
      for (i = 0; i < num_iterations; ++i) {
         // then create all the jobs...
         for (jj = 0; jj < num_directions; ++jj) {
            jobs = new thread_job(best, stepsz, jobs);
            for (ii = 0; ii < best.m_params; ++ii)
               *jobs->params.m_param[ii] = *best.m_param[ii] + ((double)std::rand()/(double)RAND_MAX - 0.5) * 10.0 * best.m_steps[ii];
         }
         // then, signal all the threads and release the mutex
         if (threads == 0) {
            for (j = jobs; j; j = jobs) {
               jobs = j->next;
               tmp1 = full_minimize(&j->params, stepsz);
               // we have the mutex... if the best score is an improvement, we edit this...
               if (tmp1 < sc) {
                  sc = tmp1;
                  best = j->params;
               }
               // we can clean up that job...
               delete j;
            }
         } else {
            pthread_cond_broadcast(cond);
            // then wait for these to finish...
            pthread_cond_wait(scond, mux);
         }
      }

      // at the end, use the best starting locations for all params
      m_params = best;

      if (threads > 0) {
         thread_status = 0;
         pthread_cond_broadcast(cond);
         pthread_mutex_unlock(mux);
         for (i = 0; i < threads; ++i)
            pthread_join(thread[i], 0);
         pthread_cond_destroy(cond);
         pthread_mutex_destroy(mux);
         pthread_cond_destroy(scond);
         delete cond;
         delete mux;
         delete scond;
         delete[] thread;
         delete[] thread_data;
      }
   }

   // report our parameter values
   void report(std::ostream& o) const
   {
      o << m_params << std::endl;
   }

   // predicts the xy values for a given r/theta/strength.  Returns the potential.
   double predict(const double* x0, double th, double r, double str, double* xy) const
   {
      return potential(x0, th, r, str, xy);
   }
   double predict(const double* x0, double th, double r, double* xy) const
   {
      return potential(x0, th, r, 1.0, xy);
   }
};

inline std::ostream& operator << (std::ostream& o, const schira::param_set& p)
{
   o << "Schira:{"
     << "a = " << p.m_a
     << "; b = " << p.m_b
     << "; psi = " << p.m_psi
     << "; shear = [1 " << p.m_shear[0][1] << "; " << p.m_shear[1][0] << " 1]"
     << "; scale = [" << p.m_scale[0] << " " << p.m_scale[1] << "]"
     << "; fc = [" << p.m_fc[0] << " " << p.m_fc[1] << "]"
     << "; sizes = [" << p.m_sizes[0] << " " << p.m_sizes[1] << " " << p.m_sizes[2] << "]"
     << "; lambda = " << p.m_lambda << '}';
   return o;
}

#endif
