////////////////////////////////////////////////////////////////////////////////////////////////////
// springs.hh
// Data structures and definitions for systems of forces that get minimized via simulation.
// by Noah Benson

#ifndef _____nben_SPRINGS_HH
#define _____nben_SPRINGS_HH

#include <vector>
#include <list>
#include <cmath>
#include <iostream>
#include <pthread.h>
#include "dhash.hh"

namespace springs
{
   using std::sqrt;

   // the whole spring system is encapsulated in this class, which tracks everything for the sim
   template <unsigned DIMS = 3> // the space class in which this system operates
   class system
   {
    public:
      class atom;
      class spring;
      class field;
      
      // first, the class that represents an object or atom or whatever is being simulated
      class atom
      {
       private:
         // we have some friends...
         friend class system::spring;
         friend class system::field;
         friend class system;

         // the mass of the atom
         double m_inertia;
         // the position of the atom
         double m_x[DIMS];
         // the velocity of the atom
         double m_v[DIMS];
         // and the acceleration
         double m_a[DIMS];
         // finally, the index
         unsigned m_idx;
         // a list of springs...
         std::list<typename system::spring*> m_springs;
         

       public:
         atom(double* x0 = 0, double inertia = 1.0)
         {
            m_inertia = inertia;
            for (unsigned i = 0; i < DIMS; ++i) {
               m_x[i] = (x0? x0[i] : 0.0);
               m_v[i] = 0;
               m_a[i] = 0;
            }
            m_idx = 0; // set in initialization
         }

         // accessors/mutators
         double inertia() const {return m_inertia;}
         void inertia(double d) {m_inertia = d;}
         const double* x() const {return m_x;}
         double* x() {return m_x;}
         const double* v() const {return m_v;}
         double* v() {return m_v;}
         const double* a() const {return m_a;}
         double* a() {return m_a;}
         unsigned index() const {return m_idx;}
         const std::list<typename system::spring*>& springs() const {return m_springs;}

         // returns the kinetic energy of the atom
         double KE() const
         {
            double d = 0;
            for (unsigned i = 0; i < DIMS; ++i)
               d += m_v[i]*m_v[i];
            return m_inertia * d;
         }

         // this function updates the spring based on the amount of time that has passed
         // and returns the square of the distance moved
         double update(double dt)
         {
            double d = 0, tmp;
            unsigned i;
            for (i = 0; i < DIMS; ++i) {
               tmp = m_x[i];
               m_x[i] += 0.5 * m_v[i] * dt;
               m_v[i] += m_a[i] * dt / m_inertia;
               m_x[i] += 0.5 * m_v[i] * dt;
               m_a[i] = 0;
               d += (m_x[i] - tmp)*(m_x[i] - tmp);
            }
            return d;
         }

         // halt the velocity
         void halt()
         {
            for (unsigned i = 0; i < DIMS; ++i)
               m_v[i] = 0;
         }

         // dampen the velocity
         void dampen(double str)
         {
            for (unsigned i = 0; i < DIMS; ++i)
               m_v[i] *= str;
         }
      };

      // next, the class that represents a single bond/spring
      class spring
      {
       private:
         // the atoms connected by the bond
         atom* m_a;
         atom* m_b;
         // the strength of the bond
         double m_k;
         // the ideal bond-length
         double m_x;

         // our friends...
         friend class system;
       public:
         spring(atom* a, atom* b, double x, double k = 1.0)
         {
            m_a = a;
            m_b = b;
            m_x = x;
            m_k = k;
         }

         // simple accessors
         inline double strength() const {return m_k;}
         inline double length() const {return m_x;}
         inline const atom* a() const {return m_a;}
         inline const atom* b() const {return m_b;}

         // this function updates the spring
         void update()
         {
            unsigned i;
            double d = 0;
            double v[DIMS];
            for (i = 0; i < DIMS; ++i) {
               v[i] = m_a->m_x[i] - m_b->m_x[i];
               d += v[i]*v[i];
            }
            d = sqrt(d);
            // and here's the force... we set v to be the force incident on a;
            // -v is the force incident on b
            for (i = 0; i < DIMS; ++i) {
               v[i] /= d;
               v[i] *= m_k * (m_x - d);
               m_a->m_a[i] += v[i];
               m_b->m_a[i] -= v[i];
            }
            // done!
         }

         // the potential of this spring
         double potential() const
         {
            double d = 0;
            for (unsigned i = 0; i < DIMS; ++i)
               d += (m_a->m_x[i] - m_b->m_x[i])*(m_a->m_x[i] - m_b->m_x[i]);
            d = std::sqrt(d);
            return 0.5 * m_k * (d - m_x)*(d - m_x);
         }

      };

      // the field class; represents (abstractly) a field
      class field
      {
       public:
         virtual ~field() {}
         // this function is called at the beginning of the simulation; it is suggested that
         // whatever subset of all the atoms the field might need be saved.
         virtual void init(const std::vector<atom*>& all_atoms, const std::vector<spring*>& all_springs) = 0;
         // this function is called once each iteration; it should update the given atom with 
         // the force field.
         virtual void update(unsigned atm_idx) = 0;
         // this function is called after the update function; it is called exactly once per field
         // and is run single-threadedly no matter how many threads are used in the program
         virtual void post_update() {}
         // this function is called at the end of the simulation
         virtual void finish() {}
         // calculate the potential energy of the field
         virtual double potential() const = 0;
      };

    private:
      // the atoms themselves
      std::vector<atom*> m_atoms;
      // and the springs
      std::vector<spring*> m_springs;
      // and any fields
      std::vector<field*> m_fields;
      // the amount of time passed
      double m_t;
      // the number of threads and the threads themselves
      unsigned m_threads;
      pthread_t** m_thread;
      // the condition we signal and its mutex and the return condition
      pthread_mutex_t m_mux;
      pthread_cond_t m_cond;
      pthread_cond_t m_done;
      // keeps track of working threads during threaded steps
      unsigned m_working;
      // also used to communicate between the threads and the front-end
      double m_dt;
      double m_result_ke;
      double m_result_d;
      // this is set to true when the worker threads need to know that they are minimizing instead
      // of integrating
      bool m_minimizing;
      // set to true when we're halting
      bool m_halt;

      // this function exists as a way to multi-thread spring updates; it updates all the springs
      // for a single atom without modifying data in any other atom (or spring)
      void update_atom_springs(unsigned i)
      {
         unsigned j;
         atom* a = m_atoms[i];
         double d;
         double v[DIMS];
         spring* s;
         typename std::list<typename system::spring*>::iterator it;
         for (it = a->m_springs.begin(); it != a->m_springs.end(); it++) {
            s = *it;
            d = 0;
            for (j = 0; j < DIMS; ++j) {
               v[j] = s->m_a->m_x[j] - s->m_b->m_x[j];
               d += v[j]*v[j];
            }
            d = sqrt(d);
            // and here's the force... we set v to be the force incident on a;
            // -v is the force incident on b
            for (j = 0; j < DIMS; ++j) {
               v[j] /= d;
               v[j] *= s->m_k * (s->m_x - d);
               if (s->m_a == a)
                  s->m_a->m_a[j] += v[j];
               else
                  s->m_b->m_a[j] -= v[j];
            }
         }
      }

      // the actual function that is run for threads
      void worker_loop(unsigned id)
      {
         // before we start the initialization, we need to wait for permission from the master thread
         pthread_mutex_lock(&m_mux);
         pthread_cond_wait(&m_cond, &m_mux);
         pthread_mutex_unlock(&m_mux);
         double d=0, ke=0;
         unsigned i, j, fields = m_fields.size();
         // we can go ahead and make our list of atoms (that we will be updating)
         unsigned wards = 0;
         unsigned* ward = new unsigned[(m_atoms.size() + m_threads - 1) / m_threads];
         for (i = id, wards = 0; i < m_atoms.size(); i += m_threads, ++wards)
            ward[wards] = i;
         // Okay, our initialization is done; now we grab the mutex and wait on a step
         pthread_mutex_lock(&m_mux);
         pthread_cond_wait(&m_cond, &m_mux);
         // the simulation has started... we loop until it halts
         while (!m_halt) {
            // if the worker threads count is not reset, do that now...
            if (m_working == 0) {
               // the first thread gets to initializa the each-step variables
               m_result_d = 0;
               m_result_ke = 0;
               m_working = m_threads;
            }
            // now, we can release the mutex and process our wards
            pthread_mutex_unlock(&m_mux);
            // first, we go through our wards and update the springs and the fields
            d = 0;
            for (i = 0; i < wards; ++i) {
               update_atom_springs(ward[i]);
               for (j = 0; j < fields; ++j)
                  m_fields[j]->update(ward[i]);
               // if we're minimizing, we need to count the total acceleration
               if (m_minimizing) {
                  for (j = 0; j < DIMS; ++j)
                     d += m_atoms[ward[i]]->a()[j]*m_atoms[ward[i]]->a()[j];
               }
            }
            // great; now wait on all workers to finish; for this, we grab the mutex and decrement
            // the number of working threads
            pthread_mutex_lock(&m_mux);
            // if we're minimizing, count up the result d
            if (m_minimizing) m_result_d += d;
            // wait on the appropriate condition...
            if (m_working == 1) {
               // all workers are finished; we reset the workers count for the next round
               m_working = m_threads;
               if (m_minimizing) {
                  // if we are minimizing, we want to wait on the main thread, which will signal
                  // back by broadcasting on m_cond
                  pthread_cond_signal(&m_done);
                  pthread_cond_wait(&m_cond, &m_mux);
               } else {
                  // okay, start the next round
                  pthread_cond_broadcast(&m_cond);
               }
            } else {
               // there's still at least one worker going; wait on it
               --m_working;
               pthread_cond_wait(&m_cond, &m_mux);
            }
            // okay, we can release the mux already...
            pthread_mutex_unlock(&m_mux);
            // now, we just need to update our wards themselves
            d = 0;
            ke = 0;
            if (m_minimizing) {
               for (i = 0; i < wards; ++i) {
                  // use m_result_d, which tells us the distance to travel
                  for (j = 0; j < DIMS; ++j) {
                     m_atoms[ward[i]]->v()[j] = m_result_d * m_atoms[ward[i]]->a()[j];
                     m_atoms[ward[i]]->a()[j] = 0;
                  }
                  m_atoms[ward[i]]->update(1.0);
               }
            } else {
               for (i = 0; i < wards; ++i) {
                  // update the atom
                  d += m_atoms[ward[i]]->update(m_dt);
                  ke += m_atoms[ward[i]]->KE();
               }
            }
            // right, all done; we relock the mutex to merge with the other workers again
            pthread_mutex_lock(&m_mux);
            m_result_d += d;
            m_result_ke += ke;
            if (m_working == 1) {
               // we signal the master thread that we're done!
               pthread_cond_signal(&m_done);
               // reset m_working for the next step
               m_working = 0;
            } else
               --m_working;
            // done with all of the updates!  Wait on the next step.
            pthread_cond_wait(&m_cond, &m_mux);
         }
         // we've finished the simulation; just release the mutex and return nicely
         pthread_mutex_unlock(&m_mux);
         delete ward;
      }
      // and the startup function for threads
      static void* thread_startup(std::pair<system*,unsigned>* p)
      {
         p->first->worker_loop(p->second);
         delete p;
         return 0;
      }

    public:
      system(unsigned nthreads = 0) 
      {
         m_t = 0;
         m_halt = false;
         m_threads = nthreads;
         m_working = 0;
         m_minimizing = false;
         pthread_cond_init(&m_cond, 0);
         pthread_cond_init(&m_done, 0);
         pthread_mutex_init(&m_mux, 0);
         if (m_threads == 0) {
            m_thread = 0;
         } else {
            pthread_mutex_lock(&m_mux);
            m_thread = new pthread_t*[nthreads];
            for (unsigned i = 0; i < nthreads; ++i) {
               m_thread[i] = new pthread_t;
               pthread_create(m_thread[i], 0, (void* (*)(void*))&thread_startup, (void*)new std::pair<system*,unsigned>(this, i));
            }
            pthread_mutex_unlock(&m_mux);
         }
      }
      ~system()
      {
         unsigned i;
         pthread_mutex_lock(&m_mux);
         // threads...
         if (m_threads) {
            m_halt = true;
            pthread_cond_broadcast(&m_cond);
            pthread_mutex_unlock(&m_mux);
            for (i = 0; i < m_threads; ++i)
               pthread_join(*m_thread[i], 0);
            pthread_mutex_lock(&m_mux);
         }
         for (i = 0; i < m_atoms.size(); ++i)
            delete m_atoms[i];
         for (i = 0; i < m_springs.size(); ++i)
            delete m_springs[i];
         for (i = 0; i < m_fields.size(); ++i)
            delete m_fields[i];
         pthread_mutex_unlock(&m_mux);
         pthread_mutex_destroy(&m_mux);
         pthread_cond_destroy(&m_cond);
         pthread_cond_destroy(&m_done);
      }

      // for adding atoms/springs/fields
      void add(atom* a) {m_atoms.push_back(a);}
      void add(atom* a, atom* b, double x, double k = 0) {m_springs.push_back(new spring(a, b, x, k));}
      void add(field* f) {m_fields.push_back(f);}
      void add(spring* s) 
      {
         m_springs.push_back(s);
         s->m_a->m_springs.push_front(s);
         s->m_b->m_springs.push_front(s);
      }

      // make springs from a distance---all atoms in the system that are within the given distance are 
      // joined by springs whose ideal length is their current distance
      unsigned make_springs(double d, double k = 1.0)
      {
         // make a distance hash for this part...
         unsigned i, j, count = 0;
         atom* a, *b;
         double mn[DIMS], mx[DIMS], ideal[DIMS];
         double dist;
         for (i = 0; i < DIMS; ++i) {
            mn[i] = 10000.0;
            mx[i] = -10000.0;
         }
         // first, find the min and max x/y values
         for (i = 0; i < m_atoms.size(); ++i) {
            a = m_atoms[i];
            for (j = 0; j < DIMS; ++j) {
               if (a->m_x[j] < mn[j]) mn[j] = a->m_x[j];
               if (a->m_x[j] > mx[j]) mx[j] = a->m_x[j];
            }
         }
         mn[0] -= (mx[0]-mn[0])*0.25;
         mn[1] -= (mx[1]-mn[1])*0.25;
         mx[0] += (mx[0]-mn[0])*0.25;
         mx[1] += (mx[1]-mn[1])*0.25;
         ideal[0] = 0.1; ideal[1] = 0.1;
         dhash<atom*,2> dh(mn, mx, ideal);
         // okay, now add all the atoms to the distance hash
         for (i = 0; i < m_atoms.size(); ++i)
            dh.insert(m_atoms[i], m_atoms[i]->m_x);
         // and then, for each atom, find the closest ones...
         std::list<typename dhash<atom*,2>::point*> l;
         typename dhash<atom*,2>::point* p;
         for (i = 0; i < m_atoms.size(); ++i) {
            a = m_atoms[i];
            dh.find(a->m_x, d, l);
            while (!l.empty()) {
               p = l.front();
               l.pop_front();
               b = p->get();
               // skip it if its the atom in question (distance is 0)
               if (b == a) continue;
               // otherwise, we want to add a spring with an ideal length of the current distance
               dist = 0;
               for (j = 0; j < DIMS; ++j) 
                  dist += (a->m_x[j] - b->m_x[j])*(a->m_x[j] - b->m_x[j]);
               dist = std::sqrt(dist);
               add(new spring(a, b, dist, k));
               ++count;
            }
         }
         // that's it!
         return count;
      }

      // for getting at atoms...
      const atom* operator [] (unsigned i) const {return m_atoms[i];}
      atom* operator [] (unsigned i) {return m_atoms[i];}
      unsigned atoms() const {return m_atoms.size();}

      // this should be called before steps in order to initialize fields
      // returns the potential energy at time 0
      double start(bool verbose = false, std::ostream& err = std::cerr)
      {
         unsigned i;
         m_t = 0;
         // give the atoms indices
         for (i = 0; i < m_atoms.size(); ++i)
            m_atoms[i]->m_idx = i;
         // init the fields...
         if (verbose) err << "Initializing fields..." << endl;
         for (i = 0; i < m_fields.size(); ++i)
            m_fields[i]->init(m_atoms, m_springs);
         // okay, now that everything is initialized, we can release the threads to do their initialization
         if (m_threads > 0)
            pthread_cond_broadcast(&m_cond);
         // and return the potential
         return PE();
      }

      // take one step (of length dt) and return the kinetic energy
      // also returns the RMSD of the step by reference if requested. 
      double step(double dt, double* rmsd = 0)
      {
         unsigned i, j;
         double d = 0, ke = 0;
         m_dt = dt;
         if (m_threads == 0) {
            // update all forces...
            for (i = 0; i < m_springs.size(); ++i)
               m_springs[i]->update();
            for (i = 0; i < m_fields.size(); ++i) {
               for (j = 0; j < m_atoms.size(); ++j)
                  m_fields[i]->update(j);
            }
            // then all atoms
            for (i = 0; i < m_atoms.size(); ++i) {
               d += m_atoms[i]->update(dt);
               ke += m_atoms[i]->KE();
            }
         } else {
            pthread_mutex_lock(&m_mux);
            pthread_cond_broadcast(&m_cond);
            // wait for the threads to finish
            pthread_cond_wait(&m_done, &m_mux);
            // that's it; we can release the mutex
            d = m_result_d;
            ke = m_result_ke;
            // and unlock the mutex
            pthread_mutex_unlock(&m_mux);
         }
         // run the post-updates single-threadedly
         for (i = 0; i < m_fields.size(); ++i)
            m_fields[i]->post_update();

         m_t += dt;
         if (rmsd)
            *rmsd = sqrt(d / (double)m_atoms.size());
         return ke;
      }

      // wipe velocities of all atoms; if the strength argument is passed,
      // then it dampens or excites the velocities by that multiple.
      void wipe(double str = 0.0)
      {
         for (unsigned i = 0; i < m_atoms.size(); ++i)
            m_atoms[i]->dampen(str);
      }

      // the amount of time that has passed
      double t() const {return m_t;}

      // calculate the current kinetic energy
      double KE() const
      {
         double ke = 0, v2;
         unsigned i, j;
         for (i = 0; i < m_atoms.size(); ++i) {
            v2 = 0;
            for (j = 0; j < DIMS; ++j)
               v2 += m_atoms[i]->v()[j]*m_atoms[i]->v()[j];
            ke += m_atoms[i]->inertia() * v2;
         }
         return ke;
      }
      // calculate the potential energy; this is obtained from the fields + all the springs
      double PE() const
      {
         // this comes from all the potentials...
         double pe = 0, tmp;
         unsigned i;
         for (i = 0; i < m_springs.size(); ++i) {
            tmp = m_springs[i]->potential();
            pe += tmp;
         }
         for (i = 0; i < m_fields.size(); ++i) {
            tmp = m_fields[i]->potential();
            pe += tmp;
         }
         return pe;
      }
      double springPE() const
      {
         // this comes from all the potentials...
         double pe = 0;
         unsigned i;
         for (i = 0; i < m_springs.size(); ++i)
            pe += m_springs[i]->potential();
         return pe;
      }
      double fieldPE() const
      {
         // this comes from all the potentials...
         double pe = 0;
         unsigned i;
         for (i = 0; i < m_fields.size(); ++i)
            pe += m_fields[i]->potential();
         return pe;
      }

      // cool the system (lower PE) across the board until net KE is <= target
      double cool(double target)
      {
         if (target <= 0) {
            wipe();
            return 0;
         }
         // first, total up our KE
         double ke = KE(), lastke;
         do {
            // note the last ke
            lastke = ke;
            // if we're done, we're done!
            if (ke <= target) return ke;
            // now, lower all velocities
            wipe(target / ke);
            // note the new ke...
            ke = KE();
         } while (lastke - ke > 0.00001);
         // if we got here, there was a failure to cool...
         return ke;
      }
      // heat the system (increase PE) across the board until the target KE is reached
      double heat(double target)
      {
         // start by noting our KE
         double ke = KE();
         if (target <= ke) return ke;
         // estimate about how much we want to change the velocity by...
         double dv = 2.0 * std::sqrt(target - ke) / (double)m_atoms.size();
         unsigned i, j, k;
         double v[DIMS];
         double d;
         while (ke < target) {
            // choose two random atoms...
            i = std::rand() % m_atoms.size();
            j = std::rand() % m_atoms.size();
            if (i == j) continue;
            // pick directions for the velocities...
            d = 0;
            for (k = 0; k < DIMS; ++k) {
               v[k] = (double)std::rand() / (double)RAND_MAX - 0.5;
               d += v[k]*v[k];
            }
            if (d < 0.0000000001) continue;
            for (k = 0; k < DIMS; ++k)
               v[k] *= dv / d;
            // subtract the energy now from the KE
            ke -= m_atoms[i]->KE();
            ke -= m_atoms[j]->KE();
            // add the new velocities
            for (k = 0; k < DIMS; ++k) {
               m_atoms[i]->v()[k] += v[k];
               m_atoms[j]->v()[k] -= v[k];
            }
            // and re update the ke
            ke += m_atoms[i]->KE();
            ke += m_atoms[j]->KE();
         }
         // that's that!
         return ke;
      }
      // private function used below by minimize---do one step of gradient descent minimization
    private:
      double minimize_step(double daring)
      {
         unsigned i, j;
         double d = 0;
         m_minimizing = true;
         // now... go through and update just the forces
         if (m_threads == 0) {
            // start by wiping all velocities
            wipe();
            // update all forces...
            for (i = 0; i < m_springs.size(); ++i)
               m_springs[i]->update();
            for (i = 0; i < m_fields.size(); ++i) {
               for (j = 0; j < m_atoms.size(); ++j)
                  m_fields[i]->update(j);
            }
            for (i = 0; i < m_atoms.size(); ++i) {
               for (j = 0; j < DIMS; ++j)
                  d += m_atoms[i]->a()[j]*m_atoms[i]->a()[j];
            }
         } else {
            pthread_mutex_lock(&m_mux);
            pthread_cond_broadcast(&m_cond);
            // wait for the threads to finish
            pthread_cond_wait(&m_done, &m_mux);
            // that's it; we can release the mutex
            d = m_result_d;
         }
         // d is the total length of the gradient vector; we want to move a <daring> distance
         // along this vector
         d = std::sqrt(d);
         m_result_d = daring;
         // now, update the atoms by giving them the proper velocities
         if (m_threads == 0) {
            atom* a;
            // then all atoms
            for (i = 0; i < m_atoms.size(); ++i) {
               a = m_atoms[i];
               for (j = 0; j < DIMS; ++j) {
                  a->v()[j] = a->a()[j] * m_result_d;
                  a->a()[j] = 0;
               }
               m_atoms[i]->update(1.0);
            }
         } else {
            pthread_cond_broadcast(&m_cond);
            // wait for the threads to finish
            pthread_cond_wait(&m_done, &m_mux);
            // and unlock the mutex
            pthread_mutex_unlock(&m_mux);
         }
         // run the post-updates single-threadedly
         for (i = 0; i < m_fields.size(); ++i)
            m_fields[i]->post_update();
         // no longer minimizing...
         m_minimizing = false;
         // return the length of the gradient vector
         return d;
      }
    public:
      // minimize a system for a number of steps (do not update the internal steps counter)
      // returns the difference in pe (new - start)
      double minimize(double daring, unsigned max_steps, 
                      double min_dke = 0.001, double max_dke2 = -0.2)
      {
         double pes[3];
         unsigned k = 2;
         double grad, dk;
         double orig_pe = PE();
         // we perform however many steps, wiping the velocity each step
         for (unsigned i = 0; i < max_steps; ++i) {
            k = (k + 1) % 3;
            grad = minimize_step(daring);
            pes[k] = PE();
            if (i > 0 && pes[k] > pes[(k+2) %3]) break;
            if (grad < 0.0001) break;
            // see if we want to stop (ie, the second deriv is < min_grad2 or the deriv < min_grad
            if (i > 2) {
               dk = pes[(k+2) % 3] - pes[(k+0) % 3];
               if (dk < min_dke) break;
               dk = (pes[(k+1) % 3] - pes[(k+2) % 3]) - (pes[(k+2) % 3] - pes[(k+0) % 3]);
               if (dk > max_dke2) break;
            }
         }
         wipe();
         return orig_pe - pes[k];
      }
   };

}

#endif
