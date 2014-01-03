////////////////////////////////////////////////////////////////////////////////////////////////////
// vdw.hh
// Specification of the van-der-waals forces
// by Noah Benson

#ifndef _____nben_VDW_HH
#define _____nben_VDW_HH

#include <cmath>
#include <vector>
#include "springs.hh"
#include "dhash.hh"

template <unsigned DIMS>
class simple_vdw:
   public springs::system<DIMS>::field
{
 public:
   // handy aliases
   typedef typename springs::system<DIMS>::atom atom_t;
   typedef typename springs::system<DIMS>::spring spring_t;
   typedef typename springs::system<DIMS> system_t;
   typedef          dhash<atom_t*, 2> dhash_t;
   typedef typename dhash<atom_t*, 2>::point point_t;
 private:
   // the atoms that matter charge-wise
   std::vector<atom_t*> m_atoms;
   // the points from the dhash get stored here...
   std::vector<point_t*> m_points;
   // for keeping track of which pairs of atoms have springs (thus are ignored by vdw)
   unsigned char* m_springs;

   // the cutoff distnace
   double m_cutoff;
   // the max VDW force
   double m_max;
   // the integration const for the potential
   double m_pot_const;
   // the integration consts for the edge potential...
   double m_epot_const[DIMS];

   // the distance hash that keeps track of everything
   dhash_t* m_hash;
   // the edge distance force cutoff
   double m_edge[DIMS];
   double m_mn[DIMS];
   double m_mx[DIMS];

   // the force and the potential...
   // the normal versions of these assume dist < cutoff
   // the edge versions do not
   inline double force(double dist) const {return (dist > m_cutoff? 0 : m_max * (2.0 / (dist / m_cutoff + 1.0) - 1.0));}
   inline double potential(double dist) const {
      return (dist > m_cutoff? 0 : -m_max * (2.0 * m_cutoff * std::log(dist + m_cutoff) - (dist + m_cutoff)) + m_pot_const);
   }
   // edge forces
   inline double edge_force(double x, unsigned i) const {
      if (x >= m_mx[i]) return -m_max;
      else if (x <= m_mn[i]) return m_max;
      else if (x > m_mx[i] - m_edge[i]) return -(m_max * (2.0 / (std::fabs(m_mx[i]-x) / m_edge[i] + 1.0) - 1.0));
      else if (x < m_mn[i] + m_edge[i]) return  (m_max * (2.0 / (std::fabs(m_mn[i]-x) / m_edge[i] + 1.0) - 1.0));
      else return 0;
   }
   inline double edge_potential(double x, unsigned i) const {
      if (x >= m_mx[i] || x <= m_mn[i]) return m_epot_const[i];
      else if (x > m_mx[i] - m_edge[i]) 
         return -m_max * (2.0 * m_edge[i] * std::log(std::fabs(m_mx[i]-x) + m_edge[i]) - (std::fabs(m_mx[i]-x) + m_edge[i])) + m_epot_const[i];
      else if (x < m_mn[i] + m_edge[i])
         return -m_max * (2.0 * m_edge[i] * std::log(std::fabs(m_mn[i]-x) + m_edge[i]) - (std::fabs(m_mn[i]-x) + m_edge[i])) + m_epot_const[i];
      else return 0;
   }
 public:
   simple_vdw(double cutoff, double max = 2.0):
      m_cutoff(cutoff),
      m_max(max)
   {
      m_hash = 0;
      m_springs = 0;
   }
   ~simple_vdw()
   {
      if (m_hash) delete m_hash;
      if (m_springs) delete m_springs;
   }

   // init must save the vectors...
   virtual void init(const std::vector<atom_t*>& all_atoms, const std::vector<spring_t*>& all_springs)
   {
      unsigned i, j;
      double* c;
      m_atoms = all_atoms;
      // if cutoff was <= 0, we guess it as half the average spring length...
      if (m_cutoff <= 0) {
         double dmu = 0;
         for (i = 0; i < all_springs.size(); ++i)
            dmu += all_springs[i]->length();
         dmu /= (double)all_springs.size();
         m_cutoff = dmu / 1.5;
      }
      // we now have both a max and a cutoff for sure; we can set the integration const
      m_pot_const = m_max * (2.0 * m_cutoff * (std::log(2.0 * m_cutoff) - 1.0));
      // next, we want to put all the atoms in the distance hash
      // to do this we first need to know the min/max size of the system
      for (i = 0; i < DIMS; ++i) {
         m_mn[i] = m_atoms[0]->x()[i];
         m_mx[i] = m_atoms[0]->x()[i];
      }
      for (i = 1; i < m_atoms.size(); ++i) {
         c = m_atoms[i]->x();
         for (j = 0; j < DIMS; ++j) {
            if (c[j] < m_mn[j]) m_mn[j] = c[j];
            if (c[j] > m_mx[j]) m_mx[j] = c[j];
         }
      }
      // okay, we have min/max; we want to add some to those
      for (i = 0; i < DIMS; ++i) {
         m_edge[i] = 0.1 * (m_mx[i] - m_mn[i]);
         m_mn[i] -= m_edge[i];
         m_mx[i] += m_edge[i];
      }
      // now create the hash...
      m_hash = new dhash_t(m_mn, m_mx);
      // add all the atoms to the hash...
      m_points.resize(m_atoms.size(), 0);
      for (i = 0; i < m_atoms.size(); ++i)
         m_points[i] = m_hash->insert(m_atoms[i], m_atoms[i]->x());
      // before we're done, we still need to set the constants of integration for the edge forces
      for (j = 0; j < DIMS; ++j)
         m_epot_const[j] = m_max * (2.0*m_edge[j] * (std::log(2.0*m_edge[j]) - 1.0));

      // okay, now, we want to make note of all the atom pairs with springs
      m_springs = new unsigned char[m_atoms.size() * m_atoms.size()];
      std::memset(m_springs, 0, sizeof(unsigned char) * m_atoms.size() * m_atoms.size());
      // and set them to true where there are springs
      for (i = 0; i < all_springs.size(); ++i) {
         m_springs[all_springs[i]->a()->index() * m_atoms.size() + all_springs[i]->b()->index()] = 1;
         m_springs[all_springs[i]->b()->index() * m_atoms.size() + all_springs[i]->a()->index()] = 1;
      }
      // that's all we need to do!
   }

   // update the acceleration for an atom...
   virtual void update(unsigned i)
   {
      unsigned j;
      const double* c;
      double* x;
      double vec[DIMS];
      double tmp;
      atom_t* a, *b;
      // we go through and do van-der-waals forces...
      std::list<point_t*> l;
      // we are given atom i
      a = m_atoms[i];
      x = a->x();
      // first, check to see if we need to add any edge repulsion forces
      for (j = 0; j < DIMS; ++j) {
         if (x[j] > m_mx[j]) x[j] = m_mx[j] - 0.0000001;
         if (x[j] < m_mn[j]) x[j] = m_mn[j] + 0.0000001;
         a->a()[j] += edge_force(x[j], j);
      }
      m_hash->find(x, m_cutoff, l);
      // we do the forces that push on a only...
      while (!l.empty()) {
         b = l.front()->get();
         l.pop_front();
         // ignore it if they're the same or they are spring-connected
         if (b == a || m_springs[a->index() * m_atoms.size() + b->index()])
            continue;
         c = b->x();
         // get the vector from the other atom to a and it's length
         tmp = 0;
         for (j = 0; j < DIMS; ++j) {
            vec[j] = x[j] - c[j];
            tmp += vec[j]*vec[j];
         }
         tmp = std::sqrt(tmp);
         // we'll want to add a vector in the direction of vec to the accel
         x = a->a();
         tmp = force(tmp) / tmp; // /tmp to normalize vec below; *force(tmp) to create the force
         // here, we actually update a's vector
         for (j = 0; j < DIMS; ++j)
            x[j] += vec[j] * tmp;
      }
      // that's it!
   }

   virtual void post_update()
   {
      // We go through and update the position of each atom in the hash.
      // Note that this is run single-threadedly no matter what
      for (unsigned i = 0; i < m_atoms.size(); ++i) {
         // move the point to wherever it is now...
         m_points[i]->move(m_atoms[i]->x());
      }
   }

   virtual double potential() const
   {
      double tot = 0, tmp;
      double* x;
      const double *c;
      unsigned i, j;
      std::list<point_t*> l;
      atom_t* a, *b;
      for (i = 0; i < m_atoms.size(); ++i) {
         a = m_atoms[i];
         x = a->x();
         m_hash->find(x, m_cutoff, l);
         // we do the forces that push on a only...
         while (!l.empty()) {
            b = l.front()->get();
            c = l.front()->coordinates();
            l.pop_front();
            if (a == b || m_springs[i * m_atoms.size() + b->index()])
               continue;
            // get the distance between the two...
            tmp = 0;
            for (j = 0; j < DIMS; ++j)
               tmp += (x[j] - c[j])*(x[j] - c[j]);
            tot += potential(std::sqrt(tmp));
            // we'll want to add a vector in the direction of vec to the accel
         }
         // we must also consider edge potential
         for (j = 0; j < DIMS; ++j)
            tot += edge_potential(x[j], j);
      }
      return tot;
   }
};

#endif
