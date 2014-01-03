////////////////////////////////////////////////////////////////////////////////////////////////////
// angle_springs.hh
// Class for handling angles between atoms and keeping them from crossing
// by Noah Benson

#ifndef _____nben_ANGLE_SPRINGS_HH
#define _____nben_ANGLE_SPRINGS_HH

#include <cmath>
#include <algorithm>
#include <vector>
#include <map>
#include "springs.hh"

namespace springs
{
   using std::sqrt;
   using std::atan2;
   using std::exp;
   using std::fabs;
   using std::vector;
   
   class angle_field:
      public system<2>::field
   {
    public:
      typedef system<2>::atom atom_t;
      typedef system<2>::spring spring_t;

      // returns the counterclockwise angle from a to b with center as the center
      inline static double ccw_angle(atom_t* center, atom_t* a, atom_t* b)
      {
         double aa = atan2(a->x()[1] - center->x()[1], a->x()[0] - center->x()[0]);
         double ba = atan2(b->x()[1] - center->x()[1], b->x()[0] - center->x()[0]);
         if (aa < 0) aa += 2.0 * M_PI;
         if (ba < 0) ba += 2.0 * M_PI;
         if (aa > ba)
            return 2.0 * M_PI - (aa - ba);
         else
            return ba - aa;
      }

      class angle_bond
      {
         // the two springs
         atom_t* m_center; // the spring at the center
         atom_t* m_ccw; // the counterclockwise spring
         atom_t* m_cw; // clockwise spring
         double m_a0; // the ideal angle from cw to ccw

       public:
         //new angle_bond(a, sprs[(i+1) % angs.size()], sprs[i], angs[(i+1) % angs.size()] - angs[i]));
         angle_bond(atom_t* center, atom_t* ccw, atom_t* cw, double a0):
            m_center(center), m_ccw(ccw), m_cw(cw), m_a0(a0)
         {
         }
         
         inline const atom_t* cw() const {return m_cw;}
         inline const atom_t* ccw() const {return m_ccw;}
         inline const atom_t* center() const {return m_center;}
         inline atom_t* cw() {return m_cw;}
         inline atom_t* ccw() {return m_ccw;}
         inline atom_t* center() {return m_center;}
         inline double a0() const {return m_a0;}

         // the current angle between the atoms
         double angle()
         {
            return angle_field::ccw_angle(m_center, m_cw, m_ccw);
         }

         friend class angle_field;
      };

    private:
      // the formula for deciding how strong spring interactions are:
      // given two springs, cw and ccw, and an ideal angle a0, a strength k, a gain q, and a
      // minimum operable distance x0
      // the potential is equal to:
      // exp(-k * a) - exp(-k * x0)
      // This goes all the way around, so for an atom with only 2 springs, there are 2 angles.
      // The derivative is always with the springs pushing away from each other
      double m_strength;
      double m_min_dist;
      // the list of all bonds
      vector<angle_bond*> m_angles;
      // a map of atom id to a list of angle bonds
      vector< vector<angle_bond*> > m_angles_by_atom;
      // and the list of atoms
      vector<atom_t*> m_atoms;

      // used when sorting things in the init function
      class sorter
      {
       public:
         vector<int> idx;
         vector<atom_t*> a;
         vector<double> v;

         inline bool operator () (int i, int j)
         {
            if (v[i] < v[j]) return true;
            else return false;
         }

         sorter(vector<atom_t*>& _a, vector<double>& _v):
            a(_a), v(_v)
         {
            unsigned i;
            for (i = 0; i < a.size(); ++i)
               idx.push_back(i);
            // sort them
            std::sort(idx.begin(), idx.end(), *this);
            // fix _a and _v based on the index
            for (i = 0; i < idx.size(); ++i) {
               _v[i] = v[idx[i]];
               _a[i] = a[idx[i]];
            }
         }
      };

    public:
      angle_field(double strength = 1.0, double min_dist = M_PI / 6)
      {
         m_strength = strength;
         m_min_dist = min_dist;
      }

      virtual void init(const std::vector<atom_t*>& all_atoms, 
                        const std::vector<spring_t*>& all_springs)
      {
         unsigned i, j;
         atom_t* a;
         vector<atom_t*> sprs;
         vector<angle_bond*> atom_angles;
         const spring_t* spr;
         vector<double> angs;
         std::list<spring_t*>::const_iterator it;
         angle_bond* ab;
         m_angles.clear();
         m_angles_by_atom.assign(all_atoms.size(), vector<angle_bond*>());
         m_atoms = all_atoms;

         // get all of the angle bonds now...
         for (i = 0; i < all_atoms.size(); ++i) {
            a = all_atoms[i];
            sprs.clear();
            angs.clear();
            for (it = a->springs().begin(); it != a->springs().end(); it++) {
               spr = *it;
               if (spr->a() == a) {
                  sprs.push_back(const_cast<atom_t*>(spr->b()));
                  // put the angle in too
                  angs.push_back(atan2(spr->b()->x()[1] - a->x()[1], 
                                       spr->b()->x()[0] - a->x()[0]));
               } else if (spr->b() == a) {
                  sprs.push_back(const_cast<atom_t*>(spr->a()));
                  angs.push_back(atan2(spr->a()->x()[1] - a->x()[1], 
                                       spr->a()->x()[0] - a->x()[0]));
               }
            }
            // okay, we have collected all springs for this atom...
            if (sprs.size() <= 2) continue;
            // sort them by angle
            sorter(sprs, angs);
            // now, go around and make angle springs for them
            atom_angles.clear();
            for (j = 0; j < angs.size(); ++j) {
               ab = new angle_bond(a, sprs[(j+1) % angs.size()], sprs[j],
                                   ccw_angle(a, sprs[j], sprs[(j+1) % angs.size()]));
               m_angles.push_back(ab);
               m_angles_by_atom[sprs[(j+1) % angs.size()]->index()].push_back(ab);
               m_angles_by_atom[sprs[j]->index()].push_back(ab);
            }
            // that's it for this atom
         }
         // and that's it for the initialization         
      }

      virtual void update(unsigned i)
      {
         // updates a single atom...
         double cwgrad[2], ccwgrad[2];
         atom_t* a = m_atoms[i];
         vector<angle_bond*>& v = m_angles_by_atom[i];
         angle_bond* ab;
         for (unsigned j = 0; j < v.size(); ++j) {
            ab = v[j];
            update(ab, cwgrad, ccwgrad);
            if (ab->m_cw == a) {
               a->a()[0] += cwgrad[0];
               a->a()[1] += cwgrad[1];
            } else {
               a->a()[0] += ccwgrad[0];
               a->a()[1] += ccwgrad[1];
            }
         }
      }

      virtual void finish()
      {
         // nothing to do here
      }

      virtual double potential() const
      {
         double sum = 0;
         for (unsigned i = 0; i < m_angles.size(); ++i) {
            sum += potential(m_angles[i]);
         }
         return 2.0 * sum;
      }

      double potential(angle_bond* ang, double* cwgrad = 0, double* ccwgrad = 0) const
      {
         double a = ang->angle();
         double x[2], xh[2], tmp;
         // the gradients are the vectors tangent to the circle pointing away from the other atom
         if (cwgrad) {
            // this gradient goes clockwise (negative direction)
            x[0] = ang->m_cw->x()[0] - ang->m_center->x()[0];
            x[1] = ang->m_cw->x()[1] - ang->m_center->x()[1];
            tmp = sqrt(x[0]*x[0] + x[1]*x[1]);
            xh[0] = x[0] / tmp;
            xh[1] = x[1] / tmp;
            // xh cross [0 0 -1] gives us the cw direction tangent to the cw point
            cwgrad[0] = -xh[1];
            cwgrad[1] = xh[0];
         }
         if (ccwgrad) {
            // this gradient goes counterclockwise (positive direction)
            x[0] = ang->m_ccw->x()[0] - ang->m_center->x()[0];
            x[1] = ang->m_ccw->x()[1] - ang->m_center->x()[1];
            tmp = sqrt(x[0]*x[0] + x[1]*x[1]);
            xh[0] = x[0] / tmp;
            xh[1] = x[1] / tmp;
            // xh cross [0 0 1] gives us the ccw direction tangent to the ccw point
            cwgrad[0] = xh[1];
            cwgrad[1] = -xh[0];
         }
         if (a > m_min_dist) return 0;
         return exp(-m_strength * a) - exp(-m_strength * m_min_dist);
      }

      double update(angle_bond* ang, double* cwgrad, double* ccwgrad) const
      {
         // start by getting the potential!
         double pot;
         pot = potential(ang, cwgrad, ccwgrad);
         // if pot is 0, we don't worry about this
         if (pot <= 0) return 0;
         // otherwise, we need to do some updating
         cwgrad[0] *= pot * m_strength;
         cwgrad[1] *= pot * m_strength;
         ccwgrad[0] *= pot * m_strength;
         ccwgrad[1] *= pot * m_strength;
         return pot;
      }

   };
}

#endif
