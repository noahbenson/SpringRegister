////////////////////////////////////////////////////////////////////////////////////////////////////
// coord_spring.hh
// Specification of a spring that just pulls something to a particular position
// by Noah Benson

#ifndef _____nben_COORD_SPRING_HH
#define _____nben_COORD_SPRING_HH

#include <cstdlib>
#include <cmath>
#include <vector>
#include "springs.hh"

template <unsigned DIMS>
class coord_spring:
   public springs::system<DIMS>::field
{
 public:
   // handy aliases
   typedef typename springs::system<DIMS>::atom atom_t;
   typedef typename springs::system<DIMS>::spring spring_t;
   typedef typename springs::system<DIMS> system_t;
 private:
   // the atom we are holding in place
   atom_t* m_atom;
   unsigned m_id;
   // the position at which we're holding it...
   double m_origin[DIMS];

 public:
   coord_spring(unsigned atomid, const double* coord)
   {
      m_id = atomid;
      std::memcpy(m_origin, coord, DIMS*sizeof(double));
   }
   ~coord_spring()
   {
   }
  
   // init--this does nothing...
   virtual void init(const std::vector<atom_t*>& all_atoms, const std::vector<spring_t*>& all_springs)
   {
      m_atom = (atom_t*)all_atoms.at(m_id);
   }

   // also does nothing...
   virtual void update(unsigned id)
   {
   }

   // single-threadedly replace the atom at the end of the step
   virtual void post_update()
   {
      // deal with this one atom here...
      for (unsigned i = 0; i < DIMS; ++i)
         m_atom->x()[i] = m_origin[i];
   }

   // the potential for this is not calculated...
   virtual double potential() const
   {
      return 0;
   }
};
#endif
