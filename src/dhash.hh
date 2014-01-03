////////////////////////////////////////////////////////////////////////////////////////////////////
// dhash.hh
// Dimensional hash, dhash
// by Noah Benson

#ifndef _____nben_DHASH_HH
#define _____nben_DHASH_HH

// since the list class is a template class, we must include it here
#include <list>
#include <set>
#include <stack>
#include <cmath>

#include <iostream>
using std::endl;
using std::cerr;

// the dimensional hash structure
template <class T, unsigned dims = 3>
class dhash
{
   class cell;

public:
   class point
   {
   private:
      double* m_x;
      T m_t;
      // our current cell
      cell* m_cell;
      // next in the cell
      point* m_next;

      inline point(const T& t): m_t(t) {m_cell = 0; m_next = 0; m_x = new double[dims];}

      // not to be used
      point(const point& p);
      point& operator = (const point& p);

      // only to be used by dhash and cell
      inline ~point()
      {
         delete[] m_x;
         if (m_next) delete m_next;
      }

      friend class dhash;
      friend class dhash::cell;
   public:
      inline const T& get() const {return m_t;}
      inline T& get() {return m_t;}

      inline const double* coordinates() const {return m_x;}

      bool move(const double* newcoords);
   };

private:
   class cell
   {
   private:
      double* m_min;
      double* m_max;
      unsigned m_index;
      point* m_nodes;
      dhash* m_hash;

      cell() {}
      ~cell()
      {
         delete[] m_min;
         delete[] m_max;
         if (m_nodes) delete m_nodes;
      }

      friend class dhash;
      friend class point;
   public:

      bool move(point* p, const double* newcoords);
   };

   double* m_sizes;
   double* m_min;
   double* m_max;
   double* m_cellsz; // ideal size of a cell
   cell* m_cell;
   unsigned m_cell_count;
   unsigned* m_cells;
   unsigned* m_sz_mult;

   // find a cell given coordinates
   inline const cell* find_cell(const double* x) const
   {
      unsigned i, idx = 0;
      for (i = 0; i < dims; ++i) {
         if (x[i] < m_min[i] || x[i] >= m_max[i]) return 0;
         idx += m_sz_mult[i] * (unsigned)((x[i] - m_min[i]) / m_cellsz[i]);
      }
      return m_cell + idx;
   }
   inline cell* find_cell(const double* x)
   {
      unsigned i, idx = 0;
      for (i = 0; i < dims; ++i) {
         if (x[i] < m_min[i] || x[i] >= m_max[i]) return 0;
         idx += m_sz_mult[i] * (unsigned)((x[i] - m_min[i]) / m_cellsz[i]);
      }
      return m_cell + idx;
   }
   // find nearby cells
   inline const cell* near(const cell* c, unsigned dim, int dist = 1) const
   {
      if ((dist < 0 && c->m_min[dim] + ((double)dist + 0.5)*m_cellsz[dim] < m_min[dim]) ||
          (dist > 0 && c->m_max[dim] + ((double)dist - 0.5)*m_cellsz[dim] >= m_max[dim]))
         return 0;
      return m_cell + ((int)c->m_index + (int)m_sz_mult[dim] * dist);
   }
   inline cell* near(const cell* c, unsigned dim, int dist = 1)
   {
      if ((dist < 0 && c->m_min[dim] + ((double)dist + 0.5)*m_cellsz[dim] < m_min[dim]) ||
          (dist > 0 && c->m_max[dim] + ((double)dist - 0.5)*m_cellsz[dim] >= m_max[dim]))
         return 0;
      return m_cell + ((int)c->m_index + (int)m_sz_mult[dim] * dist);
   }

   // norms for distance
   inline double norm(const double* a, const double* b) const
   {
      double tmp = 0;
      for (unsigned i = 0; i < dims; i++) {
         tmp += (a[i] - b[i])*(a[i] - b[i]);
      }
      return std::sqrt(tmp);
   }
   inline double norm2(const double* a, const double* b) const
   {
      double tmp = 0;
      for (unsigned i = 0; i < dims; i++) {
         tmp += (a[i] - b[i])*(a[i] - b[i]);
      }
      return tmp;
   }
   
   // decide if there is any space inside a given cell that is within distance d of the given point
   inline bool in_range(const cell* c, const double* x, double d) const
   {
      double tmp = 0;
      unsigned i;
      // first step: figure out where this is relative to the hypercube and count the number of 
      // dimensions for which it's in range
      for (i = 0; i < dims; i++) {
         if (x[i] < c->m_min[i])
            tmp += (x[i] - c->m_min[i])*(x[i] - c->m_min[i]);
         else if (x[i] > c->m_max[i]) 
            tmp += (x[i] - c->m_min[i])*(x[i] - c->m_min[i]);
         // if it's between the limits of the cube along the dimension, then this dimension contributes
         // nothing to the min distance
      }
      return tmp <= d*d;
   }
   // decide if there is any space inside a given cell that is within distance d of the given point
   inline bool in_range_2(const cell* c, const double* x, double d) const
   {
      double tmp = 0;
      unsigned i;
      // first step: figure out where this is relative to the hypercube and count the number of 
      // dimensions for which it's in range
      for (i = 0; i < dims; i++) {
         if (x[i] < c->m_min[i]) {
            tmp += (x[i] - c->m_min[i])*(x[i] - c->m_min[i]);
         } else if (x[i] > c->m_max[i]) {
            tmp += (x[i] - c->m_max[i])*(x[i] - c->m_max[i]);
         }
         // if it's between the limits of the cube along the dimension, then this dimension contributes
         // nothing to the min distance
      }
      return tmp <= d;
   }

public:
   dhash(const double* min, const double* max, const double* ideal = 0)
   {
      unsigned i, j, k;
      double tmp;

      m_min = new double[dims];
      m_max = new double[dims];
      memcpy(m_min, min, dims * sizeof(double));
      memcpy(m_max, max, dims * sizeof(double));
      m_sizes = new double[dims];
      for (i = 0; i < dims; i++) {
         m_sizes[i] = max[i] - min[i];
      }

      m_cellsz = new double[dims];
      if (ideal == 0) {
         for (i = 0; i < dims; i++)
            m_cellsz[i] = std::log(m_sizes[i]);
      } else
         memcpy(m_cellsz, ideal, sizeof(double) * dims);

      // for each dimension, we want a number as close as possible to the ideal that evenly
      // divides the space into chunks
      k = 1;
      m_cells = new unsigned[dims];
      for (i = 0; i < dims; i++) {
         // if the ideal divides the hash size, then we want the ideal to be the cell size
         // if the mod of the hash size by the ideal is close to the ideal, then we want 
         //   to round up the ideal to a divisor of the hash size
         // if the mod of the hash size by the ideal is close to the 0, then we want 
         //   to round down the ideal to a divisor of the hash size
         tmp = std::floor(m_sizes[i] / m_cellsz[i] + 0.5);
         m_cellsz[i] = m_sizes[i] / tmp;
         m_cells[i] = (unsigned)tmp;
         k *= (unsigned)tmp;
      }
      // k is now the size of the cells array
      m_cell = new cell[k];
      memset(m_cell, 0, sizeof(cell) * k);
      m_cell_count = k;

      // initialize all the cells... we need a temporary for this
      double* pos = new double[dims];
      memset(pos, 0, sizeof(double) * dims);
      for (i = 0; i < k; i++) {
         m_cell[i].m_hash = this;
         m_cell[i].m_nodes = 0;
         m_cell[i].m_index = i;
         m_cell[i].m_min = new double[dims];
         m_cell[i].m_max = new double[dims];
         for (j = 0; j < dims; j++) {
            m_cell[i].m_min[j] = pos[j] + min[j];
            m_cell[i].m_max[j] = pos[j] + min[j] + m_cellsz[j];
         }
         // this will update pos so that we count along as we go.
         for (j = 1; j <= dims; j++) {
            pos[dims - j] += m_cellsz[dims - j];
            if (pos[dims - j] + m_cellsz[dims - j]/2.0 < m_sizes[dims - j])
               break;
            else
               pos[dims - j] = 0;
         }
      }

      // next, fill out the size multipliers
      m_sz_mult = new unsigned[dims];
      k = 1;
      for (i = 1; i <= dims; i++) {
         m_sz_mult[dims - i] = k;
         k *= m_cells[dims - i];
      }
      // that should be it!

      // cleanup
      delete[] pos;
   }

   ~dhash()
   {
      delete[] m_sizes;
      delete[] m_cellsz;
      delete[] m_min;
      delete[] m_max;
      delete[] m_cells;
      delete[] m_sz_mult;
   }

   inline const double* min() const {return m_min;}
   inline const double* max() const {return m_max;}

   inline point* insert(const T& t, const double* x)
   {
      cell* c = find_cell(x);
      if (!c) return 0;
      // add the point...
      point* p = new point(t);
      p->m_cell = c;
      p->m_next = c->m_nodes;
      c->m_nodes = p;
      memcpy(p->m_x, x, sizeof(double) * dims);
      return p;
   }

   // Find a point exactly
   const point* find(const double* x) const
   {
      const cell* c = find_cell(x);
      if (!c) return 0;
      const point* p;
      unsigned i;
      for (p = c->m_nodes; p; p = p->m_next) {
         for (i = 0; i < dims; i++) {
            if (p->m_x[i] != x[i]) break;
         }
         if (i == dims) return p;
      }
      return 0;
   }
   point* find(const double* x)
   {
      cell* c = find_cell(x);
      if (!c) return 0;
      point* p;
      unsigned i;
      for (p = c->m_nodes; p; p = p->m_next) {
         for (i = 0; i < dims; i++) {
            if (p->m_x[i] != x[i]) break;
         }
         if (i == dims) return p;
      }
      return 0;
   }

   // find all points within a certain distance
   bool find(const double* x, double d, std::list<const point*>& l) const
   {
      unsigned i;
      const point *p;
      const cell* c;
      // we can use squares and norm2 to avoid sqrt function...
      d *= d;
      // first, find the base point we're operating from
      const cell* base = find_cell(x);
      if (!base) return false;
      // okay, we need to travel in a certain distance in each direction...
      // there's a -1 and +1 direction for each dim...

      // make sure we have a unique/safe mark...
      std::set<const cell*> mark;
      // push the base onto the search stack...
      std::stack<const cell*> workc;
      workc.push(base);
      mark.insert(base);
      // we start the search here...
      while (!workc.empty()) {
         // pop off the next cell to search
         base = workc.top();
         workc.pop();
         for (p = base->m_nodes; p; p = p->m_next) {
            if (norm2(p->m_x, x) < d)
               l.push_front(p);
         }
         // okay, scanned this one... now push others onto the stack
         for (i = 0; i < dims; i++) {
            // positive direction...
            c = near(base, i, 1);
            if (c && mark.find(c) == mark.end() && in_range_2(c, x, d)) {
               // mark it!
               mark.insert(c);
               // and push it onto the stack
               workc.push(c);
            }
            // negative direction...
            c = near(base, i, -1);
            if (c && mark.find(c) == mark.end() && in_range_2(c, x, d)) {
               // mark it!
               mark.insert(c);
               // and push it onto the stack
               workc.push(c);
            }
         }
         // that's it!
      }
      return true;
   }
   // same function but non-const
   // find all points within a certain distance
   bool find(const double* x, double d, std::list<point*>& l)
   {
      unsigned i;
      point *p;
      cell* c;
      // we can use squares and norm2 to avoid sqrt function...
      d *= d;
      // first, find the base point we're operating from
      cell* base = find_cell(x);
      if (!base) return false;
      // okay, we need to travel in a certain distance in each direction...
      // there's a -1 and +1 direction for each dim...

      // make sure we have a unique/safe mark...
      std::set<cell*> mark;
      // push the base onto the search stack...
      std::stack<cell*> workc;
      workc.push(base);
      mark.insert(base);
      // we start the search here...
      while (!workc.empty()) {
         // pop off the next cell to search
         base = workc.top();
         workc.pop();
         for (p = base->m_nodes; p; p = p->m_next) {
            if (norm2(p->m_x, x) < d)
               l.push_front(p);
         }
         // okay, scanned this one... now push others onto the stack
         for (i = 0; i < dims; i++) {
            // positive direction...
            c = near(base, i, 1);
            if (c && mark.find(c) == mark.end() && in_range_2(c, x, d)) {
               // mark it!
               mark.insert(c);
               // and push it onto the stack
               workc.push(c);
            }
            // negative direction...
            c = near(base, i, -1);
            if (c && mark.find(c) == mark.end() && in_range_2(c, x, d)) {
               // mark it!
               mark.insert(c);
               // and push it onto the stack
               workc.push(c);
            }
         }
         // that's it!
      }
      return true;
   }

   // remove a point
   bool remove(const point* cp)
   {
      point* p = const_cast<point*>(cp);
      cell* c = p->m_cell;
      if (c->m_dhash != this) return false;
      if (c->m_nodes == p) {
         c->m_nodes = p->m_next;
         delete p;
      } else if (!c->m_nodes) {
         return false;
      } else {
         point* tmp;
         for (tmp = c->m_nodes; tmp->m_next && tmp->m_next != p; tmp = tmp->m_next)
            ;
         if (!tmp->m_next) return false;
         tmp->m_next = p->m_next;
         delete p;
      }
      return true;
   }

private:
   // private move function
   bool move(cell* c, point* p, const double* x)
   {
      // first thing to check--can it be moved?
      cell* tmp = find_cell(x);
      if (!tmp) return false;
      // next... is it being moved within this cell?
      if (c == tmp) {
         memcpy(p->m_x, x, sizeof(double) * dims);
         return true;
      }
      // basically, just extract it then re-add it
      if (c->m_nodes == p) {
         c->m_nodes = p->m_next;
      } else if (!c->m_nodes) {
         return false;
      } else {
         point* q;
         for (q = c->m_nodes; q->m_next && q->m_next != p; q = q->m_next)
            ;
         if (!q->m_next) return false;
         q->m_next = p->m_next;
      }
      // okay, point removed... now add it...
      c = tmp;
      // add the point...
      p->m_cell = c;
      p->m_next = c->m_nodes;
      c->m_nodes = p;
      memcpy(p->m_x, x, sizeof(double) * dims);
      return true;
   }
};

// some inline funcs that were out of order...
template <class T, unsigned dims> inline bool dhash<T, dims>::point::move(const double* newcoords) 
{
   return m_cell->move(this, newcoords);
}
template <class T, unsigned dims> inline bool dhash<T, dims>::cell::move(point* p, const double* newcoords) 
{
   return m_hash->move(this, p, newcoords);
}

#endif
