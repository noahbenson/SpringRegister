////////////////////////////////////////////////////////////////////////////////////////////////////
// params.hh
// This file defines the basic parameter class, which contains information about all the parameters
// of a spring simulation.
// By Noah C. Benson

#ifndef _____nben_PARAMS_HH
#define _____nben_PARAMS_HH

#include <map>
#include <string>
#include <sstream>

namespace springs
{
   using std::string;
   using std::map;

   // classes to be declared here
   template <typename T> class typed_param;
   class untyped_param;
   template <typename T> class param;
   class param_space;
  
   namespace params_internal
   {
      // These simple functions are used in parsing below
      template <typename T> inline bool parse_internal(const string& txt, T& val)
      {
         std::istringstream ss(txt);
         if (!(ss >> val))
            return false;
         else
            return true;
      }
      template <> inline bool parse_internal<bool>(const string& txt, bool& v)
      {
         if (txt.empty() || txt == "1" || txt == "true" || txt == "TRUE" || txt == "yes") {
            v = true;
            return true;
         } else if (txt == "0" || txt == "false" || txt == "FALSE" || txt == "no") {
            v = false;
            return true;
         } else
            return false;
      }
      template <> inline bool parse_internal<string>(const string& txt, string& s)
      {
         s = txt;
         return true;
      }
   }

   // all param's must inherit from this class, but it is basically not used except in the guts
   template <typename T>
   class typed_param
   {
    protected:
      // returns the argument
      virtual const T& get() const = 0;
      // we're friends with the untyped version...
      friend class untyped_param;
      friend class param_space;
   };

   // A parameter is basically a set of information and a value...
   // This is an abstract class for representing the value generically.
   class untyped_param
   {
    private:
      string m_name;   // e.g. "journal"
      char m_abbrev; // e.g. 'j' (-j)
      string m_doc;  // short documentation
      string m_description; // long documentation

      // not to be used...
      untyped_param();

      // we're friends with param space
      friend class param_space;

    protected:
      // this may process the argument
      virtual void process()
      {
      }
      // this may post-process any argument
      virtual void post_process()
      {
      }
      // must provide a way to parse args...
      virtual void parse_arg(const string& txt) = 0;

    public:
      untyped_param(const string& name, char abbrev, const string& doc, const string& descr):
         m_name(name),
         m_abbrev(abbrev),
         m_doc(doc),
         m_description(descr)
      {
      }
      virtual ~untyped_param()
      {
      }

      // setter/copyers
      untyped_param(const untyped_param& p):
         m_name(p.m_name),
         m_abbrev(p.m_abbrev),
         m_doc(p.m_doc),
         m_description(p.m_description)
      {
      }
      untyped_param& operator = (const untyped_param& p)
      {
         m_name = p.m_name;
         m_abbrev = p.m_abbrev;
         m_doc = p.m_doc;
         m_description = p.m_description;
         return *this;
      }

      // accessors
      const string& name() const {return m_name;}
      char abbrev() const {return m_abbrev;}
      const string& doc() const {return m_doc;}
      const string& description() const {return m_description;}

      // accessor that must be overloaded: should return the text of the value (even if default
      // is the current value)
      virtual string text() const = 0;
      // must return true/false if the paramter has been parsed/not
      virtual bool provided() const = 0;

      // get the value... throws an exception if the incorrect value type is requested
      template <typename T> const T& value() const
      {
         typed_param<T>* p = dynamic_cast<typed_param<T>*>(this);
         if (p == 0)
            throw string("Incorrect value type requested of parameter");
         else
            return p->get();
      }
   };

   template <typename T>
   class param:
      public untyped_param,
      public typed_param<T>
   {
    protected:
      T m_value;   // the value of the argument...
    private:
      string m_text; // the unparsed parsed string (if any has been given)
      bool m_provided;
      
      // not to be used this way...
      param();

      // we're friends with the param_space
      friend class param_space;

    public:
      param(const string& name, char abbrev, const string& doc, const string& descr, 
            const T& default_value):
         untyped_param(name, abbrev, doc, descr),
         m_value(default_value)
      {
         m_provided = false;
      }
      virtual ~param()
      {
      }
      param& operator = (const param& p)
      {
         untyped_param::operator = (p);
         m_value = p.m_value;
         m_text = p.m_text;
         m_provided = p.m_provided;
      }
      // access the text of the argument
      string text() const
      {
         // if there is no text yet, turn the default to a string...
         if (m_text.empty()) {
            std::ostringstream ss;
            ss << m_value;
            return ss.str();
         } else {
            return m_text;
         }
      }
      bool provided() const
      {
         return m_provided;
      }

      // Since we know the type, we can provide an accessor
      const T& operator * () const
      {
         return m_value;
      }

    protected:
      // we can overload get and parse now...
      virtual const T& get() const
      {
         return m_value;
      }
      virtual void parse_arg(const string& txt)
      {
         if (!params_internal::parse_internal(txt, m_value)) {
            std::ostringstream out;
            out << "Could not parse argument " << name() << " from text \"" << txt << "\"";
            throw out.str();
         } else {
            m_text = txt;
            m_provided = true;
            process();
         }
      }
   };

   // The parameter space is the set of all possible parameters with defaults, documentation, etc.
   class param_space
   {
    public:
      // some types...
      typedef map<string, untyped_param*> param_map_t;
      typedef map<string, untyped_param*>::iterator param_iter_t;
      typedef map<string, untyped_param*>::const_iterator param_const_iter_t;
      typedef map<char, string> abbrev_map_t;
      typedef map<char, string>::iterator abbrev_iter_t;
      typedef map<char, string>::const_iterator abbrev_const_iter_t;

    private:
      // the map of parameter names to parameter data
      param_map_t m_dict;
      abbrev_map_t m_abbrev;

      // cannot be used...
      param_space();

    public:
      // default destructor and copy operator are fine, but we need to control how this is created
      param_space(untyped_param* p, ...);

      // we can get a parameter...
      const untyped_param& operator () (const string& name) const
      {
         param_const_iter_t it = m_dict.find(name);
         if (it == m_dict.end()) {
            std::ostringstream ss;
            ss << "Requested parameter name \"" << name << "\" is not found";
            throw ss.str();
         } else {
            return *it->second;
         }
      }
      const untyped_param& operator () (const char* name) const
      {
         return operator()(string(name));
      }
      
      // parse all the parameters based on arguments from main;
      // throws an exception if any of the arguments fail to parse.
      // at the end of the function call, argc and argv will represent whatever is left in
      // the arguments.  If --help or -h is encountered the text from the help message is placed in
      // the help ostream and false is returned.
      bool parse(int& argc, char**& argv, std::ostream& help);

      // print out a help message (optional argument to follow)
      void help(std::ostream& s, const string& arg = string()) const;

      // used for outputting the parameters!
      friend std::ostream& operator << (std::ostream& o, const param_space& p);
   };

}

#endif

