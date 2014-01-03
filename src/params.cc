////////////////////////////////////////////////////////////////////////////////////////////////////
// params.cc
// Implements the parameter classes defined in params.hh
// by Noah C. Benson

#include <cstdarg>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>

#include "params.hh"

namespace springs
{
   using namespace std;

   param_space::param_space(untyped_param* p, ...)
   {
      va_list v;
      va_start(v, p);
      // go through and match each...
      while (p) {
         // if we have this name or abbrev already, throw an error...
         if (m_dict.find(p->name()) != m_dict.end())
            throw (string("Argument with name ") + p->name() + " duplicated!");
         else if (m_abbrev.find(p->abbrev()) != m_abbrev.end())
            throw (string("Argument with abbrev ") + p->abbrev() + " duplicated!");
         // we add it to the dictionary and the abbrev map
         m_dict.insert(make_pair(p->name(), p));
         m_abbrev.insert(make_pair(p->abbrev(), p->name()));
         // get the next one
         p = va_arg(v, untyped_param*);
      }
      va_end(v);
   }

   // private little function for word-wrapping args
   static vector<string> wordwrap(string s, unsigned cols)
   {
      unsigned cur = cols;
      string::size_type n;
      vector<string> res;
      while (!s.empty()) {
         if (s.length() <= cols) {
            res.push_back(s);
            return res;
         }
         n = s.rfind(' ', cur);
         if (n == string::npos)
            n = s.find(' ', cur);
         if (n != string::npos) {
            res.push_back(s.substr(0, n-1));
            s = s.substr(n+1);
         } else {
            res.push_back(s);
            s.clear();
         }
      }
      return res;
   }

   void param_space::help(ostream& o, const string& s) const
   {
      param_const_iter_t it;
      unsigned maxlen = 0;
      ostringstream ss;
      if (!s.empty()) {
         // display a specific help message
         it = m_dict.find(s);
         if (it == m_dict.end())
            o << "No such option: " << s << endl;
         else if (it->second->doc().empty())
            o << "No help available on option " << s << endl;
         else
            o << it->second->doc() << endl;
      } else {
         // start with some basic message
         o << "Options: " << endl;
         for (it = m_dict.begin(); it != m_dict.end(); ++it)
            if (it->first.length() > maxlen)
               maxlen = it->first.length();
         const char* columns = getenv("COLUMNS");
         unsigned cols;
         if (!columns || !(istringstream(columns) >> cols))
            cols = 80;
         if (cols < maxlen + 12 + 20) {
            cols = 100000000;
         }
         // first, display this message...
         o << "  --" << left << setw(maxlen) << "help" << " | -h   Display this message." << endl;
         for (it = m_dict.begin(); it != m_dict.end(); ++it) {
            o << "  --" << left << setw(maxlen) << it->first << " | -" << it->second->abbrev() << "   ";
            vector<string> v = wordwrap(it->second->doc(), cols - maxlen - 12);
            // print out the help text
            for (unsigned i = 0; i < v.size(); ++i) {
               if (i > 0)
                  o << setw(maxlen + 12) << " ";
               o << v[i] << endl;
            }
         }
      }
   }

   bool param_space::parse(int& argc, char**& argv, ostream& help_out)
   {
      vector<char*> unparsed;
      string s, arg;
      string::size_type n;
      untyped_param* p;
      param<bool>* bp = 0;
      abbrev_iter_t abbrev_it;
      param_iter_t param_it;
      // we basically want to walk through each argument...
      for (--argc, ++argv; argc > 0; --argc, ++argv) {
         s = *argv;
         // we have an argument... first, sanity check...
         if (s.length() < 2) {
            unparsed.push_back(*argv);
            continue;
         } else if (s.find("--help") == 0) {
            // this means we give help text and return false
            if (s.length() > 6 && s[6] == '=')
               help(help_out, s.substr(7));
            else
               help(help_out);
            return false;
         } else if (s.find("-h") == 0) {
            if (s.length() > 2)
               help(help_out, s.substr(2));
            else
               help(help_out);
            return false;
         } else if (s[0] == '-' && s[1] == '-') {
            // starts with -- ...
            if (s == "--") {
               // this is special--no more arg parsing...
               while (argc > 1) {
                  --argc; ++argv;
                  unparsed.push_back(*argv);
               }
               break;
            } else {
               // we trim off the front and parse out the arg
               s = s.substr(2);
               n = s.find('=');
               if (n == string::npos) {
                  // this is fine as long as the argument is a boolean flag; if not, we must
                  // assume that the next arg is the followup
                  arg.clear();
               } else {
                  arg = s.substr(n+1);
                  s = s.substr(0, n);
               }
            }
         } else if (s[0] == '-') {
            // starts with - ...
            // there are a couple possibilities here... boolean argument are allowed to be 
            // grouped (-abc instead of -a -b -c), but others are allowed to have their args
            // immediately follow the abbrev (-abc instead of -a bc)
            unsigned i = 1;
            // run through as many boolean values as are in this chunk...
            bool last_boolean = false;
            do {
               abbrev_it = m_abbrev.find(s[i]);
               if (abbrev_it == m_abbrev.end()) {
                  // not found...
                  unparsed.push_back(*argv);
                  last_boolean = true;
                  break;
               }
               // now we know the argument name... see if it's a boolean
               param_it = m_dict.find(abbrev_it->second);
               bp = dynamic_cast<param<bool>*>(param_it->second);
               if (bp != 0) {
                  // it's a boolean; we can parse it and continue...
                  last_boolean = true;
                  bp->m_value = true;
                  bp->process();
               } else {
                  last_boolean = false;
                  arg = s.substr(i+1);
                  s = abbrev_it->second;
               }
            } while (last_boolean && ++i < s.length());
            // indicates that we broke because of a non-found arg
            if (last_boolean) continue;
            // we may have processed all the boolean flags already...
            if (i == s.length()) continue;
         } else {
            unparsed.push_back(*argv);
            continue;
         }
         // okay, at this point s is the name, arg is the argument; first look it up...
         param_it = m_dict.find(s);
         if (param_it == m_dict.end()) {
            // argument not found... we don't parse this one
            unparsed.push_back(*argv);
            continue;
         }
         // alright, we have a param object...
         p = param_it->second;
         // let's see if it's a boolean arg (does not need argument)
         bp = dynamic_cast<param<bool>*>(p);
         // now, if we weren't given an arg...
         if (arg.empty()) {
            if (bp != 0) {
               // ...and this is a boolean, fine.
               bp->m_value = true;
               continue;
            } else {
               // error...
               std::ostringstream ss;
               ss << "Error parsing argument: " << p->name() << " from " << *argv 
                  << ": arguments values must not be space-separated from names!" << endl;
               throw ss.str();
            }
         }
         // try parsing...
         try {
            p->parse_arg(arg);
            // if we didn't catch an error, great!
         } catch (string err) {
            // there was an error parsing... we rethrow it with some padding
            std::ostringstream ss;
            ss << "Error parsing argument: " << s << " from " << *argv << ": " << endl
               << '\t' << err << endl;
            throw ss.str();
         }
      }
      // at this point, we've exhausted all args... put them back so main can continue using the
      // unparsed arguments
      argc = unparsed.size();
      argv = new char*[argc + 1];
      for (n = 0; n < (unsigned)argc; ++n)
         argv[n] = unparsed[n];
      argv[n] = 0;
      // finally, before returning, we post-process all arguments
      for (param_it = m_dict.begin(); param_it != m_dict.end(); ++param_it)
         param_it->second->post_process();
      return true;
   }

   std::ostream& operator << (std::ostream& o, const param_space& p)
   {
      // we go through one at a time and output them
      param_space::param_const_iter_t it;
      param<string>* sp;
      bool first = true;
      for (it = p.m_dict.begin(); it != p.m_dict.end(); ++it) {
         o << (first? "{" : ",\n ") << "{\"" << it->first << "\" -> {";
         sp = dynamic_cast<param<string>*>(it->second);
         if (sp) o << '"';
         o << it->second->text();
         if (sp) o << '"';
         o << ", " << (it->second->provided()? "True" : "False") << "}}";
         first = false;
      }
      o << "}\n";
      return o;
   }

}
