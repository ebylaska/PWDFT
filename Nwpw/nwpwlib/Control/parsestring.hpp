#ifndef _PARSESTRING_H_
#define _PARSESTRING_H_
/* parsestring.hpp
   Author - Eric Bylaska
*/


#include	<vector>
#include 	<iostream>
#include	<fstream>
#include        <streambuf>
#include	<sstream>
#include	<string>
#include 	<algorithm>

using namespace std;

namespace pwdft {
using namespace pwdft;

inline string mystring_readfile(const string fname)
{
   /* load the file into string */
   std::ifstream ifile(fname);
   std::string aa((std::istreambuf_iterator<char>(ifile)),
                   std::istreambuf_iterator<char>());
   return aa;
}

inline void mystring_writefile(const string fname, const string& s)
{
    std::ofstream ofile(fname);
    ofile << s;
}

inline int mystring_contains(const string s, const string a)
{
   return (s.find(a) != std::string::npos) ;
}

const string WHITESPACE = " \n\r\t\f\v";
const string WHITESPACE2 = "/ \n\r\t\f\v";

inline string mystring_ltrim(const string& s)
{
	size_t start = s.find_first_not_of(WHITESPACE);
	return (start == string::npos) ? "" : s.substr(start);
}

inline string mystring_rtrim(const string& s)
{
	size_t end = s.find_last_not_of(WHITESPACE);
	return (end == string::npos) ? "" : s.substr(0, end + 1);
}

inline string mystring_rtrim_slash(const string& s)
{
	size_t end = s.find_last_not_of(WHITESPACE2);
	return (end == string::npos) ? "" : s.substr(0, end + 1);
}


inline string mystring_trim(const string& str)
{
    size_t first = str.find_first_not_of(' ');
    if (string::npos == first)
    {
        return str;
    }
    size_t last = str.find_last_not_of(' ');
    return str.substr(first, (last - first + 1));
}

inline vector<string> mystring_split(string s, string delimiter) {
    size_t pos_start = 0, pos_end, delim_len = delimiter.length();
    string token;
    vector<string> res;
    while ((pos_end = s.find(delimiter, pos_start)) != string::npos) {
        token = s.substr(pos_start, pos_end - pos_start);
        pos_start = pos_end + delim_len;
        res.push_back(token);
    }
    res.push_back(s.substr(pos_start));
    return res;
}

inline vector<string> mystring_split0(string s) {
   std::vector<std::string> result;
   std::istringstream iss(s);
   for(std::string s; iss >> s; )
      result.push_back(s);

   return result;
}


// fetch a vector of doubles 
inline std::vector<double> mystring_double_list(std::string s, std::string delimiter) {
   std::vector<double> result;
   std::istringstream  iss(s.substr(s.find(delimiter) + delimiter.length()));

   double x; while (iss>>x) result.push_back(x);

   return result;
}



// Make a lowercase copy of s
inline string mystring_lowercase(const string& s) {
   string lower(s);
   for(size_t i = 0; i < s.length(); ++i)
   lower[i] = tolower(lower[i]);
   return lower;
}


// Make a capitzlise copy of s
inline string mystring_capitalize(const string& s) {
   string cap(s);
   for(size_t i = 0; i < cap.length(); ++i)
      cap[i] = tolower(cap[i]);
   cap[0] = toupper(cap[0]);
   return cap;
}

inline string mystring_ireplace(const string s0, const string a, const string b)
{
   string s(s0);
   int posold = -1;
   int posnew = mystring_lowercase(s).find(mystring_lowercase(a));
   while ((posold!=posnew) && (posnew>-1))
   {
      auto pos = std::search(s.begin(), s.end(), a.begin(), a.end(), [](const char c1, const char c2){ return (std::tolower(c1)) == (std::tolower(c2));});
      if(pos == s.end())
           return "";
      auto pos2 = pos;
      //std::cout << *pos << std::endl;
      std::advance(pos2, a.size());
      s.replace(pos, pos2, b);

      posold = posnew;
      posnew = mystring_lowercase(s).find(mystring_lowercase(a));
      //std::cout << s << "posold=" << posold << " posnew=" << posnew << std::endl;
   }
   return s;
}

inline bool mystring_isfloat(const string s)
{
   std::istringstream iss(s);
   float f;
   iss >> noskipws >> f;
   return iss.eof() && !iss.fail();
}

}


#endif
