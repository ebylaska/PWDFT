#ifndef _PARSESTRING_H_
#define _PARSESTRING_H_
/* parsestring.hpp
   Author - Eric Bylaska
*/


#include <iostream>
#include <sstream>
#include <string>

using namespace std;

inline string trim(const string& str)
{
    size_t first = str.find_first_not_of(' ');
    if (string::npos == first)
    {
        return str;
    }
    size_t last = str.find_last_not_of(' ');
    return str.substr(first, (last - first + 1));
}

inline vector<string> split(string s, string delimiter) {
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

inline string lowerCase(const string& s) {
   string lower(s);
   for(size_t i = 0; i < s.length(); ++i)
   lower[i] = tolower(lower[i]);
   return lower;
}

#endif
