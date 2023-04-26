#ifndef _PARSESTRING_HPP_
#define _PARSESTRING_HPP_

#pragma once

#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <streambuf>
#include <string>
#include <vector>

namespace pwdft {

inline std::string mystring_readfile(const std::string fname) {
  /* load the file into string */
  std::ifstream ifile(fname);
  std::string aa((std::istreambuf_iterator<char>(ifile)),
                 std::istreambuf_iterator<char>());
  return aa;
}

inline void mystring_writefile(const std::string fname, const std::string &s) {
  std::ofstream ofile(fname);
  ofile << s;
}

inline int mystring_contains(const std::string s, const std::string a) {
  return (s.find(a) != std::string::npos);
}

const std::string WHITESPACE = " \n\r\t\f\v";
const std::string WHITESPACE2 = "/ \n\r\t\f\v";

inline std::string mystring_ltrim(const std::string &s) {
  size_t start = s.find_first_not_of(WHITESPACE);
  return (start == std::string::npos) ? "" : s.substr(start);
}

inline std::string mystring_rtrim(const std::string &s) {
  size_t end = s.find_last_not_of(WHITESPACE);
  return (end == std::string::npos) ? "" : s.substr(0, end + 1);
}

inline std::string mystring_rtrim_slash(const std::string &s) {
  size_t end = s.find_last_not_of(WHITESPACE2);
  return (end == std::string::npos) ? "" : s.substr(0, end + 1);
}

inline std::string mystring_trim(const std::string &str) {
  size_t first = str.find_first_not_of(' ');
  if (std::string::npos == first) {
    return str;
  }
  size_t last = str.find_last_not_of(' ');
  return str.substr(first, (last - first + 1));
}

inline std::vector<std::string> mystring_split(std::string s,
                                               std::string delimiter) {
  size_t pos_start = 0, pos_end, delim_len = delimiter.length();
  std::string token;
  std::vector<std::string> res;
  while ((pos_end = s.find(delimiter, pos_start)) != std::string::npos) {
    token = s.substr(pos_start, pos_end - pos_start);
    pos_start = pos_end + delim_len;
    res.push_back(token);
  }
  res.push_back(s.substr(pos_start));
  return res;
}

inline std::vector<std::string> mystring_split0(std::string s) {
  std::vector<std::string> result;
  std::istringstream iss(s);
  for (std::string s; iss >> s;)
    result.push_back(s);

  return result;
}

// fetch a vector of doubles
inline std::vector<double> mystring_double_list(std::string s,
                                                std::string delimiter) {
  std::vector<double> result;
  std::istringstream iss(s.substr(s.find(delimiter) + delimiter.length()));

  double x;
  while (iss >> x)
    result.push_back(x);

  return result;
}

// Make a lowercase copy of s
inline std::string mystring_lowercase(const std::string &s) {
  std::string lower(s);
  for (size_t i = 0; i < s.length(); ++i)
    lower[i] = tolower(lower[i]);
  return lower;
}

// Make a capitzlise copy of s
inline std::string mystring_capitalize(const std::string &s) {
  std::string cap(s);
  for (size_t i = 0; i < cap.length(); ++i)
    cap[i] = tolower(cap[i]);
  cap[0] = toupper(cap[0]);
  return cap;
}

inline std::string mystring_ireplace(const std::string s0, const std::string a,
                                     const std::string b) {
  std::string s(s0);
  int posold = -1;
  int posnew = mystring_lowercase(s).find(mystring_lowercase(a));
  while ((posold != posnew) && (posnew > -1)) {
    auto pos = std::search(s.begin(), s.end(), a.begin(), a.end(),
                           [](const char c1, const char c2) {
                             return (std::tolower(c1)) == (std::tolower(c2));
                           });
    if (pos == s.end())
      return "";
    auto pos2 = pos;
    // std::cout << *pos << std::endl;
    std::advance(pos2, a.size());
    s.replace(pos, pos2, b);

    posold = posnew;
    posnew = mystring_lowercase(s).find(mystring_lowercase(a));
    // std::cout << s << "posold=" << posold << " posnew=" << posnew <<
    // std::endl;
  }
  return s;
}

inline bool mystring_isfloat(const std::string s) {
  std::istringstream iss(s);
  float f;
  iss >> std::noskipws >> f;
  return iss.eof() && !iss.fail();
}

} // namespace pwdft

#endif
