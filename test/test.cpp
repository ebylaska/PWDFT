#include <iostream>
#include <sstream>
#include "json.hpp"
#include <string>


using json = nlohmann::json;

using namespace std;


string trim(const string& str)
{
    size_t first = str.find_first_not_of(' ');
    if (string::npos == first)
    {
        return str;
    }
    size_t last = str.find_last_not_of(' ');
    return str.substr(first, (last - first + 1));
}

/*
vector<string> split(const string& str, const string& delim)
{
    vector<string> tokens;
    size_t prev = 0, pos = 0;
    do
    {
        pos = str.find(delim, prev);
        if (pos == string::npos) pos = str.length();
        string token = str.substr(prev, pos-prev);
        if (!token.empty()) tokens.push_back(token);
        prev = pos + delim.length();
    }
    while (pos < str.length() && prev < str.length());
    return tokens;
}
*/

vector<string> split(string s, string delimiter) {
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

// Make a lowercase copy of s
inline string lowerCase(const string& s) {
   string lower(s);
   for(size_t i = 0; i < s.length(); ++i)
   lower[i] = tolower(lower[i]);
   return lower;
}

string ireplace(string s0,string a, string b)
{
   string s(s0);
   auto pos = std::search(s.begin(), s.end(), a.begin(), a.end(), [](const char c1, const char c2){ return (std::tolower(c1)) == (std::tolower(c2));});
   if(pos == s.end())
        return "";
   auto pos2 = pos;
   std::cout << *pos << std::endl;
   std::advance(pos2, a.size());
   s.replace(pos, pos2, b);

   std::cout << "s=" << s << std::endl;

   return s;
}

int main()
{
   json j;

   double x = 2.193939303;
   std::string s,line,nwinput;
   json rtdb;


   nwinput = "";
   while (getline(std::cin,line))
      nwinput += line + "\n";

   std::cout << "nwinput = \n" << nwinput<< std::endl << std::endl;

   //std::vector<std::string> lines = split(nwinput, "start");

   string dbname = trim(split(split(nwinput,"start")[1],"\n")[0]);
   vector<string> lines = split(nwinput,"\n");

   // Remove comments
   for (auto i = lines.begin(); i != lines.end(); ++i) 
      *i = split(*i,"#")[0];

   std::cout << "dbname = " << dbname << std::endl;
   int n = lines.size();
   int cur = 0;
   while (cur<n)
   {
      if (lowerCase(lines[cur]).find("geometry")    != std::string::npos)
      {
         std::cout << "geometry: " << lines[cur] << std::endl;
      }
      else if (lowerCase(lines[cur]).find("title")  != std::string::npos)
      {
         std::cout << "title: " << lines[cur] << std::endl;
         std::cout << "ireplace:" << ireplace(lines[cur],"MY","my") << std::endl;
      }
      else if (lowerCase(lines[cur]).find("start")  != std::string::npos)
      {
         std::cout << "start: " << lines[cur] << std::endl;
         rtdb["dbname"] = trim(split(split(lines[cur],"start")[1],"\n")[0]);
      }
      else if (lowerCase(lines[cur]).find("charge") != std::string::npos)
         std::cout << "charge: " << lines[cur] << std::endl;
      else if (lowerCase(lines[cur]).find("nwpw")   != std::string::npos)
         std::cout << "nwpw: " << lines[cur] << std::endl;
      else if (lowerCase(lines[cur]).find("driver") != std::string::npos)
         std::cout << "driver: " << lines[cur] << std::endl;
      else if (lowerCase(lines[cur]).find("task")   != std::string::npos)
         std::cout << "task: " << lines[cur] << std::endl;

      ++cur;
   }

   std::cout << "rtdb = " << rtdb.dump() << std::endl;

/*
   json k = json::parse(s);

   j["pi"] = 3.14159;
   j["eric"] = "programmer";
   j["x"] = x;

   std::cout << j.dump() << std::endl;
   std::cout << k.dump() << std::endl;
   std::cout << k["Eric"] << std::endl;
*/
}
