#include <iostream>
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


int main()
{
   json j;

   double x = 2.193939303;
   std::string s,line,nwinput;


   nwinput = "";
   while (getline(std::cin,line))
      nwinput += line + "\n";

   std::cout << "nwinput = \n" << nwinput<< std::endl << std::endl;

   //std::vector<std::string> lines = split(nwinput, "start");

   string dbname = trim(split(split(nwinput,"start")[1],"\n")[0]);

   std::cout << "dbname = " << dbname << std::endl;

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
