#include <iostream>
#include <sstream>
#include "json.hpp"
#include <string>
#include <iomanip>


using json = nlohmann::json;

using namespace std;

json periodic_table_mass = json::parse("{ \"H\"  : 1.008, \"He\" : 4.0026, \"Li\" : 7.016, \"Be\" : 9.01218, \"B\"  : 11.00931, \"C\"  : 12.0, \"N\"  : 14.00307, \"O\"  : 15.99491, \"F\"  : 18.9984, \"Ne\" : 19.99244, \"Na\" : 22.9898, \"Mg\" : 23.98504, \"Al\" : 26.98154, \"Si\" : 27.97693, \"P\"  : 30.97376, \"S\"  : 31.97207, \"Cl\" : 34.96885, \"Ar\" : 39.9624, \"K\"  : 38.96371, \"Ca\" : 39.96259, \"Sc\" : 44.95592, \"Ti\" : 45.948, \"V\"  : 50.9440, \"Cr\" : 51.9405, \"Mn\" : 54.9381, \"Fe\" : 55.9349, \"Co\" : 58.9332, \"Ni\" : 57.9353, \"Cu\" : 62.9298, \"Zn\" : 63.9291, \"Ga\" : 68.9257, \"Ge\" : 73.9219, \"As\" : 74.9216, \"Se\" : 78.9183, \"Br\" : 79.9165, \"Kr\" : 83.912, \"Rb\" : 84.9117, \"Sr\" : 87.9056, \"Y\"  : 88.9054, \"Zr\" : 89.9043, \"Nb\" : 92.9060, \"Mo\" : 97.9055, \"Tc\" : 97.9072, \"Ru\" : 101.9037, \"Rh\" : 102.9048, \"Pd\" : 105.9032, \"Ag\" : 106.90509, \"Cd\" : 113.9036, \"In\" : 114.9041, \"Sn\" : 117.9018, \"Sb\" : 120.9038, \"Te\" : 129.9067, \"I\"  : 126.9004, \"Xe\" : 131.9042, \"Cs\" : 132.9051, \"Ba\" : 137.9050, \"La\" : 138.9061, \"Ce\" : 139.9053, \"Pr\" : 140.9074, \"Nd\" : 143.9099, \"Pm\" : 144.9128, \"Sm\" : 151.9195, \"Eu\" : 152.920, \"Gd\" : 157.9241, \"Tb\" : 159.9250, \"Dy\" : 163.9288, \"Ho\" : 164.9303, \"Er\" : 165.930, \"Tm\" : 168.9344, \"Yb\" : 173.9390, \"Lu\" : 174.9409, \"Hf\" : 179.9468, \"Ta\" : 180.948, \"W\"  : 183.9510, \"Re\" : 186.9560, \"Os\" : 189.9586, \"Ir\" : 192.9633, \"Pt\" : 194.9648, \"Au\" : 196.9666, \"Hg\" : 201.9706, \"Tl\" : 204.9745, \"Pb\" : 207.9766, \"Bi\" : 208.9804, \"Po\" : 209.9829, \"At\" : 210.9875, \"Rn\" : 222.0175, \"Fr\" : 223.0198, \"Ra\" : 226.0254, \"Ac\" : 227.0278, \"Th\" : 232.0382, \"Pa\" : 231.0359, \"U\"  : 238.0508, \"Np\" : 237.0482, \"Pu\" : 244.0642, \"Am\" : 243.0614, \"Cm\" : 247.0704, \"Bk\" : 247.0703, \"Cf\" : 251.0796, \"Es\" : 252.0829, \"Fm\" : 257.0950, \"Md\" : 258.0986, \"No\" : 259.1009, \"Lr\" : 262.1100, \"Rf\" : 261.1087, \"Ha\" : 262.1138, \"Sg\" : 266.1219, \"Bh\" : 262.1229, \"Hs\" : 267.1318, \"Mt\" : 268.1388 }");



inline int mystring_contains(const string s, const string a)
{
   return (s.find(a) != std::string::npos) ;
}

string mystring_trim(const string& str)
{
    size_t first = str.find_first_not_of(' ');
    if (string::npos == first)
    {
        return str;
    }
    size_t last = str.find_last_not_of(' ');
    return str.substr(first, (last - first + 1));
}


vector<string> mystring_split(string s, string delimiter) {
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

vector<string> mystring_split0(string s) {
   std::vector<std::string> result;
   std::istringstream iss(s);
   for(std::string s; iss >> s; )
      result.push_back(s);

   return result;
}

// Make a lowercase copy of s
inline string mystring_lowercase(const string& s) {
   string lower(s);
   for(size_t i = 0; i < s.length(); ++i)
   lower[i] = tolower(lower[i]);
   return lower;
}

inline string mystring_capitalize(const string& s) {
   string cap(s);
   for(size_t i = 0; i < cap.length(); ++i)
      cap[i] = tolower(cap[i]);
   cap[0] = toupper(cap[0]);
   return cap;
}

string mystring_ireplace(const string s0, const string a, const string b)
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


class mystring : public std::string {
   std::string mystring0;

public:
   mystring() { mystring0 = "";}
   mystring(const std::string s) { mystring0 = s;}
   mystring(const mystring& s) { mystring0 = s.mystring0;}
   mystring operator = (const mystring s)    { return mystring(s); }
   mystring operator = (const std::string s) { return mystring(s); }
   mystring operator + (const std::string s ) { return mystring(mystring0 + s); }

};

json parse_geometry(int *curptr, vector<string> lines)
{
   json geomjson;

   int cur = *curptr;
   int center = 1;
   int autoz = 0;
   int autosym = 0;
   double angs_to_au = 1.0/0.52917715;
   double conv = angs_to_au;
   if (mystring_lowercase(lines[cur]).find(" au")   != std::string::npos) conv = 1.0;
   if (mystring_lowercase(lines[cur]).find(" a.u.") != std::string::npos) conv = 1.0;
   if (mystring_lowercase(lines[cur]).find(" bo")   != std::string::npos) conv = 1.0;
   if (mystring_lowercase(lines[cur]).find(" an")   != std::string::npos) conv = angs_to_au;
   if (mystring_lowercase(lines[cur]).find(" nm")   != std::string::npos) conv = 10.0*angs_to_au;
   if (mystring_lowercase(lines[cur]).find(" na")   != std::string::npos) conv = 10.0*angs_to_au;
   if (mystring_lowercase(lines[cur]).find(" pm")   != std::string::npos) conv = 0.01*angs_to_au;
   if (mystring_lowercase(lines[cur]).find(" pi")   != std::string::npos) conv = 0.01*angs_to_au;

   if (mystring_lowercase(lines[cur]).find("nocenter")      != std::string::npos)  
       center = 0;
   else if (mystring_lowercase(lines[cur]).find("center")   != std::string::npos)  
       center = 1;
   if (mystring_lowercase(lines[cur]).find("noautoz")       != std::string::npos)
      autoz  = 0;
   else if (mystring_lowercase(lines[cur]).find("autoz")    != std::string::npos)
      autoz  = 1;
   if (mystring_lowercase(lines[cur]).find("noautosym")     != std::string::npos)
      autosym  = 0;
   else if (mystring_lowercase(lines[cur]).find("autosym")  != std::string::npos)
      autosym  = 1;

   geomjson["conv"]    = conv;
   geomjson["center"]  = center;
   geomjson["autoz"]   = autoz;
   geomjson["autosym"] = autosym;

   int endcount = 1;
   vector<string> ss;
   vector<string> symbols;
   vector<double> coords;
   vector<double> velocities;
   vector<double> mass;
   ++cur;
   int nion = 0;
   double mm;
   string line;
   while (endcount>0)
   {
      line = lines[cur];
      
      ss = mystring_split0(lines[cur]);
      symbols.push_back(mystring_capitalize(ss[0]));
      coords.push_back(std::stod(ss[1])*conv);
      coords.push_back(std::stod(ss[2])*conv);
      coords.push_back(std::stod(ss[3])*conv);
      velocities.push_back(0.0);
      velocities.push_back(0.0);
      velocities.push_back(0.0);

      mm = periodic_table_mass[mystring_capitalize(ss[0])];
      if (mystring_contains(mystring_lowercase(line),"mass"))
         mm = std::stod(mystring_split0(mystring_split(line,"mass")[1])[0]);
      mass.push_back(mm);

      ++nion;
      ++cur;
      if (mystring_contains(mystring_lowercase(lines[cur]),"end"))
         --endcount;
   }

   if (center)
   {
      int ii;
      double xcm=0.0;
      double ycm=0.0;
      double zcm=0.0;
      for (ii=0; ii<nion; ++ii)
      {
         xcm += coords[3*ii];
         ycm += coords[3*ii+1];
         zcm += coords[3*ii+2];
      }
      xcm /= ((double) nion);
      ycm /= ((double) nion);
      zcm /= ((double) nion);
      for (int ii=0; ii<nion; ++ii)
      {
         coords[3*ii]   -= xcm;
         coords[3*ii+1] -= ycm;
         coords[3*ii+2] -= zcm;
      }
   }

   geomjson["symbols"]    = symbols;
   geomjson["coords"]     = coords;
   geomjson["velocities"] = velocities;
   geomjson["nion"]       = nion;
   geomjson["mass"]       = mass;

   *curptr = cur;

   return geomjson;
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

   string dbname = mystring_trim(mystring_split(mystring_split(nwinput,"start")[1],"\n")[0]);
   vector<string> lines = mystring_split(nwinput,"\n");

   // Remove comments
   for (auto i = lines.begin(); i != lines.end(); ++i) 
      *i = mystring_split(*i,"#")[0];

   std::cout << "dbname = " << dbname << std::endl;
   int n = lines.size();
   int cur = 0;
   while (cur<n)
   {
      if (mystring_lowercase(lines[cur]).find("geometry")    != std::string::npos)
      {
         std::cout << "geometry: " << lines[cur] << std::endl;
         rtdb["geometry"] = parse_geometry(&cur,lines);
         //cur = parse_geometry_endcur(cur,lines);
      }
      else if (mystring_lowercase(lines[cur]).find("title")  != std::string::npos)
      {
         rtdb["title"] = mystring_trim(mystring_ireplace(mystring_split(mystring_ireplace(lines[cur],"TITLE","title"),"title")[1],"\"",""));
      }
      else if (mystring_lowercase(lines[cur]).find("start")  != std::string::npos)
      {
         rtdb["dbname"] = mystring_trim(mystring_split(mystring_split(lines[cur],"start")[1],"\n")[0]);
      }
      else if (mystring_lowercase(lines[cur]).find("charge") != std::string::npos)
      {
         std::cout << "charge: " << lines[cur] << std::endl;
         rtdb["charge"] = std::stoi(mystring_trim(mystring_split(mystring_split(lines[cur],"charge")[1],"\n")[0]));
      }
      else if (mystring_lowercase(lines[cur]).find("nwpw")   != std::string::npos)
      {
         std::cout << "nwpw: " << lines[cur] << std::endl;
      }
      else if (mystring_lowercase(lines[cur]).find("driver") != std::string::npos)
      {
         std::cout << "driver: " << lines[cur] << std::endl;
      }
      else if (mystring_lowercase(lines[cur]).find("task")   != std::string::npos)
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
