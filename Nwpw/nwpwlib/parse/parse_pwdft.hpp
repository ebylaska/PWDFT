#ifndef _PARSE_PWDFT_HPP_
#define _PARSE_PWDFT_HPP_

#pragma once

#include <iostream>
#include <string>
#include <vector>

namespace pwdft {


extern std::string parse_nwinput(std::string);
extern std::string parse_rtdbstring(std::string);
extern int  parse_task(std::string);

extern bool parse_initialize_wvfnc(std::string,bool);
extern std::string parse_initialize_wvfnc_set(std::string,bool);

extern std::string parse_input_wavefunction_filename(std::string);
extern std::vector<std::string> parse_gen_lowlevel_rtdbstrs(std::string);
extern void parse_write(std::string);

}

#endif
