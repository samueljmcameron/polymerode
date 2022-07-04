#ifndef BEADRODPMER_INPUT_HPP
#define BEADRODPMER_INPUT_HPP

#include <string>
#include <vector>
#include <map>

namespace BeadRodPmer {
namespace input {
  std::vector<std::string> split_line(const std::string&);
  bool isInt(const std::string&,int&,std::string);
  bool isDouble(const std::string&,double&,std::string);
  void convertVariable(std::string &,
		       std::map<std::string, std::string> const&);
  
}
};
#endif
