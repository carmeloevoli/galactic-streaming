#include <iostream>

#include "grid.hpp"
#include "hello_code.hpp"

int main() {
  blabber();

  auto W = utils::Grid<double>(100, 100);
  std::cout << "grid size is " << W.get_size() << "\n";
  std::cout << "0th element " << W.get(0) << " or " << W.get_value(0, 0) << "\n";

  auto otherW = utils::Grid<double>(W);
  std::cout << " Other grid " << W.get_size() << "\n";
  return 0;
}
