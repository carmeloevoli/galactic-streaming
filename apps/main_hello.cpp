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

  {
    auto Wtest = utils::Grid<double>(2, 2, 3.5);
    std::cout << Wtest.get(0) << "\n";
    Wtest *= 2.0;
    std::cout << Wtest.get(0) << "\n";
  }

  {
    auto Wa = utils::Grid<double>(2, 2, 1);
    auto Wb = utils::Grid<double>(2, 2, 2);
    std::cout << Wa.get(0) << " " << Wb.get(0) << "\n";
    Wa += Wb;
    std::cout << Wa.get(0) << "\n";
    Wa -= Wb;
    std::cout << Wa.get(0) << "\n";
  }

  return 0;
}
