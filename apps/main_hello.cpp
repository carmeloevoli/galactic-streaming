#include "grid.hpp"
#include "hello_code.hpp"
int main() {
  blabber();

  auto W = utils::Grid2D<double>(100, 100);
  std::cout << "grid size is " << W.getSize() << "\n";
  std::cout << "0th element " << W.get(0) << " or " << W.getValue(0, 0) << "\n";
  return 0;
}
