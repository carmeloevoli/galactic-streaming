#include "galaxy.hpp"
#include "params.hpp"

int main() {
  auto params = utils::Params("");
  auto galaxy = core::Galaxy(params);
  galaxy.build_space_axis();
  galaxy.build_momentum_axis();
  galaxy.build_rigidity_axis();
  galaxy.build_initial_condition();
}