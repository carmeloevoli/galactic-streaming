#include "utils.hpp"

namespace utils {

double f_source_integral(double x, void* params) {
  double alpha = *(double*)params;
  double f = pow(x, 2. - alpha) * (sqrt(x * x + 1.) - 1.);
  return f;
}

double source_integral(double alpha, double x_max) {
  int LIMIT = 10000;
  gsl_integration_workspace* w = gsl_integration_workspace_alloc(LIMIT);
  double result, error;
  gsl_function F;
  F.function = &f_source_integral;
  F.params = &alpha;
  gsl_integration_qag(&F, 1, x_max, 0, 1e-6, LIMIT, 3, w, &result, &error);
  gsl_integration_workspace_free(w);
  return result;
}

}
