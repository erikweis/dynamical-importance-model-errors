

// g++ -O3 -I src_eigenvalues/ mpa_perco_threshold.cpp -o mpa_perco_threshold

#include "mpa_perco_threshold.hpp"


int main(int argc, char const *argv[])
{
  // Filename of the edgelist.
  std::string name = argv[1];
  // Computes the threshold.
  double thrshld = get_threshold_mpa(name);
  // Writes in screen.
  int width = 20;
  //std::cout << std::setw(width) << name    << " ";
  std::cout << thrshld;
  //std::cout << std::setw(width) << thrshld << " ";
  //std::cout << std::endl;
  // Exits successfully.
  return 0;
}
