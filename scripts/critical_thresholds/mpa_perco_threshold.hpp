
// Eigen library
#include <Eigen/Core>
#include <Eigen/SparseCore>
// Spectra library
#include <Spectra/GenEigsSolver.h>
#include <Spectra/MatOp/SparseGenMatProd.h>
// Standard template library
#include <iostream>
#include <string>
// Graph c++ library
#include "undirected_graph_t.hpp"


double get_threshold_mpa(std::string edgelist_filename)
{
  // Creates the graph object from an edgelist.
  undirected_graph_t graph(edgelist_filename);
  int nb_vertices = graph.g_prop["nb_vertices"];

  // Builds the matrix related to the Ihara-Bass determinant formula.
  Eigen::SparseMatrix<double> M(2 * nb_vertices, 2 * nb_vertices);
  M.reserve(Eigen::VectorXi::Constant(2 * nb_vertices, 1));
  

  // Builds the first block corresponding to the adjacency matrix.
  std::set< std::pair<int, int> >::iterator it1 = graph.edgelist.begin();  
  std::set< std::pair<int, int> >::iterator end = graph.edgelist.end();  
  for(int i, j; it1!=end; ++it1)
  {
    // Identifies the vertices.
    i = it1->first;
    j = it1->second;
    M.insert(i, j) = 1;
    M.insert(j, i) = 1;
  }

  // Builds the second block corresponding to the identity matrix minus the degree matrix.
  graph.degrees();
  // std::vector<double>& Vertex2Degree = graph.v_prop["degree"];
  for(int i(0), j(nb_vertices); i<nb_vertices; ++i, ++j)
  {
    M.insert(i, j) = 1 - graph.v_prop["degrees"][i];
  }

  // Builds the third block corresponding to the identity matrix.
  for(int i(nb_vertices), j(0); j<nb_vertices; ++i, ++j)
  {
    M.insert(i, j) = 1;
  }

  // Construct matrix operation object using the wrapper class SparseGenMatProd
  Spectra::SparseGenMatProd<double> op(M);
  // Construct eigen solver object, requesting the largest eigenvalue
  Spectra::GenEigsSolver< double, Spectra::LARGEST_REAL, Spectra::SparseGenMatProd<double> > eigs(&op, 1, 6);
  // Initialize and compute
  eigs.init();
  int nconv = eigs.compute();

  // Retrieve results
  Eigen::VectorXcd evalues;
  if(eigs.info() == Spectra::SUCCESSFUL)
  {
    evalues = eigs.eigenvalues();
    return 1 / std::real(evalues[0]);
  }
  else
  {
    if(eigs.info() == Spectra::NOT_COMPUTED)
      std::cerr << "NOT_COMPUTED" << std::endl;
    if(eigs.info() == Spectra::NOT_CONVERGING)
      std::cerr << "NOT_CONVERGING" << std::endl;
    if(eigs.info() == Spectra::NUMERICAL_ISSUE)
      std::cerr << "NUMERICAL_ISSUE" << std::endl;
    return -1000;
  }
}




// double get_threshold_mpa(std::string edgelist_filename)
// {
//   // Creates the graph object from an edgelist.
//   undirected_graph_t graph(edgelist_filename);
//   int nb_edges = graph.g_prop["nb_edges"];

//   // Builds the non-backtracking matrix.
//   Eigen::SparseMatrix<double> B(2 * nb_edges, 2 * nb_edges);
//   B.reserve(Eigen::VectorXi::Constant(2 * nb_edges, 1));
//   std::set< std::pair<int, int> >::iterator it1 = graph.edgelist.begin();  
//   std::set< std::pair<int, int> >::iterator it2;
//   std::set< std::pair<int, int> >::iterator end = graph.edgelist.end();  
//   for(int i, j, k, l, cnt_i(0), cnt_j; it1!=end; ++it1, cnt_i+=2)
//   {
//     // Identifies the vertices.
//     i = it1->first;
//     j = it1->second;
//     cnt_j = 0;
//     it2 = graph.edgelist.begin();
//     for(; it2!=end; ++it2, cnt_j+=2)
//     {
//       // Identifies the vertices.
//       k = it2->first;
//       l = it2->second;

//       // i -> j && k -> l
//       if((j == k) && (i != l))
//       {
//         B.insert(cnt_i, cnt_j) = 1;
//       }
//       // i -> j && l -> k
//       if((j == l) && (i != k))
//       {
//         B.insert(cnt_i, cnt_j+1) = 1;
//       }
//       // j -> i && k -> l
//       if((i == k) && (j != l))
//       {
//         B.insert(cnt_i+1, cnt_j) = 1;
//       }
//       // j -> i && l -> k
//       if((i == l) && (j != k))
//       {
//         B.insert(cnt_i+1, cnt_j+1) = 1;
//       }
//     }
//   }

//   // Construct matrix operation object using the wrapper class SparseGenMatProd
//   Spectra::SparseGenMatProd<double> op(B);
//   // Construct eigen solver object, requesting the largest eigenvalue
//   Spectra::GenEigsSolver< double, Spectra::LARGEST_REAL, Spectra::SparseGenMatProd<double> > eigs(&op, 1, 6);
//   // Initialize and compute
//   eigs.init();
//   int nconv = eigs.compute();

//   // Retrieve results
//   Eigen::VectorXcd evalues;
//   if(eigs.info() == Spectra::SUCCESSFUL)
//   {
//     evalues = eigs.eigenvalues();
//     return 1 / std::real(evalues[0]);
//   }
//   else
//   {
//     if(eigs.info() == Spectra::NOT_COMPUTED)
//       std::cerr << "NOT_COMPUTED" << std::endl;
//     if(eigs.info() == Spectra::NOT_CONVERGING)
//       std::cerr << "NOT_CONVERGING" << std::endl;
//     if(eigs.info() == Spectra::NUMERICAL_ISSUE)
//       std::cerr << "NUMERICAL_ISSUE" << std::endl;
//     return -1000;
//   }
// }
