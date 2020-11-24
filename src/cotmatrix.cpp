#include "cotmatrix.h"

// Cotagents calculated using derivations here:
// "A Cotangent Laplacian for Images as Surfaces"
void cotmatrix(
  const Eigen::MatrixXd & l,
  const Eigen::MatrixXi & F,
  Eigen::SparseMatrix<double> & L)
{
  // Add your code here
  int num_vertices = F.maxCoeff() + 1;
  L.resize(num_vertices, num_vertices);
  typedef Eigen::Triplet<double> T;
  std::vector<T> tripletList;
  tripletList.reserve(num_vertices * num_vertices);

  for(int i = 0; i < F.rows(); i++) {
    // Get edge lengths
    double l0, l1, l2;
    l0 = l(i, 0); // edge [1,2]
    l1 = l(i, 1); // edge [2,0]
    l2 = l(i, 2); // edge [0,1]
    // Get area
    double r = 0.5 * (l0 + l1 + l2);
    double area = std::sqrt(r * (r-l0) * (r-l1) * (r-l2));
    // Get cotangents
    double cot0, cot1, cot2;
    cot0 = (-1.0 * l0*l0 + l1*l1 + l2*l2) / (4.0 * area);
    cot1 = (-1.0 * l1*l1 + l0*l0 + l2*l2) / (4.0 * area);
    cot2 = (-1.0 * l2*l2 + l0*l0 + l1*l1) / (4.0 * area);

    // alpha of edge [i,j] = beta of edge [j,i]
    // L12, L21
    tripletList.push_back(T(F(i, 1), F(i,2), 0.5 * cot0));
    tripletList.push_back(T(F(i, 2), F(i,1), 0.5 * cot0));

    // L20, L02
    tripletList.push_back(T(F(i, 2), F(i,0), 0.5 * cot1));
    tripletList.push_back(T(F(i, 0), F(i,2), 0.5 * cot1));

    // L01, L10
    tripletList.push_back(T(F(i, 0), F(i,1), 0.5 * cot2));
    tripletList.push_back(T(F(i, 1), F(i,0), 0.5 * cot2));

    // L00, L11, L22
    tripletList.push_back(T(F(i, 0), F(i,0), -0.5 * cot1 - 0.5 * cot2));
    tripletList.push_back(T(F(i, 1), F(i,1), -0.5 * cot0 - 0.5 * cot2));
    tripletList.push_back(T(F(i, 2), F(i,2), -0.5 * cot0 - 0.5 * cot1));
  }
  L.setFromTriplets(tripletList.begin(), tripletList.end());
}

