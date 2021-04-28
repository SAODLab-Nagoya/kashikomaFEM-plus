#pragma once

#include <memory>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>

//constraint class
struct Dirichlet
{
  int flag[2];   ///<judge whether constraint is effective or not in each DOF
  double val[2]; ///<constrained value
};

///node class
struct Node
{
  double x[2];                       ///<coordinate
  std::shared_ptr<Dirichlet> dirich; ///<If the constraint is valid, this appears. Otherwise, this is nullptr.
};

//element class
struct Element
{
  void CalculateStiffnessMatrix(const Eigen::Matrix3f &D, std::vector<Eigen::Triplet<float>> &triplets, std::vector<Node> &nodes);

  //nodes in this element
  int ne = 3;
  //B matrix
  Eigen::Matrix<float, 3, 6> B;
  //node ID in an element
  int nodesIds[3];
};