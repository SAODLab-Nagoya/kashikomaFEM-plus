#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>

//element class
struct Element
{
  void CalculateStiffnessMatrix(const Eigen::Matrix3f &D, std::vector<Eigen::Triplet<float>> &triplets, Eigen::VectorXf nodesX, Eigen::VectorXf nodesY);

  //nodes in this element
  int ne = 3;
  //B matrix
  Eigen::Matrix<float, 3, 6> B;
  //node ID in an element
  int nodesIds[3];
};

//constraint class
struct Constraint
{
  enum Type
  {
    UX = 1 << 0,
    UY = 1 << 1,
    UXY = UX | UY
  };
  int node;
  Type type;
};