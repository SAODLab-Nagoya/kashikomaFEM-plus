#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>
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
  //UX=1, UY=2, UXY=3
  enum Type
  {
    UX = 1 << 0,
    UY = 1 << 1,
    UXY = UX | UY
  };
  //node number of constrainted
  int node;
  Type type;
};