#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include "element_constraint.h"

void SetConstraints(Eigen::SparseMatrix<float>::InnerIterator &it, int index);
void ApplyConstraints(Eigen::SparseMatrix<float> &K, const std::vector<Constraint> &constraints);

int nodesCount;
Eigen::VectorXf nodesX;
Eigen::VectorXf nodesY;
Eigen::VectorXf loads;
std::vector<Element> elements;
std::vector<Constraint> constraints;

int main(int argc, char *argv[])
{
  if (argc != 3)
  {
    std::cout << "usage: " << argv[0] << " <input file> <output file>\n";
    return 1;
  }

  std::ifstream infile(argv[1]);
  std::ofstream outfile(argv[2]);

  //input poisson's ratio & young's modulus
  std::string tmp;
  float poissonRatio, youngModulus;
  infile >> tmp >> poissonRatio >> youngModulus;

  //set D matrix
  Eigen::Matrix3f D;
  D << 1.0f, poissonRatio, 0.0f,
      poissonRatio, 1.0, 0.0f,
      0.0f, 0.0f, (1.0f - poissonRatio) / 2.0f;

  D *= youngModulus / (1.0f - pow(poissonRatio, 2.0f));

  //input node
  infile >> nodesCount;
  nodesX.resize(nodesCount);
  nodesY.resize(nodesCount);

  for (int i = 0; i < nodesCount; ++i)
  {
    infile >> nodesX[i] >> nodesY[i];
  }

  //input element
  int elementCount;
  infile >> elementCount;

  for (int i = 0; i < elementCount; ++i)
  {
    Element element;
    infile >> element.nodesIds[0] >> element.nodesIds[1] >> element.nodesIds[2];
    elements.push_back(element);
  }

  //input neumann boundary
  int constraintCount;
  infile >> constraintCount;

  for (int i = 0; i < constraintCount; ++i)
  {
    Constraint constraint;
    int type;
    infile >> constraint.node >> type;
    constraint.type = static_cast<Constraint::Type>(type);
    constraints.push_back(constraint);
  }

  //input dirichlet boundary
  loads.resize(2 * nodesCount);
  loads.setZero();

  int loadsCount;
  infile >> loadsCount;

  for (int i = 0; i < loadsCount; ++i)
  {
    int node;
    float x, y;
    infile >> node >> x >> y;
    loads[2 * node + 0] = x;
    loads[2 * node + 1] = y;
  }

  //get an element stiffness matrix
  std::vector<Eigen::Triplet<float>> triplets;
  for (std::vector<Element>::iterator it = elements.begin(); it != elements.end(); ++it)
  {
    it->CalculateStiffnessMatrix(D, triplets, nodesX, nodesY);
  }

  //set get stiffness matrix
  Eigen::SparseMatrix<float> globalK(2 * nodesCount, 2 * nodesCount);
  globalK.setFromTriplets(triplets.begin(), triplets.end());

  ApplyConstraints(globalK, constraints);

  //solver
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<float>> solver(globalK);

  Eigen::VectorXf displacements = solver.solve(loads);

  //output displacement and mises stress to vtk file
  outfile << displacements << std::endl;

  for (std::vector<Element>::iterator it = elements.begin(); it != elements.end(); ++it)
  {
    Eigen::Matrix<float, 6, 1> delta;
    delta << displacements.segment<2>(2 * it->nodesIds[0]),
        displacements.segment<2>(2 * it->nodesIds[1]),
        displacements.segment<2>(2 * it->nodesIds[2]);

    Eigen::Vector3f sigma = D * it->B * delta;
    float sigma_mises = sqrt(sigma[0] * sigma[0] - sigma[0] * sigma[1] + sigma[1] * sigma[1] + 3.0f * sigma[2] * sigma[2]);

    outfile << sigma_mises << std::endl;
  }
  return 0;
}

void SetConstraints(Eigen::SparseMatrix<float>::InnerIterator &it, int index)
{
  if (it.row() == index || it.col() == index)
  {
    it.valueRef() = it.row() == it.col() ? 1.0f : 0.0f;
  }
}

void ApplyConstraints(Eigen::SparseMatrix<float> &K, const std::vector<Constraint> &constraints)
{
  std::vector<int> indicesToConstraint;

  for (std::vector<Constraint>::const_iterator it = constraints.begin(); it != constraints.end(); ++it)
  {
    if (it->type & Constraint::UX)
    {
      indicesToConstraint.push_back(2 * it->node + 0);
    }
    if (it->type & Constraint::UY)
    {
      indicesToConstraint.push_back(2 * it->node + 1);
    }
  }

  for (int k = 0; k < K.outerSize(); ++k)
  {
    for (Eigen::SparseMatrix<float>::InnerIterator it(K, k); it; ++it)
    {
      for (std::vector<int>::iterator idit = indicesToConstraint.begin(); idit != indicesToConstraint.end(); ++idit)
      {
        SetConstraints(it, *idit);
      }
    }
  }
}