#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "element_constraint.h"

void SetConstraints(Eigen::SparseMatrix<float>::InnerIterator &it, int index);
void ApplyConstraints(Eigen::SparseMatrix<float> &K, const std::vector<Constraint> &constraints);
void output(char *outputPass, Eigen::VectorXf &displacements, std::vector<float> &sigma_mises);

//variables about node
int nodesCount;
Eigen::VectorXf nodesX;
Eigen::VectorXf nodesY;

//variables about element
int elementCount;
std::vector<Element> elements;

//variables about Boundary Condition
Eigen::VectorXf loads;
std::vector<Constraint> constraints;

int main(int argc, char *argv[])
{
  if (argc != 3)
  {
    std::cout << "usage: " << argv[0] << " <input file> <output file>\n";
    return 1;
  }

  std::ifstream infile(argv[1]);

  //input poisson's ratio & young's modulus
  float poissonRatio, youngModulus;
  infile >> poissonRatio >> youngModulus;

  //set D matrix(plane stress)
  Eigen::Matrix3f D;
  D << 1.0f, poissonRatio, 0.0f,
      poissonRatio, 1.0, 0.0f,
      0.0f, 0.0f, (1.0f - poissonRatio) / 2.0f;

  D *= youngModulus / (1.0f - pow(poissonRatio, 2.0f));

  //input node's coordinate
  infile >> nodesCount;
  nodesX.resize(nodesCount);
  nodesY.resize(nodesCount);

  for (int i = 0; i < nodesCount; ++i)
  {
    infile >> nodesX[i] >> nodesY[i];
  }

  //input element
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

  //set global stiffness matrix
  Eigen::SparseMatrix<float> globalK(2 * nodesCount, 2 * nodesCount);
  globalK.setFromTriplets(triplets.begin(), triplets.end());

  ApplyConstraints(globalK, constraints);

  Eigen::VectorXf displacements;

  //solver(deirect method)
  if (1)
  {
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<float>> solver;
    solver.compute(globalK);
    displacements = solver.solve(loads);
  }

  //solver(iterative method)
  else
  {
    Eigen::ConjugateGradient<Eigen::SparseMatrix<float>> solver;
    solver.compute(globalK);
    displacements = solver.solve(loads);
    std::cout << "#iterations:     " << solver.iterations() << std::endl;
    std::cout << "estimated error: " << solver.error() << std::endl;
  }

  //make von-Mises stress
  std::vector<float> sigma_mises;
  for (std::vector<Element>::iterator it = elements.begin(); it != elements.end(); ++it)
  {
    Eigen::Matrix<float, 6, 1> delta;
    delta << displacements.segment<2>(2 * it->nodesIds[0]),
        displacements.segment<2>(2 * it->nodesIds[1]),
        displacements.segment<2>(2 * it->nodesIds[2]);

    Eigen::Vector3f sigma = D * it->B * delta;
    sigma_mises.push_back(sqrt(sigma[0] * sigma[0] - sigma[0] * sigma[1] + sigma[1] * sigma[1] + 3.0f * sigma[2] * sigma[2]));
  }

  output(argv[2], displacements, sigma_mises);

  std::cout << "kashikoma finished" << std::endl;
  return 0;
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

void SetConstraints(Eigen::SparseMatrix<float>::InnerIterator &it, int index)
{
  if (it.row() == index || it.col() == index)
  {
    it.valueRef() = it.row() == it.col() ? 1.0f : 0.0f;
  }
}

void output(char *outputPass, Eigen::VectorXf &displacements, std::vector<float> &sigma_mises)
{
  std::ofstream outfile(outputPass);
  //header
  outfile << "# vtk DataFile Version 2.0\n"
          << "output of FEM program\n"
          << "ASCII\n\n"
          << "DATASET UNSTRUCTURED_GRID\n";

  //coordinate
  outfile << "POINTS"
          << " " << nodesCount << " "
          << "float" << std::endl;

  for (int i = 0; i < nodesCount; i++)
  {
    outfile << std::setw(16) << std::scientific << nodesX[i]
            << std::setw(16) << std::scientific << nodesY[i]
            << std::setw(16) << std::scientific << 0.0
            << std::endl;
  }
  outfile << std::endl;

  //connectivity
  outfile << "CELLS"
          << " " << elementCount << " " << (elements[0].ne + 1) * elementCount << std ::endl;

  for (int i = 0; i < elementCount; i++)
  {
    outfile << std::setw(8) << std::right << elements[i].ne;
    for (int j = 0; j < elements[i].ne; j++)
    {
      outfile << std::setw(8) << std::right << elements[i].nodesIds[j];
    }
    outfile << std::endl;
  }
  outfile << std::endl;

  //cell shape(triangle is 5,square is 9)
  outfile << "CELL_TYPES"
          << " " << elementCount << std::endl;
  for (int i = 0; i < elementCount; i++)
  {
    outfile << std::setw(8) << std::right << 5 << std::endl;
  }
  outfile << std::endl;

  //displacement
  outfile << "POINT_DATA"
          << " " << nodesCount << std::endl;
  outfile << "VECTORS displacement float" << std::endl;
  for (int i = 0; i < nodesCount; i++)
  {
    outfile << std::setw(16) << std::scientific << displacements[2 * i]
            << std::setw(16) << std::scientific << displacements[2 * i + 1]
            << std::setw(16) << std::scientific << 0.0
            << std::endl;
  }
  outfile << std::endl;

  //mises stress
  outfile << "CELL_DATA"
          << " " << elementCount << std::endl
          << "SCALARS mises_stress float" << std::endl
          << "LOOKUP_TABLE default" << std::endl;

  for (int i = 0; i < elementCount; i++)
  {
    outfile << std::setw(16) << std::scientific << sigma_mises[i]
            << std::endl;
  }

  outfile << std::endl;
}