#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "element_constraint.h"

void ApplyConstraints(Eigen::SparseMatrix<float> &K, const std::vector<Node> &nodes, Eigen::VectorXf &loads);
void output(char *outputPass, Eigen::VectorXf &displacements, std::vector<float> &sigma_mises);

//variables about node
int nodesCount;
std::vector<Node> nodes;

//variables about element
int elementCount;
std::vector<Element> elements;

//variables about Boundary Condition
Eigen::VectorXf loads;

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
  Eigen::Matrix3f De;
  De << 1.0f, poissonRatio, 0.0f,
      poissonRatio, 1.0, 0.0f,
      0.0f, 0.0f, (1.0f - poissonRatio) / 2.0f;

  De *= youngModulus / (1.0f - pow(poissonRatio, 2.0f));

  //input node's coordinate
  infile >> nodesCount;
  nodes.resize(nodesCount);

  for (int i = 0; i < nodesCount; ++i)
    infile >> nodes[i].x[0] >> nodes[i].x[1];

  //input elementf
  infile >> elementCount;

  for (int i = 0; i < elementCount; ++i)
  {
    Element element;
    infile >> element.nodesIds[0] >> element.nodesIds[1] >> element.nodesIds[2];
    elements.push_back(element);
  }

  //input dirichlet boundary
  int constraintCount;
  infile >> constraintCount;

  for (int i = 0; i < constraintCount; ++i)
  {
    int nodeID, type;
    infile >> nodeID >> type;
    if (nodes[nodeID].dirich == nullptr)
    {
      std::shared_ptr<Dirichlet> dirich(new Dirichlet);
      if (type == 1)
        dirich->flag[0] = 1;
      else if (type == 2)
        dirich->flag[1] = 1;
      else if (type == 3)
      {
        dirich->flag[0] = 1;
        dirich->flag[1] = 1;
      }
      else
      {
        std::cerr << "please put correct dirichlet value in input file" << std::endl;
        exit(1);
      }
      nodes[nodeID].dirich = dirich;
    }
  }

  //input neumann boundary
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

  //get an element stiffness matrix's triplets
  std::vector<Eigen::Triplet<float>> triplets;
  for (std::vector<Element>::iterator it = elements.begin(); it != elements.end(); ++it)
    it->CalculateStiffnessMatrix(De, triplets, nodes);

  //set global stiffness matrix
  Eigen::SparseMatrix<float> globalK(2 * nodesCount, 2 * nodesCount);
  globalK.setFromTriplets(triplets.begin(), triplets.end());

  ApplyConstraints(globalK, nodes, loads);

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

    Eigen::Vector3f sigma = De * it->B * delta;
    sigma_mises.push_back(sqrt(sigma[0] * sigma[0] - sigma[0] * sigma[1] + sigma[1] * sigma[1] + 3.0f * sigma[2] * sigma[2]));
  }

  output(argv[2], displacements, sigma_mises);

  std::cout << "kashikoma finished" << std::endl;
  return 0;
}

void ApplyConstraints(Eigen::SparseMatrix<float> &K, const std::vector<Node> &nodes, Eigen::VectorXf &loads)
{
  ///Set constraints on Kmatrix.
  Eigen::SparseMatrix<float> I(2 * nodesCount, 2 * nodesCount);
  I.setIdentity();
  Eigen::SparseMatrix<float> N(2 * nodesCount, 2 * nodesCount);
  std::vector<Eigen::Triplet<float>> triplets;
  for (int i = 0; i < nodesCount; i++)
  {
    if (nodes[i].dirich == nullptr)
    {
      for (int j = 0; j < 2; j++)
      {
        Eigen::Triplet<float> tmp(2 * i + j, 2 * i + j, 1);
        triplets.push_back(tmp);
      }
    }
    else
    {
      for (int j = 0; j < 2; j++)
      {
        if (nodes[i].dirich->flag[j] == 0)
        {
          Eigen::Triplet<float> tmp(2 * i + j, 2 * i + j, 1);
          triplets.push_back(tmp);
        }
        ///Change the external force to make it consistent with the constraint conditions.
        else
          loads[2 * i + j] = 0.0;
      }
    }
  }
  N.setFromTriplets(triplets.begin(), triplets.end());
  K = N.transpose() * K * N + (I - N);
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
    outfile << std::setw(16) << std::scientific << nodes[i].x[0]
            << std::setw(16) << std::scientific << nodes[i].x[1]
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