#include "element_constraint.h"

void Element::CalculateStiffnessMatrix(const Eigen::Matrix3f &D, std::vector<Eigen::Triplet<float>> &triplets, Eigen::VectorXf nodesX, Eigen::VectorXf nodesY)
{
  Eigen::Vector3f x, y;
  x << nodesX[nodesIds[0]], nodesX[nodesIds[1]], nodesX[nodesIds[2]];
  y << nodesY[nodesIds[0]], nodesY[nodesIds[1]], nodesY[nodesIds[2]];

  Eigen::Matrix3f C;
  C << Eigen::Vector3f(1.0f, 1.0f, 1.0f), x, y;

  Eigen::Matrix3f IC = C.inverse();

  for (int i = 0; i < 3; i++)
  {
    B(0, 2 * i + 0) = IC(1, i);
    B(0, 2 * i + 1) = 0.0f;
    B(1, 2 * i + 0) = 0.0f;
    B(1, 2 * i + 1) = IC(2, i);
    B(2, 2 * i + 0) = IC(2, i);
    B(2, 2 * i + 1) = IC(1, i);
  }
  Eigen::Matrix<float, 6, 6> K = B.transpose() * D * B * C.determinant() / 2.0f;

  for (int i = 0; i < 3; i++)
  {
    for (int j = 0; j < 3; j++)
    {
      Eigen::Triplet<float> trplt11(2 * nodesIds[i] + 0, 2 * nodesIds[j] + 0, K(2 * i + 0, 2 * j + 0));
      Eigen::Triplet<float> trplt12(2 * nodesIds[i] + 0, 2 * nodesIds[j] + 1, K(2 * i + 0, 2 * j + 1));
      Eigen::Triplet<float> trplt21(2 * nodesIds[i] + 1, 2 * nodesIds[j] + 0, K(2 * i + 1, 2 * j + 0));
      Eigen::Triplet<float> trplt22(2 * nodesIds[i] + 1, 2 * nodesIds[j] + 1, K(2 * i + 1, 2 * j + 1));

      triplets.push_back(trplt11);
      triplets.push_back(trplt12);
      triplets.push_back(trplt21);
      triplets.push_back(trplt22);
    }
  }
}