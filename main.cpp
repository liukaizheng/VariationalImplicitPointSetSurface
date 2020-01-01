#include <iostream>
#include <Eigen/Dense>
#include <vector>
#include <utils/ReadXYZ.h>
#include <internal/STLContainerTraits.h>
#include "VariationalImplicitPointSetSuface.h"

#define PrintMSG(msg) \
    std::cout << #msg << ":\t" << msg << "\n";

int main()
{
    Internal::STLVectorTraits<Eigen::Vector3d>::Vector points;
    const std::string fileName("input.xyz");
    ReadXYZ(fileName, points);
    PrintMSG(points.size());

    VariationalImplicitPointSetSurface<Eigen::Vector3d, CubicRadialKernel>  vipss(points, 0.1);
    vipss.Run();
}
