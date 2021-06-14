

#include <iostream>
#include <math.h>
#include "IntRule.h"
#include "IntRule1d.h"
#include "IntRuleQuad.h"
#include "IntRuleTetrahedron.h"
#include "IntRuleTriangle.h"
#include "Topology1d.h"
#include "TopologyTriangle.h"
#include "TopologyQuad.h"
#include "TopologyTetrahedron.h"
#include "DataTypes.h"
#include "Analysis.h"
#include "VTKGeoMesh.h"
#include "ReadGmsh.h"
#include "CompMesh.h"
#include "Poisson.h"
#include "L2Projection.h"

using std::cout;
using std::endl;
using std::cin;

int main (){
    
    GeoMesh gmesh;
    ReadGmsh read;
    read.Read(gmesh,"1element.msh");
    VTKGeoMesh plotmesh;
    plotmesh.PrintGMeshVTK(&gmesh, "1element.vtk");
    
    CompMesh cmesh(&gmesh);
    MatrixDouble perm(2,2);
    perm.setZero();
    perm.setIdentity();
    Poisson *mat1 = new Poisson(3, perm);

    MatrixDouble proj(1,1), val1(1,1), val2(1,1); 
    proj.setZero();
    val1.setZero();
    val2.setZero();
    L2Projection *bc_line_bottom = new L2Projection(0,2,proj,val1,val2);
    L2Projection *bc_line_right = new L2Projection(0,3,proj,val1,val2);
    L2Projection *bc_line_top = new L2Projection(0,4,proj,val1,val2);
    L2Projection *bc_line_left = new L2Projection(0,5,proj,val1,val2);
    L2Projection *bc_point = new L2Projection(0,6,proj,val1,val2);
    std::vector<MathStatement *> mathvec = {0,mat1,bc_line_bottom,bc_line_right,bc_line_top,bc_line_left,bc_point};
    cmesh.SetMathVec(mathvec);

    cmesh.AutoBuild();
    cmesh.Resequence();
    cmesh.Solution() (0,0) = 1.;

    plotmesh.PrintCMeshVTK(&cmesh,2, "c_mesh_1element.vtk");

//  Analysis Analysis(&cmesh);
//  Analysis.RunSimulation();

    return 0;
}
