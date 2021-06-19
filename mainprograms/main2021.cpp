
///\cond
#include <iostream>
#include <math.h>
///\endcond
#include "Topology1d.h"
#include "TopologyTriangle.h"
#include "TopologyQuad.h"
#include "TopologyTetrahedron.h"
#include "Analysis.h"
#include "VTKGeoMesh.h"
#include "ReadGmsh.h"
#include "CompMesh.h"
#include "CompElement.h"
#include "Poisson.h"
#include "L2Projection.h"

using std::cout;
using std::endl;
using std::cin;

void force(const VecDouble &co, VecDouble &result)
{
    result.resize(1);
    result[0] = 1.;
}
int main (){
    
    GeoMesh gmesh;
    ReadGmsh read;
    read.Read(gmesh,"quads.msh");
    VTKGeoMesh plotmesh;
    plotmesh.PrintGMeshVTK(&gmesh, "quads.vtk");
    
    CompMesh cmesh(&gmesh);
    MatrixDouble perm(3,3);
    perm.setZero();
    perm.setIdentity();
    Poisson *mat1 = new Poisson(1, perm);
    mat1->SetForceFunction(force);
    MatrixDouble proj(1,1), val1(1,1), val2(1,1); 
    proj.setZero();
    val1.setZero();
    val2.setZero();
    L2Projection *bc_line = new L2Projection(0,2,proj,val1,val2);
    L2Projection *bc_point = new L2Projection(0,2,proj,val1,val2);
    //L2Projection *bc_line_right = new L2Projection(0,3,proj,val1,val2);
    //L2Projection *bc_line_top = new L2Projection(0,4,proj,val1,val2);
    //L2Projection *bc_line_left = new L2Projection(0,5,proj,val1,val2);
    //L2Projection *bc_point = new L2Projection(0,6,proj,val1,val2);
    std::vector<MathStatement *> mathvec = {0,mat1,bc_point,bc_line,mat1};
    //cmesh.SetMathStatement(1, mat1);
    cmesh.SetMathVec(mathvec);
    cmesh.AutoBuild();
    cmesh.Resequence();

    //cmesh.Solution() (0,0) = 1.;
    //CompElement* cel = cmesh.GetElement(0);
    //plotmesh.PrintCMeshVTK(&cmesh,2, "c_mesh.vtk");
   
    //MatrixDouble ek(4,4), ef (4,1);

    //cel->CalcStiff(ek,ef);

//  Analysis Analysis(&cmesh);
//  Analysis.RunSimulation();

    return 0;
}
