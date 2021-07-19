
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
#include "CompElement.h"
#include "GeoElement.h"
#include "Assemble.h"
#include "PostProcessTemplate.h"
#include "PostProcess.h"

using std::cout;
using std::endl;
using std::cin;


void force(const VecDouble &x, VecDouble &res)
{
    res.resize(1);
    res[0] = 1;
}

int main (){  
    GeoMesh gmesh;
    ReadGmsh read;
    read.Read(gmesh,"triangle.msh");
    VTKGeoMesh plotmesh;
    plotmesh.PrintGMeshVTK(&gmesh, "triangle.vtk");
    CompMesh cmesh(&gmesh);
    MatrixDouble perm(3,3);
    perm.setZero();
    perm.setIdentity();
    Poisson *mat1 = new Poisson(1, perm);
    MatrixDouble proj(1,1), val1(1,1), val2(1,1); 
    proj.setZero();
    val1.setZero();
    val2.setZero();
    L2Projection *bc_line = new L2Projection(0,2,proj,val1,val2);
    L2Projection *bc_point = new L2Projection(0,1,proj,val1,val2);
    //L2Projection *bc_line_right = new L2Projection(0,3,proj,val1,val2);
    //L2Projection *bc_line_top = new L2Projection(0,4,proj,val1,val2);
    //L2Projection *bc_line_left = new L2Projection(0,5,proj,val1,val2);
    //L2Projection *bc_point = new L2Projection(0,6,proj,val1,val2);
    std::vector<MathStatement *> mathvec = {0,bc_point,bc_line,mat1};
    //cmesh.SetMathStatement(1, mat1);
    //cmesh.SetMathVec(mathvec);
    cmesh.SetMathVec(mathvec);
    cmesh.AutoBuild();
    cmesh.Resequence();  

for(auto cel:cmesh.GetElementVec())
    {
        MatrixDouble ek, ef;
        auto gel = cel->GetGeoElement();
        auto nnodes = gel->NNodes();
        VecInt nodeindices;
        IOFormat CommaInitFmt(StreamPrecision, DontAlignCols, ", ", ", ", "", ""," ", "");
        IOFormat HeavyFmt(FullPrecision, 0, ", ", "\n", "{", "}","{","}");
        gel->GetNodes(nodeindices);
        std::cout << "element index" << cel->GetIndex() << std::endl;
        std::cout << "coord = {";
        for (auto in = 0; in < nnodes; in++)
        {
            GeoNode &node = gmesh.Node(nodeindices[in]);
            std::cout << "{" << node.Co().format(CommaInitFmt) << "}";
            if(in < nnodes - 1) std::cout << ",";
        }
        std::cout << "};\n";
        cel ->CalcStiff(ek,ef);
        std::cout <<"ek = " << ek.format(HeavyFmt) << ";\n";
        std::cout <<"ef = " << ef.format(HeavyFmt) << ";\n";
    }

    Assemble assemble(&cmesh);
    auto neq = assemble.NEquations();
    MatrixDouble globmat(neq,neq),rhs(neq,1);
    assemble.Compute(globmat,rhs);

    cmesh.Solution() (0,0) = 1.;
    CompElement* cel = cmesh.GetElement(0);
    plotmesh.PrintCMeshVTK(&cmesh,2, "c_triangle.vtk");
    
    //  Analysis Analysis(&cmesh);
    //  Analysis.RunSimulation();

    return 0;

}
