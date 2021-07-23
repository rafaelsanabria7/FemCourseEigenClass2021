//
//  TestTwoDProblem.cpp MODIFICADO DO ORIGINAL
//  FemSC
//
//  Created by Eduardo Ferri on 08/17/15.
//
//
//TestOneDProblem cpp
/*
 Os testes foram preparados com um proposito educacional,
 recomenda-se que o aluno entenda a funcionalidade de cada
 teste e posteriormente use com seu cÛdigo caso a caso
 */
//      Obs: O xmax e xmin estao tomados como 4 e 0, respectivamente,
//      caso estes valores sejam alterados, editar o teste TestNodes.
//
//
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

//void exact(const VecDouble &point,VecDouble &val, MatrixDouble &deriv);

int main ()
{
    GeoMesh gmesh;
    ReadGmsh read;
    read.Read(gmesh,"TwoDmesh_triangle_05_u.msh");
    VTKGeoMesh plotmesh;
    plotmesh.PrintGMeshVTK(&gmesh, "TwoDmesh_triangle_05_u.vtk");
    CompMesh cmesh(&gmesh);
    MatrixDouble perm(3,3);
    perm.setZero();
    perm(0,0) = 1.;
    perm(1,1) = 1.;
    perm(2,2) = 1.;
    Poisson *mat1 = new Poisson(1,perm);
    mat1->SetDimension(2);

    auto force = [](const VecDouble &x, VecDouble &res)
    {
        //res[0] = 33*exp(-x[1])*cos(3.*x[1])*sin(5*x[0])-6*exp(-x[1])*sin(5.*x[0])*sin(3*x[1]);
        res[0] = (-10*(-1 + 2*x[0])*(-1 + x[1])*x[1]*cos(5*x[0])*cos(3*x[1]) + sin(5*x[0])*((-1 + x[1])*(-2*x[1] - x[0]*(4 + 33*x[1]) + x[0]*x[0]*(4 + 33*x[1]))*cos(3*x[1]) - 6*(-1 + x[0])*x[0]*(1 - 3*x[1] + x[1]*x[1])*sin(3*x[1])))/exp(x[1]);
    };
    auto exact = [](const VecDouble &x, VecDouble &Sol, MatrixDouble &DSol)
    {
        //Sol[0] = sin(5*x[0])*exp(-x[1])*cos(3*x[1]);
        Sol[0] = sin(5*x[0])*exp(-x[1])*cos(3*x[1])*(x[0]*(1-x[0])*x[1]*(1-x[1]));
        DSol(0,0) = exp(-x[1])*(-1 + x[1])*x[1]*cos(3*x[1])*(5*(-1 + x[0])*x[0]*cos(5*x[0]) + (-1 + 2*x[0])*sin(5*x[0]));
        DSol(1,0) = -exp(-x[1])*(((-1 + x[0])*x[0]*sin(5*x[0])*((1 - 3*x[1] + x[1]*x[1])*cos(3*x[1]) + 3*(-1 + x[1])*x[1]*sin(3*x[1]))));
    };
    mat1->SetForceFunction(force);
    MatrixDouble proj(1,1),val1(1,1),val2(1,1);
    proj.setZero();
    val1.setZero();
    val2.setZero();
    L2Projection *bc_linha = new L2Projection(0,2,proj,val1,val2);
    //L2Projection *bc_point = new L2Projection(0,3,proj,val1,val2);
    //bc_linha ->SetExactSolution(exact);
    //bc_point ->SetExactSolution(exact);
    mat1 ->SetExactSolution(exact);
    std::vector<MathStatement *> mathvec = {0,mat1,bc_linha}; //Verify the order of mathvec
    cmesh.SetMathVec(mathvec);
    cmesh.SetDefaultOrder(2);
    cmesh.AutoBuild();
    cmesh.Resequence(); 

    Analysis Analysis(&cmesh);
    Analysis.RunSimulation();
    PostProcessTemplate<Poisson> postprocess;

    postprocess.AppendVariable("Sol");
    postprocess.AppendVariable("DSol");
    postprocess.AppendVariable("Flux");
    postprocess.AppendVariable("Force");
    postprocess.AppendVariable("SolExact");
    postprocess.AppendVariable("DSolExact");
    
    postprocess.SetExact(exact);
    //mat1->SetExactSolution(exact);
    Analysis.PostProcessSolution("cnew_quads.vtk",postprocess);
    //plotmesh.PrintCMeshVTK(&cmesh,2,"c_TwoDmesh_triangle_05_u.vtk");

    VecDouble errvec;
    errvec = Analysis.PostProcessError(std::cout,postprocess);

    return 0;
}



