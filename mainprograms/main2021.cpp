
///\cond
#include <iostream>
#include <math.h>
#include <fstream>

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
#include "GeoElement.h"
#include "IntRule.h"

using std::cout;
using std::endl;
using std::cin;

int main(){
    IntPointData data;
    data.axes.resize(2,3);
    data.axes.setZero();
    data.detjac = 1.;
    data.dphidksi.resize(2,3);
    data.dphidksi.setZero();
    data.dphidx.resize(2,3);
    data.dphidx.setZero();
    data.gradx.resize(3,2);
    data.gradx.setZero();
    data.ksi.resize(2,1);
    data.phi.resize(3,1);
    data.weight = 1.;
    data.x.resize(3,1);
    data.phi[0] = 0.3;
    data.phi[1] = 0.3;
    data.phi[2] = 0.4;
    data.dphidksi(0,0) = -1.;
    data.dphidksi(1,0) = -1.;
    data.dphidksi(0,1) = 1.;
    data.dphidksi(1,1) = 0.;
    data.dphidksi(0,2) = 0.;
    data.dphidksi(1,2) = 1.;

    data.dphidx(0,0) = -1./sqrt(2.);
    data.dphidx(1,0) = 1./sqrt(2.) - sqrt(2.);
    data.dphidx(0,1) = 1./sqrt(2.);
    data.dphidx(1,1) = -1./sqrt(2.);
    data.dphidx(0,2) =  0.;
    data.dphidx(1,2) = sqrt(2.);

    data.gradx(0,0) = 1.;
    data.gradx(1,0) = 1.;
    data.gradx(0,1) = 0.;
    data.gradx(1,1) = 1.;

    data.axes(0,0) = 1./sqrt(2.);
    data.axes(1,0) = 1./sqrt(2.);
    data.axes(0,1) =  -1./sqrt(2.);
    data.axes(1,1) = 1./sqrt(2.);

    data.ksi[0] = 0.3;
    data.ksi[1] = 0.4;

    data.weight = 0.2;
    data.x[0] = 0.3;
    data.x[1] =0.7;

    MatrixDouble perm(3,3);
    perm.setZero();
    perm(0,0) =1.;
    perm(1,1) =1.;
    perm(2,2) =1.;

    Poisson matpoisson(1,perm);
    MatrixDouble ek(3,3), ef(3,1);
    ek.setZero();
    ef.setZero();
    matpoisson.Contribute(data, data.weight, ek, ef);
    
    std::cout << '\n' << ek << std::endl; 
    std::cout << '\n' << ef << std::endl; 

    return 0;

}
