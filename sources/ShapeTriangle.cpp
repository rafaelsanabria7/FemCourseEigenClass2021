//
//  ShapeTriangle.h
//  FemSC
//
//  Created by Philippe Devloo on 03/04/18.
//

#include "ShapeTriangle.h"
#include "Shape1d.h"

/// computes the shape functions in function of the coordinate in parameter space and orders of the shape functions (size of orders is number of sides of the element topology)
void ShapeTriangle::Shape(const VecDouble &xi, VecInt &orders, VecDouble &phi, MatrixDouble &dphi){
//  std::cout <<__PRETTY_FUNCTION__<< std::endl;
    for (int i = 0; i < orders.size(); i++)
    {
        if (orders[i] < 0) {
            std::cout << "ShapeTriangle::Shape: Invalid dimension for arguments: order\n";
            DebugStop();
        }
    }
    if (orders[0] > 1 || orders[1] > 1 || orders[2] > 1) {
        std::cout << "ShapeTriangle::Shape: Invalid dimension for arguments: order\n";
        DebugStop();
    }
    if(xi.size() != Dimension) DebugStop();
    if(orders.size() != nSides) DebugStop();

    int nshape = NShapeFunctions(orders);
    int nsides = nSides;

    if (orders[nshape - 1] > 2) {
       std::cout << "ShapeTriangle::Shape, only implemented until order = 2" << std::endl;
       DebugStop();
   }

//    std::cout << "order = " << orders << std::endl;
//    std::cout << "nshape = " << nshape << std::endl;
//    std::cout << "nsides = " << nsides << std::endl;
   
    double qsi;
    double eta;
    qsi = xi[0]; 
    eta = xi[1];
    
    // Linear order
    phi[0] =  1.- qsi - eta;
    phi[1] =  qsi;
    phi[2] =  eta;
    
    dphi(0,0) = -1.;
    dphi(1,0) = -1.;
    dphi(0,1) =  1.;
    dphi(1,1) =  0.;
    dphi(0,2) =  0.;
    dphi(1,2) =  1.; 

    int count = 3;
    int is;

    for (is = 3; is < 6; is++) 
        {
            if(orders[is] == 2)
            {
            int is1 = SideNodeLocIndex(is, 0);
            int is2 = SideNodeLocIndex(is, 1);
            phi[is] = 4. *phi[is1] * phi[is2];
            dphi(0, is) = 4. * (dphi(0, is1) * phi[is2] + phi[is1] * dphi(0, is2));
            dphi(1, is) = 4. * (dphi(1, is1) * phi[is2] + phi[is1] * dphi(1, is2));
//            std::cout << "count = " << count << std::endl;
            count++;

//            std::cout << "is1 = " << is1 << std::endl;
//            std::cout << "is2 = " << is2 << std::endl;
            }
            else if (orders[is] != 1) {DebugStop();}
        }

        if(orders[6] == 3)
    {
 //       std::cout <<__LINE__<< std::endl;
        phi[6] =  27.*phi[0]*phi[1]*phi[2];
        dphi(0, 6) =  27 * (eta - 2*eta*qsi - eta*eta);
        dphi(1, 6) =  27 * (qsi - qsi*qsi - 2*qsi*eta);
    }   
 //   std::cout << "Please implement me\n";
 //   DebugStop();   
}

/// returns the number of shape functions associated with a side
int ShapeTriangle::NShapeFunctions(int side, int order){
    switch(side) {
        case 0:
        case 1:
        case 2:
            return 1;
        case 3:
        case 4:
        case 5:
            return order-1;
        case 6:
            return 0;
    }
    
    DebugStop();
    std::cout << "ShapeTriangle::NShapeFunctions, bad parameter side " << std::endl;
    return 0;
}

/// returns the total number of shape functions
int ShapeTriangle::NShapeFunctions(VecInt &orders){
    
    int res=3;
    for(int in=3; in < orders.size(); in++) {
        res += NShapeFunctions(in, orders[in]);
    }
    
    return res;
    
}