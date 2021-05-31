//
//  ShapeQuad.cpp
//  FemSC
//
//  Created by Philippe Devloo on 03/04/18.
//

#include "DataTypes.h"
#include "tpanic.h"
#include "Shape1d.h"
#include "ShapeQuad.h"

/// computes the shape functions in function of the coordinate in parameter space and orders of the shape functions (size of orders is number of sides of the element topology)
void ShapeQuad::Shape(const VecDouble &xi, VecInt &orders, VecDouble &phi, MatrixDouble &dphi){
    for (int i = 0; i < orders.size(); i++)
    {
        if (orders[i] < 0) {
            std::cout << "ShapeQuad::Shape: Invalid dimension for arguments: order\n";
            DebugStop();
        }
    }
    if (orders[0] > 1 || orders[1] > 1 || orders[2] > 1 || orders[3] > 1) {
        std::cout << "ShapeQuad::Shape: Invalid dimension for arguments: order\n";
        DebugStop();
    }

    int nshape = NShapeFunctions(orders);

    if (orders[nshape - 1] > 2) {
        std::cout << "ShapeQuad::Shape, only implemented until order = 2" << std::endl;
        DebugStop();
    }

    if(xi.size() != Dimension) DebugStop();
    if(orders.size() != nSides) DebugStop();

    int nsides = nSides;
    orders.resize(nSides);

    double qsi;
    double eta;

    qsi = xi[0]; 
    eta = xi[1];

    phi[0] = 1./4. * (1. - qsi) * (1. - eta);
    phi[1] = 1./4. * (1. + qsi) * (1. - eta);
    phi[2] = 1./4. * (1. + qsi) * (1. + eta);
    phi[3] = 1./4. * (1. - qsi) * (1. + eta);

    dphi(0,0) = 1./4. * (-1. * (1. - eta));
    dphi(1,0) = 1./4. * (-1. * (1. - qsi));
    
    dphi(0,1) = 1./4. * (+1. * (1. - eta));
    dphi(1,1) = 1./4. * (-1. * (1. + qsi));

    dphi(0,2) = 1./4. * (+1. * (1. + eta));
    dphi(1,2) = 1./4. * (+1. * (1. + qsi));

    dphi(0,3) = 1./4. * (-1. * (1. + eta));
    dphi(1,3) = 1./4. * (+1. * (1. - qsi));

    int count  = 4;
    int is;
    for (is = 4; is < 8; is++){
        if(orders[is] == 2)
        { 
             int id = is % 4;
             int id2 = (is + 1) % 4;
             int id3 = (is + 2) % 4;
             phi[count] = 4. * phi[id] * (phi[id2] + phi[id3]);
             dphi(0,count) = 4. * (dphi(0, id) * (phi[id2] + phi[id3]) + phi[id] * (dphi(0, id2) + dphi(0, id3)));
             dphi(1,count) = 4. * (dphi(1, id) * (phi[id2] + phi[id3]) + phi[id] * (dphi(1, id2) + dphi(1, id3)));
             count++;

        }
        else if(orders[is] != 1) DebugStop();
    }
    if(orders[8] == 2)
    {
        phi[count] = 16.*phi[0]*phi[2];
        dphi(0, count) = 16.*(dphi(0, 0) * phi[2] + phi[0]*dphi(0, 2));
        dphi(1, count) = 16.*(dphi(1, 0) * phi[2] + phi[0]*dphi(1, 2));
        count++;
    }
    else if(orders[8] != 1 ) DebugStop();
    if(count != nshape) DebugStop();
//    std::cout << "Please implement me\n";
//    DebugStop();
}

/// returns the number of shape functions associated with a side
int ShapeQuad::NShapeFunctions(int side, int order){
    if(order < 1 || order >2) DebugStop();
    if(side < 4)
        return 1; //0 a 4
    else if(side < 8)
        return (order - 1); //6 a 14
    else if(side == 8)
        return ((order - 1)*(order-1));
    
    std::cout << "ShapeQuad::NShapeFunctions, bad parameter side " << side << std::endl;
    DebugStop();
    
    return 0;
}

/// returns the total number of shape functions
int ShapeQuad::NShapeFunctions(VecInt &orders){
    
    int res = 4;
    for(int in = 4; in < orders.size(); in++) {
        res += NShapeFunctions(in, orders[in]);
    }
    
    return res;
}
