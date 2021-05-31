/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "GeomTriangle.h"

GeomTriangle::GeomTriangle() {
}

GeomTriangle::~GeomTriangle() {
}

GeomTriangle::GeomTriangle(const GeomTriangle &copy) {
    fNodeIndices = copy.fNodeIndices;

}

GeomTriangle& GeomTriangle::operator=(const GeomTriangle& copy) {
    fNodeIndices = copy.fNodeIndices;

    return *this;
}

void GeomTriangle::Shape(const VecDouble& xi, VecDouble& phi, MatrixDouble& dphi) {
//    std::cout << __PRETTY_FUNCTION__<< std::endl;
//    if(xi.size() != Dimension || phi.size() != nCorners || dphi.rows() != Dimension || dphi.cols() != nCorners) DebugStop();
//    if (xi.size() <=0 || xi.size() > Dimension){DebugStop();}
    
    if (xi.size() <=0 || xi.size() > Dimension) DebugStop();
    phi.resize(nCorners);
    dphi.resize(Dimension, nCorners);

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

//    std::cout << "phi = " << phi << std::endl;
//    std::cout << "dphi = " << dphi << std::endl;
}

void GeomTriangle::X(const VecDouble &xi, MatrixDouble &NodeCo, VecDouble &x) {
//  std::cout << __PRETTY_FUNCTION__<< std::endl;
    int nrow = NodeCo.rows();
    int ncol = NodeCo.cols();
    int nnodes = NumNodes();

    if(xi.size() != Dimension) DebugStop();
    if(x.size() != NodeCo.rows()) DebugStop();
    if(NodeCo.cols() != nCorners) DebugStop();
    
    VecDouble phi;
    MatrixDouble dphi;
    x.resize(nrow);
    x.setZero();
    Shape(xi, phi, dphi);

// This part of the code follows the same logic from GeomTetrahedron lines
    for (int i = 0; i < nrow; i++) {
        for (int j = 0; j < nCorners; j++) {
            x[i] += phi[j] * NodeCo(i, j);
        }
    }
//  std::cout << "x = " << x << std::endl;
}

void GeomTriangle::GradX(const VecDouble &xi, MatrixDouble &NodeCo, VecDouble &x, MatrixDouble &gradx) {
//  std::cout << __PRETTY_FUNCTION__<< std::endl;
    if(xi.size() != Dimension) DebugStop();
    if(x.size() != NodeCo.rows()) DebugStop();
    if(NodeCo.cols() != nCorners) DebugStop();
    
    int nrow = NodeCo.rows();
    int ncol = NodeCo.cols();

    gradx.resize(nrow, Dimension);
    gradx.setZero();
    x.resize(nrow);
    x.setZero();
    
    VecDouble phi(nCorners);
    MatrixDouble dphi(Dimension, nCorners);
    Shape(xi, phi, dphi);
    X(xi,NodeCo,x);

    for (int i = 0; i < nCorners; i++) {
        for (int j = 0; j < nrow; j++) {
            gradx(j, 0) += NodeCo(j, i) * dphi(0, i);
            gradx(j, 1) += NodeCo(j, i) * dphi(1, i);     
        }
    }
//    std::cout << "GradX = " << gradx << std::endl;
//    std::cout <<__LINE__<< std::endl;
}

void GeomTriangle::SetNodes(const VecInt &nodes) {
    if(nodes.size() != nCorners) {DebugStop();}
    fNodeIndices = nodes;
}

void GeomTriangle::GetNodes(VecInt &nodes) const  {
    nodes = fNodeIndices;
}

int GeomTriangle::NodeIndex(int node) const {
    return fNodeIndices[node];
}

int GeomTriangle::NumNodes() {
    return nCorners;
}

GeoElementSide GeomTriangle::Neighbour(int side)  const {
    return fNeighbours[side];
}

void GeomTriangle::SetNeighbour(int side, const GeoElementSide &neighbour) {
    fNeighbours[side] = neighbour;
}