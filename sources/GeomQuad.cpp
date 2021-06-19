/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
#include "GeomQuad.h"

GeomQuad::GeomQuad() {
}

GeomQuad::~GeomQuad() {
}

GeomQuad::GeomQuad(const GeomQuad &copy) {
    fNodeIndices = copy.fNodeIndices;
}

GeomQuad& GeomQuad::operator=(const GeomQuad& copy) {
    fNodeIndices = copy.fNodeIndices;
    return *this;
}

void GeomQuad::Shape(const VecDouble &xi, VecDouble &phi, MatrixDouble &dphi) {
//  std::cout << __PRETTY_FUNCTION__<< std::endl;
    if (xi.size() <=0 || xi.size() > Dimension) DebugStop();
//  if(xi.size() != Dimension || phi.size() != nCorners || dphi.rows() != Dimension || dphi.cols() != nCorners) {DebugStop();}
    phi.resize(nCorners);
    dphi.resize(Dimension, nCorners);
    
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

//    std::cout << "phi = " << phi << std::endl;
//    std::cout << "dphi = " << dphi << std::endl;
}

void GeomQuad::X(const VecDouble &xi, MatrixDouble &NodeCo, VecDouble &x) {
//  std::cout <<__PRETTY_FUNCTION__<< std::endl;
    int nrow = NodeCo.rows();
    int ncol = NodeCo.cols();
    
    if(xi.size() != Dimension) DebugStop();
    if(x.size() < NodeCo.rows()) DebugStop();
    if(NodeCo.cols() != nCorners) DebugStop();
    
    VecDouble phi(nCorners);
    MatrixDouble dphi(Dimension, nCorners);
    
    if (x.size() < nrow){
        x.resize(nrow);
    }
    x.setZero();
    
    Shape(xi, phi, dphi);
// This part of the code follows the same logic from GeomTetrahedron lines
    for (int i = 0; i < Dimension; i++) {
        for (int j = 0; j < nCorners; j++) {
            x[i] += phi[j] * NodeCo(i, j);
        }
    }
//  std::cout <<__LINE__ << std::endl;
//  std::cout << "x = " << x << std::endl;
//  std::cout << "NodeCo = " << NodeCo << std::endl;
//  std::cout << "nCorners = " << nCorners << std::endl;
}

void GeomQuad::GradX(const VecDouble &xi, MatrixDouble &NodeCo, VecDouble &x, MatrixDouble &gradx) {
//  std::cout <<__PRETTY_FUNCTION__<< std::endl;
    if(xi.size() != Dimension) DebugStop();
    if(x.size() < NodeCo.rows()) DebugStop();
    if(NodeCo.cols() != nCorners) DebugStop();

    int nrow = NodeCo.rows();
    int ncol = NodeCo.cols(); 

    if(gradx.rows() < nrow || gradx.cols() < Dimension){
        gradx.resize(nrow, Dimension);
    }   
    gradx.setZero();

    if (x.size() < nrow){
        x.resize(nrow);
    }
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
//  std::cout << "GradX = " << gradx << std::endl;
}

void GeomQuad::SetNodes(const VecInt &nodes) {
    if(nodes.size() != nCorners) {DebugStop();}
    fNodeIndices = nodes;
}

void GeomQuad::GetNodes(VecInt &nodes) const{
    nodes = fNodeIndices;
}

int GeomQuad::NodeIndex(int node) const {
    return fNodeIndices[node];
}

int GeomQuad::NumNodes() {
    return nCorners;
}

GeoElementSide GeomQuad::Neighbour(int side) const {
    return fNeighbours[side];
}

void GeomQuad::SetNeighbour(int side, const GeoElementSide &neighbour) {
    fNeighbours[side] = neighbour;
}
