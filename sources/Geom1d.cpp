/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
#include "Geom1d.h"

Geom1d::Geom1d() {
}

Geom1d::~Geom1d() {
}

Geom1d::Geom1d(const Geom1d &copy) {
    fNodeIndices = copy.fNodeIndices;
}

Geom1d& Geom1d::operator=(const Geom1d& copy) {
    fNodeIndices = copy.fNodeIndices;
    return *this;
}

void Geom1d::Shape(const VecDouble &xi, VecDouble &phi, MatrixDouble &dphi) {
//  std::cout << __PRETTY_FUNCTION__<< std::endl;
    if (xi.size() <=0 || xi.size() > Dimension) DebugStop();
//  if(xi.size() != Dimension || phi.size() != nCorners || dphi.rows() != Dimension || dphi.cols() != nCorners) DebugStop();
    phi.resize(nCorners);
    dphi.resize(Dimension, nCorners);

    double qsi; 
    qsi = xi[0];

    phi[0] = (1. - qsi) * 0.5;
    phi[1] = (1. + qsi) * 0.5; 

    dphi(0, 0) = -0.5;
    dphi(0, 1) = 0.5;
//  std::cout << "dphi = " << dphi << std::endl;
}

void Geom1d::X(const VecDouble &xi, MatrixDouble &NodeCo, VecDouble &x) {  
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
//  some prints to verify the numeric value of nrow, ncol, nnodes
//  std::cout << "nrow = " << nrow << std::endl;
//  std::cout << "ncol = " << nrow << std::endl;
//  std::cout << "nnodes = " << nnodes << std::endl;  
//  std::cout <<__LINE__ << std::endl;  

    for(int i = 0; i < nrow; i++) {
       for(int j = 0; j < nCorners; j++) {  
            x[i] += NodeCo(i,j) * phi[j];
       }
    }

//    std::cout <<__LINE__ << std::endl;
//    std::cout << "x = " << x << std::endl;
//    std::cout << "NodeCo = " << NodeCo << std::endl;

// Alternatively, x can be also implemented as:    
//    for (int i = 0; i < nrow; i++) {
//        x[i] = NodeCo(i, 0)*(1. - xi[0])*0.5 + NodeCo(i, 0)*(1. + xi[0])*0.5;
//    }
}

void Geom1d::GradX(const VecDouble &xi, MatrixDouble &NodeCo, VecDouble &x, MatrixDouble &gradx) {
//  std::cout <<__PRETTY_FUNCTION__<< std::endl;
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

//  std::cout <<__LINE__<< std::endl;

    for(int i = 0; i < nCorners; i++) {
       for(int j = 0; j < nrow; j++) {
 // x[j] += NodeCo(j,i) * phi[i]; I call the function above
            gradx(j,0) += NodeCo(j,i) * dphi(0,i);         
        } 
    }
//  some prints to verify the numeric value of Grad
//  std::cout << "Grad X = " << gradx << std::endl;
}

void Geom1d::SetNodes(const VecInt &nodes) {
    if(nodes.rows() != 2) {DebugStop();}
    fNodeIndices = nodes;
}

void Geom1d::GetNodes(VecInt &nodes) const{
    nodes = fNodeIndices;
}

int Geom1d::NodeIndex(int node) const{
    return fNodeIndices[node];
}

int Geom1d::NumNodes() {
    return nCorners;    
}

GeoElementSide Geom1d::Neighbour(int side) const {
    return fNeighbours[side];
}

void Geom1d::SetNeighbour(int side, const GeoElementSide &neighbour) {
    fNeighbours[side] = neighbour;
}
