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
    if(xi.size() != Dimension || phi.size() != nCorners || dphi.rows() != Dimension || dphi.cols() != nCorners)
        DebugStop();
    double qsi = xi[0];

    phi[0] = (1. - qsi) * 0.5;
    phi[1] = (1. + qsi) * 0.5; 

    dphi(0,0) = -0.5;
    dphi(0,1) = 0.5;
}

void Geom1d::X(const VecDouble &xi, MatrixDouble &NodeCo, VecDouble &x) {
    if(xi.size() != Dimension || xi.size() != NodeCo.rows() || NodeCo.cols() != nCorners || x.rows() != NodeCo.rows()) DebugStop();
    
    int nrow = NodeCo.rows();
    x.resize(nrow);
    x.setZero();
    
    for (int i = 0; i < nrow; i++) {
        x[i] = NodeCo(i, 0)*(1 - xi[0])*0.5 + NodeCo(i, 0)*(1 + xi[0])*0.5;
    }

}

void Geom1d::GradX(const VecDouble &xi, MatrixDouble &NodeCo, VecDouble &x, MatrixDouble &gradx) {
    if(xi.size() != Dimension || xi.size() != NodeCo.rows() || NodeCo.cols() != nCorners ||  gradx.rows() != NodeCo.rows()) DebugStop();
    int nrow = NodeCo.rows();
    int ncol = NodeCo.cols();

    gradx.resize(nrow, 1);
    gradx.setZero();
    x.resize(nrow);
    x.setZero();

    VecDouble phi(2);
    MatrixDouble dphi(Dimension,2);
    Shape(xi, phi, dphi);
    for(int i = 0; i < ncol; i++) {
       for(int j = 0; j < nrow; j++) {
            x[j] += NodeCo(j,i) * phi[i];
            gradx(j,0) += NodeCo(j,i) * dphi(0,i);

        } 
    }

}

void Geom1d::SetNodes(const VecInt &nodes) {
    if(nodes.rows() != 2) DebugStop();
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
    fNeighbours[side]=neighbour;
}
