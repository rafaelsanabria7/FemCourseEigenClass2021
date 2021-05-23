/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include <iostream> 
#include "IntRuleTriangle.h"
#include "tpanic.h"

IntRuleTriangle::IntRuleTriangle(){

}

IntRuleTriangle::IntRuleTriangle(int order) {
    SetOrder(order);
}

void IntRuleTriangle::SetOrder(int order) {
    fOrder = order;
    if (order < 0 || order > MaxOrder()) DebugStop();
    switch (order)
    {
    case 0:
    case 1:
        fPoints.resize(1,2);
        fWeights.resize(1);
        fPoints(0,0) = 1./3.;
        fPoints(0,1) = 1./3.;
        fWeights[0] = 1./2.;
        break;
    
    case 2:
        fPoints.resize(3,2);
        fWeights.resize(3);

        fPoints(0,0) = 1./6.;
        fPoints(0,1) = 1./6.;

        fPoints(1,0) = 2./3.;
        fPoints(1,1) = 1./6.;

        fPoints(2,0) = 1./6.;
        fPoints(2,1) = 2./3.;

        fWeights[0] = 1./6.;
        fWeights[1] = 1./6.;
        fWeights[2] = 1./6.;
        break;

    case 3:
        fPoints.resize(4,2);
        fWeights.resize(4);

        fPoints(0,0) = 1./3.;
        fPoints(0,1) = 1./3.;

        fPoints(1,0) = 1./5.;
        fPoints(1,1) = 1./5.;

        fPoints(2,0) = 3./5.;
        fPoints(2,1) = 1./5.;
        
        fPoints(3,0) = 1./5.;
        fPoints(3,1) = 3./5.;

        fWeights[0] = -27./96.;
        fWeights[1] =  25./96.;
        fWeights[2] =  25./96.;
        fWeights[3] =  25./96.;
        break;
    
    case 4:
        fPoints.resize(6,2);
        fWeights.resize(6);

        fPoints(0,0) =  (8. - sqrt(10.) + sqrt(38. - 44. * sqrt(2. / 5.))) / 18.;
        fPoints(0,1) =  (8. - sqrt(10.) + sqrt(38. - 44. * sqrt(2. / 5.))) / 18.;

        fPoints(1,0) =  (8. - sqrt(10.) + sqrt(38. - 44. * sqrt(2. / 5.))) / 18.;
        fPoints(1,1) =  1.-2.*(8. - sqrt(10.) + sqrt(38. - 44. * sqrt(2. / 5.))) / 18.;

        fPoints(2,0) =  1.-2.*(8. - sqrt(10.) + sqrt(38. - 44. * sqrt(2. / 5.))) / 18.;
        fPoints(2,1) =  (8. - sqrt(10.) + sqrt(38. - 44. * sqrt(2. / 5.))) / 18.;
        
        fPoints(3,0) =  (8. - sqrt(10.) - sqrt(38. - 44. * sqrt(2. / 5.))) / 18.;
        fPoints(3,1) =  (8. - sqrt(10.) - sqrt(38. - 44. * sqrt(2. / 5.))) / 18.;

        fPoints(4,0) =  (8. - sqrt(10.) - sqrt(38. - 44. * sqrt(2. / 5.))) / 18.;
        fPoints(4,1) =  1.-2.*(8. - sqrt(10.) - sqrt(38. - 44. * sqrt(2. / 5.))) / 18.;

        fPoints(5,0) =  1.-2.*(8. - sqrt(10.) - sqrt(38. - 44. * sqrt(2. / 5.))) / 18.;
        fPoints(5,1) =  (8. - sqrt(10.) - sqrt(38. - 44. * sqrt(2. / 5.))) / 18.;

        fWeights[0] =  (620. + sqrt(213125. - 53320. * sqrt(10.))) / 7440.;
        fWeights[1] =  (620. + sqrt(213125. - 53320. * sqrt(10.))) / 7440.;
        fWeights[2] =  (620. + sqrt(213125. - 53320. * sqrt(10.))) / 7440.;
        fWeights[3] =  (620. - sqrt(213125. - 53320. * sqrt(10.))) / 7440.;
        fWeights[4] =  (620. - sqrt(213125. - 53320. * sqrt(10.))) / 7440.;
        fWeights[5] =  (620. - sqrt(213125. - 53320. * sqrt(10.))) / 7440.;
        break;

    case 5:
        fPoints.resize(7,2);
        fWeights.resize(7);

        fPoints(0,0) =  (6. - sqrt(15.)) / 21.;
        fPoints(0,1) =  (6. - sqrt(15.)) / 21.;

        fPoints(1,0) =  (6. - sqrt(15.)) / 21.;
        fPoints(1,1) =  1.-2.*(6. - sqrt(15.)) / 21.;

        fPoints(2,0) =  1.-2.*(6. - sqrt(15.)) / 21.;
        fPoints(2,1) =  (6. - sqrt(15.)) / 21.;
        
        fPoints(3,0) =  (6. + sqrt(15.)) / 21.;
        fPoints(3,1) =  (6. + sqrt(15.)) / 21.;

        fPoints(4,0) =  (6. + sqrt(15.)) / 21.;
        fPoints(4,1) =  1.-2.*(6. + sqrt(15.)) / 21.;

        fPoints(5,0) =  1.-2.*(6. + sqrt(15.)) / 21.;
        fPoints(5,1) =  (6. + sqrt(15.)) / 21.;

        fPoints(6,0) = 1./3.;
        fPoints(6,1) = 1./3.;

        fWeights[0] =  (155. - sqrt(15.)) / 2400.;
        fWeights[1] =  (155. - sqrt(15.)) / 2400.;
        fWeights[2] =  (155. - sqrt(15.)) / 2400.;
        fWeights[3] =  (155. + sqrt(15.)) / 2400.;
        fWeights[4] =  (155. + sqrt(15.)) / 2400.;
        fWeights[5] =  (155. + sqrt(15.)) / 2400.;
        fWeights[6] =  9./80.;
        break;

    default:
        DebugStop();
        break;
    }
}
