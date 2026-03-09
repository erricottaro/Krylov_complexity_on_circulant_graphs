#include "triangle.h"



complex<double> c_0(double theta1, double phi1, double theta2, double phi2){
    double norm = 1./sqrt(3.);

    double c1 = cos(theta1/2.);
    double s1 = sin(theta1/2.); 
    double c2 = cos(theta2/2.);
    double s2 = sin(theta2/2.);

    complex<double> phase1 = exp(phi1*I);
    complex<double> phase2 = exp(phi2*I);

    return norm*(c1 + phase1*s1*c2 + phase2*s1*s2);
}