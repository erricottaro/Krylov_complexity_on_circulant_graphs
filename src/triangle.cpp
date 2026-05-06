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

double c_plus(double theta1, double phi1, double theta2, double phi2){
    double norm = sqrt(2./3.);

    double c1 = cos(theta1/2.);
    double s1 = sin(theta1/2.); 
    double c2 = cos(theta2/2.);
    double s2 = sin(theta2/2.);

    double shift = 2*M_PI/3.;

    double c1_shifted = cos(phi1+shift);
    double c1_shifted2 = cos(phi1+2*shift);
    double c2_shifted = cos(phi2+shift);
    double c2_shifted2 = cos(phi2+2*shift);
    double cdiff_shifted = cos(phi1-phi2 + shift);
    double cdiff_shifted2 = cos(phi1-phi2 + 2*shift);

    double body = 1 + c1*s1*c2*(c1_shifted+c1_shifted2) + c1*s1*s2*(c2_shifted+c2_shifted2) + (s1*s1)*s2*c2*(cdiff_shifted+cdiff_shifted2);

    return norm*sqrt(body);
}

double a_0(double theta1, double phi1, double theta2, double phi2){
    double first = norm(c_0(theta1, phi1, theta2, phi2));
    double second = norm(c_plus(theta1, phi1, theta2, phi2));

    return 2.*first - second;
}

double b_1(double theta1, double phi1, double theta2, double phi2){
    double first = norm(c_0(theta1, phi1, theta2, phi2))*pow((2.-a_0(theta1, phi1, theta2, phi2)), 2);
    double second = norm(c_plus(theta1, phi1, theta2, phi2))*pow((1.+a_0(theta1, phi1, theta2, phi2)), 2);

    return sqrt(first+second);
}

double c_bar(double theta1, double phi1, double theta2, double phi2){
    

    return 2./9.*pow(b_1(theta1, phi1, theta2, phi2),2);    
}

