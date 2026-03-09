#include "triangle.h"
#include <iostream>
#include <fstream>

using namespace std;

int main(){

    //cout << I << endl;

    double theta_max = M_PI;
    int n_theta = 1000;

    double phi_max = 2*M_PI;
    int n_phi = 1000;

    double delta_theta = theta_max/n_theta;
    double delta_phi = phi_max/n_phi;

    for(int i=0; i<n_theta; i++){
        double theta = i*delta_theta;
        for(int j=0; i<n_phi; j++){
            double phi = j*delta_phi;
            double complexity = c_bar(theta, phi);
        }
        //cout << theta << endl;
        //double c = c_bar(theta, M_PI/2.);
    }

    return 0;
}