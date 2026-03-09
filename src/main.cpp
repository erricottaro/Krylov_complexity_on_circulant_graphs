#include "triangle.h"
#include <iostream>
#include <fstream>
#include <iomanip>

using namespace std;

int main(){

    //cout << I << endl;

    double theta_max = M_PI;
    int n_theta = 100;

    double phi_max = 2*M_PI;
    int n_phi = 100;

    double delta_theta = theta_max/n_theta;
    double delta_phi = phi_max/n_phi;

    ofstream output;

    output.open("../output/output.out");

    output << setw(20) << "theta1" << setw(20) << "phi1" << setw(20) << "theta2" << setw(20) << "phi2" << setw(20) << "C_bar" << endl;

    /*
    for(int i=0; i<n_theta; i++){
        double theta1 = i*delta_theta;
        for(int j=0; j<n_phi; j++){
            double phi1 = j*delta_phi;
            c_bar(theta1, phi1);
            output << setw(20) << theta1 << setw(20) << phi1 << setw(20) << 0. << setw(20) << 0. << endl;
        }
    }
    */
    
    for(int i=0; i<n_theta; i++){
        double theta1 = i*delta_theta;
        for(int j=0; j<n_phi; j++){
            double phi1 = j*delta_phi;
            for(int k=0; k<n_theta; k++){
                double theta2 = k*delta_theta;
                for(int l=0; l<n_phi; l++){
                    double phi2 = l*delta_phi;
                    double c_b = c_bar(theta1, phi1, theta2, phi2);
                    output << setw(20) << theta1 << setw(20) << phi1 << setw(20) << theta2 << setw(20) << phi2 << setw(20) << c_b << endl;
                }
            }
        }
    }
    

    output.close();

    return 0;
}
