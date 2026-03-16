#include "triangle.h"
#include <iostream>
#include <fstream>
#include <iomanip>

using namespace std;

int main(){

    //cout << I << endl;

    double theta_max = M_PI;
    int n_theta1 = 100;

    double phi_max = 2*M_PI;
    int n_phi1 = 100;

    int n_theta2 = 5;
    int n_phi2 = 5;

    double delta_theta1 = theta_max/n_theta1;
    double delta_phi1 = phi_max/n_phi1;

    double delta_theta2 = theta_max/n_theta2;
    double delta_phi2 = phi_max/n_phi2;

    ofstream output;

    output.open("../output/output.out");

    output << setw(20) << "theta1" << setw(20) << "phi1" << setw(20) << "theta2" << setw(20) << "phi2" << setw(20) << "C_bar" << endl;

    /*
    for(int i=0; i<n_theta1; i++){
        double theta1 = i*delta_theta1;
        for(int j=0; j<n_phi1; j++){
            double phi1 = j*delta_phi1;
            c_bar(theta1, phi1);
            output << setw(20) << theta1 << setw(20) << phi1 << setw(20) << 0. << setw(20) << 0. << endl;
        }
    }
    */
    
    for(int i=0; i<n_theta1; i++){
        double theta1 = i*delta_theta1;
        for(int j=0; j<n_phi1; j++){
            double phi1 = j*delta_phi1;
            for(int k=0; k<n_theta2; k++){
                double theta2 = k*delta_theta2;
                for(int l=0; l<n_phi2; l++){
                    double phi2 = l*delta_phi2;
                    double c_b = c_bar(theta1, phi1, theta2, phi2);
                    output << setw(20) << theta1 << setw(20) << phi1 << setw(20) << theta2 << setw(20) << phi2 << setw(20) << c_b << endl;
                }
            }
        }
    }
    

    output.close();

    return 0;
}
