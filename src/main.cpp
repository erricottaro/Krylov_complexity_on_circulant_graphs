#include "triangle.h"
#include <iostream>

using namespace std;

int main(){

    //cout << I << endl;

    double theta_max = M_PI;
    int n_theta = 10;

    double delta_theta = theta_max/n_theta;

    for(int i=0; i<n_theta; i++){
        double theta = i*delta_theta;
        //cout << theta << endl;
        //cout << c_0(theta, M_PI/2.) << endl;
    }

    return 0;
}