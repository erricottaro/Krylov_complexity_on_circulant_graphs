#include "triangle.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <random>
#include <vector>
#include <armadillo>

using namespace std;

//generate a random qu-dit and save it inside vec using generator gen
void random_qudit(arma::cx_vec&, default_random_engine&);

//angle-phase coordinates of qutrit in angles expressed in cartesian coordinates in cart 
void qutrit_XYZ(arma::vec&, arma::cx_vec);

int main(){
    default_random_engine gen;

    int n_points = 100000;
    const int dim = 3;
    //state in cartesian coordinates
    arma::cx_vec z(dim);
    //state in angle-phase coordinates
    arma::vec angle_coords(2*(dim-1));

    //generate the random state, uniformly distributed on CP^n
    //random_qudit(z, gen);

    //compute its angle-phase coordinates
    //qutrit_XYZ(angle_coords, z);
    
    /*
    cout << "First vector generated: " << endl;
    //cout << "Complex vector in CP^" << dim-1 << endl;
    for (int i = 0; i < dim; i++) {
        cout << setprecision(3) << z[i] << endl;
    }
    cout << "Norm: " << arma::norm(z) << endl;

    cout << "Angle-phase coordinates " << endl;
    for (int i = 0; i < angle_coords.size(); i++) {
        cout << setprecision(3) << angle_coords[i]/M_PI << "pi" << endl;
    }

    cout << "complexity: " << c_bar(angle_coords[0], angle_coords[1], angle_coords[2], angle_coords[3]) << endl;
    */
    
    ofstream output("../output/random.out");

    output << "n_points: " << setw(7) << n_points << endl
           << setw(7) << "complexity" << endl;
    
    for (int i = 0; i < n_points; i++) {
        //generate the random state, uniformly distributed on CP^n
        random_qudit(z, gen);
        //compute its angle-phase coordinates
        qutrit_XYZ(angle_coords, z);
        //compute complexity
        double complexity = c_bar(angle_coords[0], angle_coords[1], angle_coords[2], angle_coords[3]);
        output << setw(7) << complexity << endl;
    }

    output.close();

    return 0;
}

void random_qudit(arma::cx_vec &vec, default_random_engine &gen){

    int dim = vec.size();

    double rand_real;
    double rand_im;

    double norm=0;

    normal_distribution<double> distr(0., 1.);

    //generate complex numbers according to multivariate functions
    for (int i = 0; i < dim; i++) {
        
        rand_real = distr(gen);
        rand_im   = distr(gen);

        vec[i] = rand_real + rand_im*I;;

        norm += rand_real*rand_real + rand_im*rand_im;
    }
    
    //renormalize
    vec /= sqrt(norm);

    //compute global phase
    complex<double> v0 = vec[0];

    complex<double> phase = conj(v0)/abs(v0);

    //multiply vec by global phase

    vec*=phase;
}

void qutrit_XYZ(arma::vec &angles, arma::cx_vec cart){

    //spherical coordinates
    angles[0] = 2.*acos(abs(cart[0]));
    angles[1] = 2.*atan2(abs(cart[2]), abs(cart[1]));

    //phases
    angles[2] = arg(cart[1]);
    angles[3] = arg(cart[2]);
}