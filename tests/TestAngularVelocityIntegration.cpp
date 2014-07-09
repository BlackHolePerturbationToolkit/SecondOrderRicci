#include <vector>
#include <iostream>
#include <iomanip>
#include "Quaternions.hpp"
#include "IntegrateAngularVelocity.hpp"
using namespace std;
using namespace Quaternions;

const double alpha = 0.05;
const double T_i = 0.0;
const double T_f = 20000.0;

Quaternion R_1func(const double t) {
  return Quaternions::exp(alpha*t*Quaternions::xHat/2.)
    *Quaternions::exp((alpha*pow(t,1.2))*Quaternions::yHat/2.)
    *Quaternions::exp((alpha*pow(t,1.3))*Quaternions::zHat/2.);
}

vector<double> Omega(const double t) {
  const Quaternion R_1 = R_1func(t);
  const Quaternion Rdot_1 = (alpha/2.0) *
    (
     Quaternions::xHat*
     Quaternions::exp(alpha*t*Quaternions::xHat/2.)
     *Quaternions::exp((alpha*pow(t,1.2))*Quaternions::yHat/2.)
     *Quaternions::exp((alpha*pow(t,1.3))*Quaternions::zHat/2.)
     +
     Quaternions::exp(alpha*t*Quaternions::xHat/2.)
     *1.2*pow(t,0.2)*Quaternions::yHat
     *Quaternions::exp((alpha*pow(t,1.2))*Quaternions::yHat/2.)
     *Quaternions::exp((alpha*pow(t,1.3))*Quaternions::zHat/2.)
     +
     Quaternions::exp(alpha*t*Quaternions::xHat/2.)
     *Quaternions::exp((alpha*pow(t,1.2))*Quaternions::yHat/2.)
     *1.3*pow(t,0.3)*Quaternions::zHat
     *Quaternions::exp((alpha*pow(t,1.3))*Quaternions::zHat/2.)
     );
  return (2*Rdot_1*Quaternions::conjugate(R_1)).vec();
}

int main() {
  vector<Quaternion> Qs;
  vector<double> Ts;
  cout << setprecision(16);
  FrameFromAngularVelocity(&Omega, T_i, T_f, Qs, Ts);
  for(unsigned int i=0; i<Ts.size(); ++i) {
    cout << Ts[i]
  	 // << "\n\t" << Qs[i]
  	 // << "\t" << Qs[i].angle()
  	 // << "\n\t" << R_1func(Ts[i])
  	 // << "\n\t" << Qs[i]-R_1func(Ts[i])
  	 << "\t" << (Qs[i]-R_1func(Ts[i])).abs()
  	 // << "\t" << 1.e-14*std::sqrt(i)
  	 << endl;
  }
  cout << "Took " << Ts.size() << " steps." << endl;

  return 0;
}
