#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  VectorXd rmse(4);
  rmse<<0,0,0,0; // initialize rmse vector with all 0's

  // check validity of input data
  if(estimations.size() != ground_truth.size() || estimations.size() == 0) {
    cout<<"Invalid estimation or ground_truth data"<<endl;
    return
  }

  // iterate through the estimations and ground_truth arrays
  for(unsigned int i=0; i < estimations.size(); ++i) {

    // calculate residual
    VectorXd residual = estimations[i] - ground_truth[i];
    residual = residual.array()*residual.array(); // square it

    rmse += residual;
  }

  // calculate the mean
  rmse = rmse/estimations.size();

  // take square root
  rmse = rmse.array().sqrt();

  return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {

  MatrixXd Hj(3,4);

  //recover state parameters
  float px = x_state(0);
  float py = x_state(1);
  float vx = x_state(2);
  float vy = x_state(3);

  //preset variables to avoid repeated caculations
  float c1 = px*px+py*py;
  float c2 = sqrt(c1);
  float c3 = (c1*c2);

  //check division by 0, check absolute value of c1
  if(fabs(c1) < 0.0001) {
    cout<<"CalculateJacobian Error: Division by 0"<<endl;
    return Hj;
  }

  //compute Jacobian matrix
  Hj << (px/c2), (py/c2), 0, 0,
      -(py/c1), (px/c1), 0, 0,
      py*(vx*py - vy*px)/c3, px*(px*vy - py*vx)/c3, px/c2, py/c2;

  return Hj;
}
