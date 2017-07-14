#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

#include <iostream>

using namespace std;

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  x_ = F_ * x_;
  MatrixXd Ft = F_.transpose();
  P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Kalman Filter equations
  */

  VectorXd z_pred = H_ * x_;
  VectorXd y = z - z_pred;
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si; // kalman gain

  // new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Extended Kalman Filter equations
  */

  // cout<<"z_pred"<<endl;

  VectorXd z_pred = tools.ConvertCartesianToPolar(x_);
  VectorXd y = z - z_pred;

  MatrixXd Hj = H_; // H_ is set to Jacobian matrix by FusionEKF

  // cout<<"MatrixXd Hjt = Hj.transpose();"<<endl;

  MatrixXd Hjt = Hj.transpose();

  // cout<<"MatrixXd S = Hj * P_ * Hjt + R_;"<<endl;

  // cout<<"this is Hj: "<<Hj<<endl;
  // cout<<"this is P_: "<<P_<<endl;
  // cout<<"this is Hjt: "<<Hjt<<endl;
  // cout<<"this is R_: "<<R_<<endl;

  MatrixXd S = Hj * P_ * Hjt + R_;

  // cout<<"MatrixXd Si = S.inverse();"<<endl;

  MatrixXd Si = S.inverse();

  // cout<<"MatrixXd PHjt = P_ * Hjt;"<<endl;

  MatrixXd PHjt = P_ * Hjt;

  // cout<<"MatrixXd K = PHjt * Si;"<<endl;

  MatrixXd K = PHjt * Si; // kalman gain

  // cout<<"x_ = x_ + (K * y);"<<endl;
  // new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * Hj) * P_;

}
