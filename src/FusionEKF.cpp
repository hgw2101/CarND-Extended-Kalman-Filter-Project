#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  ekf_.x_ = VectorXd(4); //estimate vector
  ekf_.P_ = MatrixXd(4,4); //state covariance matrix, P, only need one for both laser and radar
  ekf_.F_ = MatrixXd(4,4); //state transition matrix
  ekf_.Q_ = MatrixXd(4,4); //process covariance matrix
  R_laser_ = MatrixXd(2, 2); //measurement covariance for laser
  R_radar_ = MatrixXd(3, 3); //measurement covariance for radar
  H_laser_ = MatrixXd(2, 4); //measurement matrix for laser
  Hj_ = MatrixXd(3, 4); //measurement matrix for radar (Jacobian)

  //state covariance matrix
  ekf_.P_ << 1,0,0,0,
            0,1,0,0,
            0,0,1000,0,
            0,0,0,1000;

  //state transition matrix, this will be updated with the deltaT values later
  ekf_.F_ << 1,0,1,0,
            0,1,0,1,
            0,0,1,0,
            0,0,0,1;

  //process covariance matrix, used to model acceleration related noise (noise_ax & noise_ay),
  //and dt, initialize with all 0's
  ekf_.Q_ << 0,0,0,0,
            0,0,0,0,
            0,0,0,0,
            0,0,0,0;

  // set acceleration noise (given)
  noise_ax_ = 9;
  noise_ay_ = 9;

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
              0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
              0, 0.0009, 0,
              0, 0, 0.09;

  //measurement matrix, H_laser -laser
  H_laser_ << 1,0,0,0,
              0,1,0,0;

  //measurement matrix, Hj -radar
  Hj_ << 0,0,0,0,
         0,0,0,0,
         0,0,0,0;
  // Hj_ = tools.CalculateJacobian(x_); // TODO: CAN I DO THIS in initialization???

}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    /**
    TODO:
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    // first measurement
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;

    // cout<<"2"<<endl;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      // TODO: solve the ridiculous hard multivariate equation first

    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */

      // change px and py to first measurement accordingly, leave vx and vy as 1,1 since no data for that
      ekf_.x_(0) = measurement_pack.raw_measurements_[0];
      ekf_.x_(1) = measurement_pack.raw_measurements_[1];
    }

    // cout<<"3"<<endl;

    previous_timestamp_ = measurement_pack.timestamp_;

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
   TODO:
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */

  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0; // time elapsed in seconds
  previous_timestamp_ = measurement_pack.timestamp_; // update previous_timestamp_ to current timestamp;

  // create dt variables so we can reuse in the equations below
  float dt_2 = dt * dt;
  float dt_3 = dt_2 * dt;
  float dt_4 = dt_3 * dt;

  //update the state transition matrix F according to the new elapsed time
  ekf_.F_(0,2) = dt;
  ekf_.F_(1,3) = dt;

  //update process noise covariance matrix, Q

  ekf_.Q_(0,0) = dt_4/4*noise_ax_;
  ekf_.Q_(0,2) = dt_3/2*noise_ax_;
  ekf_.Q_(1,1) = dt_4/4*noise_ay_;
  ekf_.Q_(1,3) = dt_3/2*noise_ay_;
  ekf_.Q_(2,0) = dt_3/2*noise_ax_;
  ekf_.Q_(2,2) = dt_2*noise_ax_;
  ekf_.Q_(3,1) = dt_3/2*noise_ay_;
  ekf_.Q_(3,3) = dt_2*noise_ay_;

  // cout<<"4"<<endl;
  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  // cout<<"5"<<endl;

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
    // cout<<"measurement_pack.sensor_type_ == MeasurementPackage::RADAR"<<endl;
    Hj_ = tools.CalculateJacobian(ekf_.x_);
    ekf_.H_ = Hj_;
    ekf_.R_ = R_radar_;
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);
    // cout<<"3"<<endl;
  } else {
    // cout<<"measurement_pack.sensor_type_ == MeasurementPackage::LASER"<<endl;
    // Laser updates
    ekf_.H_ = H_laser_;
    ekf_.R_ = R_laser_;
    ekf_.Update(measurement_pack.raw_measurements_);
    // cout<<"4"<<endl;
  }

  // cout<<"6"<<endl;

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
