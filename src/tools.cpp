#define _USE_MATH_DEFINES

#include <iostream>
#include "tools.h"
#include <algorithm>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;
using std::min;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */
  VectorXd rmse(4);
	rmse << 0,0,0,0;

	// check the validity of the following inputs:
	//  * the estimation vector size should not be zero
	//  * the estimation vector size should equal ground truth vector size
	if(estimations.size() != ground_truth.size()
			|| estimations.size() == 0){
		cout << "Invalid estimation or ground_truth data" << endl;
		return rmse;
	}

	//accumulate squared residuals
	for(unsigned int i=0; i < estimations.size(); ++i){

		VectorXd residual = estimations[i] - ground_truth[i];

		//coefficient-wise multiplication
		residual = residual.array()*residual.array();
		rmse += residual;
	}

	//calculate the mean
	rmse = rmse/estimations.size();

	//calculate the squared root
	rmse = rmse.array().sqrt();

	//return the result
	return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
  TODO:
    * Calculate a Jacobian here.
  */
  float px = x_state[0];
  float py = x_state[1];
  float vx = x_state[2];
  float vy = x_state[3];
  float eps = 1e-6;

  MatrixXd Hj_(3,4);
  Hj_ << 0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0;
  float pxpy_2_1 = px * px + py * py;
  float pxpy_1_2 = sqrt(pxpy_2_1);
  float pxpy_3_2 = pxpy_2_1 * pxpy_1_2;
  float pxpy_2_1_eps = max(pxpy_2_1 , eps);
  float pxpy_1_2_eps = max(pxpy_1_2 , eps);
  float pxpy_3_2_eps = max(pxpy_3_2 , eps);

  // //check division by zero
  // if(pxpy_2_1 < 0.0001){
  //   cout << "CalculateJacobian () - Error - Division by Zero" << endl;
  //   return Hj_;
  // }

  Hj_(0,0) = px / pxpy_1_2_eps;
  Hj_(0,1) = py / pxpy_1_2_eps;
  Hj_(0,2) = 0;
  Hj_(0,3) = 0;

  Hj_(1,0) = -py / pxpy_2_1_eps;
  Hj_(1,1) =  px / pxpy_2_1_eps;
  Hj_(1,2) = 0;
  Hj_(1,3) = 0;

  float temp = vx * py - vy * px;
  Hj_(2,0) = py * temp / pxpy_3_2_eps;
  Hj_(2,1) = -px * temp / pxpy_3_2_eps;
  Hj_(2,2) = Hj_(0,0);
  Hj_(2,3) = Hj_(0,1);

  return Hj_;
}

VectorXd Tools::CalculateNonlinear_h(const VectorXd& x_state) {
  /**
  TODO:
    * Calculate h(x)
  */
  float px = x_state[0];
  float py = x_state[1];
  float vx = x_state[2];
  float vy = x_state[3];
  float eps = 1e-6;
  float pxpy_1_2 = sqrt(px * px + py * py);
  VectorXd h_ =  VectorXd(3);
  h_[0] = pxpy_1_2;
  h_[1] = atan2(py, px);
  h_[2] = (px * vx + py * vy)/ max(pxpy_1_2, eps);

  return h_;
}
