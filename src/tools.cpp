#include "tools.h"
#include <cmath>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
   * TODO: Calculate the RMSE here.
   */
	VectorXd rmse(4);
   rmse << 0, 0, 0, 0;
	
	for (unsigned int i=0; i < estimations.size(); ++i) {
    VectorXd residual = estimations[i] - ground_truth[i];
    residual = residual.array()*residual.array();
    rmse += residual;
  }
  rmse = rmse/estimations.size();
  rmse = rmse.array().sqrt();
  return rmse;

}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
   * TODO:
   * Calculate a Jacobian here.
   */
	float px = x_state[0];
	float py = x_state[1];
	float vx = x_state[2];
	float vy = x_state[3];
	
	float pxpy2 = pow(px, 2) + pow(py, 2);
	if (pxpy2<0.00001){
		px = 0.001;
		py = 0.001;
		pxpy2 = pow(px, 2) + pow(py, 2);
    }
	float pxpy1 = sqrt(pxpy2);
	float pxpy3 = pxpy1*pxpy2;

	MatrixXd Hj = MatrixXd(3, 4);
	Hj << px / pxpy1, py / pxpy1, 0, 0,
		  -py / pxpy2, px / pxpy2, 0, 0,
		  py*(vx*py - vy * px) / pxpy3, px*(vy * px - vx * py) / pxpy3, px / pxpy1, py / pxpy1;
	return Hj;
}
