#ifndef QP_SIM_H_
#define QP_SIM_H_

#include <qpOASES.hpp>
#include <Eigen.h>
#include <math.h>
#include <MemoryFree.h>

#define DEBUG
#define PI 3.14159265359
#define g 9.81
#define T 0.01 //sampling time
#define dt 0.01 //sampling time
/* ================================= Operation Parameters ================================== */    
#define MODEL_SWITCH  2 //model, 1:nonlinear model, 2:linear model          
#define REF_TRAJ_SWITCH 2 //reference trajectory, 1:wv, 2:xy                
#define SWITCH_ZMP 1 //zmp constraint, 1:on 0:off                         
#define SWITCH_SlIP 0 //friction constraint, 1:on 0:off                    
#define SWTICH_INPUT 1 //input constraint, 1:on 0:off                      
#define SWITCH_STATE 0 //state constraint, 1:on 0:off                       

/* ================================= Kinematic Parameters ================================== */    
#define D 0.5 //wheel base
#define L 0.3 // width
#define h 0.6 // height of mass center

/* ================================= MPC Control Parameters ================================== */    
#define n 3 //number of states
#define m 3 //number of inputs
#define N 10 //Batch size
#define Qs 5 //state weight
#define Rs 1 //input weight

/* ================================= functions ================================== */    
void ref_trajectory_mecanum_xy(void);
inline double constrainAngle(double x);
inline double angleConv(double angle);
inline double angleDiff(double a,double b);
inline double unwrap(double previousAngle,double newAngle);

Eigen::Matrix3d A_from_ref(double dx_ref_index, double dy_ref_index, double theta_ref_index);
Eigen::Matrix3d B_from_ref(Eigen::ArrayXd theta_ref);

void Batch_formulation(void);
Eigen::Matrix3d multipleA(u_int32_t c_state, u_int32_t f_state);


/* ================================= MPC Control Parameters ================================== */    
uint32_t tTime;
uint8_t trigger = 0;
uint8_t init_case = 1;

//input constraint
float input_const_mag = 2.0;
Eigen::VectorXd thr_input(6);    

//state constraint
float state_const_mag = 0.008;
Eigen::VectorXd thr_state(6);    

//slip constraint
Eigen::VectorXd mu(n);    

//initial state : [x y theta]
Eigen::VectorXd x0(n);    

qpOASES::Options options;
qpOASES::int_t nWSR = 100;

qpOASES::real_t xOpt[n];



/* ================================= Simulation Parameters ================================== */    
Eigen::ArrayXd t;  
Eigen::ArrayXd x_ref, y_ref, theta_ref;
Eigen::ArrayXd dx_ref, dy_ref, dtheta_ref;
Eigen::ArrayXd ddx_ref, ddy_ref, ddtheta_ref;
Eigen::ArrayXXd xr(1301, n);
Eigen::ArrayXXd ur(1301, m);
Eigen::ArrayXXd ddxr(1301, n);

Eigen::MatrixXd A(n,n);
Eigen::MatrixXd B(n,m);

Eigen::VectorXd x0_tilt(n);
Eigen::VectorXd u_tilt(m);
Eigen::VectorXd x_tilt_current(n);
Eigen::VectorXd u(m);
Eigen::VectorXd x_state(n);
Eigen::VectorXd dx_state(n);
Eigen::VectorXd ddx_state(n);

uint32_t cur_state;

Eigen::MatrixXd Sx(n*N,n);
Eigen::MatrixXd Su(n*N,m*N);
Eigen::MatrixXd Qb(n*N, n*N);
Eigen::MatrixXd Rb(m*N, m*N);


//eigen library

Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> H_eigen(m*N,m*N);
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> F_eigen(n,m*N);
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> Y_eigen(n,n);

//Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> G_eigen(2,1);
//Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> wl_input_eigen(2,1);
//Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> wu_input_eigen(2,1);
//Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> wl_eigen(1,1);
//Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> wu_eigen(1,1);

qpOASES::real_t* H_qp = H_eigen.data();       
qpOASES::real_t* F_qp = F_eigen.data();       
qpOASES::real_t* Y_qp = Y_eigen.data();       
//qpOASES::real_t* G_qp = G_eigen.data();       
//qpOASES::real_t* wl_input_qp = wl_input_eigen.data();       
//qpOASES::real_t* wu_input_qp = wu_input_eigen.data();       
//qpOASES::real_t* wl_qp = wl_eigen.data();       
//qpOASES::real_t* wu_qp = wu_eigen.data();       

#endif // TURTLEBOT3_WITH_OPEN_MANIPULATOR_CONFIG_H_



