#ifndef QP_SIM_H_
#define QP_SIM_H_

#define EIGEN_RUNTIME_NO_MALLOC

#include <qpOASES.hpp>
#include <Eigen.h>
#include <math.h>
#include <MemoryFree.h>
#include <ref_t.h>
#include <ref_xr.h>
#include <ref_ur.h>
#include <ref_ddxr.h>

#define round(x) ((x)>=0?(long)((x)+0.5):(long)((x)-0.5))

#define DEBUG

#define PI 3.14159265359
#define g 9.81
#define T 0.01 //sampling time
#define dt 0.01 //sampling time
#define SAMPLING_RATE 20 //msec
// ================================= Operation Parameters ================================== //    
#define MODEL_SWITCH  2 //model, 1:nonlinear model, 2:linear model          
#define REF_TRAJ_SWITCH 2 //reference trajectory, 1:wv, 2:xy                
#define SWITCH_ZMP 1 //zmp constraint, 1:on 0:off                         
#define SWITCH_SLIP 1 //friction constraint, 1:on 0:off                    
#define SWITCH_INPUT 1 //input constraint, 1:on 0:off                      
#define SWITCH_STATE 0 //state constraint, 1:on 0:off                       

// ================================= Kinematic Parameters ================================== //    
#define D 0.5 //wheel base
#define L 0.3 // width
#define h 0.6 // height of mass center

// ================================= MPC Control Parameters ================================== //    
#define N 2 //Batch size

#define n 3 //number of states
#define m 3 //number of inputs

#define Qs 5 //state weight
#define Rs 1 //input weight

//input constraint
double input_const_mag = 2.0;
Eigen::Matrix<double, n*2, 1> thr_input;    

//state constraint
double state_const_mag = 0.008;
Eigen::Matrix<double, n*2, 1> thr_state;    

//slip constraint
Eigen::Matrix<double, n, 1> mu;    

//initial state : [x y theta]
Eigen::Matrix<double, n, 1> x0;    

double Df, Dr, Lf, Lr;

qpOASES::Options options;
qpOASES::int_t nWSR = 100;
qpOASES::returnValue rvalue;

//states//
Eigen::Matrix<double, n, 1> xr0;    
Eigen::Matrix<double, n, 1> x0_tilt;
Eigen::Matrix<double, m, 1> u_tilt;
Eigen::Matrix<double, n, 1> x_tilt_current;
Eigen::Matrix<double, m, 1> u;
Eigen::Matrix<double, n, 1> x_state;
Eigen::Matrix<double, n, 1> dx_state;
Eigen::Matrix<double, n, 1> ddx_state;
Eigen::Matrix<double, n, 1> x_tilt_m1;
Eigen::Matrix<double, n, 1> x_tilt_m1_dummy;

//batch formulation matrices//
Eigen::Matrix<double, n*N, n  > Sx;
Eigen::Matrix<double, n*N, m*N> Su;
Eigen::Matrix<double, n*N, n*N> Qb;
Eigen::Matrix<double, m*N, m*N> Rb;

Eigen::Matrix<double, m*N, m*N, Eigen::RowMajor> H_eigen;
Eigen::Matrix<double, n,   m*N, Eigen::RowMajor> F_eigen;
//Eigen::Matrix<double, n,   n  > Y_eigen;

//MPC constraints// ////qpOASES read list with rowmajor way
Eigen::Matrix<double, (n-1)*N*SWITCH_ZMP + n*N*SWITCH_STATE , m*N*(SWITCH_ZMP+SWITCH_STATE), Eigen::RowMajor> G;       
Eigen::Matrix<double, (n-1)*N*SWITCH_ZMP + n*N*SWITCH_STATE , n, Eigen::RowMajor> E;   
Eigen::Matrix<double, (n-1)*N*SWITCH_ZMP + n*N*SWITCH_STATE , 1> wl;   
Eigen::Matrix<double, (n-1)*N*SWITCH_ZMP + n*N*SWITCH_STATE , 1> wu;   
Eigen::Matrix<double, n*N , 1> wl_input;   
Eigen::Matrix<double, n*N , 1> wu_input;  

Eigen::Matrix<double, m*N, 1> Ft_x_eigen;
Eigen::Matrix<double, (n-1)*N*SWITCH_ZMP + n*N*SWITCH_STATE , 1> wl_E_x;   
Eigen::Matrix<double, (n-1)*N*SWITCH_ZMP + n*N*SWITCH_STATE , 1> wu_E_x;   

//qpOASES parameters
qpOASES::real_t* H_qp = H_eigen.data();       
qpOASES::real_t* Ft_x_qp = Ft_x_eigen.data();       
//qpOASES::real_t* Y_qp = Y_eigen.data();       
qpOASES::real_t* G_qp = G.data();
qpOASES::real_t* wl_E_x_qp = wl_E_x.data();
qpOASES::real_t* wu_E_x_qp = wu_E_x.data();
qpOASES::real_t* wl_input_qp = wl_input.data();
qpOASES::real_t* wu_input_qp = wu_input.data();

qpOASES::real_t u_tilt_qp[m];      
qpOASES::real_t cputime;


// ================================= Simulation Parameters ================================== //    
uint32_t cur_state;
uint32_t tTime;
uint8_t trigger = 0;
uint8_t init_case = 1;

// ================================= functions ================================== //    
Eigen::Matrix3d A_from_ref(double dx_ref_index, double dy_ref_index, double theta_ref_index);
Eigen::Matrix3d B_from_ref(Eigen::ArrayXd theta_ref);

void Batch_formulation(void);
Eigen::Matrix3d multipleA(u_int32_t c_state, u_int32_t f_state);

void QP_Constraints_fast(void);    

Eigen::Matrix3d A_hat(uint32_t cur_state, uint32_t final_state); 
Eigen::MatrixXd B_hat(uint32_t cur_state, uint32_t final_state);

#endif 



