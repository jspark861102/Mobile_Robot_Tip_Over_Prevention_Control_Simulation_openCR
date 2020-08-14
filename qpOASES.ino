#include <qp_sim.h>

// ============================================== setup ============================================================ //        
// ---------------------------------------------------------------------------------------------------------------- //
void setup() {    
    Serial.begin(115200);
    while(!Serial) {}
    pinMode(LED_BUILTIN, OUTPUT); 

    // ================================= parameters setting ================================== //
    thr_input.fill(input_const_mag);
    thr_state.fill(state_const_mag);
    mu << 0.46, 0.46, 0;
    //mu << 0.1, 0.1, 0;    
    
    if (REF_TRAJ_SWITCH == 1) {
        x0 << 0.1, 0.2, PI/8;
    }
    else if (REF_TRAJ_SWITCH == 2) {
        x0 << 0.0, 0.0, 0.0;
    }
    // --------------------------------------------------------------------------------------- //


    // ================================= Reference trajectory ================================== //
    //reference trajectory xr=[x_ref, y_ref, theta_ref], ur=[dx_ref, dy_ref, dtheta_ref]//
    if (REF_TRAJ_SWITCH == 1){        
    }
    else if (REF_TRAJ_SWITCH == 2) {
        //ref_trajectory_mecanum_xy(); //xr, ur, ddxr
    }    
    // --------------------------------------------------------------------------------------- //

    // ================================= initial parameters ================================== //
    Eigen::Vector3d xr0;    
    xr0 << pgm_read_float_far(&xr[0][0])/10000, pgm_read_float_far(&xr[0][1])/10000, pgm_read_float_far(&xr[0][2])/10000;
    x0_tilt = x0 - xr0;    
    x_tilt_current << x0_tilt;

    u_tilt.setZero(n,1);
    u.setZero(m,1);
    x_state.setZero(n,1);
    dx_state.setZero(n,1);
    ddx_state.setZero(n,1);
    x_tilt_m1.setZero(n,1);
    x_tilt_m1_dummy.setZero(n,1);

    cur_state = 853;
    // --------------------------------------------------------------------------------------- //    
    uint32_t t_b_start = micros();      
    Batch_formulation();
    uint32_t t_b_end = micros();
    uint32_t t_b = t_b_end - t_b_start;    
    Serial.print("ellapsed time of QP in setup = ");  Serial.print(t_b);  Serial.println("usec");
    Serial.println(" ");        
    
    //u_tilt_set, x_tilt_set    
    //DF,Dr,Lf,Lr --> let's use as a scalar D and L    
}

// ============================================== loop ============================================================ //        
// ================================================================================================================ //        
void loop() {    
    //control time
    uint32_t t_millis = millis();    
    
    // ================================= qpOASES Definition ================================== //
    //uint32_t t_q_start = micros();          
    static qpOASES::SQProblem qptest(m*N, (n-1)*N);    //define only once       
    //qptest.setOptions( options );   
    //uint32_t t_q_end = micros();
    //uint32_t t_q = t_q_end - t_q_start;
    //Serial.println(t_q);
    // --------------------------------------------------------------------------------------- //       

    if ((t_millis-tTime) >= 1000 ) { //control with sampling time 1000msec
        uint32_t t_QP_start = micros();      
    
    // ==================================== Batch formulation ================================= //           
    Batch_formulation(); //define Sx, Su, Qb, Rb, H, F, Y//
    // --------------------------------------------------------------------------------------- //       


    // ============================= -1 states for zmp constraint ============================ //           
    if (cur_state == 0 || cur_state == 1) x_tilt_m1 << x0_tilt;
    else                                  x_tilt_m1 << x_tilt_m1_dummy;    
    // --------------------------------------------------------------------------------------- //       
    
    // =================== if object is slipped, zmp bound is fluctuated===================== //
    if (SWITCH_SlIP == 0) { //when the slip is not considered
        if (ddx_state[0] > mu[0]*g) {
            Df = D * 0.5;
            Dr = D * 1.5;
        }
        else if (ddx_state[0] < -mu[0]*g) {
            Df = D * 1.5;
            Dr = D * 0.5;
        }

        if (ddx_state[1] > mu[1]*g) {
            Lf = L * 0.2;
            Lr = L * 1.8;
        }
        else if (ddx_state[1] < -mu[1]*g) {
            Lf = L * 1.8;
            Lr = L * 0.2;
        }
    }
    // ------------------------------------------------------------------------------------ //       

    // ================================= QP constraints =================================== //           
    QP_Constraints_fast(); //define G, wl_input, wu_input, wl_eigen, wu_eigen//
    // --------------------------------------------------------------------------------------- //       











        if (cur_state >= sizeof(t)/sizeof(float)) cur_state = 0; //reset when the simulation is finished
        else cur_state = ++cur_state;
        // --------------------------------------------------------------------------------------- //

        uint32_t t_QP_end = micros();
        uint32_t t_QP = t_QP_end - t_QP_start;
        Serial.print("ellapsed time of QP = ");  Serial.print(t_QP);  Serial.println("usec");
        Serial.println(" ");        
        
        // ================================= LED Check ================================== //
        if (trigger == 1) {digitalWrite(LED_BUILTIN, HIGH); trigger = 0;}
        else {digitalWrite(LED_BUILTIN, LOW); trigger = 1;}       
        // --------------------------------------------------------------------------------------- //

        tTime = t_millis;        
    }       
}

Eigen::Matrix3d A_hat(uint32_t cur_state, uint32_t final_state) {
    Eigen::Matrix3d result_Ahat; 
    
    result_Ahat = multipleA(cur_state, final_state-1);
    return result_Ahat;
}

Eigen::MatrixXd B_hat(uint32_t cur_state, uint32_t final_state) {
    uint32_t  del_fc;
    del_fc = final_state - cur_state;

    Serial.println(del_fc);

    Eigen::MatrixXd result_Bhat;     
    result_Bhat.setZero(n, m * del_fc);    
    
    //because 3dimensional matrix is unavailble with eigen library, (n*del_fc, m) is used instead of (n, m, del_fc)//
    for (int ii = del_fc; ii > 0; --ii)   { 
        if (ii == del_fc) result_Bhat.block(0, (ii-1)*m, n, m) = B_from_ref(pgm_read_float_far(&xr[final_state -1][2])/10000);  
        else result_Bhat.block(0, (ii-1)*n, n, m) = multipleA(cur_state+ii, final_state-1) * B_from_ref(pgm_read_float_far(&xr[cur_state+ii-1][2])/10000);  
    }

    Serial.println(result_Bhat.size());
    return result_Bhat;
}

void QP_Constraints_fast(void) {    
    Eigen::Matrix<double, n-1, n> P;
    P << 1/(T*T),         0, 0, 
               0,   1/(T*T), 0; 
    // ============================ backward finite difference for zmp ========================= //
    if (SWITCH_ZMP == 1) {

        //Eigen::Matrix<double, n, m*N> test;
        Eigen::MatrixXd test;
        test = B_hat(cur_state, cur_state + N);

        Serial.println(test.size());


    /*
        for (int ii = 0 ; ii < n; ++ii) {
            for (int jj = 0 ; jj < m*N; ++jj) {
                if (jj == m*N - 1)
                    Serial.println(test(ii,jj));
                else
                    Serial.print(test(ii,jj)); Serial.print(" ");    
            }
        }
        Serial.println(" ");
*/


        /*
        //G0 = {(n-1)*N, m*N}//     
        Eigen::Matrix<double, (n-1)*N, m*N> G0;
        
        Eigen::Matrix<double, n*N, m*N> Bhat; //(n ,m*N) is appeded through row line
        Bhat.setZero(n*N, m*N);
        for (int b_hat_num = 1; b_hat_num < N+1; ++b_hat_num)
            Bhat.block((b_hat_num-1)*n, 0, n, m*b_hat_num) = B_hat(cur_state, cur_state + b_hat_num);



        for (int ii = 0 ; ii < n*N; ++ii) {
            for (int jj = 0 ; jj < m*N; ++jj) {
                if (jj == m*N - 1)
                    Serial.println(Bhat(ii,jj));
                else
                    Serial.print(Bhat(ii,jj)); Serial.print(" ");    
            }
        }
        Serial.println(" ");
        */

        /*
        for (int ii = 1; ii < N+1; ++ii)   { 
            for (int jj = 1; jj < N+1; ++jj)   { 
                if (ii-jj == 0) //diagonal                    
                    G0.block((ii-1)*(n-1), (jj-1)*m, n-1, m) = P * Bhat.block((ii-1)*n, (jj-1)*m, n, m);                     
                else if (ii-jj == 1) //one lower than diagonal//
                    G0.block((ii-1)*(n-1), (jj-1)*m, n-1, m) = P * Bhat.block((ii-1)*n, (jj-1)*m, n, m) + P * -2*Bhat.block((ii-1-1)*n, (jj-1)*m, n, m);
                else if (ii-jj < 0) //upper diagonal is 0//     
                    G0.block((ii-1)*(n-1), (jj-1)*m, n-1, m) = Eigen::MatrixXd::Zero(n-1,m);
                else //lower diagonal//      
                    G0.block((ii-1)*(n-1), (jj-1)*m, n-1, m) = P * Bhat.block((ii-1)*n, (jj-1)*m, n, m) + P * (-2*Bhat.block((ii-1-1)*n, (jj-1)*m, n, m) + Bhat.block((ii-1-2)*n, (jj-1)*m, n, m));
            }        
        }


        

        for (int ii = 0 ; ii < (n-1)*N; ++ii) {
            for (int jj = 0 ; jj < m*N; ++jj) {
                if (jj == m*N - 1)
                    Serial.println(G0(ii,jj));
                else
                    Serial.print(G0(ii,jj)); Serial.print(" ");    
            }
        }
        Serial.println(" ");
        */
    


        /*
        //Bhat = np.zeros((N, n, m, N))          
        for b_hat_num in generator(1, N+1):    
            #in writing Bi(k+j), in code Bj(k+i), 
            #(1,2,3,4), 1:j, (2,3):B(:,:), 4:b_hat_num
            Bhat[b_hat_num-1, :, :, :b_hat_num] = hat.B_hat(A,B,cur_state, cur_state+b_hat_num)            

        G01 = np.zeros(((n-1)*N, m*N))
        for i in generator(1, N+1):    
            for j in generator(1, N+1):    
                if i-j == 0: #diagonal                    
                    G01[(i-1)*(n-1) : i*(n-1), (j-1)*m : j*m] = np.matmul( P, Bhat[i-1,:,:,j-1] )                    
                     
                elif i-j == 1: #one lower than diagonal                     
                    G01[(i-1)*(n-1) : i*(n-1), (j-1)*m : j*m] = np.matmul( P, Bhat[i-1,:,:,j-1] ) + np.matmul( P, -2*Bhat[i-1-1,:,:,j-1] )

                elif i-j < 0: #upper diagonal is 0                    
                    G01[(i-1)*(n-1) : i*(n-1), (j-1)*m : j*m] = 0

                else: #lower diagonal                    
                    G01[(i-1)*(n-1) : i*(n-1), (j-1)*m : j*m] = np.matmul(P, Bhat[i-1,:,:,j-1]) + np.matmul(P, -2*Bhat[i-1-1,:,:,j-1] + Bhat[i-1-2,:,:,j-1])                      
                    
        #G0 = np.append(G01,-G01, axis=0)
        G0 = G01.copy()
        */


    }

    else {        

    }















    /*
    for (int ii = 0 ; ii < n-1; ++ii) {
        for (int jj = 0 ; jj < n; ++jj) {
            if (jj == n - 1)
                Serial.println(P(ii,jj));
            else
                Serial.print(P(ii,jj)); Serial.print(" ");    
        }
    }
    */





}

void Batch_formulation(void) {    
    // ================================= Sx{n*N,n} ================================== //        
    Sx.setZero(n*N,n);    
    Sx.block(0, 0, 3, 3) = A_from_ref(pgm_read_float_far(&ur[cur_state][0])/10000, pgm_read_float_far(&ur[cur_state][1])/10000, pgm_read_float_far(&xr[cur_state][2])/10000);       
    for (int kk = 1; kk < N; ++kk)    
        Sx.block(n*kk, 0, n, n) = multipleA(cur_state, cur_state + kk);          
    // --------------------------------------------------------------------------------------- //

    // ================================= Su{n*N,m*N} ================================== //        
    Su.setZero(n*N,m*N);
    
    Eigen::Matrix<double, n,   n  > At;
    for (int ii = 1; ii < N+1; ++ii)   { 
        for (int jj = 1; jj < N+1; ++jj) {    
            if (ii - jj <0) //upper than diagonal 
                At.setZero(n,n);
            else if (ii -jj == 0) //diagonal
                At.setIdentity(n,n);
            else //ii-jj>0
                At = multipleA(cur_state+1, cur_state + (ii-jj));            
            Su.block((ii-1)*n, (jj-1)*m, n, m) = At * B_from_ref(pgm_read_float_far(&xr[cur_state + (jj-1)][2])/10000);                
        }
    }    
    // --------------------------------------------------------------------------------------- //

    // ================================= Qb{n*N,n*N} ================================== //       
    Qb.setZero(n*N, n*N);

    for (int ii = 0; ii < N; ++ii)    
        Qb.block(ii*n, ii*n, n, n) = Qs* Eigen::MatrixXd::Identity(n,n);                
    // --------------------------------------------------------------------------------------- //

    // ================================= Rb{m*N,m*N} ================================== //           
    Rb.setZero(m*N, m*N);

    for (int ii = 0; ii < N; ++ii)    
        Rb.block(ii*m, ii*m, m, m) = Rs*Eigen::MatrixXd::Identity(m,m);                
    // --------------------------------------------------------------------------------------- //
    
    // ================================= H(m*N,m*N) ================================== //   
    H_eigen.setZero(m*N, m*N);
    
    H_eigen = Su.transpose() * Qb * Su + Rb; 
    // --------------------------------------------------------------------------------------- //

    // ================================= F(n,m*N) ================================== //   
    F_eigen.setZero(n, m*N);
    
    F_eigen = Sx.transpose() * Qb * Su; 
    // --------------------------------------------------------------------------------------- //

    // ================================= Y(n,n) ================================== //   
    //Y_eigen.setZero(n, n);
    
    //Y_eigen = Sx.transpose() * Qb * Sx; 
    // --------------------------------------------------------------------------------------- //
    
    // ================================= print ================================== //       
    /*
    for (int ii = 0 ; ii < n*N; ++ii) {
        for (int jj = 0 ; jj < n; ++jj) {
            if (jj == n - 1)
                Serial.println(Sx(ii,jj));
            else
                Serial.print(Sx(ii,jj)); Serial.print(" ");    
        }
    }
    Serial.println(" ");

    
    for (int ii = 0 ; ii < n*N; ++ii) {
        for (int jj = 0 ; jj < m*N; ++jj) {
            if (jj == m*N - 1)
                Serial.println(Su(ii,jj));
            else
                Serial.print(Su(ii,jj)); Serial.print(" ");    
        }
    }
    Serial.println(" ");
    
    
    for (int ii = 0 ; ii < n*N; ++ii) {
        for (int jj = 0 ; jj < n*N; ++jj) {
            if (jj == n*N - 1)
                Serial.println(Qb(ii,jj));
            else
                Serial.print(Qb(ii,jj)); Serial.print(" ");    
        }
    }
    Serial.println(" ");


    for (int ii = 0 ; ii < m*N; ++ii) {
        for (int jj = 0 ; jj < m*N; ++jj) {
            if (jj == m*N - 1)
                Serial.println(Rb(ii,jj));
            else
                Serial.print(Rb(ii,jj)); Serial.print(" ");    
        }
    }
    Serial.println(" ");

        
    for (int ii = 0 ; ii < m*N; ++ii) {
        for (int jj = 0 ; jj < m*N; ++jj) {
            if (jj == m*N - 1)
                Serial.println(H_eigen(ii,jj));
            else
                Serial.print(H_eigen(ii,jj)); Serial.print(" ");    
        }
    }
    Serial.println(" ");

    for (int ii = 0 ; ii < n; ++ii) {
        for (int jj = 0 ; jj < m*N; ++jj) {
            if (jj == m*N - 1)
                Serial.println(F_eigen(ii,jj));
            else
                Serial.print(F_eigen(ii,jj)); Serial.print(" ");    
        }
    }
    Serial.println(" ");

    //for (int ii = 0 ; ii < n; ++ii) {
    //    for (int jj = 0 ; jj < n; ++jj) {
    //        if (jj == n - 1)
    //            Serial.println(Y_eigen(ii,jj));
    //        else
    //            Serial.print(Y_eigen(ii,jj)); Serial.print(" ");    
    //    }
    //}
    Serial.println(" ");   
    */
    // --------------------------------------------------------------------------------------- //
    
}

// ================================= multipleA ============================================================ //        
Eigen::Matrix3d multipleA(u_int32_t c_state, u_int32_t f_state) {
    Eigen::Matrix3d multiple_A; 
    
    multiple_A.setIdentity(n,n);
    for (int kk = c_state; kk <f_state+1; ++kk) {            
        //multiple_A = A_from_ref(ur(kk,0), ur(kk,1), xr(kk,2)) * multiple_A;                        
        multiple_A = A_from_ref(pgm_read_float_far(&ur[kk][0])/10000, pgm_read_float_far(&ur[kk][1])/10000, pgm_read_float_far(&xr[kk][2])/10000) * multiple_A;                        
    }
    return multiple_A;
}

// ================================= unwrap ============================================================ //        
//Normalize to [-180,180):
inline double constrainAngle(double x){
    x = fmod(x + PI,2*PI);
    if (x < 0)
        x += 2*PI;
    return x - PI;
}
// convert to [-360,360]
inline double angleConv(double angle){
    return fmod(constrainAngle(angle),2*PI);
}
inline double angleDiff(double a,double b){
    double dif = fmod(b - a + PI,2*PI);
    if (dif < 0)
        dif += 2*PI;
    return dif - PI;
}
inline double unwrap(double previousAngle,double newAngle){
    return previousAngle - angleDiff(newAngle,angleConv(previousAngle));
}

// ================================= System matrix ============================================================ //        
//for i in range(0, len(t)):    
//    A[:,:,i]= np.eye(3) + np.array([ [0, 0, (-ur[0][i]*np.sin(xr[2][i]) - ur[1][i]*np.cos(xr[2][i]))*T],
//                                    [0, 0, ( ur[0][i]*np.cos(xr[2][i]) - ur[1][i]*np.sin(xr[2][i]))*T], 
//                                    [0, 0,                                                          0] ])
//    B[:,:,i]= np.array([ [np.cos(xr[2][i])*T, -np.sin(xr[2][i])*T, 0],
//                        [np.sin(xr[2][i])*T,  np.cos(xr[2][i])*T, 0],
//                        [0,                   0,                  T] ])
// ----------------------------------------------------------------------------------------------------------- //
Eigen::Matrix3d A_from_ref(double dx_ref_index, double dy_ref_index, double theta_ref_index) {
    Eigen::Matrix3d Ar;
    
    Ar.setIdentity(n,n);
    Ar(0,2) = ( -dx_ref_index*sin(theta_ref_index) - dy_ref_index*cos(theta_ref_index) ) * T;
    Ar(1,2) = (  dx_ref_index*cos(theta_ref_index) - dy_ref_index*sin(theta_ref_index) ) * T;    

    return Ar;   
}

Eigen::Matrix3d B_from_ref(double theta_ref_index) {
    Eigen::Matrix3d Br;
    
    Br.setZero(n,m);
    Br(0,0) = cos(theta_ref_index)*T; Br(0,1) = -sin(theta_ref_index)*T;
    Br(1,0) = sin(theta_ref_index)*T; Br(1,1) =  cos(theta_ref_index)*T;
                                                                            Br(2,2) = T;    
    
    return Br;   
}





























/*
void ref_trajectory_mecanum_xy(void) {

    // ================================= reference trajectory ================================== //        
    //double* tb = LinearSpacedArray(0, 10, std::size_t(10/dt + 1) ).data();
    //Eigen::VectorXd tb = Eigen::VectorXd::Map(LinearSpacedArray(0, 10, std::size_t(10/dt + 1) ).data(), LinearSpacedArray(0, 10, std::size_t(10/dt + 1) ).size());    
    //Eigen::Matrix<double, int(10/dt + 1), 1> tb = Eigen::Matrix<double, int(10/dt + 1), 1>::Map(LinearSpacedArray(0, 10, std::size_t(10/dt + 1) ).data());    
    //Eigen::Array<double, int(10/dt + 1), 1> tb = Eigen::Array<double, int(10/dt + 1), 1>::Map(LinearSpacedArray(0, 10, std::size_t(10/dt + 1) ).data());    
    //Eigen::Array<double, int((10+3)/dt + 1), 1> t = Eigen::Array<double, int((10+3)/dt + 1), 1>::Map(LinearSpacedArray(0, (10+3), std::size_t((10+3)/dt + 1) ).data());    
    // --------------------------------------------------------------------------------------- //


    // ================================= time ================================== //
    t = Eigen::ArrayXd::LinSpaced(int((10+3)/dt + 1), 0, 13);        
    Eigen::ArrayXd tb = Eigen::ArrayXd::LinSpaced(int(10/dt + 1), 0, 10); 

    //static PROGMEM  float* ttt  = (int32_t*)t.data();                      
    // --------------------------------------------------------------------------------------- //

   
    // ================================= dx_ref ================================== //
    uint8_t num = 19;    
    uint8_t start_point = int((3*(1/T)-2*(num+1)/2)/2);
    Eigen::ArrayXd damping_profile = Eigen::ArrayXd::LinSpaced(int(0.5*num)+1, 0, 0.4737);        
    
    //Eigen::ArrayXd dx_ref(sizeof(tt)/sizeof(float));    
    dx_ref.setZero(sizeof(tt)/sizeof(float));       
    dx_ref.segment(start_point,              10        ) =  1.0*2*PI/10*(cos(2*PI*damping_profile)/2-0.5);    
    dx_ref.segment(start_point+10,           tb.size() ) = -1.0*2*PI/10*cos(2*PI*tb/10);    
    dx_ref.segment(start_point+10+tb.size(), 10        ) =  1.0*2*PI/10*(-cos(2*PI*damping_profile)/2-0.5);        
    // --------------------------------------------------------------------------------------- //

 
    // ================================= dy_ref ================================== //
    //Eigen::ArrayXd dy_ref(sizeof(tt)/sizeof(float));
    dy_ref.setZero(sizeof(tt)/sizeof(float));       
    dy_ref.segment(start_point,              10        ) =  1.0*4*PI/10*(cos(2*PI*damping_profile )/2-0.5);    
    dy_ref.segment(start_point+10,           tb.size() ) = -1.0*4*PI/10*cos(4*PI*tb/10);    
    dy_ref.segment(start_point+10+tb.size(), 10        ) =  1.0*4*PI/10*(-cos(2*PI*damping_profile )/2-0.5);    
    // --------------------------------------------------------------------------------------- //

    
    // ================================= theta_ref ================================== //
    //Eigen::ArrayXd theta_ref(sizeof(tt)/sizeof(float));
    theta_ref = dy_ref.binaryExpr(dx_ref, std::ptr_fun(::atan2));  //atan2 with array       
    theta_ref.segment(0, start_point+1) = atan2(dy_ref[start_point+1], dx_ref[start_point+1]) * Eigen::ArrayXd::Ones(start_point+1);    
    theta_ref.segment(sizeof(tt)/sizeof(float)-(start_point+1), start_point+1) = atan2(dy_ref[sizeof(tt)/sizeof(float)-(start_point+1)], dx_ref[sizeof(tt)/sizeof(float)-(start_point+1)]) * Eigen::ArrayXd::Ones(start_point+1);    
    
    //unwrap    
    for (int kk = 0; kk <sizeof(tt)/sizeof(float); ++kk)
        //unwrap(theta_ref[kk], theta_refa[kk]); //trash output
        if (theta_ref[kk] > -1.5) theta_ref[kk] = theta_ref[kk]-2*PI; //manually unwrap...
    // --------------------------------------------------------------------------------------- //

    
    // ================================= x_ref & y_ref ================================== //
    //Eigen::ArrayXd x_ref(sizeof(tt)/sizeof(float)), y_ref(sizeof(tt)/sizeof(float));
    x_ref.setZero(sizeof(tt)/sizeof(float)); y_ref.setZero(sizeof(tt)/sizeof(float));           
    for (int i = 0; i <sizeof(tt)/sizeof(float)-1; ++i) {
        x_ref[i+1] = x_ref[i] + dx_ref[i+1] * T;
        y_ref[i+1] = y_ref[i] + dy_ref[i+1] * T;
    }
    // --------------------------------------------------------------------------------------- //


    // ================================= ddx_ref & ddy_ref ================================== //
    //Eigen::ArrayXd ddx_ref(sizeof(tt)/sizeof(float)), ddy_ref(sizeof(tt)/sizeof(float));
    ddx_ref.setZero(sizeof(tt)/sizeof(float)); ddy_ref.setZero(sizeof(tt)/sizeof(float));           
    for (int i = 0; i <sizeof(tt)/sizeof(float)-1; ++i) {
        ddx_ref[i+1] = (dx_ref[i+1]-dx_ref[i])/T;
        ddy_ref[i+1] = (dy_ref[i+1]-dy_ref[i])/T;
    }
    // --------------------------------------------------------------------------------------- //


    // ================================= dtheta_ref & ddtheta_ref ================================== //
    //Eigen::ArrayXd dtheta_ref(sizeof(tt)/sizeof(float)), ddtheta_ref(sizeof(tt)/sizeof(float));
    dtheta_ref.setZero(sizeof(tt)/sizeof(float)); ddtheta_ref.setZero(sizeof(tt)/sizeof(float));           
    for (int i = 0; i <sizeof(tt)/sizeof(float)-1; ++i) {
        dtheta_ref[i+1] = (theta_ref[i+1]-theta_ref[i])/T;
        ddtheta_ref[i+1] = (dtheta_ref[i+1]-dtheta_ref[i])/T;
    }    
    // --------------------------------------------------------------------------------------- //


    // ================================= insert to global variables ================================== //
    xr << x_ref, y_ref, theta_ref;
    ur << dx_ref, dy_ref, dtheta_ref;
    ddxr << ddx_ref, ddy_ref, ddtheta_ref;
    // --------------------------------------------------------------------------------------- //


    // ================================= print ================================== //            
    for (int kk = 0 ; kk < sizeof(tt)/sizeof(float); ++kk) {
        Serial.print(xr(kk,0)); Serial.print(" ");
        Serial.print(xr(kk,1)); Serial.print(" ");
        Serial.println(xr(kk,2)); 
    }    

    for (int kk = 0 ; kk < sizeof(tt)/sizeof(float); ++kk) {
        Serial.print(x_ref(kk)); Serial.print(" ");
        Serial.print(y_ref(kk)); Serial.print(" ");
        Serial.print(theta_ref(kk)); Serial.print(" ");        

        Serial.print(dx_ref(kk)); Serial.print(" ");
        Serial.print(dy_ref(kk)); Serial.print(" ");
        Serial.print(dtheta_ref(kk)); Serial.print(" ");        

        Serial.print(ddx_ref(kk)); Serial.print(" ");
        Serial.print(ddy_ref(kk)); Serial.print(" ");
        Serial.println(ddtheta_ref(kk));
    }    
    // --------------------------------------------------------------------------------------- //    
}
*/





    
