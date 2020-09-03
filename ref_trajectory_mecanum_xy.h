// to use this function, include to qp_sim.h file
/*
void ref_trajectory_mecanum_xy(void);
inline double constrainAngle(double x);
inline double angleConv(double angle);
inline double angleDiff(double a,double b);
inline double unwrap(double previousAngle,double newAngle);
*/


// to use this function, include to main file (ino file)
/*
void ref_trajectory_mecanum_xy(void) {

    // ================================= reference trajectory =================================== //        
    //double* tb = LinearSpacedArray(0, 10, std::size_t(10/dt + 1) ).data();
    //Eigen::VectorXd tb = Eigen::VectorXd::Map(LinearSpacedArray(0, 10, std::size_t(10/dt + 1) ).data(), LinearSpacedArray(0, 10, std::size_t(10/dt + 1) ).size());    
    //Eigen::Matrix<double, int(10/dt + 1), 1> tb = Eigen::Matrix<double, int(10/dt + 1), 1>::Map(LinearSpacedArray(0, 10, std::size_t(10/dt + 1) ).data());    
    //Eigen::Array<double, int(10/dt + 1), 1> tb = Eigen::Array<double, int(10/dt + 1), 1>::Map(LinearSpacedArray(0, 10, std::size_t(10/dt + 1) ).data());    
    //Eigen::Array<double, int((10+3)/dt + 1), 1> t = Eigen::Array<double, int((10+3)/dt + 1), 1>::Map(LinearSpacedArray(0, (10+3), std::size_t((10+3)/dt + 1) ).data());    
    // ------------------------------------------------------------------------------------------ //


    // ================================= time ================================================== //
    t = Eigen::ArrayXd::LinSpaced(int((10+3)/dt + 1), 0, 13);        
    Eigen::ArrayXd tb = Eigen::ArrayXd::LinSpaced(int(10/dt + 1), 0, 10); 

    //static PROGMEM  float* ttt  = (int32_t*)t.data();                      
    // ---------------------------------------------------------------------------------------- //

   
    // ================================= dx_ref =============================================== //
    uint8_t num = 19;    
    uint8_t start_point = int((3*(1/T)-2*(num+1)/2)/2);
    Eigen::ArrayXd damping_profile = Eigen::ArrayXd::LinSpaced(int(0.5*num)+1, 0, 0.4737);        
    
    //Eigen::ArrayXd dx_ref(sizeof(tt)/sizeof(float));    
    dx_ref.setZero(sizeof(tt)/sizeof(float));       
    dx_ref.segment(start_point,              10        ) =  1.0*2*PI/10*(cos(2*PI*damping_profile)/2-0.5);    
    dx_ref.segment(start_point+10,           tb.size() ) = -1.0*2*PI/10*cos(2*PI*tb/10);    
    dx_ref.segment(start_point+10+tb.size(), 10        ) =  1.0*2*PI/10*(-cos(2*PI*damping_profile)/2-0.5);        
    // ---------------------------------------------------------------------------------------- //

 
    // ================================= dy_ref =============================================== //
    //Eigen::ArrayXd dy_ref(sizeof(tt)/sizeof(float));
    dy_ref.setZero(sizeof(tt)/sizeof(float));       
    dy_ref.segment(start_point,              10        ) =  1.0*4*PI/10*(cos(2*PI*damping_profile )/2-0.5);    
    dy_ref.segment(start_point+10,           tb.size() ) = -1.0*4*PI/10*cos(4*PI*tb/10);    
    dy_ref.segment(start_point+10+tb.size(), 10        ) =  1.0*4*PI/10*(-cos(2*PI*damping_profile )/2-0.5);    
    // ---------------------------------------------------------------------------------------- //

    
    // ================================= theta_ref ============================================ //
    //Eigen::ArrayXd theta_ref(sizeof(tt)/sizeof(float));
    theta_ref = dy_ref.binaryExpr(dx_ref, std::ptr_fun(::atan2));  //atan2 with array       
    theta_ref.segment(0, start_point+1) = atan2(dy_ref[start_point+1], dx_ref[start_point+1]) * Eigen::ArrayXd::Ones(start_point+1);    
    theta_ref.segment(sizeof(tt)/sizeof(float)-(start_point+1), start_point+1) = atan2(dy_ref[sizeof(tt)/sizeof(float)-(start_point+1)], dx_ref[sizeof(tt)/sizeof(float)-(start_point+1)]) * Eigen::ArrayXd::Ones(start_point+1);    
    
    //unwrap    
    for (int kk = 0; kk <sizeof(tt)/sizeof(float); ++kk)
        //unwrap(theta_ref[kk], theta_refa[kk]); //trash output
        if (theta_ref[kk] > -1.5) theta_ref[kk] = theta_ref[kk]-2*PI; //manually unwrap...
    // --------------------------------------------------------------------------------------- //

    
    // ================================= x_ref & y_ref ======================================= //
    //Eigen::ArrayXd x_ref(sizeof(tt)/sizeof(float)), y_ref(sizeof(tt)/sizeof(float));
    x_ref.setZero(sizeof(tt)/sizeof(float)); y_ref.setZero(sizeof(tt)/sizeof(float));           
    for (int i = 0; i <sizeof(tt)/sizeof(float)-1; ++i) {
        x_ref[i+1] = x_ref[i] + dx_ref[i+1] * T;
        y_ref[i+1] = y_ref[i] + dy_ref[i+1] * T;
    }
    // -------------------------------------------------------------------------------------- //


    // ================================= ddx_ref & ddy_ref ================================== //
    //Eigen::ArrayXd ddx_ref(sizeof(tt)/sizeof(float)), ddy_ref(sizeof(tt)/sizeof(float));
    ddx_ref.setZero(sizeof(tt)/sizeof(float)); ddy_ref.setZero(sizeof(tt)/sizeof(float));           
    for (int i = 0; i <sizeof(tt)/sizeof(float)-1; ++i) {
        ddx_ref[i+1] = (dx_ref[i+1]-dx_ref[i])/T;
        ddy_ref[i+1] = (dy_ref[i+1]-dy_ref[i])/T;
    }
    // --------------------------------------------------------------------------------------- //


    // ================================= dtheta_ref & ddtheta_ref ============================ //
    //Eigen::ArrayXd dtheta_ref(sizeof(tt)/sizeof(float)), ddtheta_ref(sizeof(tt)/sizeof(float));
    dtheta_ref.setZero(sizeof(tt)/sizeof(float)); ddtheta_ref.setZero(sizeof(tt)/sizeof(float));           
    for (int i = 0; i <sizeof(tt)/sizeof(float)-1; ++i) {
        dtheta_ref[i+1] = (theta_ref[i+1]-theta_ref[i])/T;
        ddtheta_ref[i+1] = (dtheta_ref[i+1]-dtheta_ref[i])/T;
    }    
    // --------------------------------------------------------------------------------------- //

    // ================================= insert to global variables ========================== //
    xr << x_ref, y_ref, theta_ref;
    ur << dx_ref, dy_ref, dtheta_ref;
    ddxr << ddx_ref, ddy_ref, ddtheta_ref;
    // --------------------------------------------------------------------------------------- //
}

// ================================= unwrap ================================================== //        
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
// ----------------------------------------------------------------------------------------- //
*/





    
