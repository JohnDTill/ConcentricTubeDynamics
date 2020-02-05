#include <iostream>
#include <fstream>
#include <convexoptimization.h>
#include <commonmath.h>
#include <numericalintegration.h>
#include <timemanagement.h>
using namespace ContinuumRobotLibrary;

#ifdef QT_CORE_LIB
#include <simpleplotting.h>
#endif

//Independent Parameters
const Array2d Bxy = {0, 0}; //Transverse material damping coefficients
const Array2d Bz = {0, 0};  //Axial material damping coefficients
const double alpha = 0; //BDF-alpha implicit time integration parameter
const double dt = 1e-5; //Time step
const double dt_static = 1e-3; //Time step to solve static motion prior to snap
const double snap_simulation_duration = 12e-3; //Time for dynamic simulation
const ArrayXi N =
    (ArrayXi(6)<<15,30,10,15,10,4).finished();//Spatial resolution of PDE sections
const Array2d ri = {0, 0.00062}; //Inner radii
const Array2d ro = {0.000508, 0.001055}; //Outer radii
const Array2d E = {81.97e9, 210e9}; //Youngs moduli
const Array2d rho = {6493, 8000}; //Material densities
const double beta0 = -0.112; //Inner tube transmission length (base to actuator)
const double l1 = 0.0531; //Length of outer tube
const double Ltip = 0.002; //Extension of inner tube beyond outer tube termination
const Vector3d g = -9.81*Vector3d::UnitY(); //Gravitational acceleration
const Array2d ustarx = {37.12, 8.72}; //Precurvature components in x-direction
const Array2d ustar_start =
    {l1+Ltip-0.04713, l1-0.0399}; //Precurvature start locations
const Array2d ustar_end =
    {l1+Ltip-0.0115, l1}; //Precurvature end locations
const Array2d mM = (Array2d(0.0278, 0.0714) + 0.02)/1e3; //Tracking marker masses
const Array2d roM = {0.00273, 0.00446}; //Radii of tracking markers
const Array2d tM = {0.00113, 0.00109}; //Thickness of tracking markers
double angularPosition(double t)
    { return -10*t; } //Input profile (negligible to snap motion)
double angularVelocity(double /*t*/)
    { return -10; } //Derivative of input motion profile
const double eps = 3e2; //Smoothing for dynamic friction
const double mu_prime = 0.008; //Dynamic friction coefficient

//Other parameters
const ArrayXi M = (ArrayXi(6) << 3,22,22,22,22,19).finished();
    //Number of state variables in each section
const ArrayXi P = (ArrayXi(6) << 2,11,11,11,11,9).finished();
    //Number of time derivatives in each section
const Array2d G = E/(2*(1 + 0.3)); //Shear moduli
const Array2d A = pi*(ro.square() - ri.square()); //Cross-sectional areas
const Array2d I = pi*(ro.pow(4) - ri.pow(4))/4; //Area moment of inertias
const Array2d EI = E*I; //Bending stiffnesses
const Array2d GJ = G*2*I; //Twisting stiffnesses
const Array2d rhoI = rho*I; //Mass moment of inertia
const Array2d rhoIzz = rho*2*I; //Polar mass moment of inertia
const DiagonalMatrix<double, 3> K0 (EI[0],EI[0],GJ[0]); //Tube 0 stiffness
const DiagonalMatrix<double, 3> K1 (EI[1],EI[1],GJ[1]); //Tube 1 stiffness
const DiagonalMatrix<double, 3> B0 (Bxy[0],Bxy[0],Bz[0]); //Tube 0 mat. damping
const DiagonalMatrix<double, 3> B1 (Bxy[1],Bxy[1],Bz[1]); //Tube 1 mat. damping
const DiagonalMatrix<double, 3> K0_c0B0_inv =
    DiagonalMatrix<double,3>(EI[0]+Bxy[0],EI[0]+Bxy[0],GJ[0]+Bz[0]).inverse();
const double mass_per_length = (rho*A).sum();
const double rhoA0 = rho[0]*A[0];
const DiagonalMatrix<double, 3> rhoJ0 =
    rho[0]*DiagonalMatrix<double, 3>(I[0],I[0],2*I[0]);
const double c0 = TimeManagerBdfAlpha::getC0(dt,alpha);
const Array2d IT = mM*(3*(roM.square() + ro.square()) + tM.square())/12;
const Array2d IA = mM*(roM.square() + ro.square())/2;
const DiagonalMatrix<double,3> JT (IT[0], IT[0], IA[0]);
const ArrayXd s0 = (ArrayXd(6) << beta0,
                                      0,
                         ustar_start[0],
                         ustar_start[1],
                           ustar_end[0],
                                     l1 ).finished();
const ArrayXd ds = (ArrayXd(6) << (s0(1)-s0(0))/(N(0)-1),
                                  (s0(2)-s0(1))/(N(1)-1),
                                  (s0(3)-s0(2))/(N(2)-1),
                                  (s0(4)-s0(3))/(N(3)-1),
                                  (s0(5)-s0(4))/(N(4)-1),
                                  Ltip/(N(5)-1) ).finished();
static double t = 0;

//Approximation to sign() function
double smoothSaturation(double val){
    return val/sqrt(val*val + eps*eps);
}

//S-Curve calculations
static void twoTubeOde(VectorXd& y_s_out, double s, VectorXd& y){
    //Unpack state vector
    Matrix3d R = Map<Matrix3d>(&y[3]);
    Vector3d n = Map<Vector3d>(&y[12]);
    Vector2d mb_xy = Map<Vector2d>(&y[15]);
    Vector2d mb_z = Map<Vector2d>(&y[17]);
    double theta = y[19];

    //Pack state vector derivative
    Map<Vector3d> p_s(&y_s_out[0]);
    Map<Matrix3d> R_s(&y_s_out[3]);
    Map<Vector3d> n_s(&y_s_out[12]);
    Map<Vector2d> mb_xys(&y_s_out[15]);
    Map<Vector2d> mb_zs(&y_s_out[17], 2);
    Map<VectorXd> theta_s(&y_s_out[19], 1);

    bool precurved_0 = (s > ustar_start[0] && s < ustar_end[0]);
    bool precurved_1 = (s > ustar_start[1]);
    Vector3d u_star0 = Vector3d(ustarx[0]*precurved_0, 0, 0);
    Vector3d u_star1 = Rz(theta)*Vector3d(ustarx[1]*precurved_1, 0, 0);

    Vector2d u_xy = (mb_xy+EI[0]*u_star0.head(2)+EI[1]*u_star1.head(2))/EI.sum();
    double u_z0 = mb_z(0)/GJ[0] + u_star0.z();
    double u_z1 = mb_z(1)/GJ[1] + u_star1.z();

    Vector3d mb;  mb << mb_xy, mb_z(0);
    Vector3d u;  u << u_xy, u_z0;

    //ODEs
    p_s = R.col(2);
    R_s = hat_postmultiply(R,u);
    n_s = Vector3d::Zero();
    mb_xys=(-u.cross(mb)-Vector3d::UnitZ().cross(transposeMultiply(R,n))).head(2);
    mb_zs(0) = EI[0]*(u.x()*u_star0.y() - u.y()*u_star0.x());
    mb_zs(1) = EI[1]*(u.x()*u_star1.y() - u.y()*u_star1.x());
    theta_s << u_z1 - u_z0;
}

template<int half_num_samples = 50, int integration_pts = 100>
inline Matrix2Xd getSCurve(){
    RowVectorXd theta_f = RowVectorXd::LinSpaced(2*half_num_samples,0,2*pi);
    RowVectorXd theta_b(2*half_num_samples);

    for(int i = 0; i < half_num_samples; i++){
        VectorXd yL(20);
        yL << 0,0,0, 1,0,0,0,1,0,0,0,1, VectorXd::Zero(7), theta_f(i);

        MatrixXd Y = ode4<twoTubeOde,integration_pts>(yL, l1, 0);

        double theta0 = Y(19,integration_pts-1);
        double mb_z0 = Y(17,integration_pts-1);

        theta_b(i) = theta0 - mb_z0/GJ[0]*beta0;

        //Set the upper half by symmetry
        theta_b(2*half_num_samples-1-i) = 2*pi - theta_b(i);
    }

    Matrix2Xd data(2,2*half_num_samples);
    data << theta_f, theta_b;

    return data;
}

Array3d findSnappingValues(Matrix2Xd SCurve){
    double theta_b_snap, theta_f_presnap, theta_f_postsnap;

    //Search for an inflection point in the base angle
    int snapping_index = 0;
    theta_b_snap = 0;
    for(int i = 0; i < SCurve.cols()/2; i++){
        if(SCurve(1,i) > theta_b_snap){
            theta_b_snap = SCurve(1,i);
            snapping_index = i;
        }
    }
    if(snapping_index == SCurve.cols()/2-1){
        std::cout << "Robot does not snap with given parameters." << std::endl;
        throw(1);
    }
    theta_f_presnap = SCurve(0,snapping_index);

    //Search for the corresponding point above (theta_b_snap, theta_f_presnap)
    if(theta_b_snap < 2*pi){
        for(auto i = SCurve.cols()-1; true; i--){
            if(SCurve(1,i) < theta_b_snap){
                theta_f_postsnap = SCurve(0,i);
                break;
            }
        }
    }else{
        for(int i = 0; true; i++){
            if(SCurve(1,i) > theta_b_snap - 2*pi){
                theta_f_postsnap = SCurve(0,i) + 2*pi;
            }
        }
    }

    return Array3d(theta_f_presnap, theta_f_postsnap, theta_b_snap);
}

//Internal residual calculations from ODEs
template<bool dynamic>
void lowerOde(Map<VectorXd> Ri,
              Map<VectorXd> z_out,
              VectorXd y,
              VectorXd y_s_FD,
              VectorXd z_h){
    double alpha_t = y(1);
    double mb_z = y(2);

    double u_zh = z_h(0);
    double alpha_th = z_h(1);

    double u_z = dynamic ? (mb_z - Bz[0]*u_zh)/(GJ[0] + c0*Bz[0]) : mb_z/GJ[0];

    double u_zt = dynamic ? c0*u_z + u_zh : 0;
    double alpha_tt = dynamic ? c0*alpha_t + alpha_th : 0;

    double alpha_s = u_z;
    double alpha_st = u_zt;
    double mb_zs = rhoIzz[0]*alpha_tt;

    Vector3d y_s = Vector3d(alpha_s, alpha_st, mb_zs);
    Ri = y_s - y_s_FD;

    z_out = Vector2d(u_z, alpha_t);
}

template<bool dynamic, bool precurved_0, bool precurved_1>
void upperOde(Map<VectorXd> Ri,
              Map<VectorXd> z_out,
              VectorXd y,
              VectorXd y_s_FD,
              VectorXd z_h){
    Vector4d h = y.segment<4>(3);
    Vector3d q = y.segment<3>(7);
    Vector3d w = y.segment<3>(10);
    Vector3d n = y.segment<3>(13);
    Vector2d mb_xy = y.segment<2>(16);
    double mb_z0 = y(18);
    double mb_z1 = y(19);
    double theta = y(20);
    double gamma = y(21);

    Matrix3d R = quat2rot(h);

    Vector3d ustar0 = Vector3d(ustarx[0]*precurved_0, 0, 0);
    Vector3d ustar1 = Rz(theta)*Vector3d(ustarx[1]*precurved_1, 0, 0);

    Vector3d q_h = z_h.segment<3>(0);
    Vector3d w_h = z_h.segment<3>(3);
    Vector2d u_xyh = z_h.segment<2>(6);
    double u_z0h = z_h(8);
    double u_z1h = z_h(9);
    double gamma_h = z_h(10);

    Vector2d u_xy;
    if(dynamic){
        u_xy = (mb_xy+EI[0]*ustar0.head(2)+EI[1]*ustar1.head(2)-Bxy.sum()*u_xyh)
               / (EI.sum() + c0*Bxy.sum());
    }else{
        u_xy = (mb_xy + EI[0]*ustar0.head(2) + EI[1]*ustar1.head(2)) / EI.sum();
    }
    double u_z0 = dynamic ? (mb_z0+GJ[0]*ustar0.z()-Bz[0]*u_z0h)/(GJ[0]+c0*Bz[0])
                          : mb_z0/GJ[0] + ustar0.z();
    double u_z1 = dynamic ? (mb_z1+GJ[1]*ustar1.z()-Bz[1]*u_z1h)/(GJ[1]+c0*Bz[1])
                          : mb_z1/GJ[1] + ustar1.z();

    Vector3d q_t = dynamic ? (c0*q + q_h).eval() : Vector3d::Zero();
    Vector3d w_t = dynamic ? (c0*w + w_h).eval() : Vector3d::Zero();
    Vector2d u_xyt = dynamic ? (c0*u_xy + u_xyh).eval() : Vector2d::Zero();
    double u_z0t = dynamic ? c0*u_z0 + u_z0h : 0;
    double u_z1t = dynamic ? c0*u_z1 + u_z1h : 0;
    double gamma_t = dynamic ? c0*gamma + gamma_h : 0;

    Vector3d mb;  mb << mb_xy, mb_z0+mb_z1;
    Vector3d u;  u << u_xy, u_z0;
    Vector3d u_t;  u_t << u_xyt, u_z0t;
    Vector3d mb0 = K0*(u - ustar0) + B0*u_t;
    Vector3d mb1 = K1*(u - ustar1) + B1*u_t;

    Vector3d ps = R.col(2);
    Vector4d hs = quatDerivLocal(h,u);
    Vector3d qs = -u.cross(q) + w.cross(Vector3d::UnitZ());
    Vector3d ws = u_t - u.cross(w);
    Vector3d ns = mass_per_length*( R*(w.cross(q) + q_t) - g );
    Vector2d mb_xys =
            (-u.cross(mb)-Vector3d::UnitZ().cross(transposeMultiply(R,n))
            + (rhoI).sum()*w_t
            + Vector3d(w.y(),-w.x(),0) * (rhoI[0]*w.z() + rhoI[1]*(w.z()+gamma))
            ).head(2);
    double mb_z0s = (mb0.x()*u.y() - mb0.y()*u.x()) + rhoIzz[0]*w_t.z();
    double mb_z1s = (mb1.x()*u.y() - mb1.y()*u.x()) + rhoIzz[1]*(w_t.z()+gamma_t);
    double theta_s = u_z1 - u_z0;
    double gamma_s = u_z1t - u_z0t;

    VectorXd y_s(22);
    y_s << ps, hs, qs, ws, ns, mb_xys, mb_z0s, mb_z1s, theta_s, gamma_s;
    Ri = y_s - y_s_FD;
    z_out << q, w, u_xy, u_z0, u_z1, gamma;
}

template<bool dynamic>
void tipOde(Map<VectorXd> Ri,
            Map<VectorXd> z_out,
            VectorXd y,
            VectorXd y_s_FD,
            VectorXd z_h){
    Vector4d h = y.segment<4>(3);
    Vector3d q = y.segment<3>(7);
    Vector3d w = y.segment<3>(10);
    Vector3d n = y.segment<3>(13);
    Vector3d mb = y.segment<3>(16);

    Matrix3d R = quat2rot(h);

    Vector3d u = K0_c0B0_inv*mb;

    Vector3d q_h = z_h.segment<3>(0);
    Vector3d w_h = z_h.segment<3>(3);
    Vector3d u_h = z_h.segment<3>(6);

    Vector3d q_t = dynamic ? (c0*q + q_h).eval() : Vector3d::Zero();
    Vector3d w_t = dynamic ? (c0*w + w_h).eval() : Vector3d::Zero();
    Vector3d u_t = dynamic ? (c0*u + u_h).eval() : Vector3d::Zero();

    Vector3d ps = R.col(2);
    Vector4d hs = quatDerivLocal(h,u);
    Vector3d qs = -u.cross(q) + w.cross(Vector3d::UnitZ());
    Vector3d ws = u_t - u.cross(w);
    Vector3d ns = rhoA0*( R*(w.cross(q) + q_t) - g );
    Vector3d mbs = -u.cross(mb) - Vector3d::UnitZ().cross(transposeMultiply(R,n))
                   + rhoJ0*w_t + w.cross(rhoJ0*w);

    VectorXd y_s(19);
    y_s << ps, hs, qs, ws, ns, mbs;
    Ri = y_s - y_s_FD;
    z_out << q, w, u;
}

//Objective function, grid data, and sub-functions
static MatrixXd Y[6], Z[6], Z_h[6];
void initializeGrids(){
    for(int i = 0; i < 6; i++){
        Y[i].resize(M[i],N[i]);
        Z[i].resize(P[i],N[i]-1);
        Z_h[i].resize(P[i],N[i]-1);
    }
}
static MatrixXd TIPSTATE = MatrixXd(17,1);
static MatrixXd TIPSTATE_h = MatrixXd(17,1);
static MatrixXd MIDSTATE = MatrixXd(12,1);
static MatrixXd MIDSTATE_h = MatrixXd(12,1);

void map(int grid,
         int row,
         int col,
         int num_rows,
         int num_cols,
         VectorXd& guess,
         int& index){

    Y[grid].block(row,col,num_rows,num_cols) =
            Map<MatrixXd>(&guess[index],num_rows,num_cols);

    index += num_rows*num_cols;
}

void mapGuessToGrids(VectorXd& guess){
    int index = 0;
    map(0,2,0,1,1,guess,index);
    map(0,0,1,3,N[0]-1,guess,index);
    map(1,13,0,5,1,guess,index);
    map(1,19,0,1,1,guess,index);
    map(1,0,1,22,N[1]-1,guess,index);
    map(2,0,1,22,N[2]-1,guess,index);
    map(3,0,1,22,N[3]-1,guess,index);
    map(4,0,1,22,N[4]-2,guess,index);
    map(4,0,N[4]-1,19,1,guess,index);
    map(4,20,N[4]-1,2,1,guess,index);
    map(5,0,1,19,N[5]-2,guess,index);
    map(5,0,N[5]-1,13,1,guess,index);
}

template<bool dynamic, int i>
void handleDiscontinuousCurvature(){
    double theta = Y[i+1](20,0);
    double friction_ratio = dynamic ? smoothSaturation(Y[i+1](21,0)) : 0;

    Vector2d mb_xy = Y[i+1].block<2,1>(16,0);

    double s_start = s0[i+1] - ds[i];
    double s_end = s0[i+1] + ds[i+1];
    Vector3d ustar0_before = Vector3d(ustarx[0]*(s_start >= ustar_start[0] && s_start <= ustar_end[0]), 0, 0);
    Vector3d ustar1_before = Rz(theta)*Vector3d(ustarx[1]*(s_start >= ustar_start[1] && s_start <= ustar_end[1]), 0, 0);
    Vector3d ustar0_after = Vector3d(ustarx[0]*(s_end >= ustar_start[0] && s_end <= ustar_end[0]), 0, 0);
    Vector3d ustar1_after = Rz(theta)*Vector3d(ustarx[1]*(s_end >= ustar_start[1] && s_end <= ustar_end[1]), 0, 0);

    Vector2d u_xyh_before = MIDSTATE_h.block<2,1>(4*(i-1),0);
    Vector2d u_xyh_after = MIDSTATE_h.block<2,1>(4*(i-1)+2,0);
    Vector2d u_xy_before = dynamic ? ((mb_xy + EI[0]*ustar0_before.head<2>() + EI[1]*ustar1_before.head<2>() - Bxy.sum()*u_xyh_before) / (EI.sum() + c0*Bxy.sum())).eval()
                             : (mb_xy + EI[0]*ustar0_before.head<2>() + EI[1]*ustar1_before.head<2>()) / (EI.sum());
    Vector2d u_xy_after = dynamic ? ((mb_xy + EI[0]*ustar0_after.head<2>() + EI[1]*ustar1_after.head<2>() - Bxy.sum()*u_xyh_after) / (EI.sum() + c0*Bxy.sum())).eval()
                             : (mb_xy + EI[0]*ustar0_after.head<2>() + EI[1]*ustar1_after.head<2>()) / (EI.sum());
    MIDSTATE.block<2,1>(4*(i-1),0) = u_xy_before;
    MIDSTATE.block<2,1>(4*(i-1)+2,0) = u_xy_after;

    Vector2d u_xyt_before = dynamic ? (c0*u_xy_before + u_xyh_before).eval() : Vector2d::Zero();
    Vector2d u_xyt_after = dynamic ? (c0*u_xy_before + u_xyh_before).eval() : Vector2d::Zero();

    Vector2d mb0_xy_before = EI[0]*(u_xy_before - ustar0_before.head(2)) + Bxy[0]*u_xyt_before;
    Vector2d mb0_xy_after = EI[0]*(u_xy_after - ustar0_after.head(2)) + Bxy[0]*u_xyt_after;

    Vector2d M_transferred = mb0_xy_before - mb0_xy_after;

    double tau_f = mu_prime * M_transferred.norm() * friction_ratio;

    Y[i+1](18,0) = Y[i](18,N[i]-1) - tau_f;
    Y[i+1](19,0) = Y[i](19,N[i]-1) + tau_f;
}

template<bool dynamic>
void applyStrongBCs(){
    //Actuator
    Y[0](0,0) = angularPosition(t);
    Y[0](1,0) = dynamic ? angularVelocity(t) : 0;

    //Transition to baseplate
    double alpha0 = Y[0](0,N[0]-1);
    double alpha0_t = Y[0](1,N[0]-1);
    double mb0_z0 = Y[0](2,N[0]-1);

    Vector2d mb_xy0 = Y[1].block<2,1>(16,0);
    Vector3d ustar0_0 = Vector3d(ustarx[0]*(0 >= ustar_start[0] && 0 <= ustar_end[0]), 0, 0);
    Vector3d ustar1_0 = Rz(-alpha0)*Vector3d(ustarx[1]*(0 >= ustar_start[1] && 0 <= ustar_end[1]), 0, 0);
    Vector2d u_xyh0 = TIPSTATE_h.block<2,1>(15,0);
    Vector2d u_xy0 = dynamic ? ((mb_xy0 + EI[0]*ustar0_0.head<2>() + EI[1]*ustar1_0.head<2>() - Bxy.sum()*u_xyh0) / (EI.sum() + c0*Bxy.sum())).eval()
                             : (mb_xy0 + EI[0]*ustar0_0.head<2>() + EI[1]*ustar1_0.head<2>()) / (EI.sum());
    TIPSTATE.block<2,1>(15,0) = u_xy0;
    Vector2d u_xy0t = dynamic ? (c0*u_xy0 + u_xyh0).eval() : Vector2d::Zero();
    Vector2d mb0xy0 = EI[0]*(u_xy0 - ustar0_0.head(2)) + Bxy[0]*u_xy0t;

    double gamma0 = -alpha0_t;
    double friction_ratio_0 = dynamic ? smoothSaturation(gamma0) : 0;
    double tau_f_0 = mu_prime * mb0xy0.norm() * friction_ratio_0;

    Vector3d p0 = Vector3d::Zero();
    Vector4d h0 = Vector4d(cos(alpha0/2),0,0,sin(alpha0/2));
    Vector3d q0 = Vector3d::Zero();
    Vector3d w0(0,0,alpha0_t);

    Y[1].topLeftCorner(13,1) << p0, h0, q0, w0;
    Y[1](18,0) = mb0_z0 - tau_f_0;
    Y[1](20,0) = -alpha0;
    Y[1](21,0) = -alpha0_t;

    //Points with discontinuous precurvature
    Y[2].col(0) = Y[1].col(N[1]-1);
    Y[3].col(0) = Y[2].col(N[2]-1);
    Y[4].col(0) = Y[3].col(N[3]-1);

    handleDiscontinuousCurvature<dynamic,1>();
    handleDiscontinuousCurvature<dynamic,2>();
    handleDiscontinuousCurvature<dynamic,3>();

    //Outer tube termination - accounts for tip inertias
    Matrix3d RL = quat2rot(Y[4].block<4,1>(3,N[4]-1));
    Vector3d qL = Y[4].block<3,1>(7,N[4]-1);
    Vector3d wL = Y[4].block<3,1>(10,N[4]-1);
    double theta_L = Y[4](20, N[4]-1);
    double gammaL = Y[4](21,N[4]-1);
    double wL1z = gammaL + wL.z();
    TIPSTATE.block<3,1>(0,0) = qL;
    TIPSTATE.block<3,1>(3,0) = wL;
    TIPSTATE(6,0) = wL1z;

    Vector3d qL_t = dynamic ? (c0*qL + TIPSTATE_h.block<3,1>(0,0)).eval() : Vector3d::Zero();
    Vector3d wL_t = dynamic ? (c0*wL + TIPSTATE_h.block<3,1>(3,0)).eval() : Vector3d::Zero();
    double wL1z_t = dynamic ? c0*wL1z + TIPSTATE_h(6,0) : 0;
    Vector3d aL = RL*(qL_t + wL.cross(qL));
    Vector3d FL = mM[1]*(aL - g);

    Vector2d mb_xyL = Y[4].block<2,1>(16,N[4]-1);
    Vector3d ustar0 = Vector3d(ustarx[0]*(l1 >= ustar_start[0] && l1 <= ustar_end[0]), 0, 0);
    Vector3d ustar1 = Rz(theta_L)*Vector3d(ustarx[1]*(l1 >= ustar_start[1] && l1 <= ustar_end[1]), 0, 0);
    Vector2d u_xyhL = TIPSTATE_h.block<2,1>(13,0);
    Vector2d u_xyL = dynamic ? ((mb_xyL + EI[0]*ustar0.head<2>() + EI[1]*ustar1.head<2>() - Bxy.sum()*u_xyhL) / (EI.sum() + c0*Bxy.sum())).eval()
                             : (mb_xyL + EI[0]*ustar0.head<2>() + EI[1]*ustar1.head<2>()) / (EI.sum());
    TIPSTATE.block<2,1>(13,0) = u_xyL;
    Vector2d u_xyLt = dynamic ? (c0*u_xyL + u_xyhL).eval() : Vector2d::Zero();
    Vector2d mb1xyL = EI[1]*(u_xyL - ustar1.head(2)) + Bxy[0]*u_xyLt;

    Vector2d MbxyL = IT[1]*wL_t.head(2) + IT[1]*wL1z*Vector2d(wL.y(), -wL.x()); //Tip internal moment

    Vector2d M_transferred = mb1xyL + MbxyL;

    double friction_ratio = dynamic ? smoothSaturation(gammaL) : 0;
    double tau_f = mu_prime * M_transferred.norm() * friction_ratio;

    Y[4](19,N[4]-1) = -IA[1] * wL1z_t - tau_f; //Torque on tube 1
    Y[5].block<13,1>(0,0) = Y[4].block<13,1>(0,N[4]-1);
    Y[5].block<3,1>(13,0) = Y[4].block<3,1>(13,N[4]-1) + FL;
    Y[5].block<2,1>(16,0) = mb_xyL + MbxyL;
    Y[5](18,0) = Y[4](18,N[4]-1) - tau_f;

    //Acount for inertia of second marker
    Matrix3d RT = quat2rot(Y[5].block<4,1>(3,N[5]-1));
    Vector3d qT = Y[5].block<3,1>(7,N[5]-1);
    Vector3d wT = Y[5].block<3,1>(10,N[5]-1);
    TIPSTATE.block<3,1>(7,0) = qT;
    TIPSTATE.block<3,1>(10,0) = wT;

    Vector3d qT_t = dynamic ? (c0*qT + TIPSTATE_h.block<3,1>(7,0)).eval() : Vector3d::Zero();
    Vector3d wT_t = dynamic ? (c0*wT + TIPSTATE_h.block<3,1>(10,0)).eval() : Vector3d::Zero();
    Vector3d aT = RT*(qT_t + wT.cross(qT));
    Vector3d FT = mM[0]*(aT - g);

    Y[5].block<3,1>(13,N[5]-1) = -FT;
    Y[5].block<3,1>(16,N[5]-1) = -(JT*wT_t) - wT.cross(JT*wT);
}

template<bool dynamic>
VectorXd objFunc(VectorXd guess){
    mapGuessToGrids(guess);
    applyStrongBCs<dynamic>();

    //Function pointer definition for internal residual functions
    typedef void(*InternalErrorFunction)(Map<VectorXd>,
                                         Map<VectorXd>,
                                         VectorXd,
                                         VectorXd,
                                         VectorXd);

    //List of internal residual functions for each section
    static InternalErrorFunction func[6] = {lowerOde<dynamic>,
                                            upperOde<dynamic,false,false>,
                                            upperOde<dynamic,true,false>,
                                            upperOde<dynamic,true,true>,
                                            upperOde<dynamic,false,true>,
                                            tipOde<dynamic>};

    //Loop over sections to calculate internal residual elements
    VectorXd res( (M*(N-1)).sum() );
    for(int i = 0; i < M.size(); i++){
        for(int j = 0; j < N[i]-1; j++){
            Map<VectorXd> z(&Z[i](0,j), P[i]);
            Map<VectorXd> Ri(&res((M*(N-1)).head(i).sum() + M[i]*j), M[i]);
            VectorXd ybar = (Y[i].col(j) + Y[i].col(j+1))/2;
            VectorXd y_s_FD = (Y[i].col(j+1) - Y[i].col(j)) / ds[i];
            func[i](Ri, z, ybar, y_s_FD, Z_h[i].col(j));
        }
    }

    return res;
}

int main(int, char**){
    Matrix2Xd s_curve = getSCurve();
    Array3d snapping_values = findSnappingValues(s_curve);
    //plot(s_curve.row(1), s_curve.row(0), "S-Curve",
    //     "θ<sub>b</sub> (rad)", "θ<sub>f</sub> (rad)");
	
	initializeGrids();

    //Approach the snapping point statically
    double theta_b_snap = snapping_values(2);
    VectorXd guess = VectorXd::Constant( (M*(N-1)).sum(), 1e-6 );
    guess = solveSparseLevenbergMarquardt<objFunc<false> >(guess, 1e-12, 1500, 1e-7);

    const double early_stop = 0.01;
    while(angularPosition(t) > -theta_b_snap + early_stop){
        t += dt_static;
        guess = solveSparseLevenbergMarquardt<objFunc<false> >(guess,1e-12,1500,1e-7);
        std::cout << "STATIC: t=" << t << ", theta_b=" << angularPosition(t)
                  << " of " << -theta_b_snap << std::endl;
    }

    //Start the dynamic simulation
    const int num_grids = 8;
    TimeManagerBdfAlpha grid[num_grids] = {TimeManagerBdfAlpha(Z[0], Z_h[0], dt, alpha),
                                   TimeManagerBdfAlpha(Z[1], Z_h[1], dt, alpha),
                                   TimeManagerBdfAlpha(Z[2], Z_h[2], dt, alpha),
                                   TimeManagerBdfAlpha(Z[3], Z_h[3], dt, alpha),
                                   TimeManagerBdfAlpha(Z[4], Z_h[4], dt, alpha),
                                   TimeManagerBdfAlpha(Z[5], Z_h[5], dt, alpha),
                                   TimeManagerBdfAlpha(TIPSTATE, TIPSTATE_h, dt, alpha),
                                   TimeManagerBdfAlpha(MIDSTATE, MIDSTATE_h, dt, alpha)};

    int M = static_cast<int>(snap_simulation_duration/dt);
    VectorXd theta_f(M), theta_b(M);
    MatrixXd p1(3*M,N.segment(1,4).sum()), p0(3*M,N.tail(5).sum()), RL_all(3*M,3), RT_all(3*M,3);

    for(int i = 0; i < M; i++){
        //ADVANCE
        t += dt;
        theta_b(i) = -angularPosition(t);
        for(int j = 0; j < num_grids; j++) grid[j].advanceTime();

        //SOLVE
        guess = solveSparseLevenbergMarquardt<objFunc<true> >(guess, 1e-9, 2500, 1e-7, 0.5, 1e-12, 1e-5, 1e-10);
        std::cout << "DYNAMIC: i=" << i+1 << " of " << M << std::endl;

        //PROCESS
        Matrix3d RL = quat2rot(Y[4].block<4,1>(3,N[4]-1));
        Matrix3d RT = quat2rot(Y[5].block<4,1>(3,N[5]-1));
        Vector3d diff = inv_hat(log(RL.transpose()*RT));
        theta_f(i) = Y[4](20, N[4]-1) - diff(2);
        p1.block(3*i,0,3,N.segment(1,4).sum()) << Y[1].topRows(3), Y[2].topRows(3), Y[3].topRows(3), Y[4].topRows(3);
        p0.block(3*i,0,3,N.tail(5).sum()) << p1.block(3*i,0,3,N.segment(1,4).sum()), Y[5].topRows(3);
        RL_all.block<3,3>(3*i,0) = RL;
        RT_all.block<3,3>(3*i,0) = RT;
    }

    std::fstream file;
    #define SAVE(VAR_NAME, FILE_NAME) \
        file = std::fstream(FILE_NAME, std::fstream::out); \
        file << VAR_NAME; \
        file.close();

    SAVE(p0, "p0.dat")
    SAVE(p1, "p1.dat")
    SAVE(beta0, "beta0.dat")
    SAVE(RL_all, "RL_all.dat")
    SAVE(RT_all, "RT_all.dat")
    SAVE(theta_f, "theta_f.dat")
    SAVE(theta_b, "theta_b.dat")

    #ifdef QT_CORE_LIB
    VectorXd t = VectorXd::LinSpaced(M,0,M*(dt*1e3));
    plot(t, theta_f, "Tip Angle", "t (ms)", "θ (rad)");
    #endif

    std::cout << "\nSolved tip angle trajectory:\n\ntheta_f=\n"
              << theta_f.transpose() << std::endl;

    return 0;
}
