#ifndef NUMERICALPDE_H
#define NUMERICALPDE_H

#include "eigenmatrixtypes.h"

namespace ContinuumRobotLibrary{

//A function pointer typedef for [y_s,z] = f(y)
typedef void(AutonomousOdeZFunc)(VectorXd&,Map<VectorXd>,Map<VectorXd>); //[y_s,z] = f(y)

/*! Integrate an ODE using Euler's method, while setting a z vector to solve initial conditions. */
template<AutonomousOdeZFunc ODE, int statesize, int derivsize, int N = 100>
inline void euler(
        MatrixXd& Y,
        MatrixXd& Z,
        VectorXd y0,
        double L
        ){
    Y.col(0) = y0;

    double ds = L/(N-1);

    //Euler's method
    #define y_s y0
    for(int i = 0; i < N-1; i++){
        ODE(y_s, Map<VectorXd>(&Z(0,i),derivsize), Map<VectorXd>(&Y(0,i),statesize));
        Y.col(i+1) = Y.col(i) + ds*y_s;
    }
    //Set Z at final grid point
    ODE(y_s, Map<VectorXd>(&Z(0,N-1),derivsize), Map<VectorXd>(&Y(0,N-1),statesize));
    #undef y_s
}

/*! Integrate an ODE using Euler's method, while setting a z vector to solve initial conditions. */
template<AutonomousOdeZFunc ODE, int statesize, int derivsize>
inline void euler(
		int N,
        MatrixXd& Y,
        MatrixXd& Z,
        VectorXd y0,
        double L
        ){
    Y.col(0) = y0;

    double ds = L/(N-1);

    //Euler's method
    #define y_s y0
    for(int i = 0; i < N-1; i++){
        ODE(y_s, Map<VectorXd>(&Z(0,i),derivsize), Map<VectorXd>(&Y(0,i),statesize));
        Y.col(i+1) = Y.col(i) + ds*y_s;
    }
    //Set Z at final grid point
    ODE(y_s, Map<VectorXd>(&Z(0,N-1),derivsize), Map<VectorXd>(&Y(0,N-1),statesize));
    #undef y_s
}

//A function pointer typedef for [y_s,z] = f(y,z_h)
typedef void(AutonomousPdeFunc)(VectorXd&,Map<VectorXd>,Map<VectorXd>,Map<VectorXd>);

/*! Integrate a semi-discretized PDE using Euler's method.
    The PDE function must set time-differientiated variables in the z vector,
    and the history terms z_h are calculated elsewhere. */
template<AutonomousPdeFunc PDE, int statesize, int derivsize, int N = 100>
inline void TimeBdfAlpha_SpaceEuler(
        MatrixXd& Y,
        MatrixXd& Z,
        VectorXd y0,
        double L,
        MatrixXd& Z_h
        ){
    Y.col(0) = y0;

    double ds = L/(N-1);

    //Euler's method
    #define y_s y0
    for(int i = 0; i < N-1; i++){
        PDE(y_s, Map<VectorXd>(&Z(0,i),derivsize), Map<VectorXd>(&Y(0,i),statesize), Map<VectorXd>(&Z_h(0,i),derivsize));
        Y.col(i+1) = Y.col(i) + ds*y_s;
    }
    //Set Z at final grid point
    PDE(y_s, Map<VectorXd>(&Z(0,N-1),derivsize), Map<VectorXd>(&Y(0,N-1),statesize), Map<VectorXd>(&Z_h(0,N-1),derivsize));
    #undef y_s
}

template<AutonomousPdeFunc PDE, int statesize, int derivsize>
inline void TimeBdfAlpha_SpaceEuler(
		int N,
        MatrixXd& Y,
        MatrixXd& Z,
        VectorXd y0,
        double L,
        MatrixXd& Z_h
        ){
    Y.col(0) = y0;

    double ds = L/(N-1);

    //Euler's method
    #define y_s y0
    for(int i = 0; i < N-1; i++){
        PDE(y_s, Map<VectorXd>(&Z(0,i),derivsize), Map<VectorXd>(&Y(0,i),statesize), Map<VectorXd>(&Z_h(0,i),derivsize));
        Y.col(i+1) = Y.col(i) + ds*y_s;
    }
    //Set Z at final grid point
    PDE(y_s, Map<VectorXd>(&Z(0,N-1),derivsize), Map<VectorXd>(&Y(0,N-1),statesize), Map<VectorXd>(&Z_h(0,N-1),derivsize));
    #undef y_s
}


//A function pointer typedef for [y_s,z] = f(s,y)
typedef void(OdeZFunc)(VectorXd&,Map<VectorXd>,double,Map<VectorXd>); //[y_s,z] = f(s,y)

/*! Integrate an ODE using Euler's method, while setting a z vector to solve initial conditions. */
template<OdeZFunc ODE, int statesize, int derivsize, int N = 100>
inline void euler(
        MatrixXd& Y,
        MatrixXd& Z,
        VectorXd y0,
        double s0,
        double sf
        ){
    Y.col(0) = y0;

    const int Nm1 = N-1;
    double ds = (sf-s0)/Nm1;

    //Euler's method
    #define y_s y0
    for(int i = 0; i < Nm1; i++){
        double s = s0 + ds*i;
        ODE(y_s, Map<VectorXd>(&Z(0,i),derivsize), s, Map<VectorXd>(&Y(0,i),statesize));
        Y.col(i+1) = Y.col(i) + ds*y_s;
    }
    //Set Z at final grid point
    ODE(y_s, Map<VectorXd>(&Z(0,Nm1),derivsize), sf, Map<VectorXd>(&Y(0,Nm1),statesize));
    #undef y_s
}

//A function pointer typedef for [y_s,z] = f(s,y,z_h)
typedef void(PdeFunc)(VectorXd&,Map<VectorXd>,double,Map<VectorXd>,Map<VectorXd>);

/*! Integrate a semi-discretized PDE using Euler's method.
    The PDE function must set time-differientiated variables in the z vector,
    and the history terms z_h are calculated elsewhere. */
template<PdeFunc PDE, int statesize, int derivsize, int N = 100>
inline void TimeBdfAlpha_SpaceEuler(
        MatrixXd& Y,
        MatrixXd& Z,
        VectorXd y0,
        double s0,
        double sf,
        MatrixXd& Z_h
        ){
    Y.col(0) = y0;

    const int Nm1 = N-1;
    double ds = (sf-s0)/Nm1;

    //Euler's method
    #define y_s y0
    for(int i = 0; i < Nm1; i++){
        double s = s0 + ds*i;
        PDE(y_s, Map<VectorXd>(&Z(0,i),derivsize), s, Map<VectorXd>(&Y(0,i),statesize), Map<VectorXd>(&Z_h(0,i),derivsize));
        Y.col(i+1) = Y.col(i) + ds*y_s;
    }
    //Set Z at final grid point
    PDE(y_s, Map<VectorXd>(&Z(0,Nm1),derivsize), sf, Map<VectorXd>(&Y(0,Nm1),statesize), Map<VectorXd>(&Z_h(0,Nm1),derivsize));
    #undef y_s
}


/*! Integrate an ODE using Euler's method, while setting a z vector to solve initial conditions. */
template<OdeZFunc ODE, int statesize, int derivsize, int N = 100>
inline void RK4(
        MatrixXd& Y,
        MatrixXd& Z,
        VectorXd y0,
        double s0,
        double sf
        ){
    Y.col(0) = y0;

    const int Nm1 = N-1;
    double ds = (sf-s0)/Nm1;
    double half_ds = ds/2;
    double sixth_ds = ds/6;

    //Euler's method
    #define k1 y0
    VectorXd k2(statesize), k3(statesize), k4(statesize), y_stage(statesize), throwaway(derivsize);
    for(int i = 0; i < Nm1; i++){
        double s = s0 + ds*i;
        ODE(k1, Map<VectorXd>(&Z(0,i),derivsize), s, Map<VectorXd>(&Y(0,i),statesize));
        y_stage = Y.col(i) + half_ds*k1;
        ODE(k2, Map<VectorXd>(throwaway.data(), derivsize), s+half_ds, Map<VectorXd>(y_stage.data(), statesize));
        y_stage = Y.col(i) + half_ds*k2;
        ODE(k3, Map<VectorXd>(throwaway.data(), derivsize), s+half_ds, Map<VectorXd>(y_stage.data(), statesize));
        y_stage = Y.col(i) + ds*k3;
        ODE(k4, Map<VectorXd>(throwaway.data(), derivsize), s+ds, Map<VectorXd>(y_stage.data(), statesize));
        Y.col(i+1) = Y.col(i) + sixth_ds*(k1 + 2*(k2 + k3) + k4);
    }
    //Set Z at final grid point
    ODE(k1, Map<VectorXd>(&Z(0,Nm1),derivsize), sf, Map<VectorXd>(&Y(0,Nm1),statesize));
    #undef k1
}



template<OdeZFunc ODE, int statesize, int derivsize>
inline void RK4(
        int N,
        MatrixXd& Y,
        MatrixXd& Z,
        VectorXd y0,
        double s0,
        double sf
        ){
    Y.col(0) = y0;

    const int Nm1 = N-1;
    double ds = (sf-s0)/Nm1;
    double half_ds = ds/2;
    double sixth_ds = ds/6;

    // Euler's method
    #define k1 y0
    VectorXd k2(statesize), k3(statesize), k4(statesize), y_stage(statesize), throwaway(derivsize);
    for(int i = 0; i < Nm1; i++){
        double s = s0 + ds*i;
        ODE(k1, Map<VectorXd>(&Z(0,i),derivsize), s, Map<VectorXd>(&Y(0,i),statesize));
        y_stage = Y.col(i) + half_ds*k1;
        ODE(k2, Map<VectorXd>(throwaway.data(), derivsize), s+half_ds, Map<VectorXd>(y_stage.data(), statesize));
        y_stage = Y.col(i) + half_ds*k2;
        ODE(k3, Map<VectorXd>(throwaway.data(), derivsize), s+half_ds, Map<VectorXd>(y_stage.data(), statesize));
        y_stage = Y.col(i) + ds*k3;
        ODE(k4, Map<VectorXd>(throwaway.data(), derivsize), s+ds, Map<VectorXd>(y_stage.data(), statesize));
        Y.col(i+1) = Y.col(i) + sixth_ds*(k1 + 2*(k2 + k3) + k4);
    }
    // Set Z at final grid point
    ODE(k1, Map<VectorXd>(&Z(0,Nm1),derivsize), sf, Map<VectorXd>(&Y(0,Nm1),statesize));
    #undef k1
}

template<AutonomousOdeZFunc ODE, int statesize, int derivsize>
inline void RK4(
        int N,
        MatrixXd& Y,
        MatrixXd& Z,
        VectorXd y0,
        double L
        ){
    Y.col(0) = y0;

    const int Nm1 = N-1;
    double ds = L/Nm1;
    double half_ds = ds/2;
    double sixth_ds = ds/6;

    //Euler's method
    #define k1 y0
    VectorXd k2(statesize), k3(statesize), k4(statesize), y_stage(statesize), throwaway(derivsize);
    for(int i = 0; i < Nm1; i++){
        ODE(k1, Map<VectorXd>(&Z(0,i),derivsize), Map<VectorXd>(&Y(0,i),statesize));
        y_stage = Y.col(i) + half_ds*k1;
        ODE(k2, Map<VectorXd>(throwaway.data(), derivsize), Map<VectorXd>(y_stage.data(), statesize));
        y_stage = Y.col(i) + half_ds*k2;
        ODE(k3, Map<VectorXd>(throwaway.data(), derivsize), Map<VectorXd>(y_stage.data(), statesize));
        y_stage = Y.col(i) + ds*k3;
        ODE(k4, Map<VectorXd>(throwaway.data(), derivsize), Map<VectorXd>(y_stage.data(), statesize));
        Y.col(i+1) = Y.col(i) + sixth_ds*(k1 + 2*(k2 + k3) + k4);
    }
    //Set Z at final grid point
    ODE(k1, Map<VectorXd>(&Z(0,Nm1),derivsize), Map<VectorXd>(&Y(0,Nm1),statesize));
    #undef k1
}


/*! Integrate a semi-discretized PDE using Euler's method.
    The PDE function must set time-differientiated variables in the z vector,
    and the history terms z_h are calculated elsewhere. */
template<PdeFunc PDE, int statesize, int derivsize, int N = 100>
inline void TimeBdfAlpha_SpaceRK4(
        MatrixXd& Y,
        MatrixXd& Z,
        VectorXd y0,
        double s0,
        double sf,
        MatrixXd& Z_h
        ){
    Y.col(0) = y0;

    const int Nm1 = N-1;
    double ds = (sf-s0)/Nm1;
    double half_ds = ds/2;
    double sixth_ds = ds/6;

    MatrixXd Z_h_interp = (Z_h.leftCols(Nm1) + Z_h.rightCols(Nm1))/2;

    //Euler's method
    #define k1 y0
    VectorXd k2(statesize), k3(statesize), k4(statesize), y_stage(statesize), throwaway(derivsize);
    for(int i = 0; i < Nm1; i++){
        double s = s0 + ds*i;
        PDE(k1, Map<VectorXd>(&Z(0,i),derivsize), s, Map<VectorXd>(&Y(0,i),statesize), Map<VectorXd>(&Z_h(0,i),derivsize));
        y_stage = Y.col(i) + half_ds*k1;
        PDE(k2, Map<VectorXd>(throwaway.data(), derivsize), s+half_ds, Map<VectorXd>(y_stage.data(), statesize), Map<VectorXd>(&Z_h_interp(0,i),derivsize));
        y_stage = Y.col(i) + half_ds*k2;
        PDE(k3, Map<VectorXd>(throwaway.data(), derivsize), s+half_ds, Map<VectorXd>(y_stage.data(), statesize), Map<VectorXd>(&Z_h_interp(0,i),derivsize));
        y_stage = Y.col(i) + ds*k3;
        PDE(k4, Map<VectorXd>(throwaway.data(), derivsize), s+ds, Map<VectorXd>(y_stage.data(), statesize), Map<VectorXd>(&Z_h(0,i+1),derivsize));
        Y.col(i+1) = Y.col(i) + sixth_ds*(k1 + 2*(k2+k3) + k4);
    }
    //Set Z at final grid point
    PDE(k1, Map<VectorXd>(&Z(0,Nm1),derivsize), sf, Map<VectorXd>(&Y(0,Nm1),statesize), Map<VectorXd>(&Z_h(0,Nm1),derivsize));
    #undef k1
}

// TimeBdfAlpha_SpaceRK4<lowerPDE<0>,3,2>(N[0],Y[0],Z[0],y_beta0, L[0], 0, Z_h[0]);
/*! Integrate a semi-discretized PDE using Euler's method.
    The PDE function must set time-differientiated variables in the z vector,
    and the history terms z_h are calculated elsewhere. */
	
	
template<AutonomousPdeFunc PDE, int statesize, int derivsize>
inline void TimeBdfAlpha_SpaceRK4(
		int N,
        MatrixXd& Y,
        MatrixXd& Z,
        VectorXd y0,
        double L, 
        MatrixXd& Z_h
        ){
    Y.col(0) = y0;

    const int Nm1 = N-1;
    double ds = L/Nm1;
    double half_ds = ds/2;
    double sixth_ds = ds/6;

    MatrixXd Z_h_interp = (Z_h.leftCols(Nm1) + Z_h.rightCols(Nm1))/2;

    //Euler's method
    #define k1 y0
    VectorXd k2(statesize), k3(statesize), k4(statesize), y_stage(statesize), throwaway(derivsize);
    for(int i = 0; i < Nm1; i++){
    
        PDE(k1, Map<VectorXd>(&Z(0,i),derivsize), Map<VectorXd>(&Y(0,i),statesize), Map<VectorXd>(&Z_h(0,i),derivsize));
        y_stage = Y.col(i) + half_ds*k1;
        PDE(k2, Map<VectorXd>(throwaway.data(), derivsize), Map<VectorXd>(y_stage.data(), statesize), Map<VectorXd>(&Z_h_interp(0,i),derivsize));
        y_stage = Y.col(i) + half_ds*k2;
        PDE(k3, Map<VectorXd>(throwaway.data(), derivsize), Map<VectorXd>(y_stage.data(), statesize), Map<VectorXd>(&Z_h_interp(0,i),derivsize));
        y_stage = Y.col(i) + ds*k3;
        PDE(k4, Map<VectorXd>(throwaway.data(), derivsize), Map<VectorXd>(y_stage.data(), statesize), Map<VectorXd>(&Z_h(0,i+1),derivsize));
        Y.col(i+1) = Y.col(i) + sixth_ds*(k1 + 2*(k2+k3) + k4);
    }
    //Set Z at final grid point
    PDE(k1, Map<VectorXd>(&Z(0,Nm1),derivsize), Map<VectorXd>(&Y(0,Nm1),statesize), Map<VectorXd>(&Z_h(0,Nm1),derivsize));
    #undef k1
}


}

#endif // NUMERICALPDE_H
