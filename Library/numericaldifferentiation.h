#ifndef NUMERICALDIFFERENTIATION_H
#define NUMERICALDIFFERENTIATION_H

#include "eigenmatrixtypes.h"
#include <forward_list>

namespace ContinuumRobotLibrary{

typedef VectorXd(ObjFunc)(VectorXd);

// Obtain an appropriately scaled increment for a number.
// There is a minimum increment in case the number being increment is too small, e.g. 0.
static double getIncrement(double val, double incr_scale, double incr_floor){
    double scaled = val*incr_scale;
    return (scaled > incr_floor || scaled < -incr_floor) ? scaled : incr_floor;
}

/*! Obtain the Jacobian of a function by finite differences.
 *  J is an output argument. The function value f(y) is usually already calculated,
 *  so it is an input to avoid redundant computation. */
template<ObjFunc f>
inline void getJacobian(MatrixXd& J,
                        VectorXd& y,
                        VectorXd& r,
                        double incr_scale = 1e-8,
                        double incr_floor = 1e-11){
    for(int i = 0; i < y.size(); i++){
        double temp = y(i);
        double incr = getIncrement(temp, incr_scale, incr_floor);
        y(i) += incr;
        J.col(i) = (f(y) - r) / incr; // + O(incr)
        y(i) = temp;
    }
}

template<ObjFunc f>
SparseMatrix<double> getSparseJacobian(VectorXd y,
                                       VectorXd r,
                                       double threshold = 1e-12,
                                       double incr_scale = 1e-8,
                                       double incr_floor = 1e-11){
    std::forward_list<Triplet<double> > triplet_list;

    for(int i = 0; i < y.size(); i++){
        double temp = y(i);
        double incr = getIncrement(temp, incr_scale, incr_floor);
        y(i) += incr;
        VectorXd Jcol = (f(y) - r) / incr; // + O(incr)
        y(i) = temp;

        for(int j = 0; j < r.size(); j++){
            if(Jcol(j) > threshold || Jcol(j) < -threshold)
                triplet_list.push_front( Triplet<double>(j,i,Jcol(j)) );
        }
    }

    SparseMatrix<double> J(r.size(), y.size());
    J.setFromTriplets(triplet_list.begin(), triplet_list.end());

    return J;
}

}

#endif // NUMERICALDIFFERENTIATION_H
