#include "search_eigenvalues.h"
#include <Spectra/GenEigsSolver.h>
#include <Spectra/MatOp/SparseGenMatProd.h>
#include <Spectra/MatOp/SparseGenRealShiftSolve.h>
#include <Spectra/GenEigsRealShiftSolver.h>
using namespace std;
using namespace Eigen;
using namespace Spectra;

namespace marvel{
int find_max_min_eigenvalues(const Eigen::SparseMatrix<double>& A, double& max_eigvalue, double& min_eigvalue){
    SparseGenMatProd<double> op(A);

    // Construct eigen solver object, requesting the largest three eigenvalues
    GenEigsSolver< double, LARGEST_MAGN, SparseGenMatProd<double> > eigs(&op, 1, 3);

    // Initialize and compute
    eigs.init();
    int nconv = eigs.compute();

    // Retrieve results
    Eigen::VectorXcd evalues;
    if(eigs.info() == SUCCESSFUL)
        evalues = eigs.eigenvalues();



    max_eigvalue = evalues(0).real();

    SparseGenRealShiftSolve<double> op_min(A);
    // Construct eigen solver object with shift 0
    // This will find eigenvalues that are closest to 0
    GenEigsRealShiftSolver< double, LARGEST_MAGN,
                        SparseGenRealShiftSolve<double> > eigs_min(&op_min, 1, 3, 0.0);

    eigs_min.init();
    eigs_min.compute();

    if(eigs_min.info() == SUCCESSFUL)
    {
        Eigen::VectorXd evalues = eigs_min.eigenvalues().real();
        min_eigvalue = evalues(0);
    }

    cout << "max eigvalue " << max_eigvalue << " min eigvalue " << min_eigvalue << endl;
  return 0;
}
}
