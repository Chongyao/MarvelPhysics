#ifndef MARVEL_PCG_GPU
#define MARVEL_PCG_GP



void readIAandJA(const int size_Matrix, const int size_nozeronumber, const int *IAtemp, const int *JAtemp);

template<typename DOUBLE>
class PCG {
 public:
  void function_pcg(const int Ntemp, const int NNZtemp ,const DOUBLE  *Atemp, const DOUBLE *Btemp, DOUBLE* x);
};

template class PCG<double>;
//TODO: make flaot avaiable
// template class PCG<float>;



// #include "pcg_gpu.cu"

#endif
