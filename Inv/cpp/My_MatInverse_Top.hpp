

#ifndef _MYMINV_H
#define _MYMINV_H

#include <cmath>
#include <complex>
#include "hls_stream.h"
#include <hls_math.h>
#include "hls_x_complex.h"
#include "ap_fixed.h"
using namespace std;

// number of layers, or matrix size
#define NUM_Tests 1
#define NL 32
#define RowsA NL 
#define ColsA NL 
#define RowsB NL 
#define ColsB NL 

#define BLOCK_SIZE NL/2

#define RowsColsA NL



// floating point complex data type
//typedef struct complex_float {
//  float real;
//  float imag;
//} MInv_T;

typedef float MInv_T0;
typedef float MInv_TOut0;
//
//typedef ap_fixed<18,6> MInv_T0;
//typedef ap_fixed<18,6> MInv_TOut0;
//typedef ap_fixed<36,12> MInv_TOut0;

typedef hls::x_complex<MInv_T0> MInv_T;
typedef hls::x_complex<MInv_TOut0> MInv_TOut;

typedef float MInv_T_RD0;
typedef hls::x_complex<MInv_T_RD0> MInv_T_RD;










//typedef hls::x_complex<int> MInv_T;



typedef MInv_T PROD_T;
typedef MInv_T ACCUM_T;

typedef MInv_T ADD_T;
typedef MInv_T    DIAG_T ;         // sqrt(A_minus_sum)
typedef MInv_T  OFF_DIAG_T; // sum/L










// complex square
MInv_T0 csqu(MInv_T x);
MInv_T_RD0 csquf(MInv_T_RD x);
// complex neg
MInv_T neg(MInv_T x);
// complex add
MInv_T add(MInv_T x, MInv_T y);
// complex sub
MInv_T sub(MInv_T x, MInv_T y);
// complex inv
MInv_T inv(MInv_T x);
// complex MInvtiply
MInv_T MInv(MInv_T x, MInv_T y);

MInv_T my_mac(MInv_T x, MInv_T y, MInv_T z);


void mysqrt_real(MInv_T din, MInv_T& dout);

MInv_T myconj(MInv_T A);

void mymult_complex_real(MInv_T A, MInv_T0 B, MInv_T& C);

// My_MatInverses
void My_MatInverse_Baseline(hls::stream<MInv_T>& matrixInStrm,
		hls::stream<MInv_TOut>& matrixOutStrm);


void Vitis_chol_Basic(hls::stream<MInv_T>& matrixInStrm,
		hls::stream<MInv_TOut>& matrixOutStrm);


void Vitis_chol_Alt1(hls::stream<MInv_T>& matrixInStrm,
		hls::stream<MInv_TOut>& matrixOutStrm);

void Vitis_chol_Alt2(hls::stream<MInv_T>& matrixInStrm,
		hls::stream<MInv_TOut>& matrixOutStrm);


void My_LDL(hls::stream<MInv_T>& matrixInStrm,
		hls::stream<MInv_TOut>& matrixOutStrm);









// The main function
void My_MatInverse_Top(hls::stream<MInv_T>& matrixInStrm,
		hls::stream<MInv_TOut>& matrixOutStrm);


#endif
