#include "My_MatInverse_Top.hpp"


MInv_T my_mac(MInv_T x, MInv_T y, MInv_T z){

	MInv_T0 a = z.real();
	MInv_T0 b = z.imag();
	
//	std::cout<<a<<std::endl;
	a += x.real()*y.real();
	a -= x.imag()*y.imag();
	z.real(a);
	
	b += x.real()*y.imag();
	b += x.imag()*y.real();
//	std::cout<<b<<std::endl;
	z.imag(b);
//	std::cout<<b<<std::endl;

	return(z);
}




void mysqrt_real(MInv_T din, MInv_T& dout) {
Function_cholesky_sqrt_op_complex:;
    const MInv_T0 ZERO = 0;
    MInv_T0 a = din.real();
    dout.imag(ZERO);

//    if (a < ZERO) {
//        dout.real(ZERO);
//        return (1);
//    }
    dout.real(hls::sqrt(a));
//    return (0);
}

MInv_T myconj(MInv_T A){
	MInv_T C;
	C.real() = A.real();
	C.imag() = -A.imag();
	return C;
}

void mymult_complex_real(MInv_T A, MInv_T0 B, MInv_T& C) {
Function_cholesky_prod_sum_mult_complex:;
    C.real(A.real() * B);
    C.imag(A.imag() * B);
}


//------------------
// complex inv
// z = 1./x
//------------------
//MInv_T inv(MInv_T x){
//	MInv_T z;
//	float x2 = x.real*x.real + x.imag*x.imag;
//	double invx2 = 1./x2;
//	// MInvtiply x*
//	z.real = invx2*x.real;
//	z.imag = -invx2*x.imag;
//	return(z);
//}

