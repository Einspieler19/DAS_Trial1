///////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2016 Xilinx, Inc.
// All Rights Reserved
///////////////////////////////////////////////////////////////////////////////
//   ____  ____
//  /   /\/   /
// /___/  \  /    Vendor: Xilinx
// \   \   \/     Version: 1.0
//  \   \         Application : Vivado HLS
//  /   /         Filename: basic_math.c
// /___/   /\     Timestamp: Fri Aug 19 8:08:10 PST 2016
// \   \  /  \
//  \___\/\___\
//
//Command: N/A
//Device:  Xilinx Kintex-7, US, US+; Zynq-7000, US+
//Design Name: LUInv
//Purpose:
//    This file contains basic math functions, e.g., square, add, minus, etc.
//Reference:
//    XAPP1317
///////////////////////////////////////////////////////////////////////////////

#include "My_MatInverse_Top.hpp"

//------------------
// complex square
// z = |x|^2
//------------------
MInv_T0 csqu(MInv_T x){
	MInv_T0 z = x.real()*x.real();
	z+= x.imag()*x.imag();
	return(z);
}

MInv_T_RD0 csquf(MInv_T_RD x){
	MInv_T_RD0 z = x.real()*x.real();
	z+= x.imag()*x.imag();
	return(z);
}

//------------------
// complex neg
// z = -x
//------------------
MInv_T neg(MInv_T x){
	MInv_T z;
	z.real(-x.real());
	z.imag(-x.imag());
	return(z);
}

//------------------
// complex add
// z = x + y
//------------------
MInv_T add(MInv_T x, MInv_T y){
	MInv_T z;
	z.real(x.real() + y.real());
	z.imag(x.imag() + y.imag());
	return(z);
}

//------------------
// complex subtraction
// z = x - y
//------------------
MInv_T sub(MInv_T x, MInv_T y){
	MInv_T z;
	z.real(x.real() - y.real());
	z.imag(x.imag() - y.imag());
	return(z);
}

//------------------
// complex MInvtiply
// z = x * y
//------------------
MInv_T MInv(MInv_T x, MInv_T y){
	MInv_T z;
	z.real(x.real()*y.real() - x.imag()*y.imag());
	z.imag(x.real()*y.imag() + x.imag()*y.real());
	return(z);
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

