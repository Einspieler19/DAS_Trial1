/*******************************************************************************
Vendor: Xilinx 
Associated Filename: hamming_window.h
Purpose: Vivado HLS tutorial example 
Revision History: March 1, 2013 - initial release
                : July 21, 2020 - 2020.1 Release
                                                
*******************************************************************************
Copyright (C) 2020 XILINX, Inc. 

This file contains confidential and proprietary information of Xilinx, Inc. and 
is protected under U.S. and international copyright and other intellectual 
property laws.

DISCLAIMER
This disclaimer is not a license and does not grant any rights to the materials 
distributed herewith. Except as otherwise provided in a valid license issued to 
you by Xilinx, and to the maximum extent permitted by applicable law: 
(1) THESE MATERIALS ARE MADE AVAILABLE "AS IS" AND WITH ALL FAULTS, AND XILINX 
HEREBY DISCLAIMS ALL WARRANTIES AND CONDITIONS, EXPRESS, IMPLIED, OR STATUTORY, 
INCLUDING BUT NOT LIMITED TO WARRANTIES OF MERCHANTABILITY, NON-INFRINGEMENT, OR 
FITNESS FOR ANY PARTICULAR PURPOSE; and (2) Xilinx shall not be liable (whether 
in contract or tort, including negligence, or under any other theory of 
liability) for any loss or damage of any kind or nature related to, arising under 
or in connection with these materials, including for any direct, or any indirect, 
special, incidental, or consequential loss or damage (including loss of data, 
profits, goodwill, or any type of loss or damage suffered as a result of any 
action brought by a third party) even if such damage or loss was reasonably 
foreseeable or Xilinx had been advised of the possibility of the same.

CRITICAL APPLICATIONS
Xilinx products are not designed or intended to be fail-safe, or for use in any 
application requiring fail-safe performance, such as life-support or safety 
devices or systems, Class III medical devices, nuclear facilities, applications 
related to the deployment of airbags, or any other applications that could lead 
to death, personal injury, or severe property or environmental damage 
(individually and collectively, "Critical Applications"). Customer assumes the 
sole risk and liability of any use of Xilinx products in Critical Applications, 
subject only to applicable laws and regulations governing limitations on product 
liability. 

THIS COPYRIGHT NOTICE AND DISCLAIMER MUST BE RETAINED AS PART OF THIS FILE AT 
ALL TIMES.

*******************************************************************************/
#ifndef DAS_H
#define DAS_H


#include <iostream>
#include <iomanip> 
#include <fstream>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include "constants.h"
#include <ap_int.h> // Xilinx hls的整型库
#include <ap_fixed.h>
#include "hls_math.h"


using namespace std;

#define FLOAT_DATA // Used to select error tolerance in test program

#define W_IN    48
#define IW_IN   24

typedef ap_uint<16> d_int_index;
typedef ap_int<16> d_int;
typedef ap_fixed<W_IN,IW_IN> din_angle;
typedef ap_fixed<W_IN,IW_IN> din_fix;
typedef ap_fixed<W_IN,IW_IN> d_fix;


d_fix DASTrial1(d_fix signal[SIGNALLENGTH], d_int signalLength, din_angle azimuthAngle, din_angle elevationAngle, d_int_index snapshotIndex);


#endif

