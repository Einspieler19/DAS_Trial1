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
#ifndef TEMPLATE_H
#define TEMPLATE_H

#include "DASTrial1.h"


template<typename type_angle, typename type_dout>
type_dout CalculateSin(type_angle angle)
{
    type_dout output = hls::sinf(angle);
    return output;
}

template<typename type_angle, typename type_dout>
type_dout CalculateCos(type_angle angle)
{
    type_dout output = hls::cosf(angle);
    return output;
}



template<typename din_fix, typename d_int_index, typename din_angle>
class DelayCalculator {
public:
    din_fix computeDelayTime(d_int_index elementIndex, din_angle azimuthAngle, din_angle elevationAngle) {
        din_fix sinAzimuthangle = CalculateSin(azimuthAngle);
        din_fix cosElevationangle = CalculateCos(elevationAngle);
        din_fix delayTime = ((din_fix)((elementIndex) * ELEMENTSPACING / SPEEDOFSOUND)) * sinAzimuthangle * cosElevationangle;
        return delayTime;
    }

private:
    din_fix CalculateSin(din_angle angle) {
        din_fix functionInput = angle;
        din_fix output = hls::sinf(functionInput);
        return output;
    }

    din_fix CalculateCos(din_angle angle) {
        din_fix functionInput = angle;
        din_fix output = hls::cosf(functionInput);
        return output;
    }


};

















#endif

