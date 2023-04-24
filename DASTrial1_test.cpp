/*******************************************************************************
Vendor: Xilinx 
Associated Filename: hamming_window_test.c
Purpose: Vivado HLS tutorial example 
Revision History: March 1, 2013 - initial release
                : July 21, 2020 - 2020.1 Release
                                                
*******************************************************************************
Copyright (C) 2020 XILINX, Inc. 

This file contains confidential and proprietary information of Xilinx, Inc. and 
is protected under U.S. and d_international copyright and other d_intellectual
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
Xilinx products are not designed or d_intended to be fail-safe, or for use in any
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
#include <iostream>
#include "DASTrial1.h"
#include "Template.h"
using namespace std;

////////////////////////// 确定对比精度 //////////////////////////
#ifdef FLOAT_DATA
#define ABS_ERR_THRESH 0.1
#else
#define ABS_ERR_THRESH 0.001
#endif


////////////////////////// 是否输出对比数据 //////////////////////////
#define WINDOW_FN_DEBUG 1

d_fix Calculate_multiply(d_fix a, d_fix b, d_fix c)
{
	d_fix d = a*b*c;
	return d;
}


// 延迟时间计算函数
din_fix computeDelayTime1(d_int_index elementIndex, din_angle azimuthAngle, din_angle elevationAngle) {
	d_fix sinAzimuthangle = CalculateSin<d_fix,d_fix>(azimuthAngle);
	d_fix cosElevationangle = CalculateCos<d_fix,d_fix>(elevationAngle);
	d_fix delayTime = Calculate_multiply(((d_fix)((elementIndex) * ELEMENTSPACING / SPEEDOFSOUND)), sinAzimuthangle,cosElevationangle);
	return delayTime;
}


// 生成信号函数
d_fix* generateSignal(d_fix signal[SIGNALLENGTH], d_int length, d_fix centerFreq, d_fix sampleRate) {
    for (d_int i = 0; i < length; i++) {
    	d_fix sinSignal = CalculateSin<d_fix,d_fix>((din_angle)(2 * PI  * i) * centerFreq/ sampleRate);
        signal[i] = sinSignal;
    	}
    return signal;
}

// 添加噪声函数
void addNoise(d_fix* signal, d_int length, d_fix variance) {
    // 设置随机数种子
    srand((unsigned)time(NULL));
    d_fix randNum;
    ofstream FILE;

    // Save the results to a file
    FILE.open ("result_noised.dat");
    for (d_int i = 0; i < length; i++) {
    	randNum = (d_fix)(rand() / RAND_MAX-0.5);
        signal[i] += (d_fix)hls::sqrt(variance) * randNum;
        FILE << signal[i] << endl;
    }
    FILE.close();
}

// Delay and Sum Beamformer 函数
void beamform(d_fix* signal, d_int signalLength, din_angle azimuthAngle, din_angle elevationAngle, d_fix* output) {
    d_fix* weights = new d_fix[NUMELEMENTS];
    d_fix delayTimes[NUMELEMENTS];

    // 计算各元素的延迟时间
    for (d_int i = 0; i < NUMELEMENTS; i++) {
        delayTimes[i] = computeDelayTime1(i, azimuthAngle, elevationAngle);
    std::cout << delayTimes[i]  << " ";
}
std::cout << std::endl;

    // 计算各元素的权值
    for (d_int i = 0; i < NUMELEMENTS; i++) {
        weights[i] = 1.0 / NUMELEMENTS;
    std::cout << weights[i] << " ";
}
std::cout << std::endl;
    ofstream FILE;

    // Save the results to a file
    FILE.open ("result_formed.dat");

    d_fix sum = 0.0;
    for (d_int i = 0; i < signalLength; i++) {
        d_fix signalValue = signal[i];
        d_fix arraySignal = 0.0;
        for (d_int j = 0; j < NUMELEMENTS; j++) {
            d_fix elementDelay = delayTimes[j];
            d_int delaySamples = (d_int)hls::round(elementDelay * MYSAMPLERATE);

            d_int signalIndex = i - delaySamples;
            //cout<< signalIndex << endl;
            if (signalIndex >= 0 && signalIndex < signalLength) {
                arraySignal += weights[j] * signal[signalIndex];
                cout<< weights[j] << " "<< endl;
                cout<< signal[signalIndex] <<" "<< endl;
                }
        }
        output[i] = arraySignal;
        FILE << output[i] << endl;
    }
    delete[] weights;
    FILE.close();
}

int main(int argc, char *argv[])
{
	d_fix sampleRate = MYSAMPLERATE;
	d_fix speedOfSound = SPEEDOFSOUND;
	d_fix centerFreq = CENTERFREQ;
	d_fix arrayLength = ARRAYLENGTH;
	d_int numElements = NUMELEMENTS;
	d_fix elementSpacing = ELEMENTSPACING;
	d_fix azimuthAngle = SIGNALANGLE;
	d_fix elevationAngle = 0;
	d_fix noiseVariance = NOISEVARIANCE;
	d_int signalLength = SIGNALLENGTH;
	d_fix hw_result[SIGNALLENGTH], sw_result[SIGNALLENGTH]; // HLS输出；软件输出
	d_fix signal[SIGNALLENGTH];
	d_int i;
	unsigned err_cnt = 0, check_dots = 0;
	FILE        *fp;


    d_fix * signalPod_inter = generateSignal(signal, signalLength, CENTERFREQ, MYSAMPLERATE);

    addNoise(signalPod_inter, signalLength, NOISEVARIANCE);

    // 计算 Delay and Sum Beamformer 结果

    beamform(signalPod_inter, signalLength, azimuthAngle, elevationAngle, sw_result);
    // 添加噪声

   ////////////////////////////////// 测试HLS程序 //////////////////////////////////
   for (unsigned i = 0; i < SIGNALLENGTH; i++) {
    hw_result[i] = DASTrial1(signal, signalLength, azimuthAngle, elevationAngle,i);
   }

   ////////////////////////////////// 比对+写下结果 //////////////////////////////////

   // Check results
    cout << "Checking results against a tolerance of " << ABS_ERR_THRESH << endl;
    cout << fixed << setprecision(5);
   for (unsigned i = 0; i < SIGNALLENGTH; i++) {
      float abs_err = float(hw_result[i] - sw_result[i]);
#if WINDOW_FN_DEBUG
      cout << "i = " << i << "\thw_result = " << hw_result[i];
      cout << "\t sw_result = " << sw_result[i] << endl;
#endif
      if (fabs(abs_err) > ABS_ERR_THRESH) {
         cout << "Error threshold exceeded: i = " << i;
         cout << "  Expected: " << sw_result[i];
         cout << "  Got: " << hw_result[i];
         cout << "  Delta: " << abs_err << endl;
         err_cnt++;
      }
   }
   cout << endl;

   // Print final status message
   if (err_cnt) {
      cout << "!!! TEST FAILED - " << err_cnt;
      cout << " results out of tolerance." << endl;
   } else
      cout << "Test Passed" << endl;

   // Only return 0 on success
   return err_cnt;
}







