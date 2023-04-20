/*******************************************************************************
Vendor: Xilinx 
Associated Filename: hamming_window_test.c
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
#include <iostream>
#include "DASTrial1.h"
using namespace std;

////////////////////////// 确定对比精度 //////////////////////////
#ifdef FLOAT_DATA
#define ABS_ERR_THRESH 0.01
#else
#define ABS_ERR_THRESH 0.001
#endif


////////////////////////// 是否输出对比数据 //////////////////////////
#define WINDOW_FN_DEBUG 1

// 延迟时间计算函数
double computeDelayTime1(int elementIndex, double delayAngle) {
    double delayTime = ((elementIndex - 1) * ELEMENTSPACING / SPEEDOFSOUND) * sin(delayAngle);
    return delayTime;
}



// 生成信号函数
double* generateSignal(int length, double centerFreq, double sampleRate) {
    double* signal = new double[length];
    for (int i = 0; i < length; i++) {
        signal[i] = sin(2 * PI * centerFreq * i / sampleRate);
    }
    return signal;
}

// 添加噪声函数
void addNoise(double* signal, int length, double variance) {
    // 设置随机数种子
    srand((unsigned)time(NULL));

    for (int i = 0; i < length; i++) {
        signal[i] += sqrt(variance) * (double)rand() / RAND_MAX;
    }
}

// Delay and Sum Beamformer 函数
double beamform(double* signal, int signalLength, double delayAngle) {
    double* weights = new double[NUMELEMENTS];
    double delayTimes[NUMELEMENTS];



    // 计算各元素的延迟时间
    for (int i = 0; i < NUMELEMENTS; i++) {
        delayTimes[i] = computeDelayTime1(i, delayAngle);
    std::cout << delayTimes[i]  << " ";
}
std::cout << std::endl;

    // 计算各元素的权值
    for (int i = 0; i < NUMELEMENTS; i++) {
        weights[i] = 1.0 / NUMELEMENTS;
    std::cout << weights[i] << " ";
}
std::cout << std::endl;


    double sum = 0.0;
    for (int i = 0; i < signalLength; i++) {
        double signalValue = signal[i];
        double arraySignal = 0.0;
        for (int j = 0; j < NUMELEMENTS; j++) {
            double elementDelay = delayTimes[j];
            int delaySamples = (int)round(elementDelay * MYSAMPLERATE);
            int signalIndex = i - delaySamples;
            if (signalIndex >= 0 && signalIndex < signalLength) {
                arraySignal += weights[j] * signal[signalIndex];
            }
        }
        sum += arraySignal;
    }

    delete[] weights;

    return sum / signalLength;

}

int main(int argc, char *argv[])
{
	double sampleRate = MYSAMPLERATE;
	double speedOfSound = SPEEDOFSOUND;
	double centerFreq = CENTERFREQ;
	double arrayLength = ARRAYLENGTH;
	int numElements = NUMELEMENTS;
	double elementSpacing = ELEMENTSPACING;
	double signalAngle = SIGNALANGLE;
	double noiseVariance = NOISEVARIANCE;
   int signalLength = 48000;

   double hw_result[ILENGTH], sw_result[ILENGTH]; // HLS输出；软件输出
   int i;
   unsigned err_cnt = 0, check_dots = 0;
   FILE        *fp;


    double* signal = generateSignal(signalLength, CENTERFREQ, MYSAMPLERATE);

    // 添加噪声
    addNoise(signal, signalLength, NOISEVARIANCE);

    // 计算 Delay and Sum Beamformer 结果


    for (unsigned i = 0; i < ILENGTH; i++) {
    sw_result[i] = beamform(signal, signalLength, SIGNALANGLE);
    }

   ////////////////////////////////// 测试HLS程序 //////////////////////////////////
   for (unsigned i = 0; i < ILENGTH; i++) {
	   hw_result[i] = add(signal, signalLength, SIGNALANGLE);
   }

   ////////////////////////////////// 比对+写下结果 //////////////////////////////////


   // Check results
    cout << "Checking results against a tolerance of " << ABS_ERR_THRESH << endl;
    cout << fixed << setprecision(5);
   for (unsigned i = 0; i < ILENGTH; i++) {
      float abs_err = float(hw_result[i]) - sw_result[i];
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







