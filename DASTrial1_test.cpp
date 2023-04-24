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

// 延迟时间计算函数
double computeDelayTime1(int elementIndex, double azimuthAngle, double elevationAngle) {
    double delayTime = ((elementIndex) * ELEMENTSPACING / SPEEDOFSOUND) * sin(azimuthAngle)* cos(elevationAngle);
    //cout<<delayTime<<endl;
    return delayTime;
}


// 生成信号函数
void generateSignal(double* signal, int length, double centerFreq, double sampleRate) {
	//double signal[SIGNALLENGTH];
	for (int i = 0; i < SIGNALLENGTH; i++) {
        signal[i] = sin(2 * PI * centerFreq * i / sampleRate);
    }
    //return signal;
}


d_fix* generateSignal2(d_fix signal[SIGNALLENGTH], d_int length, d_fix centerFreq, d_fix sampleRate) {
    for (d_int i = 0; i < SIGNALLENGTH; i++) {
    	d_fix sinSignal = CalculateSin<d_fix,d_fix>((din_angle)(2 * PI  * i) * centerFreq/ sampleRate);
        signal[i] = sinSignal;
        //cout<<sinSignal<<endl;
        //cout<<endl;
    	}
    return signal;
}

void transferSignal(double* signalIn, d_fix* signalOut) {
	for (d_int i = 0; i < SIGNALLENGTH; i++) {
    	signalOut[i] = (d_fix)signalIn[i];
        //cout<<sinSignal<<endl;
        //cout<<endl;
    	}
    //return signalOut;
}

// 添加噪声函数
void addNoise(double* signal, int length, double variance) {
    // 设置随机数种子
    srand((unsigned)time(NULL));
    double randNum;
    ofstream FILE;

    // Save the results to a file
    FILE.open ("result_noised.dat");
    for (int i = 0; i < length; i++) {
    	randNum = (double)rand() / RAND_MAX-0.5;
        signal[i] += sqrt(variance) * randNum;
        FILE << signal[i] << endl;
    }
    FILE.close();
}

// Delay and Sum Beamformer 函数
void beamform(double* signal, int signalLength, double azimuthAngle, double elevationAngle, double* output) {
    double* weights = new double[NUMELEMENTS];
    double delayTimes[NUMELEMENTS];

    // 计算各元素的延迟时间
    for (int i = 0; i < NUMELEMENTS; i++) {
        delayTimes[i] = computeDelayTime1(i, azimuthAngle, elevationAngle);
    //std::cout << delayTimes[i]  << " ";
}
//std::cout << std::endl;

    // 计算各元素的权值
    for (int i = 0; i < NUMELEMENTS; i++) {
        weights[i] = 1.0 / NUMELEMENTS;
    //std::cout << weights[i] << " ";
}
//std::cout << std::endl;
    ofstream FILE;

    // Save the results to a file
    FILE.open ("result_formed.dat");

    double sum = 0.0;
    for (int i = 0; i < signalLength; i++) {
        double signalValue = signal[i];
        double arraySignal = 0.0;
        for (int j = 0; j < NUMELEMENTS; j++) {
            double elementDelay = delayTimes[j];
            int delaySamples = (int)round(elementDelay * MYSAMPLERATE);
            //cout<<delaySamples<< "a" << i<<endl;
            //std::cout << std::endl;
            int delayedSnapshotIndex = i - delaySamples;
            if (delayedSnapshotIndex >= 0 && delayedSnapshotIndex < signalLength) {
                arraySignal += weights[j] * signal[delayedSnapshotIndex];
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
	double sampleRate = MYSAMPLERATE;
	double speedOfSound = SPEEDOFSOUND;
	double centerFreq = CENTERFREQ;
	double arrayLength = ARRAYLENGTH;
	int numElements = NUMELEMENTS;
	double elementSpacing = ELEMENTSPACING;
	double azimuthAngle = SIGNALANGLE;
	double elevationAngle = 0;
	double noiseVariance = NOISEVARIANCE;
	int signalLength = SIGNALLENGTH;
	double hw_result[SIGNALLENGTH], sw_result[SIGNALLENGTH]; // HLS输出；软件输出
	double signalSW[SIGNALLENGTH];
	d_fix signalHW[SIGNALLENGTH];

	int i;
	unsigned err_cnt = 0, check_dots = 0;
	FILE        *fp;


    generateSignal(signalSW, signalLength, CENTERFREQ, MYSAMPLERATE);

    addNoise(signalSW, signalLength, NOISEVARIANCE);

	transferSignal(signalSW, signalHW);



    // 计算 Delay and Sum Beamformer 结果

    beamform(signalSW, signalLength, azimuthAngle, elevationAngle, sw_result);
    // 添加噪声

   ////////////////////////////////// 测试HLS程序 //////////////////////////////////
    ofstream FILE;

    // Save the results to a file
    FILE.open ("result_HW.dat");
    for (unsigned i = 0; i < SIGNALLENGTH; i++) {
    hw_result[i] = (double)DASTrial1(signalHW, (d_int)signalLength, (din_angle)azimuthAngle, (din_angle) elevationAngle, (d_int_index)i);
    FILE <<  hw_result[i] << endl;
   }


   FILE.close();

   ////////////////////////////////// 比对+写下结果 //////////////////////////////////

   // Check results
    cout << "Checking results against a tolerance of " << ABS_ERR_THRESH << endl;
    cout << fixed << setprecision(5);
   for (unsigned i = 0; i < SIGNALLENGTH; i++) {
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


