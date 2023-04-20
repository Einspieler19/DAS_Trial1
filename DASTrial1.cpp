
#include "DASTrial1.h"


// 延迟时间计算函数
double computeDelayTime(int elementIndex, double delayAngle) {
    double delayTime = ((elementIndex - 1) * ELEMENTSPACING / SPEEDOFSOUND) * sin(delayAngle);
    return delayTime;
}

// Delay and Sum Beamformer 函数
double DASTrial1(double* signal, int signalLength, double delayAngle) {

	double sampleRate = MYSAMPLERATE;
	double speedOfSound = SPEEDOFSOUND;
	double centerFreq = CENTERFREQ;
	double arrayLength = ARRAYLENGTH;
	int numElements = NUMELEMENTS;
	double elementSpacing = ELEMENTSPACING;
	double signalAngle = SIGNALANGLE;
	double noiseVariance = NOISEVARIANCE;

    double weights[NUMELEMENTS];
    double delayTimes[NUMELEMENTS];

    // 计算各元素的延迟时间
    for (int i = 0; i < NUMELEMENTS; i++) {
        delayTimes[i] = computeDelayTime(i, delayAngle);
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

    //delete[] weights;

    double output = sum / signalLength;
    return output;
}




