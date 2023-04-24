
#include "DASTrial1.h"


// 延迟时间计算函数
double computeDelayTime(int elementIndex, double azimuthAngle, double elevationAngle) {
    double delayTime = ((elementIndex) * ELEMENTSPACING / SPEEDOFSOUND) * sin(azimuthAngle)* cos(elevationAngle);
    return delayTime;
}

// 计算 延时向量组
void computeDelayTimesVector(double* delayTimes, int numElements, double azimuthAngle, double elevationAngle) {
    for (int i = 0; i < numElements; i++) {
        delayTimes[i] = computeDelayTime(i, azimuthAngle, elevationAngle);
        //std::cout << delayTimes[i] << i<< " ";
    }
    //std::cout << std::endl;
}

// 平均权重
void computeWeights(double* weights, int numElements) {
    for (int i = 0; i < numElements; i++) {
        weights[i] = 1.0 / numElements;
        //std::cout << weights[i] << " ";
    }
    //std::cout << std::endl; // 换行
}



// 权值*信号
double mul_WeightSnapshot(double weight, double signalSnapshot) {
    return weight * signalSnapshot;
}

double delay_Snapshots(double elementDelay) {
    return (int)round(elementDelay * MYSAMPLERATE);
}

// 权值*信号
double sum_Snapshot(double* weight, double signal[SIGNALLENGTH], int snapshotIndex, double* delayTimes) {
	double summed = 0.0;
	double summedSnapshot = 0.0;
	int delayedSnapshotIndex = 0;
	// 天线index和时间片index
    for (int elementIndex = 0; elementIndex < NUMELEMENTS; elementIndex++) {
	    delayedSnapshotIndex = snapshotIndex - delay_Snapshots(delayTimes[elementIndex]);
	    if (delayedSnapshotIndex >= 0 && delayedSnapshotIndex < SIGNALLENGTH) {
                //arraySignal += weights[j] * signal[signalIndex];
                summed = mul_WeightSnapshot(weight[elementIndex], signal[delayedSnapshotIndex]);
                summedSnapshot += summed;
            }
    }
    return summedSnapshot;
}

double DASTrial1(double signal[SIGNALLENGTH], int signalLength, double azimuthAngle, double elevationAngle, int snapshotIndex) {

		double sampleRate = MYSAMPLERATE;
		double speedOfSound = SPEEDOFSOUND;
		double centerFreq = CENTERFREQ;
		double arrayLength = ARRAYLENGTH;
		int numElements = NUMELEMENTS;
		double elementSpacing = ELEMENTSPACING;
		//double azimuthAngle = SIGNALANGLE;
		double noiseVariance = NOISEVARIANCE;

	    double weights[NUMELEMENTS];
	    double delayTimes[NUMELEMENTS];
	    double summedSignal = 0.0;

	    computeWeights(weights, numElements);
	    computeDelayTimesVector(delayTimes, numElements, azimuthAngle, elevationAngle);
	    summedSignal = sum_Snapshot(weights, signal, snapshotIndex, delayTimes);

	    return summedSignal;
    }
