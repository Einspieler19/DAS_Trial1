
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
double sum_Snapshot(double* weight, double* signal, int snapshotIndex, double* delayTimes) {
	double summedSnapshot = 0.0;
	int delayedSnapshotIndex = 0;
	// 天线index和时间片index
    for (int elementIndex = 0; elementIndex < NUMELEMENTS; elementIndex++) {
	    delayedSnapshotIndex = snapshotIndex - delay_Snapshots(delayTimes[elementIndex]);
	    if (delayedSnapshotIndex >= 0 && delayedSnapshotIndex < SIGNALLENGTH) {
                //arraySignal += weights[j] * signal[signalIndex];
                summedSnapshot += mul_WeightSnapshot(weight[elementIndex], signal[delayedSnapshotIndex]);
            }
    }
    return summedSnapshot;
}

// 权值*信号
void sum_Signal(double* weight, double* signal, double* delayTimes, double* summedSignal) {

	double summedSignalAll = 0.0;
    for (int snapshotIndex = 0; snapshotIndex < SIGNALLENGTH; snapshotIndex++) {
	   summedSignal[snapshotIndex] = sum_Snapshot(weight, signal, snapshotIndex, delayTimes);
	   //summedSignalAll += summedSignal[snapshotIndex];
    }
    //return summedSignalAll;
}


// Delay and Sum Beamformer 函数
void DASTrial1(double* signal, int signalLength, double azimuthAngle, double elevationAngle, double *summedSignal) {

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

    computeWeights(weights, numElements);
    computeDelayTimesVector(delayTimes, numElements, azimuthAngle, elevationAngle);
    sum_Signal(weights, signal, delayTimes, summedSignal);

}





