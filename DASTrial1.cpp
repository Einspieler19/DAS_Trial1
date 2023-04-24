
#include "DASTrial1.h"
#include "Template.h"



// 延迟时间计算函数
din_fix computeDelayTime(d_int_index elementIndex, din_angle azimuthAngle, din_angle elevationAngle) {
	d_fix sinAzimuthangle = CalculateSin<d_fix,d_fix>(azimuthAngle);
	d_fix cosElevationangle = CalculateCos<d_fix,d_fix>(elevationAngle);
	d_fix delayTime = ((d_fix)((elementIndex) * ELEMENTSPACING / SPEEDOFSOUND)) * sinAzimuthangle * cosElevationangle;
	//cout<<delayTime<<endl;
	return delayTime;
}

// 计算 延时向量组
void computeDelayTimesVector(d_fix* delayTimes, d_int_index numElements, din_angle azimuthAngle, din_angle elevationAngle) {
    for (int i = 0; i < numElements; i++) {
        delayTimes[i] = computeDelayTime(i, azimuthAngle, elevationAngle);
        //std::cout << delayTimes[i] << i<< " ";
    }
    //std::cout << std::endl;
}

// 平均权重
void computeWeights(d_fix* weights, d_int_index numElements) {
    for (int i = 0; i < numElements; i++) {
        weights[i] = 1.0 / numElements;
        //std::cout << weights[i] << " ";
    }
    //std::cout << std::endl; // 换行
}


// 权值*信号
d_fix mul_WeightSnapshot(d_fix weight, d_fix signalSnapshot) {
    return weight * signalSnapshot;
}

d_fix delay_Snapshots(d_fix elementDelay) {
    return (d_int)hls::round(elementDelay * MYSAMPLERATE);
}

// 权值*信号
d_fix sum_Snapshot(d_fix* weight, d_fix signal[SIGNALLENGTH], d_int_index snapshotIndex, d_fix* delayTimes) {
	d_fix multed = 0.0;
	d_fix summedSnapshot = 0.0;
	d_int_index delayedSnapshotIndex = 0;
	// 天线index和时间片index
    for (d_int_index elementIndex = 0; elementIndex < NUMELEMENTS; elementIndex++) {
    	d_int_index delays = delay_Snapshots(delayTimes[elementIndex]);
	    delayedSnapshotIndex = snapshotIndex - delays;
	    //cout<<delays<< "b " << snapshotIndex<< endl;
	    //std::cout << std::endl;
	    if (delayedSnapshotIndex >= 0 && delayedSnapshotIndex < SIGNALLENGTH) {
                //arraySignal += weights[j] * signal[signalIndex];
                multed = mul_WeightSnapshot(weight[elementIndex], signal[delayedSnapshotIndex]);
                //cout<<multed<< "multed " <<endl;
                //std::cout << std::endl;
                summedSnapshot += multed;
                //cout<<delayedSnapshotIndex<< "b " <<endl;
                //std::cout << std::endl;
            }
    }
    return summedSnapshot;
}

d_fix DASTrial1(d_fix signal[SIGNALLENGTH], d_int signalLength, din_angle azimuthAngle, din_angle elevationAngle, d_int_index snapshotIndex) {

		d_fix sampleRate = MYSAMPLERATE;
		d_fix speedOfSound = SPEEDOFSOUND;
		d_fix centerFreq = CENTERFREQ;
		d_fix arrayLength = ARRAYLENGTH;
		d_int numElements = NUMELEMENTS;
		d_fix elementSpacing = ELEMENTSPACING;
		//d_fix azimuthAngle = SIGNALANGLE;
		d_fix noiseVariance = NOISEVARIANCE;

	    d_fix weights[NUMELEMENTS];
	    d_fix delayTimes[NUMELEMENTS];
	    d_fix summedSignal = 0.0;

	    computeWeights(weights, numElements);
	    computeDelayTimesVector(delayTimes, numElements, azimuthAngle, elevationAngle);
	    summedSignal = sum_Snapshot(weights, signal, snapshotIndex, delayTimes);
        //cout<<summedSignal<< "summedSignal " <<endl;
        //std::cout << std::endl;
	    return summedSignal;
    }
