/*******************************************************************
Describes: 	routines for detection/score of wake, sleep and REM
   Author:	Jango
  Version: 	v1.1
     Date: 	2017-9-19
Copyright:	EEGSmart
********************************************************************/
#include "sleepWakeDetection.h"



#define MAXHOUR				(20)
#define SECONDS_OF_EPOCH 	(30)
#define EPOCH_OF_HOUR		(120)
#define NFFT				(256)
#define RAW_SAMPLE_FS       (256)
//#define FFT50				(50)
#define NaN					(-1)
#define ORIGINAL_THRESHOLD	(3)
#define ORIGINAL_THRESHOLD2	(3)
#define REM					(5)
#define WAKE				(4)
#define LIGHT				(3)
#define DEEP				(2)
#define WINDOWS				(5)
#define MAXENPOCH 			(MAXHOUR*EPOCH_OF_HOUR)
#define MAXSECOND			(SECONDS_OF_EPOCH*EPOCH_OF_HOUR*MAXHOUR)

typedef struct{
	int totalEpoch;
	int score[MAXENPOCH];
	int epochValid[MAXENPOCH];
	float epoch[MAXENPOCH];
	float epoch2[MAXENPOCH];
}outAASM;


static void rtRemoveDC(int *rawdata)
{
	int i;
	int j;
	double sum0[SECONDS_OF_EPOCH] = {0};
	double sum;

	sum = 0;

	for(i=0; i<SECONDS_OF_EPOCH; i++)
	{
		for(j=0; j<RAW_SAMPLE_FS; j++)
			sum0[i] += 1.0*rawdata[i*RAW_SAMPLE_FS+j];
		sum0[i] /= RAW_SAMPLE_FS;
	}

	for(i=0; i<SECONDS_OF_EPOCH; i++)
		sum += sum0[i];
	sum /= SECONDS_OF_EPOCH;

	for(i=0; i<SECONDS_OF_EPOCH*RAW_SAMPLE_FS; i++)
		rawdata[i] -= (int)sum;

	return;
}


static void judge(float *epochFilted2, float *th, float *epoch2Filted2, float *th2, int *epochValid, int numEpoch, int *score)
{
	int i;
	int haveSlept;

	haveSlept = 0;
	for(i=0; i<numEpoch; i++)
	{
		if(epochValid[i]==NaN)
			score[i] = NaN;
		else
		{
			if(epochFilted2[i]<=th[i])
			{
				if(haveSlept==0)
					score[i] = WAKE;
				else
					score[i] = REM;
			}
			else
			{
				if(epoch2Filted2[i]<=th2[i])
					score[i] = DEEP;
				else
					score[i] = LIGHT;

				haveSlept = 1;
			}
		}
	}

	while(score[i-1]==REM)
		score[i---1] = WAKE;

	return;
}


static void rtClearRaw(int *rawdata, float threshold, int *indClear)
{
	int i;
	int j;
	int offset;
	int index;

	index = 0;
	for(i=0; i<SECONDS_OF_EPOCH; i++)
	{
		indClear[index] = 1;
		for(j=0; j<RAW_SAMPLE_FS; j++)
		{
			offset = i*RAW_SAMPLE_FS + j;
			if((rawdata[offset]>threshold)||(rawdata[offset]<-1*threshold))
			{
				indClear[index] = NaN;
				break;
			}
		}
		index++;
	}

	return;
}


static float mean(float *data, int num)
{
	int i;
	float sum;
	float result;

	sum = 0;
	for(i=0; i<num; i++)
		sum += data[i];
	result = sum/num;

	return result;
}

static float ratioFreqBand(float *absFft)
{
	int numAlpha;
	int numBeta;
	const int alphaLow = 9;
	const int alphaHigh = 14;
	const int betaLow = 17;
	const int betaHigh = 22;
	float meanAlpha;
	float meanBeta;
	float ratioBand;

	numAlpha = alphaHigh - alphaLow + 1;
	numBeta = betaHigh - betaLow + 1;
	meanAlpha = mean(&absFft[alphaLow],numAlpha);
	meanBeta = mean(&absFft[betaLow],numBeta);
	ratioBand = meanAlpha/meanBeta;

	return ratioBand;
}

static float ratioFreqBand2(float *absFft)
{
	int numAlpha;
	int numDelta;
	const int alphaLow = 9;
	const int alphaHigh = 14;
	const int deltaLow = 1;
	const int deltaHigh = 3;
	float meanAlpha;
	float meanDelta;
	float ratioBand;

	numAlpha = alphaHigh - alphaLow + 1;
	numDelta = deltaHigh - deltaLow + 1;
	meanAlpha = mean(&absFft[alphaLow],numAlpha);
	meanDelta = mean(&absFft[deltaLow],numDelta);
	ratioBand = meanAlpha/meanDelta;

	return ratioBand;
}

static void averagingFilter(float *in, int length, int window, float *out)
{
	int i;
	int j;
	int midpoint;
	float tmp;

	midpoint = window/2 + 1;

	for(i=0; i<midpoint-1; i++)
		out[i] = in[i];

	for(i=midpoint-1; i<length-midpoint+1; i++)
	{
		tmp = 0;
		for(j=i-midpoint; j<i-midpoint+window; j++)
			tmp += in[j];
		out[i] = tmp/window;
	}

	for(i=length-midpoint+1; i<length; i++)
		out[i] = in[i];

	return;
}


static void thSleepWake(float *ratio, int length, float *th, float gain)
{
	int i;
	const int window = 5;
	const int startPoint = 1;
	int node[2] = {0,startPoint-1};
	int flag;

	const float coef1 = 0.5;
	const float coef2 = 0.5;

	for(i=0; i<startPoint; i++)
		th[i] = gain;

	flag = 0;
	for(i=startPoint; i<length; i++)
	{
		if(((th[i-1]>=ratio[i-1])&&(th[i-1]<=ratio[i]))||((th[i-1]<=ratio[i-1])&&(th[i-1]>=ratio[i])))
		{
			flag = 1;
			node[0] = node[1];
			node[1] = i;

			if(node[1]-node[0]>=window)
				th[i] = coef1*mean(&ratio[node[0]],node[1]-node[0]) + coef2*ratio[i];
			else
				th[i] = th[i-1];
		}
		else
		{
			if(flag==0)
				th[i] = th[i-1];
			else
				th[i] = coef1*mean(&ratio[node[0]],node[1]-node[0]) + coef2*mean(&ratio[node[1]],i-node[1]+1);
		}
	}

	return;
}

static void rtEpochCal(float *absFft, int *indClear, float *epoch, float *epoch2, int *epochValid, int *totalEpoch)
{
	int j;
	int tmp;
	int numValid;

	tmp = 0;
	numValid = 0;
	*epoch = 0;
	*epoch2 = 0;
	for(j=0; j<SECONDS_OF_EPOCH; j++)
	{
		tmp += indClear[j];
		if(indClear[j]!=NaN)
		{
			*epoch += ratioFreqBand(absFft+j*NFFT);
			*epoch2 += ratioFreqBand2(absFft+j*NFFT);
			numValid++;
		}
	}

	if(tmp==SECONDS_OF_EPOCH*NaN)
		*epochValid = NaN;
	else
		*epochValid = 1;

	if(numValid==0)
	{
		if(totalEpoch[0]>0)
			*epoch = *(epoch-1);
		else
			*epoch = ORIGINAL_THRESHOLD;
	}
	else
		*epoch /= numValid;

	totalEpoch[0]++;

	return;
}

int outaasm_size_of_working_space(){
	return sizeof(outAASM);
}

unsigned char* initAASM()
{
	int size = outaasm_size_of_working_space();
	unsigned char * working_space = (unsigned char *)malloc(sizeof(unsigned char)*size);
	memset(working_space, 0, size);

	outAASM* outaasm = (outAASM*)working_space;
	outaasm->totalEpoch = 0;

	return working_space;
}

int* getTotalEpoch(unsigned char *working_space)
{
	outAASM* outaasm = (outAASM*)working_space;

	return outaasm->score;
}


void rtEpoch(int* rawdata, float *fftdata, unsigned char *working_space)
{
	int indClear[SECONDS_OF_EPOCH];
	int index;
	const int rawTh = 500;

	outAASM* outaasm = (outAASM*)working_space;

	int *totalEpoch;
	int *epochValid;
	float *epoch;
	float *epoch2;

	totalEpoch = &outaasm->totalEpoch;
	epochValid = outaasm->epochValid;
	epoch = outaasm->epoch;
	epoch2 = outaasm->epoch2;

	index = outaasm->totalEpoch;

	rtRemoveDC(rawdata);
	rtClearRaw(rawdata, rawTh, indClear);
	rtEpochCal(fftdata, indClear, &epoch[index], &epoch2[index], &epochValid[index], totalEpoch);

	return;
}


void nrtAASM(unsigned char *working_space)
{
	const int window = 5;
	float *th;
	float *th2;
	float *epochFilted;
	float *epochFilted2;
	float *epoch2Filted;
	float *epoch2Filted2;

	outAASM* outaasm = (outAASM*)working_space;

	int totalEpoch;
	int *score;
	int *epochValid;
	float *epoch;
	float *epoch2;

	totalEpoch = outaasm->totalEpoch;
	epochValid = outaasm->epochValid;
	score = outaasm->score;
	epoch = outaasm->epoch;
	epoch2 = outaasm->epoch2;

	th = (float*)malloc(totalEpoch*sizeof(float));
	th2 = (float*)malloc(totalEpoch*sizeof(float));
	epochFilted = (float*)malloc(totalEpoch*sizeof(float));
	epochFilted2 = (float*)malloc(totalEpoch*sizeof(float));
	epoch2Filted = (float*)malloc(totalEpoch*sizeof(float));
	epoch2Filted2 = (float*)malloc(totalEpoch*sizeof(float));

	averagingFilter(epoch, totalEpoch, window, epochFilted);
	averagingFilter(epochFilted, totalEpoch, window, epochFilted2);
	thSleepWake(epochFilted2, totalEpoch, th, ORIGINAL_THRESHOLD);
	averagingFilter(epoch2, totalEpoch, window, epoch2Filted);
	averagingFilter(epoch2Filted, totalEpoch, window, epoch2Filted2);
	thSleepWake(epoch2Filted2, totalEpoch, th2, ORIGINAL_THRESHOLD2);
	judge(epochFilted2, th, epoch2Filted2, th2, epochValid, totalEpoch, score);

	free(th);
	free(th2);
	free(epochFilted);
	free(epochFilted2);
	free(epoch2Filted);
	free(epoch2Filted2);

	return;
}

static void combine_deal(int *bp_style_vector, int *bp_move_vector, int style_thirty_seconds_len, int raw_thirty_seconds_len, int *score)
{
	int index = 0;
	const int bp_style_num = 0;
	const int bp_move_num = 0;
	
	for(index=0; index<raw_thirty_seconds_len; index++)
	{
		if(index<=style_thirty_seconds_len && bp_style_vector[index]>bp_style_num && score[index]!= NaN)
		{
			score[index] = WAKE;
		}
		
		if(index<=style_thirty_seconds_len && bp_style_vector[index]==bp_style_num && bp_move_vector[index]==bp_move_num && score[index]== LIGHT)
		{
			score[index] = LIGHT;
		}
		
		if(index<=style_thirty_seconds_len && bp_style_vector[index]==bp_style_num && bp_move_vector[index]>bp_move_num && (score[index]== LIGHT || score[index]== DEEP))
		{
			score[index] = LIGHT;
		}
		
		if(index<=style_thirty_seconds_len && score[index]== NaN) 
		{
			score[index] = NaN;
		} 
	
	}
	
	return;
}

static void get_precise_body_position_style(int* bodyposition, int* bodymove, int *bodypositionVector, int *bodymoveVector,  int style_totalSecond)
{
	
	int style_thirty_seconds_num = style_totalSecond/(2*SECONDS_OF_EPOCH);
	int index = 0;
	int i = 0, j = 0;
	
	memset(bodypositionVector, 0, sizeof(int)*style_thirty_seconds_num);
	memset(bodymoveVector, 0, sizeof(int)*style_thirty_seconds_num);
	
	for(index=0; index<style_totalSecond; index++)
	{
		if(bodyposition[index]==0 || bodyposition[index]==7)
		{
			if(index==0)
			{
				bodyposition[index] = 3;
			}
			else
			{
				bodyposition[index] =bodyposition[index-1];
			}
			
		} 
		
		bodyposition[index] = bodyposition[index]; 
	}
	
	for(i=0; i<style_thirty_seconds_num; i++)
	{
		for(j=0; j<2*SECONDS_OF_EPOCH; j++)
		{
			if(bodyposition[i*2*SECONDS_OF_EPOCH+j]==5 || bodyposition[i*2*SECONDS_OF_EPOCH+j]==6)
			{
				bodypositionVector[i]+=1;
			}
			if(bodymove[i*2*SECONDS_OF_EPOCH+j]>0)
			{
				bodymoveVector[i]+=1;
			}
	    }	 
	}
	
}


void eeg_detect_algorithm(int *bodyposition, int *bodymove, int style_totalSecond, unsigned char *working_space)
{																				
	 outAASM* outaasm = (outAASM*)working_space;

	 int style_thirty_seconds_len = style_totalSecond/(2*SECONDS_OF_EPOCH);
	 int raw_thirty_seconds_len = outaasm->totalEpoch;
	 int *bodypositionVector;
	 int *bodymoveVector;
	 int *score;

	 score = outaasm->score;
	 
	 bodypositionVector = (int*)malloc(sizeof(int)*style_thirty_seconds_len);
	 bodymoveVector = (int*)malloc(sizeof(int)*style_thirty_seconds_len);
	 
	 get_precise_body_position_style(bodyposition, bodymove, bodypositionVector, bodymoveVector, style_totalSecond);
	 combine_deal(bodypositionVector, bodymoveVector, style_thirty_seconds_len, raw_thirty_seconds_len, score);
	 
	 return;
}
