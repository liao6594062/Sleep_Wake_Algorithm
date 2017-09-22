/*******************************************************************
Describes: 	main function for test to detection/score of wake, sleep and REM
   Author:	Jango
  Version: 	v1.1
     Date: 	2017-9-19
Copyright:	EEGSmart
********************************************************************/
#include "readcsv.h"
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


FILE* fp;

int main()
{
	int			m,n;
	int 		i;
	int 		numEpoch;
	int 		totalEpoch;
	int         style_totalSecond;
	int*		rawdata;
	int*		fftdata;
	int*        bodyposition;
	int*        bodymove;
	float*		fftdata2;

	char*		rawfile = "C:\\workspace\\data\\eegs\\rawdata2.csv";
	char*		fftfile = "C:\\workspace\\data\\eegs\\fftdata2.csv";
	char*       bpnfile = "C:\\workspace\\data\\eegs\\bodyposition2.csv";
	char*       bmvfile = "C:\\workspace\\data\\eegs\\bodymove2.csv";

	int *score;

	unsigned char* working_space = initAASM();

	n = get_col_of_raw(rawfile);
	m = get_row_of_raw(rawfile,n);
	rawdata = (int*)malloc(n*m*sizeof(int));
	get_raw_from_csv(rawfile,rawdata,m,n);

	n = get_col_of_raw(fftfile);
	m = get_row_of_raw(fftfile,n);
	fftdata = (int*)malloc(n*m*sizeof(int));
	get_raw_from_csv(fftfile,fftdata,m,n);

	fftdata2 = (float*)malloc(n*m*sizeof(float));

	for(i=0; i<m*n; i++)
		fftdata2[i] = fftdata[i];
		
	free(fftdata);

	totalEpoch = m/SECONDS_OF_EPOCH;
	numEpoch = 0;

	while(1)
	{
		rtEpoch(&rawdata[SECONDS_OF_EPOCH*NFFT*numEpoch], &fftdata2[SECONDS_OF_EPOCH*NFFT*numEpoch], working_space);

		numEpoch++;

		if(numEpoch==totalEpoch)
			break;
	}

	nrtAASM(working_space);

	n = get_col_of_raw(bpnfile);
	m = get_row_of_raw(bpnfile,n);
	bodyposition = (int*)malloc(n*m*sizeof(int));
	get_raw_from_csv(bpnfile,bodyposition,m,n);
	
	n = get_col_of_raw(bmvfile);
	m = get_row_of_raw(bmvfile,n);
	bodymove = (int*)malloc(n*m*sizeof(int));
	get_raw_from_csv(bmvfile,bodymove,m,n);
	style_totalSecond = m;
	
    eeg_detect_algorithm(bodyposition, bodymove, style_totalSecond, working_space);
    
    score = getTotalEpoch(working_space);

	fp = fopen("rtscore.txt", "w+");
	for(i=0; i<totalEpoch; i++)
	{
		fprintf(fp,"%d\n",score[i]);
	}
	fclose(fp);
	
	free(rawdata);
	free(bodyposition);
	free(bodymove);
	free(fftdata2);

	printf("\nBingo.");

	return 0;
}

