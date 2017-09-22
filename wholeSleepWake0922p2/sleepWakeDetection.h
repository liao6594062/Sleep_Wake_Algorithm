/*******************************************************************
Describes: 	header for detection/score of wake, sleep and REM
   Author:	Jango
  Version: 	v1.0
     Date: 	2017-9-19
Copyright:	EEGSmart
********************************************************************/
#ifndef __SLEEPWAKEDETECTION_H__
#define __SLEEPWAKEDETECTION_H__
#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <string.h>


#ifdef __cplusplus

extern "C" {

#endif
//void aasm(int* rawdata, float *fftdata, int totalSecond, char *score);
unsigned char* initAASM();
int* getTotalEpoch(unsigned char *working_space);
void rtEpoch(int* rawdata, float *fftdata, unsigned char *working_space);
/* @输入参数
 * 			rawdata: 数组，例如总长度为（256*3600秒*20小时），原始时域信号，采样率256Hz
 * 			fftdata: 数组，长度与rawdata相同，例如总长度为（256*3600秒*20小时），原始时域信号每秒256点FFT后结果的绝对值平方
            notFirst: 整数，初始值为0，之后每隔30s值加1 
 * *
 * @输出参数
 * 			epoch：数组，例如总长度为2400，区分是否睡觉的门限值变化曲线 
 *          epoch2：数组，总长度与epoch相同，例如总长度为2400，区分深睡和浅睡的门限值变化曲线
            epochValid：数组，总长度与epoch相同，判断每30s是否有效，1为有效，-1则为无效 
 * @return
 * 			无
 */
 
void nrtAASM(unsigned char *working_space);
/* @输入参数
 * 			epoch：数组，例如总长度为2400，区分是否睡觉的门限值变化曲线 
 *          epoch2：数组，总长度与epoch相同，例如总长度为2400，区分深睡和浅睡的门限值变化曲线
            totalEpoch: 整数， 表示总共有多少个30s
			epochValid：数组，总长度与epoch相同，判断每30s是否有效，1为有效，-1则为无效
 * *
 * @输出参数
 * 			score：数组，总长度与epoch长度相同，脑电分期结果 
 *         
 * @return
 * 			无
 */

void eeg_detect_algorithm(int *bodyposition, int *bodymove, int style_totalSecond, unsigned char *working_space);

/* @输入参数
 * 			bodyposition：数组，记录的整夜体位变化值 
 *          bodymove：数组，记录的整夜体动变化值 
            style_totalSecond: 整数， bodyposition对应的元素个数 
			raw_totalSecond：数组，总共记录脑电的秒数 
 * *
 * @输出参数
 * 			score：数组，总长度与epoch长度相同，脑电与体位融合后的脑电分期结果 
 *         
 * @return
 * 			无
 */
#ifdef __cplusplus

}

#endif /* end of __cplusplus */
 
#endif
