#include <stdio.h>
#include <vector>

#include <abycore/aby/abyparty.h>
//#include <../ENCRYPTO_utils/src/ENCRYPTO_utils/timer.h> // hack {:

#define BUILD_TIMING 0
#define EXEC_TIMING 0

void collectTiming();
void collectCommunication();
void PrintSumTimings(int x);
void ClearSumTimings();
void PrintSumCommunication();
