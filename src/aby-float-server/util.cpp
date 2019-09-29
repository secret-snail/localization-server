#include <stdio.h>
#include <iostream>
#include <vector>

#include <../ENCRYPTO_utils/src/ENCRYPTO_utils/timer.h>  // hack {:
#include <abycore/aby/abyparty.h>
#include <jlog.h>
#include "util.h"

#define BUILD_TIMING 0
#define EXEC_TIMING 0

aby_timings sumTimes[P_LAST - P_FIRST + 1];
aby_comm sumSend[P_LAST - P_FIRST + 1];
aby_comm sumRecv[P_LAST - P_FIRST + 1];

void collectTiming() {
  for (int phase = P_FIRST; phase <= P_LAST; phase++) {
    sumTimes[phase].timing += m_tTimes[phase].timing;
  }
}

void collectCommunication() {
  for (int phase = P_FIRST; phase <= P_LAST; phase++) {
    sumSend[phase].totalcomm += m_tSend[phase].totalcomm;
    sumRecv[phase].totalcomm += m_tRecv[phase].totalcomm;
  }
}

void PrintSumTimings(int x) {
  std::string unit = " ms";
  std::cout << "Sum Total Timings: " << std::endl;
  MSG("SeNtInAl,xy,%s,%s,%d,%g\n", __FUNCTION__, "gngd_total", x,
      sumTimes[P_TOTAL].timing);
  MSG("SeNtInAl,xy,%s,%s,%d,%g\n", __FUNCTION__, "gngd_init", x,
      sumTimes[P_INIT].timing);
  MSG("SeNtInAl,xy,%s,%s,%d,%g\n", __FUNCTION__, "gngd_circuit", x,
      sumTimes[P_CIRCUIT].timing);
  MSG("SeNtInAl,xy,%s,%s,%d,%g\n", __FUNCTION__, "gngd_network", x,
      sumTimes[P_NETWORK].timing);
  MSG("SeNtInAl,xy,%s,%s,%d,%g\n", __FUNCTION__, "gngd_baseot", x,
      sumTimes[P_BASE_OT].timing);
  MSG("SeNtInAl,xy,%s,%s,%d,%g\n", __FUNCTION__, "gngd_setup", x,
      sumTimes[P_SETUP].timing);
  MSG("SeNtInAl,xy,%s,%s,%d,%g\n", __FUNCTION__, "gngd_otext", x,
      sumTimes[P_OT_EXT].timing);
  MSG("SeNtInAl,xy,%s,%s,%d,%g\n", __FUNCTION__, "gngd_garble", x,
      sumTimes[P_GARBLE].timing);
  MSG("SeNtInAl,xy,%s,%s,%d,%g\n", __FUNCTION__, "gngd_online", x,
      sumTimes[P_ONLINE].timing);
  for (int phase = P_FIRST; phase <= P_LAST; phase++) {
    sumTimes[phase].timing = 0;
  }
}

void ClearSumTimings() {
  for (int phase = P_FIRST; phase <= P_LAST; phase++) {
    sumTimes[phase].timing = 0;
  }
}

void PrintSumCommunication() {
  std::string unit = " bytes";
  std::cout << "Communication: " << std::endl;
  std::cout << "Total Sent / Rcv\t" << m_tSend[P_TOTAL].totalcomm << " " << unit
            << " / " << m_tRecv[P_TOTAL].totalcomm << unit << std::endl;
  std::cout << "BaseOTs Sent / Rcv\t" << m_tSend[P_BASE_OT].totalcomm << " "
            << unit << " / " << m_tRecv[P_BASE_OT].totalcomm << unit
            << std::endl;
  std::cout << "Setup Sent / Rcv\t" << m_tSend[P_SETUP].totalcomm << " " << unit
            << " / " << m_tRecv[P_SETUP].totalcomm << unit << std::endl;
  std::cout << "OTExtension Sent / Rcv\t" << m_tSend[P_OT_EXT].totalcomm << " "
            << unit << " / " << m_tRecv[P_OT_EXT].totalcomm << unit
            << std::endl;
  std::cout << "Garbling Sent / Rcv\t" << m_tSend[P_GARBLE].totalcomm << " "
            << unit << " / " << m_tRecv[P_GARBLE].totalcomm << unit
            << std::endl;
  std::cout << "Online Sent / Rcv\t" << m_tSend[P_ONLINE].totalcomm << " "
            << unit << " / " << m_tRecv[P_ONLINE].totalcomm << unit
            << std::endl;
}
