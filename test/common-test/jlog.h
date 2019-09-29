#ifndef __JLOG_H
#define __JLOG_H

////////////////////////
// Logging Intrinsics //
////////////////////////
#include <stdint.h>

#define RESET "\033[0m"
#define BLACK "\033[30m"              /* Black */
#define RED "\033[31m"                /* Red */
#define GREEN "\033[32m"              /* Green */
#define YELLOW "\033[33m"             /* Yellow */
#define BLUE "\033[34m"               /* Blue */
#define MAGENTA "\033[35m"            /* Magenta */
#define CYAN "\033[36m"               /* Cyan */
#define WHITE "\033[37m"              /* White */
#define BOLDBLACK "\033[1m\033[30m"   /* Bold Black */
#define BOLDRED "\033[1m\033[31m"     /* Bold Red */
#define BOLDGREEN "\033[1m\033[32m"   /* Bold Green */
#define BOLDYELLOW "\033[1m\033[33m"  /* Bold Yellow */
#define BOLDBLUE "\033[1m\033[34m"    /* Bold Blue */
#define BOLDMAGENTA "\033[1m\033[35m" /* Bold Magenta */
#define BOLDCYAN "\033[1m\033[36m"    /* Bold Cyan */
#define BOLDWHITE "\033[1m\033[37m"   /* Bold White */

#define MSG(...)         \
  do {                   \
    printf(__VA_ARGS__); \
  } while (0)

#define MSGF(enable, ...)                                               \
  do {                                                                  \
    if (_verbosity & enable) {                                          \
      printf(MAGENTA "[" BOLDMAGENTA "%s" MAGENTA "] " RESET, #enable); \
      printf(__VA_ARGS__);                                              \
    }                                                                   \
  } while (0)

#define MSGO(...)               \
  do {                          \
    static uint8_t printed = 0; \
    if (!printed) {             \
      printf(__VA_ARGS__);      \
      printed = 1;              \
    }                           \
  } while (0)

#define DEBUG_MSG(...)                                                        \
  do {                                                                        \
    printf(BLUE "(%s:%d - %s)" RESET "\n", __FILE__, __LINE__, __FUNCTION__); \
    printf(__VA_ARGS__);                                                      \
  } while (0)

#define WARNING_MSG(...)                                                  \
  do {                                                                    \
    printf(YELLOW "WARNING! (%s:%d - %s)" RESET "\n", __FILE__, __LINE__, \
           __FUNCTION__);                                                 \
    printf(__VA_ARGS__);                                                  \
  } while (0)

#define ERROR_MSG(...)                                                        \
  do {                                                                        \
    fprintf(stderr, RED "ERROR! (%s:%d - %s)" RESET "\n", __FILE__, __LINE__, \
            __FUNCTION__);                                                    \
    fprintf(stderr, __VA_ARGS__);                                             \
  } while (0)

/////////////////////////
// CSV Plotting Tools ///
/////////////////////////

#include <time.h>

// TIC TOC only counts process time, it will not count sleeps or io waits
#define CLOCK(ID_) static clock_t tic##ID_;

#define TIC(ID_) tic##ID_ = clock();

#define TOC(ID_)                                                               \
  do {                                                                         \
    clock_t toc##ID_ = clock();                                                \
    if (tic##ID_ > toc##ID_) {                                                 \
      MSG(BLUE "%s" RESET " - " RED "Timer %s Overflowed" RESET "\n",          \
          __FUNCTION__, #ID_);                                                 \
    } else {                                                                   \
      MSG(BLUE "%s" RESET " - " GREEN "Timer %s: %g s" RESET "\n",             \
          __FUNCTION__, #ID_, (toc##ID_ - tic##ID_) / (double)CLOCKS_PER_SEC); \
    }                                                                          \
  } while (0)

#define TOC_CSV_XY(ID_, X_)                                                    \
  do {                                                                         \
    clock_t toc##ID_ = clock();                                                \
    if (tic##ID_ > toc##ID_) {                                                 \
      ERROR_MSG("SeNtInAl,error,%s,%s,%d,overflow\n", __FUNCTION__, #ID_, X_); \
    } else {                                                                   \
      MSG("SeNtInAl,xy,%s,%s,%d,%g\n", __FUNCTION__, #ID_, X_,                 \
          (toc##ID_ - tic##ID_) / (double)CLOCKS_PER_SEC);                     \
    }                                                                          \
  } while (0)

// will stack ID's with same XLABEL_ on top
#define TOC_CSV_BAR(ID_, XLABEL_)                                         \
  do {                                                                    \
    clock_t toc##ID_ = clock();                                           \
    if (tic##ID_ > toc##ID_) {                                            \
      ERROR_MSG("SeNtInAl,error,%s,%s,%s,overflow\n", __FUNCTION__, #ID_, \
                XLABEL_);                                                 \
    } else {                                                              \
      MSG("SeNtInAl,bar,%s,%s,%s,%g\n", __FUNCTION__, #ID_, XLABEL_,      \
          (toc##ID_ - tic##ID_) / (double)CLOCKS_PER_SEC);                \
    }                                                                     \
  } while (0)

#define TOC_CSV_GROUPED_BAR(ID_, XGROUP_)                                 \
  do {                                                                    \
    clock_t toc##ID_ = clock();                                           \
    if (tic##ID_ > toc##ID_) {                                            \
      ERROR_MSG("SeNtInAl,error,%s,%s,%s,overflow\n", __FUNCTION__, #ID_, \
                XGROUP_);                                                 \
    } else {                                                              \
      MSG("SeNtInAl,bar,%s,%s,%s,%g\n", __FUNCTION__, #ID_, XGROUP_,      \
          (toc##ID_ - tic##ID_) / (double)CLOCKS_PER_SEC);                \
    }                                                                     \
  } while (0)

#define TOC_CSV_BOX(ID_)                                                     \
  do {                                                                       \
    clock_t toc##ID_ = clock();                                              \
    if (tic##ID_ > toc##ID_) {                                               \
      ERROR_MSG("SeNtInAl,error,%s,%s,no-x,overflow\n", __FUNCTION__, #ID_); \
    } else {                                                                 \
      MSG("SeNtInAl,box,%s,%s,no-x,%g\n", __FUNCTION__, #ID_,                \
          (toc##ID_ - tic##ID_) / (double)CLOCKS_PER_SEC);                   \
    }                                                                        \
  } while (0)

// Wall clock measures absolute time

#include <chrono>
#include <ctime>
#include <ratio>
using namespace std::chrono;

#define WALL_CLOCK(ID_) static time_point<high_resolution_clock> tic##ID_;

#define WALL_TIC(ID_) tic##ID_ = high_resolution_clock::now();

#define WALL_TOC(ID_)                                                          \
  do {                                                                         \
    time_point<high_resolution_clock> toc##ID_ = high_resolution_clock::now(); \
    if (tic##ID_ > toc##ID_) {                                                 \
      MSG(BLUE "%s" RESET " - " RED "Timer %s Overflowed" RESET "\n",          \
          __FUNCTION__, #ID_);                                                 \
    } else {                                                                   \
      MSG(BLUE "%s" RESET " - " GREEN "Timer %s: %g s" RESET "\n",             \
          __FUNCTION__, #ID_,                                                  \
          std::chrono::duration_cast<std::chrono::microseconds>(toc##ID_ -     \
                                                                tic##ID_)      \
                  .count() /                                                   \
              1000000.0);                                                      \
    }                                                                          \
  } while (0)

#define WALL_TOC_CSV_XY(ID_, X_)                                               \
  do {                                                                         \
    time_point<high_resolution_clock> toc##ID_ = high_resolution_clock::now(); \
    if (tic##ID_ > toc##ID_) {                                                 \
      ERROR_MSG("SeNtInAl,error,%s,%s,%d,overflow\n", __FUNCTION__, #ID_, X_); \
    } else {                                                                   \
      MSG("SeNtInAl,xy,%s,%s,%d,%g\n", __FUNCTION__, #ID_, X_,                 \
          std::chrono::duration_cast<std::chrono::microseconds>(toc##ID_ -     \
                                                                tic##ID_)      \
                  .count() /                                                   \
              1000000.0);                                                      \
    }                                                                          \
  } while (0)

#define WALL_TOC_CSV_XY_NORMALIZED(ID_, X_, DIV_)                              \
  do {                                                                         \
    time_point<high_resolution_clock> toc##ID_ = high_resolution_clock::now(); \
    if (tic##ID_ > toc##ID_) {                                                 \
      ERROR_MSG("SeNtInAl,error,%s,%s,%d,overflow\n", __FUNCTION__, #ID_, X_); \
    } else {                                                                   \
      MSG("SeNtInAl,xy,%s,%s_normalized,%d,%g\n", __FUNCTION__, #ID_, X_,      \
          std::chrono::duration_cast<std::chrono::microseconds>(toc##ID_ -     \
                                                                tic##ID_)      \
                  .count() /                                                   \
              (1000000.0 * DIV_));                                             \
    }                                                                          \
  } while (0)

// will stack ID's with same XLABEL on top
#define WALL_TOC_CSV_BAR(ID_, XLABEL_)                                         \
  do {                                                                         \
    time_point<high_resolution_clock> toc##ID_ = high_resolution_clock::now(); \
    if (tic##ID_ > toc##ID_) {                                                 \
      ERROR_MSG("SeNtInAl,error,%s,%s,%s,overflow\n", __FUNCTION__, #ID_,      \
                XLABEL_);                                                      \
    } else {                                                                   \
      MSG("SeNtInAl,bar,%s,%s,%s,%g\n", __FUNCTION__, #ID_, XLABEL_,           \
          std::chrono::duration_cast<std::chrono::microseconds>(toc##ID_ -     \
                                                                tic##ID_)      \
                  .count() /                                                   \
              1000000.0);                                                      \
    }                                                                          \
  } while (0)

#define WALL_TOC_CSV_GROUPED_BAR(ID_, XGROUP_)                                 \
  do {                                                                         \
    time_point<high_resolution_clock> toc##ID_ = high_resolution_clock::now(); \
    if (tic##ID_ > toc##ID_) {                                                 \
      ERROR_MSG("SeNtInAl,error,%s,%s,%s,overflow\n", __FUNCTION__, #ID_,      \
                XGROUP_);                                                      \
    } else {                                                                   \
      MSG("SeNtInAl,grouped_bar,%s,%s,%s,%g\n", __FUNCTION__, #ID_, XGROUP_,   \
          std::chrono::duration_cast<std::chrono::microseconds>(toc##ID_ -     \
                                                                tic##ID_)      \
                  .count() /                                                   \
              1000000.0);                                                      \
    }                                                                          \
  } while (0)

#define WALL_TOC_CSV_BOX(ID_)                                                \
  do {                                                                       \
    static time_point<high_resolution_clock> toc##ID_ =                      \
        high_resolution_clock::now();                                        \
    if (tic##ID_ > toc##ID_) {                                               \
      ERROR_MSG("SeNtInAl,error,%s,%s,no-x,overflow\n", __FUNCTION__, #ID_); \
    } else {                                                                 \
      MSG("SeNtInAl,box,%s,%s,no-x,%g\n", __FUNCTION__, #ID_,                \
          std::chrono::duration_cast<std::chrono::microseconds>(toc##ID_ -   \
                                                                tic##ID_)    \
                  .count() /                                                 \
              1000000.0);                                                    \
    }                                                                        \
  } while (0)

// Counter
#define COUNTER(ID_) static uint32_t count##ID_ = 0;

#define COUNTER_RESET(ID_, VAL_) \
  do {                           \
    count##ID_ = VAL_;           \
  } while (0)

#define COUNTER_INC(ID_, INC_BY_)                                         \
  do {                                                                    \
    static uint32_t tmp = count##ID_ + INC_BY_;                           \
    if (count##ID_ > tmp) {                                               \
      printf("SeNtInAl,error,%s,%s,no-x,overflow\n", __FUNCTION__, #ID_); \
    }                                                                     \
    count##ID_ = tmp;                                                     \
  } while (0)

#define COUNTER_PRINT(ID_, X_)                                          \
  do {                                                                  \
    MSG(BLUE "%s" RESET " - " GREEN "Counter %s: \t%d,\t%d" RESET "\n", \
        __FUNCTION__, #ID_, X_, count##ID_);                            \
  } while (0)

#define COUNTER_CSV_XY(ID_, X_)                                              \
  do {                                                                       \
    printf("SeNtInAl,xy,%s,%s,%d,%d\n", __FUNCTION__, #ID_, X_, count##ID_); \
  } while (0)

#define COUNTER_CSV_BOX(ID_)                                                \
  do {                                                                      \
    printf("SeNtInAl,box,%s,%s,no-x,%d\n", __FUNCTION__, #ID_, count##ID_); \
  } while (0)

// Histogram
#define HISTOCSV(ID_, VAL_)                                                 \
  do {                                                                      \
    printf("SeNtInAl,histogram,%s,%s,no-x,%d\n", __FUNCTION__, #ID_, VAL_); \
  } while (0)

#endif
