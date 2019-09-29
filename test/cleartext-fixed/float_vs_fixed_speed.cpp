#include <jlog.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string>

#include <printutil.h>
#include <util.h>

#include <fixed_point_emp.h>
#include "emp-sh2pc/emp-sh2pc.h"
#include "emp-tool/emp-tool.h"

using namespace std;

#define TEST_SZ 1000 * 10

#ifdef TEST_ADD
#define TEST_OP +
#define START()        \
  WALL_CLOCK(EMP_ADD); \
  WALL_TIC(EMP_ADD);
#define STOP(LOGSTR_) WALL_TOC_CSV_GROUPED_BAR(EMP_ADD, LOGSTR_);
#endif

#ifdef TEST_MUL
#define TEST_OP *
#define START()        \
  WALL_CLOCK(EMP_MUL); \
  WALL_TIC(EMP_MUL);
#define STOP(LOGSTR_) WALL_TOC_CSV_GROUPED_BAR(EMP_MUL, LOGSTR_);
#endif

int main(int argc, char** argv) {
  if (argc != 2) {
    MSG("Usage: %s <party number>\n", argv[0]);
    return 1;
  }
  int party = atoi(argv[1]) + 1;
  int port = 8080;
  cout << "party: " << party << " port: " << port << endl;
  NetIO* io = new NetIO(party == ALICE ? nullptr : "127.0.0.1", port);

  {
    MSG("\nBenchmarking EMP float\n");

    setup_semi_honest(io, party);

    Float* s_testa =
        static_cast<Float*>(operator new[](TEST_SZ * sizeof(Float)));
    Float* s_testb =
        static_cast<Float*>(operator new[](TEST_SZ * sizeof(Float)));
    Float* s_testc =
        static_cast<Float*>(operator new[](TEST_SZ * sizeof(Float)));
    for (int i = 0; i < TEST_SZ; i++) {
      s_testa[i] = Float(1.0f, ALICE);
      s_testb[i] = Float(2.0f, ALICE);
    }

    START();
    auto start = clock_start();
    for (int i = 0; i < TEST_SZ; i++) {
      s_testc[i] = s_testa[i] TEST_OP s_testb[i];
    }
    double interval = time_from(start);
    STOP("32 bit float");

    uint64_t numand = CircuitExecution::circ_exec->num_and();
    cout << "number of and gates: " << numand << endl;
    cout << "garbling speed : " << numand / interval
         << " million gate per second\n";

    // for(int i=0; i<TEST_SZ; i++) {
    //     cout << s_testc[i].reveal<double>() << " ";
    // }
    // cout << endl;

    delete[] s_testa;
    delete[] s_testb;
    delete[] s_testc;
    cout << "deleted\n";
    finalize_semi_honest();
  }

  {
    MSG("\nBenchmarking EMP Integer (32 bit)\n");

    setup_semi_honest(io, party);

    Integer s_testa[TEST_SZ];
    Integer s_testb[TEST_SZ];
    Integer s_testc[TEST_SZ];
    for (int i = 0; i < TEST_SZ; i++) {
      s_testa[i] = Integer(32, 1, ALICE);
      s_testb[i] = Integer(32, 2, ALICE);
    }

    START();
    auto start = clock_start();
    for (int i = 0; i < TEST_SZ; i++) {
      s_testc[i] = s_testa[i] TEST_OP s_testb[i];
    }
    double interval = time_from(start);
    STOP("32 bit int");

    uint64_t numand = CircuitExecution::circ_exec->num_and();
    cout << "number of and gates: " << numand << endl;
    cout << "garbling speed : " << numand / interval
         << " million gate per second\n";

    // for(int i=0; i<TEST_SZ; i++) {
    //     cout << s_testc[i].reveal<int>() << " ";
    // }
    // cout << endl;

    finalize_semi_honest();
  }

  {
    MSG("\nBenchmarking EMP Integer (64 bit)\n");

    setup_semi_honest(io, party);

    Integer s_testa[TEST_SZ];
    Integer s_testb[TEST_SZ];
    Integer s_testc[TEST_SZ];
    for (int i = 0; i < TEST_SZ; i++) {
      s_testa[i] = Integer(64, 1, ALICE);
      s_testb[i] = Integer(64, 2, ALICE);
    }

    START();
    auto start = clock_start();
    for (int i = 0; i < TEST_SZ; i++) {
      s_testc[i] = s_testa[i] TEST_OP s_testb[i];
    }
    double interval = time_from(start);
    STOP("64 bit int");

    uint64_t numand = CircuitExecution::circ_exec->num_and();
    cout << "number of and gates: " << numand << endl;
    cout << "garbling speed : " << numand / interval
         << " million gate per second\n";

    // for(int i=0; i<TEST_SZ; i++) {
    //     cout << s_testc[i].reveal<int>() << " ";
    // }
    // cout << endl;

    finalize_semi_honest();
  }

  {
    MSG("\nBenchmarking EMP fixed point (32 bits)\n");

    setup_semi_honest(io, party);

    fixed_point_emp<16, 16> s_testa[TEST_SZ];
    fixed_point_emp<16, 16> s_testb[TEST_SZ];
    fixed_point_emp<16, 16> s_testc[TEST_SZ];
    for (int i = 0; i < TEST_SZ; i++) {
      s_testa[i] = 1;
      s_testb[i] = 2;
    }

    START()
    auto start = clock_start();
    for (int i = 0; i < TEST_SZ; i++) {
      s_testc[i] = s_testa[i] TEST_OP s_testb[i];
    }
    double interval = time_from(start);
    STOP("32 bit fixed");

    uint64_t numand = CircuitExecution::circ_exec->num_and();
    cout << "number of and gates: " << numand << endl;
    cout << "garbling speed : " << numand / interval
         << " million gate per second\n";

    // for(int i=0; i<TEST_SZ; i++) {
    //     cout << s_testc[i].reveal<int>() << " ";
    // }
    // cout << endl;

    finalize_semi_honest();
  }

  {
    MSG("\nBenchmarking EMP fixed point (64 bits)\n");

    setup_semi_honest(io, party);

    fixed_point_emp<32, 32> s_testa[TEST_SZ];
    fixed_point_emp<32, 32> s_testb[TEST_SZ];
    fixed_point_emp<32, 32> s_testc[TEST_SZ];
    for (int i = 0; i < TEST_SZ; i++) {
      s_testa[i] = 1;
      s_testb[i] = 2;
    }

    START();
    auto start = clock_start();
    for (int i = 0; i < TEST_SZ; i++) {
      s_testc[i] = s_testa[i] TEST_OP s_testb[i];
    }
    double interval = time_from(start);
    STOP("64 bit fixed");

    uint64_t numand = CircuitExecution::circ_exec->num_and();
    cout << "number of and gates: " << numand << endl;
    cout << "garbling speed : " << numand / interval
         << " million gate per second\n";

    // for(int i=0; i<TEST_SZ; i++) {
    //     cout << s_testc[i].reveal<int>() << " ";
    // }
    // cout << endl;

    finalize_semi_honest();
  }
  delete io;
  return 0;
}
