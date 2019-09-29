#pragma once
#include <stdio.h>
#include <iostream>
#include <string>

#include "abycore/aby/abyparty.h"
#include "abycore/circuit/arithmeticcircuits.h"
#include "abycore/circuit/booleancircuits.h"
#include "abycore/circuit/circuit.h"
#include "abycore/sharing/sharing.h"

void BuildFabsCircuit(share* a, BooleanCircuit* c);

void BuildSignCircuit(share* a, share* b, BooleanCircuit* c);

share* BuildPythagCircuit(share* a, share* b, BooleanCircuit* c,
                          uint32_t bitlen);

// uint32_t arguments must be from PutSharedOUTGate()
// then calling get_clear_value() and circuit->Reset();
// This function builds/executes multiple circuits.
void RunSvd(uint32_t** a, int nRows, int nCols, uint32_t* w, uint32_t** v,
            BooleanCircuit* c, ABYParty* p, e_role role);

// Wrapper function around RunSvd which
// takes shares instead of secret shares.
// Creates dummy circuits to build secret shares.
void BuildAndRunSvd(
    share*** s_a, int nRows, int nCols, share** s_w, share*** s_v,
    BooleanCircuit* c, ABYParty* party, e_role role,
    std::function<void()> toRawShares = []() {},
    std::function<void()> toShareObjects = []() {});

uint32_t test_fabs_circuit(e_role role, const std::string& address,
                           uint16_t port, seclvl seclvl, uint32_t nthreads,
                           e_mt_gen_alg mt_alg, e_sharing sharing, float a);

uint32_t test_fmax_circuit(e_role role, const std::string& address,
                           uint16_t port, seclvl seclvl, uint32_t nthreads,
                           e_mt_gen_alg mt_alg, e_sharing sharing, float a,
                           float b);

uint32_t test_negative_circuit(e_role role, const std::string& address,
                               uint16_t port, seclvl seclvl, uint32_t nthreads,
                               e_mt_gen_alg mt_alg, e_sharing sharing, float a);
