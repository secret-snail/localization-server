#include "abycore/aby/abyparty.h"
#include "abycore/circuit/arithmeticcircuits.h"
#include "abycore/circuit/booleancircuits.h"
#include "abycore/circuit/circuit.h"
#include "abycore/sharing/sharing.h"

void BuildInvertCircuit(
    share** _s_in[], int m, int n, share** s_res[], BooleanCircuit* c,
    ABYParty* party, e_role role, std::function<void()> toRawShares = []() {},
    std::function<void()> toShareObjects = []() {});
