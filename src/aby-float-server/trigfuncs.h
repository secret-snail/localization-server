#include <jlog.h>
#include <stdio.h>
#include <iostream>
#include <string>

#include "abycore/aby/abyparty.h"
#include "abycore/circuit/arithmeticcircuits.h"
#include "abycore/circuit/booleancircuits.h"
#include "abycore/circuit/circuit.h"
#include "abycore/sharing/sharing.h"

using namespace std;

share* BuildSinCircuit(share* theta, BooleanCircuit* c);
share* BuildCosCircuit(share* theta, BooleanCircuit* c);
