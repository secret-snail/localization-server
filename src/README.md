# Localization Source

This directory contains source code of the localization libraries for each MPC
protocol, ABY and EMP, and some code shared between them.

Each library implements all sub-operations for localization (e.g. point
projection or singular value decomposition) in a separate source and header
files. The main optimization loops for Gauss Newton or Levenburg-Marquardt
algorithms are located in the `gaussnewtonlocalization.cpp` and
`lmlocalization.cpp` files (for both ABY and EMP). The EMP library additionally
contains `server-lm.cpp`, a simple network server for the Levenburg-Marquardt
algorithm used by the [robot](https://github.com/secret-snail/snail) described
in section 7.4 of the [paper](https://arxiv.org/pdf/2403.14916.pdf).

Each library, ABY or EMP, can be included independently in other CMake projects
using `target_link_libraries(... EmpFloatLocalization)`. See the `test`
directory for examples, which also demonstrate how to invoke localization. If
one just wants to invoke EMP LM-based localization from over the network, see
the server `emp-float-server/server-lm.cpp`, the client
`test/emp-float/client.cpp` and an
[example demo](https://github.com/secret-snail/snail/blob/main/src/emp_client.cpp).