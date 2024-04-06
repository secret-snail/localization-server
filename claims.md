# Main Results and Claims

The results and accompanying claims are enumerated below.

1. Figure 2 compares the runtime of localization using two MPC libraries, EMP
and ABY. Claim: EMP is better suited to localization than ABY.
2. Figure 3 measures the runtime of arithmetic operations on various data
representations. Claim: 64 bit fixed point multiplication is slower than 32 bit
floating.
3. Figure 4 counts the number of arithmetic operations performed during
localization. Claim: Multiplication is the dominant operation.
4. Figure 6 compares the time to localize using data oblivious (DO) vs. single
iteration localization (SIL). Claim: SIL localizes faster than DO with the
Levenburg Marquardt (LM) optimization algorithm being faster than Gauss Newton
(GN).
5. Figure 7 measures runtime and network IO for different localization
configurations at large input sizes. Claim: LM better scales to large input
sizes than GN, making LM the better approach of the two.

## Main Result 1: Figure 2

The EMP MPC library is better suited to localization than the ABY MPC library.
Open the file `plots/emp_vs_aby.pdf` and compare it to Figure 2 in the paper.
The EMP LM and EMP GN bars should be smaller than the ABY bars.

## Main Result 2: Figure 3

64 bit fixed point multiplication is slower than 32 bit floating point. Open
`plots/emp_float_vs_fixed_benchmark_add_mul.pdf` and verify the size of these
two bars noting the y axis for multiplication is is on the right.

## Main Result 3: Figure 4

Multiplication is the dominant operation in LM-based localization.
Open `plots/emp_arith_ops.pdf` and verify the green multiplication bars are
larger than all other bars at any of the measured number of features.

## Main Result 4: Figure 6

Single iteration localization (SIL) localizes faster than the data oblivious
(DO) adaptation with LM being the faster optimization algorithm. Open
`plots/loopleak_vs_dataobl.pdf` and verify the DO lines are (much) higher than
the SIL lines. Also check the purple LM SIL line is the lowest in the figure
for all the numbers of features.

## Main Result 5: Figure 7

LM better scales to large input size i.e. large numbers of input features.
Open `plots/emp_float_runtime_long.pdf` and verify the GN lines are higher than
their respective LM lines for the same latency. Next open `plots/netio.pdf`
and verify the LM lines are below their respective GN lines.