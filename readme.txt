Curvep is a heuristic algorithm for dose-response curve processing.
It assumes monotonicity, but does not fit curve to observed responses (as isotonic regression).
Instead, it minimizes the number of corrections that need to be made to restore monotonicity.
---------------

Curvep accepts *.hts and *.htsx formats, see "input_demo*"-files
Example runs: 
	curvep input_demo1 -HTSX=uM -RNG=-100 -THR=15 -MXDV=5
	curvep input_demo2 -RNG=100 -THR=15 -MXDV=5

Key parameters are 
	baseline noise threshold (-THR=) - points with smaller response are suppressed
	maximum allowed deviation from the  monotonicity (-MXDV=) - affects adjacent points
	range of response values (-RNG=) - can be negative (e.g., inhibition) or positive (e.g., agonism)

----------------
Run Curvep without parameters in command line to print out more detailed help:

#CurveP V5.12 - qHTS curve processing and noise-reduction
Sedykh et al, Environ Health Perspect (2011), 119, p.364

Usage:     curvep filename [flags]
Allowed inputs are tab-delimited text files of following formats:
.hts file: 1st line - #samples and #test.concentrations (e.g., 100  8)
2nd line - column headers, 'ID(optional) ExclPoints 0.6nM_% 3uM_% 3mM_%'..
Following lines: ids(optional), mask, responses at each test concentration
.htsx file: 1st line - same as in .hts file
2nd line - column headers, 'ID(optional) ExclPoints Conc1 Conc2..Resp1 Resp2..'
Following lines: ids(optional), mask, test conc-s, responses
Mask is either in '{1.2, 1.4}' format or a binary string
('00010010' or '0 0 0 1 0 0 1 0') were 1's denote invalidated points
NB: Mask in output file is extended to include corrected points
NB: Conc.-based output values are in log10(M).

Optional Flags:
'-HTSX=' - input units: 'nM' 'uM'<def.> 'M' 'lgM' etc.
'-OUT=' output file name; '-LOG' detailed log; '-ALT' altern.fingerprint
'-MXDV=' max.deviation <def. 5.0>; '-THR=' baseline threshold <def.15.0>
'-RNG=' max.response <def. -100> (can be positive)
'-CRO=' carryover threshold <def. 80% of max.response>, not used if 0
'-BLFX' baseline shift correction mode, off by def.
'-USHAPE=' min.#points for u-shape <def.4>
'-BSHIFT=' min.#flat points to detect baseline-shift <def.3>
'-BYHI','-BYLO'<def> favors corrections based on hi or lo conc-s
'-XPLAX' allow extrapolation beyond test conc. boundaries
'-NA=' dummy value for missing data <def. -999>
'-TAG=' optional marker to add to column headers of output file
