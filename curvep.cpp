/* curvep.cpp Nov 1 2011 - 
Console application (Linux, Windows). 

Based on dataset-class
TODO:
detect baseline for the baselineshift better??

detect outlier by a consensus of trials? (to see how many times each point get removed, depending on the starting position)
handle carryovers with decreasing/constant signal in a similar way as with increasing? (i.e. detect minimal conc., and erase before that, not the entire curve)???


DONE: 

(history list of recent changes)
5.12	April 28 2014	fix for the rescue of 5.10 (it can overlook small spikes if overall st.dev is small)
5.11	April 26 2014	fix help for new USHAPE default
5.10	April 26 2014	fix for overritten WARNING flags
						rescue fix for flat curves with low variance and significant signal
						default for IGNORED_N_USHAPE readjusted to 4
						
5.03	April  9 2014	fixed wrong indexation in Baseline shift part
5.02	March  7 2014	wAUC calculation fixed
5.01	March  6 2014	Fixed wrong handling of -HTSX= and -XPLAX options

5.00
		Feb 25 2014		fix with remote folder execution
						expanded output, -XTINF removed and made as default output

						added support for log10 input and explicit base: nM, uM <def.>, mM, M

4.11-4.12
		Sept 9 2013		fix in baseline detection for rising carryover or baseline shift

4.10	August 12 2013	fix in carryover handling, added "INVERSE" warning if curve monotonicity contradicts the signal scale

4.01-4.03
		August 11 2013	small fix in baseline shift handling to suppress signals below baseline
						added restrictions on u-shapes (by total #spikes, #corrections, #support for retained part of u-shape)

4.00	August 8 2013	u-shape analysis reorganized to find pivot region (where monotonicity changes), 
						-MONINT key is removed, instead -USHAPE= is used for both and set to 3 as default,
						-BYHI and -BYLO(new def.) flags now control, which set of corrections to prefer if several of equal size (from high or low concentrations)
						added to output #corrected points as a last column
						added blip handling at an early stage (this can also handle weak carry-over cases)

3.43-3.44
		August 6 2013	skipped outlier detection if curve is constant (xRange == 0)
						fixed 'set TrialBest(0, nCols);' in outlier detection
						modified U_SHAPE analysis, checks it only if (fabs(xWrk) > alwdDeviation) (was '> 0' in V3.4; '> thresholdHTS' in V3.3)
3.40-3.42
		July 26 2013	overhauled U_SHAPE treatment, fix in curve flattening
						overhauled outlier detection (now generates a set of corrections from every point, then
						keeps minimum set (if several equal, keeps the one generated from the highest conc.)

3.30	July 25 2013	hardcoded parameters made as command-line keys
						-USHAPE=2 -MONINT=2 -BSHIFT=3
						-CRO= was allowed to be exactly 0 to avoid carryover check

3.20	July 23 2013	fix in POD reporting and lgConc conversion from input (units are treated as is)

3.11	June 11 2013	fix in fingerprint calculation
3.10	June 7 2013		fix in POD calculation, update in default file extension handling

3.04	June 3-6 2013	added alternative format support (array of conc, array of responses)
						added key for: carry over, baseline shift corrections, missing values
						modified output format
						new provisional metrics are printed via undocumented key -XTINF (ExtendedMode)

2.01	Jan 2 2012		added handling for spikes on the otherwise constant-response curve
						(big spikes could affect the average value used for rewriting all points)

2.00	Nov 18 2011		stricter out-of-range handling was added 
						CurveP calculation adjusted, alternative fingerprint mode was added "-ALT"
						(changes can be found by "11.18.2011" tag)
						help system and log messages were expanded
						'-TRSH=' was replaced to '-THR=' so as to correspond to the publication

1.01	Nov 2 2011		Basic version completed

*/

#include "core.h"
#include "qsar.h"

#define Version		"5.12"
#define COMMENT		"#"
#define	HTS_FILE	".hts"
#define	HTSX_FILE	".htsx"


UNSIGNED_1B_TYPE IGNORED_N_USHAPE = 4,		//for u-shape curves, min.#points to avoid flattening
				MIN_N_BASEP = 3;			//min. #baseline points needed to detect baseline shift
//bool ExtendedMode = false;				//false is temporary, 2013
bool TrustHighConc = false,				//assume more noise at high conc-s.
StrictExtrapolation = true;				//for POD, ECxx etc, used in Impute()
SIGNED_2B_TYPE Log10Input = 0;
UNSIGNED_1B_TYPE CONC_BASE = 6;			//uM are assumed

REALNUM_TYPE Impute(apvector<REALNUM_TYPE> &lgC, apvector<REALNUM_TYPE> &V, REALNUM_TYPE L, SIGNED_4B_TYPE JUNK = -999, bool strict = true)
{//returns estimated log conc., at which response level L is reached.
	SIGNED_4B_TYPE z, s = 0, e = lgC.length(); 
	while ( (lgC[--e] == JUNK) );
	while (lgC[s] == JUNK) s++;
	
	if (fabs(V[e]) < L) if (strict) return JUNK; 	
	if (fabs(V[s]) >= L) if (strict) return lgC[s];

	REALNUM_TYPE xT = JUNK;
	z = s;
	while (z++ < e)
	{
		if (lgC[z] == JUNK) continue;
					
		if ( (fabs(V[z]) >= L) || (z == e) )
		if ( fabs(V[z]) > fabs(V[s]) )
		{
			xT  = L - fabs(V[s]);
			xT /= fabs(V[z] - V[s]);
			xT *= lgC[z] - lgC[s];
			xT += lgC[s];
			break;
		}
		s = z;
	}

	return xT;
}

void handleHTSdata (STRING_TYPE inf, STRING_TYPE outf, STRING_TYPE tag, bool ifSuppl, REALNUM_TYPE thresholdHTS, REALNUM_TYPE fullRange, REALNUM_TYPE alwdDeviation, REALNUM_TYPE crOver = 0, bool fixBase = false, bool AltFingerprint = false, SIGNED_4B_TYPE DUMV = -999)
{	
	FILETYPE_IN inHTS(inf.c_str());	
	if (inHTS.eof() || inHTS.fail())
	{
		cout << "Cannot open '" << inf << "' file" << endl;
		return;
	}
	
	bool ifHTSX = false, ifIDs = false;
	STRING_TYPE xS = inf, bs = outf;
	xS.tolowercase();
	ifHTSX = CheckStrEnding(xS, HTSX_FILE);
	if (ifHTSX) bs += HTSX_FILE;	else	bs += HTS_FILE;

	SIGNED_4B_TYPE r, c, f, v, nRows, nCols, cIndex = 1;
	inHTS >> nRows; //#rows i.e. data entries
	inHTS >> nCols; //#columns, i.e. not counting first column of strings
	xS.getlinewithtabs(inHTS); //read the rest of 1st line

	QSAR QQ;
	apvector<STRING_TYPE> fields;
	apvector<REALNUM_TYPE> Conc(nCols), lgConc(nCols), HTS(nCols);
	STRING_TYPE BadPoints, Warn;

	xS.getlinewithtabs(inHTS);
	xS.parse_string();
	SplitString(xS, "\t", fields);
	cIndex = fields.length() - nCols*(1 + UNSIGNED_1B_TYPE(ifHTSX));
	if ( (cIndex <  1) || (cIndex > 2) )
	{
		cout << "wrong format in '" << inf << "' file" << endl;
		return;
	}
	else
		ifIDs = (cIndex == 2);
	
	FILETYPE_OUT outHTS(bs.c_str());
	if (ifSuppl)
	{
		GetTimeStamp(xS);
		outHTS << "#CurveP V" << Version << "; running on " << xS << endl;
		outHTS << "#Input: " << inf << endl;
		outHTS << "# Theoretical range of activity = 0 .. " << fullRange << endl;
		outHTS << "# Applied threshold (THR) = " << thresholdHTS << endl;
		outHTS << "# Limit of acceptable deviation (MXDV) = " << alwdDeviation << endl;
		outHTS << "# Carry-over threshold (CRO) = " << crOver << " (NB: also used for potency check of u-shapes)" << endl;
		outHTS << "# min.#points to detect u-shape (USHAPE) =  " << (SIGNED_4B_TYPE)IGNORED_N_USHAPE << endl;
		outHTS << "# min.#flat points to detect baseline-shift (BSHIFT) =  " << (SIGNED_4B_TYPE)MIN_N_BASEP << endl;
		outHTS << "# if several smallest sets of corrections found, based on " << (TrustHighConc ? "highest" : "lowest") << " conc. will be used" << endl;

		if (fixBase)
			outHTS << "# Baseline shift detection mode is ON" << endl;

		if (!StrictExtrapolation)
			outHTS << "# Extrapoation mode beyond test conc. range is ON" << endl;

		outHTS << "# Fill-in for missing  values = " << DUMV << endl;
		if (AltFingerprint) 
			outHTS << "# Alternative fingerprint mode is ON (uses baseline and 3 bins instead of 4-bin scheme)" << endl;
		outHTS << "# NB: points are adjusted on assumption of curve monotonicity." << endl;
		outHTS << "# NB: concentration-based parameters are output as log10 of Molar conc." << endl;
		outHTS << "#####################################################################" << endl;		
	}

	//print out header of the output file, prepare conc.
	//NB: tag will normally be '_assayname'
	outHTS << nRows << TAB << nCols << endl << fields[0] << TAB;
	if (ifIDs) outHTS << fields[1] << TAB;
	outHTS.precision(9);

	if (ifHTSX)
	{//output as original input		
		for (c = 0; c < nCols; c++) outHTS << "Conc" << c << tag << TAB;
		for (c = 0; c < nCols; c++) outHTS << "Resp" << c << tag << TAB;
	}
	else		
		for (c = 0; c < nCols; c++)
		{	
			xS = fields[c + cIndex];
			Conc[c] = atof(xS.c_str());
			
			if (xS.find('n') == INVALID) 
			{				
				if (xS.find('u') == INVALID)
				{
					if (xS.find("mM") == INVALID)
						Conc[c] *= 0.001; //convert milli to M
				}
				else
					Conc[c] *= 0.000001; //convert micro to M

			}
			else 
				Conc[c] *= 0.000000001; //convert nano to M
			outHTS << xS << tag << TAB;		//output original entries
		}	

	
	apvector<UNSIGNED_1B_TYPE> rlevels;
	rlevels.resize(12);
	rlevels[0] = 1;	rlevels[1] = 5; rlevels[2] = 10; rlevels[3] = 20; rlevels[4] = 25; rlevels[5] = 50;
	rlevels[6] = 75; rlevels[7] = 80; rlevels[8] = 90; rlevels[9] = 95; rlevels[10] = 99; rlevels[11] = 100; 
	
	outHTS << "CurvP" << tag << TAB << "log10CurvP" << tag << TAB;
	for (c = 0; c < rlevels.length(); c++) outHTS << "C@" << (SIGNED_4B_TYPE)rlevels[c] << "%" << tag << TAB;
	for (c = 0; c < rlevels.length(); c++) outHTS << "EC" << (SIGNED_4B_TYPE)rlevels[c] << tag << TAB;
	outHTS << "Emax" << tag << TAB << "EC50_slope" << tag << TAB << "wConc" << tag << TAB << "wResponse" << tag << TAB << "wAUC" << tag << TAB << "POD" << tag << TAB;
	outHTS << "Remarks" << tag << TAB << "#Corrections"  << tag << endl;

	for (r = 0; r < nRows; r++)
	{
		if (inHTS.eof())
		{
			cout << "Unexpected end in '" << inf << "' file" << endl;
			break;
		}

		xS.getlinewithtabs(inHTS);
		xS.parse_string();
		SplitString(xS, "\t", fields);
		if (ifIDs) 
		{	
			outHTS << fields[0] << TAB; //output current id
			BadPoints = fields[1]; 
		}
		else 
			BadPoints = fields[0];

		if (ifHTSX)
			for (c = 0; c < nCols; c++)	
			{				
				Conc[c] = atof(fields[c + cIndex].c_str());
				HTS[c] = atof(fields[c + cIndex + nCols].c_str());
				if (Conc[c] == DUMV) continue;

				if (Log10Input ==  0) 
					Conc[c] *= pow((REALNUM_TYPE)10, -CONC_BASE); //converts to moles
				else
					Conc[c] = pow(10, Conc[c]*Log10Input + CONC_BASE);	//log10 or -log10 				
			}
		else
			for (c = 0; c < nCols; c++)	HTS[c] = atof(fields[c + cIndex].c_str());

		for (c = 0; c < nCols; c++)	if (Conc[c] == DUMV) lgConc[c] = DUMV; else lgConc[c] = log10(Conc[c]);

		//read the missing points
		set Baddies;
		for (c = 0; c < nCols; c++)	if ((Conc[c] == DUMV) || (HTS[c] == DUMV)) Baddies.PutInSet( c );

		BadPoints.parse_string();
		BadPoints.replace(BLANK, "");
		BadPoints.replace("null", "");		 
		if ( BadPoints.length() )
		{
			if (BadPoints[0] == '{')
			{//{1.2,1.10} mask format (it is old, kept for compatibility)
				BadPoints += ',';
				c = BadPoints.find('.');
				while (c > INVALID)
				{
					BadPoints[c] = '_';
					c++;
					f = BadPoints.find(',');
					if (f == INVALID) break;
					xS = BadPoints.substr(c, f - c);
					c = atoi( xS.c_str() );
					if (c == ZERO) break;
					Baddies.PutInSet( --c );
					BadPoints[f] = '_';
					c = BadPoints.find('.');
				}
			}
			else
			{//unrecognizable mask formats will be ignored (which amounts to "all points are correct" mask)				
				//00..001000 format assumed
				f = nCols - BadPoints.length();
				if (f > 0) 
				{
					xS = "";
					for (f = c = 0; c < nCols; c++)	
					{//fix incomplete masks that skip unused test conc. (defined by DUMMY VALUES)
						if (HTS[c] == DUMV) { xS += "1"; continue; }						
						if ( (f < BadPoints.length())&&(('0'==BadPoints[f])||('1'==BadPoints[f])) ) xS += BadPoints[f++]; else xS += "0";
					}
					BadPoints = xS;
				}
				for (c = 0; c < BadPoints.length(); c++)	if (BadPoints[c] == '1')	Baddies.PutInSet(c);	
			}
		}//if ( BadPoints.length() )

		//------- process curve		
		for (c = 0; c < nCols; c++)
		{//corrects out of range values
			if ( HTS[c] == DUMV )	{ Baddies.PutInSet(c); continue; };
			if ( HTS[c] > max(0, fullRange) )	HTS[c] = max(0, fullRange);
			if ( HTS[c] < min(0, fullRange) )	HTS[c] = min(0, fullRange);
			if ( fabs( HTS[c] ) < thresholdHTS )	HTS[c] = ZERO;
		}

		set origBaddies(Baddies); //save original mask for later comparison

		//------- the curve should be monotonous, therefore, need to detect "violating" points
		apvector<REALNUM_TYPE> tVals(nCols);		
		for (f = c = 0; c < nCols; c++) if (!Baddies.IsInSet(c))	tVals[f++] = HTS[c];
		tVals.resize(f);
		
		REALNUM_TYPE xRange = 0, xWrk, tdff;
		if (nCols == Baddies.Size()) 
		{
			Warn = "NO_VALID_POINTS";
			for (c = 0; c < nCols; c++) HTS[c] = 0;	
			goto AHEAD;
		}

		f = 0; while ( Baddies.IsInSet(f) ) f++;
		c = nCols; while ( Baddies.IsInSet(--c) );
		xRange = HTS[f] - HTS[c];
		Warn = "";
		
		if (fabs(HTS[f]) > 0)
		{//blips can affect xRange, check next two conc.
			v = f;
			while ( Baddies.IsInSet(++v) ); 
			if (v < nCols) if (HTS[v] == 0) while ( Baddies.IsInSet(++v) );
			if (v < nCols) if (HTS[v] == 0) { xRange -= HTS[f]; HTS[f] = 0; } //erase starting point (likely blip)
		}

		if (fabs(xRange) < thresholdHTS)
		{//redo the slope-direction analysis differently
			xRange = 0;		//set the curve type as constant
			if (QQ.stdev(tVals) > alwdDeviation)
			{
				apvector<SIGNED_2B_TYPE> Ivals(nCols, 0); 	//first interval Ivals[0] is fake, it is from infinite dilution to the first conc. tested
				for (xWrk = v = c = 0; c < nCols; c++)
				{
					if ( Baddies.IsInSet(c) ) continue;					
					if (fabs(HTS[c] - xWrk) > alwdDeviation)
					{
						if (HTS[c] < xWrk) { v++; Ivals[c] =  1; } //decreasing
						if (HTS[c] > xWrk) { v--; Ivals[c] = -1; } //increasing
					}
					xWrk = HTS[c];
				}
				
				//flat or U-shaped curves					
				SIGNED_2B_TYPE cP, xP, vP, mP = 0; //work vars and optimal pivot value
				SIGNED_2B_TYPE errP = nCols, bP = 0; //for optimal pivot, store #corrections and pivot index
				for (f = c = 0; c < nCols; c++)						
				{
					if (Ivals[c] == 0) continue;
					f += Ivals[c]; v -= Ivals[c];
					cP = abs(v) + abs(f); //reflects the width of the spike
					if (cP < mP) continue;
					xP = vP = c;
					if (fullRange*(f - v) > 0)	xP++; 
					else {	while ((++vP < nCols)&&(Ivals[vP] == 0));	xP = nCols - vP--; }		
					if ((cP == mP) && (errP < xP)) continue;
					mP = cP; bP = vP; errP = xP;
				}					
				
				for (xP = vP = v = f = c = 0; c < nCols; c++)
				{
					if (Ivals[c] == 0)	continue;
					if (Ivals[c]*Ivals[f] < 0) v++; //count number of spikes (changes in monotonicity)
					if (c > bP) vP += Ivals[c]; else xP += Ivals[c];
					f = c;
				}
				if ((bP+1) == errP) f = abs(vP); else f = abs(xP);
				xWrk = fabs(HTS[bP]); if (crOver == 0)  xWrk = -1;
				if ( (mP < IGNORED_N_USHAPE) || (bP == 0) || (v > mP) || ((f<<1) < IGNORED_N_USHAPE) || ((xWrk < crOver)&&(errP > mP)) )
				{ Warn = " NOISY"; if (xWrk > crOver) Warn += " CHECK"; }
				else
				{//treat U-shape, its pivot-point represents peak/plateau
					Warn = " U_SHAPE"; if (errP > mP) Warn += " CHECK";
					if (fullRange > 0) xRange = -1; else xRange = 1;
					//invalidate wrong part to aid outlier-detection
					if (++bP == errP) {c = 0; f = bP; v = vP; }	else { c = bP; f = nCols; v = xP; }
					for (; c < f; c++) Baddies.PutInSet(c);
				}
			} //if (QQ.stdev(tVals) > alwdDeviation) 
			else
			{//04.26.2014 fix
				xWrk = QQ.meanV(tVals);
				if (xWrk >= thresholdHTS)
				{//const curve with low variance and signif signal: can be baseline shift or carry over, do not erase
					for (c = 0; c < nCols; c++) 
					{
						if (Conc[c] == DUMV) continue; 
						if (fabs(HTS[c] - xWrk) > alwdDeviation) { Baddies.PutInSet(c); HTS[c] = xWrk; }
					}
					goto AHEAD;
				}				
			}
			
		}//if (fabs(xRange) < thresholdHTS)

		//Now xRange stores a range of the curve, if < 0 then it's rising
		if (xRange == 0)
		{//if still flat, make it so			
			for (c = 0; c < nCols; c++) 
			{
				if (Conc[c] == DUMV) continue; 
				if (fabs(HTS[c]) > 0) { Baddies.PutInSet(c); HTS[c] = 0; }
			}
			goto AHEAD;
		}
		else
		{//Detect a minimum set of violating points
			set TrialBest(0, nCols);
			UNSIGNED_2B_TYPE tbSize = nCols, bdSize = Baddies.Size();
			
			for (v = 0; v < nCols; v++)
			{ //v seeds the starting point to trust as "true"
				if ( Baddies.IsInSet(v) ) continue;
				f = v;
				set Trial(Baddies);
				for (c = f+1; c < nCols; c++)
				{//detecting discrepancies from v to the right
					if ( Baddies.IsInSet(c) ) continue;				
					tdff = HTS[c] - HTS[f];
					if ( ((xRange * tdff) > 0) && (fabs(tdff) > alwdDeviation) )
					{//bad point
						Trial.PutInSet(c); continue;
					}				
					f = c;
				}//for c

				f = c = v;			
				while (c-- > 0)
				{//detecting discrepancies on the curve from v to left
					if ( Baddies.IsInSet(c) ) continue;				
					tdff = HTS[f] - HTS[c];
					if ( ((xRange * tdff) > 0) && (fabs(tdff) > alwdDeviation) )
					{//bad point
						Trial.PutInSet(c); continue;
					}				
					f = c;
				}//while c

				f = Trial.Size();
				if (tbSize < f) continue;
				if ( (tbSize > f) || TrustHighConc ) { TrialBest = Trial; tbSize = f; } //update set
				if (!TrustHighConc)		if (tbSize == bdSize) break; //optimum reached
			}//for v
			
			Baddies = TrialBest;	//output mask
		} //else..if (xRange == 0)

		//------- replace bad points with appropriate extrapolations
		for (c = 0; c < nCols; c++)
		{
			if (Conc[c] == DUMV) continue;
			if ( Baddies.IsInSet(c) )
			{
				f = v = c; 
				while ( Baddies.IsInSet(--v) );
				while ( Baddies.IsInSet(++f) );
				if ( v == INVALID )
				{
					if (f < nCols) HTS[c] = HTS[f]; else { Warn += " NO_VALID_POINTS"; HTS[c] = ZERO; }	//all points are invalidated, set to baseline
				}
				else
				{
					if ( f == nCols )
						HTS[c] = HTS[v];	
					else
					{
						HTS[c]	 = lgConc[c] - lgConc[v];
						HTS[c]	/= lgConc[f] - lgConc[v];
						HTS[c]	*= HTS[f] - HTS[v];
						HTS[c]	+= HTS[v];
					}
				}				
			}				
		}//for c

AHEAD:
		for (f = v = c = 0; c < nCols; c++)
		{
			if (HTS[c] == 0) continue;
			if (HTS[c] == DUMV) continue;
			f++;
			if (Baddies.IsInSet(c)) continue;			
			v++;
		}
		if ((v == 1) && (f > 1)) Warn += " SINGLE_POINT_ACT";

		tdff = xRange * fullRange;
		if (tdff > 0) Warn += " INVERSE";

		//------- detection & handling of carry-over and related problems
		if ( (crOver > 0) && (v > 0) )
		{
			c = 0; while (Conc[c] == DUMV) c++;
			f = nCols; while ( (Conc[--f] == DUMV) );
			xWrk = fabs(HTS[c]);						
			if (xWrk > 0)
			{
					if (tdff > 0)
					{//decrease in signal, unconditional carryover if inhibitor or potency-conditional if agonist
						if ( (fullRange < 0) || (xWrk < crOver) )
						{
							Warn += " CARRY_OVER";
							while (c <= f) if (Conc[c] == DUMV) c++; else { if (HTS[c] == 0) break; Baddies.PutInSet(c); HTS[c++] = 0; }
						}
						else Warn += " CHECK";
					}
					else //increase or constant signal					
						if (xWrk < crOver)										
						{
							if (xRange == 0)
							{//constant, likely carryover or baseline shift
								Warn += " CARRY_OVER? BASE_SHIFT?";
								while (c <= f) if (Conc[c] == DUMV) c++; else { if (HTS[c] == 0) break; Baddies.PutInSet(c); HTS[c++] = 0; }
							}
							else
							{//increasing, can be potent active, carry over or baseline shift
								tVals.resize(1);
								tVals[0] = HTS[c];
								f = ++c;
								//advanced baseline detection, 09.09.2013
								for (v = 0; f < nCols; f++)
								{
									if (Conc[f] == DUMV) continue;
									c = tVals.length();
									tVals.resize( ++c );
									tVals[c - 1] = HTS[f];
									if (QQ.stdev(tVals) < alwdDeviation) v = c;									
								}
								tVals.resize( v );
								if (v < MIN_N_BASEP)
									Warn += " CHECK";
								else
								{//baseline shift, adjust									
									Warn += " BASE_SHIFT?";
									xWrk = QQ.meanV(tVals);									
									for (f = 0; f < nCols; f++)
									{
										if (Conc[f] == DUMV) continue;
										v--;
										if (fixBase) 
										{
											HTS[f] -= xWrk;		//fix as base-shift
											if ( fabs( HTS[f] ) < thresholdHTS )	HTS[f] = 0;
										}
										else 
										{											
											HTS[f] = 0;
											Baddies.PutInSet(f); //fix as carry-over
											if (v == 0) break;
										}										
									}																		
								}								
							}
						} //if (xWrk < crOver)	
						else Warn += " TOO_POTENT";					
			} //if (xWrk > 0)
		}//if ( (crOver > 0) && (v > 0) )
		

		//------- print out updated mask, and adjusted data
		xS = "";
		for (c = 0; c < nCols; c++) if (Baddies.IsInSet(c))	xS += "1 "; else xS += "0 ";
		xS.parse_string();
		outHTS << xS << TAB;
		if (ifHTSX)	for (c = 0; c < nCols; c++) outHTS << Conc[c] << TAB;	//print concentrations (extended format only)
		for (c = 0; c < nCols; c++)	outHTS << HTS[c] << TAB;	//print corrected responses


		REALNUM_TYPE Emax = 0, POD = DUMV, slope = 0;
		apvector<REALNUM_TYPE> CAs(rlevels.length(), DUMV), ECs(rlevels.length(), DUMV);

		//------- curve fingerprint
		UNSIGNED_4B_TYPE fn = 0;
		xRange	 = fabs(fullRange);		//cannot use local range, because a particular entry may have no activity
		if (AltFingerprint) xRange	/= 3; else xRange	/= 4;		
		for(c = 0; c < nCols; c++)
		{
			if ( Conc[c] == DUMV) continue; //skips "invalid" bits

			UNSIGNED_1B_TYPE lfn = 0;
			if ( HTS[c] == DUMV) tdff = 0; else tdff = fabs(HTS[c]);

			if (AltFingerprint)
			{//alternative scheme: 0 - for baseline, 1 for <= 1/3 of fullRange, 2 for <= 2/3, 3 for >2/3				
				if ( tdff > 0) lfn++;
				if ( tdff > xRange) lfn++;
				if ( tdff > 2*xRange) lfn++;
			}
			else
			{//basic scheme, 0 for <= 1/4, 1 for <= 2/4, 2 for <=3/4, 3 for >3/4
				if ( tdff > xRange) lfn++;
				if ( tdff > 2*xRange) lfn++;
				if ( tdff > 3*xRange) lfn++;
			}
			
			fn <<= 2;
			fn += lfn;			
		}

		//------- print curve fingerprint and its log10
		outHTS << fn << TAB; 		
		outHTS << (fn ? log10((long double)fn) : -1) << TAB;

		//-------- print C@.. and EC..		
		if (fn > 0)
		{//some signal
			for(Emax = c = 0; c < nCols; c++)
			{
				if ( Conc[c] == DUMV) continue;
				if ( (Emax > HTS[c])^(fullRange > 0) ) Emax= HTS[c];
			}
			slope = Emax/100;
			xWrk = fabs(slope);
			for(c = 0; c < rlevels.length(); c++)
			{
				CAs[c] = Impute(lgConc, HTS, rlevels[c], DUMV, StrictExtrapolation);
				ECs[c] = Impute(lgConc, HTS,  xWrk*rlevels[c], DUMV, StrictExtrapolation);
			}
			POD = Impute(lgConc, HTS, thresholdHTS, DUMV, StrictExtrapolation);
			slope *= 50/(ECs[6] - ECs[4]);
		}

		for (c = 0; c < rlevels.length(); c++) outHTS << CAs[c] << TAB;
		for (c = 0; c < rlevels.length(); c++) outHTS << ECs[c] << TAB;

		//------- addl.metrics
		xWrk = tdff = xRange = 0;
		for (c = 0; c < nCols; c++)
		{
			if (Conc[c] == DUMV) continue;
			tdff += lgConc[c]; 
			xWrk += lgConc[c]*HTS[c];
			xRange += HTS[c];
		}

		outHTS << Emax << TAB << slope  << TAB <<((xRange == 0) ? 0 : xWrk/xRange) << TAB << xWrk/tdff << TAB;	//Emax, slope, wConc, wResponse

		if (POD == DUMV)
			xWrk = 0;
		else
		{//---- wAUC
			f = nCols; while ( (Conc[--f] == DUMV) );
			c = 0; while (Conc[c] == DUMV) c++;
			xRange = 2*(lgConc[c] - lgConc[f]);	//scaling for wAUC
			for (xWrk =0, v = c+1; v <= f; v++)
			{
				if (Conc[v] == DUMV) continue;
				xWrk += (lgConc[v] - lgConc[c])*(HTS[v] + HTS[c]); //NB: all HTS values should be of the same sign at this point
				c = v;
			}			
			xWrk *= POD/xRange;
		}
		outHTS << xWrk << TAB << POD << TAB;	//print wAUC, POD

		//print remarks
		Warn.parse_string();
		if (Warn.length() == 0) Warn = "OK";
		Warn.replace(" ", "|");
		outHTS << Warn << TAB << (Baddies.Size() - origBaddies.Size()) << endl;
	}//for r

	inHTS.close();
	outHTS.close();
	cout << "Done." << endl;
}

int main(int argc, char* argv[])
{
	REALNUM_TYPE THR = 15, MXDV = 5, RANGE = -100, CRO = 80;	//will be used for noise-reduction
	SIGNED_4B_TYPE DUMMYVAL = -999;

	bool DetailedOutput = false;
	bool AltcurveP = false;
	bool BSHIFT = false;

	if (argc < 2)
	{//print terse help
		cout << endl << "#CurveP V" << Version << " - qHTS curve processing and noise-reduction" << endl;
		cout << "Sedykh et al, Environ Health Perspect (2011), 119, p.364" << endl << endl;

		cout << "Usage:     curvep filename [flags]" << endl;
		cout << "Allowed inputs are tab-delimited text files of following formats:" << endl;
		cout << ".hts file: 1st line - #samples and #test.concentrations (e.g., 100  8)" << endl;		
		cout << "2nd line - column headers, 'ID(optional) ExclPoints 0.6nM_% 3uM_% 3mM_%'.." << endl;
		cout << "Following lines: ids(optional), mask, responses at each test concentration" << endl;
		cout << ".htsx file: 1st line - same as in .hts file" << endl;
		cout << "2nd line - column headers, 'ID(optional) ExclPoints Conc1 Conc2..Resp1 Resp2..'" << endl;
		cout << "Following lines: ids(optional), mask, test conc-s, responses" << endl;
		cout << "Mask is either in '{1.2, 1.4}' format or a binary string" << endl;
		cout << "('00010010' or '0 0 0 1 0 0 1 0') were 1's denote invalidated points" << endl;
		cout << "NB: Mask in output file is extended to include corrected points" << endl;
		cout << "NB: Conc.-based output values are in log10(M)." << endl << endl;

		cout << "Optional Flags:" << endl;
		cout << "'-HTSX=' - input units: 'nM' 'uM'<def.> 'M' 'lgM' etc." << endl;
		cout << "'-OUT=' output file name; '-LOG' detailed log; '-ALT' altern.fingerprint" << endl;
		cout << "'-MXDV=' max.deviation <def. 5.0>; '-THR=' baseline threshold <def.15.0>" << endl;
		cout << "'-RNG=' max.response <def. " << RANGE << "> (can be positive)" << endl;
		cout << "'-CRO=' carryover threshold <def. 80% of max.response>, not used if 0" << endl;
		cout << "'-BLFX' baseline shift correction mode, off by def." << endl;
		cout << "'-USHAPE=' min.#points for u-shape <def.4>" << endl;
		cout << "'-BSHIFT=' min.#flat points to detect baseline-shift <def.3>" << endl;
		cout << "'-BYHI','-BYLO'<def> favors corrections based on hi or lo conc-s" << endl;
		//cout << "'-XTINF' extended output with additional metrics" << endl;
		cout << "'-XPLAX' allow extrapolation beyond test conc. boundaries" << endl;
		cout << "'-NA=' dummy value for missing data <def. -999>" << endl;
		cout << "'-TAG=' optional marker to add to column headers of output file" << endl << endl;		
		return 0;
	}

	STRING_TYPE stInput = argv[1], stOutput = argv[1], stArg, stJ, stTAG = "";
	stJ = stInput;
	stJ.tolowercase();
	
	if (CheckStrEnding(stJ, HTS_FILE) || CheckStrEnding(stJ, HTSX_FILE))	//if right file extension is there, adjust default output file name		
		CutStrEnding(stOutput);		
	else
	{//if no extension then try to add appropriate one
		stInput += HTS_FILE;	
		FILETYPE_IN fiTest(stInput.c_str());
		if (fiTest.eof() || fiTest.fail())
		{//revert to extended format
			CutStrEnding(stInput);
			stInput += HTSX_FILE;
		}
		fiTest.close();
	}
	
	stOutput += "-mdf";
		
	REALNUM_TYPE rtX;
	SIGNED_4B_TYPE nArg = 1, intU;
	
	while (argc > ++nArg) 	
	{ 
		stArg = argv[nArg];
		stArg.parse_string();
		stJ = stArg;
		stArg.touppercase();
		
		intU = stArg.find("-OUT="); //output file
		if (intU >= 0)
		{//use original case!
			stOutput = stJ.substr(intU + 5, stJ.length());
			if (CheckStrEnding(stOutput, HTS_FILE) || CheckStrEnding(stOutput, HTSX_FILE))
				CutStrEnding(stOutput);
		}
		
		intU = stArg.find("-TAG="); //tag label
		if (intU >= 0)	//use original case!
		{
			stTAG = stJ.substr(intU + 5, stJ.length());
			stTAG.parse_string();
		}

		if (stArg.find("-LOG") >= 0) DetailedOutput = true;
		if (stArg.find("-ALT") >= 0) AltcurveP = true;
	
		intU = stArg.find("-MXDV=");
		if (intU == 0)
		{
			rtX = atof( stArg.substr(intU + 6, stJ.length()).c_str() );
			if (rtX > 0)	MXDV  = rtX;
		}

		intU = stArg.find("-THR=");
		if (intU == 0)
		{
			rtX = atof( stArg.substr(intU + 5, stJ.length()).c_str() );
			if (rtX > 0)	THR  = rtX;
		}
		
		intU = stArg.find("-RNG=");
		if (intU == 0)
		{
			rtX = atof( stArg.substr(intU + 5, stJ.length()).c_str() );
			if (fabs(rtX) > 0)	RANGE  = rtX;
		}

		intU = stArg.find("-NA=");
		if (intU == 0)	DUMMYVAL = atoi( stArg.substr(intU + 4, stJ.length()).c_str() );

		intU = stArg.find("-CRO=");
		if (intU == 0)
		{
			rtX = atof( stArg.substr(intU + 5, stJ.length()).c_str() );
			CRO  = rtX;
		}

		intU = stArg.find("-HTSX=");
		if (intU == 0)
		{
			STRING_TYPE inp = stJ.substr(intU + 6, stJ.length());
			
			if (inp.substr(0, 2) == "lg") { Log10Input = 1; inp.replace("lg", ""); };
			if (inp.substr(0, 1) == "p") { Log10Input = -1; inp.replace("p", ""); };
			if (inp == "nM") CONC_BASE = 9;
			if (inp == "uM") CONC_BASE = 6;	//default
			if (inp == "mM") CONC_BASE = 3;
			if (inp == "M") CONC_BASE = 0;
		}

		if (0 == stArg.find("-XPLAX"))	StrictExtrapolation =false;

		if (0 == stArg.find("-BLFX"))	BSHIFT = true;
		//if (0 == stArg.find("-XTINF"))	ExtendedMode= true; //temporary, 06.06.2013	

		intU = stArg.find("-USHAPE=");
		if (intU == 0)	IGNORED_N_USHAPE = atoi( stArg.substr(intU + 8, stJ.length()).c_str() );

		intU = stArg.find("-BSHIFT=");
		if (intU == 0) MIN_N_BASEP	= atoi( stArg.substr(intU + 8, stJ.length()).c_str() );

		if (stArg.find("-BYHI") == 0) TrustHighConc = true;
		if (stArg.find("-BYLO") == 0) TrustHighConc = false;				
	}//while (argc > ++nArg) loop

	if (stTAG.length() > 0)	
		if (stTAG[0] != '_')
			stTAG = '_' + stTAG;

	handleHTSdata(stInput, stOutput, stTAG, DetailedOutput, THR, RANGE, MXDV, CRO, BSHIFT, AltcurveP, DUMMYVAL);
}

