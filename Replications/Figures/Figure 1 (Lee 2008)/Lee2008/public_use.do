#delimit;
set more off;
set mem 1G;

log using public_use.log, replace;

program drop _all;

************************************************************************************************;
/*FIGURE ONE*/;
************************************************************************************************;

use figure1.dta, clear;

list;

twoway (connected col1 year, ms(s) mc(gs0) msize(medium) mfcolor(gs0) mlcolor(gs0) lpattern(solid) lcolor(gs0)) || 
				(connected col2 year, ms(sh) mc(gs0) msize(medium) mfcolor(none) mlcolor(gs0) lpattern(dash) lcolor(gs0)) || 
				(connected col3 year, ms(x) mc(gs0) msize(medium) mfcolor(none) mlcolor(gs0) lpattern(dash) lcolor(gs0)),
				title("Figure 1", color(gs0)) xtitle("Year") ytitle("Proportion Winning Election") legend(label(1 "Incumbent Party") label(2 
				"Winning Candidate") label(3 "Runner-up Candidate")) yscale(r(0 1)) xscale(r(1948 1998)) ylabel(0(.1)1) 
				xlabel(1948(10)1998) saving("./final/Fig1.gph", replace); 


************************************************************************************************;
/*FIGURES TWO AND THREE*/;
************************************************************************************************;

use individ_final.dta, clear;

graph twoway (scatter mmyoutcomenext difshare, ms(O) mc(gs0) msize(small)) (line mpmyoutcomenext difshare, sort lcolor(gs0)) 
				if difshare>-0.251 & difshare<0.251 & use==1, xline(0, lcolor(gs0)) title("Figure 2a", color(gs0)) xtitle("Democratic Vote Share Margin of Victory, Election t") 
				ytitle("Probability of Winning, Election t+1") yscale(r(0 1)) ylabel(0(.1)1) xscale(r(-0.25 0.25)) xlabel(-0.25(.05)0.25) 
				legend(label(1 "Local Average") label(2 "Logit Fit")) saving("./final/Fig2a.gph", replace);

graph twoway (scatter mofficeexp difshare, ms(O) mc(gs0) msize(small)) (line mpofficeexp difshare, sort lcolor(gs0)) 
				if difshare>-0.251 & difshare<0.251 & use==1, xline(0, lcolor(gs0))
				title("Figure 2b", color(gs0)) xtitle("Democratic Vote Share Margin of Victory, Election t") ytitle("No. of Past Victories as of Election t")
				yscale(r(0 5)) ylabel(0(.5)5) xscale(r(-0.25 0.25)) xlabel(-0.25(.05)0.25) 
				legend(label(1 "Local Average") label(2 "Polynomial Fit")) saving("./final/Fig2b.gph", replace);

graph twoway (scatter mrunagain difshare, ms(O) mc(gs0) msize(small)) (line mprunagain difshare, sort lcolor(gs0)) 
				if difshare>-0.251 & difshare<0.251 & use==1, xline(0, lcolor(gs0))
				title("Figure 3a", color(gs0)) xtitle("Democratic Vote Share Margin of Victory, Election t") 
				ytitle("Probability of Candidacy, Election t+1") 
				yscale(r(0 1)) ylabel(0(.1)1) xscale(r(-0.25 0.25)) xlabel(-0.25(.05)0.25)
				legend(label(1 "Local Average") label(2 "Logit Fit")) saving("./final/Fig3a.gph", replace);

graph twoway (scatter melectexp difshare, ms(O) mc(gs0) msize(small)) (line mpelectexp difshare, sort lcolor(gs0)) 
				if difshare>-0.251 & difshare<0.251 & use==1, xline(0, lcolor(gs0)) 
				title("Figure 3b", color(gs0)) xtitle("Democratic Vote Share Margin of Victory, Election t") 
				ytitle("No. of Past Attempts as of Election t") 
				yscale(r(0 5)) ylabel(0(.5)5) xscale(r(-0.25 0.25)) xlabel(-0.25(.05)0.25)
				legend(label(1 "Local Average") label(2 "Polynomial Fit")) saving("./final/Fig3b.gph", replace);

sort difgrp use;
quietly by difgrp: keep if _n==1;

outsheet mpmyoutcomenext mmyoutcomenext mprunagain mrunagain mpofficeexp mofficeexp mpelectexp melectexp difgrp 
			difshare using fig23 if difgrp>=150 & difgrp<=250, replace comma;
			

************************************************************************************************;
/*FIGURES FOUR AND FIVE*/;
************************************************************************************************;

use group_final, clear;

graph twoway (scatter mdemsharenext difdemshare, ms(O) mc(gs0) msize(small)) (line mpdemsharenext difdemshare, sort lcolor(gs0)) if 
				difdemshare>-0.251 & difdemshare<0.251 & use==1, 
				xline(0, lcolor(gs0)) title("Figure 4a", color(gs0)) xtitle("Democratic Vote Share Margin of Victory, Election t") 
				ytitle("Vote Share, Election t+1") 
				yscale(r(0.30 0.70)) ylabel(0.3(.05)0.7) xscale(r(-0.25 0.25)) xlabel(-0.25(.05)0.25)
				legend(label(1 "Local Average") label(2 "Polynomial Fit")) saving("./final/Fig4a.gph", replace);

graph twoway (scatter mdemshareprev difdemshare, ms(O) mc(gs0) msize(small)) (line mpdemshareprev difdemshare, sort lcolor(gs0)) 
				if difdemshare>-0.251 & difdemshare<0.251 & use==1, 
				xline(0, lcolor(gs0)) title("Figure 4b", color(gs0)) xtitle("Democratic Vote Share Margin of Victory, Election t") 
				ytitle("Vote Share, Election t-1") 
				yscale(r(0.30 0.70)) ylabel(0.3(.05)0.7) xscale(r(-0.25 0.25)) xlabel(-0.25(.05)0.25)
				legend(label(1 "Local Average") label(2 "Polynomial Fit")) saving("./final/Fig4b.gph", replace);

graph twoway (scatter mdemwinnext difdemshare, ms(O) mc(gs0) msize(small)) (line mpdemwinnext difdemshare, sort lcolor(gs0)) 
				if difdemshare>-0.251 & difdemshare<0.251 & use==1, 
				xline(0, lcolor(gs0)) title("Figure 5a", color(gs0)) xtitle("Democratic Vote Share Margin of Victory, Election t") 
				ytitle("Probability of Victory, Election t+1") 
				yscale(r(0 1)) ylabel(0(.1)1) xscale(r(-0.25 0.25)) xlabel(-0.25(.05)0.25)
				legend(label(1 "Local Average") label(2 "Logit Fit")) saving("./final/Fig5a.gph", replace);

graph twoway (scatter mdemwinprev difdemshare, ms(O) mc(gs0) msize(small)) (line mpdemwinprev difdemshare, sort lcolor(gs0)) 
				if difdemshare>-0.251 & difdemshare<0.251 & use==1, 
				xline(0, lcolor(gs0)) title("Figure 5b", color(gs0)) xtitle("Democratic Vote Share Margin of Victory, Election t") 
				ytitle("Probability of Victory, Election t+1") 
				yscale(r(0 1)) ylabel(0(.1)1) xscale(r(-0.25 0.25)) xlabel(-0.25(.05)0.25)
				legend(label(1 "Local Average") label(2 "Logit Fit")) saving("./final/Fig5b.gph", replace);

sort difgrp use;
quietly by difgrp: keep if _n==1;

outsheet mpdemsharenext mdemsharenext mpdemwinnext mdemwinnext mpdemshareprev mdemshareprev mpdemwinprev mdemwinprev 
			difgrp difdemshare using fig45 if difgrp>=150 & difgrp<=250, replace comma;
			
************************************************************************************************;
/*TABLE ONE*/;
************************************************************************************************;

use table_one_final, clear;

outsheet Variable Col_1 Col_2 Col_3 Col_4 Col_5 Col_6 Col_7 Col_8
	 using "./final/table1.csv", comma replace;



************************************************************************************************;
/*TABLE TWO*/;
************************************************************************************************;
program drop _all;

program define fitpoly;
	regress `1' difdemshare difdemshare2 difdemshare3 difdemshare4 rdifdemshare
	rdifdemshare2 rdifdemshare3 rdifdemshare4 right `2' if use==1, cluster(statedisdec);
end;

use table_two_final, clear;

/*column one*/;
fitpoly demsharenext;
estimates store c1;

/*column two*/;
fitpoly demsharenext "demshareprev demwinprev";
estimates store c2;

/*column three*/;
fitpoly demsharenext "demofficeexp othofficeexp";
estimates store c3;

/*column four*/;
fitpoly demsharenext "demelectexp othelectexp";
estimates store c4;

/*column five*/;
fitpoly demsharenext "demshareprev demwinprev demofficeexp othofficeexp demelectexp othelectexp";
estimates store c5;

reg demsharenext demshareprev demwinprev demofficeexp othofficeexp demelectexp othelectexp if use==1;
predict demsharenextres if use==1, residuals;

/*column six*/;
fitpoly demsharenextres;
estimates store c6;

/*column seven*/;
fitpoly difdemsharenext "demwinprev demofficeexp othofficeexp demelectexp othelectexp";
estimates store c7;

/*column eight*/;
fitpoly demshareprev "demwinprev demofficeexp othofficeexp demelectexp othelectexp";
estimates store c8;

keep _est_c1 _est_c2 _est_c3 _est_c4 _est_c5 _est_c6 _est_c7 _est_c8;

esttab c1 c2 c3 c4 c5 c6 c7 c8 using "./final/table2.rtf", replace b(3) se(3) noconstant nostar noobs 
		title({\b "Table 2"} {\i "Effect of winning an election on subsequent party electoral success"}) 
		mtitles("Vote share t+1" "Vote share t+1" "Vote share t+1"
		"Vote share t+1" "Vote share t+1" "Res. vote share t+1" "1st dif. vote share, t+1" 
		"Vote share t-1") nodepvars coeflabels(right "Victory, election t" demshareprev 
		"Dem. vote share, t-1" demwinprev "Dem. win, t-1" demofficeexp "Dem. political experience"
		othofficeexp "Opp. political experience" demelectexp "Dem. electoral experience" 
		othelectexp "Opp. electoral experience") drop(difdemshare* rdifdemshare* _cons)
		addnotes("Note: Details of data processing in Appendix A. N=6558 in all regressions.");

log close;
clear;
