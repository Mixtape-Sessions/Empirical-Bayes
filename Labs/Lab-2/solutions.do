

*******************************************
**** Chris Walters
**** 10/2023
**** This program applies empirical Bayes methods
**** to study employment discrimination using data
**** from Kline, Rose, and Walters (QJE 2022)
********************************************


	********************
	*** Basic setup ****
	********************
	
		*Set up stata
		clear
		cap set more off
		cap set trace off
		cap log close
		set seed 1028
    /* cd "Labs/Lab-2" */
		
		*Switches
		local estimate_gaps=0
		local multiple_test=1
		
	**************************************
	*** Load and format data
	***************************************
		
		****Load data
		insheet using "krw_data.csv", names comma clear
		
		****Summarize variables in this data set
		sum	
		
		*****Overall treatment effects -- compare robust and job-clustered SE
		reg callback white, r
		reg callback white, cluster(job)
	
		****Collapse to data set of firm-specific gaps and SEs
		gen callback_white=callback if white==1
		gen callback_black=callback if white==0
		bys job: egen ybar_white=mean(callback_white)
		bys job: egen ybar_black=mean(callback_black)
		bys job: gen n_j=_N
		bys job: keep if _n==1
		gen thetahat_jf=ybar_white - ybar_black
		drop if thetahat_jf==.
		bys firm: egen thetahat_f = mean(thetahat_jf)
		bys firm: egen std=sd(thetahat_jf)
		bys firm: gen J=_N
		gen s_f=sqrt((1/J)*(std^2))
		bys firm: keep if _n==1
		keep thetahat_f s_f
		
		***Save temporary file to load in either exercise
		tempfile basic_estimates
		save "`basic_estimates'", replace
		
	**************************************
	*** Estimate distribution of firm race gaps
	***************************************
		if `estimate_gaps'==1 {
		
			***Load estimates
			use "`basic_estimates'", clear

			****Estimate variance of mixing distribution
			qui sum thetahat_f
			local mu=r(mean)
			local F=r(N)
			gen v=(thetahat_f - (`mu'))^2 - (((`F'-1)/`F')*(s_f^2))
			qui sum v
			local sigma=sqrt(r(mean))
			disp "Estimated SD of firm gaps: `sigma'"

			
			****Linear shrinkage estimates
			gen lambda=((`sigma')^2)/(((`sigma')^2)+(s_f^2))
			gen thetastar=lambda*thetahat_f + (1-lambda)*(`mu')

			****Histogram of estimates vs. linear shrinkage posteriors
			graph twoway hist thetahat_f,  color(navy) fintensity(0) width(0.005) ///
						|| hist thetastar, width(0.005) barwidth(0.0025) color(navy) ///
						legend(lab(1 "Unbiased estimates") lab(2 "Linear shrinkage estimates") order(1 2)) ///
						xtitle("White/Black contact gap") ytitle("Number of firms") ///
						scheme(s1color)
						
			****Execute R command for nonparametric deconvolution
			gen id_f=_n
			preserve
			keep id_f thetahat_f s_f
			order id_f thetahat_f s_f
			outsheet using "disc_estimates.csv", nonames comma replace
			rsource using "decon_example.R", rpath("/usr/local/bin/R") roptions(`"--vanilla"')
			restore
		
			****Histogram of estimates with nonparametric mixing distribution estimate
			preserve
			gen sample1=1
			tempfile sample
			save "`sample'"
			insheet using "decon_estimates.csv", nonames comma clear
			ren v1 supp_theta
			ren v2 g_theta
			gen sample1=0
			append using "`sample'"
			graph twoway hist thetahat if sample1==1,  color(navy) fintensity(0) width(0.005) density ///
						|| line g_theta supp_theta if sample1==0,  lcolor(maroon) ///
						legend(lab(1 "Unbiased estimates") lab(2 "Nonparametric prior") order(1 2)) ///
						xtitle("White/Black contact gap") ytitle("Density") ///
						scheme(s1color)	
			restore
			
			*****Plot of linear shrinkage vs. non-parametric posterior means
		
				*Bring in non-parametric posteriors from R
				tempfile mergefile
				save "`mergefile'", replace
				insheet using "postmean_estimates.csv", nonames comma clear
				ren v1 id_f
				ren v2 thetastar_nonpar
				merge 1:1 id_f using "`mergefile'"
				tab _merge
				drop _merge
				
				*Plot linear shrinkage vs. nonparametric posteriors
				graph twoway scatter thetastar thetastar_nonpar, color(navy) ///
					ytitle("Linear shrinkage estimate") xtitle("Nonparametric posterior mean") ///
					scheme(s1color)

				
			preserve
			gen sample1=1
			tempfile sample
			save "`sample'"
			insheet using "decon_estimates.csv", nonames comma clear
			ren v1 supp_theta
			ren v2 g_theta
			gen sample1=0
			append using "`sample'"
			graph twoway hist thetahat if sample1==1,  color(navy) fintensity(0) width(0.005) density ///
						|| line g_theta supp_theta if sample1==0,  lcolor(maroon) ///
						legend(lab(1 "Unbiased estimates") lab(2 "Nonparametric prior") order(1 2)) ///
						xtitle("White/Black contact gap") ytitle("Density") ///
						scheme(s1color)	
			restore
			
    }
	
	**************************************
	*** Multiple testing: Classify firm discrimination while controlling FDR
	***************************************
		if `multiple_test'==1 {
		
			***Load estimates
			use "`basic_estimates'", clear
		
			*Choose cutoff lambda for calculating bound on share of true nulls (pi0)
			local lambda=0.5
			
			*Choose False Discovery Rate control threshold
			local FDR=0.05
		
			*Compute one-tailed tests of H0: theta=0 vs. HA: theta>0
			gen z_f=thetahat_f/s_f
			gen p_f=1-normal(z_f)
						
			*Bound pi_0
			gen c=(p_f*(p_f>`lambda'))/(1-`lambda')
			qui sum c
			local pi_0=r(mean)
			disp "Bound on pi_0: `pi_0'"
				
			*Compute q-values
			sort p_f
			gen F_p=_n/_N
			gen q_f=(p_f*(`pi_0'))/F_p
				
			*Count firms with q-vals below FDR threshold, and find p-val cutoff
			count if q_f<`FDR'
			local N_firms=r(N)
			sum p_f if q_f<`FDR', detail
			local p_cutoff=r(max)
			
			*Plot histogram of p-vals
			local pi0=round(`pi_0',0.01)
			local lambda_loc=`lambda'+0.06
			local pi0_loc=`pi0'+0.2
			local p_loc=`p_cutoff'+0.175
			graph twoway hist p_f, color(navy) fintensity(0) width(0.05) /// 
					legend(off) ///
					xtitle("P-value from test of no disc. against Black applicants") ytitle("Density") ///
					scheme(s1color) 
					
				
					
			*Plot histogram of p-vals with lambda, pi0, and p-val cutoff for FDR control
			graph twoway hist p_f, color(navy) fintensity(0) width(0.05) /// 
					|| function y=`pi_0', range(0 1) lcolor(maroon) lpattern(dash) ///
					xline(`lambda', lcolor(black) lpattern(dash)) ///
					legend(off) ///
					xtitle("P-value from test of no disc. against Black applicants") ytitle("Density") ///
					scheme(s1color) ///
					text(3.5 `lambda_loc' "{&lambda} = 0`lambda'") ///
					text(`pi0_loc' 0.725 "{&pi}{subscript:0} = 0`pi0'")	
				
			*Plot histogram of p-vals with lambda, pi0, and p-val cutoff for FDR control
			local pi0=round(`pi_0',0.01)
			local lambda_loc=`lambda'+0.06
			local pi0_loc=`pi0'+0.2
			local p_loc=`p_cutoff'+0.175
			graph twoway hist p_f, color(navy) fintensity(0) width(0.05) /// 
					|| function y=`pi_0', range(0 1) lcolor(maroon) lpattern(dash) ///
					xline(`lambda', lcolor(black) lpattern(dash)) ///
					xline(`p_cutoff', lcolor(black) lpattern(dash)) ///
					legend(off) ///
					xtitle("P-value from test of no disc. against Black applicants") ytitle("Density") ///
					scheme(s1color) ///
					text(3.5 `lambda_loc' "{&lambda} = 0`lambda'") ///
					text(`pi0_loc' 0.725 "{&pi}{subscript:0} = 0`pi0'") ///
					text(5.5 `p_loc' "`N_firms' firms with q-vals < 0`FDR'")

		}
	
	
