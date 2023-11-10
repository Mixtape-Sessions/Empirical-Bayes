

*******************************************
**** Chris Walters
**** 10/2023
**** This program estimates the model parameters using empirical Bayes methods
********************************************

	********************
	*** Basic setup ****
	********************
	
		clear
		cap set more off
		cap set trace off
		cap log close
		set seed 1028
		/* cd "Labs/Lab-1/" */
	
		****Load data
		insheet using "vam_example_data.csv", names comma clear
		ren y Y
		ren x X
		ren d D
		ren theta_d theta_D

    *Number of schools
    local J=50
	
		****Generate school attendance dummies and sample sizes
			sum D, detail
			local J=r(max)
			foreach j of numlist 1/`J' {
				gen D_`j'=(D==`j')
				count if D==`j'
				gen N`j'=r(N)
			}
	
		******Value-added models
			
			***Fit VAM regressions
			
				**No controls
				
					*Run regression
					reg Y D_*, r nocons
					
					*Collect estimates and standard errors
					foreach j of numlist 1/`J' {
						gen thetahat_unc`j'=_b[D_`j']
						gen s_unc`j'=_se[D_`j']
					}
				
				**With control for X

					*Run regression
					reg Y D_* X, r nocons
					
					*Collect estimates and standard errors
					foreach j of numlist 1/`J' {
						gen thetahat_con`j'=_b[D_`j']
						gen s_con`j'=_se[D_`j']
					}
					

			***Collapse to school-level data set with true and estimated VA
			foreach j of numlist 1/`J' {
				sum theta_D if D==`j'
				gen theta`j'=r(mean)
				
			}
			drop theta_D
			keep theta* s_* N*
			duplicates drop
			gen i=1
			reshape long theta thetahat_unc s_unc thetahat_con s_con N, i(i) j(j)
			drop i

		
		******Estimate hyperparameters for each model
		foreach m in unc con {
			
			*Mean
			sum thetahat_`m'
			local muhat_`m'=r(mean)
			
			*Standard dev
			gen c=(thetahat_`m'-(`muhat_`m''))^2 - s_`m'^2
			sum c
			local sigmahat_`m'=sqrt(r(mean))
			drop c
				
			*Display results
			disp "Estimated hyperparameters for model `m': mean = `muhat_`m'', SD = `sigmahat_`m''"
			
		}
			
		****Compute linear shrinkage estimates for the controlled model
		foreach m in unc con {
		
			*Shrinkage factor
			gen lambda_`m'=((`sigmahat_`m'')^2)/(((`sigmahat_`m'')^2)+(s_`m'^2))
			
			*Linear shrinkage formula
			gen thetastar_`m'=lambda_`m'*thetahat_`m'+(1-lambda_`m')*(`muhat_`m'')
		
		}
	
	
	***************************************
	*** Assess performance of estimators **
	***************************************
	
		*****Compare dispersion of unbiased estimates, estimated prior, shrunk posteriors, and true parameters
		foreach m in unc con {
		
			*Unbiased estimates
			sum thetahat_`m'
			local SD_1=r(sd)
			
			*Shrunk posteriors
			sum thetastar_`m'
			local SD_2=r(sd)
			
			*True parameters
			sum theta
			local SD_3=r(sd)
			
			*Report comparison
			disp "Model `m': SD estimates=`SD_1', SD prior =`sigmahat_`m'', SD posteriors=`SD_2', SD truth = `SD_3'"

		
		}
		
		****Compare MSE of each estimator
		foreach e in thetahat_unc thetastar_unc thetahat_con thetastar_con {
		
			gen MSE_`e'=(`e'-theta)^2
		
		}
		gen MSE_diff = MSE_thetahat_con - MSE_thetastar_con
		sum MSE*, detail
		
		
		****For controlled model, graph true VA against unbiased estimates and linear shrinkage estimates, along w/45-deg line
		graph twoway scatter theta thetahat_con, mcolor(navy) ///
			|| scatter theta thetastar_con, mcolor(maroon) ///
			|| function y=x, range(thetahat_con) lcolor(black) lpattern(dash) ///
			xtitle("Estimated value-added") ytitle("True value-added") ///
			legend(lab(1 "Unbiased estimates") lab(2 "Linear shrinkage estimates")) ///
			scheme(s1color)
			
		****Histograms of unbiased estimates, prior, posteriors
		local N=_N
		local width=0.06
		local bwidth=`width'/1.5
    
    * For the figure, the true sd of theta
    local sigma_theta = 0.2
		local min = -4*`sigma_theta'
		local max = 4*`sigma_theta'
    
		graph twoway hist thetahat_con, freq width(`width') color(navy) fintensity(0) ///
			|| hist thetastar_con, freq width(`width') barwidth(`bwidth') color(maroon) ///
			|| function y = (`N'*`width')*normalden(x,`muhat_con',`sigmahat_con'), range(`min' `max') lcolor(black) ///
			xtitle("Value-added") ytitle("Schools (frequency)") ///
			legend(lab(1 "Unbiased estimates") lab(2 "Linear shrinkage estimates") lab(3 "Prior distribution") order(1 3 2)) ///
			scheme(s1color)
	
	
	
	
	
	
	
	
	
	
	
