*******************************************
**** Chris Walters
**** 10/2023
**** This program simulates data from a school value-added model
********************************************


/* Parameter info and reasonable defaults: 
**Sample dimensions
  *Number of students
  local n_students=2500
  
  *Number of schools
  local n_schools=50

**Distribution of school-level parameters
  *Parametric form for VA distribution (normal by default, or lognormal)
  local va_model="normal"

  *Mean and SD of school value-added
  local mu_theta=0
  local sigma_theta=0.2
  
  *Dispersion of school mean utilities
  local sigma_delta=1
  
  *Dispersion of school utility coefficients on student covariate X
  local sigma_gamma=0.5
  
**Student-level parameters
  *Standard dev of X
  local sigma_x=1

  *Coefficient on X in outcome equation
  local beta=1
  
  *Residual variance of outcome
  local sigma_y=1 
*/

cap program drop simulate_data
program define simulate_data
  version 13
  syntax , n_students(integer) n_schools(integer) sigma_y(real) sigma_x(real) beta(real) va_model(string) mu_theta(real) sigma_theta(real) sigma_delta(real) sigma_gamma(real) 

  
  ******Draw school-level parameters
  
    *Generate data set of schools
    qui set obs `n_schools'
    qui gen j =_n
  
    *Value-added
    qui gen theta=`mu_theta'+`sigma_theta'*invnormal(uniform())
    if "`va_model'"=="lognormal" {
      qui replace theta=exp(theta)
    }
    
    *Utility parameters
    qui gen delta=`sigma_delta'*invnormal(uniform())
    qui gen gamma=`sigma_gamma'*invnormal(uniform())
    
  ****Generate student-level data
  
    *Reshape wide
    qui gen i=1
    qui reshape wide theta delta gamma, i(i) j(j)
    drop i
    
    *Create student-level obs
    qui expand `n_students'
    gen student_id=_n
    
    *Generate student covariate 
    qui gen X=`sigma_x'*invnormal(uniform())
    
    *Code school attended by each student
    qui gen school_id=.
    qui gen Umax=.
    qui gen theta_D=.
    foreach j of numlist 1/`n_schools' {
      qui gen U`j'=delta`j'+gamma`j'*X - log(-log(uniform()))
      qui replace school_id=`j' if U`j'>Umax | Umax==.
      qui replace Umax=U`j' if school_id==`j'
      qui replace theta_D=theta`j' if school_id==`j'
    }
    
    *Generate student outcome
    qui gen Y=theta_D + `beta'*X + `sigma_y'*invnormal(uniform())
    
    *Generate school attendance dummies
    foreach j of numlist 1/`n_schools' {
      qui gen D_`j'=(school_id==`j')
      label var D_`j' "Dummy if student is in school `j'"
      qui count if school_id==`j'
      qui gen N`j'=r(N)
    }

  ****Clean-up output
  order school_id student_id X theta_D 
  keep  student_id school_id Y X theta* D_* N*
end 

simulate_data, ///
  /// Number of students and schools
  n_students(`n_students') n_schools(`n_schools') ///
  /// standard deviation of outcome residual
  sigma_y(`sigma_y') ///
  /// Distribution of true value-added theta
  va_model("normal") mu_theta(0) sigma_theta(`sigma_theta') ///
  /// Dispersion of school mean utilities
  sigma_delta(1) ///
  /// Dispersion of school utility coefficients on student covariate X
  sigma_gamma(0.5) ///
  /// Covariate X
  sigma_x(1) beta(1)

outsheet using "vam_example_data.csv", names comma replace


