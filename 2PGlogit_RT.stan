functions{
real h_alpha(real alpha1, real alpha2, real eta){
if (eta > 0){
	if (alpha1 > 0)
		return (exp(alpha1 * eta) - 1) / alpha1;
	else if (alpha1 == 0)
		return eta;
	else //alpha1 < 0
		return -log(1 - alpha1 * eta) / alpha1;
}
else{ 	if (alpha2 > 0)
		return -(exp(-alpha2 * eta) - 1) / alpha2;
	else if (alpha2 == 0)
		return eta;
	else //alpha2 < 0
		return log(1 + alpha2 * eta) / alpha2;
}
}
}

data{
int<lower=0> n_student;
int<lower=0> n_item;
int<lower=0,upper=1> Y[n_student,n_item];
real logT[n_student,n_item];
}

parameters{
  row_vector[2] theta_tao[n_student];
  real<lower=0> sigma_tao;
  cholesky_factor_corr[2] L_Omega_person;
  vector<lower=0,upper=3>[n_item] a;
  vector[n_item] b;
  vector<lower=-1>[n_item] alpha1;
  vector<lower=-1>[n_item] alpha2;
  vector[n_item] lambda;
  vector<lower=0>[n_item] sigma_tj;
}

transformed parameters{
  vector[n_student] theta;
  vector[n_student] tao;
  corr_matrix[2] Omega_person;
  cov_matrix[2] Sigma_person;


  Omega_person = multiply_lower_tri_self_transpose(L_Omega_person);
  Sigma_person = quad_form_diag(Omega_person,[1,sigma_tao]);
  theta = to_vector(theta_tao[,1]);
  tao = to_vector(theta_tao[,2]);
}

model{
 sigma_tao ~ cauchy(0,5);
 L_Omega_person ~ lkj_corr_cholesky( 1);
 theta_tao ~ multi_normal( [0,0], Sigma_person);
 a ~ uniform(0,3);
 b ~  normal(0,10);
 lambda ~  normal(0,10);
 alpha1 ~ normal(0,0.5); 
 alpha2 ~ normal(0,0.5); 
 sigma_tj ~ cauchy(0,5);
  for(i in 1:n_student){
	  for (j in 1:n_item){
        Y[i,j] ~ bernoulli_logit( h_alpha(alpha1[j], alpha2[j], a[j]*(theta[i] - b[j])));
        logT[i,j] ~  normal(  lambda[j] - tao[i] , sigma_tj[j]);
	}
 }
}

generated quantities{
  vector[n_item] log_lik_Y[n_student];
  vector[n_item] log_lik_logT[n_student];

for (i in 1: n_student){
	for (j in 1: n_item){
	  log_lik_Y[i,j] = bernoulli_logit_lpmf(Y[i,j] | h_alpha(alpha1[j], alpha2[j], a[j]*(theta[i] - b[j])));
        log_lik_logT[i,j] = normal_lpdf(logT[i,j] | lambda[j] - tao[i] , sigma_tj[j]);
	  }
   }
}