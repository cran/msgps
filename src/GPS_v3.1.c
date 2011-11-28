//Ｒ用 Ｃ
//#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <R.h>
#include <Rdefines.h>
#include <R_ext/Parse.h>
#include <R_ext/BLAS.h>


//////////////////////////////////////////////////
//BLAS: y := alpha*A*x + beta*y
long int gl_incx = 1;
long int gl_incy = 1;
double alphaOne = 1.0;
double alphaminusOne = -1.0;
double betaOne = 1.0;
double betaZero = 0.0;
double Zero = 0.0;
char TRANSN = 'N';
char TRANST = 'T';
char UPPER = 'U';


void _dgemv(long int m, long int n, double *A, long int lda, double *x, double *ans){
	dgemv_(&TRANSN, &m, &n, &alphaOne, A, &lda, x, &gl_incx, &betaOne, ans, &gl_incy);
}

SEXP Rdgemv01(SEXP ex_A, SEXP ex_x, SEXP ex_y)
{
	SEXP ans;
	int m = INTEGER(GET_DIM(ex_A))[0];
	int n = INTEGER(GET_DIM(ex_A))[1];
	double *A;
	double *x;
	double *y;
	int lda = m;
	long int i;
	
	A = REAL(ex_A);
	x = REAL(ex_x);
	y = REAL(ex_y);
	PROTECT(ans = allocVector(REALSXP, m));
	for( i=0; i<m; i++){
		REAL(ans)[i] = y[i];
	}
	//dgemv_(&TRANSN, &m, &n, &alphaOne, A, &lda, x, &gl_incx, &betaOne, REAL(ans), &gl_incy);
	_dgemv(m, n, A, lda, x, REAL(ans));
	UNPROTECT(1);
	return(ans);
}



// GPS
// 入力：データ(y,X)，step幅（delta_t），ステップ数（STEP），
// 出力：係数の推定値行列(betahat_matrix)，自由度ベクトル（df_zou）
SEXP gps(SEXP ex_X, SEXP ex_y, SEXP ex_delta_t, SEXP ex_penalty_type, SEXP ex_para, SEXP ex_STEP, SEXP ex_pmax, SEXP ex_weight,SEXP ex_standardize_vec)
{
	long int p = INTEGER(GET_DIM(ex_X))[1];
	long int N = INTEGER(GET_DIM(ex_X))[0];
	long int NN =N*N;
	long int STEP = INTEGER(ex_STEP)[0];
	long int i,i1,j, step0;
	int one=1;
	int flag_St;
	int flag_XTxjstar=0;
	long int jstar=0;
	long int jstar0=0;
	//long int beta_vecminus2_index;
	long int beta_vecminus2_flag;
	long int sign_gt_jstar_int;
	//long int STEP_adj_int;
	long int jsort;
	long int pmax=INTEGER(ex_pmax)[0];
	long int jstar_sort=0;
	
	int penalty_type_int = INTEGER(ex_penalty_type)[0];
	int check_int = 20;
	long int check2_int = 0;
	
	
	double alphaN2 = 2.0 / (double) N;
	double alphabetaminusN2;
	double alpha_minusdelta_t;
	//double jstar_max;
	double sign_gt_jstar;
	double temp;
	double beta_jstar_previous=0.0;
	//double increment_covpenalty;
	//double increment_covpenalty_minus;
	//double zeros_p_vec[p];
	
	
	//double *RSS;             
	
	SEXP betahat_index_vec;
	SEXP betahat_index_vec_adj;
	SEXP beta_vec;
	SEXP g_t;
	SEXP absg_t;          
	SEXP g_t0;
	SEXP P_t;
	SEXP Lambda_t;
	SEXP absLambda_t;          
	SEXP S_t;
	SEXP x_jstar;
	SEXP STEP_adj;
	SEXP RSS_vec;
	SEXP increment_covpenalty_vec_power;
	SEXP selected_variable_index_vec;
	SEXP selected_variable_index_sort_vec;
	SEXP XTXNA_adj;
	SEXP XTxjstar;
	SEXP ans;
	SEXP penalty_type_str;  
	SEXP tuning;  
	SEXP tuning_stand;  
	SEXP sum_lambda_t;
	
	
	
	
	
	PROTECT(betahat_index_vec = allocVector(INTSXP, STEP+1));
	PROTECT(betahat_index_vec_adj = allocVector(INTSXP, STEP+1));
	PROTECT(beta_vec = allocVector(REALSXP,p));
	PROTECT(g_t = allocVector(REALSXP,p));
	PROTECT(absg_t = allocVector(REALSXP,p));
	PROTECT(g_t0 = allocVector(REALSXP,N));
	PROTECT(P_t = allocVector(REALSXP,p));
	PROTECT(Lambda_t = allocVector(REALSXP,p));
	PROTECT(absLambda_t = allocVector(REALSXP,p));
	PROTECT(S_t = allocVector(INTSXP,p));
	PROTECT(x_jstar = allocVector(REALSXP,N));
	PROTECT(RSS_vec = allocVector(REALSXP, STEP+1));
	PROTECT(STEP_adj = allocVector(INTSXP,1));
	PROTECT(increment_covpenalty_vec_power = allocVector(REALSXP,STEP+1));
	PROTECT(selected_variable_index_vec = allocVector(INTSXP,p));
	PROTECT(selected_variable_index_sort_vec = allocVector(INTSXP,p));
	PROTECT(XTXNA_adj = allocMatrix(REALSXP,p,pmax));
	PROTECT(XTxjstar = allocVector(REALSXP,p));
	PROTECT(sum_lambda_t = allocVector(REALSXP,STEP+1));
	PROTECT(tuning = allocVector(REALSXP,STEP+1));
	PROTECT(tuning_stand = allocVector(REALSXP,STEP+1));
	PROTECT(ans = allocVector(VECSXP, 8)); 


	//初期化
	
	
	for( i=0; i<p; i++){
		REAL(P_t)[i] = 1.0; //lassoとして設定
		REAL(beta_vec)[i] = 0.0;
		REAL(g_t)[i] = 0.0;
		REAL(Lambda_t)[i] = 0.0;
		REAL(absLambda_t)[i] = 0.0;
		REAL(absg_t)[i] = 0.0;
		REAL(XTxjstar)[i] = 0.0;
		INTEGER(S_t)[i] = 0;
		INTEGER(selected_variable_index_vec)[i] =0;
		INTEGER(selected_variable_index_sort_vec)[i]=0;
		
		
	}


		//p_tの初期値を与える
		if(penalty_type_int == 0){ //enet
			for( i=0; i<p; i++){
				REAL(P_t)[i] = REAL(ex_para)[0] * fabs(REAL(beta_vec)[i]) + (1.0 - REAL(ex_para)[0]); 
			}
		}
		
		if(penalty_type_int == 1){ //genet
			for( i=0; i<p; i++){
				REAL(P_t)[i] = (1.0 - REAL(ex_para)[0]) / ( (1.0 - REAL(ex_para)[0])* fabs(REAL(beta_vec)[i]) +  REAL(ex_para)[0] ); 
			}
		}
		
		if(penalty_type_int == 2){ //alasso
			for( i=0; i<p; i++){
				REAL(P_t)[i] = REAL(ex_weight)[i]; 
			}
		}		
		
	
	
	for( i=0; i<STEP+1; i++){
		REAL(RSS_vec)[i] = 0.0;
		INTEGER(betahat_index_vec)[i] = 0;
		REAL(increment_covpenalty_vec_power)[i]=0.0;
		INTEGER(betahat_index_vec_adj)[i] = 0;          
		REAL(sum_lambda_t)[i]=0.0;
		REAL(tuning)[i]=0.0;
		REAL(tuning_stand)[i]=0.0;
	}          
	
	//generalized elastic netのtuningの初期値
	if(penalty_type_int == 1){ //genet
		REAL(tuning)[0] = -log(REAL(ex_para)[0]);
		REAL(tuning_stand)[0] = -log(REAL(ex_para)[0]);
	}


	
	for( i=0; i<N; i++){
		REAL(x_jstar)[i] = 0.0;
		REAL(g_t0)[i] = REAL(ex_y)[i];
	}
	
	
	for( i=0; i<p*pmax; i++){
		REAL(XTXNA_adj)[i] = 0.0;
	}
	
	//RSSの計算
	REAL(RSS_vec)[0] = dnrm2_(&N,REAL(ex_y),&one);
	REAL(RSS_vec)[0] = REAL(RSS_vec)[0]*REAL(RSS_vec)[0];
	
	
	//g_tの定義
	F77_CALL(dgemv)(&TRANST, &N, &p, &alphaN2, REAL(ex_X), &N, REAL(ex_y), &one, &betaZero, REAL(g_t), &one);
	
	//jstarの初期値
	jstar=0;
	sign_gt_jstar=1.0;
	
	//jsortの初期値
	jsort=0;
	
	
	////////////////////
	////////GPS//////
	//////////////////
	for( step0=1; step0<STEP+1; step0++){
		
		//g_t <- 2/N * t(X) %*% (y - X %*% beta_t )を計算
		//daxpy_(&TRANSN, &N, &p, &alphaminusOne, REAL(ex_X), &N, REAL(beta_vec), &one, &betaOne, REAL(g_t0), &one);
		//dgemv_(&TRANSN, &N, &p, &alphaminusOne, REAL(ex_X), &N, REAL(beta_vec), &one, &betaOne, REAL(g_t0), &one);

		//p_tの計算
		if(penalty_type_int == 0 && REAL(ex_para)[0] != 0.0){ //enet
			REAL(P_t)[jstar] = REAL(ex_para)[0] * fabs(REAL(beta_vec)[jstar]) + (1.0 - REAL(ex_para)[0]); 
		}
		
		if(penalty_type_int == 1){ //genet
			REAL(P_t)[jstar] = (1.0 - REAL(ex_para)[0]) / ( (1.0 - REAL(ex_para)[0])* fabs(REAL(beta_vec)[jstar]) +  REAL(ex_para)[0] ); 
		}
		
		//if(penalty_type_int == 2){ //alasso
			//for( i=0; i<p; i++){
			//	REAL(P_t)[i] = REAL(ex_weight)[i]; 
			//}
		//}		
		
		
		alphabetaminusN2 = -1.0 * alphaN2 * REAL(ex_delta_t)[0] * sign_gt_jstar;
		if(flag_XTxjstar==1){
			F77_CALL(dgemv)(&TRANST, &N, &p, &alphaOne, REAL(ex_X), &N, REAL(x_jstar), &one, &betaZero, REAL(XTxjstar), &one);
			for( i=0; i<p; i++){			
				REAL(XTXNA_adj)[i+jstar_sort*p] = REAL(XTxjstar)[i];
			}
		}else{
			for( i=0; i<p; i++){
				REAL(XTxjstar)[i] = REAL(XTXNA_adj)[i+jstar_sort*p];
			}
		}
		
		F77_CALL(daxpy)(&p, &alphabetaminusN2, REAL(XTxjstar), &one, REAL(g_t), &one);		
		//dgemv_(&TRANST, &N, &p, &alphabetaminusN2, REAL(ex_X), &N, REAL(x_jstar), &one, &betaOne, REAL(g_t), &one);
		
		
		

		flag_St =0;		  
		for( i=0; i<p; i++){
		//	REAL(absg_t)[i] = fabs(REAL(g_t)[i]);
			REAL(Lambda_t)[i] = REAL(g_t)[i] / REAL(P_t)[i];
			REAL(absLambda_t)[i] = fabs(REAL(Lambda_t)[i]);
			REAL(absg_t)[i] = fabs(REAL(g_t)[i]);
			if(REAL(Lambda_t)[i]*REAL(beta_vec)[i]< 0  && REAL(absg_t)[i] > 2.0 * REAL(ex_delta_t)[0] / (double) N ){
			    REAL(absLambda_t)[i] = REAL(absLambda_t)[i] + 1000000000000.0;
		        flag_St =1;		  
			}
			if(REAL(absg_t)[i] < 2.0 * REAL(ex_delta_t)[0] / (double) N ){
				REAL(absLambda_t)[i] = 0.0;
				flag_St =1;		  
			}
		}
		
		
		jstar = F77_CALL(idamax)(&p,REAL(absLambda_t),&one);
		jstar0 = F77_CALL(idamax)(&p,REAL(Lambda_t),&one);
		jstar = jstar-1;
		jstar0 = jstar0-1;
		
		
		
		flag_XTxjstar=0;
		if(INTEGER(selected_variable_index_vec)[jstar]==0){
			INTEGER(selected_variable_index_vec)[jstar] = jsort+1;
			INTEGER(selected_variable_index_sort_vec)[jsort] = jstar+1;
			jsort=jsort+1;
			flag_XTxjstar=1;
		}
		
		jstar_sort=INTEGER(selected_variable_index_vec)[jstar]-1;
		
		
		INTEGER(betahat_index_vec_adj)[step0] = INTEGER(selected_variable_index_vec)[jstar];
		
		//beta_t[jstar] =   beta_t[jstar] + delta_t * sign(g_t[jstar])
		if(REAL(Lambda_t)[jstar]>0){
			sign_gt_jstar=1.0;
			sign_gt_jstar_int = 1;
		}
		else {
			sign_gt_jstar=-1.0;
			sign_gt_jstar_int = -1;
		}
		beta_jstar_previous = REAL(beta_vec)[jstar];
		REAL(beta_vec)[jstar] = REAL(beta_vec)[jstar] + REAL(ex_delta_t)[0] * sign_gt_jstar;
		
        INTEGER(betahat_index_vec)[step0] = (jstar+1) * sign_gt_jstar_int;
		
		//cov_penaltyのための係数計算
		REAL(increment_covpenalty_vec_power)[step0] = 1.0 /  REAL(absg_t)[jstar];
		
		
		//g_t0をy-Xbetaと定義
		//g_t0=g_t0 - delta_beta*x[jstar]*sign
		for( i=0; i<N; i++){
			REAL(x_jstar)[i] = REAL(ex_X)[i +jstar*N];
		}
		alpha_minusdelta_t = -1.0*REAL(ex_delta_t)[0]* sign_gt_jstar;
		F77_CALL(daxpy)(&N,&alpha_minusdelta_t,REAL(x_jstar),&one,REAL(g_t0), &one);
		REAL(RSS_vec)[step0] = dnrm2_(&N,REAL(g_t0),&one);
		REAL(RSS_vec)[step0] = REAL(RSS_vec)[step0]*REAL(RSS_vec)[step0];		
		
		//tuning parameterの計算
		if(penalty_type_int == 0){ //enet
			temp = REAL(ex_para)[0] * (REAL(beta_vec)[jstar]*REAL(beta_vec)[jstar] - beta_jstar_previous*beta_jstar_previous) / 2.0 + (1.0 - REAL(ex_para)[0]) * ( fabs(REAL(beta_vec)[jstar]) - fabs(beta_jstar_previous) ); 
		}
		
		if(penalty_type_int == 1){ //genet
			temp = log( (1.0 - REAL(ex_para)[0])* fabs(REAL(beta_vec)[jstar]) +  REAL(ex_para)[0] ) - log( (1.0 - REAL(ex_para)[0])* fabs(beta_jstar_previous) +  REAL(ex_para)[0] ); 
		}
		
		if(penalty_type_int == 2){ //alasso
			temp = REAL(ex_weight)[jstar]*fabs(REAL(beta_vec)[jstar]) -  REAL(ex_weight)[jstar]*fabs(beta_jstar_previous); 
		}		
		
		REAL(tuning)[step0] = REAL(tuning)[step0-1] + temp; 
		REAL(tuning_stand)[step0] = REAL(tuning_stand)[step0-1] + REAL(ex_standardize_vec)[jstar] * temp; 

		//dfが妥当かcheck!
		beta_vecminus2_flag = 0;
		//beta_vecminus2_index = step0-1;
		//if(fabs(REAL(df_vec)[step0] - REAL(df_vec)[step0-1]) > REAL(ex_delta_t)[0] * 100.0) beta_vecminus2_flag=1;
		//収束してるのかcheck!
		check2_int = 0;
		for( i=0; i<p; i++){
			if(REAL(absg_t)[i] < 2.0 * REAL(ex_delta_t)[0] / (double) N) check2_int=check2_int+1;
		}
		if(check2_int==p) beta_vecminus2_flag=1;
		
		//収束してるのかcheck! ver2
		REAL(sum_lambda_t)[step0] = F77_CALL(dasum)(&p,REAL(Lambda_t),&one);
		//if(REAL(sum_lambda_t)[step0] < REAL(tol)[1]) beta_vecminus2_flag=1;
		//さらにちぇっっく
		if(step0>check_int){
			for( i=step0-check_int; i<step0; i++){
				if(fabs(REAL(sum_lambda_t)[step0] - REAL(sum_lambda_t)[i]) == 0.0 ){
					beta_vecminus2_flag=1;
				}
			}
		}
		//pmaxに達したのかcheck
		if(jsort==pmax+1) beta_vecminus2_flag=1;
		if(beta_vecminus2_flag==1) break;
		
		
	}	
	
	
	
	
	INTEGER(STEP_adj)[0] = step0;  
	
	
	
	SET_VECTOR_ELT(ans, 0, betahat_index_vec);
	SET_VECTOR_ELT(ans, 1, STEP_adj);	
	SET_VECTOR_ELT(ans, 2, RSS_vec);	
	SET_VECTOR_ELT(ans, 3, increment_covpenalty_vec_power);
	SET_VECTOR_ELT(ans, 4, selected_variable_index_sort_vec);
	SET_VECTOR_ELT(ans, 5, betahat_index_vec_adj);
	SET_VECTOR_ELT(ans, 6, tuning);
	SET_VECTOR_ELT(ans, 7, tuning_stand);
	UNPROTECT(22);
	
	return(ans);
	
	
}





// CVGPS
// 入力：データ(y,X)，step幅（delta_t），ステップ数（STEP），
// 出力：係数の推定値行列(betahat_matrix)，自由度ベクトル（df_zou）
SEXP cvgps(SEXP ex_X, SEXP ex_y, SEXP ex_delta_t, SEXP ex_penalty_type, SEXP ex_para, SEXP ex_STEP,  SEXP ex_pmax, SEXP ex_weight, SEXP ex_X_test, SEXP ex_y_test, SEXP ex_standardize_vec)
{
	long int p = INTEGER(GET_DIM(ex_X))[1];
	long int N = INTEGER(GET_DIM(ex_X))[0];
	long int Ntest = INTEGER(GET_DIM(ex_X_test))[0];
	long int NN =N*N;
	long int NNtest = Ntest*Ntest;
	long int STEP = INTEGER(ex_STEP)[0];
	long int i,j, step0;
	int one=1;
	int flag_St;
	int flag_XTxjstar=0;
	long int jstar=0;
	long int jstar0=0;
	//long int beta_vecminus2_index;
	long int beta_vecminus2_flag;
	long int sign_gt_jstar_int;
	//long int STEP_adj_int;
	long int jsort;
	long int pmax=INTEGER(ex_pmax)[0];
	long int jstar_sort=0;

	int penalty_type_int = INTEGER(ex_penalty_type)[0];
	int check_int = 20;
	long int check2_int=0;
	
	double alphaN2 = 2.0 / (double) N;
	double alphabetaminusN2;
	double alpha_minusdelta_t;
	//double jstar_max;
	double sign_gt_jstar;
	double temp;
	//double increment_covpenalty;
	//double increment_covpenalty_minus;
	//double zeros_p_vec[p];
	//double sum_lambda_t=0.0;
	
	
	//double *RSS;             
	
	SEXP betahat_index_vec;
	SEXP betahat_index_vec_adj;
	SEXP beta_vec;
	SEXP g_t;
	SEXP absg_t;          
	SEXP g_t0;
	SEXP g_t0_test;
	SEXP P_t;
	SEXP Lambda_t;
	SEXP absLambda_t;          
	SEXP S_t;
	SEXP x_jstar;
	SEXP x_jstar_test;
	SEXP STEP_adj;
	SEXP RSS_vec;
	SEXP RSS_vec_test;
	SEXP increment_covpenalty_vec;
	SEXP selected_variable_index_vec;
	SEXP selected_variable_index_sort_vec;
	SEXP XTXNA_adj;
	SEXP XTxjstar;
	SEXP sum_lambda_t;
	SEXP tuning;  
	SEXP tuning_stand;  
	SEXP tuning_sum;  
	SEXP tuning_sum_stand;  
	SEXP ans;
	
	
	
	
	
	PROTECT(betahat_index_vec = allocVector(INTSXP, STEP+1));
	PROTECT(betahat_index_vec_adj = allocVector(INTSXP, STEP+1));
	PROTECT(beta_vec = allocVector(REALSXP,p));
	PROTECT(g_t = allocVector(REALSXP,p));
	PROTECT(absg_t = allocVector(REALSXP,p));
	PROTECT(g_t0 = allocVector(REALSXP,N));
	PROTECT(g_t0_test = allocVector(REALSXP,Ntest));
	PROTECT(P_t = allocVector(REALSXP,p));
	PROTECT(Lambda_t = allocVector(REALSXP,p));
	PROTECT(absLambda_t = allocVector(REALSXP,p));
	PROTECT(S_t = allocVector(INTSXP,p));
	PROTECT(x_jstar = allocVector(REALSXP,N));
	PROTECT(x_jstar_test = allocVector(REALSXP,Ntest));
	PROTECT(RSS_vec = allocVector(REALSXP, STEP+1));
	PROTECT(RSS_vec_test = allocVector(REALSXP, STEP+1));
	PROTECT(STEP_adj = allocVector(INTSXP,1));
	PROTECT(increment_covpenalty_vec = allocVector(REALSXP,STEP+1));
	PROTECT(selected_variable_index_vec = allocVector(INTSXP,p));
	PROTECT(selected_variable_index_sort_vec = allocVector(INTSXP,p));
	PROTECT(XTXNA_adj = allocMatrix(REALSXP,p,pmax));
	PROTECT(XTxjstar = allocVector(REALSXP,p));
	PROTECT(sum_lambda_t = allocVector(REALSXP,STEP+1));
	PROTECT(tuning = allocVector(REALSXP,STEP+1));
	PROTECT(tuning_stand = allocVector(REALSXP,STEP+1));
	PROTECT(tuning_sum = allocVector(REALSXP,STEP+1));
	PROTECT(tuning_sum_stand = allocVector(REALSXP,STEP+1));
	PROTECT(ans = allocVector(VECSXP, 7));
	
	
	
	
	//初期化
	
	
	for( i=0; i<p; i++){
		REAL(P_t)[i] = 1.0;  //lassoとして設定
		REAL(beta_vec)[i] = 0.0;
		REAL(g_t)[i] = 0.0;
		REAL(Lambda_t)[i] = 0.0;
		REAL(absLambda_t)[i] = 0.0;
		REAL(absg_t)[i] = 0.0;
		REAL(XTxjstar)[i] = 0.0;
		INTEGER(S_t)[i] = 0;
		INTEGER(selected_variable_index_vec)[i] =0;
		INTEGER(selected_variable_index_sort_vec)[i]=0;
		
		
	}
	
//alassoのP_tをあらかじめ計算（計算時間短縮）
	if(penalty_type_int == 2){ //alasso
		for( i=0; i<p; i++){
			REAL(P_t)[i] = REAL(ex_weight)[i]; 
		}
	}		

	
	
	for( i=0; i<STEP+1; i++){
		REAL(RSS_vec)[i] = 0.0;
		REAL(RSS_vec_test)[i] = 0.0;
		INTEGER(betahat_index_vec)[i] = 0;
		REAL(increment_covpenalty_vec)[i]=0.0;
		REAL(sum_lambda_t)[i]=0.0;
		INTEGER(betahat_index_vec_adj)[i] = 0;          
		REAL(tuning)[i]=0.0;
		REAL(tuning_stand)[i]=0.0;
		REAL(tuning_sum)[i]=0.0;
		REAL(tuning_sum_stand)[i]=0.0;
	}          
	
	for( i=0; i<N; i++){
		REAL(x_jstar)[i] = 0.0;
		REAL(g_t0)[i] = REAL(ex_y)[i];
	}
	
	
	for( i=0; i<Ntest; i++){
		REAL(g_t0_test)[i] = REAL(ex_y_test)[i];
		REAL(x_jstar_test)[i] = 0.0;
	}
	
	
	for( i=0; i<p*pmax; i++){
		REAL(XTXNA_adj)[i] = 0.0;
	}
	
	//RSSの計算
	REAL(RSS_vec)[0] = dnrm2_(&N,REAL(ex_y),&one);
	REAL(RSS_vec)[0] = REAL(RSS_vec)[0]*REAL(RSS_vec)[0];
	REAL(RSS_vec_test)[0] = dnrm2_(&Ntest,REAL(ex_y_test),&one);
	REAL(RSS_vec_test)[0] = REAL(RSS_vec_test)[0]*REAL(RSS_vec_test)[0];
	
	//g_tの定義
	F77_CALL(dgemv)(&TRANST, &N, &p, &alphaN2, REAL(ex_X), &N, REAL(ex_y), &one, &betaZero, REAL(g_t), &one);
	
	//jstarの初期値
	jstar=0;
	sign_gt_jstar=1.0;
	
	//jsortの初期値
	jsort=0;
	
	
	////////////////////
	////////GPS//////
	//////////////////
	for( step0=1; step0<STEP+1; step0++){
		
		//p_tの計算
		if(penalty_type_int == 0 && REAL(ex_para)[0] != 0.0){ //enet
			for( i=0; i<p; i++){
				REAL(P_t)[i] = REAL(ex_para)[0] * fabs(REAL(beta_vec)[i]) + (1.0 - REAL(ex_para)[0]); 
			}
		}
		
		if(penalty_type_int == 1){ //genet
			for( i=0; i<p; i++){
				REAL(P_t)[i] = (1.0 - REAL(ex_para)[0]) / ( (1.0 - REAL(ex_para)[0])* fabs(REAL(beta_vec)[i]) +  REAL(ex_para)[0] ); 
			}
		}
		
		if(penalty_type_int == 2){ //alasso
			//for( i=0; i<p; i++){
			//	REAL(P_t)[i] = REAL(ex_weight)[i]; 
			//}
		}		


		
		//g_t <- 2/N * t(X) %*% (y - X %*% beta_t )を計算
		//daxpy_(&TRANSN, &N, &p, &alphaminusOne, REAL(ex_X), &N, REAL(beta_vec), &one, &betaOne, REAL(g_t0), &one);
		//dgemv_(&TRANSN, &N, &p, &alphaminusOne, REAL(ex_X), &N, REAL(beta_vec), &one, &betaOne, REAL(g_t0), &one);
		alphabetaminusN2 = -1.0 * alphaN2 * REAL(ex_delta_t)[0] * sign_gt_jstar;
		if(flag_XTxjstar==1){
			F77_CALL(dgemv)(&TRANST, &N, &p, &alphaOne, REAL(ex_X), &N, REAL(x_jstar), &one, &betaZero, REAL(XTxjstar), &one);
			for( i=0; i<p; i++){			
				REAL(XTXNA_adj)[i+jstar_sort*p] = REAL(XTxjstar)[i];
			}
		}else{
			for( i=0; i<p; i++){
				REAL(XTxjstar)[i] = REAL(XTXNA_adj)[i+jstar_sort*p];
			}
		}
		
		F77_CALL(daxpy)(&p, &alphabetaminusN2, REAL(XTxjstar), &one, REAL(g_t), &one);		
		//dgemv_(&TRANST, &N, &p, &alphabetaminusN2, REAL(ex_X), &N, REAL(x_jstar), &one, &betaOne, REAL(g_t), &one);
		
		
		flag_St =0;		  
		for( i=0; i<p; i++){
			//	REAL(absg_t)[i] = fabs(REAL(g_t)[i]);
			REAL(Lambda_t)[i] = REAL(g_t)[i] / REAL(P_t)[i];
			REAL(absLambda_t)[i] = fabs(REAL(Lambda_t)[i]);
			REAL(absg_t)[i] = fabs(REAL(g_t)[i]);
			if(REAL(Lambda_t)[i]*REAL(beta_vec)[i]< 0  && REAL(absg_t)[i] > 2.0 * REAL(ex_delta_t)[0] / (double) N ){
			    REAL(absLambda_t)[i] = REAL(absLambda_t)[i] + 1000000000000.0;
		        flag_St =1;		  
			}
			if(REAL(absg_t)[i] < 2.0 * REAL(ex_delta_t)[0] / (double) N ){
				REAL(absLambda_t)[i] = 0.0;
				flag_St =1;		  
			}
		}
		
		
		jstar = F77_CALL(idamax)(&p,REAL(absLambda_t),&one);
		jstar0 = F77_CALL(idamax)(&p,REAL(Lambda_t),&one);
		jstar = jstar-1;
		jstar0 = jstar0-1;
		
		
		
		flag_XTxjstar=0;
		if(INTEGER(selected_variable_index_vec)[jstar]==0){
			INTEGER(selected_variable_index_vec)[jstar] = jsort+1;
			INTEGER(selected_variable_index_sort_vec)[jsort] = jstar+1;
			jsort=jsort+1;
			flag_XTxjstar=1;
		}
		
		jstar_sort=INTEGER(selected_variable_index_vec)[jstar]-1;
		
		
		INTEGER(betahat_index_vec_adj)[step0] = INTEGER(selected_variable_index_vec)[jstar];
		
		//beta_t[jstar] =   beta_t[jstar] + delta_t * sign(g_t[jstar])
		if(REAL(Lambda_t)[jstar]>0){
			sign_gt_jstar=1.0;
			sign_gt_jstar_int = 1;
		}
		else {
			sign_gt_jstar=-1.0;
			sign_gt_jstar_int = -1;
		}			
		REAL(beta_vec)[jstar] = REAL(beta_vec)[jstar] + REAL(ex_delta_t)[0] * sign_gt_jstar;
		
        INTEGER(betahat_index_vec)[step0] = (jstar+1) * sign_gt_jstar_int;
		
		//cov_penaltyのための係数計算
		REAL(increment_covpenalty_vec)[step0] = 2.0 * REAL(ex_delta_t)[0] / ( (double) N * fabs(REAL(g_t)[jstar]) );
		
		//g_t0をy-Xbetaと定義
		//g_t0=g_t0 - delta_beta*x[jstar]*sign
		for( i=0; i<N; i++){
			REAL(x_jstar)[i] = REAL(ex_X)[i +jstar*N];
		}
		
		for( i=0; i<Ntest; i++){
			REAL(x_jstar_test)[i] = REAL(ex_X_test)[i +jstar*Ntest];
		}
		
		alpha_minusdelta_t = -1.0*REAL(ex_delta_t)[0]* sign_gt_jstar;
		F77_CALL(daxpy)(&N,&alpha_minusdelta_t,REAL(x_jstar),&one,REAL(g_t0), &one);
		F77_CALL(daxpy)(&Ntest,&alpha_minusdelta_t,REAL(x_jstar_test),&one,REAL(g_t0_test), &one);
		REAL(RSS_vec)[step0] = dnrm2_(&N,REAL(g_t0),&one);
		REAL(RSS_vec)[step0] = REAL(RSS_vec)[step0]*REAL(RSS_vec)[step0];		
		REAL(RSS_vec_test)[step0] = dnrm2_(&Ntest,REAL(g_t0_test),&one);
		REAL(RSS_vec_test)[step0] = REAL(RSS_vec_test)[step0]*REAL(RSS_vec_test)[step0];	
		
		
		
		//tuning parameterの計算
		if(penalty_type_int == 0){ //enet
			for( i=0; i<p; i++){
				temp = REAL(ex_para)[0] * REAL(beta_vec)[i]*REAL(beta_vec)[i] / 2.0 + (1.0 - REAL(ex_para)[0])*fabs(REAL(beta_vec)[i]); 
				REAL(tuning)[step0] = REAL(tuning)[step0] + temp; 
				REAL(tuning_stand)[step0] = REAL(tuning_stand)[step0] + REAL(ex_standardize_vec)[i] * temp; 
			}
		}
		
		if(penalty_type_int == 1){ //genet
			for( i=0; i<p; i++){
				temp = log( (1.0 - REAL(ex_para)[0])* fabs(REAL(beta_vec)[i]) +  REAL(ex_para)[0] ) - log(REAL(ex_para)[0]); 
				REAL(tuning)[step0] = REAL(tuning)[step0] + temp;
				REAL(tuning_stand)[step0] = REAL(tuning_stand)[step0] + REAL(ex_standardize_vec)[i] * temp; 
			}
		}
		
		if(penalty_type_int == 2){ //alasso
			for( i=0; i<p; i++){
				temp = REAL(ex_weight)[i]*fabs(REAL(beta_vec)[i]); 
				REAL(tuning)[step0] = REAL(tuning)[step0] + temp;
				REAL(tuning_stand)[step0] = REAL(tuning_stand)[step0] + REAL(ex_standardize_vec)[i] * temp; 
			}
		}		

			//tuning_sum
		for( i=0; i<p; i++){
			temp = fabs(REAL(beta_vec)[i]);
			REAL(tuning_sum)[step0] = REAL(tuning_sum)[step0] + temp; 
			REAL(tuning_sum_stand)[step0] = REAL(tuning_sum_stand)[step0] + REAL(ex_standardize_vec)[i] * temp;
		}
			


		//dfが妥当かcheck!
		beta_vecminus2_flag = 0;
		//beta_vecminus2_index = step0-1;
		//if(fabs(REAL(df_vec)[step0] - REAL(df_vec)[step0-1]) > REAL(ex_delta_t)[0] * 100.0) beta_vecminus2_flag=1;
		//収束してるのかcheck!
		check2_int = 0;
		for( i=0; i<p; i++){
			if(REAL(absg_t)[i] < 2.0 * REAL(ex_delta_t)[0] / (double) N) check2_int=check2_int+1;
		}
		if(check2_int==p) beta_vecminus2_flag=1;
		
		//収束してるのかcheck! ver2
		REAL(sum_lambda_t)[step0] = F77_CALL(dasum)(&p,REAL(Lambda_t),&one);
		//if(REAL(sum_lambda_t)[step0] < REAL(tol)[1]) beta_vecminus2_flag=1;
		//さらにちぇっっく
		if(step0>check_int){
			for( i=step0-check_int; i<step0; i++){
				if(fabs(REAL(sum_lambda_t)[step0] - REAL(sum_lambda_t)[i]) == 0.0 ){
					beta_vecminus2_flag=1;
				}
			}
		}
		//pmaxに達したのかcheck
		if(jsort==pmax+1) beta_vecminus2_flag=1;
		if(beta_vecminus2_flag==1) break;
		
		
	}	
	
	
	
	
	INTEGER(STEP_adj)[0] = step0;  
	
	
	
	SET_VECTOR_ELT(ans, 0, betahat_index_vec);
	SET_VECTOR_ELT(ans, 1, STEP_adj);	
	SET_VECTOR_ELT(ans, 2, RSS_vec_test);	
	SET_VECTOR_ELT(ans, 3, tuning);
	SET_VECTOR_ELT(ans, 4, tuning_stand);
	SET_VECTOR_ELT(ans, 5, tuning_sum);
	SET_VECTOR_ELT(ans, 6, tuning_sum_stand);
	UNPROTECT(27);
	
	return(ans);
	
	
}


SEXP findtuning(SEXP tuning_kouho, SEXP tuning_eval, SEXP RSSvec)
{
	int STEP = LENGTH(tuning_kouho);
	long int i;
	long int j;
	long int j0;
	int one=1;
	SEXP ans;
	SEXP comparison;
	PROTECT(comparison = allocVector(REALSXP,STEP));
	PROTECT(ans = allocVector(REALSXP,STEP));
	
	for( i=0; i<STEP; i++){
		REAL(ans)[i] = 0.0;
		REAL(comparison)[i] = 0.0;
	}
	
	for( i=0; i<STEP; i++){
		for( j=0; j<STEP; j++){
			REAL(comparison)[j] = REAL(tuning_kouho)[j] - REAL(tuning_eval)[i] ;
			REAL(comparison)[j] = 1.0 / (1.0 + fabs(REAL(comparison)[j]));
		}
		j0 = F77_CALL(idamax)(&STEP,REAL(comparison),&one);
		REAL(ans)[i] = REAL(RSSvec)[j0-1];
	}
			
			
	UNPROTECT(2);
	return(ans);
}


SEXP findtuning2(SEXP tuning_candidate, SEXP tuning)
{
	long int STEP_candidate = LENGTH(tuning_candidate);
	long int STEP = LENGTH(tuning);
	long int i;
	long int j;
	long int j0;
	int one=1;
	SEXP ans;
	SEXP comparison;
	PROTECT(comparison = allocVector(REALSXP,STEP_candidate));
	PROTECT(ans = allocVector(INTSXP,STEP));
	
	for( i=0; i<STEP_candidate; i++){
		REAL(comparison)[i] = 0.0;
	}
	
	for( i=0; i<STEP; i++){
		INTEGER(ans)[i] = 1;
	}

	
	
	for( i=0; i<STEP; i++){
		for( j=0; j<STEP_candidate; j++){
			REAL(comparison)[j] = REAL(tuning)[i] - REAL(tuning_candidate)[j] ;
			REAL(comparison)[j] = 1.0 / (1.0 + fabs(REAL(comparison)[j]));
		}
		j0 = F77_CALL(idamax)(&STEP_candidate,REAL(comparison),&one);
		INTEGER(ans)[i] = j0;
			}
			
			
		UNPROTECT(2);
	return(ans);
}



SEXP DFNAIVE(SEXP ex_X, SEXP ex_y, SEXP ex_betahat_index_vec, SEXP ex_STEP_adj, SEXP ex_increment_vec){
	
	
	long int step0;
	long int jstar;
	long int i;
	long int N=length(ex_y);
	long int NN=N*N;
	long int STEP_adj = INTEGER(ex_STEP_adj)[0];
	int one=1;
	
	double increment_covpenalty;
	double increment_covpenalty_minus;
	
	SEXP cov_penalty_matrix;
	SEXP xx_matrix;
    SEXP xx_cov_matrix;
    SEXP x_cov_vec;
    SEXP df_vec;
    SEXP x_jstar;
	
	PROTECT(cov_penalty_matrix = allocMatrix(REALSXP, N,N)); 
	PROTECT(xx_matrix = allocMatrix(REALSXP, N,N)); 
	PROTECT(xx_cov_matrix = allocMatrix(REALSXP, N,N)); 
	PROTECT(x_cov_vec = allocVector(REALSXP, N)); 
	PROTECT(df_vec = allocVector(REALSXP, STEP_adj+1)); 
	PROTECT(x_jstar = allocVector(REALSXP,N));
	
	for( i=0; i<STEP_adj+1; i++){
		REAL(df_vec)[i] = 0.0;
	}
	
	for( i=0; i<NN; i++){
		REAL(cov_penalty_matrix)[i] = 0.0;
		REAL(xx_matrix)[i] = 0.0;
		REAL(xx_cov_matrix)[i] = 0.0;
	}
	
	
	for( i=0; i<N; i++){
		REAL(x_cov_vec)[i] = 0.0;
		REAL(x_jstar)[i] = 0.0;
	}	    	
	
	
	//covariance penalty の計算
	//\frac{\mathrm{cov}(\hat{\bm{y}}(t+\Delta t),\bm{y})}{\tau^2} =\frac{\mathrm{cov}(\hat{\bm{y}}(t) ,\bm{y})}{\tau^2}  + \frac{2\Delta t}{N|g_{k}(t) |}  \ \bm{x}_{k} \bm{x}_{k}^T  \left\{ I- \frac{\mathrm{cov}(\hat{\bm{y}}(t) ,\bm{y})}{\tau^2} \right\}
	//初期化
	for( step0=1; step0<STEP_adj+1; step0++){
		F77_CALL(dscal)(&N,&Zero,REAL(x_cov_vec),&one);
		F77_CALL(dscal)(&NN,&Zero,REAL(xx_cov_matrix),&one);
		F77_CALL(dscal)(&NN,&Zero,REAL(xx_matrix),&one);   
		
		
		jstar = abs(INTEGER(ex_betahat_index_vec)[step0])-1;
		for( i=0; i<N; i++){
			REAL(x_jstar)[i] = REAL(ex_X)[i +jstar*N];
		}   		
		
		increment_covpenalty = REAL(ex_increment_vec)[step0];
		increment_covpenalty_minus = -increment_covpenalty;
		//	  for( i=0; i<N; i++){
		//		REAL(x_jstar)[i] = REAL(ex_X)[i +jstar*N];
		//	}
		// cov * x
		F77_CALL(dsymv)( &UPPER, &N, &alphaOne, REAL(cov_penalty_matrix),  &N, REAL(x_jstar), &one, &betaZero, REAL(x_cov_vec) ,&one);
		// xx* cov
		F77_CALL(dger)(&N, &N,  &increment_covpenalty_minus, REAL(x_jstar), &one,  REAL(x_cov_vec), &one, REAL(xx_cov_matrix), &N);
		// xx* (I-cov)
		F77_CALL(dger)(&N, &N,  &increment_covpenalty, REAL(x_jstar), &one,  REAL(x_jstar), &one, REAL(xx_cov_matrix), &N);
		
		// covの計算
		F77_CALL(daxpy)(&NN, &alphaOne, REAL(xx_cov_matrix),&one,REAL(cov_penalty_matrix),&one);
		
		
		//dfの計算
		//REAL(df_vec)[step0] = REAL(df_vec)[step0-1];
		for( i=0; i<N; i++){
			REAL(df_vec)[step0] = REAL(df_vec)[step0] + REAL(cov_penalty_matrix)[i*N + i];
		}
		
	}
	
	
	
	UNPROTECT(6);
	
	return(df_vec);
	
}




SEXP DFNAIVE2(SEXP ex_X_selected, SEXP ex_y, SEXP ex_betahat_index_vec_adj, SEXP ex_STEP_adj, SEXP ex_increment_vec){
	
	
	long int step0;
	long int jstar;
	long int i;
	long int N=length(ex_y);
	long int NN=N*N;
	long int STEP_adj = INTEGER(ex_STEP_adj)[0];
	int one=1;
	
	double increment_covpenalty;
	double increment_covpenalty_minus;
	
	SEXP Bt;
	SEXP xx_matrix;
    SEXP xx_cov_matrix;
    SEXP x_cov_vec;
    SEXP df_vec;
    SEXP x_jstar;
	
	PROTECT(Bt = allocMatrix(REALSXP, N,N)); 
	PROTECT(xx_matrix = allocMatrix(REALSXP, N,N)); 
	PROTECT(xx_cov_matrix = allocMatrix(REALSXP, N,N)); 
	PROTECT(x_cov_vec = allocVector(REALSXP, N)); 
	PROTECT(df_vec = allocVector(REALSXP, STEP_adj+1)); 
	PROTECT(x_jstar = allocVector(REALSXP,N));
	
	for( i=0; i<STEP_adj+1; i++){
		REAL(df_vec)[i] = 0.0;
	}
	
	for( i=0; i<NN; i++){
		REAL(Bt)[i] = 0.0;
		REAL(xx_matrix)[i] = 0.0;
		REAL(xx_cov_matrix)[i] = 0.0;
	}
	
	
	for( i=0; i<N; i++){
		REAL(x_cov_vec)[i] = 0.0;
		REAL(x_jstar)[i] = 0.0;
	}	
	
	for( i=0; i<N; i++){
		REAL(Bt)[i + i*N] = 1.0;
	}
	
	
	//covariance penalty の計算
	//\frac{\mathrm{cov}(\hat{\bm{y}}(t+\Delta t),\bm{y})}{\tau^2} =\frac{\mathrm{cov}(\hat{\bm{y}}(t) ,\bm{y})}{\tau^2}  + \frac{2\Delta t}{N|g_{k}(t) |}  \ \bm{x}_{k} \bm{x}_{k}^T  \left\{ I- \frac{\mathrm{cov}(\hat{\bm{y}}(t) ,\bm{y})}{\tau^2} \right\}
	//初期化
	for( step0=1; step0<STEP_adj+1; step0++){
		F77_CALL(dscal)(&N,&Zero,REAL(x_cov_vec),&one);
		F77_CALL(dscal)(&NN,&Zero,REAL(xx_cov_matrix),&one);
		F77_CALL(dscal)(&NN,&Zero,REAL(xx_matrix),&one);   
		
		
		jstar = abs(INTEGER(ex_betahat_index_vec_adj)[step0])-1;
		for( i=0; i<N; i++){
			REAL(x_jstar)[i] = REAL(ex_X_selected)[i +jstar*N];
		}   		
		
		increment_covpenalty = REAL(ex_increment_vec)[step0];
		increment_covpenalty_minus = -increment_covpenalty;
		//	  for( i=0; i<N; i++){
		//		REAL(x_jstar)[i] = REAL(ex_X)[i +jstar*N];
		//	}
		// cov * x
		//			  F77_CALL(dsymv)( &UPPER, &N, &alphaOne, REAL(Bt),  &N, REAL(x_jstar), &one, &betaZero, REAL(x_cov_vec) ,&one);
		F77_CALL(dgemv)( &TRANST, &N,&N, &alphaOne, REAL(Bt),  &N, REAL(x_jstar), &one, &betaZero, REAL(x_cov_vec) ,&one);
		// cov - xx* cov
		F77_CALL(dger)(&N, &N,  &increment_covpenalty_minus, REAL(x_jstar), &one,  REAL(x_cov_vec), &one, REAL(Bt), &N);
		
		
		//dfの計算
		//REAL(df_vec)[step0] = REAL(df_vec)[step0-1];
		for( i=0; i<N; i++){
			REAL(df_vec)[step0] = REAL(df_vec)[step0] - REAL(Bt)[i*N + i];
		}
		REAL(df_vec)[step0] = REAL(df_vec)[step0] + (double) N;
		
	}
	
	
	
	UNPROTECT(6);
	
	return(df_vec);
	
}



SEXP DFMODIFIED(SEXP ex_qr_X, SEXP ex_y, SEXP ex_betahat_index_vec_adj, SEXP ex_STEP_adj, SEXP ex_increment_vec, SEXP ex_selected_variable_index){
	
	
	long int q_N = INTEGER(GET_DIM(ex_qr_X))[0];
	long int q_Nq_N = q_N*q_N;
	long int step0;
	long int jstar;
	long int i;
	//	long int N=length(ex_y);
	//	long int NN=N*N;
	long int STEP_adj = INTEGER(ex_STEP_adj)[0];
	int one=1;
	
	double increment_covpenalty;
	double increment_covpenalty_minus;
	
	SEXP Bt;
	SEXP rr_matrix;
    SEXP rr_cov_matrix;
    SEXP r_cov_vec;
    SEXP df_vec;
    SEXP r_jstar;
	
	PROTECT(Bt = allocMatrix(REALSXP, q_N,q_N)); 
	PROTECT(rr_matrix = allocMatrix(REALSXP, q_N,q_N)); 
	PROTECT(rr_cov_matrix = allocMatrix(REALSXP, q_N,q_N)); 
	PROTECT(r_cov_vec = allocVector(REALSXP, q_N)); 
	PROTECT(df_vec = allocVector(REALSXP, STEP_adj+1)); 
	PROTECT(r_jstar = allocVector(REALSXP,q_N));
	
	
	
	//covariance penalty の計算
	//\frac{\mathrm{cov}(\hat{\bm{y}}(t+\Delta t),\bm{y})}{\tau^2} =\frac{\mathrm{cov}(\hat{\bm{y}}(t) ,\bm{y})}{\tau^2}  + \frac{2\Delta t}{N|g_{k}(t) |}  \ \bm{x}_{k} \bm{x}_{k}^T  \left\{ I- \frac{\mathrm{cov}(\hat{\bm{y}}(t) ,\bm{y})}{\tau^2} \right\}
	//初期化
	
	
	
	for( i=0; i<STEP_adj+1; i++){
		REAL(df_vec)[i] = 0.0;
	}
	
	for( i=0; i<q_Nq_N; i++){
		REAL(Bt)[i] = 0.0;
		REAL(rr_matrix)[i] = 0.0;
		REAL(rr_cov_matrix)[i] = 0.0;
	}
	
	
	for( i=0; i<q_N; i++){
		REAL(Bt)[i + i*q_N] = 1.0;
	}	    	
	
	for( i=0; i<q_N; i++){
		REAL(r_cov_vec)[i] = 0.0;
		REAL(r_jstar)[i] = 0.0;
	}
	
	
	//START!!
	
	for( step0=1; step0<STEP_adj+1; step0++){
		F77_CALL(dscal)(&q_N,&Zero,REAL(r_cov_vec),&one);
		//			dscal_(&q_Nq_N,&Zero,REAL(rr_cov_matrix),&one);
		//			dscal_(&q_Nq_N,&Zero,REAL(rr_matrix),&one);   
		
		
		
		jstar = INTEGER(ex_betahat_index_vec_adj)[step0]-1;
		for( i=0; i<q_N; i++){
			REAL(r_jstar)[i] = REAL(ex_qr_X)[i +jstar*q_N];
		}   		
		
		increment_covpenalty = REAL(ex_increment_vec)[step0];
		increment_covpenalty_minus = -increment_covpenalty;
		//	  for( i=0; i<N; i++){
		//		REAL(x_jstar)[i] = REAL(ex_X)[i +jstar*N];
		//	}
		// cov * r
		
		//F77_CALL(dsymv)( &UPPER, &q_N, &alphaOne, REAL(Bt),  &q_N, REAL(r_jstar), &one, &betaZero, REAL(r_cov_vec) ,&one);
		F77_CALL(dgemv)( &TRANST, &q_N, &q_N, &alphaOne, REAL(Bt),  &q_N, REAL(r_jstar), &one, &betaZero, REAL(r_cov_vec) ,&one);
		// Bt - r(cov*r)^T
		F77_CALL(dger)(&q_N, &q_N,  &increment_covpenalty_minus, REAL(r_jstar), &one,  REAL(r_cov_vec), &one, REAL(Bt), &q_N);
		//	          // rr* (I-cov)
		//	          dger_(&q_N, &q_N,  &increment_covpenalty, REAL(r_jstar), &one,  REAL(r_jstar), &one, REAL(rr_cov_matrix), &q_N);
		//	          
		//			  // covの計算
		//			  daxpy_(&q_Nq_N, &alphaOne, REAL(rr_cov_matrix),&one,REAL(Bt),&one);
		
		
		//dfの計算
		//REAL(df_vec)[step0] = REAL(df_vec)[step0-1];
		for( i=0; i<q_N; i++){
			REAL(df_vec)[step0] = REAL(df_vec)[step0] - REAL(Bt)[i*q_N + i];
		}
		
		REAL(df_vec)[step0] = REAL(df_vec)[step0] + (double) q_N;
		
		
	}
	
	
	
	UNPROTECT(6);
	
	return(df_vec);
	
}







SEXP DFMODIFIED2(SEXP ex_qr_X, SEXP ex_y, SEXP ex_betahat_index_vec_adj, SEXP ex_STEP_adj, SEXP ex_increment_vec, SEXP ex_selected_variable_index){
	
	
	long int q_N = INTEGER(GET_DIM(ex_qr_X))[0];
	long int q_Nq_N = q_N*q_N;
	long int step0;
	long int jstar;
	long int i,j;
	//	long int N=length(ex_y);
	//	long int NN=N*N;
	long int STEP_adj = INTEGER(ex_STEP_adj)[0];
	int one=1;
	long int q_N_adj = 0;
	int flag_BtBt0;
	
	double increment_covpenalty;
	double increment_covpenalty_minus;
	
	SEXP Bt;
	SEXP Bt0; //けちな計算を行うための行列
	SEXP rr_matrix;
    SEXP rr_cov_matrix;
    SEXP r_cov_vec;
    SEXP df_vec;
    SEXP r_jstar;
	
	PROTECT(Bt = allocMatrix(REALSXP, q_N,q_N)); 
	PROTECT(Bt0 = allocMatrix(REALSXP, q_N, q_N)); 
	PROTECT(rr_matrix = allocMatrix(REALSXP, q_N,q_N)); 
	PROTECT(rr_cov_matrix = allocMatrix(REALSXP, q_N,q_N)); 
	PROTECT(r_cov_vec = allocVector(REALSXP, q_N)); 
	PROTECT(df_vec = allocVector(REALSXP, STEP_adj+1)); 
	PROTECT(r_jstar = allocVector(REALSXP,q_N));
	
	
	
	//covariance penalty の計算
	//\frac{\mathrm{cov}(\hat{\bm{y}}(t+\Delta t),\bm{y})}{\tau^2} =\frac{\mathrm{cov}(\hat{\bm{y}}(t) ,\bm{y})}{\tau^2}  + \frac{2\Delta t}{N|g_{k}(t) |}  \ \bm{x}_{k} \bm{x}_{k}^T  \left\{ I- \frac{\mathrm{cov}(\hat{\bm{y}}(t) ,\bm{y})}{\tau^2} \right\}
	//初期化
	
	
	
	for( i=0; i<STEP_adj+1; i++){
		REAL(df_vec)[i] = 0.0;
	}
	
	for( i=0; i<q_Nq_N; i++){
		REAL(Bt)[i] = 0.0;
		REAL(Bt0)[i] = 0.0;
		REAL(rr_matrix)[i] = 0.0;
		REAL(rr_cov_matrix)[i] = 0.0;
	}
	
	
	for( i=0; i<q_N; i++){
		REAL(Bt)[i + i*q_N] = 1.0;
	}	    	
	
	for( i=0; i<q_N; i++){
		REAL(r_cov_vec)[i] = 0.0;
		REAL(r_jstar)[i] = 0.0;
	}
	
	
	//START!!
	
	for( step0=1; step0<STEP_adj+1; step0++){
		F77_CALL(dscal)(&q_N,&Zero,REAL(r_cov_vec),&one);
		
		//			dscal_(&q_Nq_N,&Zero,REAL(rr_cov_matrix),&one);
		//			dscal_(&q_Nq_N,&Zero,REAL(rr_matrix),&one);   
		
		
		flag_BtBt0 = 0;
		jstar = INTEGER(ex_betahat_index_vec_adj)[step0]-1;
		if(q_N_adj  == 	jstar && q_N_adj < q_N){
			q_N_adj = q_N_adj+1;
			flag_BtBt0=1;
			if(q_N_adj==1){
				REAL(Bt0)[0]=1.0;
			}else {
				for( i=0; i<(q_N_adj-1); i++){
					for( j=0; j<(q_N_adj-1); j++){
						REAL(Bt)[j+i*q_N] = REAL(Bt0)[j + i*(q_N_adj-1)];
					}
				}
				//dscal_(&q_Nq_N,&Zero,REAL(Bt0),&one);   
				for( i=0; i<q_N_adj; i++){
					for( j=0; j<q_N_adj; j++){
						REAL(Bt0)[j+i*q_N_adj] = REAL(Bt)[j+i*q_N];
					}
				}
			}
		 	
		}
		
		
		for( i=0; i<q_N; i++){
			REAL(r_jstar)[i] = REAL(ex_qr_X)[i +jstar*q_N];
		}   		
		
		increment_covpenalty = REAL(ex_increment_vec)[step0];
		increment_covpenalty_minus = -increment_covpenalty;
		//	  for( i=0; i<N; i++){
		//		REAL(x_jstar)[i] = REAL(ex_X)[i +jstar*N];
		//	}
		
		
		//つじつまがあうよう行列作成
		
		// cov * r
		F77_CALL(dgemv)( &TRANST, &q_N_adj, &q_N_adj, &alphaOne, REAL(Bt0),  &q_N_adj, REAL(r_jstar), &one, &betaZero, REAL(r_cov_vec) ,&one);
		// Bt - r(cov*r)^T
		F77_CALL(dger)(&q_N_adj, &q_N_adj,  &increment_covpenalty_minus, REAL(r_jstar), &one,  REAL(r_cov_vec), &one, REAL(Bt0), &q_N_adj);
		
		
		
		
		//dfの計算
		//REAL(df_vec)[step0] = REAL(df_vec)[step0-1];
		for( i=0; i<q_N_adj; i++){
			REAL(df_vec)[step0] = REAL(df_vec)[step0] - REAL(Bt0)[i*q_N_adj + i];
		}
		
		REAL(df_vec)[step0] = REAL(df_vec)[step0] + (double) q_N_adj;
		
		
	}
	
	
	
	UNPROTECT(7);
	
	return(df_vec);
	
}


SEXP DFMODIFIED3(SEXP ex_betahat_index_vec_adj, SEXP ex_increment_vec, SEXP ex_p_adj, SEXP ex_STEP_adj){
	
	
	long int p_adj = INTEGER(ex_p_adj)[0];
	long int STEP_adj = INTEGER(ex_STEP_adj)[0];
	
	long int i;
	long int j;
	
    SEXP df0;
    SEXP df_vec;
	
	PROTECT(df0 = allocVector(REALSXP, p_adj)); 
	PROTECT(df_vec = allocVector(REALSXP, STEP_adj)); 

	for( i=0; i<STEP_adj; i++){
		REAL(df_vec)[i] = 0.0;
	}
	
	for( i=0; i<p_adj; i++){
		REAL(df0)[i] = 1.0;
	}
	
	for(i=1; i<STEP_adj; i++){
		REAL(df0)[INTEGER(ex_betahat_index_vec_adj)[i]-1] = REAL(df0)[INTEGER(ex_betahat_index_vec_adj)[i]-1] * REAL(ex_increment_vec)[i];
		for(j=0; j<p_adj; j++){
			REAL(df_vec)[i] = REAL(df_vec)[i] + REAL(df0)[j];
		}
		REAL(df_vec)[i] = (double) p_adj - REAL(df_vec)[i];
	}

	
	UNPROTECT(2);
	
	return(df_vec);
	
}



SEXP betaOUT(SEXP ex_betahat_index_vec, SEXP ex_step, SEXP ex_p, SEXP ex_delta_t, SEXP ex_step_adj, SEXP ex_stepmax){
	long int p = INTEGER(ex_p)[0];
	long int n_tuning = length(ex_step);
	long int stepmax = INTEGER(ex_stepmax)[0];
	long int i;
	long int j;
	long int j1;
	long int flag_is=0;
	long int i0;
	int sign_int;
	SEXP beta;
	SEXP betaMATRIX;
	double sign;
	PROTECT(beta = allocVector(REALSXP,p));
	PROTECT(betaMATRIX = allocMatrix(REALSXP,p,n_tuning));
	//PROTECT(sign = allocVector(REALSXP,1));
	
	for( i=0; i<p; i++){
		REAL(beta)[i] = 0.0;
	}	
	
	
	for( i=0; i<p*n_tuning; i++){
		REAL(betaMATRIX)[i] = 0.0;
	}	
	
	
	for( i=1; i<stepmax+1; i++){
		i0 = INTEGER(ex_betahat_index_vec)[i];
		if(i0>0){
			sign = 1.0;
			sign_int = 1;
		}
		else {
		    sign = -1.0;
			sign_int = -1;
		}
		i0=i0*sign_int;
		REAL(beta)[i0-1] = REAL(beta)[i0-1] + sign * REAL(ex_delta_t)[0];
		if(INTEGER(ex_step_adj)[i-1] > 0){
			for( j1=0; j1<INTEGER(ex_step_adj)[i-1]; j1++){
				for( j=0; j<p; j++){
					REAL(betaMATRIX)[flag_is*p + j] = REAL(beta)[j];
				}
				flag_is=flag_is + 1;
			}
		}
	}
    
    UNPROTECT(2);
    
    return(betaMATRIX);
	
	
	
}


SEXP betaOUT_MATRIX(SEXP ex_betahat_index_vec, SEXP ex_p, SEXP ex_delta_t){
	long int STEP=length(ex_betahat_index_vec);
	long int p = INTEGER(ex_p)[0];
	long int i;
	long int j;
	long int i0;
	int sign_int;
	SEXP beta;
	SEXP betaMATRIX;
	double sign;
	PROTECT(beta = allocVector(REALSXP,p));
	PROTECT(betaMATRIX = allocMatrix(REALSXP,p,STEP));
	//PROTECT(sign = allocVector(REALSXP,1));
	
	for( i=0; i<p; i++){
		REAL(beta)[i] = 0.0;
	}	
	
	for( i=0; i<p*STEP; i++){
		REAL(betaMATRIX)[i] = 0.0;
	}
	
	for( i=1; i<STEP; i++){
		i0 = INTEGER(ex_betahat_index_vec)[i];
		if(i0>0){
			sign = 1.0;
			sign_int = 1;
		}
		else {
		    sign = -1.0;
			sign_int = -1;
		}
		i0=i0*sign_int;
		REAL(beta)[i0-1] = REAL(beta)[i0-1] + sign * REAL(ex_delta_t)[0];
		for( j=0; j<p; j++){
			REAL(betaMATRIX)[i*p + j] = REAL(beta)[j];
		}
	}
    
    UNPROTECT(2);
    
    return(betaMATRIX);
	
	
	
}



SEXP tOUT(SEXP ex_betahat_index_vec, SEXP ex_p, SEXP ex_delta_t){
	long int STEP=length(ex_betahat_index_vec);
	long int p = INTEGER(ex_p)[0];
	long int i;
	long int i0;
	int sign_int;
	SEXP beta;
	SEXP tout;
	double sign;
	PROTECT(beta = allocVector(REALSXP,p));
	PROTECT(tout = allocVector(REALSXP,STEP));
	//PROTECT(sign = allocVector(REALSXP,1));
	
	for( i=0; i<p; i++){
		REAL(beta)[i] = 0.0;
	}	
	for( i=0; i<STEP; i++){
		REAL(tout)[i] = 0.0;
	}	
	
	
	
	
	for( i=1; i<STEP; i++){
		i0 = INTEGER(ex_betahat_index_vec)[i];
		if(i0>0){
			sign = 1.0;
			sign_int = 1;
		}
		else {
		    sign = -1.0;
			sign_int = -1;
		}
		i0=i0*sign_int;
		REAL(beta)[i0-1] = REAL(beta)[i0-1] + sign * REAL(ex_delta_t)[0];
		if(REAL(beta)[i0-1] * sign < 0.0){
			REAL(tout)[i] = REAL(tout)[i-1] - REAL(ex_delta_t)[0];
		}else{
			REAL(tout)[i] = REAL(tout)[i-1] +  REAL(ex_delta_t)[0];
		}
		
	}
    
    UNPROTECT(2);
    
    return(tout);
}
