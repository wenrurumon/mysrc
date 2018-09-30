#include <RcppArmadillo.h>
#include <unistd.h>
#include <math.h>


inline double covsel_obj(const arma::mat & s, const arma::mat & x, const arma::mat & z, const double & lambda){
  return( trace(s*x) -log(det(x)) + lambda*accu(abs(z)) );
}


double norm_fro(const arma::mat & x){
  return( sqrt(trace( (x.t() * x) )) );
}


// [[Rcpp::export]]
arma::mat shrinkage(const arma::mat & a, double kappa) {
  /*
  Rcpp::Rcout << "kappa " << kappa << std::endl
              << "a.max " << a.max() << std::endl
              << "a.min " << a.min() << std::endl;
  */
  arma::mat s1 = clamp(  a-kappa, 0,  std::max( a.max()-kappa, kappa));
  arma::mat s2 = clamp( -a-kappa, 0,  std::max(-a.min()-kappa, kappa));
  s1 = s1 - s2;
  return( s1 );
}

// [[Rcpp::export]]
arma::mat covsel(const arma::mat & d, const int max_iter=1000, double alpha=1.2,
           double lambda=0.0001, double rho=0.5){
  double ABSTOL = 1e-5;
  double RELTOL = 1e-3;
  arma::mat s = cov(d);
  unsigned int n = s.n_rows;
  arma::mat x = arma::zeros(n,n);
  arma::mat z = arma::zeros(n,n);
  arma::mat u = arma::zeros(n,n);

  double  r_norm = 0, s_norm = 0,  eps_pri = 0, eps_dual = 0;
  //double objval;
  for( int i = 0 ; i < max_iter && r_norm <= eps_pri && s_norm <= eps_dual; i++){
    arma::vec es;
    arma::mat Q;
    eig_sym(es,Q,rho*(z-u) - s);
    arma::vec xi = (es +sqrt(pow(es,2)+4*rho))/(2*rho);
    x = Q * diagmat(xi) * Q.t();

    arma::mat zold = z;
    arma::mat x_hat = alpha*x + (1-alpha)*zold;
    z = shrinkage( x_hat +u, lambda/rho);

    u = u + (x_hat - z);

    //objval =  covsel_obj(s,x,z,lambda);
    r_norm = norm_fro( x - z);
    s_norm = -rho*norm_fro(z - zold);

    eps_pri = n*ABSTOL + RELTOL*std::max(norm_fro(x), norm_fro(z) );
    eps_dual = n*ABSTOL + rho*RELTOL*norm_fro(u);
    /*
    Rcpp::Rcout << "i " << i  << std::endl;
    Rcpp::Rcout << "r_norm " << r_norm <<std::endl;
    Rcpp::Rcout << "s_norm " << s_norm <<std::endl;
    Rcpp::Rcout << "eps_pri " << eps_pri <<std::endl;
    Rcpp::Rcout << "eps_dual " << eps_dual <<std::endl;
    Rcpp::Rcout << "objval " << objval <<std::endl;
    */
  }
  return x;

}

void factor(const arma::mat & A, arma::mat & L, const double & rho){
  unsigned m = A.n_rows;
  unsigned n = A.n_cols;

  if ( m >=n ){
    chol(L, A.t()*A+rho*arma::eye(n,n), "lower");
  }else{
    chol(L, arma::eye(m,m)+1/rho*(A*A.t()),"lower");
  }
}

double norm_1(const arma::mat & x){
  return( arma::as_scalar( max(sum(abs(x) ) ) ) );
}

double norm(const arma::mat & x){
  return( arma::as_scalar( max(svd(x) ) ) );
}

double objective_lasso(const arma::mat & A, const arma::vec & b, const arma::vec & x,
                       const arma::vec & z, double lambda=0.0001){
  return( 0.5*sum(pow(A*x - b, 2))+ lambda*norm_1(z) );
}


// [[Rcpp::export]]
arma::vec lasso(const arma::mat & A, const arma::vec & b, const int max_iter=1000, double alpha=1.2,
          double lambda=0.1, double rho=1){
  double ABSTOL = 1e-4;
  double RELTOL = 1e-2;
  unsigned int m = A.n_rows;
  unsigned int n = A.n_cols;

  arma::mat At = A.t();
  arma::mat Atb = At*b;
  arma::mat L;

  factor(A,L,rho);
  arma::mat U = L.t();

  arma::vec x = arma::zeros(n);
  arma::vec z = arma::zeros(n);
  arma::vec u = arma::zeros(n);

  arma::vec q;
  double  r_norm = 0, s_norm = 0, eps_pri = 0, eps_dual = 0;
  //double objval;
  double sqrt_n = sqrt(n);
  double rho_2 = rho*rho;
  for( int i = 0 ; i < max_iter && r_norm <= eps_pri && s_norm <= eps_dual; i++){
    q = Atb + rho *(z-u);

    if( m >= n ){
      solve(x,U,solve(L,q));
    }else{
      x = q/rho-(At*(solve(U,solve(L, A*q))))/rho_2;
    }

    arma::vec zold = z;
    arma::vec x_hat = alpha*x + (1-alpha)*zold;
    z = shrinkage(x_hat+u,lambda/rho);

    u = u+(x_hat -z);
    //objval = objective_lasso(A,b,x,z,lambda);

    r_norm = norm( x - z);
    s_norm = -rho*norm(z - zold);

    eps_pri = sqrt_n*ABSTOL + RELTOL*std::max(norm(x), norm(-z) );
    eps_dual = sqrt_n*ABSTOL + rho*RELTOL*norm(u);
    /*
     Rcpp::Rcout << "i " << i  << std::endl;
     Rcpp::Rcout << "r_norm " << r_norm <<std::endl;
     Rcpp::Rcout << "s_norm " << s_norm <<std::endl;
     Rcpp::Rcout << "eps_pri " << eps_pri <<std::endl;
     Rcpp::Rcout << "eps_dual " << eps_dual <<std::endl;
     Rcpp::Rcout << "objval " << objval <<std::endl;
     */
  }
  return z;

}




// [[Rcpp::export]]
arma::mat ca_mean(const arma::mat & x){
  arma::mat m = sum(x)/x.n_rows;
  return(m);
}


// [[Rcpp::export]]
arma::mat ca_var(const arma::mat & x){
  return(var(x));
}



arma::mat center(const arma::mat & x){
  return( x - repmat(sum(x)/x.n_rows,x.n_rows,1) );
}

void center_replace( arma::mat & x){
  x = x - repmat(sum(x)/x.n_rows,x.n_rows,1);
}

void scale_replace( arma::mat & x){
  x = x / repmat(sqrt(var(x)),x.n_rows,1);
}

arma::mat scale( arma::mat & x){
  return( x / repmat(sqrt(var(x)),x.n_rows,1) );
}

// [[Rcpp::export]]
arma::mat center_scale(const arma::mat & x){
  arma::mat m = sum(x)/x.n_rows;
  return( (x -repmat(m,x.n_rows,1)) / repmat(sqrt(var(x)),x.n_rows,1) ); ////TODO: add check if any var = 0
}

void center_scale_replace(arma::mat & x){
  arma::mat m = sum(x)/x.n_rows;
  x = (x -repmat(m,x.n_rows,1)) / repmat(sqrt(var(x)),x.n_rows,1) ; ////TODO: add check if any var = 0
}



// [[Rcpp::export]]
arma::mat initial_sem_nofix( const arma::mat & data, const arma::uvec & variable,
                double rho = 1.0, double lambda = 0.1, double prec=1e-5, int max_iter = 1000){
  unsigned int n = data.n_rows;
  unsigned int d = data.n_cols - 1;

  int num_scr;
  if( d>n){
    if( n <=3 ){
      num_scr = n;
    }else{
      num_scr = std::ceil( n/log(n) );
    }
  }else{
    if(d<=4){
      num_scr = d;
    }else{
      num_scr = std::ceil( sqrt(d) );
    }
  }

  arma::mat td = center_scale(data) ;
  arma::mat coef_mat( td.n_cols, td.n_cols );
  coef_mat.fill(0.0);


  for( unsigned int index = 0; index < variable.n_rows; index++){
    unsigned int idx = variable(index) - 1;
    arma::mat W = td;
    W.shed_col( idx  );
    arma::mat Y = td.col(idx );
    arma::vec xy0 = W.t() * Y;
    arma::uvec order = sort_index(xy0,"descend");
    arma::uvec idx_scr = order.subvec(0,num_scr-1);
    W = W.cols(idx_scr);
    idx_scr.elem( find( idx_scr>=idx )  ) += 1;

    arma::vec beta = lasso(W,Y);
    //Rcpp::Rcout << "beta " << beta << std::endl;
    for(unsigned int i = 0; i < idx_scr.n_rows; i++ ){
      coef_mat(idx_scr(i),idx) = beta(i);
    }
    //Rcpp::Rcout << "idx_scr " << idx_scr << std::endl;
  }

  return( coef_mat.cols(variable-1).t() );
}


// [[Rcpp::export]]
arma::mat initial_sem(const arma::mat & data, const arma::uvec &  fixed, const arma::uvec & variable,
                double rho = 1.0, double lambda = 0.1, double prec=1e-5, int max_iter = 1000){
  arma::mat x = data.cols( fixed - 1);
  arma::mat xt = x.t();
  arma::mat y = data.cols(variable - 1);

  arma::mat y_hat = x* covsel(xt*x)*xt*y;

  return( initial_sem_nofix( join_rows(y_hat, x),  variable) );
}

arma::mat initial_sem_stable(const arma::uvec &  phenos, const arma::uvec &  fixed, const arma::mat & data,
                double rho = 1, int max_iter = 100, double lambda=0.1){
  unsigned int n = data.n_rows;
  unsigned int d = data.n_cols - 1;
  //unsigned int maxdf = (n<d)?d:n;
  arma::mat data_var = ca_var(data);
  arma::mat csx = center_scale(data);

  int num_scr1, num_scr2;
  if( d>n){
    if( n <=3 ){
      num_scr1 = n;
      num_scr2 = n;
    }else{
      num_scr1 = std::ceil( n/log(n) );
      num_scr2 = n-1;
    }
  }else{
    if(d<=4){
      num_scr1 = d;
      num_scr2 = d;
    }else{
      num_scr1 = std::ceil( sqrt(d) );
      num_scr2 = std::ceil( d/log(d) );
    }
  }

  for( unsigned int index = 0; index < phenos.n_rows; index++){
    arma::mat x0 = csx;
    x0.shed_col(phenos(index) - 1);
    arma::mat y = csx.unsafe_col(phenos(index) - 1);
    arma::vec xy0 = x0.t() * y;
    arma::mat beta(d-1,1);
    arma::uvec order0 = sort_index(xy0,"descend");
    arma::uvec idx_scr = order0;
    //unsigned int num_scr = idx_scr.n_rows;
    arma::uvec idx_scr1 = order0.subvec(0,num_scr1-1);
    arma::uvec idx_scr2 = order0.subvec(0,num_scr2-1);
    arma::mat x1 = x0.cols(idx_scr);
    arma::mat xx = x1.t() * x1;
    //double gamma = sum(abs(xx)).max();
    arma::vec beta0;
    beta0.fill(0);


  }
  return(csx);
}

void lasso_ladm_scr(const arma::mat & xy0, const arma::mat & x0, const arma::mat & xx0, double & T,
                    arma::mat & xy, arma::mat & x, arma::mat & xx, double rho, int max_ite,
                    double prec, const arma::vec & beta, const arma::uvec & idx_scr,
                    double lambda, bool flag ){
  arma::vec beta0, beta1;
  unsigned int n = idx_scr.n_cols;
  beta0 = beta.rows(idx_scr);


  arma::uvec idx_a(n), idx_i(n), idx_a1(n), idx_i1(n);
  idx_i.fill(1);
  idx_i.elem( find(beta0!=0)).zeros();

  for( unsigned int j = 0,size_a = 0; j < n; j++){
    if ( beta0[j] != 0){
      idx_a[size_a] = j;
      size_a++;
    }
  }

  x = x0.cols(idx_scr);
  xy = xy0.cols(idx_scr);
  if(flag) {
    xx = xx0;
    T = accu(abs(xx))*0.8;
  }
}

void lineaization_lasso(arma::vec & beta0, arma::vec & beta1, arma::vec & beta_tild, arma::vec & g,
                        arma::uvec & idx_a1, arma::uvec & idx_i1, double T,
                        double threshold, int dim){
  int size_a1 = 0;
  for( int j = 0; j < dim; j++){
    beta_tild(j) = beta0(j) - g(j)/T;
    if( fabs(beta_tild[j]) <= threshold ){
      beta1(j) = 0;
      idx_i1(j) = 1;
    }else{
      idx_a1(size_a1) = j;
      size_a1++;
      idx_i1(j) = 0;
      if( beta_tild(j) > 0 ){
        beta1(j) = beta_tild(j) - threshold;
      }else{
        beta1(j) = beta_tild(j) + threshold;
      }
    }
  }
}
