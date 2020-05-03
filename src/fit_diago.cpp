// from gaston::lmm.diago, modify by haohao.

#include <Rcpp.h>
#include <RcppEigen.h>

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;
using namespace Eigen;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::SparseMatrix;
using Eigen::Lower;

typedef Map<MatrixXd> Map_MatrixXd;

// cette fonction badibulgue x en l'utilisant pour calculer SD
// dans le cas où on n'a plus besoin de x une fois qu'on a son inverse ça le
// fait Cette fonction n'utilise que le triangle supérieur de x... et ne remplit
// que le triangle supérieur de y !! avec eps = 0 calcule l'inverse avec eps
// petit... (1e-6 ?) pseudo inverse
void blocki(Eigen::MatrixXd &x, int x1, int n, Eigen::MatrixXd &y, int y1,
            double &log_det, double &det, double eps) {
  if (n == 1) {
    double d = (std::abs(x(x1, x1)) < eps) ? 0 : x(x1, x1);
    y(y1, y1) = (d == 0) ? 0 : 1 / d;
    det = d;
    log_det = log(d);
    return;
  }

  int m1 = n / 2;
  int m2 = n - m1;

  Block<MatrixXd> A = x.block(x1, x1, m1, m1);
  Block<MatrixXd> B = x.block(x1, x1 + m1, m1, m2);

  Block<MatrixXd> TL = y.block(y1, y1, m1, m1);
  Block<MatrixXd> BL = y.block(y1 + m1, y1, m2, m1);
  Block<MatrixXd> TR = y.block(y1, y1 + m1, m1, m2);
  Block<MatrixXd> BR = y.block(y1 + m1, y1 + m1, m2, m2);

  // BR = inverse(D)
  double log_detD, detD;
  blocki(x, x1 + m1, m2, y, y1 + m1, log_detD, detD, eps);

  // BL = inverse(D)*Bt
  BL.noalias() = BR.selfadjointView<Upper>() * B.transpose();

  // le bloc A est écrasé par SD
  A.triangularView<Upper>() -= B * BL; // on ne calcule que le triangle
                                       // supérieur car on n'utilise pas l'autre

  // TL = inverse(SD)
  double log_detSD, detSD;
  blocki(x, x1, m1, y, y1, log_detSD, detSD, eps);

  TR.noalias() = TL.selfadjointView<Upper>() * (-BL.transpose());
  // on sait que cette matrice doit être symmétrique : on ne fait que la moitié
  // des calculs
  BR.triangularView<Upper>() -= BL * TR;

  log_det = log_detD + log_detSD;
  det = detD * detSD;
}

inline void sym_inverse(Eigen::MatrixXd &X, Eigen::MatrixXd &Y, double &log_det,
                        double &det, double eps) {
  blocki(X, 0, X.rows(), Y, 0, log_det, det, eps);
  Y.triangularView<Lower>() = Y.transpose(); // symétriser
}

// structure qui hérite fun pour l'optimisation par Brent
// Normalement l'utilisation doit être
// MATRIX = MatrixXd, VECTOR = VectorXd, scalar_t = double
// MATRIX = MatrixXf, VECTOR = VectorXf, scalar_t = float
template <typename MATRIX, typename VECTOR, typename scalar_t>
class diag_likelihood {
public:
  int p, n, r;
  const MATRIX Y;
  MATRIX X;
  const MATRIX Sigma;
  VECTOR P0y;
  scalar_t v;
  MATRIX XViX_i;
  VECTOR Deltab;
  scalar_t d, log_d;
  VECTOR V0b, V0bi;
  MATRIX ViX, XViX, xtx;
  scalar_t yP0y;
  diag_likelihood(int p, const MATRIX &Y, const MATRIX &X, const VECTOR &Sigma)
      : p(p), n(Sigma.rows()), r(X.cols()), Y(Y), X(X), Sigma(Sigma) {
    Deltab = Sigma.bottomRows(n - p) - VECTOR::Ones(n - p);
    XViX_i = MATRIX(r, r);
  };

  void update(scalar_t h2) {
    V0b = h2 * Sigma.bottomRows(n - p) + (1 - h2) * VECTOR::Ones(n - p);
    V0bi = V0b.cwiseInverse();
    ViX = V0bi.asDiagonal() * X.bottomRows(n - p); //      V_{0b}^{-1} X_b
    XViX = X.bottomRows(n - p).transpose() * ViX;  // X_b' V_{0b}^{-1} X_b
    sym_inverse(XViX, XViX_i, log_d, d, 1e-5);

    P0y = V0bi.asDiagonal() * Y.bottomRows(n - p) - ViX * (XViX_i * (ViX.transpose() * Y.bottomRows(n - p)));
    yP0y = P0y.dot(Y.bottomRows(n - p).col(0));
    v = yP0y / (n - r - p);
  }

  scalar_t f(scalar_t h2) {
    update(h2);
    return 0.5 * (V0b.array().log().sum() + log_d + (n - r - p) * log(yP0y) + (n - r - p) * (1 - log((scalar_t)(n - r - p)))); // avec le terme constant
  }

  scalar_t Brent_fmax(scalar_t ax, scalar_t bx, scalar_t tol) {
    /*  c is the squared inverse of the golden ratio */
    const scalar_t c = (3. - sqrt(5.)) * .5;

    /* Local variables */
    scalar_t a, b, d, e, p, q, r, u, v, w, x;
    scalar_t t2, fu, fv, fw, fx, xm, eps, tol1, tol3;

    /*  eps is approximately the square root of the relative machine precision.
     */
    eps = DBL_EPSILON;
    tol1 = eps + 1.; /* the smallest 1.000... > 1 */
    eps = sqrt(eps);

    a = ax;
    b = bx;
    v = a + c * (b - a);
    w = v;
    x = v;

    d = 0.; /* -Wall */
    e = 0.;
    fx = f(x);
    fv = fx;
    fw = fx;
    tol3 = tol / 3.;

    /*  main loop starts here ----------------------------------- */

    for (;;) {
      xm = (a + b) * .5;
      tol1 = eps * fabs(x) + tol3;
      t2 = tol1 * 2.;

      /* check stopping criterion */

      if (fabs(x - xm) <= t2 - (b - a) * .5)
        break;
      p = 0.;
      q = 0.;
      r = 0.;
      if (fabs(e) > tol1) { /* fit parabola */
        r = (x - w) * (fx - fv);
        q = (x - v) * (fx - fw);
        p = (x - v) * q - (x - w) * r;
        q = (q - r) * 2.;
        if (q > 0.)
          p = -p;
        else
          q = -q;
        r = e;
        e = d;
      }

      if (fabs(p) >= fabs(q * .5 * r) || p <= q * (a - x) ||
          p >= q * (b - x)) { /* a golden-section step */
        if (x < xm)
          e = b - x;
        else
          e = a - x;
        d = c * e;
      } else { /* a parabolic-interpolation step */
        d = p / q;
        u = x + d;
        /* f must not be evaluated too close to ax or bx */
        if (u - a < t2 || b - u < t2) {
          d = tol1;
          if (x >= xm)
            d = -d;
        }
      }

      /* f must not be evaluated too close to x */
      if (fabs(d) >= tol1)
        u = x + d;
      else if (d > 0.)
        u = x + tol1;
      else
        u = x - tol1;

      fu = f(u);
      /*  update  a, b, v, w, and x */

      if (fu <= fx) {
        if (u < x)
          b = x;
        else
          a = x;
        v = w;
        w = x;
        x = u;
        fv = fw;
        fw = fx;
        fx = fu;
      } else {
        if (u < x)
          a = u;
        else
          b = u;
        if (fu <= fw || w == x) {
          v = w;
          fv = fw;
          w = u;
          fw = fu;
        } else if (fu <= fv || v == x || v == w) {
          v = u;
          fv = fu;
        }
      }
    }
    /* end of main loop */

    return x;
  }

  // *********** CALCUL DES BLUPS ************************
  // Attention P0y n'est que (P0y)b, les n-p dernières composantes ! (les p
  // premières sont nulles)
  void blup(scalar_t h2, VECTOR &beta, VECTOR &omega) {

    VECTOR Sigmab = Sigma.bottomRows(n - p);
    omega = h2 * Sigmab.asDiagonal() * P0y;

    VECTOR z = Y;
    z.tail(n - p) -= omega + (1 - h2) * P0y;
    // Xb' Xb
    // MATRIX xtx( MATRIX(r,r).setZero().selfadjointView<Lower>().rankUpdate(
    // X.bottomRows(n-p).transpose() ));
    xtx = MATRIX(r, r).setZero();
    SelfAdjointView<MATRIX, Lower> xtx_sa(xtx);
    xtx_sa.rankUpdate(X.bottomRows(n - p).transpose());
    MATRIX xtx0(xtx_sa);
    MATRIX xtxi(r, r); // et son inverse
    scalar_t d1, ld1;
    sym_inverse(xtx0, xtxi, d1, ld1, 1e-5); // détruit xtx0

    beta = VECTOR(r + p);
    beta.topRows(r) =
        xtxi * X.bottomRows(n - p).transpose() * z.bottomRows(n - p);
    beta.bottomRows(p) = z.topRows(p) - X.topRows(p) * beta.topRows(r);
  }
};

void min_max_h2(NumericVector Sigma, double &min_h2, double &max_h2) {
  int n = Sigma.size();

  for (int i = 0; i < n; i++) {
    double s = Sigma[i];
    if (s > 1) {
      double u = 1 / (1 - s) + 1e-6;
      min_h2 = (min_h2 > u) ? min_h2 : u;
    } else if (s < 1) {
      double u = 1 / (1 - s) - 1e-6;
      max_h2 = (max_h2 < u) ? max_h2 : u;
    }
  }
}

// [[Rcpp::export]]
List fit_diago_brent(NumericVector Y, NumericMatrix X, IntegerVector p_,
                     NumericVector Sigma, NumericMatrix U, double min_h2,
                     double max_h2, double tol, double verbose) {
  Map_MatrixXd y0(as<Map<MatrixXd>>(Y));
  Map_MatrixXd x0(as<Map<MatrixXd>>(X));
  Map<VectorXd> sigma(as<Map<VectorXd>>(Sigma));
  Map_MatrixXd u(as<Map<MatrixXd>>(U));

  MatrixXd x = u.transpose() * x0;
  MatrixXd y = u.transpose() * y0;

  List R;

  min_max_h2(Sigma, min_h2, max_h2);
  if (verbose)
    Rcout << "Optimization in interval [" << min_h2 << ", " << max_h2 << "]"
          << std::endl;

  for (int i = 0; i < p_.length(); i++) {
    int p = p_(i);
    if (verbose)
      Rcout << "Optimizing with p = " << p << "\n";

    // likelihood maximization
    double h2 = min_h2;
    diag_likelihood<MatrixXd, VectorXd, double> A(p, y, x, sigma);

    h2 = A.Brent_fmax(min_h2, max_h2, tol);

    // calcul blups transféré dans la classe diag_likelihood !...
    VectorXd beta, omega;
    A.blup(h2, beta, omega);
    double s2 = (1 - h2) * A.v, tau = h2 * A.v;

    List L;
    L["sigma2"] = s2;
    L["tau"] = tau;

    if (p_.length() > 1) {
      R.push_back(L);
    } else {
      R = L;
    }
  }
  return R;
}

//computes X'WX where W is diagonal (input w as vector)
MatrixXd xtwx(const MatrixXd& xx, const MatrixXd& ww) {
  const int n(xx.cols());
  MatrixXd AtWA(MatrixXd(n, n).setZero().
    selfadjointView<Lower>().rankUpdate(xx.adjoint() * ww.asDiagonal()));
  return (AtWA);
}

//computes X'SX where S is not diagonal (input ss as matrix)
MatrixXd xtsx(const MatrixXd& xx, const MatrixXd& ss) {
  const int n(xx.cols());
  MatrixXd AtSA(MatrixXd(n, n).setZero().
    selfadjointView<Lower>().rankUpdate(xx.adjoint() * ss));
  return (AtSA);
}

MatrixXd xtx(const MatrixXd& xx) {
  const int n(xx.cols());
  MatrixXd AtA(MatrixXd(n, n).setZero().
    selfadjointView<Lower>().rankUpdate(xx.adjoint()));
  return (AtA);
}

MatrixXd xxt(const MatrixXd& xx) {
  const int m(xx.rows());
  MatrixXd AtA(MatrixXd(m, m).setZero().
    selfadjointView<Lower>().rankUpdate(xx));
  return (AtA);
}


MatrixXd conjugate_gradient(const MatrixXd& A, const VectorXd& b, int maxit, double tol)
{

  const int n(A.cols());
  VectorXd x(n);
  VectorXd r(n);
  VectorXd p(n);
  VectorXd Ap(n);
  x.fill(0);
  
  double rsold;
  double rsnew;
  double alpha;
 
  r = b; 
  p = r;
  rsold = r.squaredNorm();
  
  for (int i = 0; i < maxit; i++) {
    Ap = A * p;
    alpha = rsold / (p.transpose() * Ap);
    x = x + (alpha * p.array()).matrix();
    r = r - (alpha * Ap.array()).matrix();
    rsnew = r.squaredNorm();
    if (sqrt(rsnew) < tol) {
      break;
    }
    p = r + ((rsnew / rsold) * p.array()).matrix();
    rsold = rsnew;
  }
  return(x);
}

// [[Rcpp::export]]
SEXP crossprodcpp(SEXP X)
{
  try {
    using Eigen::Map;
    using Eigen::MatrixXd;
    using Eigen::Lower;
    const Eigen::Map<MatrixXd> A(as<Map<MatrixXd> >(X));
    const int n(A.cols());
    MatrixXd AtA(MatrixXd(n, n).setZero().
    selfadjointView<Lower>().rankUpdate(A.adjoint()));
    return wrap(AtA);
  } catch (std::exception &ex) {
    forward_exception_to_r(ex);
  } catch (...) {
    ::Rf_error("C++ exception (unknown reason)");
  }
  return R_NilValue; //-Wall
}

// [[Rcpp::export]]
SEXP geninv(SEXP GG)
{
  try {
    using Eigen::Map;
    using Eigen::MatrixXd;
    using Eigen::Lower;
    const Eigen::Map<MatrixXd> G(as<Map<MatrixXd> >(GG));
    const int n(G.rows());
    const int m(G.cols());
    const int mn(std::min(n, m));
    
    bool transp(false);
    double tol(1.0e-10);
    MatrixXd A(MatrixXd(mn, mn));
    MatrixXd L(MatrixXd(mn, mn).setZero());
    
    
    
    if (n < m) {
      transp = true;
      A = xxt(G);
    } else {
      A = xtx(G);
    }

    int r = 0;
    for (int k = 0; k < mn; k++) {
      r++;
      
      if (r == 1) {
        L.block(k, r - 1, mn - k, 1) = A.block(k, k, mn - k, 1);
      } else {
        L.block(k, r - 1, mn - k, 1) = A.block(k, k, mn - k, 1) - 
                L.block(k, 0, mn - k, r - 1) * L.block(k, 0, 1, r - 1).adjoint();
      }
      
      if (L(k, r - 1) > tol) {
        L.block(k, r - 1, 1, 1) = L.block(k, r - 1, 1, 1).array().sqrt();
        if (k + 1 < mn) {
          L.block(k + 1, r - 1, mn - k - 1, 1) = L.block(k + 1, r - 1, mn - k - 1, 1) / L(k, r - 1);
        }
      } else {
        r--;
      }
    }

    MatrixXd M(MatrixXd(r, r));
    M = xtx(L.block(0, 0, mn, r)).inverse();

    MatrixXd Y(MatrixXd(m, n));
    
    if (transp) {
      Y = G.adjoint() * L.block(0, 0, mn, r) * M * M * L.block(0, 0, mn, r).adjoint();
    } else {
      Y = L.block(0, 0, mn, r) * M * M * L.block(0, 0, mn, r).adjoint() * G.adjoint();
    }

    return wrap(Y);
  } catch (std::exception &ex) {
    forward_exception_to_r(ex);
  } catch (...) {
    ::Rf_error("C++ exception (unknown reason)");
  }
  return R_NilValue; //-Wall
}

/*** R
lmm.diago.brent <- function(Y, X = matrix(1, nrow = length(Y)), eigenK, p = 0, min_h2 = 0, max_h2 = 1, verbose = TRUE, tol = .Machine$double.eps ^ 0.25) {
    if (any(is.na(Y))) {
        stop('Missing data in Y.')
    }
    if (!is.null(X)) {
        if (length(Y) < (ncol(X) + max(p))) {
            stop('The total number of covariables and PCs as fixed effects cannot exceed the number of observations.')
        }
        if (length(Y) != nrow(X)) {
            stop('Length of Y and the number of rows of X differ.')
        }
    } else if (length(Y) < max(p)) {
        stop('The total number of covariables and PCs as fixed effects cannot exceed the number of observations.')
    }
    Sigma <- eigenK$values
    if (length(Y) != length(Sigma)) {
        stop('Length of Y and number of eigenvalues differ.')
    }
    w <- which(Sigma < 1e-6)
    Sigma[w] <- 1e-6
        
    U <- eigenK$vectors
            
    return(fit_diago_brent(Y, X, p, Sigma, U, min_h2, max_h2, tol, verbose))
}
*/
