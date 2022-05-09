#include<cgLASSO.h>

using namespace Rcpp;
using namespace arma;
using namespace std;
using namespace cgLASSO;
using namespace quic;
using namespace workingparam;

// cgLASSO with uniform L1 penalty, for comparason rather than productive use. 

// [[Rcpp::export]]
List cgLASSO_path(arma::mat X,
              arma::mat Y,
              arma::vec lambdas,
              arma::vec xis,
              int diag_penalty,
              int max_iter,
              double eps,
              int s_max_condition,
              int obj_counter_max,
              int verbose)
{
  gcgWorkingParam Worker(X,Y);
  List results = cgLASSO::cgLASSO_path<gcgWorkingParam>(Worker, lambdas, 
                xis,
                diag_penalty,max_iter,eps,
                s_max_condition,obj_counter_max,verbose);
  return results;
}
