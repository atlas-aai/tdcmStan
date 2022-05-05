data {
  int<lower=1> I;                    // number of items
  int<lower=1> J;                    // number of respondents
  int<lower=1> N;                    // number of observations
  int<lower=1> C;                    // number of classes
  int<lower=1> A;                    // number of attributes
  int<lower=1,upper=I> ii[N, 2];     // item for obs n
  int<lower=1,upper=J> jj[N, 2];     // respondent for obs n
  int<lower=0,upper=1> y[N, 2];      // score for obs n
  int<lower=1,upper=N> s[J, 2];      // starting row for j
  int<lower=1,upper=I> l[J, 2];      // number of items for j
  matrix[C,A] Alpha;                 // attribute pattern for each C
}
parameters {
  simplex[C] tau[C];
  simplex[C] Vc;
  real l_0;
  real<lower=0> l_1;
}
transformed parameters {
  matrix[I,C] pi;

  pi[1,1] = inv_logit(l_0);
  pi[2,1] = inv_logit(l_0);
  pi[3,1] = inv_logit(l_0);
  pi[4,1] = inv_logit(l_0);
  pi[5,1] = inv_logit(l_0);
  pi[6,1] = inv_logit(l_0);
  pi[7,1] = inv_logit(l_0);
  pi[8,1] = inv_logit(l_0);
  pi[9,1] = inv_logit(l_0);
  pi[10,1] = inv_logit(l_0);
  pi[11,1] = inv_logit(l_0);
  pi[12,1] = inv_logit(l_0);
  pi[13,1] = inv_logit(l_0);
  pi[14,1] = inv_logit(l_0);
  pi[15,1] = inv_logit(l_0);
  pi[16,1] = inv_logit(l_0);
  pi[17,1] = inv_logit(l_0);
  pi[18,1] = inv_logit(l_0);
  pi[19,1] = inv_logit(l_0);
  pi[20,1] = inv_logit(l_0);
  pi[21,1] = inv_logit(l_0);
  pi[22,1] = inv_logit(l_0);
  pi[23,1] = inv_logit(l_0);
  pi[24,1] = inv_logit(l_0);
  pi[1,2] = inv_logit(l_0+l_1);
  pi[2,2] = inv_logit(l_0+l_1);
  pi[3,2] = inv_logit(l_0+l_1);
  pi[4,2] = inv_logit(l_0+l_1);
  pi[5,2] = inv_logit(l_0+l_1);
  pi[6,2] = inv_logit(l_0+l_1);
  pi[7,2] = inv_logit(l_0+l_1);
  pi[8,2] = inv_logit(l_0+l_1);
  pi[9,2] = inv_logit(l_0+l_1);
  pi[10,2] = inv_logit(l_0+l_1);
  pi[11,2] = inv_logit(l_0+l_1);
  pi[12,2] = inv_logit(l_0+l_1);
  pi[13,2] = inv_logit(l_0+l_1);
  pi[14,2] = inv_logit(l_0+l_1);
  pi[15,2] = inv_logit(l_0+l_1);
  pi[16,2] = inv_logit(l_0+l_1);
  pi[17,2] = inv_logit(l_0+l_1);
  pi[18,2] = inv_logit(l_0+l_1);
  pi[19,2] = inv_logit(l_0+l_1);
  pi[20,2] = inv_logit(l_0+l_1);
  pi[21,2] = inv_logit(l_0+l_1);
  pi[22,2] = inv_logit(l_0+l_1);
  pi[23,2] = inv_logit(l_0+l_1);
  pi[24,2] = inv_logit(l_0+l_1);
}
model {
  real ps[C, C];

  // Priors
  l_0 ~ normal(0, 2);
  l_1 ~ lognormal(0, 1);

  // Likelihood
  for (j in 1:J) {
    vector[C] tmp;
    for (c1 in 1:C) {
      for (c2 in 1:C) {
        real log_items[l[j, 1]];
        for (m in 1:l[j, 1]) {
          int i = ii[s[j, 1] + m - 1, 1];
          log_items[m] = y[s[j, 1] + m - 1, 1] * log(pi[i,c1]) + (1 - y[s[j, 1] + m - 1, 1]) * log(1 - pi[i,c1]) + y[s[j, 1] + m - 1, 2] * log(pi[i,c2]) + (1 - y[s[j, 1] + m - 1, 2]) * log(1 - pi[i,c2]);
        }
        ps[c1, c2] = log(Vc[c1]) + log(tau[c1, c2]) + sum(log_items);
      }
      tmp[c1] = log_sum_exp(ps[c1,]);
    }
    target += log_sum_exp(tmp);
  }
}
generated quantities {
  vector[J] log_lik;
  matrix[C, C] prob_transition_class[J];
  matrix[A, 2] prob_resp_attr[J];

  // Likelihood
  for (j in 1:J) {
    vector[C] tmp;
    real ps[C, C];
    for (c1 in 1:C) {
      for (c2 in 1:C) {
        real log_items[l[j, 1]];
        for (m in 1:l[j, 1]) {
          int i = ii[s[j, 1] + m - 1, 1];
          log_items[m] = y[s[j, 1] + m - 1, 1] * log(pi[i,c1]) + (1 - y[s[j, 1] + m - 1, 1]) * log(1 - pi[i,c1]) + y[s[j, 1] + m - 1, 2] * log(pi[i,c2]) + (1 - y[s[j, 1] + m - 1, 2]) * log(1 - pi[i,c2]);
        }
        ps[c1, c2] = log(Vc[c1]) + log(tau[c1, c2]) + sum(log_items);
      }
      tmp[c1] = log_sum_exp(ps[c1,]);
    }
    log_lik[j] = log_sum_exp(tmp);
  }

  // latent class probabilities
  for (j in 1:J) {
    vector[C] tmp;
    matrix[C, C] prob_joint;
    for (c1 in 1:C) {
      for (c2 in 1:C) {
        real log_items[l[j, 1]];
        for (m in 1:l[j, 1]) {
          int i = ii[s[j, 1] + m - 1, 1];
          log_items[m] = y[s[j, 1] + m - 1, 1] * log(pi[i,c1]) + (1 - y[s[j, 1] + m - 1, 1]) * log(1 - pi[i,c1]) + y[s[j, 1] + m - 1, 2] * log(pi[i,c2]) + (1 - y[s[j, 1] + m - 1, 2]) * log(1 - pi[i,c2]);
        }
        prob_joint[c1, c2] = log(Vc[c1]) + log(tau[c1, c2]) + sum(log_items);
      }
    }
    prob_transition_class[j] = exp(prob_joint) / sum(exp(prob_joint));
  }

  for (j in 1:J) {
    for (a in 1:A) {
      vector[C] prob_attr_class_t1;
      vector[C] prob_attr_class_t2;
      for (c in 1:C) {
        prob_attr_class_t1[c] = sum(prob_transition_class[j,c,]) * Alpha[c,a];
        prob_attr_class_t2[c] = sum(prob_transition_class[j,,c]) * Alpha[c,a];
      }
      prob_resp_attr[j,a,1] = sum(prob_attr_class_t1);
      prob_resp_attr[j,a,2] = sum(prob_attr_class_t2);
    }
  }
}
