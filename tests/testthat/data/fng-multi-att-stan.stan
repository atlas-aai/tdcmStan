data {
  int<lower=1> I;                      // number of items
  int<lower=1> J;                      // number of respondents
  int<lower=1> N;                      // number of observations
  int<lower=1> C;                      // number of classes
  int<lower=1> A;                      // number of attributes
  array[N, 2] int<lower=1,upper=I> ii; // item for obs n
  array[N, 2] int<lower=0,upper=1> y;  // score for obs n
  array[J, 2] int<lower=1,upper=N> s;  // starting row for j
  array[J, 2] int<lower=1,upper=I> l;  // number of items for j
  matrix[C,A] Alpha;                   // attribute pattern for each C
}
parameters {
  array[C] simplex[C] tau;
  simplex[C] Vc;
  real l_0;
  real<lower=0> l_1;
  real<lower=-1 * fmin(l4_11, l4_12)> l4_212;
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
  pi[1,2] = inv_logit(l_0+l_1);
  pi[2,2] = inv_logit(l_0+l_1);
  pi[3,2] = inv_logit(l_0+l_1);
  pi[4,2] = inv_logit(l_0+l_1);
  pi[5,2] = inv_logit(l_0);
  pi[6,2] = inv_logit(l_0);
  pi[7,2] = inv_logit(l_0);
  pi[8,2] = inv_logit(l_0);
  pi[1,3] = inv_logit(l_0);
  pi[2,3] = inv_logit(l_0);
  pi[3,3] = inv_logit(l_0);
  pi[4,3] = inv_logit(l_0+l_1);
  pi[5,3] = inv_logit(l_0+l_1);
  pi[6,3] = inv_logit(l_0+l_1);
  pi[7,3] = inv_logit(l_0+l_1);
  pi[8,3] = inv_logit(l_0+l_1);
  pi[1,4] = inv_logit(l_0+l_1);
  pi[2,4] = inv_logit(l_0+l_1);
  pi[3,4] = inv_logit(l_0+l_1);
  pi[4,4] = inv_logit(l_0+l_1+l_1+l4_212);
  pi[5,4] = inv_logit(l_0+l_1);
  pi[6,4] = inv_logit(l_0+l_1);
  pi[7,4] = inv_logit(l_0+l_1);
  pi[8,4] = inv_logit(l_0+l_1);
}
model {
  array[C, C] real ps;

  // Priors
  l_0 ~ normal(0, 2);
  l_1 ~ lognormal(0, 1);
  l4_212 ~ normal(0, 2);

  // Likelihood
  for (j in 1:J) {
    vector[C] tmp;
    for (c1 in 1:C) {
      for (c2 in 1:C) {
        array[l[j, 1]] real log_items;
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
  array[J] matrix[C, C] prob_transition_class;
  array[J] matrix[A, 2] prob_resp_attr;

  // Likelihood
  for (j in 1:J) {
    vector[C] tmp;
    array[C, C] real ps;
    for (c1 in 1:C) {
      for (c2 in 1:C) {
        array[l[j, 1]] real log_items;
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
        array[l[j, 1]] real log_items;
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
