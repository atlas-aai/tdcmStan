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
  real l1_0;
  real l2_0;
  real l3_0;
  real l4_0;
  real l5_0;
  real l6_0;
  real l7_0;
  real l8_0;
  real l9_0;
  real l10_0;
  real l11_0;
  real l12_0;
  real l13_0;
  real l14_0;
  real l15_0;
  real l16_0;
  real l17_0;
  real l18_0;
  real l19_0;
  real l20_0;
  real l21_0;
  real l22_0;
  real l23_0;
  real l24_0;
  real<lower=0> l1_11;
  real<lower=0> l2_11;
  real<lower=0> l3_11;
  real<lower=0> l4_11;
  real<lower=0> l5_11;
  real<lower=0> l6_11;
  real<lower=0> l7_11;
  real<lower=0> l8_11;
  real<lower=0> l9_12;
  real<lower=0> l10_12;
  real<lower=0> l11_12;
  real<lower=0> l12_12;
  real<lower=0> l13_12;
  real<lower=0> l14_12;
  real<lower=0> l15_12;
  real<lower=0> l16_12;
  real<lower=0> l17_13;
  real<lower=0> l18_13;
  real<lower=0> l19_13;
  real<lower=0> l20_13;
  real<lower=0> l21_13;
  real<lower=0> l22_13;
  real<lower=0> l23_13;
  real<lower=0> l24_13;
}
transformed parameters {
  matrix[I,C] pi;

  pi[1,1] = inv_logit(l1_0);
  pi[2,1] = inv_logit(l2_0);
  pi[3,1] = inv_logit(l3_0);
  pi[4,1] = inv_logit(l4_0);
  pi[5,1] = inv_logit(l5_0);
  pi[6,1] = inv_logit(l6_0);
  pi[7,1] = inv_logit(l7_0);
  pi[8,1] = inv_logit(l8_0);
  pi[9,1] = inv_logit(l9_0);
  pi[10,1] = inv_logit(l10_0);
  pi[11,1] = inv_logit(l11_0);
  pi[12,1] = inv_logit(l12_0);
  pi[13,1] = inv_logit(l13_0);
  pi[14,1] = inv_logit(l14_0);
  pi[15,1] = inv_logit(l15_0);
  pi[16,1] = inv_logit(l16_0);
  pi[17,1] = inv_logit(l17_0);
  pi[18,1] = inv_logit(l18_0);
  pi[19,1] = inv_logit(l19_0);
  pi[20,1] = inv_logit(l20_0);
  pi[21,1] = inv_logit(l21_0);
  pi[22,1] = inv_logit(l22_0);
  pi[23,1] = inv_logit(l23_0);
  pi[24,1] = inv_logit(l24_0);
  pi[1,2] = inv_logit(l1_0+l1_11);
  pi[2,2] = inv_logit(l2_0+l2_11);
  pi[3,2] = inv_logit(l3_0+l3_11);
  pi[4,2] = inv_logit(l4_0+l4_11);
  pi[5,2] = inv_logit(l5_0+l5_11);
  pi[6,2] = inv_logit(l6_0+l6_11);
  pi[7,2] = inv_logit(l7_0+l7_11);
  pi[8,2] = inv_logit(l8_0+l8_11);
  pi[9,2] = inv_logit(l9_0);
  pi[10,2] = inv_logit(l10_0);
  pi[11,2] = inv_logit(l11_0);
  pi[12,2] = inv_logit(l12_0);
  pi[13,2] = inv_logit(l13_0);
  pi[14,2] = inv_logit(l14_0);
  pi[15,2] = inv_logit(l15_0);
  pi[16,2] = inv_logit(l16_0);
  pi[17,2] = inv_logit(l17_0);
  pi[18,2] = inv_logit(l18_0);
  pi[19,2] = inv_logit(l19_0);
  pi[20,2] = inv_logit(l20_0);
  pi[21,2] = inv_logit(l21_0);
  pi[22,2] = inv_logit(l22_0);
  pi[23,2] = inv_logit(l23_0);
  pi[24,2] = inv_logit(l24_0);
  pi[1,3] = inv_logit(l1_0);
  pi[2,3] = inv_logit(l2_0);
  pi[3,3] = inv_logit(l3_0);
  pi[4,3] = inv_logit(l4_0);
  pi[5,3] = inv_logit(l5_0);
  pi[6,3] = inv_logit(l6_0);
  pi[7,3] = inv_logit(l7_0);
  pi[8,3] = inv_logit(l8_0);
  pi[9,3] = inv_logit(l9_0+l9_12);
  pi[10,3] = inv_logit(l10_0+l10_12);
  pi[11,3] = inv_logit(l11_0+l11_12);
  pi[12,3] = inv_logit(l12_0+l12_12);
  pi[13,3] = inv_logit(l13_0+l13_12);
  pi[14,3] = inv_logit(l14_0+l14_12);
  pi[15,3] = inv_logit(l15_0+l15_12);
  pi[16,3] = inv_logit(l16_0+l16_12);
  pi[17,3] = inv_logit(l17_0);
  pi[18,3] = inv_logit(l18_0);
  pi[19,3] = inv_logit(l19_0);
  pi[20,3] = inv_logit(l20_0);
  pi[21,3] = inv_logit(l21_0);
  pi[22,3] = inv_logit(l22_0);
  pi[23,3] = inv_logit(l23_0);
  pi[24,3] = inv_logit(l24_0);
  pi[1,4] = inv_logit(l1_0);
  pi[2,4] = inv_logit(l2_0);
  pi[3,4] = inv_logit(l3_0);
  pi[4,4] = inv_logit(l4_0);
  pi[5,4] = inv_logit(l5_0);
  pi[6,4] = inv_logit(l6_0);
  pi[7,4] = inv_logit(l7_0);
  pi[8,4] = inv_logit(l8_0);
  pi[9,4] = inv_logit(l9_0);
  pi[10,4] = inv_logit(l10_0);
  pi[11,4] = inv_logit(l11_0);
  pi[12,4] = inv_logit(l12_0);
  pi[13,4] = inv_logit(l13_0);
  pi[14,4] = inv_logit(l14_0);
  pi[15,4] = inv_logit(l15_0);
  pi[16,4] = inv_logit(l16_0);
  pi[17,4] = inv_logit(l17_0+l17_13);
  pi[18,4] = inv_logit(l18_0+l18_13);
  pi[19,4] = inv_logit(l19_0+l19_13);
  pi[20,4] = inv_logit(l20_0+l20_13);
  pi[21,4] = inv_logit(l21_0+l21_13);
  pi[22,4] = inv_logit(l22_0+l22_13);
  pi[23,4] = inv_logit(l23_0+l23_13);
  pi[24,4] = inv_logit(l24_0+l24_13);
  pi[1,5] = inv_logit(l1_0+l1_11);
  pi[2,5] = inv_logit(l2_0+l2_11);
  pi[3,5] = inv_logit(l3_0+l3_11);
  pi[4,5] = inv_logit(l4_0+l4_11);
  pi[5,5] = inv_logit(l5_0+l5_11);
  pi[6,5] = inv_logit(l6_0+l6_11);
  pi[7,5] = inv_logit(l7_0+l7_11);
  pi[8,5] = inv_logit(l8_0+l8_11);
  pi[9,5] = inv_logit(l9_0+l9_12);
  pi[10,5] = inv_logit(l10_0+l10_12);
  pi[11,5] = inv_logit(l11_0+l11_12);
  pi[12,5] = inv_logit(l12_0+l12_12);
  pi[13,5] = inv_logit(l13_0+l13_12);
  pi[14,5] = inv_logit(l14_0+l14_12);
  pi[15,5] = inv_logit(l15_0+l15_12);
  pi[16,5] = inv_logit(l16_0+l16_12);
  pi[17,5] = inv_logit(l17_0);
  pi[18,5] = inv_logit(l18_0);
  pi[19,5] = inv_logit(l19_0);
  pi[20,5] = inv_logit(l20_0);
  pi[21,5] = inv_logit(l21_0);
  pi[22,5] = inv_logit(l22_0);
  pi[23,5] = inv_logit(l23_0);
  pi[24,5] = inv_logit(l24_0);
  pi[1,6] = inv_logit(l1_0+l1_11);
  pi[2,6] = inv_logit(l2_0+l2_11);
  pi[3,6] = inv_logit(l3_0+l3_11);
  pi[4,6] = inv_logit(l4_0+l4_11);
  pi[5,6] = inv_logit(l5_0+l5_11);
  pi[6,6] = inv_logit(l6_0+l6_11);
  pi[7,6] = inv_logit(l7_0+l7_11);
  pi[8,6] = inv_logit(l8_0+l8_11);
  pi[9,6] = inv_logit(l9_0);
  pi[10,6] = inv_logit(l10_0);
  pi[11,6] = inv_logit(l11_0);
  pi[12,6] = inv_logit(l12_0);
  pi[13,6] = inv_logit(l13_0);
  pi[14,6] = inv_logit(l14_0);
  pi[15,6] = inv_logit(l15_0);
  pi[16,6] = inv_logit(l16_0);
  pi[17,6] = inv_logit(l17_0+l17_13);
  pi[18,6] = inv_logit(l18_0+l18_13);
  pi[19,6] = inv_logit(l19_0+l19_13);
  pi[20,6] = inv_logit(l20_0+l20_13);
  pi[21,6] = inv_logit(l21_0+l21_13);
  pi[22,6] = inv_logit(l22_0+l22_13);
  pi[23,6] = inv_logit(l23_0+l23_13);
  pi[24,6] = inv_logit(l24_0+l24_13);
  pi[1,7] = inv_logit(l1_0);
  pi[2,7] = inv_logit(l2_0);
  pi[3,7] = inv_logit(l3_0);
  pi[4,7] = inv_logit(l4_0);
  pi[5,7] = inv_logit(l5_0);
  pi[6,7] = inv_logit(l6_0);
  pi[7,7] = inv_logit(l7_0);
  pi[8,7] = inv_logit(l8_0);
  pi[9,7] = inv_logit(l9_0+l9_12);
  pi[10,7] = inv_logit(l10_0+l10_12);
  pi[11,7] = inv_logit(l11_0+l11_12);
  pi[12,7] = inv_logit(l12_0+l12_12);
  pi[13,7] = inv_logit(l13_0+l13_12);
  pi[14,7] = inv_logit(l14_0+l14_12);
  pi[15,7] = inv_logit(l15_0+l15_12);
  pi[16,7] = inv_logit(l16_0+l16_12);
  pi[17,7] = inv_logit(l17_0+l17_13);
  pi[18,7] = inv_logit(l18_0+l18_13);
  pi[19,7] = inv_logit(l19_0+l19_13);
  pi[20,7] = inv_logit(l20_0+l20_13);
  pi[21,7] = inv_logit(l21_0+l21_13);
  pi[22,7] = inv_logit(l22_0+l22_13);
  pi[23,7] = inv_logit(l23_0+l23_13);
  pi[24,7] = inv_logit(l24_0+l24_13);
  pi[1,8] = inv_logit(l1_0+l1_11);
  pi[2,8] = inv_logit(l2_0+l2_11);
  pi[3,8] = inv_logit(l3_0+l3_11);
  pi[4,8] = inv_logit(l4_0+l4_11);
  pi[5,8] = inv_logit(l5_0+l5_11);
  pi[6,8] = inv_logit(l6_0+l6_11);
  pi[7,8] = inv_logit(l7_0+l7_11);
  pi[8,8] = inv_logit(l8_0+l8_11);
  pi[9,8] = inv_logit(l9_0+l9_12);
  pi[10,8] = inv_logit(l10_0+l10_12);
  pi[11,8] = inv_logit(l11_0+l11_12);
  pi[12,8] = inv_logit(l12_0+l12_12);
  pi[13,8] = inv_logit(l13_0+l13_12);
  pi[14,8] = inv_logit(l14_0+l14_12);
  pi[15,8] = inv_logit(l15_0+l15_12);
  pi[16,8] = inv_logit(l16_0+l16_12);
  pi[17,8] = inv_logit(l17_0+l17_13);
  pi[18,8] = inv_logit(l18_0+l18_13);
  pi[19,8] = inv_logit(l19_0+l19_13);
  pi[20,8] = inv_logit(l20_0+l20_13);
  pi[21,8] = inv_logit(l21_0+l21_13);
  pi[22,8] = inv_logit(l22_0+l22_13);
  pi[23,8] = inv_logit(l23_0+l23_13);
  pi[24,8] = inv_logit(l24_0+l24_13);
}
model {
  real ps[C, C];

  // Priors
  l1_0 ~ normal(0, 2);
  l2_0 ~ normal(0, 2);
  l3_0 ~ normal(0, 2);
  l4_0 ~ normal(0, 2);
  l5_0 ~ normal(0, 2);
  l6_0 ~ normal(0, 2);
  l7_0 ~ normal(0, 2);
  l8_0 ~ normal(0, 2);
  l9_0 ~ normal(0, 2);
  l10_0 ~ normal(0, 2);
  l11_0 ~ normal(0, 2);
  l12_0 ~ normal(0, 2);
  l13_0 ~ normal(0, 2);
  l14_0 ~ normal(0, 2);
  l15_0 ~ normal(0, 2);
  l16_0 ~ normal(0, 2);
  l17_0 ~ normal(0, 2);
  l18_0 ~ normal(0, 2);
  l19_0 ~ normal(0, 2);
  l20_0 ~ normal(0, 2);
  l21_0 ~ normal(0, 2);
  l22_0 ~ normal(0, 2);
  l23_0 ~ normal(0, 2);
  l24_0 ~ normal(0, 2);
  l1_11 ~ lognormal(0, 1);
  l2_11 ~ lognormal(0, 1);
  l3_11 ~ lognormal(0, 1);
  l4_11 ~ lognormal(0, 1);
  l5_11 ~ lognormal(0, 1);
  l6_11 ~ lognormal(0, 1);
  l7_11 ~ lognormal(0, 1);
  l8_11 ~ lognormal(0, 1);
  l9_12 ~ lognormal(0, 1);
  l10_12 ~ lognormal(0, 1);
  l11_12 ~ lognormal(0, 1);
  l12_12 ~ lognormal(0, 1);
  l13_12 ~ lognormal(0, 1);
  l14_12 ~ lognormal(0, 1);
  l15_12 ~ lognormal(0, 1);
  l16_12 ~ lognormal(0, 1);
  l17_13 ~ lognormal(0, 1);
  l18_13 ~ lognormal(0, 1);
  l19_13 ~ lognormal(0, 1);
  l20_13 ~ lognormal(0, 1);
  l21_13 ~ lognormal(0, 1);
  l22_13 ~ lognormal(0, 1);
  l23_13 ~ lognormal(0, 1);
  l24_13 ~ lognormal(0, 1);

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
