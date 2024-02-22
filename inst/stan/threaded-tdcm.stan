functions {
  real minmax (real x) {
    if (x < .01) {
      return 0.01;
    }

    if (x > 0.99) {
      return 0.99;
    }

    return x;
  }

  vector sum_probs(vector beta, vector theta, array[] real xr, array[] int xi) {
    int Z = num_elements(xi);
    int ys = xi[Z - 1];
    int iis = xi[Z];

    array[iis] int y1 = xi[1:iis];
    array[iis] int ii1 = xi[(iis + 1):(2 * iis)];
    array[iis] int jj1 = xi[((2 * iis) + 1):(3 * iis)];
    array[ys] int s1 = xi[((3 * iis) + 1):((3 * iis) + ys)];
    array[ys] int l1 = xi[((3 * iis) + ys + 1):((3 * iis) + (2 * ys))];

    array[iis] int y2 = xi[((3 * iis) + (2 * ys) + 1):((4 * iis) + (2 * ys))];
    array[iis] int ii2 = xi[((4 * iis) + (2 * ys) + 1):((5 * iis) + (2 * ys))];
    array[iis] int jj2 = xi[((5 * iis) + (2 * ys) + 1):((6 * iis) + (2 * ys))];
    array[ys] int s2 = xi[((6 * iis) + (2 * ys) + 1):((6 * iis) + (3 * ys))];
    array[ys] int l2 = xi[((6 * iis) + (3 * ys) + 1):((6 * iis) + (4 * ys))];

    int I = xi[Z - 7];
    int N = xi[Z - 6];
    int C = xi[Z - 5];
    int A = xi[Z - 4];
    int J = xi[Z - 3];
    int M = xi[Z - 2];

    vector[C] Vc = beta[1:C];
    vector[I * C] pic = beta[(C + 1):(C + (I * C))];
    vector[C * C] tauc = beta[(C + (I * C) + 1):(C + (I * C) + (C * C))];

    matrix[C, C] ps;
    matrix[C, C] tau_c;
    real person = 0;

    matrix[I, C] pi_c;
    for(c in 1:C) {
      for(i in 1:I) {
        int ic = i + ((c-1) * I);
        pi_c[i, c] = pic[ic];
      }
    }

    for(c1 in 1:C) {
      for(c2 in 1:C) {
        int cc = c1 + ((c2 - 1) * C);
        tau_c[c1, c2] = tauc[cc];
      }
    }

    // Likelihood
    for (j in 1:J) {
      vector[C] tmp;
      for (c1 in 1:C) {
        for (c2 in 1:C) {
          array[l1[j]] real log_items;
          for (m in 1:l1[j]) {
            int i = ii1[s1[j] + m - 1];
            log_items[m] = y1[s1[j] + m - 1] * log(pi_c[i,c1]) + (1 - y1[s1[j] + m - 1]) * log(1 - pi_c[i,c1]) + y2[s1[j] + m - 1] * log(pi_c[i,c2]) + (1 - y2[s1[j] + m - 1]) * log(1 - pi_c[i,c2]);
          }
          ps[c1, c2] = log(Vc[c1]) + log(tau_c[c1, c2]) + sum(log_items);
        }
        tmp[c1] = log_sum_exp(ps[c1,]);
      }
      person += log_sum_exp(tmp);
    }

    return [person]';
  }

  vector person_loglik(vector beta, vector theta, array[] real xr, array[] int xi) {
    int Z = num_elements(xi);
    int ys = xi[Z - 1];
    int iis = xi[Z];

    array[iis] int y1 = xi[1:iis];
    array[iis] int ii1 = xi[(iis + 1):(2 * iis)];
    array[iis] int jj1 = xi[((2 * iis) + 1):(3 * iis)];
    array[ys] int s1 = xi[((3 * iis) + 1):((3 * iis) + ys)];
    array[ys] int l1 = xi[((3 * iis) + ys + 1):((3 * iis) + (2 * ys))];

    array[iis] int y2 = xi[((3 * iis) + (2 * ys) + 1):((4 * iis) + (2 * ys))];
    array[iis] int ii2 = xi[((4 * iis) + (2 * ys) + 1):((5 * iis) + (2 * ys))];
    array[iis] int jj2 = xi[((5 * iis) + (2 * ys) + 1):((6 * iis) + (2 * ys))];
    array[ys] int s2 = xi[((6 * iis) + (2 * ys) + 1):((6 * iis) + (3 * ys))];
    array[ys] int l2 = xi[((6 * iis) + (3 * ys) + 1):((6 * iis) + (4 * ys))];

    int I = xi[Z - 7];
    int N = xi[Z - 6];
    int C = xi[Z - 5];
    int A = xi[Z - 4];
    int J = xi[Z - 3];
    int M = xi[Z - 2];

    vector[C] Vc = beta[1:C];
    vector[I * C] pic = beta[(C + 1):(C + (I * C))];
    vector[C * C] tauc = beta[(C + (I * C) + 1):(C + (I * C) + (C * C))];

    matrix[C, C] ps;
    matrix[C, C] tau_c;
    vector[J] person;

    matrix[I, C] pi_c;
    for(c in 1:C) {
      for(i in 1:I) {
        int ic = i + ((c-1) * I);
        pi_c[i, c] = pic[ic];
      }
    }

    for(c1 in 1:C) {
      for(c2 in 1:C) {
        int cc = c1 + ((c2 - 1) * C);
        tau_c[c1, c2] = tauc[cc];
      }
    }

    // Likelihood
    for (j in 1:J) {
      vector[C] tmp;
      for (c1 in 1:C) {
        for (c2 in 1:C) {
          array[l1[j]] real log_items;
          for (m in 1:l1[j]) {
            int i = ii1[s1[j] + m - 1];
            log_items[m] = y1[s1[j] + m - 1] * log(pi_c[i,c1]) + (1 - y1[s1[j] + m - 1]) * log(1 - pi_c[i,c1]) + y2[s1[j] + m - 1] * log(pi_c[i,c2]) + (1 - y2[s1[j] + m - 1]) * log(1 - pi_c[i,c2]);
          }
          ps[c1, c2] = log(Vc[c1]) + log(tau_c[c1, c2]) + sum(log_items);
        }
        tmp[c1] = log_sum_exp(ps[c1,]);
      }
      person[j] = log_sum_exp(tmp);
    }

    return person;
  }

  vector resp_transition(vector beta, vector theta, array[] real xr, array[] int xi) {
    int Z = num_elements(xi);
    int ys = xi[Z - 1];
    int iis = xi[Z];

    array[ii2] int y1 = xi[1:iis];
    array[ii2] int ii1 = xi[(iis + 1):(2 * iis)];
    array[ii2] int jj1 = xi[((2 * iis) + 1):(3 * iis)];
    array[ys] int s1 = xi[((3 * iis) + 1):((3 * iis) + ys)];
    array[ys] int l1 = xi[((3 * iis) + ys + 1):((3 * iis) + (2 * ys))];

    array[ii2] int y2 = xi[((3 * iis) + (2 * ys) + 1):((4 * iis) + (2 * ys))];
    array[ii2] int ii2 = xi[((4 * iis) + (2 * ys) + 1):((5 * iis) + (2 * ys))];
    array[ii2] int jj2 = xi[((5 * iis) + (2 * ys) + 1):((6 * iis) + (2 * ys))];
    array[ys] int s2 = xi[((6 * iis) + (2 * ys) + 1):((6 * iis) + (3 * ys))];
    array[ys] int l2 = xi[((6 * iis) + (3 * ys) + 1):((6 * iis) + (4 * ys))];

    int I = xi[Z - 7];
    int N = xi[Z - 6];
    int C = xi[Z - 5];
    int A = xi[Z - 4];
    int J = xi[Z - 3];
    int M = xi[Z - 2];

    vector[C] Vc = beta[1:C];
    vector[I * C] pic = beta[(C + 1):(C + (I * C))];
    vector[C * C] tauc = beta[(C + (I * C) + 1):(C + (I * C) + (C * C))];

    matrix[C, C] ps;
    matrix[C, C] tau_c;
    matrix[C, C] prob_transition_class[J];

    vector[J * C * C] person;

    matrix[I, C] pi_c;
    for(c in 1:C) {
      for(i in 1:I) {
        int ic = i + ((c-1) * I);
        pi_c[i, c] = pic[ic];
      }
    }

    for(c1 in 1:C) {
      for(c2 in 1:C) {
        int cc = c1 + ((c2 - 1) * C);
        tau_c[c1, c2] = tauc[cc];
      }
    }

    // latent class probabilities
    for (j in 1:J) {
      vector[C] tmp;
      matrix[C, C] prob_joint;
      for (c1 in 1:C) {
        for (c2 in 1:C) {
          array[l1[j]] real log_items;
          for (m in 1:l1[j]) {
            int i = ii1[s1[j] + m - 1];
            log_items[m] = y1[s1[j] + m - 1] * log(pi_c[i,c1]) + (1 - y1[s1[j] + m - 1]) * log(1 - pi_c[i,c1]) + y2[s1[j] + m - 1] * log(pi_c[i,c2]) + (1 - y2[s1[j] + m - 1]) * log(1 - pi_c[i,c2]);
          }
          prob_joint[c1, c2] = log(Vc[c1]) + log(tau_c[c1, c2]) + sum(log_items);
        }
      }
      prob_transition_class[j] = exp(prob_joint) / sum(exp(prob_joint));

      for (c1 in 1:C) {
        for (c2 in 1:C) {
          person[((c1 - 1) * 8) + c2 + ((j - 1) * C * C)] = prob_transition_class[j,c1,c2];
        }
      }
    }
    return person;
  }
}
data {
  int<lower=1> I;                       // number of items
  int<lower=1> J;                       // number of respondents
  int<lower=1> N;                       // number of observations
  int<lower=1> C;                       // number of classes
  int<lower=1> A;                       // number of attributes
  array[N, 2] int<lower=1,upper=I> ii;  // item for obs n
  array[N, 2] int<lower=1,upper=J> jj;  // respondent for obs n
  array[N, 2] int<lower=0,upper=1> y;   // score for obs n
  array[J, 2] int<lower=1,upper=N> s;   // starting row for j
  array[J, 2] int<lower=1,upper=I> l;   // number of items for j
  matrix[C,A] Alpha;                    // attribute pattern for each C
  int<lower=1> n_shards;                // the number of shards to split the data into
}
transformed data {
  int ys = num_elements(s) / 2 / n_shards;
  int iis = num_elements(ii) / 2 / n_shards;

  int M = iis;

  array[n_shards, (4 * ys) + (6 * iis) + 8] int xi;

  // an empty set of per-shard parameters
  vector[0] theta[n_shards];

  array[n_shards,1] real xr;
  for(kk in 1:n_shards) {
    xr[kk, 1] = 1.0;
  }

  // split into shards
  for (i in 1:n_shards) {
    int ylower;
    int yupper;
    int iilower;
    int iiupper;

    ylower = ((i - 1) * ys) + 1;
    yupper = i * ys;
    iilower = ((i - 1) * iis) + 1;
    iiupper = i * iis;

    xi[i, 1:iis] = y[iilower:iiupper, 1];
    xi[i, (iis + 1):(iis + iis)] = ii[iilower:iiupper, 1];
    xi[i, ((2 * iis) + 1):((2 * iis) + iis)] = jj[iilower:iiupper, 1];
    xi[i, ((3 * iis) + 1):((3 * iis) + ys)] = s[1:ys, 1];
    xi[i, ((3 * iis) + ys + 1):((3 * iis) + (2 * ys))] = l[ylower:yupper, 1];
    xi[i, ((3 * iis) + (2 * ys) + 1):((4 * iis) + (2 * ys))] = y[iilower:iiupper, 2];
    xi[i, ((4 * iis) + (2 * ys) + 1):((5 * iis) + (2 * ys))] = ii[iilower:iiupper, 2];
    xi[i, ((5 * iis) + (2 * ys) + 1):((6 * iis) + (2 * ys))] = jj[iilower:iiupper, 2];
    xi[i, ((6 * iis) + (2 * ys) + 1):((6 * iis) + (3 * ys))] = s[1:ys, 2];
    xi[i, ((6 * iis) + (3 * ys) + 1):((6 * iis) + (4 * ys))] = l[ylower:yupper, 2];
    xi[i, ((6 * iis) + (4 * ys) + 1)] = I;
    xi[i, ((6 * iis) + (4 * ys) + 2)] = N / n_shards;
    xi[i, ((6 * iis) + (4 * ys) + 3)] = C;
    xi[i, ((6 * iis) + (4 * ys) + 4)] = A;
    xi[i, ((6 * iis) + (4 * ys) + 5)] = J / n_shards;
    xi[i, ((6 * iis) + (4 * ys) + 6)] = iis;
    xi[i, ((6 * iis) + (4 * ys) + 7)] = ys;
    xi[i, ((6 * iis) + (4 * ys) + 8)] = iis;
  }
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
  real<lower=0> l1_11;
  real<lower=0> l2_12;
  real<lower=0> l3_11;
  real<lower=0> l4_12;
  real<lower=0> l5_11;
  real<lower=0> l5_12;
  real<lower=0> l6_11;
  real<lower=0> l6_12;
  real<lower=-1 * fmin(l5_11, l5_12)> l5_212;
  real<lower=-1 * fmin(l6_11, l6_12)> l6_212;
}
transformed parameters {
  matrix[I,C] pi;

  pi[1,1] = inv_logit(l1_0);
  pi[2,1] = inv_logit(l2_0);
  pi[3,1] = inv_logit(l3_0);
  pi[4,1] = inv_logit(l4_0);
  pi[5,1] = inv_logit(l5_0);
  pi[6,1] = inv_logit(l6_0);
  pi[1,2] = inv_logit(l1_0+l1_11);
  pi[2,2] = inv_logit(l2_0);
  pi[3,2] = inv_logit(l3_0+l3_11);
  pi[4,2] = inv_logit(l4_0);
  pi[5,2] = inv_logit(l5_0+l5_11);
  pi[6,2] = inv_logit(l6_0+l6_11);
  pi[1,3] = inv_logit(l1_0);
  pi[2,3] = inv_logit(l2_0+l2_12);
  pi[3,3] = inv_logit(l3_0);
  pi[4,3] = inv_logit(l4_0+l4_12);
  pi[5,3] = inv_logit(l5_0+l5_12);
  pi[6,3] = inv_logit(l6_0+l6_12);
  pi[1,4] = inv_logit(l1_0+l1_11);
  pi[2,4] = inv_logit(l2_0+l2_12);
  pi[3,4] = inv_logit(l3_0+l3_11);
  pi[4,4] = inv_logit(l4_0+l4_12);
  pi[5,4] = inv_logit(l5_0+l5_11+l5_12+l5_212);
  pi[6,4] = inv_logit(l6_0+l6_11+l6_12+l6_212);

  array[I * C] real pic;
  for(c in 1:C) {
    for(i in 1:I) {
      int ic = i + ((c - 1) * I);
      pic[ic] = pi[i, c];
    }
  }

  array[C * C] real tauc;
  for(c1 in 1:C) {
    for(c2 in 1:C) {
      int cc = c2 + ((c1 - 1) * C);
      tauc[cc] = tau[c2, c1];
    }
  }

  // a set of shared parameters
  vector[C + (I * C) + (C * C)] beta;
  beta[1:C] = Vc[1:C];
  beta[(C + 1):(C + (I * C))] = to_vector(pic[1:(I * C)]);
  beta[(C + (I * C) + 1):(C + (I * C) + (C * C))] = to_vector(tauc[1:(C * C)]);
}
model {
  array[C, C] real ps;

  // Priors
  l1_0 ~ normal(0, 2);
  l2_0 ~ normal(0, 2);
  l3_0 ~ normal(0, 2);
  l4_0 ~ normal(0, 2);
  l5_0 ~ normal(0, 2);
  l6_0 ~ normal(0, 2);
  l1_11 ~ lognormal(0, 1);
  l2_12 ~ lognormal(0, 1);
  l3_11 ~ lognormal(0, 1);
  l4_12 ~ lognormal(0, 1);
  l5_11 ~ lognormal(0, 1);
  l5_12 ~ lognormal(0, 1);
  l6_11 ~ lognormal(0, 1);
  l6_12 ~ lognormal(0, 1);
  l5_212 ~ normal(0, 2);
  l6_212 ~ normal(0, 2);

  target += sum(map_rect(sum_probs, beta, theta, xr, xi));
}
generated quantities {
  vector[J] log_lik;
  matrix[C, C] format_prob_transition_class[J];
  vector[J*C*C] prob_transition_class;
  matrix[A, 2] prob_resp_attr[J];

  log_lik = map_rect(person_loglik, beta, theta, xr, xi);

  prob_transition_class = map_rect(resp_transition, beta, theta, xr, xi);

  for (j in 1:J) {
    for (c1 in 1:C) {
      for (c2 in 1:C) {
        int iter = ((c1 - 1) * 8) + c2 + ((j - 1) * C * C);
        format_prob_transition_class[j,c1,c2] = prob_transition_class[iter];
      }
    }
  }

  for (j in 1:J) {
    for (a in 1:A) {
      vector[C] prob_attr_class_t1;
      vector[C] prob_attr_class_t2;
      for (c in 1:C) {
        prob_attr_class_t1[c] = sum(format_prob_transition_class[j,c,]) * Alpha[c,a];
        prob_attr_class_t2[c] = sum(format_prob_transition_class[j,,c]) * Alpha[c,a];
      }
      prob_resp_attr[j,a,1] = sum(prob_attr_class_t1);
      prob_resp_attr[j,a,2] = sum(prob_attr_class_t2);
    }
  }
}
