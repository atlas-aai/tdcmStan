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
    array[J] matrix[C, C] prob_transition_class;

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
  array[n_shards] vector[0] theta;

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
  array[C] simplex[C] tau;
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

  target += sum(map_rect(sum_probs, beta, theta, xr, xi));
}
generated quantities {
  vector[J] log_lik;
  array[J] matrix[C, C] format_prob_transition_class;
  vector[J*C*C] prob_transition_class;
  array[J] matrix[A, 2] prob_resp_attr;

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
