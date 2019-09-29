#include <jlog.h>
#include <privacyconf.h>
#include <stdio.h>
#include <svd.h>
#include <iostream>
#include <string>

#include "emp-sh2pc/emp-sh2pc.h"
#include "emp-tool/emp-tool.h"

#include <math.h>
#include <util.h>

static int iminarg1, iminarg2;
#define IMIN(a, b)                 \
  (iminarg1 = (a), iminarg2 = (b), \
   (iminarg1 < (iminarg2) ? (iminarg1) : iminarg2))

using namespace emp;
using namespace std;

// EMP NOTE - new value = falseValue.If(bool, trueValue)

void BuildSignCircuit(Float* a, Float* b) {
  // WARNING:
  // cannot simply set msb of a to msb of b because
  // of the case when b is 0.0
  Float zero_gate = Float(0.0, PUBLIC);

  // share* cmp = c->PutFPGate(b, zero_gate, CMP);
  Bit bgtzero = zero_gate.less_than(*b);
  // share* tmp = c->PutINVGate(cmp);
  // a->set_wire_id(31, tmp->get_wire_id(0));
  a->value[31] = !bgtzero;
}

Float BuildPythagCircuit(Float* a, Float* b) {
  Float aa = a->sqr();
  Float bb = b->sqr();
  Float sum = aa + bb;
  Float res = sum.sqrt();
  return res;
}

int BuildSvdCircuit(Float** a, int nRows, int nCols, Float* w, Float** v) {
  Float zero_gate = Float(0.0, PUBLIC);
  // Float one_gate = Float(1.0, PUBLIC);
  // Float two_gate = Float(2.0, PUBLIC);

  int flag, i, its, j, jj, k, l, nm;
  Float anorm = Float(0.0, PUBLIC);
  Float c = Float(0.0, PUBLIC);
  Float f = Float(0.0, PUBLIC);
  Float g = Float(0.0, PUBLIC);
  Float h = Float(0.0, PUBLIC);
  Float s = Float(0.0, PUBLIC);
  Float scale = Float(0.0, PUBLIC);
  Float x = Float(0.0, PUBLIC);
  Float y = Float(0.0, PUBLIC);
  Float z = Float(0.0, PUBLIC);
  Float* rv1 = (Float*)malloc(sizeof(Float) * nCols);
  if (rv1 == NULL) {
    printf("BuildSvdCircuit(): Unable to allocate vector\n");
    return (-1);
  }

  // Householder reduction to bidiagonal form
  // g = scale = anorm = 0.0;
  for (i = 0; i < nCols; i++) {
    l = i + 1;
    rv1[i] = scale * g;
    g = Float(0.0, PUBLIC);
    s = Float(0.0, PUBLIC);
    scale = Float(0.0, PUBLIC);
    if (i < nRows) {
      for (k = i; k < nRows; k++) {
        scale = scale + a[k][i].abs();
      }
      /*if (scale)*/ {
        Bit ifscale = scale.equal(zero_gate);
        for (k = i; k < nRows; k++) {
          Float temp = a[k][i] / scale;
          a[k][i] = temp.If(ifscale, a[k][i]);
          s = s + a[k][i].sqr();
          // s = s + (a[k][i] * a[k][i]);
        }
        // if (i == 2) return 0; // good here
        // if (i==2) cout << "a[i][i] " << a[i][i].reveal<double>() << endl;
        // good
        f = a[i][i];  // no .If okay
        Float temp = s.sqrt();
        BuildSignCircuit(&temp, &f);
        temp = -temp;
        g = temp.If(ifscale, g);
        h = f * g - s;  // no .If okay
        temp = f - g;
        a[i][i] = temp.If(ifscale, a[i][i]);
        // if (i==2) cout << "f " << f.reveal<double>() << endl; good
        // if (i==2) cout << "g " << g.reveal<double>() << endl; good
        // if (i==2) cout << "h " << h.reveal<double>() << endl; good
        // if (i==2) cout << "a[i][i] " << a[i][i].reveal<double>() << endl;
        // good
        for (j = l; j < nCols; j++) {
          s = Float(0.0);
          for (k = i; k < nRows; k++) {
            // if (i==2) cout << "a[k][i] " << a[k][i].reveal<double>() << endl;
            // if (i==2) cout << "a[k][j] " << a[k][j].reveal<double>() << endl;
            s = s + (a[k][i] * a[k][j]);
            // if (i==2) cout << "mid s " << s.reveal<double>() << endl;
          }
          f = s / h;
          for (k = i; k < nRows; k++) {
            Float temp = a[k][j] + (f * a[k][i]);
            a[k][j] = temp.If(ifscale, a[k][j]);
          }
        }
        for (k = i; k < nRows; k++) {
          Float temp = a[k][i] * scale;
          a[k][i] = temp.If(ifscale, a[k][i]);
        }
      }
    }
    // if (i == 2) return 0; // bad here
    //  discrepancy in A here because of underflow?
    w[i] = scale * g;
    g = Float(0.0, PUBLIC);
    s = Float(0.0, PUBLIC);
    scale = Float(0.0, PUBLIC);
    if (i < nRows && i != nCols - 1) {
      for (k = l; k < nCols; k++)
        scale = scale + a[i][k].abs();
      /*if (scale)*/ {
        Bit ifscale = scale.equal(zero_gate);
        for (k = l; k < nCols; k++) {
          Float temp = a[i][k] / scale;
          a[i][k] = temp.If(ifscale, a[i][k]);
          s = s + a[i][k].sqr();
          // s = s + (a[i][k] * a[i][k]);
        }
        f = a[i][l];  // no .If okay
        Float temp = s.sqrt();
        BuildSignCircuit(&temp, &f);
        temp = -temp;
        g = temp.If(ifscale, g);
        temp = f * g - s;
        h = temp.If(ifscale, Float(1.0, PUBLIC));  // h is used to set rv1 if
            // scale==0, set to 1 for noop
        temp = f - g;
        a[i][l] = temp.If(ifscale, a[i][l]);
        for (k = l; k < nCols; k++)
          rv1[k] = a[i][k] / h;  // no mux because h is one ifscale
        for (j = l; j < nRows; j++) {
          for (s = Float(0.0, PUBLIC), k = l; k < nCols; k++)
            s = s + (a[j][k] * a[i][k]);
          for (k = l; k < nCols; k++) {
            Float temp = a[j][k] + (s * rv1[k]);
            a[j][k] = temp.If(ifscale, a[j][k]);
          }
        }
        for (k = l; k < nCols; k++) {
          Float temp = a[i][k] * scale;
          a[i][k] = temp.If(ifscale, a[i][k]);
        }
      }
    }
    // anorm = FMAX(anorm, (fabs(w[i]) + fabs(rv1[i])));
    Float temp = w[i].abs() + rv1[i].abs();
    Bit cmp = anorm.less_than(temp);
    anorm = anorm.If(cmp, temp);

    printf(".");
    fflush(stdout);
  }

  // accumulation of right hand transformations
  for (i = nCols - 1; i >= 0; i--) {
    if (i < nCols - 1) {
      /*if (g)*/ {
        Bit ifg = g.equal(zero_gate);
        for (j = l; j < nCols; j++) {
          Float temp = (a[i][j] / a[i][l]) / g;
          v[j][i] = temp.If(ifg, v[j][i]);
        }
        for (j = l; j < nCols; j++) {
          s = Float(0.0, PUBLIC);
          for (k = l; k < nCols; k++)
            s = s + (a[i][k] * v[k][j]);
          for (k = l; k < nCols; k++) {
            Float temp = v[k][j] + (s * v[k][i]);
            v[k][j] = temp.If(ifg, v[k][j]);
          }
        }
      }
      for (j = l; j < nCols; j++) {
        v[i][j] = Float(0.0, PUBLIC);
        v[j][i] = Float(0.0, PUBLIC);
      }
    }
    v[i][i] = Float(1.0, PUBLIC);
    g = rv1[i];
    l = i;
    printf(":");
    fflush(stdout);
  }

  // accumulation of left hand transformations
  for (i = IMIN(nRows, nCols) - 1; i >= 0; i--) {
    l = i + 1;
    g = w[i];
    for (j = l; j < nCols; j++)
      a[i][j] = Float(0.0, PUBLIC);
    /*if (g)*/ {
      Bit ifg = g.equal(zero_gate);
      g = Float(1.0, PUBLIC) / g;
      for (j = l; j < nCols; j++) {
        s = Float(0.0, PUBLIC);
        for (k = l; k < nRows; k++) {
          s = s + (a[k][i] * a[k][j]);  // no .If okay
        }
        f = (s / a[i][i]) * g;  // no .If okay
        for (k = i; k < nRows; k++) {
          Float temp = a[k][j] + (f * a[k][i]);
          a[k][j] = temp.If(ifg, a[k][j]);
        }
      }
      for (j = i; j < nRows; j++) {
        Float temp = a[j][i] * g;
        a[j][i] = temp.If(ifg, Float(0.0, PUBLIC));
      }
      /*} else {  else  - done in for loop above with mux */
      // for (j = i; j < nRows; j++) {
      //     a[j][i] = a[j][i].If(ifg, zero_gate);
      // }
    }
    a[i][i] = a[i][i] + Float(1.0, PUBLIC);
    printf("|");
    fflush(stdout);
  }
  // NOTE - ignore the fact that last row may be different here.
  // Cause is floats running out of precision

  // Diagonalization of the bidiagonal form: loop over singular
  //  values and over allowed iterations
#if PPL_FLOW == PPL_FLOW_DO
  for (k = nCols - 1; k >= 0; k--) {
    Bit converged = Bit(false, PUBLIC);
    for (its = 0; its < 30; its++) {
      cout << k << flush;
      Bit flag = Bit(true, PUBLIC);
      Bit l_found = Bit(false, PUBLIC);
      Integer priv_l = Integer(32, 0, PUBLIC);

      for (l = k; l >= 0; l--) {  // test for splitting
        nm = l - 1;               // note rv1[0] is always zero

        Integer this_l = Integer(32, l, PUBLIC);

        Bit temp = (rv1[l] + anorm).equal(anorm);
        flag = flag & !(temp & !l_found);
        priv_l = priv_l.If(temp & (!l_found), this_l);
        l_found = l_found | temp;

        if (nm >= 0) {
          temp = (w[nm].abs() + anorm).equal(anorm);  // TODO need abs here?
          priv_l = priv_l.If(temp & (!l_found), this_l);
          l_found = l_found | temp;
        }
      }

      Integer priv_nm = priv_l - Integer(32, 1, PUBLIC);
      {  // if(flag) - cancellation of rv1(l), note l is always >= 1
        c = Float(0.0, PUBLIC);         // no If needed
        s = Float(1.0, PUBLIC);         // no If needed
        for (int i = 1; i <= k; i++) {  // start at lowest l (which is 1)
          Bit igel = Integer(32, i, PUBLIC).geq(priv_l);
          Bit modify = flag & igel & (!converged);

          f = s * rv1[i];
          rv1[i] = rv1[i].If(modify, c * rv1[i]);

          Bit breaked = (f.abs() + anorm).equal(anorm);
          modify = modify & (!breaked);

          g = w[i];
          h = BuildPythagCircuit(&f, &g);
          w[i] = w[i].If(modify, h);
          h = Float(1.0, PUBLIC) / h;
          c = c.If(modify, g * h);
          s = s.If(modify, -f * h);
          for (j = 0; j < nRows; j++) {
            for (int nm_finder = 0; nm_finder < k; nm_finder++) {
              Bit found_nm = priv_nm.equal(Integer(32, nm_finder, PUBLIC));
              Bit mafnm = modify & found_nm;
              y = a[j][nm_finder];
              z = a[j][i];
              a[j][nm_finder] = a[j][nm_finder].If(mafnm, y * c + z * s);
              a[j][i] = a[j][i].If(mafnm, z * c - y * s);
            }
          }
        }
      }
      z = w[k];
      Bit newconverged = priv_l.equal(Integer(32, k, PUBLIC));
      /*if (l == k)*/ {     // convergence
        /*if (z < 0.0)*/ {  // singular value is made nonnegative
          Bit ifzlt = z.less_than(zero_gate);
          Bit modify = ifzlt & newconverged & (!converged);
          w[k] = w[k].If(modify, -z);
          for (j = 0; j < nCols; j++) {
            v[j][k] = v[j][k].If(modify, -v[j][k]);
          }
        }
      }
      converged = converged | newconverged;
      if (k == 0)
        break;  // if k == 0 we are done, code below segfaults from nm

      // shift from bottom 2-by-2 minor
      for (int l_finder = 0; l_finder < k; l_finder++) {
        Bit found_l = priv_l.equal(Integer(32, l_finder, PUBLIC));
        x = x.If(found_l & (!converged), w[l_finder]);
      }
      nm = k - 1;
      y = w[nm];
      g = rv1[nm];
      h = rv1[k];
      f = ((y - z) * (y + z) + (g - h) * (g + h)) /
          (Float(2.0, PUBLIC) * h * y);
      Float one_gate = Float(1.0, PUBLIC);
      g = BuildPythagCircuit(&f, &one_gate);
      BuildSignCircuit(&g, &f);  // modifies g
      f = ((x - z) * (x + z) + h * ((y / (f + g)) - h)) / x;
      // Next QR transformation
      c = Float(1.0, PUBLIC);
      s = Float(1.0, PUBLIC);

      for (j = 0; j <= nm; j++) {
        Integer priv_j = Integer(32, j, PUBLIC);
        Bit jgtl = priv_j.geq(priv_l);
        Bit modify = jgtl & (!converged);

        i = j + 1;
        g = rv1[i];
        y = w[i];
        h = s * g;
        g = c * g;
        z = BuildPythagCircuit(&f, &h);
        rv1[j] = rv1[j].If(modify, z);
        c = c.If(modify, f / z);
        s = s.If(modify, h / z);
        f = f.If(modify, x * c + g * s);
        g = g * c - x * s;
        h = y * s;
        y = y * c;
        for (jj = 0; jj < nCols; jj++) {
          x = x.If(modify, v[jj][j]);
          // z = z.If(modify, v[jj][i]);
          // x = v[jj][j];
          z = v[jj][i];
          // if (k==1 && modify.reveal<bool>()) {
          //     printf("v[%d][%d] = %d\n", jj, j, x.reveal<double>());
          //     printf("v[%d][%d] = %d\n", jj, i, z.reveal<double>());
          // }
          v[jj][j] = v[jj][j].If(modify, x * c + z * s);
          v[jj][i] = v[jj][i].If(modify, z * c - x * s);
        }
        z = BuildPythagCircuit(&f, &h);
        w[j] = w[j].If(modify, z);  // rotation can be artibitrary if z=0
        /*if (z)*/ {
          Bit ifz = !(z.equal(zero_gate));
          Bit ifzam = ifz & modify;
          z = Float(1.0, PUBLIC) / z;  // no If needed
          c = c.If(ifzam, f * z);
          s = s.If(ifzam, h * z);
        }
        f = f.If(modify, c * g + s * y);
        x = x.If(modify, c * y - s * g);

        for (jj = 0; jj < nRows; jj++) {
          y = a[jj][j];
          z = a[jj][i];
          a[jj][j] = a[jj][j].If(modify, y * c + z * s);
          a[jj][i] = a[jj][i].If(modify, z * c - y * s);
        }
      }
      for (int l_finder = 0; l_finder < k; l_finder++) {
        Bit found_l = priv_l.equal(Integer(32, l_finder, PUBLIC));
        rv1[l_finder] =
            rv1[l_finder].If(found_l & (!converged), Float(0.0, PUBLIC));
      }
      rv1[k] = rv1[k].If(!converged, f);
      w[k] = w[k].If(!converged, x);
    }
    printf("/");
    fflush(stdout);
  }

#elif PPL_FLOW == PPL_FLOW_LOOP_LEAK
  for (k = nCols - 1; k >= 0; k--) {
    for (its = 0; its < 30; its++) {
      cout << k;
      flag = 1;
      for (l = k; l >= 0; l--) {  // test for splitting
        nm = l - 1;               // note rv1(1) is always zero
        Bit temp = (rv1[l] + anorm).equal(anorm);
        if (temp.reveal<bool>(PUBLIC)) {
          flag = 0;
          break;
        }
        assert(nm >= 0);  // sanity check, should never happen since rv1[0] = 0
        temp = (w[nm].abs() + anorm).equal(anorm);
        if (temp.reveal<bool>(PUBLIC)) {
          break;
        }
      }

      if (flag) {  // cancellation of rv1(l), if l > 1
        c = Float(0.0, PUBLIC);
        s = Float(1.0, PUBLIC);
        for (i = l; i <= k; i++) {
          f = s * rv1[i];
          rv1[i] = c * rv1[i];
          Bit temp = (f.abs() + anorm).equal(anorm);
          if (temp.reveal<bool>(PUBLIC)) {
            break;
          }
          g = w[i];
          h = BuildPythagCircuit(&f, &g);
          w[i] = h;
          h = Float(1.0, PUBLIC) / h;
          c = g * h;
          s = -f * h;
          for (j = 0; j < nRows; j++) {
            y = a[j][nm];
            z = a[j][i];
            a[j][nm] = y * c + z * s;
            a[j][i] = z * c - y * s;
          }
        }
      }
      z = w[k];
      if (l == k) {         // convergence
        /*if (z < 0.0)*/ {  // singular value is made nonnegative
          Bit ifzlt = z.less_than(zero_gate);
          Float nz = -z;
          w[k] = w[k].If(ifzlt, nz);
          for (j = 0; j < nCols; j++) {
            Float nv = -v[j][k];
            v[j][k] = v[j][k].If(ifzlt, nv);
          }
        }
        break;
      }
      if (its == 29)
        printf("no convergence in 30 svdcmp iterations\n");
      // shift from bottom 2-by-2 minor
      x = w[l];
      nm = k - 1;
      y = w[nm];
      g = rv1[nm];
      h = rv1[k];
      f = ((y - z) * (y + z) + (g - h) * (g + h)) /
          (Float(2.0, PUBLIC) * h * y);
      Float one_gate = Float(1.0, PUBLIC);
      g = BuildPythagCircuit(&f, &one_gate);
      BuildSignCircuit(&g, &f);  // modifies g
      f = ((x - z) * (x + z) + h * ((y / (f + g)) - h)) / x;
      // Next QR transformation
      c = Float(1.0, PUBLIC);
      s = Float(1.0, PUBLIC);
      for (j = l; j <= nm; j++) {
        i = j + 1;
        g = rv1[i];
        y = w[i];
        h = s * g;
        g = c * g;
        z = BuildPythagCircuit(&f, &h);
        rv1[j] = z;
        c = f / z;
        s = h / z;
        f = x * c + g * s;
        g = g * c - x * s;
        h = y * s;
        y = y * c;
        for (jj = 0; jj < nCols; jj++) {
          x = v[jj][j];
          z = v[jj][i];
          v[jj][j] = x * c + z * s;
          v[jj][i] = z * c - x * s;
        }
        z = BuildPythagCircuit(&f, &h);
        w[j] = z;  // rotation can be artibitrary if z=0
        /*if (z)*/ {
          Bit ifz = z.equal(zero_gate);
          Float temp = Float(1.0, PUBLIC) / z;
          z = temp.If(ifz, z);
          temp = f * z;
          c = temp.If(ifz, c);
          temp = h * z;
          s = temp.If(ifz, s);
        }
        f = c * g + s * y;
        x = c * y - s * g;
        for (jj = 0; jj < nRows; jj++) {
          y = a[jj][j];
          z = a[jj][i];
          a[jj][j] = y * c + z * s;
          a[jj][i] = z * c - y * s;
        }
      }
      rv1[l] = Float(0.0, PUBLIC);
      rv1[k] = f;
      w[k] = x;
    }
    printf("/");
    fflush(stdout);
  }

#elif PPL_FLOW == PPL_FLOW_SiSL
  for (k = nCols - 1; k >= 0; k--) {
    for (its = 0; its < 2; its++) {
      cout << k;
      l = 0;
      nm = -1;
      flag = 0;

      z = w[k];
      if (l == k) {         // convergence
        /*if (z < 0.0)*/ {  // singular value is made nonnegative
          Bit ifzlt = z.less_than(zero_gate);
          Float nz = -z;
          w[k] = w[k].If(ifzlt, nz);
          for (j = 0; j < nCols; j++) {
            Float nv = -v[j][k];
            v[j][k] = v[j][k].If(ifzlt, nv);
          }
        }
        break;
      }

      if (k == 0)
        break;

      // shift from bottom 2-by-2 minor
      x = w[l];
      nm = k - 1;
      y = w[nm];
      g = rv1[nm];
      h = rv1[k];
      f = ((y - z) * (y + z) + (g - h) * (g + h)) /
          (Float(2.0, PUBLIC) * h * y);
      Float one_gate = Float(1.0, PUBLIC);
      g = BuildPythagCircuit(&f, &one_gate);
      BuildSignCircuit(&g, &f);  // modifies g
      f = ((x - z) * (x + z) + h * ((y / (f + g)) - h)) / x;
      // Next QR transformation
      c = Float(1.0, PUBLIC);
      s = Float(1.0, PUBLIC);
      for (j = l; j <= nm; j++) {
        i = j + 1;
        g = rv1[i];
        y = w[i];
        h = s * g;
        g = c * g;
        z = BuildPythagCircuit(&f, &h);
        rv1[j] = z;
        c = f / z;
        s = h / z;
        f = x * c + g * s;
        g = g * c - x * s;
        h = y * s;
        y = y * c;
        for (jj = 0; jj < nCols; jj++) {
          x = v[jj][j];
          z = v[jj][i];
          v[jj][j] = x * c + z * s;
          v[jj][i] = z * c - x * s;
        }
        z = BuildPythagCircuit(&f, &h);
        w[j] = z;  // rotation can be artibitrary if z=0
        /*if (z)*/ {
          Bit ifz = z.equal(zero_gate);
          Float temp = Float(1.0, PUBLIC) / z;
          z = temp.If(ifz, z);
          temp = f * z;
          c = temp.If(ifz, c);
          temp = h * z;
          s = temp.If(ifz, s);
        }
        f = c * g + s * y;
        x = c * y - s * g;
        for (jj = 0; jj < nRows; jj++) {
          y = a[jj][j];
          z = a[jj][i];
          a[jj][j] = y * c + z * s;
          a[jj][i] = z * c - y * s;
        }
      }
      rv1[l] = Float(0.0, PUBLIC);
      rv1[k] = f;
      w[k] = x;
    }
    printf("/");
    fflush(stdout);
  }
#endif
  printf("\n");

  free(rv1);
  return 0;
}
