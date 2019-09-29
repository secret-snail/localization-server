//#include <opencv2/opencv.hpp>
#include <opencv2/calib3d/calib3d.hpp>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>

#include <jlog.h>
#include <stdio.h>
#include <hoff_features.hpp>
#include <iostream>
#include <string>

#include <gaussnewtonlocalization.h>
#include <invert.h>
#include <lmlocalization.h>
#include <matmult.h>
#include <printutil.h>
#include <projectpoints.h>
#include <rodrigues.h>
#include <svd.h>
#include <trigfuncs.h>
#include <twonormsq.h>
#include <util.h>

#include <cleartext-ref/gaussnewtonlocalization.hpp>
#include <cleartext-ref/invert.hpp>
#include <cleartext-ref/lmlocalization.hpp>
#include <cleartext-ref/matmult.hpp>
#include <cleartext-ref/projectpoints.hpp>
#include <cleartext-ref/rodrigues.hpp>
#include <cleartext-ref/svd.hpp>
#include <cleartext-ref/trigfuncs.hpp>
#include <cleartext-ref/twonormsq.hpp>

#include <math.h>

using namespace emp;
using namespace std;

typedef fixed_point<int64_t, 32> baset;

int main(int argc, char** argv) {

  float f_float = 715;
  float cx_float = 354;
  float cy_float = 245;
  float _cM_float[] = {f_float, 0, cx_float, 0, f_float, cy_float, 0, 0, 1};
  baset f = 715;
  baset cx = 354;
  baset cy = 245;
  baset _cM[] = {f, 0, cx, 0, f, cy, 0, 0, 1};

  vector<cv::Point_<baset>> imagePoints;
  vector<cv::Point3_<baset>> objectPoints;
  Hoffs2DPoints(imagePoints);
  Hoffs3DPoints(objectPoints);

  MSG("\n\ntesting trigfuncs - cleartext\n");
  {
    float lulz = .2;
    float lulzsin = sin(lulz);
    float lulzcos = cos(lulz);

    cout << "libc      sin: " << lulzsin << " cos: " << lulzcos << endl;
    cout << "custom    sin: " << mysin(lulz) << " cos: " << mycos(lulz) << endl;
  }
  {
    baset lulz = .2;
    cout << "baset     sin: " << mysin(lulz) << " cos: " << mycos(lulz) << endl;
  }

  MSG("\n\ntesting sqrt - cleartext\n");
  {
    baset i = .2;
    baset o = sqrt(i);
    cout << "baset sqrt" << o << endl;
  }

  MSG("\n\ntesting matmult - cleartext\n");
  {
    float A[] = {1, 1, 1, 1,
                 1,  // 3x5
                 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    float B[] = {2, 2, 2,
                 2,  // 5x4
                 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2};
    float C[3 * 4];  // 3x4
    matmult(A, 3, 5, B, 5, 4, C);
    printVector("float matmult", C, 3 * 4);
  }
  {
    baset A[] = {1, 1, 1, 1,
                 1,  // 3x5
                 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    baset B[] = {2, 2, 2,
                 2,  // 5x4
                 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2};
    baset C[3 * 4];  // 3x4
    matmult(A, 3, 5, B, 5, 4, C);
    printVector("baset matmult", C, 3 * 4);
  }

  MSG("\n\ntesting matmult2DwTranspose - cleartext\n");
  {
    float** A = new float*[3];  // 3x5
    for (int p = 0; p < 3; p++) {
      A[p] = new float[5];
      for (int pp = 0; pp < 5; pp++) {
        A[p][pp] = 1;
      }
    }
    float** B = new float*[3];  // 3x5
    for (int p = 0; p < 3; p++) {
      B[p] = new float[5];
      for (int pp = 0; pp < 5; pp++) {
        B[p][pp] = 1;
      }
    }
    float** C = new float*[3];  // 3x3
    for (int p = 0; p < 3; p++)
      C[p] = new float[3];
    matmult2DwTranspose(A, 3, 5, false, B, 3, 5, true, C);
    printMatrix("float matmult", C, 3, 3);
  }
  {
    baset** A = new baset*[3];  // 3x5
    for (int p = 0; p < 3; p++) {
      A[p] = new baset[5];
      for (int pp = 0; pp < 5; pp++) {
        A[p][pp] = 1;
      }
    }
    baset** B = new baset*[3];  // 3x5
    for (int p = 0; p < 3; p++) {
      B[p] = new baset[5];
      for (int pp = 0; pp < 5; pp++) {
        B[p][pp] = 1;
      }
    }
    baset** C = new baset*[3];  // 3x3
    for (int p = 0; p < 3; p++)
      C[p] = new baset[3];
    matmult2DwTranspose(A, 3, 5, false, B, 3, 5, true, C);
    printMatrix("baset matmult", C, 3, 3);
  }

  MSG("\n\ntesting rodrigues - cleartext\n");
  {
    float r[3] = {.2, .4, .2};  // 3x1
    float R[3 * 4] = {};
    rodrigues(r, R);
    printVector("float rodrigues", R, 3 * 4);
  }
  {
    baset r[3] = {.2, .4, .2};  // 3x1
    baset R[3 * 4] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    rodrigues(r, R);
    printVector("baset rodrigues", R, 3 * 4);
  }

  MSG("\n\ntesting projectpoints - cleartext\n");
  {
    float P[] = {0, 0, 2, 0, 0, 10, 2, 0, 10, 2, 6, 6, 6, 2, 2, 1, 1, 1, 1, 1};
    float x[] = {1.4888731,  -0.58786786, 0.71628469,
                 0.93423641, 2.983731,    18.328608};
    float projected[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    projectPoints(P, x, _cM_float, projected, 5);
    printVector("float projection", projected, 2 * 5);
  }
  {
    baset P[] = {0, 0, 2, 0, 0, 10, 2, 0, 10, 2, 6, 6, 6, 2, 2, 1, 1, 1, 1, 1};
    baset x[] = {1.4888731,  -0.58786786, 0.71628469,
                 0.93423641, 2.983731,    18.328608};
    baset projected[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    projectPoints(P, x, _cM, projected, 5);
    printVector("baset projection", projected, 2 * 5);
  }

  MSG("testing svd - sign subfunc - cleartext\n");
  {
    float a = 2.40, b = 0;
    cout << "float sign: " << mysign(a, b) << endl;
  }
  {
    baset a = 2.40, b = 0;
    cout << "baset sign: " << mysign(a, b) << endl;
  }

  MSG("testing svd - pythag subfunc - cleartext\n");
  {
    float a = -2.40, b = 96.0;
    cout << "float pythag: " << mypythag(a, b) << endl;
  }
  {
    baset a = -2.40, b = 96.0;
    cout << "baset pythag: " << mypythag(a, b) << endl;
  }

  MSG("\n\ntesting svd - cleartext\n");
  {
    int m = 5, n = 4;
    // svd overwrites with u matrix
    float** in = new float*[m];
    for (int i = 0; i < m; i++) {
      in[i] = new float[n]();
      for (int j = 0; j < n; j++) {
        in[i][j] = i + j;
      }
    }
    in[0][0] = 0;
    in[1][0] = 0;
    in[2][0] = 0;
    in[3][0] = 0;
    in[4][0] = 0;

    float* w = new float[n]();  // aka sigma, only diag
    float** v = new float*[n];  // nxn
    for (int i = 0; i < n; i++)
      v[i] = new float[n]();

    svdcmp(in, m, n, w, v);

    printVector("w", w, n);
    printMatrix("v", v, n, n);
    printMatrix("a", in, m, n);
    for (int i = 0; i < m; i++) {
      delete[] in[i];
    }
    delete[] in;
    delete[] w;
    for (int i = 0; i < n; i++) {
      delete[] v[i];
    }
    delete[] v;
  }
  {
    int m = 5, n = 4;
    // svd overwrites with u matrix
    baset** in = new baset*[m];
    for (int i = 0; i < m; i++) {
      in[i] = new baset[n]();
      for (int j = 0; j < n; j++) {
        in[i][j] = i + j;
      }
    }
    in[0][0] = 0;
    in[1][0] = 0;
    in[2][0] = 0;
    in[3][0] = 0;
    in[4][0] = 0;

    baset* w = new baset[n]();  // aka sigma, only diag
    baset** v = new baset*[n];  // nxn
    for (int i = 0; i < n; i++)
      v[i] = new baset[n]();

    svdcmp(in, m, n, w, v);

    printVector("w", w, n);
    printMatrix("v", v, n, n);
    printMatrix("a", in, m, n);
    for (int i = 0; i < m; i++) {
      delete[] in[i];
    }
    delete[] in;
    delete[] w;
    for (int i = 0; i < n; i++) {
      delete[] v[i];
    }
    delete[] v;
  }

  //    MSG("\ntesting myinvert\n");
  //    //int m=7,n=3;
  //    //float one_d_M[] = {3, 0, 1, // toy example
  //    //                   0, 0, 0,
  //    //                   0, 4, 0,
  //    //                   0, 0, 0,
  //    //                   2, 0, 0,
  //    //                   0, 0, 3,
  //    //                   0, 0, 0};
  //    int m=12,n=6;
  //    float one_d_M[] = {-9.1552734, 149.53613, -59.509277, 18.310547,
  //    0, 3.0517578,
  //                       -164.79492, -64.086914, -47.302246,
  //                       0, 19.836426, 1.5258789,
  //                       27.46582, 77.819824, 50.354004, 22.888184,
  //                       0, 1.5258789, -56.45752, -42.724609, -53.405762,
  //                       0, 22.888184, 3.0517578,
  //                       45.776367, 79.345703, 79.345703, 21.362305, 0, 0,
  //                       -32.806396, -13.73291, -20.599365,
  //                       0, 24.414062, 4.5776367, -21.362305, 108.3374,
  //                       -93.078613, 16.784668, 0, 3.0517578, -152.58789,
  //                       -45.776367, -10.681152, 0, 19.836426, 3.0517578,
  //                       6.1035156, 39.672852, 3.0517578, 21.362305, 0, 0,
  //                       -39.672852, -21.362305, -16.784668,
  //                       0, 22.888184, 1.5258789,
  //                       24.414062, 48.828125, 27.46582, 24.414062, 0, 0,
  //                       -13.73291, 6.1035156, 12.207031,
  //                       0, 22.888184, 1.5258789};
  //
  //    cv::Mat cvM = cv::Mat(m, n, cv::DataType<float>::type, one_d_M);
  //    cv::Mat ores;
  //    invert(cvM, ores, cv::DECOMP_SVD);
  //    cout << "opencv result:\n" << ores << endl;
  //
  //    float** M = new float*[m];
  //    for(int i=0; i<m; i++) {
  //        M[i] = new float[n];
  //        for (int j=0; j<n; j++) {
  //            M[i][j] = one_d_M[i*n+j];
  //        }
  //    }
  //
  //    // res is nxm (not mxn)
  //    float** res = new float*[n];
  //    for(int i=0; i<n; i++) {
  //        res[i] = new float[m];
  //    }
  //
  //    myinvert(M, m, n, res);
  //    printMatrix("cleartext result", res, n, m);
  //
  //    cout << "testing invert\n";
  //    test_invert_circuit(party, io, M, m, n, res);
  //    printMatrix("result", res, n, m);
  //
  //
  //
  //    cout << "\ntesting 2 norm squared - cleartext\n";
  //    float tvect[] = { 3.1, 5, 2, 4, 1 };
  //    cout << twonormsq(tvect, 5) << endl;
  //
  //    cout << "testing 2 norm squared\n";
  //    test_twonormsq_circuit(party, io, tvect, 5);
  //
  //
  //
  //    cout << "\ntesting gauss newton - opencv\n";
  //    distCoeffs = cv::Mat::zeros(4,1,cv::DataType<float>::type);
  //    rvec = cv::Mat::zeros(3,1,cv::DataType<float>::type);
  //    tvec = cv::Mat::zeros(3,1,cv::DataType<float>::type);
  //    CLOCK(opencv);
  //    TIC(opencv);
  //    // OpenCV PnP method
  //    cv::solvePnP(objectPoints, imagePoints, cameraMatrix, distCoeffs, rvec,
  //    tvec, false, cv::SOLVEPNP_ITERATIVE); TOC(opencv); cout << "opencv
  //    result:" << endl; cout << rvec << endl; cout << tvec << endl << endl;
  //
  //#ifdef PPL_GN
  //    {
  //        cout << "testing gauss newton - cleartext\n";
  //        float rt[6] = {0, 0, 0, 0, 0, 0}; // initial guess
  //        CLOCK(cleartext);
  //        TIC(cleartext);
  //        gaussNewton<float>( objectPoints, imagePoints, f, cx, cy, rt);
  //        TOC(cleartext);
  //        printVector("[rotation; translation]", rt, 6);
  //
  //        cout << "testing gauss newton\n";
  //        float srt[6] = {0, 0, 0, 0, 0, 0}; // initial guess
  //        test_gaussnewton_circuit(party, io,
  //            objectPoints, imagePoints, f, cx, cy, srt);
  //        printVector("[rotation; translation]", srt, 6);
  //    }
  //#endif
  //
  //#ifdef PPL_LM
  //    {
  //        cout << "testing lm - cleartext\n";
  //        float rt[6] = {0, 0, 0, 0, 0, 0}; // initial guess
  //        CLOCK(cleartext);
  //        TIC(cleartext);
  //        lm<float>( objectPoints, imagePoints, f, cx, cy, rt);
  //        TOC(cleartext);
  //        printVector("[rotation; translation]", rt, 6);
  //
  //        cout << "testing lm\n";
  //        float srt[6] = {0, 0, 0, 0, 0, 0}; // initial guess
  //        test_lm_circuit(party, io,
  //            objectPoints, imagePoints, f, cx, cy, srt);
  //        printVector("[rotation; translation]", srt, 6);
  //    }
  //#endif

  return 0;
}
