// Brief: Call matlab engine to generate SG golay
//        Shows how to put & get variable from MATLAB

// NOTE: engOpen("matlab -nosplash") will open matlab
//       without start splash screen

// Matlab include
#include <iostream>
#include <stdlib.h>
#include <engine.h>
#include <string.h>

int main()
{
    std::cout << ".... Starting MATLAB engine ..." << std::endl;
    Engine *meng;
    if (!(meng = engOpen("matlab -nosplash"))) {
        std::cerr << "Can't start MATLAB engine" << std::endl;
        return EXIT_FAILURE;
    } else {
        std::cout << ".... MATLAB engine started ..." << std::endl;
    }

    const int N = 4;   // order of polynomial fit
    const int F = 11;  // window size

    // create N, F, send to MATLAB
    mxArray *mN = NULL;
    mN = mxCreateDoubleScalar(N);
    engPutVariable(meng, "N", mN);

    mxArray *mF = NULL;
    mF = mxCreateDoubleScalar(F);
    engPutVariable(meng, "F", mF);

    // call sgolay & retrieve return value
    double* ret;
    engEvalString(meng, "[b,g] = sgolay(N,F)");
    mxArray *mRet;
    mRet = engGetVariable(meng, "g");
    if (mRet == NULL) {
        std::cerr << "No g has been created" << std::endl;
    } else {
        ret = mxGetPr(mRet);
        std::cout << "g = " << ret[0] << " " << ret[5] << std::endl;
        std::cout << "row = " << mxGetM(mRet) << std::endl;
        std::cout << "col = " << mxGetN(mRet) << std::endl;
    }

    // save & print coefficients
    double coef_smooth[F];
    double coef_derivative[F];

    for (size_t i = 0; i < F; i++) {
        coef_smooth[i] = ret[i];
        coef_derivative[i] = ret[i+F];
        std::cout << "s = " << coef_smooth[i] << " d = " << coef_derivative[i] << std::endl;
    }

    // close engine & clean up
    std::cout << "Done" << std::endl;
    mxDestroyArray(mN);
    mxDestroyArray(mF);
    mxDestroyArray(mRet);
    engClose(meng);

    return EXIT_SUCCESS;
}