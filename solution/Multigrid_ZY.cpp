//-----------Written by ZhangYu HIT-------------------
#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <ctime>

#include <Eigen/Dense>
#include <Eigen/Sparse>

using namespace Eigen;
using namespace std;

int get_NodeCondition(const int ii, const int ij, const int N2);

inline int coordinate(const int x, const int y, const int NX)
{
    return (x + y * NX);
}

template <typename T>
inline T p2(const T x)
{
    return (x * x);
}

namespace Basic
{
    const int N = 80 + 1;
    const double MAX_ERR = 1e-6;

    // ---------------------------------------
    const bool uni_mesh = 1;

    double S_x = 1.0 / 0.98;
    double S_y = 0.98;

    //--------------------Node Spaces-----------------
    vector<double> X_delta(N - 1, 0);
    vector<double> Y_delta(N - 1, 0);

    //--------------------Node Positions--------------
    vector<double> X_position(N, 0);
    vector<double> Y_position(N, 0);

}

namespace Solve_equ
{
    using Basic::N;

    double L2_error;

    // ----------------- FD Coeffieient Matrix
    //MatrixXd A = MatrixXd::Constant(N * N, N * N, 0);

    SparseMatrix<double, RowMajor> A(N *N, N *N);

    SparseMatrix<double, RowMajor> A1(p2(N - 2), p2(N - 2));

    // -----------------[A][X] = [B]
    VectorXd B = VectorXd::Constant(N * N, 0);

    VectorXd B1 = VectorXd::Constant(p2(N - 2), 0);

    VectorXd X = VectorXd::Constant(N * N, 1);

    VectorXd X1 = VectorXd::Constant(p2(N - 2), 1);

    VectorXd phi_ana = VectorXd::Constant(N * N, 0);
}

namespace MultiGrids
{

}

void get_Bacis_Mesh(const bool uni_mesh, const int N, double S_x, double S_y,
                    vector<double> &X_delta, vector<double> &Y_delta,
                    vector<double> &X_position, vector<double> &Y_position)
{
    // using namespace Basic;

    if (uni_mesh != 1)
    {
        S_x = 1.0 / S_x;
        S_y = 1.0 / S_y;
    }
    else
    {
        S_x = S_y = 1.0;
    }

    double x0, y0;

    if (uni_mesh == 1)
    {
        x0 = 1.0 / (N - 1);
        y0 = 1.0 / (N - 1);
    }
    else
    {
        x0 = 1.0 * (1 - S_x) / (1 - pow(S_x, N - 1));
        y0 = 1.0 * (1 - S_y) / (1 - pow(S_y, N - 1));
    }

    X_delta[0] = x0;
    Y_delta[0] = y0;

    cout << "-----------N = " << N << "   deltax   deltay----------" << endl;
    for (int i = 1; i < N - 1; ++i)
    {
        X_delta[i] = X_delta[i - 1] * S_x;
        Y_delta[i] = Y_delta[i - 1] * S_y;
        cout << X_delta[i] << " " << Y_delta[i] << endl;
    }

    cout << "-----------N = " << N << "    positionx   positiony----------" << endl;
    for (int i = 1; i < N; ++i)
    {
        X_position[i] = X_position[i - 1] + X_delta[i - 1];
        Y_position[i] = Y_position[i - 1] + Y_delta[i - 1];
        cout << X_position[i] << " " << Y_position[i] << endl;
    }
}

inline double anaSolu(const double X_, const double Y_)
{
    //double ana = (5000000*p2(Y_) + 5000000*p2(X_ - 1) - 100000)*exp(-50*p2(Y_) - 50*p2(1 - X_));
    double ana = 500 * exp(-50 * (p2(1.0 - X_) + p2(Y_))) + 100 * X_ * (1 - Y_);
    return ana;
}

void get_Analytic_Solution()
{

    using Basic::N;
    using Basic::X_position;
    using Basic::Y_position;

    using Solve_equ::phi_ana;

    for (int j = 0; j < N; ++j)
    {
        for (int i = 0; i < N; ++i)
        {
            const int index = coordinate(i, j, N);
            const double X_ = X_position[i];
            const double Y_ = Y_position[j];
            //phi_ana(index) = 500 * exp(-50 * (p2(1.0 - X_) + p2(Y_))) + 100 * X_ * (1 - Y_);
            phi_ana(index) = anaSolu(X_, Y_);
        }
    }
}

int get_NodeCondition(const int ii, const int ij, const int N2)
{

    int NodeCondition;
    // ------------6 2 5
    // ------------3 0 1
    // ------------7 4 8
    if (ii != 0 && ij != 0 && ii != N2 - 1 && ij != N2 - 1)
    {
        NodeCondition = 0;
    }
    else if (ij == 0 && ii == 0)
    {
        //NodeCondition = "LeftBottom";
        NodeCondition = 7;
    }
    else if (ij == 0 && ii != 0 && ii != N2 - 1)
    {
        //NodeCondition = "Bottom";
        NodeCondition = 4;
    }
    else if (ij == 0 && ii == N2 - 1)
    {
        //NodeCondition = "RightBottom";
        NodeCondition = 8;
    }
    else if (ii == 0 && ij != 0 && ij != N2 - 1)
    {
        //NodeCondition = "Left";
        NodeCondition = 3;
    }
    else if (ii == N2 - 1 && ij != 0 && ij != N2 - 1)
    {
        //NodeCondition = "Right";
        NodeCondition = 1;
    }
    else if (ii == N2 - 1 && ij == N2 - 1)
    {
        //NodeCondition = "RightTop";
        NodeCondition = 5;
    }
    else if (ij == N2 - 1 && ii != 0 && ii != N2 - 1)
    {
        //NodeCondition = "Top";
        NodeCondition = 2;
    }
    else if (ii == 0 && ij == N2 - 1)
    {
        //NodeCondition = "LeftTop";
        NodeCondition = 6;
    }
    return NodeCondition;
}

inline double SRight(const double X_, const double Y_)
{
    double S = (5000000 * p2(Y_) + 5000000 * p2(X_ - 1) - 100000) * exp(-50 * p2(Y_) - 50 * p2(1 - X_));
    return S;
}

void get_Right_term_resi(const int N, const double *X_delta, const double *Y_delta,
                         VectorXd &r, VectorXd &B)
{
    const int N2 = N - 2;

    for (int i = 0; i < p2(N2); ++i)
    {

        int ii = i % ((int)N2);
        int ij = i / ((int)N2);

        int NodeCondition = get_NodeCondition(ii, ij, N2);

        double dx1 = X_delta[1 + ii - 1];
        double dx2 = X_delta[1 + ii];
        double dy1 = Y_delta[1 + ij - 1];
        double dy2 = Y_delta[1 + ij];

        double B_N, B_S, B_W, B_E;

        switch (NodeCondition)
        {
        case 1:
            B_E = 0;
            B(i) = r(i) - 2.0 * B_E / (dx2 * (dx1 + dx2)); //---E
            break;
        case 2:
            B_N = 0;
            B(i) = r(i) - 2.0 * B_N / (dy2 * (dy1 + dy2)); //---N
            break;
        case 3:
            B_W = 0;
            B(i) = r(i) - 2.0 * B_W / (dx1 * (dx1 + dx2)); //---W
            break;
        case 4:
            B_S = 0;
            B(i) = r(i) - 2.0 * B_S / (dy1 * (dy1 + dy2)); //---S
            break;
        case 5:
            B_N = 0;
            B_E = 0;
            B(i) = r(i) - (2.0 * B_E / (dx2 * (dx1 + dx2)) + 2.0 * B_N / (dy2 * (dy1 + dy2))); //---NE
            break;
        case 6:
            B_N = 0;
            B_W = 0;
            B(i) = r(i) - (2.0 * B_W / (dx1 * (dx1 + dx2)) + 2.0 * B_N / (dy2 * (dy1 + dy2))); //---NW
            break;
        case 7:
            B_S = 0;
            B_W = 0;
            B(i) = r(i) - (2.0 * B_W / (dx1 * (dx1 + dx2)) + 2.0 * B_S / (dy1 * (dy1 + dy2))); //---SW
            break;
        case 8:
            B_S = 0;
            B_E = 0;
            B(i) = r(i) - (2.0 * B_E / (dx2 * (dx1 + dx2)) + 2.0 * B_S / (dy1 * (dy1 + dy2))); //---SE
            break;
        default:
            B(i) = r(i);
            break;
        }
    }
}

void get_Right_term_interior(const int N, const double *X_delta, const double *Y_delta,
                             vector<double> &X_position, vector<double> &Y_position,
                             VectorXd &B1)
{
    const int N2 = N - 2;

    for (int i = 0; i < p2(N2); ++i)
    {

        int ii = i % ((int)N2);
        int ij = i / ((int)N2);

        int NodeCondition = get_NodeCondition(ii, ij, N2);

        double dx1 = X_delta[1 + ii - 1];
        double dx2 = X_delta[1 + ii];
        double dy1 = Y_delta[1 + ij - 1];
        double dy2 = Y_delta[1 + ij];

        double *X_ = &X_position[0] + ii + 1;
        double *Y_ = &Y_position[0] + ij + 1;

        double B_N, B_S, B_W, B_E;

        B1(i) = SRight(*X_, *Y_);
        switch (NodeCondition)
        {
        case 1:
            B_E = anaSolu(*(X_ + 1), *Y_);
            B1(i) -= 2.0 * B_E / (dx2 * (dx1 + dx2)); //---E
            break;
        case 2:
            B_N = anaSolu(*X_, *(Y_ + 1));
            B1(i) -= 2.0 * B_N / (dy2 * (dy1 + dy2)); //---N
            break;
        case 3:
            B_W = anaSolu(*(X_ - 1), *Y_);
            B1(i) -= 2.0 * B_W / (dx1 * (dx1 + dx2)); //---W
            break;
        case 4:
            B_S = anaSolu(*X_, *(Y_ - 1));
            B1(i) -= 2.0 * B_S / (dy1 * (dy1 + dy2)); //---S
            break;
        case 5:
            B_N = anaSolu(*X_, *(Y_ + 1));
            B_E = anaSolu(*(X_ + 1), *Y_);
            B1(i) -= (2.0 * B_E / (dx2 * (dx1 + dx2)) + 2.0 * B_N / (dy2 * (dy1 + dy2))); //---NE
            break;
        case 6:
            B_N = anaSolu(*X_, *(Y_ + 1));
            B_W = anaSolu(*(X_ - 1), *Y_);
            B1(i) -= (2.0 * B_W / (dx1 * (dx1 + dx2)) + 2.0 * B_N / (dy2 * (dy1 + dy2))); //---NW
            break;
        case 7:
            B_S = anaSolu(*X_, *(Y_ - 1));
            B_W = anaSolu(*(X_ - 1), *Y_);
            B1(i) -= (2.0 * B_W / (dx1 * (dx1 + dx2)) + 2.0 * B_S / (dy1 * (dy1 + dy2))); //---SW
            break;
        case 8:
            B_S = anaSolu(*X_, *(Y_ - 1));
            B_E = anaSolu(*(X_ + 1), *Y_);
            B1(i) -= (2.0 * B_E / (dx2 * (dx1 + dx2)) + 2.0 * B_S / (dy1 * (dy1 + dy2))); //---SE
            break;
        default:
            //B1(i) = SRight(*X_, *Y_);
            break;
        }
    }
    //cout << B1 << endl;
}

void get_iter_coeff(VectorXd &A_c, const double *X_delta, const double *Y_delta, const int N)
{
    int N2 = N - 2;
    int N22 = p2(N2);

    for (int i = 0; i < N22; ++i)
    {

        int ii = i % ((int)N2);
        int ij = i / ((int)N2);

        int im5 = i * 5;

        double dx1 = X_delta[1 + ii - 1];
        double dx2 = X_delta[1 + ii];
        double dy1 = Y_delta[1 + ij - 1];
        double dy2 = Y_delta[1 + ij];

        int NodeCondition = get_NodeCondition(ii, ij, N2);
        A_c(im5) = -(2.0 / (dx1 * dx2) + 2.0 / (dy1 * dy2));

        switch (NodeCondition)
        {
            // center - left - bottom - right - top
        case 0:
            A_c(im5 + 1) = 2.0 / (dx1 * (dx1 + dx2)); //---W
            A_c(im5 + 2) = 2.0 / (dy1 * (dy1 + dy2)); //---S
            A_c(im5 + 3) = 2.0 / (dx2 * (dx1 + dx2)); //---E
            A_c(im5 + 4) = 2.0 / (dy2 * (dy1 + dy2)); //---N
            break;
        case 1:
            A_c(im5 + 1) = 2.0 / (dx1 * (dx1 + dx2)); //---W
            A_c(im5 + 2) = 2.0 / (dy1 * (dy1 + dy2)); //---S
            A_c(im5 + 4) = 2.0 / (dy2 * (dy1 + dy2)); //---N
            break;
        case 2:
            A_c(im5 + 1) = 2.0 / (dx1 * (dx1 + dx2)); //---W
            A_c(im5 + 2) = 2.0 / (dy1 * (dy1 + dy2)); //---S
            A_c(im5 + 3) = 2.0 / (dx2 * (dx1 + dx2)); //---E
            break;
        case 3:
            A_c(im5 + 2) = 2.0 / (dy1 * (dy1 + dy2)); //---S
            A_c(im5 + 3) = 2.0 / (dx2 * (dx1 + dx2)); //---E
            A_c(im5 + 4) = 2.0 / (dy2 * (dy1 + dy2)); //---N
            break;
        case 4:
            A_c(im5 + 1) = 2.0 / (dx1 * (dx1 + dx2)); //---W
            A_c(im5 + 3) = 2.0 / (dx2 * (dx1 + dx2)); //---E
            A_c(im5 + 4) = 2.0 / (dy2 * (dy1 + dy2)); //---N
            break;
        case 5:
            A_c(im5 + 1) = 2.0 / (dx1 * (dx1 + dx2)); //---W
            A_c(im5 + 2) = 2.0 / (dy1 * (dy1 + dy2)); //---S
            break;
        case 6:
            A_c(im5 + 2) = 2.0 / (dy1 * (dy1 + dy2)); //---S
            A_c(im5 + 3) = 2.0 / (dx2 * (dx1 + dx2)); //---E
            break;
        case 7:
            A_c(im5 + 3) = 2.0 / (dx2 * (dx1 + dx2)); //---E
            A_c(im5 + 4) = 2.0 / (dy2 * (dy1 + dy2)); //---N
            break;
        case 8:
            A_c(im5 + 1) = 2.0 / (dx1 * (dx1 + dx2)); //---W
            A_c(im5 + 4) = 2.0 / (dy2 * (dy1 + dy2)); //---N
            break;

        default:
            break;
        }
    }
}

void Gauss_Seidel(VectorXd &X, VectorXd &Ac, VectorXd &B, const int N)
{
    for (int j = 0; j < N; ++j)
    {
        for (int i = 0; i < N; ++i)
        {
            const int index = coordinate(i, j, N);

            const int i_L = coordinate(i - 1, j, N);
            const int i_B = coordinate(i, j - 1, N);
            const int i_R = coordinate(i + 1, j, N);
            const int i_T = coordinate(i, j + 1, N);

            const int im5 = 5 * index;

            int NodeCondition = get_NodeCondition(i, j, N);

            switch (NodeCondition)
            {
            case 0:
                X(index) = 1.0 / Ac(im5) * (B(index) - (Ac(im5 + 1) * X(i_L) + Ac(im5 + 2) * X(i_B) + Ac(im5 + 3) * X(i_R) + Ac(im5 + 4) * X(i_T)));
                break;
            case 1:
                X(index) = 1.0 / Ac(im5) * (B(index) - (Ac(im5 + 1) * X(i_L) + Ac(im5 + 2) * X(i_B) + Ac(im5 + 4) * X(i_T)));
                break;
            case 2:
                X(index) = 1.0 / Ac(im5) * (B(index) - (Ac(im5 + 1) * X(i_L) + Ac(im5 + 2) * X(i_B) + Ac(im5 + 3) * X(i_R)));
                break;
            case 3:
                X(index) = 1.0 / Ac(im5) * (B(index) - (Ac(im5 + 2) * X(i_B) + Ac(im5 + 3) * X(i_R) + Ac(im5 + 4) * X(i_T)));
                break;
            case 4:
                X(index) = 1.0 / Ac(im5) * (B(index) - (Ac(im5 + 1) * X(i_L) + Ac(im5 + 3) * X(i_R) + Ac(im5 + 4) * X(i_T)));
                break;
            case 5:
                X(index) = 1.0 / Ac(im5) * (B(index) - (Ac(im5 + 1) * X(i_L) + Ac(im5 + 2) * X(i_B)));
                break;
            case 6:
                X(index) = 1.0 / Ac(im5) * (B(index) - (Ac(im5 + 2) * X(i_B) + Ac(im5 + 3) * X(i_R)));
                break;
            case 7:
                X(index) = 1.0 / Ac(im5) * (B(index) - (Ac(im5 + 3) * X(i_R) + Ac(im5 + 4) * X(i_T)));
                break;
            case 8:
                X(index) = 1.0 / Ac(im5) * (B(index) - (Ac(im5 + 1) * X(i_L) + Ac(im5 + 4) * X(i_T)));
                break;

            default:
                break;
            }
        }
    }
}

void from_X1_to_X(VectorXd &X, VectorXd &X1, int N)
{
    using Solve_equ::phi_ana;

    for (int j = 0; j < N; ++j)
    {
        for (int i = 0; i < N; ++i)
        {

            const int index = coordinate(i, j, N);

            int index_in;

            int NodeCondition = get_NodeCondition(i, j, N);

            switch (NodeCondition)
            {
            case 0:
                index_in = coordinate(i - 1, j - 1, N - 2);
                X(index) = X1(index_in);
                break;
            //case 1:
            //    X(index) = phi_ana(index);
            //    break;
            //case 2:
            //    X(index) = phi_ana(index);
            //    break;
            //case 3:
            //    X(index) = phi_ana(index);
            //    break;
            //case 4:
            //    X(index) = phi_ana(index);
            //    break;
            //case 5:
            //    X(index) = phi_ana(index);
            //    break;
            //case 6:
            //    X(index) = phi_ana(index);
            //    break;
            //case 7:
            //    X(index) = phi_ana(index);
            //    break;
            //case 8:
            //    X(index) = phi_ana(index);
            //    break;
            default:
                X(index) = phi_ana(index);
                break;
            }
        }
    }
}

void Gauss_Seidel_Iteration(VectorXd &X, VectorXd &X1, VectorXd &B,
                            vector<double> &X_delta, vector<double> &Y_delta,
                            vector<double> &X_position, vector<double> &Y_position, const int N)
{
    // ********** The Gauss Seidel Iteration don't include Boundary nodes
    const int N2 = N - 2;
    // FD coefficient----------------
    VectorXd A_coeff = VectorXd::Constant(p2(N2) * 5, 0);

    get_iter_coeff(A_coeff, &X_delta[0], &Y_delta[0], N);

    get_Right_term_interior(N, &X_delta[0], &Y_delta[0], X_position, Y_position, B);

    //Gauss_Seidel(X1, A_coeff, B, N);

    double err_sum = 1.0;
    int iter_num = 0;

    while (err_sum > 1e-8)
    {
        VectorXd X_old = X1;
        Gauss_Seidel(X1, A_coeff, B, N2);

        VectorXd error = X1 - X_old;
        err_sum = error.lpNorm<2>();

        ++iter_num;
    }

    cout << "Iter steps is " << iter_num << ", "
         << "iter max err is " << err_sum << endl;

    from_X1_to_X(X, X1, N);
}

void residual(VectorXd &rh, const int N,
              const MatrixXd &Ah_coeff, const VectorXd &Xh, const VectorXd &Bh)
{
    for (int j = 0; j < N; ++j)
    {
        for (int i = 0; i < N; ++i)
        {
            const int index = coordinate(i, j, N);

            const int i_L = coordinate(i - 1, j, N);
            const int i_B = coordinate(i, j - 1, N);
            const int i_R = coordinate(i + 1, j, N);
            const int i_T = coordinate(i, j + 1, N);

            const int im5 = 5 * index;

            int NodeCondition = get_NodeCondition(i, j, N);

            switch (NodeCondition)
            {
            case 0:
                rh(index) = Bh(index) - (Ah_coeff(im5) * Xh(index) + Ah_coeff(im5 + 1) * Xh(i_L) + Ah_coeff(im5 + 2) * Xh(i_B) + Ah_coeff(im5 + 3) * Xh(i_R) + Ah_coeff(im5 + 4) * Xh(i_T));
                break;
            case 1:
                rh(index) = Bh(index) - (Ah_coeff(im5) * Xh(index) + Ah_coeff(im5 + 1) * Xh(i_L) + Ah_coeff(im5 + 2) * Xh(i_B) + Ah_coeff(im5 + 4) * Xh(i_T));
                break;
            case 2:
                rh(index) = Bh(index) - (Ah_coeff(im5) * Xh(index) + Ah_coeff(im5 + 1) * Xh(i_L) + Ah_coeff(im5 + 2) * Xh(i_B) + Ah_coeff(im5 + 3) * Xh(i_R));
                break;
            case 3:
                rh(index) = Bh(index) - (Ah_coeff(im5) * Xh(index) + Ah_coeff(im5 + 2) * Xh(i_B) + Ah_coeff(im5 + 3) * Xh(i_R) + Ah_coeff(im5 + 4) * Xh(i_T));
                break;
            case 4:
                rh(index) = Bh(index) - (Ah_coeff(im5) * Xh(index) + Ah_coeff(im5 + 1) * Xh(i_L) + Ah_coeff(im5 + 3) * Xh(i_R) + Ah_coeff(im5 + 4) * Xh(i_T));
                break;
            case 5:
                rh(index) = Bh(index) - (Ah_coeff(im5) * Xh(index) + Ah_coeff(im5 + 1) * Xh(i_L) + Ah_coeff(im5 + 2) * Xh(i_B));
                break;
            case 6:
                rh(index) = Bh(index) - (Ah_coeff(im5) * Xh(index) + Ah_coeff(im5 + 2) * Xh(i_B) + Ah_coeff(im5 + 3) * Xh(i_R));
                break;
            case 7:
                rh(index) = Bh(index) - (Ah_coeff(im5) * Xh(index) + Ah_coeff(im5 + 3) * Xh(i_R) + Ah_coeff(im5 + 4) * Xh(i_T));
                break;
            case 8:
                rh(index) = Bh(index) - (Ah_coeff(im5) * Xh(index) + Ah_coeff(im5 + 1) * Xh(i_L) + Ah_coeff(im5 + 4) * Xh(i_T));
                break;

            default:
                break;
            }
        }
    }
}

void Restriction(VectorXd &r1, const VectorXd &r0, const int N1, const int N0)
{

    // ------------------- simply injection
    for (int i = 0; i < N1; ++i)
    {
        for (int j = 0; j < N1; ++j)
        {
            const int index1 = coordinate(i, j, N1);
            const int index0 = coordinate(i * 2 + 1, j * 2 + 1, N0);
            r1(index1) = r0(index0);
        }
    }
}

void Prolongation(const VectorXd &eh1, VectorXd &rh0, const int N1, const int N0,
                  vector<double> &X_weight, vector<double> &Y_weight)
{
    // Bilinear interpolation -----------------------
    for (int i = 0; i < N1 - 1; ++i)
    {
        for (int j = 0; j < N1 - 1; ++j)
        {
            const int ip1 = i + 1;
            const int jp1 = j + 1;
            const int im2p1 = i * 2 + 1;
            const int im2p2 = i * 2 + 2;
            const int jm2p1 = j * 2 + 1;
            const int jm2p2 = j * 2 + 2;

            const int index1 = coordinate(i, j, N1);
            const int index1_i = coordinate(ip1, j, N1);
            const int index1_j = coordinate(i, jp1, N1);
            const int index1_ij = coordinate(ip1, jp1, N1);

            const int index0 = coordinate(im2p1, jm2p1, N0);
            const int index0_i = coordinate(im2p2, jm2p1, N0);
            const int index0_j = coordinate(im2p1, jm2p2, N0);
            const int index0_ij = coordinate(im2p2, jm2p2, N0);

            rh0(index0) = eh1(index1);
            rh0(index0_i) = X_weight[im2p2] * eh1(index1) + X_weight[im2p2 + 1] * eh1(index1_i);
            rh0(index0_j) = Y_weight[jm2p2] * eh1(index1) + Y_weight[jm2p2 + 1] * eh1(index1_j);
            rh0(index0_ij) = Y_weight[jm2p2] * (X_weight[im2p2] * eh1(index1) + X_weight[im2p2 + 1] * eh1(index1_i)) + Y_weight[jm2p2 + 1] * (X_weight[im2p2] * eh1(index1_j) + X_weight[im2p2 + 1] * eh1(index1_ij));
        }
    }

    const int N1_1 = N1 - 1;

    for (int i = 0; i < N1_1; ++i)
    {
        const int ip1 = i + 1;
        const int im2p1 = i * 2 + 1;
        const int im2p2 = i * 2 + 2;

        const int index1 = coordinate(i, N1_1, N1);
        const int index1_1 = coordinate(ip1, N1_1, N1);

        const int index0 = coordinate(im2p1, N1_1 * 2 + 1, N0);
        const int index0_1 = coordinate(im2p2, N1_1 * 2 + 1, N0);

        rh0(index0) = eh1(index1);
        rh0(index0_1) = X_weight[im2p2] * eh1(index1) + X_weight[im2p2 + 1] * eh1(index1_1);

        const int index0_pj = coordinate(im2p1, N1_1 * 2 + 2, N0);
        const int index0_1pj = coordinate(im2p2, N1_1 * 2 + 2, N0);
        rh0(index0_pj) = Y_weight[N1_1 * 2 + 2] * rh0(index0) + Y_weight[N1_1 * 2 + 3] * 0;
        rh0(index0_1pj) = Y_weight[N1_1 * 2 + 2] * rh0(index0_1) + Y_weight[N1_1 * 2 + 3] * 0;

        const int index0_0mj = coordinate(im2p1, 1, N0);
        const int index0_01mj = coordinate(im2p2, 1, N0);
        rh0(coordinate(im2p1, 0, N0)) = Y_weight[1] * rh0(coordinate(im2p1, 1, N0)) + Y_weight[0] * 0;
        rh0(coordinate(im2p2, 0, N0)) = Y_weight[1] * rh0(coordinate(im2p2, 1, N0)) + Y_weight[0] * 0;
    }

    for (int j = 0; j < N1_1; ++j)
    {
        const int jp1 = j + 1;
        const int jm2p1 = j * 2 + 1;
        const int jm2p2 = j * 2 + 2;

        const int index1 = coordinate(N1_1, j, N1);
        const int index1_1 = coordinate(N1_1, jp1, N1);

        const int index0 = coordinate(N1_1 * 2 + 1, jm2p1, N0);
        const int index0_1 = coordinate(N1_1 * 2 + 1, jm2p2, N0);

        rh0(index0) = eh1(index1);
        rh0(index0_1) = Y_weight[jm2p2] * eh1(index1) + Y_weight[jm2p2 + 1] * eh1(index1_1);

        const int index0_pi = coordinate(N1_1 * 2 + 2, jm2p1, N0);
        const int index0_1pi = coordinate(N1_1 * 2 + 2, jm2p2, N0);
        rh0(index0_pi) = X_weight[N1_1 * 2 + 2] * rh0(index0) + X_weight[N1_1 * 2 + 3] * 0;
        rh0(index0_1pi) = X_weight[N1_1 * 2 + 2] * rh0(index0_1) + X_weight[N1_1 * 2 + 3] * 0;

        const int index0_0mi = coordinate(1, jm2p1, N0);
        const int index0_01mi = coordinate(1, jm2p2, N0);
        rh0(coordinate(0, jm2p1, N0)) = X_weight[1] * rh0(coordinate(1, jm2p1, N0)) + X_weight[0] * 0;
        rh0(coordinate(0, jm2p2, N0)) = X_weight[1] * rh0(coordinate(1, jm2p2, N0)) + X_weight[0] * 0;
    }

    rh0(coordinate(2 * N1_1 + 1, 2 * N1_1 + 1, N0)) = eh1(coordinate(N1_1, N1_1, N1));

    rh0(coordinate(2 * N1_1 + 1, 0, N0)) = Y_weight[1] * rh0(coordinate(2 * N1_1 + 1, 1, N0)) + Y_weight[0] * 0;
    rh0(coordinate(2 * N1_1 + 1, 2 * N1_1 + 2, N0)) = Y_weight[N1_1 * 2 + 2] * rh0(coordinate(2 * N1_1 + 1, 2 * N1_1 + 1, N0)) + Y_weight[N1_1 * 2 + 3] * 0;

    rh0(coordinate(0, 2 * N1_1 + 1, N0)) = X_weight[1] * rh0(coordinate(1, 2 * N1_1 + 1, N0)) + Y_weight[0] * 0;
    rh0(coordinate(2 * N1_1 + 2, 2 * N1_1 + 1, N0)) = X_weight[N1_1 * 2 + 2] * rh0(coordinate(2 * N1_1 + 1, 2 * N1_1 + 1, N0)) + Y_weight[N1_1 * 2 + 3] * 0;

    rh0(coordinate(0, 0, N0)) = X_weight[1] * rh0(coordinate(1, 0, N0));
    rh0(coordinate(0, N0 - 1, N0)) = X_weight[1] * rh0(coordinate(1, N0 - 1, N0));
    rh0(coordinate(N0 - 1, 0, N0)) = X_weight[N0 - 1] * rh0(coordinate(N0 - 2, 0, N0));
    rh0(coordinate(N0 - 1, N0 - 1, N0)) = X_weight[N0 - 1] * rh0(coordinate(N0 - 2, N0 - 1, N0));
}

void get_Weight_Interpolation(vector<double> &weight, const int N_c,
                              const vector<double> &X_f, const vector<double> &X_c)
{
    for (int i = 0; i < N_c; ++i)
    {
        const int i2 = 2 * i;
        //const double inv_delta = 1.0 / X_c[i];
        const double inv_delta = 1.0 / (X_f[i2 + 1] + X_f[i2]);

        weight[i2] = X_f[i2 + 1] * inv_delta;
        weight[i2 + 1] = X_f[i2 + 0] * inv_delta;
    }
}
void get_resi_mesh(const int N1, vector<double> &X1_delta, vector<double> &Y1_delta,
                   vector<double> &X1_position, vector<double> &Y1_position,
                   vector<double> &X0_delta, vector<double> &Y0_delta,
                   vector<double> &X0_position, vector<double> &Y0_position)
{
    for(int i = 0; i < N1 - 1; ++i){
        X1_delta[i] = X0_delta[2 * i] + X0_delta[2 * i + 1];
        Y1_delta[i] = Y0_delta[2 * i] + Y0_delta[2 * i + 1];
    }

    X1_position[0] = X0_position[0];
    Y1_position[0] = Y0_position[0];
    for (int i = 1; i < N1; ++i)
    {
        X1_position[i] = X1_position[i - 1] + X1_delta[i - 1];
        Y1_position[i] = Y1_position[i - 1] + Y1_delta[i - 1];
        // cout << X_position[i] << " " << Y_position[i] << endl;
    }
}

void MultiGrid_Iter(const int N, VectorXd &X, VectorXd &Xh0, VectorXd &Bh0,
                    const bool unimesh, const double S_x, const double S_y,
                    vector<double> &Xh0_delta, vector<double> &Yh0_delta,
                    vector<double> &Xh0_position, vector<double> &Yh0_position)
{
    const int N0 = N - 2;
    const int N1 = (N0 - 1) / 2;
    const int N2 = (N1 - 1) / 2;
    const int N3 = (N2 - 1) / 2;

    //--------------------Node Spaces-----------------
    vector<double> Xh1_delta(N1 + 2 - 1, 0);
    vector<double> Yh1_delta(N1 + 2 - 1, 0);
    vector<double> Xh2_delta(N2 + 2 - 1, 0);
    vector<double> Yh2_delta(N2 + 2 - 1, 0);
    vector<double> Xh3_delta(N3 + 2 - 1, 0);
    vector<double> Yh3_delta(N3 + 2 - 1, 0);

    vector<double> X01_weight(N0 + 2 - 1, 0);
    vector<double> Y01_weight(N0 + 2 - 1, 0);
    vector<double> X12_weight(N1 + 2 - 1, 0);
    vector<double> Y12_weight(N1 + 2 - 1, 0);
    vector<double> X23_weight(N2 + 2 - 1, 0);
    vector<double> Y23_weight(N2 + 2 - 1, 0);

    //--------------------Node Positions--------------
    vector<double> Xh1_position(N1 + 2, 0);
    vector<double> Yh1_position(N1 + 2, 0);
    vector<double> Xh2_position(N2 + 2, 0);
    vector<double> Yh2_position(N2 + 2, 0);
    vector<double> Xh3_position(N3 + 2, 0);
    vector<double> Yh3_position(N3 + 2, 0);

    get_resi_mesh(N1 + 2, Xh1_delta, Yh1_delta, Xh1_position, Yh1_position, Xh0_delta, Yh0_delta, Xh0_position, Yh0_position);
    get_resi_mesh(N2 + 2, Xh2_delta, Yh2_delta, Xh2_position, Yh2_position, Xh1_delta, Yh1_delta, Xh1_position, Yh1_position);
    get_resi_mesh(N3 + 2, Xh3_delta, Yh3_delta, Xh3_position, Yh3_position, Xh2_delta, Yh2_delta, Xh2_position, Yh2_position);
    // get_Bacis_Mesh(unimesh, N1 + 2, S_x, S_y, Xh1_delta, Yh1_delta, Xh1_position, Yh1_position);
    // get_Bacis_Mesh(unimesh, N2 + 2, S_x, S_y, Xh2_delta, Yh2_delta, Xh2_position, Yh2_position);
    // get_Bacis_Mesh(unimesh, N3 + 2, S_x, S_y, Xh3_delta, Yh3_delta, Xh3_position, Yh3_position);

    get_Weight_Interpolation(X01_weight, N1 + 1, Xh0_delta, Xh1_delta);
    get_Weight_Interpolation(Y01_weight, N1 + 1, Yh0_delta, Yh1_delta);
    get_Weight_Interpolation(X12_weight, N2 + 1, Xh1_delta, Xh2_delta);
    get_Weight_Interpolation(Y12_weight, N2 + 1, Yh1_delta, Yh2_delta);
    get_Weight_Interpolation(X23_weight, N3 + 1, Xh2_delta, Xh3_delta);
    get_Weight_Interpolation(Y23_weight, N3 + 1, Yh2_delta, Yh3_delta);

    VectorXd rh0 = VectorXd::Constant(p2(N0), 0);
    VectorXd eh0 = VectorXd::Constant(p2(N0), 0);
    VectorXd rh1 = VectorXd::Constant(p2(N1), 0);
    VectorXd eh1 = VectorXd::Constant(p2(N1), 0);
    VectorXd rh2 = VectorXd::Constant(p2(N2), 0);
    VectorXd eh2 = VectorXd::Constant(p2(N2), 0);
    VectorXd rh3 = VectorXd::Constant(p2(N3), 0);
    VectorXd eh3 = VectorXd::Constant(p2(N3), 0);

    VectorXd Bh1 = VectorXd::Constant(p2(N1), 1);
    VectorXd Bh2 = VectorXd::Constant(p2(N2), 1);
    VectorXd Bh3 = VectorXd::Constant(p2(N3), 1);

    VectorXd Ah0_coeff = VectorXd::Constant(p2(N0) * 5, 0);
    VectorXd Ah1_coeff = VectorXd::Constant(p2(N1) * 5, 0);
    VectorXd Ah2_coeff = VectorXd::Constant(p2(N2) * 5, 0);
    VectorXd Ah3_coeff = VectorXd::Constant(p2(N3) * 5, 0);

    get_iter_coeff(Ah0_coeff, &Xh0_delta[0], &Yh0_delta[0], N0 + 2);
    get_Right_term_interior(N0 + 2, &Xh0_delta[0], &Yh0_delta[0], Xh0_position, Yh0_position, Bh0);

    get_iter_coeff(Ah1_coeff, &Xh1_delta[0], &Yh1_delta[0], N1 + 2);
    get_iter_coeff(Ah2_coeff, &Xh2_delta[0], &Yh2_delta[0], N2 + 2);
    get_iter_coeff(Ah3_coeff, &Xh3_delta[0], &Yh3_delta[0], N3 + 2);

    double err_sum = 1.0;
    int iter_num = 0;
    VectorXd X_old = Xh0;
    VectorXd error = Xh0;

    std::ofstream err_out("err.txt");
    for (int step = 0; step < 20000; ++step)
    //for (int step = 0; step < 1; ++step)
    {
        X_old = Xh0;
        eh1.setZero();
        eh2.setZero();
        eh3.setZero();

        for (int steph0 = 0; steph0 < 5; ++steph0)
        {
            Gauss_Seidel(Xh0, Ah0_coeff, Bh0, N0);
            error = Xh0 - X_old;
            err_sum = error.lpNorm<2>();
        }

        if (err_sum < 1e-6)
        {
            cout << "Iterations of MG = " << step << " x 5" << endl;
            break;
        }

        err_out << step << " " << err_sum << endl;

        residual(rh0, N0, Ah0_coeff, Xh0, Bh0);

        Restriction(rh1, rh0, N1, N0);
        get_Right_term_resi(N1 + 2, &Xh1_delta[0], &Yh1_delta[0], rh1, Bh1);

        for (int steph0 = 0; steph0 < 10; ++steph0)
        {
            Gauss_Seidel(eh1, Ah1_coeff, rh1, N1);
        }

        residual(rh1, N1, Ah1_coeff, eh1, Bh1);
        Restriction(rh2, rh1, N2, N1);
        get_Right_term_resi(N2 + 2, &Xh2_delta[0], &Yh2_delta[0], rh2, Bh2);
        for (int steph0 = 0; steph0 < 5; ++steph0)
        {
            Gauss_Seidel(eh2, Ah2_coeff, Bh2, N2);
        }

        residual(rh2, N2, Ah2_coeff, eh2, Bh2);
        Restriction(rh3, rh2, N3, N2);
        get_Right_term_resi(N3 + 2, &Xh3_delta[0], &Yh3_delta[0], rh3, Bh3);
        for (int steph0 = 0; steph0 < 20; ++steph0)
        {
            Gauss_Seidel(eh3, Ah3_coeff, Bh3, N3);
        }
        //   cout << eh3 << endl;

        //   -------------------------------------------------------
        Prolongation(eh3, rh2, N3, N2, X23_weight, Y23_weight);
        eh2 += rh2;
        for (int steph0 = 0; steph0 < 5; ++steph0)
        {
            Gauss_Seidel(eh2, Ah2_coeff, Bh2, N2);
        }

        Prolongation(eh2, rh1, N2, N1, X12_weight, Y12_weight);
        eh1 += rh1;
        for (int steph0 = 0; steph0 < 5; ++steph0)
        {
            Gauss_Seidel(eh1, Ah1_coeff, Bh1, N1);
        }

        Prolongation(eh1, rh0, N1, N0, X01_weight, Y01_weight);
        Xh0 += rh0;
    }

    err_out.close();

    rh0.setZero();
    VectorXd ph1 = VectorXd::Constant(p2(N1), 500);
    Prolongation(ph1, rh0, N1, N0, X01_weight, Y01_weight);

    cout << "------------------------------------------------------------------ll" << endl;
    from_X1_to_X(X, Xh0, N0 + 2);
}

void out_tec()
{

    using Basic::N;

    using Basic::X_position;
    using Basic::Y_position;

    using Solve_equ::phi_ana;
    using Solve_equ::X;

    VectorXd error_num2ana = VectorXd::Constant(p2(N), 0);
    VectorXd error_num2ana_L2 = VectorXd::Constant(p2(N), 0);

    std::ostringstream name;
    //name << "Phi_" << N << "_.dat";
    if (Basic::uni_mesh == 1)
    {
        name << "UniMesh_Phi_" << N << "_.dat";
    }
    else
    {
        name << "Un_uniMesh_Phi_" << N << "_.dat";
    }
    std::ofstream out(name.str().c_str());
    out << "Title= \"Poisson_\"\n"
        << "VARIABLES = \"X\", \"Y\", \"phi_num\", \"phi_ana\", \"error\", \"rela_error\" \n";
    out << "ZONE T= \"BOX\",I=" << N - 2 << ",J=" << N - 2 << ", F = POINT" << endl;
    for (int j = 1; j < N - 1; ++j)
    {
        for (int i = 1; i < N - 1; ++i)
        {

            const int index = coordinate(i, j, N);

            const int NodeCondition = get_NodeCondition(i, j, N);

            switch (NodeCondition)
            {
            case 0:
                //error_num2ana(index) = sqrt(p2(X(index) - phi_ana(index)));
                error_num2ana = (X - phi_ana).cwiseAbs();
                error_num2ana_L2(index) = error_num2ana(index) / phi_ana(index);
                break;
            default:
                break;
            }

            out << X_position[i] << " " << Y_position[j] << " "
                << X(index) << " "
                << phi_ana(index) << " "
                << error_num2ana(index) << " "
                << error_num2ana_L2(index) << " "
                //<< phi_ana(index) << " "
                << endl;
        }
    }
    out.close();
}

int main()
{
    cout << "Hello wo!" << endl;

    get_Bacis_Mesh(Basic::uni_mesh, Basic::N, Basic::S_x, Basic::S_y,
                   Basic::X_delta, Basic::Y_delta, Basic::X_position, Basic::Y_position);
    get_Analytic_Solution();

    clock_t start = clock();
    //// --------------- Solve AX = B directly -------------
    //solve_AX_B();

    // Gauss_Seidel_Iteration(Solve_equ::X, Solve_equ::X1, Solve_equ::B1,
    //                        Basic::X_delta, Basic::Y_delta, Basic::X_position, Basic::Y_position, Basic::N);

    MultiGrid_Iter(Basic::N, Solve_equ::X, Solve_equ::X1, Solve_equ::B1,
                   Basic::uni_mesh, Basic::S_x, Basic::S_y,
                   Basic::X_delta, Basic::Y_delta,
                   Basic::X_position, Basic::Y_position);

    clock_t end = clock();

    Solve_equ::L2_error = ((Solve_equ::phi_ana - Solve_equ::X).lpNorm<1>()) / Solve_equ::phi_ana.lpNorm<1>();

    cout << "N = " << Basic::N << ", ";
    cout << "Time : " << double(end - start) << " ms ,";
    cout << "L1 error = " << Solve_equ::L2_error << endl;

    std::ofstream time_out("Steepest_descent.txt", ios::app);
    time_out << "N = " << Basic::N << ", ";
    time_out << "Time : " << double(end - start) << " ms ,";
    time_out << "L1 error = " << Solve_equ::L2_error << endl;
    time_out.close();

    out_tec();
}
