
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

int coordinate(const int x, const int y, const int N){
    return(x + y * N);
}

void Gaussian_Elimination(double* A, double* B, double* res, const int N){      
    for(int j = 0; j < N - 1; ++j){                                             
        for(int i = j + 1; i < N; ++i){                                         
            double f_eli = A[coordinate(j, i, N)] / A[coordinate(j, j, N)];               
                                                                                
            B[i] = B[i] - f_eli * B[j];                                         
            for(int k = 0; k < N; ++k){                                         
                A[coordinate(k, i, N)] = A[coordinate(k, i, N)] - f_eli * A[coordinate(k, j, N)];
            }                                                                   
        }                                                                       
    }                                                                           
                                                                                
    for(int j = N - 1; j > -1; --j){                                            
        res[j] = B[j];                                                          
        for(int i = j + 1; i < N; ++i){                                         
            if(i != j){                                                         
                res[j] -= A[coordinate(i, j, N)] * res[i];                           
            }                                                                   
        }                                                                       
        res[j] = res[j] / A[coordinate(j, j, N)];                                    
    }                                                                           
}       

int main(){
    const double L = 1.0;
    const int N = 21;

    // ******* Exponential mesh**************************
    const double S = 1.5;
    const double x0 = L * (1 - S) / (1 - pow(S, N - 1));

    //// ******* Uniform mesh**************************
    //const double S = 1.00;
    //const double x0 = L / (N - 1);

    printf("The first delta_x is : x0 = %f\n", x0);

    double* Coefficient_Matrix = (double*) malloc(N * N * sizeof(double));

    double* Right_Term = (double*) malloc((N    ) * sizeof(double));
    double* X_Posi     = (double*) malloc((N    ) * sizeof(double));
    double* Delta_x    = (double*) malloc((N - 1) * sizeof(double)); 

    double* phi_num    = (double*) malloc((N    ) * sizeof(double));
    double* phi_ana    = (double*) malloc((N    ) * sizeof(double));

    //**********We need get delta_x, position of each node,***********
    Delta_x[0] = x0; 
    phi_ana[0] = 0;
    X_Posi[0] = 0;
    Right_Term[0] = 2 * X_Posi[0] - 1; 

    for(int i = 1; i < N; ++i){

        if(i < N - 1){
            Delta_x[i] = S * Delta_x[i - 1];

        }

        X_Posi[i] = Delta_x[i - 1] + X_Posi[i - 1];

        phi_ana[i] = 1.0 / 3.0 * X_Posi[i] * X_Posi[i] * X_Posi[i] 
                   - 1.0 / 2.0 * X_Posi[i] * X_Posi[i] 
                   + 7.0 / 6.0 * X_Posi[i];

        Right_Term[i] = 2 * X_Posi[i] - 1; 
    }


    //*********** Constract Coeffieient Matrix of Dirichlet Boundary condiction
    for(int j = 1; j < N - 1; ++j){
        for(int i = 0; i < N ; ++i){
            const int coordinate_ = coordinate(i, j, N);
            Coefficient_Matrix[coordinate_] = 0;
            if(i == j){
                Coefficient_Matrix[coordinate_] = -(1.0 + S);
            }
            if(i + 1 == j){
                Coefficient_Matrix[coordinate_] = S;
            }
            if(i - 1 == j){
                Coefficient_Matrix[coordinate_] = 1.0;
            }
        }
    }
    for(int i = 0; i < N; ++i){
        const int coordinate_T = coordinate(i, 0    , N);
        const int coordinate_B = coordinate(i, N - 1, N);
        Coefficient_Matrix[coordinate_T] = 0;
        Coefficient_Matrix[coordinate_B] = 0;
    }
    Coefficient_Matrix[coordinate(0    , 0    , N)] = 1;
    Coefficient_Matrix[coordinate(N - 1, N - 1, N)] = 1;
    //***End of constracting Coeffieient Matrix of Dirichlet Boundary condiction

    //***** Reconstract Right hand term of linear algebraic equations
    for(int i = 1; i < N - 1; ++i){
        Right_Term[i] *= Delta_x[i - 1] * Delta_x[i - 1] * S * (1 + S) / 2.0;
    }
    Right_Term[0    ] = 0;
    Right_Term[N - 1] = 1;
    //**End of reconstracting Right hand term of linear algebraic equations

    Gaussian_Elimination(Coefficient_Matrix, Right_Term, phi_num, N);


    //*** Output results*********************
    // std::ostringstream name;
	// name << N << "_nodes_results" << ".dat";
    // std::ofstream fout(name.str().c_str( ) );

    // double* abs_err  = (double*) malloc((N) * sizeof(double));
    // double* rela_err  = (double*) malloc((N) * sizeof(double));

    // fout << "i Position Right_Term Analysitic Numerical Abs_error Relative_error" << endl;
    // for(int i = 0; i < N; ++i){
    //     abs_err[i] = fabs(phi_ana[i] - phi_num[i]);
    //     rela_err[i] = sqrt((phi_ana[i] - phi_num[i]) * (phi_ana[i] - phi_num[i]) / phi_ana[i]);
    //     //fout << std::setw(4) << i << " "
    //     fout << i << " "
    //          << std::fixed  << std::setprecision(10)
    //          << X_Posi[i] << " "
    //          << Right_Term[i] << " "
    //          << phi_ana[i] << " "
    //          << phi_num[i] << " "
    //          << abs_err[i] << " "
    //          << rela_err[i] << " "
    //          << endl;
    // }
    // fout.close();

    // free(abs_err);
    // free(rela_err);
    free(Coefficient_Matrix);
    free(Right_Term );
    free(X_Posi     );
    free(Delta_x    ); 
    free(phi_num    );
    free(phi_ana    );
}
