#include <stdlib.h>
#include <stdio.h>
#include <math.h>

// Performs the Jacobi eigenvalue algorithm until the largest absolute
// off-diagonal element is smaller than Epsilon.
// The matrices S and U^T of size NxN are stored in row major memory layout.
// The runtime of the algorithm is in O(N^2) per iteration.
static int Jacobi(double *S, double *UT, int N, double Epsilon) {
    // Set U^T to identity matrix
    for(int K = 0; K < N*N; ++K) {
        UT[K] = 0.0;
    }
    for(int K = 0; K < N; ++K) {
        UT[K*N + K] = 1.0;
    }
    
    for(int It = 0, I, J;; ++It) {
        // Seek largest (absolute) off-diagonal element
        // Initialize with lower right most off-diagonal
        I = N - 2;
        J = N - 1;
        for(int R = 0; R < N - 2; ++R) {
            for(int C = R + 1; C < N; ++C) {
                if(fabs(S[R*N + C]) > fabs(S[I*N + J])) {
                    I = R;
                    J = C;
                }
            }
        }
        
        // Done if largest element is smaller than epsilon
        if(fabs(S[I*N + J]) < Epsilon) {
            return It;
        }
        
        // Compute sine and cosine (pick larger root because why not)
        double Sii = S[I*N + I];
        double Sjj = S[J*N + J];
        double Sij = S[I*N + J];
        double Cot2A = (Sjj - Sii)/(2.0*Sij);
        double TanA = -Cot2A + sqrt(Cot2A*Cot2A + 1.0);
        double CosA = 1.0/sqrt(1.0 + TanA*TanA);
        double SinA = TanA*CosA;
        double CSq = CosA*CosA;
        double SSq = SinA*SinA;
        double SC = SinA*CosA;
        
        // S <- R*S*R^T
        S[I*N + I] = CSq*Sii - 2*SC*Sij + SSq*Sjj;
        S[J*N + J] = SSq*Sii + 2*SC*Sij + CSq*Sjj;
        S[I*N + J] = 0.0;
        S[J*N + I] = 0.0;
        for(int K = 0; K < N; ++K) {
            if(K != I && K != J) {
                double Sik = S[I*N + K];
                double Sjk = S[J*N + K];
                S[I*N + K] = S[K*N + I] = CosA*Sik - SinA*Sjk;
                S[J*N + K] = S[K*N + J] = SinA*Sik + CosA*Sjk;
            }
        }
        
        // U^T <- U^T*R^T (U^T will contain the eigenvectors)
        for(int K = 0; K < N; ++K) {
            double UTki = UT[K*N + I];
            double UTkj = UT[K*N + J];
            UT[K*N + I] = CosA*UTki - SinA*UTkj;
            UT[K*N + J] = CosA*UTkj + SinA*UTki;
        }
    }
}

static void PrintMatrix(char *Title, double *M, int N) {
    printf("%s\n", Title);
    for(int I = 0; I < N; ++I) {
        for(int J = 0; J < N; ++J) {
            printf("%.9f ", M[I*N + J]);
        }
        printf("\n");
    }
    printf("\n");
}

// Computes M = A*B
static void MatrixMultiply(double *M, double *A, double *B, int N) {
    for(int R = 0; R < N; ++R) {
        for(int C = 0; C < N; ++C) {
            double D = 0.0;
            for(int I = 0; I < N; ++I) {
                D += A[R*N + I]*B[I*N + C];
            }
            M[R*N + C] = D;
        }
    }
}

// Computes M = A*B^T
static void MatrixMultiplyBT(double *M, double *A, double *B, int N) {
    for(int R = 0; R < N; ++R) {
        for(int C = 0; C < N; ++C) {
            double D = 0.0;
            for(int I = 0; I < N; ++I) {
                D += A[R*N + I]*B[C*N + I];
            }
            M[R*N + C] = D;
        }
    }
}

// Computes M = A^T*B
static void MatrixMultiplyAT(double *M, double *A, double *B, int N) {
    for(int R = 0; R < N; ++R) {
        for(int C = 0; C < N; ++C) {
            double D = 0.0;
            for(int I = 0; I < N; ++I) {
                D += A[I*N + R]*B[I*N + C];
            }
            M[R*N + C] = D;
        }
    }
}

static void MatrixCopy(double *D, double *S, int N) {
    for(int I = 0; I < N*N; ++I) {
        D[I] = S[I];
    }
}

static void MatrixSetZero(double *M, int N) {
    for(int I = 0; I < N*N; ++I) {
        M[I] = 0.0;
    }
}

static void MatrixSetIdentity(double *M, int N) {
    MatrixSetZero(M, N);
    for(int I = 0; I < N; ++I) {
        M[I*N + I] = 1.0;
    }
}

static double MatrixMaxError(double *M, double *A, int N) {
    double MaxError = 0.0;
    for(int I = 0; I < N*N; ++I) {
        double Error = fabs(M[I] - A[I]);
        if(MaxError < Error) {
            MaxError = Error;
        }
    }
    return MaxError;
}

// Returns G(I, J, Angle)
static void MatrixGivens(double *G, int I, int J, double Angle, int N) {
    MatrixSetIdentity(G, N);
    
    double CosA = cos(Angle);
    double SinA = sin(Angle);
    G[I*N + I] = CosA;
    G[I*N + J] = -SinA;
    G[J*N + I] = SinA;
    G[J*N + J] = CosA;
}

static void TestJacobi(double *S, double *Sp, double *UT, double *UTSp, double *UTSpU, int N) {
    MatrixCopy(Sp, S, N);
    //PrintMatrix("Matrix S:", S, N);
    
    int NumIterations = Jacobi(Sp, UT, N, 1e-12);
    printf("Number of iterations done: %d\n", NumIterations);
    
    //PrintMatrix("Matrix S':", Sp, N);
    //PrintMatrix("Matrix U^T:", UT, N);
    
    MatrixMultiply(UTSp, UT, Sp, N);
    MatrixMultiplyBT(UTSpU, UTSp, UT, N);
    //PrintMatrix("Matrix U^T*S'*U", UTSpU, N);
    printf("MaxError(U^T*S'*U - S): %.12f\n", MatrixMaxError(UTSpU, S, N));
}

static void TestRandom() {
#define MaxN (20)
#define C (MaxN*MaxN)
    double S[C], G[C], Tmp[C], U[C], UT[C], Sp[C], UTSp[C], UTSpU[C];
    
    for(int TestIt = 0; TestIt < 100; ++TestIt) {
        int N = (rand() % (MaxN - 2) + 2);
        
        MatrixSetZero(Sp, N);
        MatrixSetIdentity(U, N);
        
        for(int I = 0; I < N; ++I) {
            Sp[I*N + I] = (double)rand()/((double)rand() + 1.0);
        }
        for(int It = 0; It < 100; ++It) {
            int I = rand()%N;
            int J;
            do {
                J = rand()%N;
            } while(I == J);
            if(J < I) {
                int TmpI = I;
                I = J;
                J = TmpI;
            }
            double Angle = (double)rand()/((double)rand() + 1.0);
            MatrixGivens(G, I, J, Angle, N);
            
            //printf("%d, I=%d, J=%d, A=%f", It, I, J, Angle);
            //PrintMatrix("", G, N);
            
            MatrixCopy(Tmp, U, N);
            MatrixMultiply(U, G, Tmp, N);
        }
        
        MatrixMultiply(Tmp, Sp, U, N);
        MatrixMultiplyAT(S, U, Tmp, N);
        
        printf("N: %d\n", N);
        //PrintMatrix("Original S':", Sp, N);
        //PrintMatrix("Original U:", U, N);
        
        TestJacobi(S, Sp, UT, UTSp, UTSpU, N);
    }
}

static void TestWikipediaExample() {
    int N = 4;
    double S[] = {
        4, -30, 60, -35,
        -30, 300, -675, 420,
        60, -675, 1620, -1050,
        -35, 420, -1050, 700,
    };
    double Sp[16], UT[16];
    double UTSp[16];
    double UTSpU[16];
    TestJacobi(S, Sp, UT, UTSp, UTSpU, N);
}

int main(void) {
    TestRandom();
    return 0;
}
