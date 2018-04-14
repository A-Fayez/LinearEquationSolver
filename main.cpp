#include <iostream>
#include <cmath>
#define N 10 // Number of variables
#define M 11 // Number of columns for matrix (A|b)
using namespace std;
double x[N];
double adj[N][N];
double inv[N][N];
double A[N][N];
double b[N];
double G[N][M]; // That matrix represents A|b matrix
 // result of multiplying matrices
// global to initialize its values by 0
void multiply(double a1[N][N], double a2[N],double x[N]) // multiply two matrices
{

    for (int i=0;i<N;i++)
    {
        for (int j=0;j<N;j++)
        {
            x[i]+=(a1[i][j]*a2[j]) ;
        }
    }
}

double det(double arr[N][N] ,int n) {//det of 2x2 matrix
    double determ = 0;//determ variable to store determinant
    double subarr[N][N];
    if (n==2)
    {
    determ+=(arr[0][0]*arr[1][1]) - arr[1][0]*arr[0][1];
    return determ;
    }
    for (int k=0;k<n;k++)
        {
            int subi=0;
            for (int i=1;i<n;i++)
            {
                int subj=0;
                for (int j=0;j<n;j++)
                {
                    if (j==k) continue;

                    subarr[subi][subj] = arr[i][j];
                    subj++;
                }
                subi++;
            }
        determ=determ+ ( pow(-1 ,k) * arr[0][k] * det(subarr, n-1));
        }

     return determ;
}

void factors(double A[N][N], double temp[N][N], int k, int m, int n) // function for getting cofactors of elment A[k][m]
{
    int i = 0, j = 0;
    //for each element of the matrix
    for (int row = 0; row < n; row++)
    {
        for (int col = 0; col < n; col++)
        {
            if (row != k && col != m) //skipping the row and col of the element A[row][col]
            {
                temp[i][j++] = A[row][col];

                if (j == n - 1)
                {
                    j = 0;
                    i++;
                }
            }
        }
    }
}


void adjoint(double A[N][N],double adj[N][N])
{
    // temp is used to store cofactors of A[i][j]
    double temp[N][N];
    for (int i=0; i<N; i++)
    {
        for (int j=0; j<N; j++)
        {
            factors(A, temp, i, j, N); // temp is filled
            adj[j][i] = (pow(-1,i+j))*(det(temp, N-1)); //Transpose and multiply by sub matrix
        }
    }
}


void inverse(double A[N][N],double inv[N][N]) {// filling the inv matrix when passing adj matrix
    double d = det(A,N);
    for (int i=0;i<N;i++) {
        for (int j=0;j<N;j++)  {
            inv[i][j] = (adj[i][j]/d) ;
        }
    }

}

void swap(double mat[N][N], int row1, int row2,int col)
{

    for (int i = 0; i < col; i++)
    {
        int temp = mat[row1][i];
        mat[row1][i] = mat[row2][i];
        mat[row2][i] = temp;
    }
}
void swapG(double mat[N][M], int row1, int row2,int col)
{
    for (int i = 0; i < col; i++)
    {
        int temp = mat[row1][i];
        mat[row1][i] = mat[row2][i];
        mat[row2][i] = temp;
    }
}
int rankof(double mat[N][N]) // used online help here
{
    int rank = N;

    for (int row = 0; row < rank; row++)
    {
        // Before we visit current row 'row', we make
        // sure that mat[row][0],....mat[row][row-1]
        // are 0.

        // Diagonal element is not zero
        if (mat[row][row])
        {
           for (int col = 0; col < N; col++)
           {
               if (col != row)
               {
                 // This makes all entries of current
                 // column as 0 except entry 'mat[row][row]'
                 double mult = (double)mat[col][row] /
                                       mat[row][row];
                 for (int i = 0; i < rank; i++)
                   mat[col][i] -= mult * mat[row][i];
              }
           }
        }

        // Diagonal element is already zero. Two cases
        // arise:
        // 1) If there is a row below it with non-zero
        //    entry, then swap this row with that row
        //    and process that row
        // 2) If all elements in current column below
        //    mat[r][row] are 0, then remvoe this column
        //    by swapping it with last column and
        //    reducing number of columns by 1.
        else
        {
            bool reduce = true;

            /* Find the non-zero element in current
                column  */
            for (int i = row + 1; i < N;  i++)
            {
                // Swap the row with non-zero element
                // with this row.
                if (mat[i][row])
                {
                    swap(mat, row, i, rank);
                    reduce = false;
                    break ;
                }
            }

            // If we did not find any row with non-zero
            // element in current columnm, then all
            // values in this column are 0.
            if (reduce)
            {
                // Reduce number of columns
                rank--;

                // Copy the last column here
                for (int i = 0; i < N; i ++)
                    mat[i][row] = mat[i][rank];
            }

            // Process this row again
            row--;
        }

       // Uncomment these lines to see intermediate results
       // display(mat, R, C);
       // printf("\n");
    }
    return rank;
}

int rankofG(double mat[N][M]) // used online help here
{
    int rank = M;

    for (int row = 0; row < rank; row++)
    {
        // Before we visit current row 'row', we make
        // sure that mat[row][0],....mat[row][row-1]
        // are 0.

        // Diagonal element is not zero
        if (mat[row][row])
        {
           for (int col = 0; col < M; col++)
           {
               if (col != row)
               {
                 // This makes all entries of current
                 // column as 0 except entry 'mat[row][row]'
                 double mult = (double)mat[col][row] /
                                       mat[row][row];
                 for (int i = 0; i < rank; i++)
                   mat[col][i] -= mult * mat[row][i];
              }
           }
        }

        // Diagonal element is already zero. Two cases
        // arise:
        // 1) If there is a row below it with non-zero
        //    entry, then swap this row with that row
        //    and process that row
        // 2) If all elements in current column below
        //    mat[r][row] are 0, then remvoe this column
        //    by swapping it with last column and
        //    reducing number of columns by 1.
        else
        {
            bool reduce = true;

            /* Find the non-zero element in current
                column  */
            for (int i = row + 1; i < M;  i++)
            {
                // Swap the row with non-zero element
                // with this row.
                if (mat[i][row])
                {
                    swapG(mat, row, i, rank);
                    reduce = false;
                    break ;
                }
            }

            // If we did not find any row with non-zero
            // element in current columnm, then all
            // values in this column are 0.
            if (reduce)
            {
                // Reduce number of columns
                rank--;

                // Copy the last column here
                for (int i = 0; i < M; i ++)
                    mat[i][row] = mat[i][rank];
            }

            // Process this row again
            row--;
        }

       // Uncomment these lines to see intermediate results
       // display(mat, R, C);
       // printf("\n");
    }
    return rank;
}

int main()
{
    cout << "For solving system of linear equations as follows   " << endl;
    cout << "a1X1+a2X2+....+a10X10 = z1" << endl;
    cout << "As Ax = b please enter the matrix A which represents cofactors of the 10 variables " << endl << endl;

   for (int i=0;i<N;i++) // entering elements of matrix A
   {
       for (int j=0;j<N;j++)
       {
           cin >> A[i][j];
       }
   }

   cout <<" Enter matrix b which represents Z(n) (the resulting factors )" <<endl << endl;

   for (int i=0;i<N;i++) // entering elements of matrix B
   {
      cin>>b[i];
   }
cout << endl << endl;
cout << "The x vector which reprsesent the system solution is" << endl;

// Making matrix G which represents A|b
for (int i=0;i<N;i++) // entering elements of matrix A
   {
       for (int j=0;j<N+1;j++)
       {
           if (j==N) G[i][N] = b[i];
           else G[i][j] = A[i][j];
       }
   }
adjoint(A,adj); //adj matrix is filled
inverse(A,inv); //inv matrix is filled
multiply(inv,b,x); //x vector is filled
if (  rankof(A)==rankofG(G) && rankof(A)==N    )
{
cout << " The System has one solution and it is " << endl;
for (int i=0;i<N;i++) cout << x[i] << endl;
}


else if (rankof(A)==rankofG(G) && rankof(A)<N ) cout << " The System has infinite solutions " << endl;
else if (rankof(A)<rankofG(G)) cout << "The system has no solution" << endl;

        return 0;
}


