#include <iostream>
#include <vector>

using namespace std;

vector<double> solveLinearSystem(vector<vector<double>>& A) {
    int n = A.size();
    vector<double> x(n, 0.0);
    vector<double> z(n, 0.0);
    vector<vector<double>> L(n, vector<double>(n, 0.0));
    vector<vector<double>> U(n, vector<double>(n, 0.0));

    // Step 1
    L[0][0] = A[0][0];
    U[0][1] = A[0][1] / L[0][0];
    z[0] = A[0][n]/L[0][0];

    // Step 2
    for (int i = 1; i <= n - 2; i++) {
        L[i][i-1] = A[i][i-1];
        z[i] = A[i][n] -L[i][i-1]*U[i-1][i] ;
        L[i][i] = A[i][i] - L[i][i - 1] * U[i - 1][i];
        U[i][i + 1] = A[i][i + 1] / L[i][i];
    }

    // Step 3
    L[n - 1][n - 2] = A[n - 1][n - 2];
    L[n - 1][n - 1] = A[n - 1][n - 1] - L[n - 1][n - 2] * U[n - 2][n - 1];
    z[n-1] = (A[n-1][n] -L[n-1][n-2]*z[n-2] )/ L[n-1][n-1];

    x[n-1] = z[n-1] ;
    // Step 5
    for (int i = 1; i < n; i++) {
        z[i] = (A[i][n] - L[i][i - 1] * z[i - 1]) / L[i][i];
    }

    // Step 6
    x[n - 1] = z[n - 1];
    for (int i = n - 2; i >= 0; i--) {
        x[i] = z[i] - U[i][i + 1] * x[i + 1];
    }

    return x;
}

int main() {
    int n;
    cout << "Enter the dimension n: ";
    cin >> n;

    vector<vector<double>> A(n, vector<double>(n + 1, 0.0));

    cout << "Enter the entries of matrix A:" << endl;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n + 1; j++) {
            cin >> A[i][j];
        }
    }

    vector<double> solution = solveLinearSystem(A);

    cout << "Solution: ";
    for (double x : solution) {
        cout << x << " ";
    }
    cout << endl;

    return 0;
}
