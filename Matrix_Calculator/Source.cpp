#include<iostream>
#include<fstream>
#include<vector>
#include<ctime>
using namespace std;
ifstream file("D:\\data.txt");

// -i*4+2^4+7*i+3E2-12/4+5 = 318+3*i
// 5+(4*(12-5)-30/(4+6)+4)+7-(15+2^5+(4-1)^(7-5)-13)+17 = 15

const double eps = 1e-6;

template<typename T1, typename T2> void operator +=(T1& A, const T2& B)
{
	A = A + B;
}
template<typename T1, typename T2> void operator -=(T1& A, const T2& B)
{
	A = A - B;
}
template<typename T1, typename T2> void operator *=(T1& A, const T2& B)
{
	A = A * B;
}
template<typename T1, typename T2> void operator /=(T1& A, const T2& B)
{
	A = A / B;
}
template<typename T> T operator -(const T& A)
{
	return -1 * A;
}
template<typename T> T operator -(const T& A, const T& B)
{
	return A + (-B);
}

void timer()
{
	int L[12] = { 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };
	const char* X[5], * W[7] = { "MONDAY", "TUESDAY", "WEDNESDAY", "THURSDAY", "FRIDAY", "SATURDAY", "SUNDAY" };
	int T[7];
	T[0] = time(NULL);
	T[1] = T[0] / 60;
	T[0] %= 60;
	T[2] = T[1] / 60 + 3;
	T[1] %= 60;
	T[3] = T[2] / 24;
	T[2] %= 24;
	T[6] = (T[3] + 3) % 7;
	T[5] = 1970;
	while (T[3] >= 1461)
	{
		T[3] -= 1461;
		T[5] += 4;
	}
	T[4] = 0;
	while (T[3] >= L[T[4]])
	{
		T[3] -= L[T[4]];
		T[4]++;
		T[4] %= 12;
		if (T[4] == 0)
		{
			T[5]++;
			if (T[5] % 4 == 0)
				L[1] = 29;
			else
				L[1] = 28;
		}
	}
	T[3]++;
	T[4]++;
	for (int i = 0; i < 5; i++)
		if (T[i] < 10)
			X[i] = "0";
		else
			X[i] = "";
	cout << X[3] << T[3] << "." << X[4] << T[4] << "." << T[5] << " | "
		<< X[2] << T[2] << ":" << X[1] << T[1] << ":" << X[0] << T[0] << " | " << W[T[6]];
}

class complex
{
private:
	int out(double X)
	{
		bool b = false;
		int i, j, M[7], N = 0, S = 11;
		if (X < 0)
			X = -X;
		if (X > 1)
			while (X >= 10)
			{
				X /= 10;
				N++;
			}
		else if (X != 0)
			while (X < 1)
			{
				X *= 10;
				N--;
			}
		for (i = 0; i < 7; i++)
		{
			M[i] = X;
			X -= M[i];
			X *= 10;
		}
		if (M[6] >= 5)
			for (i = 5; i >= 0; i--)
			{
				M[i]++;
				if (M[i] != 10)
					break;
				else
					M[i] = 0;
			}
		if (M[0] == 0)
		{
			M[0] = 1;
			N++;
		}
		if (N >= 6)
		{
			cout << M[0];
			S--;
			for (i = 1; i < 6; i++)
			{
				b = true;
				for (j = i; j < 6; j++)
					if (M[j] != 0)
						b = false;
				if (b)
					break;
				if (i == 1)
				{
					cout << ".";
					S--;
				}
				cout << M[i];
				S--;
			}
			cout << "E" << N;
			S--;
			for (i = 1; i <= N; i *= 10)
				S--;
		}
		else if (N >= 0)
		{
			for (i = 0; i < 6; i++)
			{
				cout << M[i];
				S--;
				if (i >= N)
				{
					b = true;
					for (j = i + 1; j < 6; j++)
						if (M[j] != 0)
							b = false;
				}
				if (b)
					break;
				if (i == N)
				{
					cout << ".";
					S--;
				}
			}
		}
		else if (N >= -6)
		{
			if (M[7 + N] >= 5 && N != -1)
				for (i = 6 + N; i >= 0; i--)
				{
					M[i]++;
					if (M[i] != 10)
						break;
					else
						M[i] = 0;
				}
			if (M[0] == 0)
			{
				M[0] = 1;
				N++;
				for (i = 1; i < 7; i++)
					M[i] = 0;
			}
			cout << "0.";
			S -= 2;
			for (i = -1; i > N; i--)
			{
				cout << "0";
				S--;
			}
			for (i = 0; i < 7 + N; i++)
			{
				cout << M[i];
				S--;
				b = true;
				for (j = i + 1; j < 7 + N; j++)
					if (M[j] != 0)
						b = false;
				if (b)
					break;
			}
		}
		else
		{
			cout << "0";
			S--;
		}
		return S;
	}
public:
	double Re, Im;
	complex(double Re, double Im)
	{
		this->Re = Re;
		this->Im = Im;
	}
	complex(double Re)
	{
		this->Re = Re;
		Im = 0;
	}
	complex()
	{
		Re = 0;
		Im = 0;
	}
	void output()
	{
		double X = Re, Y = Im;
		if (X > -eps && X < eps)
			X = 0;
		if (Y > -eps && Y < eps)
			Y = 0;
		if (Y == 0)
			if (X == 0)
				cout << 0;
			else
				cout << X;
		else
			if (X == 0)
				cout << Y << "*i";
			else
			{
				cout << X;
				if (Y > 0)
					cout << "+";
				cout << Y << "*i";
			}
	}
	void output_tab()
	{
		int S = 1;
		double X = Re, Y = Im;
		if (X > -eps && X < eps)
			X = 0;
		if (Y > -eps && Y < eps)
			Y = 0;
		if (Y == 0)
			if (X == 0)
			{
				cout << " 0";
				S += 24;
			}
			else
			{
				if (X > 0)
					cout << " ";
				else
					cout << "-";
				S += out(X) + 14;
			}
		else
			if (X == 0)
			{
				if (Y > 0)
					cout << " ";
				else
					cout << "-";
				S += 12 + out(Y);
				cout << "*i";
			}
			else
			{
				if (X > 0)
					cout << " ";
				else
					cout << "-";
				S += out(X);
				if (Y > 0)
					cout << "+";
				else
					cout << "-";
				S += out(Y);
				cout << "*i";
			}
		for (int i = 0; i < S; i++)
			cout << " ";
	}
	void output_endl()
	{
		output_tab();
		cout << "\n\n";
	}
	double mod_2()
	{
		return Re * Re + Im * Im;
	}
	double mod()
	{
		return sqrt(mod_2());
	}
	friend complex operator ~(complex Z)
	{
		return complex(Z.Re, -Z.Im);
	}
	friend complex operator +(complex Z1, complex Z2)
	{
		return complex(Z1.Re + Z2.Re, Z1.Im + Z2.Im);
	}
	friend complex operator -(complex Z1, complex Z2)
	{
		return complex(Z1.Re - Z2.Re, Z1.Im - Z2.Im);
	}
	friend complex operator *(complex Z1, complex Z2)
	{
		return complex(Z1.Re * Z2.Re - Z1.Im * Z2.Im, Z1.Re * Z2.Im + Z1.Im * Z2.Re);
	}
	friend complex operator /(complex Z1, complex Z2)
	{
		double X = Z2.mod_2();
		return complex((Z1.Re * Z2.Re + Z1.Im * Z2.Im) / X, (Z1.Im * Z2.Re - Z1.Re * Z2.Im) / X);
	}
	complex exp(complex Z)
	{
		complex W = *this;
		int i, S, X = 0, Y = 1;
		if (Z.Re < 0)
		{
			Z.Re = -Z.Re;
			W = 1 / W;
		}
		vector<complex> D({ W });
		for (i = 0; Y <= Z.Re; i++)
		{
			D.push_back(D[i] * D[i]);
			Y *= 2;
		}
		W = 1;
		for (i = D.size() - 1; i >= 0; i--)
		{
			S = X + Y;
			if (S <= Z.Re)
			{
				X = S;
				W *= D[i];
			}
			Y /= 2;
		}
		return W;
	}
	complex operator ^=(complex Z)
	{
		return *this = exp(Z);
	}
	friend bool operator ==(complex Z1, complex Z2)
	{
		return ((Z1 - Z2).mod() < eps);
	}
	friend bool operator !=(complex Z1, complex Z2)
	{
		return ((Z1 - Z2).mod() > eps);
	}
	friend bool operator >(complex Z1, complex Z2)
	{
		return (Z1.mod() > Z2.mod() + eps);
	}
	friend bool operator <(complex Z1, complex Z2)
	{
		return (Z1.mod() + eps < Z2.mod());
	}
};

class matrix
{
private:
	void Rayleigh(complex Z)
	{
		bool X = true;
		int i;
		complex W;
		matrix A, B, E = exp(0);
		for (i = 0; i < table.size(); i++)
			A.table.push_back({ 1 });
		i = 0;
		while (X)
		{
			i++;
			W = Z;
			B = (*this - Z * E).exp(-1) * A;
			A = B / B.F_norm();
			Z = (~A * *this * A) % (~A * A);
			if (Z == W)
				break;
			if (i > 10)
				X = false;
		}
		if (X)
		{
			eigenvalue.push_back(Z);
			eigenvector.push_back(A);
		}
	}
public:
	vector<vector<complex>> table;
	vector<complex> eigenvalue;
	vector<matrix> eigenvector;
	matrix() {}
	void input()
	{
		bool X;
		int M, N, i, j;
		double Re, Im = 0;
		vector<complex> row;
		table.clear();
		file >> X >> M >> N;
		for (i = 0; i < M; i++)
		{
			for (j = 0; j < N; j++)
			{
				file >> Re;
				if (X)
					file >> Im;
				row.push_back(complex(Re, Im));
			}
			table.push_back(row);
			row.clear();
		}
	}
	void output()
	{
		cout << " ============================================  M  A  T  R  I  X  ============================================\n\n";
		for (int i = 0; i < table.size(); i++)
		{
			cout << " ";
			for (int j = 0; j < table[0].size(); j++)
				table[i][j].output_tab();
			cout << "\n\n";
		}
		cout << " ============================================================================================================\n\n";
	}
	void output_close()
	{
		cout << " ============================================  M  A  T  R  I  X  ============================================\n\n";
		for (int i = 0; i < table.size(); i++)
		{
			cout << " ";
			for (int j = 0; j < table[0].size(); j++)
			{
				cout << " ";
				table[i][j].output();
			}
			cout << "\n\n";
		}
		cout << " ============================================================================================================\n\n";
	}
	friend matrix operator &(matrix A)
	{
		vector<complex> row;
		matrix B;
		for (int j = 0; j < A.table[0].size(); j++)
		{
			for (int i = 0; i < A.table.size(); i++)
				row.push_back(A.table[i][j]);
			B.table.push_back(row);
			row.clear();
		}
		return B;
	}
	friend matrix operator ~(matrix A)
	{
		for (int i = 0; i < A.table.size(); i++)
			for (int j = 0; j < A.table[0].size(); j++)
				A.table[i][j] = ~A.table[i][j];
		return &A;
	}
	friend matrix operator +(matrix A, matrix B)
	{
		for (int i = 0; i < A.table.size(); i++)
			for (int j = 0; j < A.table[0].size(); j++)
				A.table[i][j] += B.table[i][j];
		return A;
	}
	friend matrix operator -(matrix A, matrix B)
	{
		for (int i = 0; i < A.table.size(); i++)
			for (int j = 0; j < A.table[0].size(); j++)
				A.table[i][j] -= B.table[i][j];
		return A;
	}
	friend matrix operator *(matrix A, matrix B)
	{
		complex cell;
		vector<complex> row;
		matrix C;
		B = &B;
		for (int i = 0; i < A.table.size(); i++)
		{
			for (int j = 0; j < B.table.size(); j++)
			{
				cell = 0;
				for (int k = 0; k < B.table[0].size(); k++)
					cell += A.table[i][k] * B.table[j][k];
				row.push_back(cell);
			}
			C.table.push_back(row);
			row.clear();
		}
		return C;
	}
	friend matrix operator *(matrix A, complex Z)
	{
		for (int i = 0; i < A.table.size(); i++)
			for (int j = 0; j < A.table[0].size(); j++)
				A.table[i][j] *= Z;
		return A;
	}
	friend matrix operator *(complex Z, matrix A)
	{
		return A * Z;
	}
	friend matrix operator /(matrix A, complex Z)
	{
		for (int i = 0; i < A.table.size(); i++)
			for (int j = 0; j < A.table[0].size(); j++)
				A.table[i][j] /= Z;
		return A;
	}
	friend complex operator %(complex Z, matrix A)
	{
		return Z / A.table[0][0];
	}
	friend complex operator %(matrix A, complex Z)
	{
		return A.table[0][0] / Z;
	}
	friend complex operator %(matrix A, matrix B)
	{
		return A.table[0][0] / B.table[0][0];
	}
	matrix exp(int N)
	{
		int i, j, k, S, X = 0, Y = 1;
		matrix A = *this, E = *this;
		for (i = 0; i < E.table.size(); i++)
			for (j = 0; j < E.table.size(); j++)
				if (i == j)
					E.table[i][j] = 1;
				else
					E.table[i][j] = 0;
		if (N < 0)
		{
			N = -N;
			complex cell;
			vector<complex> row;
			A = concat_hor(A, E);
			for (i = 0; i < A.table.size(); i++)
			{
				if (A.table[i][i] == 0)
					for (j = i + 1; j < A.table.size(); j++)
						if (A.table[j][i] != 0)
						{
							row = A.table[i];
							A.table[i] = A.table[j];
							A.table[j] = row;
							break;
						}
				cell = A.table[i][i];
				for (j = 0; j < A.table[0].size(); j++)
					A.table[i][j] /= cell;
				for (j = 0; j < A.table.size(); j++)
					if (i != j)
					{
						cell = A.table[j][i];
						for (k = 0; k < A.table[0].size(); k++)
							A.table[j][k] -= cell * A.table[i][k];
					}
			}
			A = &A;
			A.table.erase(A.table.begin(), A.table.begin() + A.table[0].size());
			A = &A;
		}
		vector<matrix> D({ A });
		for (i = 0; Y <= N; i++)
		{
			D.push_back(D[i] * D[i]);
			Y *= 2;
		}
		A = E;
		for (i = D.size() - 1; i >= 0; i--)
		{
			S = X + Y;
			if (S <= N)
			{
				X = S;
				A *= D[i];
			}
			Y /= 2;
		}
		return A;
	}
	matrix operator ^=(int N)
	{
		return *this = exp(N);
	}
	friend bool operator ==(matrix A, matrix B)
	{
		if (A.table.size() != B.table.size())
			return false;
		if (A.table[0].size() != B.table[0].size())
			return false;
		for (int i = 0; i < A.table.size(); i++)
			for (int j = 0; j < A.table[0].size(); j++)
				if (A.table[i][j] != B.table[i][j])
					return false;
		return true;
	}
	friend bool operator !=(matrix A, matrix B)
	{
		return !(A == B);
	}
	friend matrix concat_ver(matrix A, matrix B)
	{
		for (int i = 0; i < B.table.size(); i++)
			A.table.push_back(B.table[i]);
		return A;
	}
	friend matrix concat_hor(matrix A, matrix B)
	{
		return &concat_ver(&A, &B);
	}
	complex det()
	{
		if (table.size() != table[0].size())
			return 0;
		if (table.size() == 1)
			return table[0][0];
		if (table.size() == 2)
			return table[0][0] * table[1][1] - table[0][1] * table[1][0];
		if (table.size() > 2)
		{
			complex cell, det;
			for (int j = 0; j < table.size(); j++)
			{
				matrix minor = *this;
				minor.table.erase(minor.table.begin());
				minor = &minor;
				minor.table.erase(minor.table.begin() + j);
				minor = &minor;
				cell = table[0][j] * minor.det();
				if (j % 2 == 0)
					det += cell;
				else
					det -= cell;
			}
			return det;
		}
	}
	int rank()
	{
		int i, j, rk = 0;
		complex cell;
		vector<complex> row;
		matrix A = *this;
		while (true)
		{
			if (A.table[0].size() == 1)
				A = &A;
			if (A.table.size() == 1)
			{
				for (i = 0; i < A.table[0].size(); i++)
					if (A.table[0][i] != 0)
					{
						rk++;
						break;
					}
				break;
			}
			if (A.table[0][0] == 0)
				for (i = 1; i < A.table.size(); i++)
					if (A.table[i][0] != 0)
					{
						row = A.table[0];
						A.table[0] = A.table[i];
						A.table[i] = row;
						break;
					}
			if (A.table[0][0] == 0)
			{
				A = &A;
				A.table.erase(A.table.begin());
				A = &A;
			}
			else
			{
				cell = A.table[0][0];
				for (i = 0; i < A.table[0].size(); i++)
					A.table[0][i] /= cell;
				for (i = 1; i < A.table.size(); i++)
					for (j = 0; j < A.table[0].size(); j++)
						A.table[i][j] -= A.table[i][0] * A.table[0][j];
				A.table.erase(A.table.begin());
				A = &A;
				A.table.erase(A.table.begin());
				A = &A;
				rk++;
			}
		}
		return rk;
	}
	double F_norm()
	{
		double norm = 0;
		for (int i = 0; i < table.size(); i++)
			for (int j = 0; j < table[0].size(); j++)
				norm += table[i][j].mod_2();
		return sqrt(norm);
	}
	void grid(int N)
	{
		if (table.size() == table[0].size())
		{
			eigenvalue.clear();
			eigenvector.clear();
			int i, j;
			double X = F_norm();
			complex Z;
			for (i = -N; i <= N; i++)
				for (j = -N; j <= N; j++)
				{
					Z = complex(X * i / N, X * j / N);
					if (Z.mod() < X + eps)
						Rayleigh(Z);
				}
			i = 0;
			while (i < eigenvalue.size())
			{
				Z = eigenvalue[i];
				j = i + 1;
				while (j < eigenvalue.size())
					if (Z == eigenvalue[j])
					{
						eigenvalue.erase(eigenvalue.begin() + j);
						eigenvector.erase(eigenvector.begin() + j);
					}
					else
						j++;
				i++;
			}
		}
	}
	void info()
	{
		cout << " ===============================================  I  N  F  O  ===============================================\n\n";
		cout << "  SIZE " << table.size() << " x " << table[0].size() << "     ";
		cout << "  RANK " << rank() << "\n\n";
		cout << "  FROBENIUS NORM  ";
		complex(F_norm()).output_tab();
		cout << "\n\n";
		if (table.size() == table[0].size())
		{
			cout << "  DETERMINANT     ";
			det().output_tab();
			cout << "\n\n";
		}
		cout << " ============================================================================================================\n\n";
	}
	void eigen()
	{
		if (eigenvalue.size() != 0)
		{
			cout << " ======================================= EIGENVALUES AND EIGENVECTORS =======================================\n\n";
			for (int i = 0; i < eigenvalue.size(); i++)
			{
				cout << "  [" << i + 1 << "] ";
				eigenvalue[i].output_tab();
				cout << "\n\n";
				for (int j = 0; j < eigenvector[0].table.size(); j++)
				{
					cout << "  ";
					eigenvector[i].table[j][0].output_tab();
					cout << "\n\n";
				}
				cout << " ============================================================================================================\n\n";
			}
		}
	}
	friend matrix Kramer(matrix A, matrix B)
	{
		matrix X, Y;
		complex D = A.det();
		for (int i = 0; i < A.table.size(); i++)
		{
			Y = A;
			for (int j = 0; j < A.table.size(); j++)
				Y.table[j][i] = B.table[j][0];
			X.table.push_back({ Y.det() / D });
		}
		return X;
	}
	friend matrix Gauss(matrix A, matrix B)
	{
		int i, j, k, l = 0;
		complex cell;
		vector<complex> row;
		A = concat_hor(A, B);
		for (k = 0; k < A.table.size(); k++)
		{
			if (A.table[k][l] == 0)
				for (i = k + 1; i < A.table.size(); i++)
					if (A.table[i][l] != 0)
					{
						row = A.table[k];
						A.table[k] = A.table[i];
						A.table[i] = row;
						break;
					}
			if (A.table[k][l] == 0)
				l++;
			else
			{
				cell = A.table[k][l];
				for (i = 0; i < A.table[0].size(); i++)
					A.table[k][i] /= cell;
				for (i = 0; i < A.table.size(); i++)
					for (j = 0; j < A.table[0].size(); j++)
						if (i != k)
							A.table[i][j] -= A.table[i][l] * A.table[k][j];
			}
			l++;
		}
		return A;
	}
};

complex computing(char text[], vector<complex> computed, int N)
{
	bool b = true;
	int i, j, next, L = 0;
	double Re;
	vector<int> I, C, S;
	vector<complex> R;
	while (text[L] != '\0')
		L++;
	for (i = 0; i < L; i++)
		if (text[i] == '#')
			C.push_back(i);
		else if (text[i] == 'i')
			I.push_back(i);
		else if (text[i] == '+' || text[i] == '-' || text[i] == '*' || text[i] == '/'
			|| text[i] == '^' || text[i] == 'E')
			S.push_back(i);
	S.push_back(L);
	if (text[S[0]] == '+' || text[S[0]] == '-')
	{
		for (i = 0; i < S[0]; i++)
			if (text[i] != ' ')
			{
				b = false;
				break;
			}
	}
	else
		b = false;
	if (b)
		R.push_back(0);
	else
		S.emplace(S.begin(), -1);
	for (i = 0; i < S.size() - 1; i++)
	{
		b = true;
		ofstream output("D:\\buffer.txt");
		next = i + 1;
		for (j = S[i] + 1; j < S[next]; j++)
			if (text[j] == '#')
			{
				R.push_back(computed[N]);
				N++;
				b = false;
				break;
			}
			else if (text[j] == 'i')
			{
				R.push_back(complex(0, 1));
				b = false;
				break;
			}
			else
				output << text[j];
		output.close();
		if (b)
		{
			ifstream input("D:\\buffer.txt");
			input >> Re;
			R.push_back(Re);
			input.close();
		}
	}
	if (S[0] == -1)
		S.erase(S.begin());
	S.pop_back();
	while (true)
	{
		if (S.size() == 0)
			break;
		j = 0;
		for (i = 0; i < S.size(); i++)
			if ((text[S[i]] == '*' || text[S[i]] == '/')
				&& text[S[j]] != '*' && text[S[j]] != '/' && text[S[j]] != '^' && text[S[j]] != 'E'
				|| text[S[i]] == '^' && text[S[j]] != '^' && text[S[j]] != 'E'
				|| text[S[i]] == 'E' && text[S[j]] != 'E')
				j = i;
		next = j + 1;
		if (text[S[j]] == '+')
			R[j] += R[next];
		else if (text[S[j]] == '-')
			R[j] -= R[next];
		else if (text[S[j]] == '*')
			R[j] *= R[next];
		else if (text[S[j]] == '/')
			R[j] /= R[next];
		else if (text[S[j]] == '^')
			R[j] ^= R[next];
		else if (text[S[j]] == 'E')
			R[j] *= complex(10).exp(R[next]);
		R.erase(R.begin() + j + 1);
		S.erase(S.begin() + j);
	}
	return R[0];
}

void read(char text[])
{
	vector<complex> computed;
	int i, before, inside, left = 0, right = 0, length = 0;
	while (text[length] != '\0')
		length++;
	while (true)
	{
		before = 0;
		for (i = 0; i < length; i++)
		{
			if (text[i] == '#')
				before++;
			else if (text[i] == '(')
				left = i;
			else if (text[i] == ')')
			{
				right = i;
				break;
			}
		}
		if (left < right)
		{
			text[left] = '#';
			text[right] = ' ';
			char next[10000] = "";
			inside = 0;
			for (i = left + 1; i < right; i++)
			{
				if (text[i] == '#')
					inside++;
				next[i - left - 1] = text[i];
				text[i] = ' ';
			}
			before -= inside;
			computed.push_back(computing(next, computed, before));
			computed.erase(computed.begin() + before, computed.begin() + before + inside);
			left = right;
		}
		else
			break;
	}
	cout << "  = ";
	computing(text, computed, 0).output_endl();
}

void calculation(vector<matrix> X)
{
	/*matrix A = Gauss(X[9], X[8]);
	A.output_close();*/

	matrix A = X[8] * X[7];
	X[8].output();
	X[7].output();
	X[6].output();
	A.output();
}

int main()
{
	cout << "\n";
	cout << "  C  A  L  C  U  L  A  T  O  R                                             ";
	timer();
	cout << "\n\n";
	cout << " ============================================================================================================\n\n";
	matrix A;
	vector<matrix> X;
	while (!file.eof())
	{
		A.input();
		X.emplace(X.begin(), A);
	}
	file.close();
	calculation(X);
	const int length = 10000;
	while (true)
	{
		cout << "  ";
		char text[length] = "";
		cin.getline(text, length);
		if (text[0] == 'e' && text[1] == 'x' && text[2] == 'i' && text[3] == 't' && text[4] == '\0')
			break;
		else
			read(text);
	}
	return 0;
}