#include <fstream>
#include <string>
#include <memory>
#include <vector>

std::string s;
int cur;

class EWrongInputData {
};

enum ConnectionClass {
	CONNECTION_CLASS_EDGE,
	CONNECTION_CLASS_SEQUENTIAL,
	CONNECTION_CLASS_PARALLEL
};

struct Connection {
	virtual ~Connection() {}
	virtual ConnectionClass GetConnectionClass() const = 0;
};

typedef std::auto_ptr<Connection> ConnectionP;

struct EdgeConnection : public Connection {
	double p;
	EdgeConnection(double p1)
		: p(p1)
	{
	}
	ConnectionClass GetConnectionClass() const
	{
		return CONNECTION_CLASS_EDGE;
	}
};

struct SequentialConnection : public Connection {
	ConnectionP c1, c2;
	SequentialConnection(ConnectionP _c1, ConnectionP _c2)
		: c1(_c1), c2(_c2)
	{
	}
	ConnectionClass GetConnectionClass() const
	{
		return CONNECTION_CLASS_SEQUENTIAL;
	}
};

struct ParallelConnection : public Connection {
	ConnectionP c1, c2;
	ParallelConnection(ConnectionP _c1, ConnectionP _c2)
		: c1(_c1), c2(_c2)
	{
	}
	ConnectionClass GetConnectionClass() const
	{
		return CONNECTION_CLASS_PARALLEL;
	}
};

class EPolynomDivideError {
};

struct Polynom {
	std::vector<double> coeff;
	Polynom operator*(Polynom p)
	{
		Polynom res;
		res.coeff.resize(coeff.size()+p.coeff.size()-1);
		for (int i = 0; i < coeff.size(); ++i) {
			for (int j = 0; j < p.coeff.size(); ++j) {
				res.coeff[i+j] += coeff[i]*p.coeff[j];
			}
		}
		return res;
	}
	Polynom operator*(double a)
	{
		Polynom res;
		for (int i = 0; i < coeff.size(); ++i)
			res.coeff.push_back(coeff[i]*a);
		return res;
	}
	Polynom operator+(Polynom p)
	{
		Polynom res;
		if (coeff.size() > p.coeff.size())
			res.coeff.resize(coeff.size());
		else
			res.coeff.resize(p.coeff.size());
		for (int i = 0; i < res.coeff.size(); ++i) {
			if (i < coeff.size())
				res.coeff[i] += coeff[i];
			if (i < p.coeff.size())
				res.coeff[i] += p.coeff[i];
		}
		return res;
	}
	void Divide(int k)
	{
		for (int i = 0; i < k; ++i) {
			if (coeff[i] != 0)
				throw EPolynomDivideError();
		}
		int n = coeff.size();
		for (int i = 0; i < n-k; ++i)
			coeff[i] = coeff[i+k];
		coeff.resize(n-k);
	}
};

void miss_spaces()
{
	while (cur < s.length() && s[cur] == ' ')
		++cur;
}

bool is_digit(char ch)
{
	return ch >= '0' && ch <= '9';
}

ConnectionP ParseConnection()
{
	miss_spaces();
	if (cur == s.length())
		throw EWrongInputData();
	if (s[cur] == '(') {
		++cur;
		ConnectionP c1 = ParseConnection();
		while (1) {
			ConnectionP c2 = ParseConnection();
			c1 = ConnectionP(new SequentialConnection(c1, c2));
			miss_spaces();
			if (cur < s.length() && s[cur] == ')') {
				++cur;
				break;
			}
		}
		return c1;
	}
	else if (s[cur] == '[') {
		++cur;
		ConnectionP c1 = ParseConnection();
		while (1) {
			ConnectionP c2 = ParseConnection();
			c1 = ConnectionP(new ParallelConnection(c1, c2));
			miss_spaces();
			if (cur < s.length() && s[cur] == ']') {
				++cur;
				break;
			}
		}
		return c1;
	}
	else if (is_digit(s[cur])) {
		double num = 0.0;
		while (cur < s.length() && is_digit(s[cur])) {
			int digit = (int)(s[cur]-'0');
			num = num*10 + digit;
			++cur;
		}
		if (cur < s.length() && s[cur] == '.') {
			++cur;
			if (!(cur < s.length() && is_digit(s[cur])))
				throw EWrongInputData();
			double z = 1.0;
			while (cur < s.length() && is_digit(s[cur])) {
				int digit = (int)(s[cur]-'0');
				z /= 10;
				num += z*digit;
				++cur;
			}
		}
		if (num > 1)
			throw EWrongInputData();
		return ConnectionP(new EdgeConnection(num));
	}
	else
		throw EWrongInputData();
}

void probability_of_connectivity(ConnectionP &c, double &A, double &B)
{
	double A1, B1, A2, B2;
	switch (c->GetConnectionClass()) {
		case CONNECTION_CLASS_EDGE:
			A = static_cast<EdgeConnection*>(c.get())->p;
			B = 1-A;
			break;
		case CONNECTION_CLASS_SEQUENTIAL:
			probability_of_connectivity(static_cast<SequentialConnection*>(c.get())->c1, A1, B1);
			probability_of_connectivity(static_cast<SequentialConnection*>(c.get())->c2, A2, B2);
			A = A1*A2;
			B = A1*B2 + B1*A2;
			break;
		case CONNECTION_CLASS_PARALLEL:
			probability_of_connectivity(static_cast<ParallelConnection*>(c.get())->c1, A1, B1);
			probability_of_connectivity(static_cast<ParallelConnection*>(c.get())->c2, A2, B2);
			A = A1*A2 + B1*A2 + A1*B2;
			B = B1*B2;
			break;
	}
}

int conj(int a, int b)
{
	return a*b;
}

int disj(int a, int b)
{
	if (a == 0 && b == 0)
		return 0;
	else
		return 1;
}

void moments(ConnectionP &c, double (&m1)[2], double (&m2)[2], double (&P)[2])
{
	double m1_1[2], m1_2[2];
	double m2_1[2], m2_2[2];
	double P1[2], P2[2];
	switch (c->GetConnectionClass()) {
		case CONNECTION_CLASS_EDGE:
			P[1] = static_cast<EdgeConnection*>(c.get())->p;
			P[0] = 1-P[1];
			m1[0] = 2;
			m1[1] = 1;
			m2[0] = 4;
			m2[1] = 1;
			break;
		case CONNECTION_CLASS_SEQUENTIAL:
			moments(static_cast<SequentialConnection*>(c.get())->c1, m1_1, m2_1, P1);
			moments(static_cast<SequentialConnection*>(c.get())->c2, m1_2, m2_2, P2);
			P[1] = P1[1]*P2[1];
			P[0] = 1-P[1];
			for (int a = 0; a <= 1; ++a) {
				m1[a] = 0.0;
				for (int a1 = 0; a1 <= 1; ++a1) {
					for (int a2 = 0; a2 <= 1; ++a2) {
						if (conj(a1, a2) == a)
							m1[a] += P1[a1]*P2[a2]/P[a]*(m1_1[a1]+m1_2[a2]-1);							
					}
				}
				m2[a] = 0.0;
				for (int a1 = 0; a1 <= 1; ++a1) {
					for (int a2 = 0; a2 <= 1; ++a2) {
						if (conj(a1, a2) == a)
							m2[a] += P1[a1]*P2[a2]/P[a]*(m2_1[a1]+m2_2[a2]+1+2*m1_1[a1]*m1_2[a2]-2*m1_1[a1]-2*m1_2[a2]);
					}
				}
			}
			break;
		case CONNECTION_CLASS_PARALLEL:
			moments(static_cast<ParallelConnection*>(c.get())->c1, m1_1, m2_1, P1);
			moments(static_cast<ParallelConnection*>(c.get())->c2, m1_2, m2_2, P2);
			P[0] = P1[0]*P2[0];
			P[1] = 1-P[0];
			for (int a = 0; a <= 1; ++a) {
				m1[a] = 0.0;
				for (int a1 = 0; a1 <= 1; ++a1) {
					for (int a2 = 0; a2 <= 1; ++a2) {
						if (disj(a1, a2) == a)
							m1[a] += P1[a1]*P2[a2]/P[a]*(m1_1[a1]+m1_2[a2]-2+conj(a1, a2));
					}
				}
				m2[a] = 0.0;
				for (int a1 = 0; a1 <= 1; ++a1) {
					for (int a2 = 0; a2 <= 1; ++a2) {
						if (disj(a1, a2) == a)
							m2[a] += P1[a1]*P2[a2]/P[a]*(m2_1[a1]+m2_2[a2]+4+2*m1_1[a1]*m1_2[a2]-4*m1_1[a1]-4*m1_2[a2]+(2*m1_1[a1]+2*m1_2[a2]-3)*conj(a1, a2));
					}
				}
			}
			break;
	}
}

void generating_function(ConnectionP &c, Polynom (&gf)[2], double (&P)[2])
{
	Polynom gf1[2], gf2[2];
	double P1[2], P2[2];
	switch (c->GetConnectionClass()) {
		case CONNECTION_CLASS_EDGE:
			P[1] = static_cast<EdgeConnection*>(c.get())->p;
			P[0] = 1-P[1];
			gf[1].coeff.resize(2);
			gf[1].coeff[1] = 1.0;
			gf[0].coeff.resize(3);
			gf[0].coeff[2] = 1.0;
			break;
		case CONNECTION_CLASS_SEQUENTIAL:
			generating_function(static_cast<SequentialConnection*>(c.get())->c1, gf1, P1);
			generating_function(static_cast<SequentialConnection*>(c.get())->c2, gf2, P2);
			P[1] = P1[1]*P2[1];
			P[0] = 1-P[1];
			for (int a = 0; a <= 1; ++a) {
				gf[a].coeff.resize(0);
				for (int a1 = 0; a1 <= 1; ++a1) {
					for (int a2 = 0; a2 <= 1; ++a2) {
						if (conj(a1, a2) == a) {
							Polynom p = gf1[a1]*gf2[a2];
							p.Divide(1);
							gf[a] = gf[a] + p*(P1[a1]*P2[a2]/P[a]);
						}
					}
				}
			}
			break;
		case CONNECTION_CLASS_PARALLEL:
			generating_function(static_cast<ParallelConnection*>(c.get())->c1, gf1, P1);
			generating_function(static_cast<ParallelConnection*>(c.get())->c2, gf2, P2);
			P[0] = P1[0]*P2[0];
			P[1] = 1-P[0];
			for (int a = 0; a <= 1; ++a) {
				gf[a].coeff.resize(0);
				for (int a1 = 0; a1 <= 1; ++a1) {
					for (int a2 = 0; a2 <= 1; ++a2) {
						if (disj(a1, a2) == a) {
							Polynom p = gf1[a1]*gf2[a2];
							p.Divide(2-conj(a1, a2));
							gf[a] = gf[a] + p*(P1[a1]*P2[a2]/P[a]);
						}
					}
				}
			}
			break;
	}
}

int main()
{
	std::ifstream fin("input.txt");
	std::ofstream fout("output.txt");
	getline(fin, s);

	cur = 0;
	try {
		ConnectionP c = ParseConnection();
		miss_spaces();
		if (cur < s.length())
			throw EWrongInputData();

		double A, tmp;
		probability_of_connectivity(c, A, tmp);
		double m1[2], m2[2], P[2];
		moments(c, m1, m2, P);
		double M = P[0]*m1[0] + P[1]*m1[1];
		double X = P[0]*m2[0] + P[1]*m2[1];
		double D = X - M*M;
		Polynom gf[2];
		generating_function(c, gf, P);
		Polynom GF = gf[0]*P[0] + gf[1]*P[1];

		fout << "Вероятность связности = " << A << "\n";
		fout << "Мат. ожидание числа компонент связности = " << M << "\n";
		fout << "Дисперсия числа компонент связности = " << D << "\n";
		fout << "Производящая функция числа компонент связности = ";
		for (int i = 0; i < GF.coeff.size(); ++i) {
			if (i > 0)
				fout << " + ";
			fout << GF.coeff[i] << "z^" << i;
		}
	}
	catch(EWrongInputData) {
		fout << "Ошибка во входном файле, символ " << cur+1;
	}

	return 0;
}