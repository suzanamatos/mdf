//#include <QCoreApplication>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <stdio.h>

using namespace Eigen;
using namespace std;

typedef SparseMatrix<double> SpMatrixD;
typedef Triplet<double> T;

enum P
{
    CONSTANT,
    SINSIN,
    POISSON_2
};

double pFunction(P _typeP, double _constant, const Vector2d &_x)
{
    switch (_typeP) {
    case POISSON_2:
    case CONSTANT:
        return _constant;
        break;
    case SINSIN:
        return -sin(_x(0))*sin(_x(1))*_constant;
        break;
    default:
        break;
    }

    return 0;
}

double calculateError(const VectorXd &_calculated, const VectorXd &_analytical)
{
    return (_calculated-_analytical).norm()/_analytical.norm();
}
double calculateError(double _calculated, double _analytical)
{
    return fabs(_calculated-_analytical)/fabs(_analytical);
}

double calculateAnalytical(P _typeP, double _constant, const Vector2d &_x, double L)
{
    switch (_typeP) {
    case CONSTANT:
        return _constant;
        break;
    case SINSIN:
        return 0.5*sin(_x(0))*sin(_x(1));
        break;
//    case MAZZIA:
//        return cos(M_PI);
//        break;
    case POISSON_2:
    {
        //Assumindo domínio 2x2
        //resultado pegue do documento que o prof. Creto mandou
        double value = 0;
        double x = _x(0);
        double y = _x(1);
        double H = L;
        unsigned int max_m = 5, max_n = 5;

        for(unsigned int n = 2; n <= max_n; n++)
        {
            double piN = M_PI*n;
            for(unsigned int m = 2; m <= max_m; m++)
            {
                double piM = M_PI*m;

                //O B JA TEM O H E L = 2
                double B = 4.0*H*H*L*L*_constant*(cos(piM)-1.0)*(cos(piN)-1.0);
                cout << " B " << B << endl;
                double dividend = B*sin((piN*x)/L) * sin((piM*y)/H);
                double divisor =M_PI*M_PI*piM*piN*(H*H*n*n+L*L*m*m);


                value -= dividend/divisor;
            }
        }
        return value;
        break;
    }
    default:
        break;
    }

    return 0;
}

void savePositions(const std::stringstream &_label, unsigned int _n_elem, double _h)
{
    std::stringstream outputPosFileName;
    outputPosFileName.str("");
    outputPosFileName<<"posX" << _label.str();
    cout << outputPosFileName.str()<<endl;
    ofstream posFile(outputPosFileName.str(), std::ofstream::out);

    if (!posFile.is_open())
    {
        cerr<<"[MAIN] o arquivo de resultado não foi aberto"<<endl;
        exit(1);
    }


    for(unsigned int i = 0; i < _n_elem+1; i++)
    {
        posFile << setprecision(15) << _h*i <<endl;
    }

    posFile.close();
}

void saveDisplacements(const std::stringstream &_label, unsigned int _n_elem, const VectorXd &_x)
{
    std::stringstream outputFileName;
    outputFileName.str("");
    outputFileName<<"resultado" << _label.str();
    cout << outputFileName.str()<<endl;
    ofstream resultFile(outputFileName.str(), std::ofstream::out);


    if (!resultFile.is_open())
    {
        cerr<<"[MAIN] o arquivo de resultado não foi aberto"<<endl;
        exit(1);
    }


    //fazendo um pós processamento do resultado para que o
    VectorXd newX = VectorXd::Zero((_n_elem+1)*(_n_elem+1));
    unsigned int count = 0;
    for(unsigned int i = 1; i < _n_elem; i++)
    {
        for(unsigned int j = 1; j < _n_elem; j++)
        {
            newX(i*(_n_elem+1)+j) = _x(count++);
        }
    }

    for(unsigned int i = 0; i < newX.rows(); i++)
    {
        resultFile << setprecision(15) << newX(i) <<endl;
    }

}

int main(int argc, char *argv[])
{

    //entrada
    double L = M_PI; //tamanho do domínio quadrado
    unsigned int n_elem = 10; //numero de elementos em cada direção (SEMPRE PAR)
    double p = 1;   //constante para a função p
    P typeP = P::SINSIN;

    if(argc>1)
        n_elem = atoi(argv[1]);

    //-------------------------------------------------------------------
    //criando as matrizes/vetores
    MatrixXd A = MatrixXd::Constant(n_elem+1, n_elem+1,-1);
    unsigned int size = (n_elem-1)*(n_elem-1);
    VectorXd f = VectorXd::Zero(size);
    VectorXd analytical = VectorXd::Zero(size);
    SpMatrixD K;
    K.resize(size, size);
    K.setZero();
    vector<T> coeffK;

    double h = L/n_elem;


    //-------------------------------------------------------------------
    std::stringstream label;
    label.str("");
    label << "_L-"<<L<<"_N-"<<n_elem<<"_p-"<<p<<"_tipo-"<<int(typeP)<<".txt";
//    savePositions(label, n_elem, h);


    //-------------------------------------------------------------------
    //preenchendo a matriz A que contem o nº das equacoes
    unsigned int count = 0;
    for(unsigned int i = 1; i < n_elem; i++)
    {
        for(unsigned int j = 1; j < n_elem; j++)
        {
            A(i,j) = count++;
        }
    }

    //-------------------------------------------------------------------
    //preenchendo o sistema
    for(unsigned int i = 1; i < n_elem; i++)
    {
        for(unsigned int j = 1; j < n_elem; j++)
        {
            if(A(i-1, j) >= 0)
                coeffK.push_back(T(A(i, j), A(i-1, j), 1));

            if(A(i, j-1) >= 0)
                coeffK.push_back(T(A(i, j), A(i, j-1), 1));

            if(A(i+1, j) >= 0)
                coeffK.push_back(T(A(i, j), A(i+1, j), 1));

            if(A(i, j+1) >= 0)
                coeffK.push_back(T(A(i, j), A(i, j+1), 1));

            coeffK.push_back(T(A(i, j), A(i, j), -4));

            Vector2d pos(h*i, h*j);
            f(A(i, j)) = h*h * pFunction(typeP, p, pos);
            analytical(A(i, j)) = calculateAnalytical(typeP, p, pos, L);
        }
    }

    K.setFromTriplets(coeffK.begin(), coeffK.end());
    coeffK.clear();

    //MatrixXd K0=MatrixXd(K);


    //-------------------------------------------------------------------
    //resolvendo o sistema
    SparseLU<SpMatrixD> lu;
    K.makeCompressed();
    lu.analyzePattern(K);
    lu.factorize(K);

    if(lu.info()!=0)
    {
        std::cout<<"\t\t\t\tnão conseguiu fatorar a matriz de !!!"<<std::endl;
        exit(1);
    }
    VectorXd x = lu.solve(f);


    //testar esses valores: https://www.ufsj.edu.br/portal2-repositorio/File/nepomuceno/mn/21MN_EDO4.pdf
    //esse aqui tmb: https://phkonzen.github.io/notas/MatematicaNumerica/cap_edp_sec_Poisson.html
    //testados e passados!


    cout <<"acabou processo"<< endl;

    saveDisplacements(label, n_elem, x);

    //-------------------------------------------------------------------
    std::stringstream outputFileName;
    outputFileName.str("");
    outputFileName<<"info" << label.str();
    ofstream resultInfoFile(outputFileName.str(), std::ofstream::out);

    if (!resultInfoFile.is_open())
    {
        cerr<<"[MAIN] o arquivo de info não foi aberto"<<endl;
        exit(1);
    }

    cout <<"Adicionando no arquivo"<< endl;

    /*resultInfoFile << setprecision(15) << "Residual: "
                   << (K0 * x - f).norm() << endl;*/
    resultInfoFile << setprecision(15) << "Erro relativo: "
                   << calculateError(x, analytical) << endl;

    unsigned int middle = (int)ceil(size/2);

    resultInfoFile << "\n\nponto central: "
                   << h*ceil(n_elem/2) << endl;
    resultInfoFile << "ponto central calculado: "
                   << x(middle) << endl;
    resultInfoFile << "ponto central analitico: "
                   << analytical(middle) << endl;
    resultInfoFile << "erro do ponto central: "
                   << calculateError(x(middle), analytical(middle)) << endl;

    cout <<"Fim"<< endl;

    return 0;
}
