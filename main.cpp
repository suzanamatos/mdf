//#include <QCoreApplication>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <limits>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include "timeSpent.h"

using namespace Eigen;
using namespace std;

typedef SparseMatrix<double> SpMatrixD;
typedef Triplet<double> T;

enum P
{
    CONSTANT,
    SINSIN,
    POISSON_2,
    POISSON_EXTRA
};

double pFunction(P _typeP, double _constant, const Vector2d &_x)
{
    switch (_typeP) {
    case POISSON_EXTRA:
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
    case POISSON_EXTRA:
    case CONSTANT:
        cout << "Não implementado resultado analítico" <<endl;
        exit(1);
        return _constant; //ta errado!!!???
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
    outputPosFileName<<"posXY" << _label.str();
    cout << outputPosFileName.str()<<endl;
    ofstream posFile(outputPosFileName.str(), std::ofstream::out);

    if (!posFile.is_open())
    {
        cerr<<"[MAIN] o arquivo de resultado não foi aberto"<<endl;
        exit(1);
    }


    for(unsigned int i = 0; i < _n_elem+1; i++)
    {
        for(unsigned int j = 0; j < _n_elem+1; j++)
//        {
            posFile << setprecision(15) << _h*i << "\t"<< _h*j <<endl;
//            posFile << setprecision(15) << _h*j <<endl;
//        }
    }

    posFile.close();
}

void savePositions(const std::stringstream &_label,
                   double _L, double _D, unsigned int _n_elemL, unsigned int _n_elemD,
                   unsigned int _n_elemL_part, unsigned int _n_elemD_part,
                   double _hx, double _hy)
{
    std::stringstream outputPosFileName;
    outputPosFileName.str("");
    outputPosFileName<<"extra_posXY" << _label.str();
    cout << outputPosFileName.str()<<endl;
    ofstream posFile(outputPosFileName.str(), std::ofstream::out);

    if (!posFile.is_open())
    {
        cerr<<"[MAIN] o arquivo de resultado não foi aberto"<<endl;
        exit(1);
    }
//    MatrixXd *A = new MatrixXd(_n_elemD+1, _n_elemL+1);
//    A->setConstant(-1);

//    unsigned axis = 0;
    for(unsigned int i = 0; i < _n_elemD_part; i++)
    {
        for(unsigned int j = _n_elemL_part/2; j <= _n_elemL; j++)
        {
            Vector2d position = Vector2d(j*_hx, _D - (i*_hy));

            posFile << setprecision(15) << position(0) << "\t"<< position(1) <<endl;

//            (*A)(i,j) = position(axis);
        }
    }
    for(unsigned int j = 0; j <= _n_elemL; j++)
    {
        Vector2d position = Vector2d(j*_hx, 0);

        posFile << setprecision(15) << position(0) << "\t"<< position(1) <<endl;

//        (*A)(_n_elemD_part, j) = position(axis);
    }
    for(unsigned int i = _n_elemD_part+1; i <= _n_elemD; i++)
    {
        for(unsigned int j = 0; j <= _n_elemL_part; j++)
        {
            Vector2d position = Vector2d(j*_hx, _D - (i*_hy));

            posFile << setprecision(15) << position(0) << "\t"<< position(1) <<endl;

//            (*A)(i,j) = position(axis);
        }
    }
//    cout << "Matrix positions X\n" << *A << endl;
//    delete A;

//    for(unsigned int i = 0; i < _n_elemL+1; i++)
//    {
//        Vector2d position = Vector2d(i * _hx, 0.0);

//        posFile << setprecision(15) << position(0) << "\t"<< position(1) <<endl;
//    }

//    for(unsigned int i = 0; i < _n_elemL_part+1; i++)
//    {
//        for(unsigned int j = 1; j < _n_elemD_part+1; j++)
//        {

//            Vector2d position = Vector2d(i * _hx, -(j * _hy));
//            posFile << setprecision(15) << position(0) << "\t"<< position(1) <<endl;


//            position = Vector2d((3.0*_L/2.0) - (i * _hx), j * _hy);
//            posFile << setprecision(15) << position(0) << "\t"<< position(1) <<endl;
//        }
//    }

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

//    cout << "valor calculado\n"<<_x << endl;


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

//void saveDisplacements(const std::stringstream &_label,
//                       double _L, double _D, unsigned int _n_elemL, unsigned int _n_elemD,
//                       unsigned int _n_elemL_part, unsigned int _n_elemD_part,
//                       const VectorXd &_x)
//{
//    std::stringstream outputFileName;
//    outputFileName.str("");
//    outputFileName<<"extra_resultado" << _label.str();
//    cout << outputFileName.str()<<endl;
//    ofstream resultFile(outputFileName.str(), std::ofstream::out);


//    if (!resultFile.is_open())
//    {
//        cerr<<"[MAIN] o arquivo de resultado não foi aberto"<<endl;
//        exit(1);
//    }


//    //fazendo um pós processamento do resultado para que o
////    cout << "size " << (2*(_n_elemL_part+1)*((_n_elemD_part+1)) - ceil((_n_elemL_part+1)/2.)) << endl;
//    VectorXd newX = VectorXd::Zero(2*(_n_elemL_part+1)*((_n_elemD_part+1)) - ceil((_n_elemL_part+1)/2.));
//    unsigned int count = 0;
//    unsigned int i = 0;
//    for(i = 1; i < _n_elemD_part; i++)
//    {
//        for(unsigned int j = 1; j < _n_elemL_part; j++)
//        {
////            cout << "1 POSICAO " << i*(_n_elemL_part+1)+j <<  endl;
//            newX(i*(_n_elemL_part+1)+j) = _x(count++);
//        }
//    }
//    for(unsigned int j = _n_elemL_part/2 + 1; j < _n_elemL_part; j++)
//    {
////        cout << "2 POSICAO " << i*(_n_elemL_part+1)+j <<  endl;
//        newX(i*(_n_elemL_part+1)+j) = _x(count++);
//    }
//    for(i = _n_elemD_part +1; i < _n_elemD; i++)
//    {
//        for(unsigned int j = 1; j < _n_elemL_part; j++)
//        {
////            cout << "3 POSICAO " << i*(_n_elemL_part+1)+j + (_n_elemL_part)/2.  <<  endl;
//            newX(i*(_n_elemL_part+1)+j + (_n_elemL_part)/2.) = _x(count++);
//        }
//    }

//    for(unsigned int i = 0; i < newX.rows(); i++)
//    {
//        resultFile << setprecision(15) << newX(i) <<endl;
//    }

//}

void saveDisplacements(const std::stringstream &_label, double _L, double _D,
                       double _hx, double _hy, unsigned int _n_elemL, unsigned int _n_elemD,
                       unsigned int _n_elemL_part, unsigned int _n_elemD_part,
                       const VectorXd &_x)
{
    std::stringstream outputFileName;
    outputFileName.str("");
    outputFileName<<"extra_resultado" << _label.str();
    cout << outputFileName.str()<<endl;
    ofstream resultFile(outputFileName.str(), std::ofstream::out);


    if (!resultFile.is_open())
    {
        cerr<<"[MAIN] o arquivo de resultado não foi aberto"<<endl;
        exit(1);
    }


    MatrixXd *A = new MatrixXd(_n_elemD+1, _n_elemL+1);
    A->setConstant(std::numeric_limits<double>::infinity());

    unsigned int count = 0;
    for(unsigned int i = 0; i < _n_elemD_part; i++)
    {
        for(unsigned int j = _n_elemL_part/2; j <= _n_elemL; j++)
        {
            if ( (i == 0) || (j == _n_elemL_part/2) || (j == _n_elemL))
                (*A)(i,j) = 0;
            else
                (*A)(i,j) = _x[count++];
        }
    }
    cout << "_n_elemL_part/2 " << _n_elemL_part/2 << endl;
    for(unsigned int j = 0; j <= _n_elemL; j++)
    {
        if ((j <= _n_elemL_part/2) || (j >=_n_elemL_part))
            (*A)(_n_elemD_part,j) = 0;
        else
            (*A)(_n_elemD_part, j) = _x[count++];
    }
    for(unsigned int i = _n_elemD_part+1; i <= _n_elemD; i++)
    {
        for(unsigned int j = 0; j <= _n_elemL_part; j++)
        {
            if ( (i == _n_elemD) || (j == 0) || (j == _n_elemL_part))
                (*A)(i,j) = 0;
            else
                (*A)(i,j) = _x[count++];
        }
    }

    cout << "resultado" << endl << *A << endl;


    resultFile << setprecision(15) << _L << "\t" << _D << "\t" << A->rows() << "\t" << A->cols() << "\t" << _hx << "\t" << _hy  <<endl;
    for(unsigned int i = 0; i < A->rows(); i++)
    {
        for(unsigned int j = 0; j < A->cols(); j++)
        {
            resultFile << setprecision(15) << (*A)(i,j) << "\t";
        }
        resultFile <<endl;
    }




    delete A;
}

void quadModel(unsigned int n_elem, double L, double p, P typeP)
{
    //-------------------------------------------------------------------
    //criando as matrizes/vetores
    MatrixXd *A = new MatrixXd(n_elem+1, n_elem+1);
    A->setConstant(-1);
    unsigned int size = (n_elem-1)*(n_elem-1);
    VectorXd f = VectorXd::Zero(size);
    SpMatrixD K;
    K.resize(size, size);
    K.reserve(size*5);
    VectorXd analytical = VectorXd::Zero(size);
    vector<T> coeffK;
    coeffK.reserve(size*5);

    cout << "total memory (resrvado) " << size*5*sizeof(T)/(1e6) << endl;

    double h = L/n_elem;


    //-------------------------------------------------------------------
    std::stringstream label;
    label.str("");
    label << "_L-"<<L<<"_N-"<<n_elem<<"_p-"<<p<<"_tipo-"<<int(typeP)<<"COEFFK_Cholesky.txt";
    savePositions(label, n_elem, h);

    //exit(1);


    cout << "Criando matriz auxiliar" << endl;
    TimeSpent criandoATime("Criada matriz auxiliar | TEMPO: ");
    criandoATime.startCount();
    //-------------------------------------------------------------------
    //preenchendo a matriz A que contem o nº das equacoes
    unsigned int count = 0;
    for(unsigned int i = 1; i < n_elem; i++)
    {
        for(unsigned int j = 1; j < n_elem; j++)
        {
            (*A)(i,j) = count++;
        }
    }
    criandoATime.endCount();
    criandoATime.print();


    cout << "Preenchendo o sistema" << endl;
    TimeSpent preenchendoSistemaTime("Preenchido o sistema | TEMPO: ");
    preenchendoSistemaTime.startCount();
    //-------------------------------------------------------------------
    //preenchendo o sistema
//    Eigen::initParallel();
//    Eigen::setNbThreads(2);
//    #pragma omp parallel for num_threads(2)
    for(unsigned int j = 1; j < n_elem; j++)
    {
        for(unsigned int i = 1; i < n_elem; i++)
        {
            if((*A)(i-1, j) >= 0)
                coeffK.push_back(T((*A)(i, j), (*A)(i-1, j), 1));

            if((*A)(i, j-1) >= 0)
                coeffK.push_back(T((*A)(i, j), (*A)(i, j-1), 1));

            if((*A)(i+1, j) >= 0)
                coeffK.push_back(T((*A)(i, j), (*A)(i+1, j), 1));

            if((*A)(i, j+1) >= 0)
                coeffK.push_back(T((*A)(i, j), (*A)(i, j+1), 1));

            coeffK.push_back(T((*A)(i, j), (*A)(i, j), -4));

            Vector2d pos(h*i, h*j);

            f.coeffRef((*A)(i, j), 0) = h*h * pFunction(typeP, p, pos);
            analytical((*A)(i, j)) = calculateAnalytical(typeP, p, pos, L);
        }
    }

    cout << "total memory (real) " << coeffK.size()*sizeof(T)/(1e6) << endl;
    delete A;
    K.setFromTriplets(coeffK.begin(), coeffK.end());
    coeffK.clear();
    K.makeCompressed();
    preenchendoSistemaTime.endCount();
    preenchendoSistemaTime.print();

//    cout << "\n\nK"<<endl;
//    cout << MatrixXd(K)<<endl;

//    cout << "\n\nf"<<endl;
//    cout << MatrixXd(f)<<endl;


    cout << "Resolvendo o sistema" << endl;
    TimeSpent resolvendoSistemaTime("Resolvido o sistema | TEMPO: ");
    resolvendoSistemaTime.startCount();
    //-------------------------------------------------------------------
    //resolvendo o sistema
    SimplicialCholesky<SpMatrixD> chol(K);
    VectorXd x = chol.solve(f);
    resolvendoSistemaTime.endCount();
    resolvendoSistemaTime.print();
//    cout <<"\n\n"<< MatrixXd(x)<<endl;


    //testar esses valores: https://www.ufsj.edu.br/portal2-repositorio/File/nepomuceno/mn/21MN_EDO4.pdf
    //esse aqui tmb: https://phkonzen.github.io/notas/MatematicaNumerica/cap_edp_sec_Poisson.html
    //testados e passados!


    cout <<"acabou processo"<< endl;

    cout << "Salvando o resultado" << endl;
    TimeSpent salvandoTime("Salvo o resultado | TEMPO: ");
    salvandoTime.startCount();
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

    salvandoTime.endCount();
    salvandoTime.print();
    cout <<"Fim"<< endl;
}

void ladderModel(unsigned int _n_elemL_part, unsigned int _n_elemD_part, double L, double D, double p, P typeP)
{
    unsigned int n_elemL = _n_elemL_part + _n_elemL_part/2;
    unsigned int n_elemD = 2*_n_elemD_part;


    //-------------------------------------------------------------------
    //criando as matrizes/vetores
    MatrixXd *A = new MatrixXd(n_elemD+1, n_elemL+1);
    A->setConstant(-1);
    unsigned int size = 2*(_n_elemL_part-1)*(_n_elemD_part-1) + _n_elemL_part/2 - 1;
    VectorXd f = VectorXd::Zero(size);
    SpMatrixD K;
    K.resize(size, size);
    VectorXd analytical = VectorXd::Zero(size);
    vector<T> coeffK;


    double hx = L/_n_elemL_part;
    double hy = D/_n_elemD_part;

    cout << "_n_elemL_part " << _n_elemL_part << endl;
    cout << "_n_elemD_part " << _n_elemD_part << endl;
    cout << "n_elemL " << n_elemL << endl;
    cout << "n_elemD " << n_elemD << endl;
    cout << "hx " << hx << endl;
    cout << "hy " << hy << endl;
    cout << "tamanho do sistema " << size << endl;


    //-------------------------------------------------------------------
    std::stringstream label;
    label.str("");
    label << "_L-"<<L<< "_D-"<<D<<"_Nl-"<<_n_elemL_part<<"_Nd-"<<_n_elemD_part<<"_p-"<<p<<"_tipo-"<<int(typeP)<<"COEFFK_Cholesky.txt";

//    savePositions(label, L, D, n_elemL, n_elemD, _n_elemL_part, _n_elemD_part, hx, hy);


    cout << "Criando matriz auxiliar" << endl;
    TimeSpent criandoATime("Criada matriz auxiliar | TEMPO: ");
    criandoATime.startCount();
    //-------------------------------------------------------------------
    //preenchendo a matriz A que contem o nº das equacoes
    //notar que é diferente da posição, uma vez que que a matriz
    //linha -- Y
    //coluna -- X
    unsigned int count = 0;
    for(unsigned int i = 1; i < _n_elemD_part; i++)
    {
        for(unsigned int j = _n_elemL_part/2 + 1; j < n_elemL; j++)
        {
            (*A)(i,j) = count++;
        }
    }
    for(unsigned int j = _n_elemL_part/2 + 1; j < _n_elemL_part; j++)
    {
        (*A)(_n_elemD_part, j) = count++;
    }
    for(unsigned int i = _n_elemD_part +1; i < n_elemD; i++)
    {
        for(unsigned int j = 1; j < _n_elemL_part; j++)
        {
            (*A)(i,j) = count++;
        }
    }

    cout << "A" << endl << *A << endl;
    criandoATime.endCount();
    criandoATime.print();


    cout << "Preenchendo o sistema" << endl;
    TimeSpent preenchendoSistemaTime("Preenchido o sistema | TEMPO: ");
    preenchendoSistemaTime.startCount();
    //-------------------------------------------------------------------
    //preenchendo o sistema
//    Eigen::initParallel();
//    Eigen::setNbThreads(2);
//    #pragma omp parallel for num_threads(2)
    for(unsigned int j = 1; j < n_elemL; j++)
    {
        for(unsigned int i = 1; i < n_elemD; i++)
        {
            if((*A)(i, j) >= 0)
            {
                coeffK.push_back(T((*A)(i, j), (*A)(i, j), -4));
//                cout << "Linha e coluna da matriz K  " <<(*A)(i, j)<< endl;
            }
            else continue;

            if((*A)(i-1, j) >= 0)
            {
                coeffK.push_back(T((*A)(i, j), (*A)(i-1, j), 1));
//                cout << "Linha da matriz K  " <<(*A)(i, j)<< endl;
//                cout << "Coluna da matriz K  " <<(*A)(i-1, j)<< endl;
            }

            if((*A)(i, j-1) >= 0)
            {
                coeffK.push_back(T((*A)(i, j), (*A)(i, j-1), 1));
//                cout << "Linha da matriz K  " <<(*A)(i, j)<< endl;
//                cout << "Coluna da matriz K  " <<(*A)(i, j-1)<< endl;
            }

            if((*A)(i+1, j) >= 0)
            {
                coeffK.push_back(T((*A)(i, j), (*A)(i+1, j), 1));
//                cout << "Linha da matriz K  " <<(*A)(i, j)<< endl;
//                cout << "Coluna da matriz K  " <<(*A)(i+1, j)<< endl;
            }

            if((*A)(i, j+1) >= 0)
            {
                coeffK.push_back(T((*A)(i, j), (*A)(i, j+1), 1));
//                cout << "Linha da matriz K  " <<(*A)(i, j)<< endl;
//                cout << "Coluna da matriz K  " <<(*A)(i, j+1)<< endl;
            }


            Vector2d pos(hx*j, D - (hy*i));
//            cout << "posição" <<pos<< endl<<endl<<endl;

            f.coeffRef((*A)(i, j), 0) = hx*hy * pFunction(typeP, p, pos);
        }
    }

    delete A;
    K.setFromTriplets(coeffK.begin(), coeffK.end());
    coeffK.clear();
    K.makeCompressed();
    preenchendoSistemaTime.endCount();
    preenchendoSistemaTime.print();

//    cout << "\n\nK " <<endl;
//    cout << MatrixXd(K)<<endl;

//    cout << "\n\nf"<<endl;
//    cout << MatrixXd(f)<<endl;


    cout << "Resolvendo o sistema" << endl;
    TimeSpent resolvendoSistemaTime("Resolvido o sistema | TEMPO: ");
    resolvendoSistemaTime.startCount();
    //-------------------------------------------------------------------
    //resolvendo o sistema
    SimplicialCholesky<SpMatrixD> chol(K);
    VectorXd x = chol.solve(f);
    resolvendoSistemaTime.endCount();
    resolvendoSistemaTime.print();
    cout <<"\n\n"<< MatrixXd(x)<<endl;



    cout <<"acabou processo"<< endl;

    cout << "Salvando o resultado" << endl;
    TimeSpent salvandoTime("Salvo o resultado | TEMPO: ");
    salvandoTime.startCount();
    saveDisplacements(label, L, D, hx, hy, n_elemL, n_elemD, _n_elemL_part, _n_elemD_part, x);

    cout <<"Fim"<< endl;
}

int main(int argc, char *argv[])
{
//    unsigned int n_elem =10; //numero de elementos em cada direção (SEMPRE PAR)
//    double L = M_PI; //tamanho do domínio quadrado
//    double p = 1;   //constante para a função p
//    P typeP = P::SINSIN;
//    if(argc>1)
//        n_elem = atoi(argv[1]);
//    quadModel(n_elem, L, p, typeP);

    unsigned int n_elemL = 6; //numero de elementos em cada direção X
    unsigned int n_elemD = 6; //numero de elementos em cada direção Y
    double L = 2; //tamanho do domínio quadrado
    double D = 2; //tamanho do domínio quadrado
    double p = 1;   //constante para a função p
    P typeP = P::POISSON_EXTRA;
    if(argc>1)
        n_elemL = n_elemD= atoi(argv[1]);
    ladderModel(n_elemL, n_elemD, L, D, p, typeP);

    return 0;
}
