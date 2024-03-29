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
    PATCH_TEST,
    SINSIN,
    POISSON_2,
    POISSON_EXTRA
};


enum BORDER_TYPE
{
    PROG_X_REGR_Y,
    REGR_X_REGR_Y,
    CENT_X_REGR_Y,
    PROG_X_CENT_Y,
    REGR_X_CENT_Y,
    CENT_X_CENT_Y,
    PROG_X_PROG_Y,
    REGR_X_PROG_Y,
    CENT_X_PROG_Y
};

struct XY
{
    XY() : x (-1), y(-1), borderType(-1){}
    XY(int _x, int _y, unsigned int _borderType) : x (_x),y(_y), borderType(_borderType){}
    int x;
    int y;
    int borderType;
};

typedef Matrix<XY, Dynamic, Dynamic> NumEquationMatrix;

double pFunction(P _typeP, double _constant, const Vector2d &_x = Vector2d(0,0))
{
    switch (_typeP) {
    case PATCH_TEST:
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



double calculateConditionNumber(const MatrixXd &A)
{
//    cout << "\n-------------------------\n\n----------------------\n" <<  setprecision(15) << A <<endl;


    VectorXcd eivals = A.eigenvalues();
    cout << "autovalores \n" << eivals <<endl;

    std::vector<double> myvector;
    for(auto i = 0; i < eivals.rows(); i++)
        myvector.push_back( fabs(eivals[i].real()));

    std::sort (myvector.begin(), myvector.end());

    double conditionNumber = fabs(myvector.back())/fabs(myvector.front());

    cout << "Condition Number: " << conditionNumber << endl;

    return conditionNumber;
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
    case PATCH_TEST:
        return _x(0)+_x(1);
        break;
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


#define VEUSZ 1
void saveDisplacements(const std::stringstream &_label, double _L, double _D,
                       double _hx, double _hy, unsigned int _n_elemL, unsigned int _n_elemD,
                       unsigned int _n_elemL_part, unsigned int _n_elemD_part,
                       const VectorXd &_x)
{
    std::stringstream outputFileName;
    outputFileName.str("");
    #if VEUSZ
        outputFileName<<"WITHOUTTOPextra_resultado" << _label.str();
    #else
        outputFileName<<"extra_resultado" << _label.str();
    #endif
    cout << outputFileName.str()<<endl;
    ofstream resultFile(outputFileName.str(), std::ofstream::out);


    if (!resultFile.is_open())
    {
        cerr<<"[MAIN] o arquivo de resultado não foi aberto"<<endl;
        exit(1);
    }


    MatrixXd *A = new MatrixXd(_n_elemD+1, _n_elemL+1);

    A->setConstant(0);

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
//    cout << "_n_elemL_part/2 " << _n_elemL_part/2 << endl;
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

//    cout << "resultado" << endl << *A << endl;


    #if VEUSZ
    for(unsigned int j = 0; j < A->cols(); j++)
    {
        resultFile << setprecision(15) << ", " << j*_hy;
    }
    resultFile <<endl;
    #else
    resultFile << setprecision(15) << _L << "\t" << _D << "\t" << A->rows() << "\t" << A->cols() << "\t" << _hx << "\t" << _hy  <<endl;
    #endif
    for(unsigned int i = 0; i < A->rows(); i++)
    {
        #if VEUSZ
        resultFile << setprecision(15) << - (i*_hx-_D);
        #endif
        for(unsigned int j = 0; j < A->cols(); j++)
        {
           #if VEUSZ
           resultFile << ","<< setprecision(15) << (*A)(i,j);
           #else
           resultFile << setprecision(15) << (*A)(i,j) << "\t";
           #endif
        }
        resultFile <<endl;
    }




    delete A;
}

void saveDisplacementsHaunch(const std::stringstream &_label, double _L, double _H,
                           double _hx, double _hy, unsigned int _n_elemL, unsigned int _n_elemH,
                           NumEquationMatrix *_A, const VectorXd &_x)
{
    std::stringstream outputFileName;
    outputFileName.str("");
    #if VEUSZ
        outputFileName<<"_haunch_resultado_VEUSZ" << _label.str();
    #else
        outputFileName<<"_haunch_resultado" << _label.str();
    #endif
    cout << outputFileName.str()<<endl;
    ofstream resultFile(outputFileName.str(), std::ofstream::out);


    if (!resultFile.is_open())
    {
        cerr<<"[MAIN] o arquivo de resultado não foi aberto"<<endl;
        exit(1);
    }


    resultFile << "\"" << "Label" << "\"\t"
               << "\"" << "x_pos" << "\"\t"
               << "\"" << "y_pos" << "\"\t"
               << "\"" << "x_defl" << "\"\t"
               << "\"" << "y_defl" << "\"\t\n";



    for(unsigned int i = 0; i <= _n_elemH; i++) // linhas de A -> y no sistema de coordenadas
    {
        for(unsigned int j = 0; j <= _n_elemL; j++)// colunas de A -> x no sistema de coordenadas
        {
            if( (*_A)(i,j).x != -1  )
            {
                resultFile << (*_A)(i,j).x << "\t" << j*_hx << "\t" <<_H-(i*_hy) << "\t"
                           <<(_x((*_A)(i,j).x)) <<"\t";
                if((*_A)(i,j).y != -1)
                {
                    resultFile  <<  (_x((*_A)(i,j).y)) << endl;
                }
                else
                {
                    resultFile  <<  0 << endl;
                }

            }
        }
    }

}


void quadModel(unsigned int n_elem, double L, double p, P typeP)
{
    //-------------------------------------------------------------------
    //criando as matrizes/vetores
    MatrixXi *A = new MatrixXi(n_elem+1, n_elem+1);
    A->setConstant(-1);
    unsigned int size = (n_elem-1)*(n_elem-1);
//    VectorXd f = VectorXd::Zero(size);
    SpMatrixD K, f;
    K.resize(size, size);
    K.reserve(size*5);
    f.resize(size, 1);
    f.reserve(size);
    VectorXd analytical = VectorXd::Zero(size);
    vector<T> coeffK, coeffF;
    coeffK.reserve(size*5);
    coeffF.reserve(size);

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

//    cout << *A << endl<< endl<< endl;

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
            else if(typeP == P::PATCH_TEST)
            {
                coeffF.push_back(T((*A)(i, j), 0, -(h*(i-1) + h*j) ));
            }

            if((*A)(i, j-1) >= 0)
                coeffK.push_back(T((*A)(i, j), (*A)(i, j-1), 1));
            else if(typeP == P::PATCH_TEST)
            {
                coeffF.push_back(T((*A)(i, j), 0, -(h*(i) + h*(j-1))  ));
            }

            if((*A)(i+1, j) >= 0)
                coeffK.push_back(T((*A)(i, j), (*A)(i+1, j), 1));
            else if(typeP == P::PATCH_TEST)
            {
                coeffF.push_back(T((*A)(i, j), 0, -(h*(i+1) + h*j) ));
            }

            if((*A)(i, j+1) >= 0)
                coeffK.push_back(T((*A)(i, j), (*A)(i, j+1), 1));
            else if(typeP == P::PATCH_TEST)
            {
                coeffF.push_back(T((*A)(i, j), 0, -(h*(i) + h*(j+1))  ));
            }

            coeffK.push_back(T((*A)(i, j), (*A)(i, j), -4));

            Vector2d pos(h*i, h*j);

            if(typeP != P::PATCH_TEST)
                coeffF.push_back( T((*A)(i, j), 0, h*h * pFunction(typeP, p, pos)) );

            analytical((*A)(i, j)) = calculateAnalytical(typeP, p, pos, L);
        }
    }

    cout << "total memory (real) " << coeffK.size()*sizeof(T)/(1e6) << endl;
    delete A;
    K.setFromTriplets(coeffK.begin(), coeffK.end());
    coeffK.clear();
    K.makeCompressed();

    calculateConditionNumber(MatrixXd(K));


    f.setFromTriplets(coeffF.begin(), coeffF.end());
    coeffF.clear();
    f.makeCompressed();
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

    cout << "Calculado  \t\tAnalitico" << std::fixed << std::setprecision(15)<< endl;
    for(unsigned int i = 0; i < size; i++)
    {
        cout << x(i) << "  \t\t" << analytical(i) << "\t\t" << calculateError(x(i), analytical(i))  << endl;
    }
    cout << "\n\n\nerro global  " << std::scientific << calculateError(x, analytical) <<endl<<endl<<endl;

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

//    cout << "A" << endl << *A << endl;
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
//    cout <<"\n\n"<< MatrixXd(x)<<endl;



    cout <<"acabou processo"<< endl;

    cout << "Salvando o resultado" << endl;
    TimeSpent salvandoTime("Salvo o resultado | TEMPO: ");
    salvandoTime.startCount();
    saveDisplacements(label, L, D, hx, hy, n_elemL, n_elemD, _n_elemL_part, _n_elemD_part, x);
    salvandoTime.endCount();
    salvandoTime.print();

    cout <<"Fim"<< endl;
}

void haunchedBeamModel(double L, double l_haunch, double l_base, double H,
                       double h_haunch, double b, double E, double nu, double p, P typeP,
                       unsigned int n_xb)
{
    //valores calculados para plane strain
    double lambda_planeStrain = (nu*E)/((1+nu)*(1-(2*nu)));
    double mu = E/(2*(1+nu));

    //lambda para plane stress
    double lambda = (2*lambda_planeStrain*mu)/(lambda_planeStrain + 2*mu);

    double hx = l_base/(n_xb - 1);
    double hy = (hx*h_haunch)/(l_haunch - l_base);
    unsigned int n_elemL = ceil(L/hx);
    unsigned int n_elemH = ceil(H/hy);

    unsigned int j_h_min = ceil(((n_xb -1)*l_haunch)/l_base);
    unsigned int j_h_max = n_elemL - j_h_min;

    //testando se vai dar certo
    cout << "tamanho total L: " << hx*n_elemL << " = " << L << endl;
    cout << "tamanho total H: " << hy*n_elemH << " = " << H << endl;
    if( fabs(L - hx*n_elemL) >= 0.0001  )
    {
        cerr << "hx não comporta o L" << endl;
        exit(1);
    }
    if( fabs(H - hy*n_elemH) >= 0.0001  )
    {
        cerr << "hy não comporta o H" << endl;
        exit(1);
    }



    //-------------------------------------------------------------------
    //matriz que contem as equacoes associadas aos nós
    NumEquationMatrix *A = new NumEquationMatrix(n_elemH+1, n_elemL+1);
    //-1 não existe nó

    cout << "Criando matriz auxiliar" << endl;
    TimeSpent criandoATime("Criada matriz auxiliar | TEMPO: ");
    criandoATime.startCount();
    //-------------------------------------------------------------------
    //preenchendo a matriz A que contem o nº das equacoes
    //os valores vão ser só os pares (os impartes serão da outra direção)
    //quem não tem equação são os pontos presos, isso é, os "peszinhos" da mísula
    //  o lado da esquerda é zerado em x e em y
    //  o lado da direita é zerado somente na direção y (não subir)
    //linha -- Y
    //coluna -- X
    unsigned int count = 0;

    //quantidade das linhas da matriz que fazem parte do retangulo de cima
    unsigned int n_elemH_top = ceil((n_elemH * (H - h_haunch))/H);
    //criando as equações do retangulo de cima
    for(unsigned int i = 0; i <= n_elemH_top; i++)
    {
        for(unsigned int j = 0; j <= n_elemL; j++)
        {
            unsigned int borderType;
            //definindo os tipos de borda
            if(j == 0) // j na ponta esquerda
            {
                if(i == 0) //i no topo
                {
                    borderType = BORDER_TYPE::PROG_X_REGR_Y;
                }
                else //i em qualquer outro ponto
                {
                    borderType = BORDER_TYPE::PROG_X_CENT_Y;
                }
            }
            else if(j == n_elemL) //j na ponta direita
            {
                if(i == 0)//i no topo
                {
                    borderType = BORDER_TYPE::REGR_X_REGR_Y;
                }
                else //i em qualquer outro ponto
                {
                    borderType = BORDER_TYPE::REGR_X_CENT_Y;
                }
            }
            else //j é qualquer outro sem ser os da ponta
            {
//                double j_hx = j*hx;
                if(i == 0)//i no topo
                {
                    borderType = BORDER_TYPE::CENT_X_REGR_Y;
                }
                else if( (i == n_elemH_top) && ( j >= j_h_min ) && ( j <= j_h_max ) )
//                    if( (i == n_elemH_top) && ( j_hx >= l_haunch ) && ( j_hx <= (L - l_haunch)) )
                {
                    borderType = BORDER_TYPE::CENT_X_PROG_Y;
                }
                else
                {
                    borderType = BORDER_TYPE::CENT_X_CENT_Y;
                }
            }

            (*A)(i,j) = XY(count, count+1, borderType);
            count += 2;
        }
    }


    //criando a mísula da esquerda
    //parte quadrada
    unsigned int i, j;
    for(i = n_elemH_top +1; i <= n_elemH; i++)
    {
        for(j = 0; j < n_xb-1; j++)
        {
            unsigned int borderType;
            if(i == n_elemH)
            {
                if(j == 0)
                {
                    borderType = BORDER_TYPE::PROG_X_PROG_Y;
                }
                else
                {
                    borderType = BORDER_TYPE::CENT_X_PROG_Y;
                }
            }
            else
            {
                if(j == 0)
                {
                    borderType = BORDER_TYPE::PROG_X_CENT_Y;
                }
                else
                {
                    borderType = BORDER_TYPE::CENT_X_CENT_Y;
                }
            }

            if( !((j == (n_xb-1)/2) && (i == n_elemH)))
            {
                (*A)(i,j) = XY(count, count+1, borderType);
                count += 2;
            }
        }
    }
    //parte triangular
    unsigned int current_j_h_min = j_h_min;
    for(i = n_elemH_top +1; i <= n_elemH; i++)
    {
        for(j = n_xb - 1; j < current_j_h_min; j++)
        {
            unsigned int borderType = BORDER_TYPE::CENT_X_CENT_Y;
            if(j == current_j_h_min-1)
            {
                borderType  = BORDER_TYPE::REGR_X_PROG_Y;
            }

            (*A)(i,j) = XY(count, count+1, borderType);
            count += 2;
        }
        current_j_h_min--;
    }


    //criando a mísula da direita
    //parte quadrada
    for(i = n_elemH_top +1; i <= n_elemH; i++)
    {
        for(j = n_elemL - (n_xb -1) +1; j <= n_elemL; j++)
        {
            unsigned int borderType;
            if(i == n_elemH)
            {
                if(j == n_elemL)
                {
                    borderType = BORDER_TYPE::REGR_X_PROG_Y;
                }
                else
                {
                    borderType = BORDER_TYPE::CENT_X_PROG_Y;
                }
            }
            else
            {
                if(j == n_elemL)
                {
                    borderType = BORDER_TYPE::REGR_X_CENT_Y;
                }
                else
                {
                    borderType = BORDER_TYPE::CENT_X_CENT_Y;
                }
            }

            if( ((j == n_elemL - (n_xb-1)/2) && (i == n_elemH)))
            {
                (*A)(i,j) = XY(count, -1, borderType);
                count += 1;
//                (*A)(i,j) = XY(-1, -1, -1);
            }
            else
            {
                (*A)(i,j) = XY(count, count+1, borderType);
                count += 2;
            }
        }
    }
    //parte triangular
    unsigned int current_j_h_max = j_h_max;
    unsigned int max_j = n_elemL - (n_xb -1);
    for(i = n_elemH_top +1; i <= n_elemH; i++)
    {
        for(j = current_j_h_max+1; j <= max_j; j++)
        {
            unsigned int borderType = BORDER_TYPE::CENT_X_CENT_Y;
            if(j == current_j_h_max+1)
            {
                borderType  = BORDER_TYPE::PROG_X_PROG_Y;
            }

            (*A)(i,j) = XY(count, count+1, borderType);
            count += 2;
        }
        current_j_h_max++;
    }



    criandoATime.endCount();
    criandoATime.print();


//    //printando a matriz A
//    cout << "-----------------------x"<<endl;
//    for(unsigned int i = 0; i <= n_elemH; i++)
//    {
//        for(unsigned int j = 0; j <= n_elemL; j++)
//        {
//            cout << (*A)(i,j).x << " \t ";
//        }
//        cout  << endl;
//    }
//    cout << "-----------------------y"<<endl;
//    for(unsigned int i = 0; i <= n_elemH; i++)
//    {
//        for(unsigned int j = 0; j <= n_elemL; j++)
//        {
//            cout << (*A)(i,j).y << " \t ";
//        }
//        cout  << endl;
//    }
//    cout << "-----------------------tipo"<<endl;
//    for(unsigned int i = 0; i <= n_elemH; i++)
//    {
//        for(unsigned int j = 0; j <= n_elemL; j++)
//        {
//            cout << (*A)(i,j).borderType << " \t ";
//        }
//        cout  << endl;
//    }


    //-------------------------------------------------------------------
    std::stringstream label;
    label.str("");
    label << "_L-"<<L<<"_l-h-"<<l_haunch<<"_l-b-"<<l_base<< "_H-"<<H<<"_h-h-"<<h_haunch << "_b-" << b
          << "_E-" << E << "_nu-" << nu << "_n-xb-" << n_xb
          <<"_p-"<<p<<"_tipo-"<<int(typeP)<<"COEFFK_LU.txt";

    cout << "count " << count << endl;
    VectorXd f = VectorXd::Zero(count);
    SpMatrixD K;
    K.resize(count, count);
    vector<T> coeffK;


    double hx_2 = hx*hx;
    double hy_2 = hy*hy;

    cout << "Preenchendo o sistema" << endl;
    TimeSpent preenchendoSistemaTime("Preenchido o sistema | TEMPO: ");
    preenchendoSistemaTime.startCount();
    //-------------------------------------------------------------------
    //preenchendo o sistema
    double forcePerNode = p/n_elemL;

    for(unsigned int j = 0; j <= n_elemL; j++) // colunas de A -> x no sistema de coordenadas
    {

        for(unsigned int i = 0; i <= n_elemH; i++)// linhas de A -> y no sistema de coordenadas
        {

            switch ((*A)(i, j).borderType)
            {
                case BORDER_TYPE::PROG_X_REGR_Y:
                {
                    /*
                         |0|-(3)-(4)
                          |  \
                         (1) (5)
                          |
                         (2)
                    */

                    f.coeffRef((*A)(i, j).y, 0) =  forcePerNode/2.;

                    //adicionando as informações da direção X e Y
                    coeffK.push_back(T((*A)(i, j).x, (*A)(i, j).x, (hy_2*lambda-hx*hy*lambda+2*hy_2*mu-hx*hy*mu+hx_2*mu)/(hx_2*hy_2)));
                    coeffK.push_back(T((*A)(i, j).y, (*A)(i, j).y, -(hx*hy*lambda-hx_2*lambda-hy_2*mu+hx*hy*mu-2*hx_2*mu)/(hx_2*hy_2)));

                    if((*A)(i+1, j).x >= 0) //baixo e meio 1
                    {
                        coeffK.push_back(T((*A)(i, j).x, (*A)(i+1, j).x, (hy*lambda+hy*mu-2*hx*mu)/(hx*hy_2)));
                        coeffK.push_back(T((*A)(i, j).y, (*A)(i+1, j).y, (hy*lambda-2*hx*lambda+hy*mu-4*hx*mu)/(hx*hy_2)));
                    }
                    else
                    {
                        cout << "[PROG_X_REGR_Y] errou em ((*A)(i+1, j).x >= 0) " << i +1 << " - " << j +1 << endl;
                    }

                    if((*A)(i+2, j).x >= 0) //baixo baixo e meio 2
                    {
                        coeffK.push_back(T((*A)(i, j).x, (*A)(i+2, j).x, mu/hy_2));
                        coeffK.push_back(T((*A)(i, j).y, (*A)(i+2, j).y, (lambda+2*mu)/hy_2));
                    }
                    else
                    {
                        cout << "[PROG_X_REGR_Y] errou em ((*A)(i+2, j).x >= 0) " << i +1 << " - " << j +1 << endl;
                    }

                    if((*A)(i, j+1).x >= 0) //meio e direita 3
                    {
                        coeffK.push_back(T((*A)(i, j).x, (*A)(i, j+1).x, -(2*hy*lambda-hx*lambda+4*hy*mu-hx*mu)/(hx_2*hy)));
                        coeffK.push_back(T((*A)(i, j).y, (*A)(i, j+1).y, (hx*lambda-2*hy*mu+hx*mu)/(hx_2*hy)));
                    }
                    else
                    {
                        cout << "[PROG_X_REGR_Y] errou em ((*A)(i, j+1).x >= 0) " << i +1 << " - " << j +1 << endl;
                    }

                    if((*A)(i, j+2).x >= 0) //meio e direita direita 4
                    {
                        coeffK.push_back(T((*A)(i, j).x, (*A)(i, j+2).x, (lambda+2*mu)/hx_2));
                        coeffK.push_back(T((*A)(i, j).y, (*A)(i, j+2).y, mu/hx_2));
                    }
                    else
                    {
                        cout << "[PROG_X_REGR_Y] errou em ((*A)(i, j+2).x >= 0) " << i +1 << " - " << j +1 << endl;
                    }

                    if((*A)(i+1, j+1).x >= 0) //diagonal baixo e direita 5
                    {
                        coeffK.push_back(T((*A)(i, j).x, (*A)(i+1, j+1).x, -(lambda + mu)/(hx*hy)));
                        coeffK.push_back(T((*A)(i, j).y, (*A)(i+1, j+1).y, -(lambda + mu)/(hx*hy)));
                    }
                    else
                    {
                        cout << "[PROG_X_REGR_Y] errou em ((*A)(i+1, j+1).x >= 0) " << i +1 << " - " << j +1 << endl;
                    }
                    break;
                }
                case BORDER_TYPE::REGR_X_REGR_Y:
                {
                    /*
                    (4)-(3)-|0|
                          /  |
                        (5) (1)
                             |
                            (2)
                    */


                    f.coeffRef((*A)(i, j).y, 0) =  forcePerNode/2.;

                    //adicionando as informações da direção X e Y
                    coeffK.push_back(T((*A)(i, j).x, (*A)(i, j).x, (hy_2*lambda+hx*hy*lambda+2*hy_2*mu+hx*hy*mu+hx_2*mu)/(hx_2*hy_2)));
                    coeffK.push_back(T((*A)(i, j).y, (*A)(i, j).y, (hx*hy*lambda+hx_2*lambda+hy_2*mu+hx*hy*mu+2*hx_2*mu)/(hx_2*hy_2)));

                    if((*A)(i+1, j).x >= 0) //baixo e meio 1
                    {
                        coeffK.push_back(T((*A)(i, j).x, (*A)(i+1, j).x, -(hy*lambda+hy*mu+2*hx*mu)/(hx*hy_2)));
                        coeffK.push_back(T((*A)(i, j).y, (*A)(i+1, j).y, -(hy*lambda+2*hx*lambda+hy*mu+4*hx*mu)/(hx*hy_2)));
                    }
                    else
                    {
                        cout << "[REGR_X_REGR_Y] errou em ((*A)(i+1, j).x >= 0) " << i +1 << " - " << j +1 << endl;
                    }

                    if((*A)(i+2, j).x >= 0) //baixo baixo e meio 1
                    {
                        coeffK.push_back(T((*A)(i, j).x, (*A)(i+2, j).x, mu/hy_2));
                        coeffK.push_back(T((*A)(i, j).y, (*A)(i+2, j).y, (lambda+2*mu)/hy_2));
                    }
                    else
                    {
                        cout << "[REGR_X_REGR_Y] errou em ((*A)(i+2, j).x >= 0) " << i +1 << " - " << j +1 << endl;
                    }

                    if((*A)(i, j-1).x >= 0) //meio e esquerda 3
                    {
                        coeffK.push_back(T((*A)(i, j).x, (*A)(i, j-1).x, -(2*hy*lambda+hx*lambda+4*hy*mu+hx*mu)/(hx_2*hy)));
                        coeffK.push_back(T((*A)(i, j).y, (*A)(i, j-1).y, -(hx*lambda+2*hy*mu+hx*mu)/(hx_2*hy)));
                    }
                    else
                    {
                        cout << "[REGR_X_REGR_Y] errou em ((*A)(i, j-1).x >= 0) " << i +1 << " - " << j +1 << endl;
                    }

                    if((*A)(i, j-2).x >= 0) //meio e esquerda esquerda 4
                    {
                        coeffK.push_back(T((*A)(i, j).x, (*A)(i, j-2).x, (lambda+2*mu)/hx_2));
                        coeffK.push_back(T((*A)(i, j).y, (*A)(i, j-2).y, mu/hx_2));
                    }
                    else
                    {
                        cout << "[REGR_X_REGR_Y] errou em ((*A)(i, j-2).x >= 0) " << i +1 << " - " << j +1 << endl;
                    }


                    if((*A)(i+1, j-1).x >= 0) //diagonal baixo e esquerda 5
                    {
                        coeffK.push_back(T((*A)(i, j).x, (*A)(i+1, j-1).x, (lambda + mu)/(hx*hy)));
                        coeffK.push_back(T((*A)(i, j).y, (*A)(i+1, j-1).y, (lambda + mu)/(hx*hy)));
                    }
                    else
                    {
                        cout << "[REGR_X_REGR_Y] errou em ((*A)(i+1, j-1).x >= 0) " << i +1 << " - " << j +1 << endl;
                    }
                    break;
                }
                case BORDER_TYPE::CENT_X_REGR_Y:
                {
                    /*
                       (3)-|0|-(4)
                         /  |  \
                       (5) (1) (6)
                            |
                           (2)
                    */

                    f.coeffRef((*A)(i, j).y, 0) =  forcePerNode;

                    //adicionando as informações da direção X e Y
                    coeffK.push_back(T((*A)(i, j).x, (*A)(i, j).x, -(2*hy_2*lambda+4*hy_2*mu-hx_2*mu)/(hx_2*hy_2)));
                    coeffK.push_back(T((*A)(i, j).y, (*A)(i, j).y, (hx_2*lambda-2*hy_2*mu+2*hx_2*mu)/(hx_2*hy_2)));

                    if((*A)(i+1, j).x >= 0) //baixo e meio 1
                    {
                        coeffK.push_back(T((*A)(i, j).x, (*A)(i+1, j).x, -2*mu/hy_2));
                        coeffK.push_back(T((*A)(i, j).y, (*A)(i+1, j).y, -2*(lambda+2*mu)/hy_2));
                    }
                    else
                    {
                        cout << "[CENT_X_REGR_Y] errou em ((*A)(i+1, j).x >= 0) " << i +1 << " - " << j +1 << endl;
                    }

                    if((*A)(i+2, j).x >= 0) //baixo baixo e meio 2
                    {
                        coeffK.push_back(T((*A)(i, j).x, (*A)(i+2, j).x, mu/hy_2));
                        coeffK.push_back(T((*A)(i, j).y, (*A)(i+2, j).y, (lambda+2*mu)/hy_2));
                    }
                    else
                    {
                        cout << "[CENT_X_REGR_Y] errou em ((*A)(i+2, j).x >= 0) " << i +1 << " - " << j +1 << endl;
                    }


                    if((*A)(i, j-1).x >= 0) //meio e esquerda 3
                    {
                        coeffK.push_back(T((*A)(i, j).x, (*A)(i, j-1).x, (lambda+2*mu)/hx_2 - (lambda+mu)/(2*hx*hy)));
                        coeffK.push_back(T((*A)(i, j).y, (*A)(i, j-1).y, mu/hx_2 - (lambda+mu)/(2*hx*hy)));
                    }
                    else
                    {
                        cout << "[CENT_X_REGR_Y] errou em ((*A)(i, j-1).x >= 0) " << i +1 << " - " << j +1 << endl;
                    }

                    if((*A)(i, j+1).x >= 0) //meio e direita 4
                    {
                        coeffK.push_back(T((*A)(i, j).x, (*A)(i, j+1).x, (lambda+mu)/(2*hx*hy) + (lambda+2*mu)/hx_2));
                        coeffK.push_back(T((*A)(i, j).y, (*A)(i, j+1).y, (lambda+mu)/(2*hx*hy) + mu/hx_2));
                    }
                    else
                    {
                        cout << "[CENT_X_REGR_Y] errou em ((*A)(i, j+1).x >= 0) " << i +1 << " - " << j +1 << endl;
                    }


                    if((*A)(i+1, j-1).x >= 0) //diagonal baixo e esquerda 5
                    {
                        coeffK.push_back(T((*A)(i, j).x, (*A)(i+1, j-1).x, (lambda + mu)/(2*hx*hy)));
                        coeffK.push_back(T((*A)(i, j).y, (*A)(i+1, j-1).y, (lambda + mu)/(2*hx*hy)));
                    }
                    else
                    {
                        cout << "[CENT_X_REGR_Y] errou em ((*A)(i+1, j-1).x >= 0) " << i +1 << " - " << j +1 << endl;
                    }

                    if((*A)(i+1, j+1).x >= 0) //diagonal baixo e direita 6
                    {
                        coeffK.push_back(T((*A)(i, j).x, (*A)(i+1, j+1).x, -(lambda + mu)/(2*hx*hy)));
                        coeffK.push_back(T((*A)(i, j).y, (*A)(i+1, j+1).y, -(lambda + mu)/(2*hx*hy)));
                    }
                    else
                    {
                        cout << "[CENT_X_REGR_Y] errou em ((*A)(i+1, j+1).x >= 0) " << i +1 << " - " << j +1 << endl;
                    }
                    break;
                }
                case BORDER_TYPE::PROG_X_CENT_Y:
                {
                    /*
                       (3) (5)
                        |  /
                       |0|-(1)-(2)
                        |  \
                       (4) (6)
                    */


                    //adicionando as informações da direção X e Y
                    coeffK.push_back(T((*A)(i, j).x, (*A)(i, j).x, (hy_2*lambda+2*hy_2*mu-2*hx_2*mu)/(hx_2*hy_2)));
                    coeffK.push_back(T((*A)(i, j).y, (*A)(i, j).y, -(2*hx_2*lambda-hy_2*mu+4*hx_2*mu)/(hx_2*hy_2)));


                    if((*A)(i, j+1).x >= 0) //meio e direita 1
                    {
                        coeffK.push_back(T((*A)(i, j).x, (*A)(i, j+1).x, -2*(lambda+2*mu)/hx_2));
                        coeffK.push_back(T((*A)(i, j).y, (*A)(i, j+1).y, -2*mu/hx_2));
                    }
                    else
                    {
                        cout << "[PROG_X_CENT_Y] errou em ((*A)(i, j+1).x >= 0) " << i +1 << " - " << j +1 << endl;
                    }


                    if((*A)(i, j+2).x >= 0) //meio e direita direita 2
                    {
                        coeffK.push_back(T((*A)(i, j).x, (*A)(i, j+2).x, (lambda+2*mu)/hx_2));
                        coeffK.push_back(T((*A)(i, j).y, (*A)(i, j+2).y, mu/hx_2));
                    }
                    else
                    {
                        cout << "[PROG_X_CENT_Y] errou em ((*A)(i, j+2).x >= 0) " << i +1 << " - " << j +1 << endl;
                    }


                    if((*A)(i-1, j).x >= 0) //cima e meio 3
                    {
                        coeffK.push_back(T((*A)(i, j).x, (*A)(i-1, j).x, mu/hy_2 -(lambda+mu)/(2*hx*hy)));
                        coeffK.push_back(T((*A)(i, j).y, (*A)(i-1, j).y, (lambda+2*mu)/hy_2 - (lambda+mu)/(2*hx*hy)));
                    }
                    else
                    {
                        cout << "[PROG_X_CENT_Y] errou em ((*A)(i-1, j).x >= 0) " << i +1 << " - " << j +1 << endl;
                    }

                    if((*A)(i+1, j).x >= 0) //baixo e meio 4
                    {
                        coeffK.push_back(T((*A)(i, j).x, (*A)(i+1, j).x, mu/hy_2 + (lambda+mu)/(2*hx*hy)));
                        coeffK.push_back(T((*A)(i, j).y, (*A)(i+1, j).y, (lambda+2*mu)/hy_2 + (lambda+mu)/(2*hx*hy)));
                    }
                    else
                    {
                        cout << "[PROG_X_CENT_Y] errou em ((*A)(i+1, j).x >= 0) " << i +1 << " - " << j +1 << endl;
                    }

                    if((*A)(i-1, j+1).x >= 0) //diagonal cima e direita 5
                    {
                        coeffK.push_back(T((*A)(i, j).x, (*A)(i-1, j+1).x, (lambda + mu)/(2*hx*hy)));
                        coeffK.push_back(T((*A)(i, j).y, (*A)(i-1, j+1).y, (lambda + mu)/(2*hx*hy)));
                    }
                    else
                    {
                        cout << "[PROG_X_CENT_Y] errou em ((*A)(i-1, j+1).x >= 0) " << i +1 << " - " << j +1 << endl;
                    }

                    if((*A)(i+1, j+1).x >= 0) //diagonal baixo e direita 6
                    {
                        coeffK.push_back(T((*A)(i, j).x, (*A)(i+1, j+1).x, -(lambda + mu)/(2*hx*hy)));
                        coeffK.push_back(T((*A)(i, j).y, (*A)(i+1, j+1).y, -(lambda + mu)/(2*hx*hy)));
                    }
                    else
                    {
                        cout << "[PROG_X_CENT_Y] errou em ((*A)(i+1, j+1).x >= 0) " << i +1 << " - " << j +1 << endl;
                    }
                    break;
                }
                case BORDER_TYPE::REGR_X_CENT_Y:
                {
                    /*
                        (5) (3)
                          \  |
                    (2)-(1)-|0|
                          /  |
                        (6) (4)
                    */


                    //adicionando as informações da direção X e Y
                    coeffK.push_back(T((*A)(i, j).x, (*A)(i, j).x, (hy_2*lambda+2*hy_2*mu-2*hx_2*mu)/(hx_2*hy_2)));
                    coeffK.push_back(T((*A)(i, j).y, (*A)(i, j).y, -(2*hx_2*lambda-hy_2*mu+4*hx_2*mu)/(hx_2*hy_2)));


                    if((*A)(i, j-1).x >= 0) //meio e esquerda 1
                    {
                        coeffK.push_back(T((*A)(i, j).x, (*A)(i, j-1).x, -2*(lambda+2*mu)/hx_2));
                    }
                    else
                    {
                        cout << "[REGR_X_CENT_Y] errou em ((*A)(i, j-1).x >= 0) " << i +1 << " - " << j +1 << endl;
                    }
                    if((*A)(i, j-1).y >= 0) //meio e esquerda 1
                    {
                        coeffK.push_back(T((*A)(i, j).y, (*A)(i, j-1).y, -2*mu/hx_2));
                    }
                    else
                    {
                        cout << "[REGR_X_CENT_Y] errou em ((*A)(i, j-1).y >= 0) " << i +1 << " - " << j +1 << endl;
                    }

                    if((*A)(i, j-2).x >= 0) //meio e esquerda esquerda 2
                    {
                        coeffK.push_back(T((*A)(i, j).x, (*A)(i, j-2).x, (lambda+2*mu)/hx_2));
                    }
                    else
                    {
                        cout << "[REGR_X_CENT_Y] errou em ((*A)(i, j-2).x >= 0) " << i +1 << " - " << j +1 << endl;
                    }
                    if((*A)(i, j-2).y >= 0) //meio e esquerda esquerda 2
                    {
                        coeffK.push_back(T((*A)(i, j).y, (*A)(i, j-2).y, mu/hx_2));
                    }
                    else
                    {
                        cout << "[REGR_X_CENT_Y] errou em ((*A)(i, j-2).y >= 0) " << i +1 << " - " << j +1 << endl;
                    }

                    if((*A)(i-1, j).x >= 0) //cima e meio 3
                    {
                        coeffK.push_back(T((*A)(i, j).x, (*A)(i-1, j).x, mu/hy_2 + (lambda+mu)/(2*hx*hy)));
                    }
                    else
                    {
                        cout << "[REGR_X_CENT_Y] errou em ((*A)(i-1, j).x >= 0) " << i +1 << " - " << j +1 << endl;
                    }
                    if((*A)(i-1, j).y >= 0) //cima e meio 3
                    {
                        coeffK.push_back(T((*A)(i, j).y, (*A)(i-1, j).y, (lambda+2*mu)/hy_2 + (lambda+mu)/(2*hx*hy)));
                    }
                    else
                    {
                        cout << "[REGR_X_CENT_Y] errou em ((*A)(i-1, j).y >= 0) " << i +1 << " - " << j +1 << endl;
                    }

                    if((*A)(i+1, j).x >= 0) //baixo e meio 4
                    {
                        coeffK.push_back(T((*A)(i, j).x, (*A)(i+1, j).x, mu/hy_2 - (lambda+mu)/(2*hx*hy)));
                    }
                    else
                    {
                        cout << "[REGR_X_CENT_Y] errou em ((*A)(i+1, j).x >= 0) " << i +1 << " - " << j +1 << endl;
                    }
                    if((*A)(i+1, j).y >= 0) //baixo e meio 4
                    {
                        coeffK.push_back(T((*A)(i, j).y, (*A)(i+1, j).y, (lambda+2*mu)/hy_2 - (lambda+mu)/(2*hx*hy)));
                    }
                    else
                    {
                        cout << "[REGR_X_CENT_Y] errou em ((*A)(i+1, j).y >= 0) " << i +1 << " - " << j +1 << endl;
                    }

                    if((*A)(i-1, j-1).x >= 0) //diagonal cima e esquerda 5
                    {
                        coeffK.push_back(T((*A)(i, j).x, (*A)(i-1, j-1).x, -(lambda + mu)/(2*hx*hy)));
                    }
                    else
                    {
                        cout << "[REGR_X_CENT_Y] errou em ((*A)(i-1, j-1).x >= 0) " << i +1 << " - " << j +1 << endl;
                    }
                    if((*A)(i-1, j-1).y >= 0) //diagonal cima e esquerda 5
                    {
                        coeffK.push_back(T((*A)(i, j).y, (*A)(i-1, j-1).y, -(lambda + mu)/(2*hx*hy)));
                    }
                    else
                    {
                        cout << "[REGR_X_CENT_Y] errou em ((*A)(i-1, j-1).y >= 0) " << i +1 << " - " << j +1 << endl;
                    }

                    if((*A)(i+1, j-1).x >= 0) //diagonal baixo e esquerda 6
                    {
                        coeffK.push_back(T((*A)(i, j).x, (*A)(i+1, j-1).x, (lambda + mu)/(2*hx*hy)));
                    }
                    else
                    {
                        cout << "[REGR_X_CENT_Y] errou em ((*A)(i+1, j-1).x >= 0) " << i +1 << " - " << j +1 << endl;
                    }
                    if((*A)(i+1, j-1).y >= 0) //diagonal baixo e esquerda 6
                    {
                        coeffK.push_back(T((*A)(i, j).y, (*A)(i+1, j-1).y, (lambda + mu)/(2*hx*hy)));
                    }
                    else
                    {
                        cout << "[REGR_X_CENT_Y] errou em ((*A)(i+1, j-1).y >= 0) " << i +1 << " - " << j +1 << endl;
                    }

                    break;
                }
                case BORDER_TYPE::CENT_X_CENT_Y:
                {
                    /*
                       (5) (1) (7)
                         \  |  /
                       (3)-|0|-(4)
                         /  |  \
                       (6) (2) (8)
                    */


                    //adicionando as informações da direção X e Y
                    coeffK.push_back(T((*A)(i, j).x, (*A)(i, j).x, -2*(hy_2*lambda+2*hy_2*mu+hx_2*mu)/(hx_2*hy_2)));
                    coeffK.push_back(T((*A)(i, j).y, (*A)(i, j).y, -2*(hx_2*lambda+hy_2*mu+2*hx_2*mu)/(hx_2*hy_2)));

                    if((*A)(i-1, j).x >= 0) //cima e meio 1
                    {
                        coeffK.push_back(T((*A)(i, j).x, (*A)(i-1, j).x, mu/hy_2));
                    }
                    else
                    {
                        cout << "[CENT_X_CENT_Y] errou em ((*A)(i-1, j).x >= 0) " << i +1 << " - " << j +1 << endl;
                    }
                    if((*A)(i-1, j).y >= 0) //cima e meio 1
                    {
                        coeffK.push_back(T((*A)(i, j).y, (*A)(i-1, j).y, (lambda+2*mu)/hy_2));
                    }
                    else
                    {
                        cout << "[CENT_X_CENT_Y] errou em ((*A)(i-1, j).y >= 0) " << i +1 << " - " << j +1 << endl;
                    }

                    if((*A)(i+1, j).x >= 0) //baixo e meio 2
                    {
                        coeffK.push_back(T((*A)(i, j).x, (*A)(i+1, j).x, mu/hy_2));
                    }
                    else
                    {
                        cout << "[CENT_X_CENT_Y] errou em ((*A)(i+1, j).x >= 0) " << i +1 << " - " << j +1 << endl;
                    }
                    if((*A)(i+1, j).y >= 0) //baixo e meio 2
                    {
                        coeffK.push_back(T((*A)(i, j).y, (*A)(i+1, j).y, (lambda+2*mu)/hy_2));
                    }
                    else
                    {
                        cout << "[CENT_X_CENT_Y] errou em ((*A)(i+1, j).y >= 0) " << i +1 << " - " << j +1 << endl;
                    }

                    if((*A)(i, j-1).x >= 0) //meio e esquerda 3
                    {
                        coeffK.push_back(T((*A)(i, j).x, (*A)(i, j-1).x, (lambda+2*mu)/hx_2));
                    }
                    else
                    {
                        cout << "[CENT_X_CENT_Y] errou em ((*A)(i, j-1).x >= 0) " << i +1 << " - " << j +1 << endl;
                    }
                    if((*A)(i, j-1).y >= 0) //meio e esquerda 3
                    {
                        coeffK.push_back(T((*A)(i, j).y, (*A)(i, j-1).y, mu/hx_2));
                    }
                    else
                    {
                        cout << "[CENT_X_CENT_Y] errou em ((*A)(i, j-1).y >= 0) " << i +1 << " - " << j +1 << endl;
                    }

                    if((*A)(i, j+1).x >= 0) //meio e direita 4
                    {
                        coeffK.push_back(T((*A)(i, j).x, (*A)(i, j+1).x, (lambda+2*mu)/hx_2));
                    }
                    else
                    {
                        cout << "[CENT_X_CENT_Y] errou em ((*A)(i, j+1).x >= 0) " << i +1 << " - " << j +1 << endl;
                    }
                    if((*A)(i, j+1).y >= 0) //meio e direita 4
                    {
                        coeffK.push_back(T((*A)(i, j).y, (*A)(i, j+1).y, mu/hx_2));
                    }
                    else
                    {
                        cout << "[CENT_X_CENT_Y] errou em ((*A)(i, j+1).y >= 0) " << i +1 << " - " << j +1 << endl;
                    }

                    if((*A)(i-1, j-1).y >= 0) //diagonal cima e esquerda 5
                    {
                        coeffK.push_back(T((*A)(i, j).x, (*A)(i-1, j-1).y, -(lambda + mu)/(4*hx*hy)));
                    }
                    else
                    {
                        cout << "[CENT_X_CENT_Y] errou em ((*A)(i-1, j-1).x >= 0) " << i +1 << " - " << j +1 << endl;
                    }
                    if((*A)(i-1, j-1).x >= 0) //diagonal cima e esquerda 5
                    {
                        coeffK.push_back(T((*A)(i, j).y, (*A)(i-1, j-1).x, -(lambda + mu)/(4*hx*hy)));
                    }
                    else
                    {
                        cout << "[CENT_X_CENT_Y] errou em ((*A)(i-1, j-1).y >= 0) " << i +1 << " - " << j +1 << endl;
                    }

                    if((*A)(i+1, j-1).y >= 0) //diagonal baixo e esquerda 6
                    {
                        coeffK.push_back(T((*A)(i, j).x, (*A)(i+1, j-1).y, (lambda + mu)/(4*hx*hy)));
                    }
                    else
                    {
                        cout << "[CENT_X_CENT_Y] errou em ((*A)(i+1, j-1).x >= 0) " << i +1 << " - " << j +1 << endl;
                    }
                    if((*A)(i+1, j-1).x >= 0) //diagonal baixo e esquerda 6
                    {
                        coeffK.push_back(T((*A)(i, j).y, (*A)(i+1, j-1).x, (lambda + mu)/(4*hx*hy)));
                    }
                    else
                    {
                        cout << "[CENT_X_CENT_Y] errou em ((*A)(i+1, j-1).y >= 0) " << i +1 << " - " << j +1 << endl;
                    }

                    if((*A)(i-1, j+1).y >= 0) //diagonal cima e direita 7
                    {
                        coeffK.push_back(T((*A)(i, j).x, (*A)(i-1, j+1).y, (lambda + mu)/(4*hx*hy)));
                    }
                    else
                    {
                        cout << "[CENT_X_CENT_Y] errou em ((*A)(i-1, j+1).x >= 0) " << i +1 << " - " << j +1 << endl;
                    }
                    if((*A)(i-1, j+1).x >= 0) //diagonal cima e direita 7
                    {
                        coeffK.push_back(T((*A)(i, j).y, (*A)(i-1, j+1).x, (lambda + mu)/(4*hx*hy)));
                    }
                    else
                    {
                        cout << "[CENT_X_CENT_Y] errou em ((*A)(i-1, j+1).x >= 0) " << i +1 << " - " << j +1 << endl;
                    }

                    if((*A)(i+1, j+1).y >= 0) //diagonal baixo e direita 8
                    {
                        coeffK.push_back(T((*A)(i, j).x, (*A)(i+1, j+1).y, -(lambda + mu)/(4*hx*hy)));
                    }
                    else
                    {
                        cout << "[CENT_X_CENT_Y] errou em ((*A)(i+1, j+1).x >= 0) " << i +1 << " - " << j +1 << endl;
                    }
                    if((*A)(i+1, j+1).x >= 0) //diagonal baixo e direita 8
                    {
                        coeffK.push_back(T((*A)(i, j).y, (*A)(i+1, j+1).x, -(lambda + mu)/(4*hx*hy)));
                    }
                    else
                    {
                        cout << "[CENT_X_CENT_Y] errou em ((*A)(i+1, j+1).y >= 0) " << i +1 << " - " << j +1 << endl;
                    }
                    break;
                }
                case BORDER_TYPE::PROG_X_PROG_Y:
                {
                    /*
                       (2)
                        |
                       (1) (5)
                        |  /
                       |0|-(3)-(4)
                    */


                    //adicionando as informações da direção X e Y
                    coeffK.push_back(T((*A)(i, j).x, (*A)(i, j).x, (hy_2*lambda+hx*hy*lambda+2*hy_2*mu+hx*hy*mu+hx_2*mu)/(hx_2*hy_2)));
                    coeffK.push_back(T((*A)(i, j).y, (*A)(i, j).y, (hx*hy*lambda+hx_2*lambda+hy_2*mu+hx*hy*mu+2*hx_2*mu)/(hx_2*hy_2)));

                    if((*A)(i-1, j).x >= 0) //cima e meio 1
                    {
                        coeffK.push_back(T((*A)(i, j).x, (*A)(i-1, j).x, -(hy*lambda+hy*mu+2*hx*mu)/(hx*hy_2)));
                    }
                    else
                    {
                        cout << "[PROG_X_PROG_Y] errou em ((*A)(i-1, j).x >= 0) " << i +1 << " - " << j +1 << endl;
                    }
                    if((*A)(i-1, j).y >= 0) //cima e meio 1
                    {
                        coeffK.push_back(T((*A)(i, j).y, (*A)(i-1, j).y, -(hy*lambda+2*hx*lambda+hy*mu+4*hx*mu)/(hx*hy_2)));
                    }
                    else
                    {
                        cout << "[PROG_X_PROG_Y] errou em ((*A)(i-1, j).y >= 0) " << i +1 << " - " << j +1 << endl;
                    }


                    if((*A)(i-2, j).x >= 0) //cima cima e meio 2
                    {
                        coeffK.push_back(T((*A)(i, j).x, (*A)(i-2, j).x, mu/hy_2));
                    }
                    else
                    {
                        cout << "[PROG_X_PROG_Y] errou em ((*A)(i-2, j).x >= 0) " << i +1 << " - " << j +1 << endl;
                    }
                    if((*A)(i-2, j).y >= 0) //cima cima e meio 2
                    {
                        coeffK.push_back(T((*A)(i, j).y, (*A)(i-2, j).y, (lambda+2*mu)/hy_2));
                    }
                    else
                    {
                        cout << "[PROG_X_PROG_Y] errou em ((*A)(i-2, j).y >= 0) " << i +1 << " - " << j +1 << endl;
                    }


                    if((*A)(i, j+1).x >= 0) //meio e direita 3
                    {
                        coeffK.push_back(T((*A)(i, j).x, (*A)(i, j+1).x, -(2*hy*lambda+hx*lambda+4*hy*mu+hx*mu)/(hx_2*hy)));
                    }
                    else
                    {
                        cout << "[PROG_X_PROG_Y] errou em ((*A)(i, j+1).x >= 0) " << i +1 << " - " << j +1 << endl;
                    }
                    if((*A)(i, j+1).y >= 0) //meio e direita 3
                    {
                        coeffK.push_back(T((*A)(i, j).y, (*A)(i, j+1).y, -(hx*lambda+2*hy*mu+hx*mu)/(hx_2*hy)));
                    }
                    else
                    {
                        cout << "[PROG_X_PROG_Y] errou em ((*A)(i, j+1).y >= 0) " << i +1 << " - " << j +1 << endl;
                    }


                    if((*A)(i, j+2).x >= 0) //meio e direita direita 4
                    {
                        coeffK.push_back(T((*A)(i, j).x, (*A)(i, j+2).x, (lambda+2*mu)/hx_2));
                    }
                    else
                    {
                        cout << "[PROG_X_PROG_Y] errou em ((*A)(i, j+2).x >= 0) " << i +1 << " - " << j +1 << endl;
                    }
                    if((*A)(i, j+2).y >= 0) //meio e direita direita 4
                    {
                        coeffK.push_back(T((*A)(i, j).y, (*A)(i, j+2).y, mu/hx_2));
                    }
                    else
                    {
                        cout << "[PROG_X_PROG_Y] errou em ((*A)(i, j+2).y >= 0) " << i +1 << " - " << j +1 << endl;
                    }

                    if((*A)(i-1, j+1).x >= 0) //diagonal cima e direita 5
                    {
                        coeffK.push_back(T((*A)(i, j).x, (*A)(i-1, j+1).x, (lambda + mu)/(hx*hy)));
                    }
                    else
                    {
                        cout << "[PROG_X_PROG_Y] errou em ((*A)(i-1, j+1).x >= 0) " << i +1 << " - " << j +1 << endl;
                    }
                    if((*A)(i-1, j+1).y >= 0) //diagonal cima e direita 5
                    {
                        coeffK.push_back(T((*A)(i, j).y, (*A)(i-1, j+1).y, (lambda + mu)/(hx*hy)));
                    }
                    else
                    {
                        cout << "[PROG_X_PROG_Y] errou em ((*A)(i-1, j+1).y >= 0) " << i +1 << " - " << j +1 << endl;
                    }
                    break;
                }
                case BORDER_TYPE::REGR_X_PROG_Y:
                {
                    /*
                            (2)
                             |
                        (5) (1)
                          \  |
                    (4)-(3)-|0|
                    */


                    //adicionando as informações da direção X e Y
                    coeffK.push_back(T((*A)(i, j).x, (*A)(i, j).x, (hy_2*lambda-hx*hy*lambda+2*hy_2*mu-hx*hy*mu+hx_2*mu)/(hx_2*hy_2)));
                    coeffK.push_back(T((*A)(i, j).y, (*A)(i, j).y, -(hx*hy*lambda-hx_2*lambda-hy_2*mu+hx*hy*mu-2*hx_2*mu)/(hx_2*hy_2)));

                    if((*A)(i-1, j).x >= 0) //cima e meio 1
                    {
                        coeffK.push_back(T((*A)(i, j).x, (*A)(i-1, j).x, (hy*lambda+hy*mu-2*hx*mu)/(hx*hy_2)));
                    }
                    else
                    {
                        cout << "[REGR_X_PROG_Y] errou em ((*A)(i-1, j).x >= 0) " << i +1 << " - " << j +1 << endl;
                    }
                    if((*A)(i-1, j).y >= 0) //cima e meio 1
                    {
                        coeffK.push_back(T((*A)(i, j).y, (*A)(i-1, j).y, (hy*lambda-2*hx*lambda+hy*mu-4*hx*mu)/(hx*hy_2)));
                    }
                    else
                    {
                        cout << "[REGR_X_PROG_Y] errou em ((*A)(i-1, j).y >= 0) " << i +1 << " - " << j +1 << endl;
                    }

                    if((*A)(i-2, j).x >= 0) //cima cima e meio 2
                    {
                        coeffK.push_back(T((*A)(i, j).x, (*A)(i-2, j).x, mu/hy_2));
                    }
                    else
                    {
                        cout << "[REGR_X_PROG_Y] errou em ((*A)(i-2, j).x >= 0) " << i +1 << " - " << j +1 << endl;
                    }                    
                    if((*A)(i-2, j).y >= 0) //cima cima e meio 2
                    {
                        coeffK.push_back(T((*A)(i, j).y, (*A)(i-2, j).y, (lambda+2*mu)/hy_2));
                    }
                    else
                    {
                        cout << "[REGR_X_PROG_Y] errou em ((*A)(i-2, j).y >= 0) " << i +1 << " - " << j +1 << endl;
                    }

                    if((*A)(i, j-1).x >= 0) //meio e esquerda 3
                    {
                        coeffK.push_back(T((*A)(i, j).x, (*A)(i, j-1).x, -(2*hy*lambda-hx*lambda+4*hy*mu-hx*mu)/(hx_2*hy)));
                    }
                    else
                    {
                        cout << "[REGR_X_PROG_Y] errou em ((*A)(i, j-1).x >= 0) " << i +1 << " - " << j +1 << endl;
                    }                    
                    if((*A)(i, j-1).y >= 0) //meio e esquerda 3
                    {
                        coeffK.push_back(T((*A)(i, j).y, (*A)(i, j-1).y, (hx*lambda-2*hy*mu+hx*mu)/(hx_2*hy)));
                    }
                    else
                    {
                        cout << "[REGR_X_PROG_Y] errou em ((*A)(i, j-1).y >= 0) " << i +1 << " - " << j +1 << endl;
                    }


                    if((*A)(i, j-2).x >= 0) //meio e esquerda esquerda 4
                    {
                        coeffK.push_back(T((*A)(i, j).x, (*A)(i, j-2).x, (lambda+2*mu)/hx_2));
                    }
                    else
                    {
                        cout << "[REGR_X_PROG_Y] errou em ((*A)(i, j-2).x >= 0) " << i +1 << " - " << j +1 << endl;
                    }                    
                    if((*A)(i, j-2).y >= 0) //meio e esquerda esquerda 4
                    {
                        coeffK.push_back(T((*A)(i, j).y, (*A)(i, j-2).y, mu/hx_2));
                    }
                    else
                    {
                        cout << "[REGR_X_PROG_Y] errou em ((*A)(i, j-2).y >= 0) " << i +1 << " - " << j +1 << endl;
                    }


                    if((*A)(i-1, j-1).x >= 0) //diagonal cima e esquerda 5
                    {
                        coeffK.push_back(T((*A)(i, j).x, (*A)(i-1, j-1).x, -(lambda + mu)/(hx*hy)));
                    }
                    else
                    {
                        cout << "[REGR_X_PROG_Y] errou em ((*A)(i-1, j-1).x >= 0) " << i +1 << " - " << j +1 << endl;
                    }
                    if((*A)(i-1, j-1).y >= 0) //diagonal cima e esquerda 5
                    {
                        coeffK.push_back(T((*A)(i, j).y, (*A)(i-1, j-1).y, -(lambda + mu)/(hx*hy)));
                    }
                    else
                    {
                        cout << "[REGR_X_PROG_Y] errou em ((*A)(i-1, j-1).y >= 0) " << i +1 << " - " << j +1 << endl;
                    }

                    break;
                }
                case BORDER_TYPE::CENT_X_PROG_Y:
                {
                    /*
                           (2)
                            |
                       (5) (1) (6)
                         \  |  /
                       (3)-|0|-(4)
                    */
                    //pensando no problema que virá caso o ponto seja preso em X e solto em Y, é preciso verificar se pode colocar em y os dados


                    bool setYAxis = ((*A)(i, j).y != -1);
                    //adicionando as informações da direção X e verificando se pode incluir na direção Y
    //              cout << "Linha e coluna da matriz K  " <<(*A)(i, j).x<< endl;
                    coeffK.push_back(T((*A)(i, j).x, (*A)(i, j).x, -(2*hy_2*lambda+4*hy_2*mu-hx_2*mu)/(hx_2*hy_2)));
                    if(setYAxis) coeffK.push_back(T((*A)(i, j).y, (*A)(i, j).y, (hx_2*lambda-2*hy_2*mu+2*hx_2*mu)/(hx_2*hy_2)));

                    if((*A)(i-1, j).x >= 0) //cima e meio 1
                    {
                        coeffK.push_back(T((*A)(i, j).x, (*A)(i-1, j).x, -2*mu/hy_2));
                    }
                    else
                    {
                        cout << "[CENT_X_PROG_Y] errou em ((*A)(i-1, j).x >= 0) " << i +1 << " - " << j +1 << endl;
                    }
                    if( ((*A)(i-1, j).y >= 0) && setYAxis ) //cima e meio 1
                    {
                        coeffK.push_back(T((*A)(i, j).y, (*A)(i-1, j).y, -2*(lambda+2*mu)/hy_2));
                    }
                    else
                    {
                        cout << "[CENT_X_PROG_Y] errou em ((*A)(i-1, j).y >= 0) " << i +1 << " - " << j +1 << endl;
                    }

                    if((*A)(i-2, j).x >= 0) //cima cima e meio 2
                    {
                        coeffK.push_back(T((*A)(i, j).x, (*A)(i-2, j).x, mu/hy_2));
                    }
                    else
                    {
                        cout << "[CENT_X_PROG_Y] errou em ((*A)(i-2, j).x >= 0) " << i +1 << " - " << j +1 << endl;
                    }
                    if(((*A)(i-2, j).y >= 0) && setYAxis) //cima cima e meio 2
                    {
                        coeffK.push_back(T((*A)(i, j).y, (*A)(i-2, j).y, (lambda+2*mu)/hy_2));
                    }
                    else
                    {
                        cout << "[CENT_X_PROG_Y] errou em ((*A)(i-2, j).y >= 0) " << i +1 << " - " << j +1 << endl;
                    }

                    if((*A)(i, j-1).x >= 0) //meio e esquerda 3
                    {
                        coeffK.push_back(T((*A)(i, j).x, (*A)(i, j-1).x, (lambda+mu)/(2*hx*hy) + (lambda+2*mu)/hx_2 ));
                    }
                    else
                    {
                        cout << "[CENT_X_PROG_Y] errou em ((*A)(i, j-1).x >= 0) " << i +1 << " - " << j +1 << endl;
                    }
                    if(((*A)(i, j-1).y >= 0) && setYAxis )  //meio e esquerda 3
                    {
                        coeffK.push_back(T((*A)(i, j).y, (*A)(i, j-1).y, (lambda+mu)/(2*hx*hy) + mu/hx_2));
                    }
                    else
                    {
                        cout << "[CENT_X_PROG_Y] errou em ((*A)(i, j-1).y >= 0) " << i +1 << " - " << j +1 << endl;
                    }

                    if((*A)(i, j+1).x >= 0) //meio e direita 4
                    {
                        coeffK.push_back(T((*A)(i, j).x, (*A)(i, j+1).x, -(lambda+mu)/(2*hx*hy) + (lambda+2*mu)/hx_2));
                    }
                    else
                    {
                        cout << "[CENT_X_PROG_Y] errou em ((*A)(i, j+1).x >= 0) " << i +1 << " - " << j +1 << endl;
                    }
                    if(((*A)(i, j+1).y >= 0) && setYAxis )  //meio e direita 4
                    {
                        coeffK.push_back(T((*A)(i, j).y, (*A)(i, j+1).y, -(lambda+mu)/(2*hx*hy) + mu/hx_2));
                    }
                    else
                    {
                        cout << "[CENT_X_PROG_Y] errou em ((*A)(i, j+1).y >= 0) " << i +1 << " - " << j +1 << endl;
                    }

                    if((*A)(i-1, j-1).x >= 0) //diagonal cima e esquerda 5
                    {
                        coeffK.push_back(T((*A)(i, j).x, (*A)(i-1, j-1).x, -(lambda + mu)/(2*hx*hy)));
                    }
                    else
                    {
                        cout << "[CENT_X_PROG_Y] errou em ((*A)(i-1, j-1).x >= 0) " << i +1 << " - " << j +1 << endl;
                    }
                    if(((*A)(i-1, j-1).y >= 0) && setYAxis )  //diagonal cima e esquerda 5
                    {
                        coeffK.push_back(T((*A)(i, j).y, (*A)(i-1, j-1).y, -(lambda + mu)/(2*hx*hy)));
                    }
                    else
                    {
                        cout << "[CENT_X_PROG_Y] errou em ((*A)(i-1, j-1).y >= 0) " << i +1 << " - " << j +1 << endl;
                    }

                    if((*A)(i-1, j+1).x >= 0) //diagonal cima e direita 6
                    {
                        coeffK.push_back(T((*A)(i, j).x, (*A)(i-1, j+1).x, (lambda + mu)/(2*hx*hy)));
                    }
                    else
                    {
                        cout << "[CENT_X_PROG_Y] errou em ((*A)(i-1, j+1).x >= 0) " << i +1 << " - " << j +1 << endl;
                    }
                    if(((*A)(i-1, j+1).y >= 0 ) && setYAxis ) //diagonal cima e direita 6
                    {
                        coeffK.push_back(T((*A)(i, j).y, (*A)(i-1, j+1).y, (lambda + mu)/(2*hx*hy)));
                    }
                    else
                    {
                        cout << "[CENT_X_PROG_Y] errou em ((*A)(i-1, j+1).y >= 0) " << i +1 << " - " << j +1 << endl;
                    }

                    break;
                }
                default:
    //                cout << "Não tem tipo" << endl;
                    break;
            }



        }
    }


//    cout << " hx = " << hx << " ;  hy = " << hy << endl;
//    cout << " lambda = " << lambda << " ;  mu = " << mu << endl << endl;
    for(unsigned int i = 0; i< coeffK.size(); i++)
    {
        if(coeffK[i].row()==-1 || coeffK[i].col()==-1)
            cout << "K[" << coeffK[i].row() << ", "<< coeffK[i].col() << "] = " << coeffK[i].value() << endl;
    }


    K.setFromTriplets(coeffK.begin(), coeffK.end());
    coeffK.clear();
    K.makeCompressed();
    preenchendoSistemaTime.endCount();
    preenchendoSistemaTime.print();

//    cout << "\n\nK " <<endl;
//    cout << MatrixXd(K)<<endl;

//    cout << "\n\nf"<<endl;
//    cout << MatrixXd(f)<<endl;

//    cout << "Determinante de K " << MatrixXd(K).determinant() << endl;

    cout << "Resolvendo o sistema" << endl;
    TimeSpent resolvendoSistemaTime("Resolvido o sistema | TEMPO: ");
    resolvendoSistemaTime.startCount();
    //-------------------------------------------------------------------
    //resolvendo o sistema por cholesky
//    SimplicialCholesky<SpMatrixD> chol(K);
//    VectorXd x = chol.solve(f);
//    resolvendoSistemaTime.endCount();
//    resolvendoSistemaTime.print();

    //resolvendo o sistema por LU
    SparseLU<SpMatrixD> lu;
    K.makeCompressed();
    lu.analyzePattern(K);
    lu.factorize(K);
    if(lu.info()!=Success)
    {
        std::cout<<"\t\t\t\tnão conseguiu fatorar a matriz de !!!"<<std::endl;
        exit(1);
    }
    cout << "!Determinante " << lu.logAbsDeterminant() << endl;
    VectorXd x = lu.solve(-f);

//    bool isSymmetric = false;
//    lu.isSymmetric(isSymmetric);
//    if (isSymmetric) cout << "é simetrico " << endl;
//    else cout << "NAO é simetrico " << endl;
//    cout <<"\n\n"<< MatrixXd(f)<<endl;
    cout <<"\n\n"<< MatrixXd(x)<<endl;

    //--------------------------------------------------------------------
//    cout << "Verificadndo o condicionamento da matriz" <<endl;
//    VectorXcd eivals = MatrixXd(K).eigenvalues();
////    cout << "autovalores \n" << eivals <<endl;

//    std::vector<double> myvector;
//    for(auto i = 0; i < eivals.rows(); i++)
//        myvector.push_back(eivals[i].real());
//    cout << "Adicionados e indo ordenar"<< endl;

//    std::sort (myvector.begin(), myvector.end());
////    for(auto i = 0; i < eivals.rows(); i++)
////        cout << myvector[i] << " " << endl;
//    cout << "Primeiro autovalor " << (myvector.front()) << endl;
//    cout << "Ultimo autovalor " << myvector.back() << endl;

//    cout <<"Numero de condicionamento:  "<< fabs(myvector.back())/fabs(myvector.front()) << endl;

    //--------------------------------------------------------------------




    cout <<"acabou processo"<< endl;

    cout << "Salvando o resultado" << endl;
    TimeSpent salvandoTime("Salvo o resultado | TEMPO: ");
    salvandoTime.startCount();
    saveDisplacementsHaunch(label, L, H, hx, hy, n_elemL, n_elemH, A, x);
    salvandoTime.endCount();
    salvandoTime.print();

    delete A;

    cout <<"Fim"<< endl;
}

int main(int argc, char *argv[])
{
    unsigned int n_elem = 24; //numero de elementos em cada direção (SEMPRE PAR) 6 12 24 48
    double L = 2; //tamanho do domínio quadrado
    double p = 0;   //constante para a função p
    P typeP = P::PATCH_TEST;
    if(argc>1)
        n_elem = atoi(argv[1]);
    quadModel(n_elem, L, p, typeP);



//    unsigned int n_elem =10; //numero de elementos em cada direção (SEMPRE PAR)
//    double L = M_PI; //tamanho do domínio quadrado
//    double p = 1;   //constante para a função p
//    P typeP = P::SINSIN;
//    if(argc>1)
//        n_elem = atoi(argv[1]);
//    quadModel(n_elem, L, p, typeP);





//    unsigned int n_elemL = 2000; //numero de elementos em cada direção X
//    unsigned int n_elemD = 2000; //numero de elementos em cada direção Y
//    double L = 2; //tamanho do domínio quadrado
//    double D = 2; //tamanho do domínio quadrado
//    double p = 1;   //constante para a função p
//    P typeP = P::POISSON_EXTRA;
//    if(argc>1)
//        n_elemL = n_elemD= atoi(argv[1]);
    //ladderModel(n_elemL, n_elemD, L, D, p, typeP);




//    double L = 12; //largura total (m)
//    double l_haunch = 0.5; //tamanho da misula
//    double l_base = 0.1; //tamanho da base reta da misula
//    double H = 0.6; //altura total
//    double h_haunch = 0.2; //tamanho da misula
//    double b = 0.3; //thickness
//    double p = -5000;   //constante para a função p
//    double E = 2.21*pow(10, 10); //concreto
//    double nu = 0.15;
//    P typeP = P::CONSTANT;
//    unsigned int n_xb = 3; //numero de nós na base (número impar!)


//    haunchedBeamModel(L, l_haunch, l_base, H,
//                      h_haunch, b, E, nu, p, typeP, n_xb);

    return 0;
}
