//
// Created by bryan on 18/03/18.
//

#ifndef TAREA03_LU_H
#define TAREA03_LU_H

#include <vector>
#include <iostream>
#include "Matrix.hpp"
#include "Exception.hpp"
#include "QR.h"
#include "LUDoolittle.hpp"
#include "LUCrout.hpp"

namespace anpi {

/**
 * Este metodo resuelve un sistema de ecuaciones utilizando la descomposicion LU
 * @param A Matriz que contiene los coeficientes del sistema de ecuaciones
 * @param x Vector que contiene las incognitas del sistema de ecuaciones
 * @param b Vector que contiene los resultados del sistema de ecuaciones
 */
    template<typename T>
    void solveLU(Matrix<T> A, std::vector <T> &x, const std::vector <T> &b) {
        if (x.size() != A.rows()) { throw anpi::Exception("Wrong size of vector x"); }
        if (b.size() != A.rows()) { throw anpi::Exception("Wrong size of vector b"); }
        if (singular(A)) { throw anpi::Exception("A is a singular matrix"); }
        anpi::Matrix<T> Lu(A.rows(), A.cols(), DoNotInitialize);
        std::vector<size_t > permut(A.rows());
        for(int i=0;i<A.rows();++i){
            permut[i]=size_t(i);
        }
        luCrout(A, Lu,permut);
        for(int i=0;i<A.rows();++i){
            std::swap(b[i],b[permut[i]]);
        }
        std::vector <T> xtmp(A.rows(), T(0));
        std::vector <T> y(A.rows(), T(0)); //vector y que se obtiene de Ly=b donde y = Ux
        y[0] = b[0] / Lu(0, 0);//primer termino yi que se obtiene por la sustitucion hacia adelante
        T restador = T(); //factor que se obtiene de la sumatoria de Lijyj el cual se usa en yi=(1/Lii)(bi-restador)
        for (int i = 1; i < y.size(); ++i) {
            for (int j = 0; j < i; ++j) {
                restador += Lu(i, j) * y[j];//se obtiene el restador de Lijyj
            }
            y[i] = (b[i] - restador) / Lu(i, i);//se obtiene el yi de la iteracion actual
            restador = T();//se reinicia el restador a 0
        }
        for (int i = 0; i < Lu.rows(); ++i) {
            Lu(i, i) = T(
                    1);//como las matrices L y U se encuentran en una misma matriz LU, se hace la identidad de esta en 1 para utilizar a U
        }
        x[x.size() - 1] = y[y.size() - 1] /
                          Lu(Lu.rows() - 1, Lu.cols() - 1);//ultimo xi que se obtiene por la sustitucion hacia atras
        for (int i = x.size() - 2; i > -1; --i) {
            for (int j = i + 1; j < x.size(); ++j) {
                restador += Lu(i, j) * x[j];//este restador es para la ecuacion xi = (1/Uii)(yi-restador)
            }
            x[i] = (y[i] - restador) / Lu(i, i);//se obtiene el xi de la iteracion actual
            restador = T();//se reinicia el restador para la siguiente iteracion
        }
    }


    /**
 * Metodo que calcula la inversa de una matriz utilizando la descomposicion LU
 * @param A: Matriz original
 * @param Ai: Inversa de A
 */
    template<typename T>
    void invert(const Matrix<T>& A, Matrix<T>& Ai){
        if (singular(A)) { throw anpi::Exception("A is a singular matrix"); }
        for (int j = 0; j < A.cols(); ++j) {
            std::vector<T> b(A.rows(), T(0));
            std::vector<T> x(A.rows(), T(0));
            b[j] = 1;
            //std::cout<<x[0]<<"  "<<x[1]<<"  "<<x[2]<<std::endl<<std::endl;
            solveLU(A, x, b);
            //std::cout<<x[0]<<"  "<<x[1]<<"  "<<x[2]<<std::endl<<std::endl;
            for (int i = 0; i < A.rows(); ++i) {
                Ai(i, j) = x[i];
            }
        }
    }
}
#endif //TAREA03_LU_H
