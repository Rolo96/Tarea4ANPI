//
// Created by bryan on 18/03/18.
//

#ifndef TAREA03_QR_H
#define TAREA03_QR_H

#include <vector>
#include "Matrix.hpp"
#include "Exception.hpp"
#include "AuxMethods.h"

namespace anpi{

    /**
 * Metodo que aplica el algoritmo de Householder para encontrar
 * las matrices Q y R de una matriz dada
 * @param A Matriz a la cual se le hace el analisis
 * @param Q Matriz Q resultante
 * @param R Matriz R resultante (triangular superior)
 */
    template <typename T>
    void qr(const Matrix<T>& A, Matrix<T>& Q,Matrix<T>& R ){
        R = A;
        anpi::Matrix<T> Qtmp(A.rows(),A.rows(),T(0));
        for(size_t j = A.rows(); j > 1; --j){
            anpi::Matrix<T> I(j,j,T(0));
            std::vector<T> e(j,0); e[0] = 1;
            std::vector<T> u(j);
            std::vector<T> a(j);
            anpi::Matrix<T> uu(u.size(),u.size(),anpi::DoNotInitialize);
            T unorma;

            for(size_t fila = 0; fila < j; ++fila){ //iteraciones para crear la matriz identidad y formar el vector a
                I(fila, fila) = 1; //diagonal igual a 1
                a[fila] = R(fila+abs(j-A.rows()),abs(j-A.rows()));//se crea el vector a desde el valor R(j,1), en la primera iteracion R=A
            }
            u = sign(a[0])*(norm(a)*e) - a;//se encuentra el vector u

            for(size_t i = 0; i<u.size(); ++i){ //iteraciones para crear la matriz uu
                for (size_t j = 0; j < u.size(); ++j) {
                    uu(i,j) = u[i]*u[j]; //multiplicacion de (u)*(u transpuesto)
                }
            }
            unorma = norm(u); //se encuentra la norma de u

            for(size_t i = 0; i < I.rows(); ++i) {//iteraciones para encontrar Qi
                for(size_t jj = 0; jj< I.cols(); ++jj) {
                    Qtmp(i+abs(j-A.rows()),jj+abs(j-A.rows())) = I(i,jj) - (T(2)/(unorma*unorma))*(uu(i,jj));//aplicando lo indicado en el metodo Householder
                }
            }
            R = Qtmp*R; //modificacion de R (cada Qi por la matriz R anterior, la primera R es igual a A)
            if (j != A.rows())//se pregunta si estamos en la primera iteracion
                Q = Q*(trans(Qtmp));//si no se esta entonces Q toma su valor anterior por el nuevo Qi tranpuesta
            else
                Q = trans(Qtmp); //si es la primera iteracion Q es igual a Qi transpuesta
            Qtmp.fill(T(0)); //se ponen todos los valores del siguiente Qi en 0
            for(size_t fila = 0; fila < Qtmp.rows(); ++fila) { //se hace la diagonal del Qi siguiente en 1
                Qtmp(fila, fila) = 1;
            }
        }
    }

    /**
 * Este metodo resuelve un sistema de ecuaciones utilizando la descomposicion QR
 * @param A Matriz que contiene los coeficientes del sistema de ecuaciones
 * @param x Vector de las incognitas que seran determinadas con este metodo
 * @param b Vector que contiene los resultados de cada ecuacion del sistema
 */
    template <typename T>
    void solveQR(const Matrix<T>& A, std::vector<T>& x, const std::vector<T> b){
        if (x.size() != A.rows()) { throw anpi::Exception("Wrong size of vector x"); }
        if (b.size() != A.rows()) { throw anpi::Exception("Wrong size of vector b"); }
        if (singular(A)) { throw anpi::Exception("A is a singular matrix"); }
        anpi::Matrix<T> Q(A.rows(), A.cols(), DoNotInitialize);
        anpi::Matrix<T> R(A.rows(), A.cols(), DoNotInitialize);
        qr(A, Q, R);
        std::vector<T> Qb(Q.rows());
        Qb = trans(Q) * b;
        T restador = T();
        x[x.size() - 1] = Qb[Qb.size() - 1] /
                          R(R.rows() - 1, R.cols() - 1);
        for (int i = x.size() - 2; i > -1; --i) {
            for (int j = i + 1; j < x.size(); ++j) {
                restador += R(i, j) * x[j];
            }
            x[i] = (Qb[i] - restador) / R(i, i);
            restador = T();
        }
    }
}

#endif //TAREA03_QR_H
