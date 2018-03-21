//
// Created by bryan on 18/03/18.
//

#ifndef TAREA03_AUXMETHODS_H
#define TAREA03_AUXMETHODS_H

#include <vector>
#include <cmath>
#include "Matrix.hpp"
#include "Exception.hpp"

namespace anpi{

    template<typename T>
    void pivot(anpi::Matrix<T>& a,anpi::Matrix<T>& lu,int diag,std::vector<size_t>& permut){
        T bigger=fabs(a(diag,diag));
        int pivotRow=diag;
        for(int i=diag+1;i<a.rows();++i){
            if(fabs(a(i,diag))>bigger){bigger=fabs(a(i,diag));pivotRow=i;}
        }
        if(pivotRow!=diag){
            std::swap(permut[pivotRow],permut[diag]);
            for(int i=0;i<a.cols();++i){
                std::swap(a(diag,i),a(pivotRow,i));
                std::swap(lu(diag,i),lu(pivotRow,i));
            }
        }
    }

    /**
 * Metodo para sobrecargar el operador "-" (resta) para dos vectores dados
 * @param a primer vector
 * @param b segundo vector
 * @return resultado de restar el vector a - vector b
 */
    template<typename T>
    std::vector<T> operator-(const std::vector<T>& a, const std::vector<T>& b) {
        if(a.size() != b.size()){ throw anpi::Exception("Vectors with different sizes");}
        std::vector<T> c(a.size(), T(0));
        for (unsigned long i = 0; i < a.size(); ++i) {
            c.at(i) = a.at(i) - b.at(i);
        }
        return c;
    }

    /**
 * Metodo para sobrecargar el operador "*" (multiplicacion) de un escalar por un vector (debe hacerse en ese orden)
 * @param cons escalar para la multiplicacion
 * @param a vector que sera multiplicado por el escalar
 * @return resultado de multiplicar el escalar * el vector
 */
    template<class T>
    std::vector<T> operator*(const T cons, const std::vector<T>& a) {
        std::vector<T> c(a.size(),T(0));
        for(unsigned long i=0;i<a.size();++i){
            c.at(i)=a.at(i)*cons;
        }
        return c;
    }


    template<typename T>
    Matrix<T> trans(Matrix<T> A) {
        for(unsigned long i=0;i<A.rows();++i){
            for(unsigned long j=0;j<A.cols();++j){
                std::swap(A(i,j),A(j,i));
            }
        }
        return A;
    }

    /**
 * Metodo utilizado para encontrar la norma de un vector dado
 * @param a vector al cual se le va a encontrar su norma
 * @return norma del vector
 */
    template<typename T>
    T norm(std::vector<T> a){
        T norma=T();
        for(int i=0;i<a.size();++i){
            norma+=a[i]*a[i];
        }
        return sqrt(norma);
    }
/**
 * Metodo para obtener el sigo de un numero dado
 * @param numero numero para determinarle el signo
 * @return devulve 1 si numero>=0 y -1 si numero<0
 */
    template<class T>
    T sign(T numero){
        if(numero<0)return -1;
        return 1;
    }

    /**
 * Metodo que obtiene el determinante de una matriz e indica si es igual a 0
 * @param A Matriz a la que se le obtiene el determinante
 * @return 1(true) si el determiante es 0, caso contrario 0(false)
 */
    template <typename T>
    bool singular(const anpi::Matrix<T>& A){
        anpi::Matrix<T> Atmp(A); //Como la matriz se modifica se crea una temporal para mantener los valores de A
        T det = T(); //variable que almacena el valor del determinante
        det = Atmp(0,0); // se toma el primer termino de la matriz como determinante inicial
        for(int k = 0; k < Atmp.rows()-1; ++k){
            for (int i = k+1; i < Atmp.rows() ; ++i) {
                for (int j = k+1; j < A.rows(); ++j) {
                    if(Atmp(k,k) != T(0))//se valida que no suceda division entre cero
                        Atmp(i,j) = (Atmp(k,k)*Atmp(i,j) - Atmp(k,j)*Atmp(i,k))/Atmp(k,k);//se modifica la matriz Atmp con nuevos valores
                    else
                        break;//si existe entonces se sale de la iteracion
                }
            }
            det = det*Atmp(k+1,k+1);//el determinante toma el valor anterior multiplicado por un nuevo valor de Atmp
        }
        if(det == T(0))//si el determinante es 0
            return true;//retorna un 1
        return false;//si no, retorna un 0
    }
}

#endif //TAREA03_AUXMETHODS_H
