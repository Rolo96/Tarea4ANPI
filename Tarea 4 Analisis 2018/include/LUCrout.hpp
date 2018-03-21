/**
 * Copyright (C) 2018
 * Área Académica de Ingeniería en Computadoras, ITCR, Costa Rica
 *
 * This file is part of the numerical analysis lecture CE3102 at TEC
 *
 * @Author: 
 * @Date  : 03.03.2018
 */

#include <cmath>
#include <limits>
#include <functional>

#include "Exception.hpp"
#include "Matrix.hpp"

#ifndef ANPI_LU_CROUT_HPP
#define ANPI_LU_CROUT_HPP

namespace anpi {

  /**
   * Auxiliary method used to debug LU decomposition.
   *
   * It separates a packed LU matrix into the lower triangular matrix
   * L and the upper triangular matrix U, such that the diagonal of U
   * is composed by 1's.
   */
  template<typename T>
  void unpackCrout(const Matrix<T>& LU,
                   Matrix<T>& L,
                   Matrix<T>& U) {
    if(LU.rows()!=LU.cols()){throw anpi::Exception("LU is not a square matrix");}
    if(L.rows()!=L.cols()){throw anpi::Exception("L is not a square matrix");}
    if(U.rows()!=U.cols()){throw anpi::Exception("U is not a square matrix");}
    if(((LU.rows()!=L.rows())&&(LU.rows()!=U.rows()))||
       ((LU.cols()!=L.cols())&&(LU.cols()!=U.cols())))
    {throw anpi::Exception("LU,L and U don't have the same size");}

    for(unsigned int i=0;i<LU.rows();++i){
      for (unsigned int j = 0; j < LU.cols(); ++j) {
        if(i<j){///When it is in the upper triangular matrix of LU
          L(i,j)=0;
          U(i,j)=LU(i,j);
        }else if(i>j){///When it is in the lower triangular matrix LU
          U(i,j)=0;
          L(i,j)=LU(i,j);
        }else{///When it is on the diagonal of LU
          L(i,j)=LU(i,j);
          U(i,j)=T(1);
        }
      }
    }
  }
  
  /**
   * Decompose the matrix A into a lower triangular matrix L and an
   * upper triangular matrix U.  The matrices L and U are packed into
   * a single matrix LU.  
   *
   * Crout's way of packing assumes a diagonal of
   * 1's in the U matrix.
   *
   * @param[in] A a square matrix 
   * @param[out] LU matrix encoding the L and U matrices
   * @param[out] permut permutation vector, holding the indices of the
   *             original matrix falling into the corresponding element.
   *             For example if permut[5]==3 holds, then the fifth row
   *             of the LU decomposition in fact is dealing with the third
   *             row of the original matrix.
   *
   * @throws anpi::Exception if matrix cannot be decomposed, or input
   *         matrix is not square.
   */
  template<typename T>
  void luCrout(const Matrix<T>& A,
               Matrix<T>& LU,
               std::vector<size_t>& permut) {

    if(A.rows()!=A.cols()){throw anpi::Exception("A is not a square matrix");}
    if(LU.rows()!=LU.cols()){throw anpi::Exception("LU is'n a square matrix");}
    if((A.rows()!=LU.rows()) || (A.cols()!=LU.cols())){throw anpi::Exception("A and LU are not the same size");}
    if((permut.size()!=A.rows())){throw anpi::Exception("Invalid size of the permut vector");}
    ///The copy is made because it has to be edited
    Matrix<T> ACopy=A;
    pivot(ACopy,LU,0,permut);
    T l0 = LU(0, 0) = ACopy(0, 0);
    T n = ACopy.rows();
    for (int i = 1; i < n; ++i) {
      LU(i, 0) = ACopy(i, 0);
      LU(0, i) = ACopy(0, i) / l0;
    }
    T summation = T(0);
    for (int j = 1; j < n - 1; ++j) {
      pivot(ACopy,LU,j,permut);
      for (int i = j; i < n; ++i) {
        for (int k = 0;
             k < j; ++k) {
          summation += LU(i, k) * LU(k, j);
        }
        LU(i, j) = ACopy(i, j) - summation;
        summation = T(0);
      }
      summation = T(0);
      for (int k = j + 1; k < n; ++k) {
        for (int i = 0;
             i < j; ++i) {
          summation += LU(j, i) * LU(i, k);
        }
        LU(j, k) = (ACopy(j, k) - summation) / LU(j, j);
        summation = T(0);
      }
      summation = T(0);
    }
    summation = T(0);
    for (int k = 0; k < n; ++k) {
      summation += LU(n - 1, k) * LU(k, n - 1);
    }
    LU(n - 1, n - 1) = ACopy(n - 1, n - 1) - summation;
  }
}
  
#endif

