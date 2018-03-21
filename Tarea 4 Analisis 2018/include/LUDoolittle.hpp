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
#include <algorithm>

#include "Exception.hpp"
#include "Matrix.hpp"

#ifndef ANPI_LU_DOOLITTLE_HPP
#define ANPI_LU_DOOLITTLE_HPP

namespace anpi {

  /**
   * Auxiliary method used to debug LU decomposition.
   *
   * It separates a packed LU matrix into the lower triangular matrix
   * L and the upper triangular matrix U, such that the diagonal of L
   * is composed by 1's.
   */
  template<typename T>
  void unpackDoolittle(const Matrix<T>& LU,
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
          L(i,j)=T(1);
          U(i,j)=LU(i,j);
        }
      }
    }
  }
  
  /**
   * Decompose the matrix A into a lower triangular matrix L and an
   * upper triangular matrix U.  The matrices L and U are packed into
   * a single matrix LU. 
   *
   * The L matrix will have in the Doolittle's LU decomposition a
   * diagonal of 1's
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
  void luDoolittle(const Matrix<T>& A,
                   Matrix<T>& LU,
                   std::vector<size_t>& permut) {
    if(A.rows()!=A.cols()){throw anpi::Exception("A is not a square matrix");}
    if(LU.rows()!=LU.cols()){throw anpi::Exception("LU is'n a square matrix");}
    if((A.rows()!=LU.rows()) || (A.cols()!=LU.cols())){throw anpi::Exception("A and LU are not the same size");}
    if((permut.size()!=A.rows())){throw anpi::Exception("Invalid size of the permut vector");}
    ///The copy is made because it has to be edited
    Matrix<T> ACopy=A;
    ///The next three cycles, make the pivot,
    /// then the one with j as the iteration
    /// variable goes row by row, and the one
    /// with k as the iteration variable goes
    /// colum by colum computing the elimination
    /// factor and doing the elimination process
    /// while the factors and the pivot row are
    /// stored in LU
    for(unsigned int i=0;i<ACopy.rows()-1;++i){
      pivot(ACopy,LU,i,permut);
      for(unsigned int j=i+1;j<ACopy.rows();++j){
        T factor = ACopy(j,i)/ACopy(i,i);
        for(unsigned int k=i;k<ACopy.cols();++k){
          LU(i,k)=ACopy(i,k);
          if(j>k) {
            LU(j, k) = ACopy(j,k)/ACopy(k,k);
          }
          ACopy(j,k)=ACopy(j,k)-factor*ACopy(i,k);
        }
      }
    }
    LU(LU.rows()-1,LU.cols()-1)=ACopy(ACopy.rows()-1,ACopy.cols()-1);
  }
}
  
#endif

