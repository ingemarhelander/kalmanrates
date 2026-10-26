#include <iomanip>
#include <fstream>
#include <math.h>

#include "boost/numeric/ublas/matrix.hpp"
#include "boost/numeric/ublas/io.hpp"
#include "boost/lexical_cast.hpp"
//#include "XiMath/boostHelper.h"
#include "XiMath/sparseMatrix.h"
//#include "XiMath/iomatlab.h"

using namespace boost::numeric::ublas;
using namespace XiMath;
using namespace XiFin;

#ifndef KALMANLOGLIKELIHOODAD_H
#define KALMANLOGLIKELIHOODAD_H

class ExcKal
{
public:
	ExcKal(const std::string& msg) : msg_(msg) {}
	virtual void print() const {std::cerr << msg_ << std::endl;}
	std::string msg_;
};

template <class T>
class AutoDiff
{
public:
  AutoDiff() {}
  virtual ~AutoDiff() {}

  virtual void count(std::size_t& rMax, std::size_t& dMax) const {}
  virtual void fAD(vector<T>& d, vector<std::size_t>& col, vector<std::size_t>& row, vector<std::size_t>& rInd, std::size_t& nr, std::size_t& nd) {}

  T testAD()
  { // Compute the norm of the difference of Jacobian between automatic differentiation and numeric differentiation
    std::size_t nIn;
    std::size_t nOut;
    countSize(nIn, nOut);

    matrix<T> Jacobian(nIn, nOut);
    std::size_t nc=0;
    AutoDiff<T>* ad = createLocalObj(nc); // Creates local copies of Ind vectors and matrices
//std::cout << "Local obj ok " << ad << std::endl;
    if (ad == NULL)
      throw ExcKal(std::string("AutoDiff: Object not created "));

    // Automatic differentiation
    for (std::size_t j=0 ; j<nOut ; j++)
    {
      std::size_t rMax = 0;
      std::size_t dMax = 0;
      ad->count(rMax, dMax);
//std::cout << "count ok" << std::endl;

      // Initialize computations for automatic differentiation

      vector<T> d(dMax, 0.0);
      vector<std::size_t> col(dMax);
      vector<std::size_t> row(rMax+1);
      vector<std::size_t> rInd(rMax);
      std::size_t nr=0;
      std::size_t nd=0;

//std::cout << "before fAD" << std::endl;
      ad->fAD(d, col, row, rInd, nr, nd);
//std::cout << "fAD ok" << std::endl;

      row(nr) = nd; // Final column

      if (nr > rMax)
        std::cout << "AD: Incorrect number of rows:" << nr << ">" << rMax << std::endl;
      if (nd > dMax)
        std::cout << "AD: Incorrect number of data:" << nd << ">" << dMax << std::endl;

      vector<T> dOdx(nc, 0.0);

      dOdx(nIn+j) = 1;

      for (std::size_t i=nr ; i-- > 0 ; )
      {
        const T v1 = dOdx(rInd(i));
        if (v1 != 0)
        {
    //std::cerr << i << " " << rInd(i) << " : " << row(i) << " ->" << row(i+1) << " mult " << v1 << std::endl; 
          for (std::size_t j=row(i) ; j<row(i+1) ; j++)
          {  
            T& v2 = dOdx(col(j));
    //std::cerr << xn(rInd(i)) << " " << rInd(i) << " : "  << xn(col(j)) << " = " << v2 << " + " << v1*d(j) << " = " << v2 + v1*d(j) << std::endl; 
            if (v2==0)
              v2 = v1*d(j);
            else
              v2 += v1*d(j);
          }
        }
      }
//std::cout << "AD ok" << std::endl;

      vector<T> grad(nIn);
      ad->testADgrad(dOdx, grad);
      for (std::size_t i=0 ; i<nIn ; i++)
        Jacobian(i,j) = grad(i);
    }

    ad->freeLocalObj(); // Deletes local copies of Ind vectors and matrices
    delete ad;
//std::cout << "free ok" << std::endl;

    matrix<T> JacobianNum(nIn, nOut);
    const T epsilon = 1E-6;
    gradNum(epsilon, JacobianNum);
//std::cout << "num ok" << std::endl;

    T norm = 0;
    for (std::size_t i=0 ; i<nIn ; i++)
      for (std::size_t j=0 ; j<nOut ; j++)
      {
        T diff = abs(Jacobian(i,j) - JacobianNum(i,j));
        if (diff > norm)
          norm = diff;
        //if (abs(diff) > 1E-3)
        //  std::cout << i << " " << j << " : " << Jacobian(i,j) << " " << JacobianNum(i,j) << std::endl;
      }
    //std::cout << "norm = " << norm << std::endl;
    //std::ofstream of("Jac.m");
    //write2matlab(of, "Jacobian", Jacobian);
    //std::ofstream ofNum("JacNum.m");
    //write2matlab(ofNum, "JacobianNum", JacobianNum);

    return norm;
  }
protected:
  template<class S>
  inline void checkDim(const std::string& strFunc, const std::string& strVec, const vector<S>& a, const std::size_t m) const
  {
    if (a.size() != m)
      throw ExcKal(std::string("AutoDiff:") + strFunc + std::string(" : Incorrect number of elements in ") + strVec);
  }
  template<class S>
  inline void checkDim(const std::string& strFunc, const std::string& strMat, const matrix<S>& A, const std::size_t m, const std::size_t n) const
  {
    if (A.size1() != m)
      throw ExcKal(std::string("AutoDiff:") + strFunc + std::string(" : Incorrect number of rows in ") + strMat);
    if (A.size2() != n)
      throw ExcKal(std::string("AutoDiff:") + strFunc + std::string(" : Incorrect number of columns in ") + strMat);
  }
  template<class S>
  inline void checkDim(const std::string& strFunc, const std::string& strMat, const SparseMatrix<S>& A, const std::size_t m, const std::size_t n) const
  {
    if (A.size1() != m)
      throw ExcKal(std::string("AutoDiff:") + strFunc + std::string(" : Incorrect number of rows in ") + strMat);
    if (A.size2() != n)
      throw ExcKal(std::string("AutoDiff:") + strFunc + std::string(" : Incorrect number of columns in ") + strMat);
  }

  virtual void countSize(std::size_t& nIn, std::size_t& nOut)
  {
    nIn = 0;
    nOut = 0;
  }
  virtual AutoDiff<T>* createLocalObj(std::size_t& nc)
  {
    nc = 0;
    return NULL;
  }
  virtual void freeLocalObj()
  { // Free objects created by createLocalObj
  }
  virtual void testADgrad(const vector<T>& dOdx, vector<T>& grad) // Compute gradient with respect to input arguments
  { }

  virtual void gradNum(const T epsilon, matrix<T>& Jacobian) {}

};

template <class T>
class AD_a_cumsumBx : public AutoDiff<T>
{ // y(x) = cumsum(Bx) / a
public:
  AD_a_cumsumBx(const T _a, const matrix<T>* Bxp, matrix<T>* Yxp, 
                const matrix<std::size_t>* BxIndp, const matrix<std::size_t>* YxIndp) 
  : a(_a), Bx(*Bxp), Yx(*Yxp), BxInd(*BxIndp), YxInd(*YxIndp) 
  {
    // Check dimensions
    std::size_t m = Bx.size1();
    std::size_t n = Bx.size2();
    std::string strFunc("AD_a_cumsumBx");
    checkDim(strFunc, "Bx", Bx, m, n);
    checkDim(strFunc, "Yx", Yx, m+1, n);
    checkDim(strFunc, "BxInd", BxInd, m, n);
    checkDim(strFunc, "YxInd", YxInd, m+1, n);
  }
  virtual ~AD_a_cumsumBx() {}

  virtual void count(std::size_t& rMax, std::size_t& dMax) const 
  {
    std::size_t m = Bx.size1();
    std::size_t n = Bx.size2();
    rMax += m*n; 
    dMax += 2*m*n; 
  }
  virtual void f() 
  {
    std::size_t m = Bx.size1();
    std::size_t n = Bx.size2();

    for (std::size_t j=0 ; j<n ; j++)
      Yx(0,j) = 0.0;

    for (std::size_t i=0 ; i<m ; i++)
      for (std::size_t j=0 ; j<n ; j++)
        Yx(i+1,j) = Yx(i,j) + Bx(i,j)/a;
  }
  virtual void fAD(vector<T>& d, vector<std::size_t>& col, vector<std::size_t>& row, vector<std::size_t>& rInd, std::size_t& nr, std::size_t& nd) 
  {
    std::size_t m = Bx.size1();
    std::size_t n = Bx.size2();

    for (std::size_t j=0 ; j<n ; j++)
      Yx(0,j) = 0.0;

    for (std::size_t i=0 ; i<m ; i++)
      for (std::size_t j=0 ; j<n ; j++)
      {
        rInd(nr) = YxInd(i+1,j);
        row(nr++) = nd;
        Yx(i+1,j) = Yx(i,j) + Bx(i,j)/a;
        d(nd) = 1;                                  // d Yx(i+1,j) / d Yx(i,j) 
        col(nd++) = YxInd(i,j);
        d(nd) = 1/a;                            // d Yx(i+1,j) / d Bx(i,j) 
        col(nd++) = BxInd(i,j);
      }
  }
protected:
  virtual void countSize(std::size_t& nIn, std::size_t& nOut)
  {
    nIn = Bx.size1() * Bx.size2();
    nOut = Yx.size1() * Yx.size2();
  }
  virtual AutoDiff<T>* createLocalObj(std::size_t& nc)
  { // Deallocated in freeLocalObj
    matrix<std::size_t>* BxInd_= new matrix<std::size_t>(Bx.size1(), Bx.size2()); 
    matrix<std::size_t>* YxInd_= new matrix<std::size_t>(Yx.size1(), Yx.size2()); 
    nc = 0;
    for (std::size_t i=0 ; i<Bx.size1() ; i++)
      for (std::size_t j=0 ; j<Bx.size2() ; j++)
        (*BxInd_)(i,j) = nc++;
    for (std::size_t i=0 ; i<Yx.size1() ; i++)
      for (std::size_t j=0 ; j<Yx.size2() ; j++)
        (*YxInd_)(i,j) = nc++;
    std::size_t nIn;
    std::size_t nOut;
    countSize(nIn, nOut);
    if (nc != nIn+nOut)
      throw ExcKal(std::string("AutoDiff: : Incorrect size in createLocalObj"));

    return new AD_a_cumsumBx(a, &Bx, &Yx, BxInd_, YxInd_);
  }

  virtual void freeLocalObj()
  { // Free objects created by createLocalObj
    delete &BxInd;
    delete &YxInd;
  }

  virtual void testADgrad(const vector<T>& dOdx, vector<T>& grad) // Compute gradient with respect to input arguments
  { 
    std::size_t k=0;
    for (std::size_t i=0 ; i<Bx.size1() ; i++)
      for (std::size_t j=0 ; j<Bx.size2() ; j++)
        grad(k++) = dOdx(BxInd(i,j));
  }

  virtual void gradNum(const T epsilon, matrix<T>& Jacobian)
  {
    std::size_t nIn;
    std::size_t nOut;
    countSize(nIn, nOut);
    Jacobian.resize(nIn, nOut);
    f();
    matrix<T> Yx_(Yx);
    matrix<T> Bx_(Bx);

    std::size_t ii=0;
    for (std::size_t i=0 ; i<Bx.size1() ; i++)
      for (std::size_t j=0 ; j<Bx.size2() ; j++)
      {        
        Bx_(i,j) += epsilon;
        AD_a_cumsumBx<T> ad(a, &Bx_, &Yx_, &BxInd, &YxInd);
        ad.f();
        std::size_t jj=0;
        for (std::size_t iii=0 ; iii<Yx.size1() ; iii++)
          for (std::size_t jjj=0 ; jjj<Yx.size2() ; jjj++)
            Jacobian(ii, jj++) = (Yx_(iii,jjj) - Yx(iii,jjj))/epsilon;
        ii++;
        Bx_(i,j) = Bx(i,j);
      if (jj != nOut)
        throw ExcKal(std::string("AutoDiff: : Incorrect size in gradNum"));
      }
    if (ii != nIn)
      throw ExcKal(std::string("AutoDiff: : Incorrect size in gradNum"));
  }

  const T a;
  const matrix<T>& Bx;
  matrix<T>& Yx;
  const matrix<std::size_t>& BxInd;
  const matrix<std::size_t>& YxInd;
};

template <class T>
class AD_adx_Bs_cx : public AutoDiff<T>
{ // y(x) = diag(a(x)) * Bsparse * c(x)
public:
  AD_adx_Bs_cx(const vector<T>* adxp, const SparseMatrix<T>* Bsp, const vector<T>* cxp, vector<T>* yxp, 
    const vector<std::size_t>* adxIndp, const vector<std::size_t>* cxIndp, const vector<std::size_t>* yxIndp) 
    : adx(*adxp), Bs(*Bsp), cx(*cxp), yx(*yxp), adxInd(*adxIndp), cxInd(*cxIndp), yxInd(*yxIndp) 
  {
    // Check dimensions
    std::size_t m = Bs.size1();
    std::size_t n = Bs.size2();
    std::string strFunc("AD_adx_Bs_cx");
    checkDim(strFunc, "adx", adx, m);
    checkDim(strFunc, "Bs", Bs, m, n);
    checkDim(strFunc, "cx", cx, n);
    checkDim(strFunc, "yx", yx, m);
    checkDim(strFunc, "adxInd", adxInd, m);
    checkDim(strFunc, "cxInd", cxInd, n);
    checkDim(strFunc, "yxInd", yxInd, m);
  }
  virtual ~AD_adx_Bs_cx() {}

  virtual void count(std::size_t& rMax, std::size_t& dMax) const 
  {
    const std::size_t m = Bs.size1();
    const std::size_t n = Bs.size2();
    const vector<T>& Bv = Bs.val();
    rMax += m;
    dMax += 2*Bv.size();
  }
  virtual void f() 
  {
    const std::size_t m = Bs.size1();

    const vector<std::size_t>& Br = Bs.row();
    const vector<std::size_t>& Bc = Bs.col();
    const vector<T>& Bv = Bs.val();

    for (std::size_t i=0 ; i<m ; i++)
    {
      yx(i) = 0.0;
      for (std::size_t k=Br(i) ; k<Br(i+1) ; k++)
        yx(i) += adx(i)*Bv(k)*cx(Bc(k));
    }
  }
  virtual void fAD(vector<T>& d, vector<std::size_t>& col, vector<std::size_t>& row, vector<std::size_t>& rInd, std::size_t& nr, std::size_t& nd) 
  {
    const std::size_t m = Bs.size1();

    const vector<std::size_t>& Br = Bs.row();
    const vector<std::size_t>& Bc = Bs.col();
    const vector<T>& Bv = Bs.val();

    for (std::size_t i=0 ; i<m ; i++)
    {
      rInd(nr) = yxInd(i);
      row(nr++) = nd;
      yx(i) = 0.0;
      for (std::size_t k=Br(i) ; k<Br(i+1) ; k++)
      {
        yx(i) += adx(i)*Bv(k)*cx(Bc(k));
        d(nd) = Bv(k)*cx(Bc(k));                // d yx(i) / d adx(i) 
        col(nd++) = adxInd(i);
        d(nd) = Bv(k)*adx(i);                   // d yx(i) / d cx(Bc(k)) 
        col(nd++) = cxInd(Bc(k));
      }
    }
  }
protected:
  virtual void countSize(std::size_t& nIn, std::size_t& nOut)
  {
    nIn = adx.size() + cx.size();
    nOut = yx.size();
  }

  virtual AutoDiff<T>* createLocalObj(std::size_t& nc)
  { // Deallocated in freeLocalObj
    vector<std::size_t>* adxInd_= new vector<std::size_t>(adx.size()); 
    vector<std::size_t>* cxInd_= new vector<std::size_t>(cx.size()); 
    vector<std::size_t>* yxInd_= new vector<std::size_t>(yx.size()); 
    nc = 0;
    for (std::size_t i=0 ; i<adx.size() ; i++)
      (*adxInd_)(i) = nc++;
    for (std::size_t i=0 ; i<cx.size() ; i++)
      (*cxInd_)(i) = nc++;
    for (std::size_t i=0 ; i<yx.size() ; i++)
      (*yxInd_)(i) = nc++;
    std::size_t nIn;
    std::size_t nOut;
    countSize(nIn, nOut);
    if (nc != nIn+nOut)
      throw ExcKal(std::string("AutoDiff: : Incorrect size in createLocalObj"));

    return new AD_adx_Bs_cx(&adx, &Bs, &cx, &yx, adxInd_, cxInd_, yxInd_);
  }

  virtual void freeLocalObj()
  { // Free objects created by createLocalObj
    delete &adxInd;
    delete &cxInd;
    delete &yxInd;
  }

  virtual void testADgrad(const vector<T>& dOdx, vector<T>& grad) // Compute gradient with respect to input arguments
  { 
    std::size_t k=0;
    for (std::size_t i=0 ; i<adx.size() ; i++)
      grad(k++) = dOdx(adxInd(i));
    for (std::size_t i=0 ; i<cx.size() ; i++)
      grad(k++) = dOdx(cxInd(i));
  }

  virtual void gradNum(const T epsilon, matrix<T>& Jacobian)
  {
    std::size_t nIn;
    std::size_t nOut;
    countSize(nIn, nOut);
    Jacobian.resize(nIn, nOut);
    f();
    vector<T> adx_(adx);
    vector<T> cx_(cx);
    vector<T> yx_(yx);

    std::size_t ii=0;
    for (std::size_t i=0 ; i<adx.size() ; i++)
    {        
      adx_(i) += epsilon;
      AD_adx_Bs_cx<T> ad(&adx_, &Bs, &cx_, &yx_, &adxInd, &cxInd, &yxInd);
      ad.f();
      std::size_t jj=0;
      for (std::size_t iii=0 ; iii<yx.size() ; iii++)
        Jacobian(ii, jj++) = (yx_(iii) - yx(iii))/epsilon;
      ii++;
      adx_(i) = adx(i);
    if (jj != nOut)
      throw ExcKal(std::string("AutoDiff: : Incorrect size in gradNum"));
    }
    for (std::size_t i=0 ; i<cx.size() ; i++)
    {        
      cx_(i) += epsilon;
      AD_adx_Bs_cx<T> ad(&adx_, &Bs, &cx_, &yx_, &adxInd, &cxInd, &yxInd);
      ad.f();
      std::size_t jj=0;
      for (std::size_t iii=0 ; iii<yx.size() ; iii++)
        Jacobian(ii, jj++) = (yx_(iii) - yx(iii))/epsilon;
      ii++;
      cx_(i) = cx(i);
    if (jj != nOut)
      throw ExcKal(std::string("AutoDiff: : Incorrect size in gradNum"));
    }
    if (ii != nIn)
      throw ExcKal(std::string("AutoDiff: : Incorrect size in gradNum"));
  }

  const vector<T>& adx;
  const SparseMatrix<T>& Bs;
  const vector<T>& cx;
  vector<T>& yx;
  const vector<std::size_t>& adxInd;
  const vector<std::size_t>& cxInd;
  const vector<std::size_t>& yxInd;
};

template <class T>
class AD_As_bdx : public AutoDiff<T>
{ // Y(x) = Asparse * diag(bd(x))
public:
  AD_As_bdx(const SparseMatrix<T>* Asp, const vector<T>* bdxp, matrix<T>* Yxp, 
    const vector<std::size_t>* bdxIndp, const matrix<std::size_t>* YxIndp) 
    : As(*Asp), bdx(*bdxp), Yx(*Yxp), bdxInd(*bdxIndp), YxInd(*YxIndp) 
  {
    // Check dimensions
    std::size_t m = As.size1();
    std::size_t n = As.size2();
    std::string strFunc("AD_As_bdx");
    checkDim(strFunc, "As", As, m, n);
    checkDim(strFunc, "bdx", bdx, n);
    checkDim(strFunc, "Yx", Yx, m, n);
    checkDim(strFunc, "bdxInd", bdxInd, n);
    checkDim(strFunc, "YxInd", YxInd, m, n);
  }
  virtual ~AD_As_bdx() {}

  virtual void count(std::size_t& rMax, std::size_t& dMax) const 
  {
    std::size_t m = As.size1();
    std::size_t n = As.size2();
    const vector<T>& Av = As.val();
    rMax += m*n; 
    dMax += n*Av.size(); 
  }
  virtual void f() 
  {
    std::size_t m = As.size1();
    std::size_t n = As.size2();

    const vector<std::size_t>& Ar = As.row();
    const vector<std::size_t>& Ac = As.col();
    const vector<T>& Av = As.val();

    for (std::size_t i=0 ; i<m ; i++)
      for (std::size_t j=0 ; j<n ; j++)
      {
        Yx(i,j) = 0.0;
        for (std::size_t k=Ar(i) ; k<Ar(i+1) ; k++)
          if (Ac(k) == j)
            Yx(i,j) += Av(k)*bdx(Ac(k));
      }
  }
  virtual void fAD(vector<T>& d, vector<std::size_t>& col, vector<std::size_t>& row, vector<std::size_t>& rInd, std::size_t& nr, std::size_t& nd) 
  {
    std::size_t m = As.size1();
    std::size_t n = As.size2();

    const vector<std::size_t>& Ar = As.row();
    const vector<std::size_t>& Ac = As.col();
    const vector<T>& Av = As.val();

    for (std::size_t i=0 ; i<m ; i++)
      for (std::size_t j=0 ; j<n ; j++)
      {
        rInd(nr) = YxInd(i,j);
        row(nr++) = nd;
        Yx(i,j) = 0.0;
        for (std::size_t k=Ar(i) ; k<Ar(i+1) ; k++)
          if (Ac(k) == j)
          {
            Yx(i,j) += Av(k)*bdx(Ac(k));
            d(nd) = Av(k);                      // d Yx(i,j) / d bdx(Ac(k)) 
            col(nd++) = bdxInd(Ac(k));
          }
      }
  }

protected:
  virtual void countSize(std::size_t& nIn, std::size_t& nOut)
  {
    nIn = bdx.size();
    nOut = Yx.size1() * Yx.size2();
  }

  virtual AutoDiff<T>* createLocalObj(std::size_t& nc)
  { // Deallocated in freeLocalObj
    vector<std::size_t>* bdxInd_= new vector<std::size_t>(bdx.size()); 
    matrix<std::size_t>* YxInd_= new matrix<std::size_t>(Yx.size1(), Yx.size2()); 
    nc = 0;
    for (std::size_t i=0 ; i<bdx.size() ; i++)
        (*bdxInd_)(i) = nc++;
    for (std::size_t i=0 ; i<Yx.size1() ; i++)
      for (std::size_t j=0 ; j<Yx.size2() ; j++)
        (*YxInd_)(i,j) = nc++;
    std::size_t nIn;
    std::size_t nOut;
    countSize(nIn, nOut);
    if (nc != nIn+nOut)
      throw ExcKal(std::string("AutoDiff: : Incorrect size in createLocalObj"));

    return new AD_As_bdx(&As, &bdx, &Yx, bdxInd_, YxInd_);
  }

  virtual void freeLocalObj()
  { // Free objects created by createLocalObj
    delete &bdxInd;
    delete &YxInd;
  }

  virtual void testADgrad(const vector<T>& dOdx, vector<T>& grad) // Compute gradient with respect to input arguments
  { 
    std::size_t k=0;
    for (std::size_t i=0 ; i<bdx.size() ; i++)
      grad(k++) = dOdx(bdxInd(i));
  }

  virtual void gradNum(const T epsilon, matrix<T>& Jacobian)
  {
    std::size_t nIn;
    std::size_t nOut;
    countSize(nIn, nOut);
    Jacobian.resize(nIn, nOut);
    f();
    vector<T> bdx_(bdx);
    matrix<T> Yx_(Yx);

    std::size_t ii=0;
    for (std::size_t i=0 ; i<bdx.size() ; i++)
    {        
      bdx_(i) += epsilon;
      AD_As_bdx<T> ad(&As, &bdx_, &Yx_, &bdxInd, &YxInd);
      ad.f();
      std::size_t jj=0;
      for (std::size_t iii=0 ; iii<Yx.size1() ; iii++)
        for (std::size_t jjj=0 ; jjj<Yx.size2() ; jjj++)
          Jacobian(ii, jj++) = (Yx_(iii,jjj) - Yx(iii,jjj))/epsilon;
      ii++;
      bdx_(i) = bdx(i);
    if (jj != nOut)
      throw ExcKal(std::string("AutoDiff: : Incorrect size in gradNum"));
    }
    if (ii != nIn)
      throw ExcKal(std::string("AutoDiff: : Incorrect size in gradNum"));
  }

  const SparseMatrix<T>& As;
  const vector<T>& bdx;
  matrix<T>& Yx;
  const vector<std::size_t>& bdxInd;
  const matrix<std::size_t>& YxInd;
};

template <class T>
class AD_As_Bx : public AutoDiff<T>
{ // Y(x) = Asparse * B(x)
public:
  AD_As_Bx(const SparseMatrix<T>* Asp, const matrix<T>* Bxp, matrix<T>* Yxp, 
    const matrix<std::size_t>* BxIndp, const matrix<std::size_t>* YxIndp) 
    : As(*Asp), Bx(*Bxp), Yx(*Yxp), BxInd(*BxIndp), YxInd(*YxIndp) 
  {
    // Check dimensions
    std::size_t m = As.size1();
    std::size_t K = As.size2();
    std::size_t n = Bx.size2();
    std::string strFunc("AD_As_Bx");
    checkDim(strFunc, "As", As, m, K);
    checkDim(strFunc, "Bx", Bx, K, n);
    checkDim(strFunc, "Yx", Yx, m, n);
    checkDim(strFunc, "BxInd", BxInd, K, n);
    checkDim(strFunc, "YxInd", YxInd, m, n);
  }
  virtual ~AD_As_Bx() {}

  virtual void count(std::size_t& rMax, std::size_t& dMax) const 
  {
    std::size_t m = As.size1();
    std::size_t n = As.size2();
    const vector<T>& Av = As.val();
    rMax += m*n; 
    dMax += n*Av.size(); 
  }
  virtual void f() 
  {
    std::size_t m = As.size1();
    std::size_t n = Bx.size2();

    const vector<std::size_t>& Ar = As.row();
    const vector<std::size_t>& Ac = As.col();
    const vector<T>& Av = As.val();

    for (std::size_t i=0 ; i<m ; i++)
      for (std::size_t j=0 ; j<n ; j++)
      {
        Yx(i,j) = 0.0;
        for (std::size_t k=Ar(i) ; k<Ar(i+1) ; k++)
          Yx(i,j) += Av(k)*Bx(Ac(k), j);
      }
  }
  virtual void fAD(vector<T>& d, vector<std::size_t>& col, vector<std::size_t>& row, vector<std::size_t>& rInd, std::size_t& nr, std::size_t& nd) 
  {
    std::size_t m = As.size1();
    std::size_t n = Bx.size2();

    const vector<std::size_t>& Ar = As.row();
    const vector<std::size_t>& Ac = As.col();
    const vector<T>& Av = As.val();

    for (std::size_t i=0 ; i<m ; i++)
      for (std::size_t j=0 ; j<n ; j++)
      {
        rInd(nr) = YxInd(i,j);
        row(nr++) = nd;
        Yx(i,j) = 0.0;
        for (std::size_t k=Ar(i) ; k<Ar(i+1) ; k++)
        {
          Yx(i,j) += Av(k)*Bx(Ac(k), j);
          d(nd) = Av(k);                      // d Yx(i,j) / d Bx(Ac(k), j) 
          col(nd++) = BxInd(Ac(k), j);
        }
      }
  }
protected:
  virtual void countSize(std::size_t& nIn, std::size_t& nOut)
  {
    nIn = Bx.size1() * Bx.size2();
    nOut = Yx.size1() * Yx.size2();
  }

  virtual AutoDiff<T>* createLocalObj(std::size_t& nc)
  { // Deallocated in freeLocalObj
    matrix<std::size_t>* BxInd_= new matrix<std::size_t>(Bx.size1(), Bx.size2()); 
    matrix<std::size_t>* YxInd_= new matrix<std::size_t>(Yx.size1(), Yx.size2()); 
    nc = 0;
    for (std::size_t i=0 ; i<Bx.size1() ; i++)
      for (std::size_t j=0 ; j<Bx.size2() ; j++)
        (*BxInd_)(i,j) = nc++;
    for (std::size_t i=0 ; i<Yx.size1() ; i++)
      for (std::size_t j=0 ; j<Yx.size2() ; j++)
        (*YxInd_)(i,j) = nc++;
    std::size_t nIn;
    std::size_t nOut;
    countSize(nIn, nOut);
    if (nc != nIn+nOut)
      throw ExcKal(std::string("AutoDiff: : Incorrect size in createLocalObj"));

    return new AD_As_Bx(&As, &Bx, &Yx, BxInd_, YxInd_);
  }

  virtual void freeLocalObj()
  { // Free objects created by createLocalObj
    delete &BxInd;
    delete &YxInd;
  }

  virtual void testADgrad(const vector<T>& dOdx, vector<T>& grad) // Compute gradient with respect to input arguments
  { 
    std::size_t k=0;
    for (std::size_t i=0 ; i<Bx.size1() ; i++)
      for (std::size_t j=0 ; j<Bx.size2() ; j++)
        grad(k++) = dOdx(BxInd(i,j));
  }

  virtual void gradNum(const T epsilon, matrix<T>& Jacobian)
  {
    std::size_t nIn;
    std::size_t nOut;
    countSize(nIn, nOut);
    Jacobian.resize(nIn, nOut);
    f();
    matrix<T> Bx_(Bx);
    matrix<T> Yx_(Yx);

    std::size_t ii=0;
    for (std::size_t i=0 ; i<Bx.size1() ; i++)
      for (std::size_t j=0 ; j<Bx.size2() ; j++)
      {        
        Bx_(i,j) += epsilon;
        AD_As_Bx<T> ad(&As, &Bx_, &Yx_, &BxInd, &YxInd);
        ad.f();
        std::size_t jj=0;
        for (std::size_t iii=0 ; iii<Yx.size1() ; iii++)
          for (std::size_t jjj=0 ; jjj<Yx.size2() ; jjj++)
            Jacobian(ii, jj++) = (Yx_(iii,jjj) - Yx(iii,jjj))/epsilon;
        ii++;
        Bx_(i,j) = Bx(i,j);
        if (jj != nOut)
          throw ExcKal(std::string("AutoDiff: : Incorrect size in gradNum"));
      }
    if (ii != nIn)
      throw ExcKal(std::string("AutoDiff: : Incorrect size in gradNum"));
  }

  const SparseMatrix<T>& As;
  const matrix<T>& Bx;
  matrix<T>& Yx;
  const matrix<std::size_t>& BxInd;
  const matrix<std::size_t>& YxInd;
};

template <class T>
class AD_Ax_BsT : public AutoDiff<T>
{ // Y(x) = A(x) * transpose(Bsparse) - note that Bs is not the transpose
public:
  AD_Ax_BsT(const matrix<T>* Axp, const SparseMatrix<T>* Bsp, matrix<T>* Yxp, 
    const matrix<std::size_t>* AxIndp, const matrix<std::size_t>* YxIndp) 
    : Ax(*Axp), Bs(*Bsp), Yx(*Yxp), AxInd(*AxIndp), YxInd(*YxIndp) 
  {
    // Check dimensions
    std::size_t m = Ax.size1();
    std::size_t K = Ax.size2();
    std::size_t n = Bs.size1();
    std::string strFunc("AD_Ax_BsT");
    checkDim(strFunc, "Ax", Ax, m, K);
    checkDim(strFunc, "Bs", Bs, n, K);
    checkDim(strFunc, "Yx", Yx, m, n);
    checkDim(strFunc, "AxInd", AxInd, m, K);
    checkDim(strFunc, "YxInd", YxInd, m, n);
  }
  virtual ~AD_Ax_BsT() {}

  virtual void count(std::size_t& rMax, std::size_t& dMax) const 
  {
    std::size_t m = Ax.size1();
    std::size_t K = Ax.size2();
    std::size_t n = Bs.size1();
    const vector<T>& Bv = Bs.val();
    rMax += m*n; 
    dMax += m*Bv.size(); 
  }
  virtual void f() 
  {
    std::size_t m = Ax.size1();
    std::size_t K = Ax.size2();
    std::size_t n = Bs.size1();

    const vector<std::size_t>& Br = Bs.row();
    const vector<std::size_t>& Bc = Bs.col();
    const vector<T>& Bv = Bs.val();

    for (std::size_t i=0 ; i<m ; i++)
      for (std::size_t j=0 ; j<n ; j++)
      {
        Yx(i,j) = 0.0;
        for (std::size_t k=Br(j) ; k<Br(j+1) ; k++)
          Yx(i,j) += Ax(i, Bc(k))*Bv(k);
      }
  }
  virtual void fAD(vector<T>& d, vector<std::size_t>& col, vector<std::size_t>& row, vector<std::size_t>& rInd, std::size_t& nr, std::size_t& nd) 
  {
    std::size_t m = Ax.size1();
    std::size_t K = Ax.size2();
    std::size_t n = Bs.size1();

    const vector<std::size_t>& Br = Bs.row();
    const vector<std::size_t>& Bc = Bs.col();
    const vector<T>& Bv = Bs.val();

    for (std::size_t i=0 ; i<m ; i++)
      for (std::size_t j=0 ; j<n ; j++)
      {
        rInd(nr) = YxInd(i,j);
        row(nr++) = nd;
        Yx(i,j) = 0.0;
        for (std::size_t k=Br(j) ; k<Br(j+1) ; k++)
        {
          Yx(i,j) += Ax(i, Bc(k))*Bv(k);
          d(nd) = Bv(k);                        // d Yx(i,j) / d Ax(i, Bc(k)) 
          col(nd++) = AxInd(i, Bc(k));
        }
      }
  }
protected:
  virtual void countSize(std::size_t& nIn, std::size_t& nOut)
  {
    nIn = Ax.size1() * Ax.size2();
    nOut = Yx.size1() * Yx.size2();
  }

  virtual AutoDiff<T>* createLocalObj(std::size_t& nc)
  { // Deallocated in freeLocalObj
    matrix<std::size_t>* AxInd_= new matrix<std::size_t>(Ax.size1(), Ax.size2()); 
    matrix<std::size_t>* YxInd_= new matrix<std::size_t>(Yx.size1(), Yx.size2()); 
    nc = 0;
    for (std::size_t i=0 ; i<Ax.size1() ; i++)
      for (std::size_t j=0 ; j<Ax.size2() ; j++)
        (*AxInd_)(i,j) = nc++;
    for (std::size_t i=0 ; i<Yx.size1() ; i++)
      for (std::size_t j=0 ; j<Yx.size2() ; j++)
        (*YxInd_)(i,j) = nc++;
    std::size_t nIn;
    std::size_t nOut;
    countSize(nIn, nOut);
    if (nc != nIn+nOut)
      throw ExcKal(std::string("AutoDiff: : Incorrect size in createLocalObj"));

    return new AD_Ax_BsT(&Ax, &Bs, &Yx, AxInd_, YxInd_);
  }

  virtual void freeLocalObj()
  { // Free objects created by createLocalObj
    delete &AxInd;
    delete &YxInd;
  }

  virtual void testADgrad(const vector<T>& dOdx, vector<T>& grad) // Compute gradient with respect to input arguments
  { 
    std::size_t k=0;
    for (std::size_t i=0 ; i<Ax.size1() ; i++)
      for (std::size_t j=0 ; j<Ax.size2() ; j++)
        grad(k++) = dOdx(AxInd(i,j));
  }

  virtual void gradNum(const T epsilon, matrix<T>& Jacobian)
  {
    std::size_t nIn;
    std::size_t nOut;
    countSize(nIn, nOut);
    Jacobian.resize(nIn, nOut);
    f();
    matrix<T> Ax_(Ax);
    matrix<T> Yx_(Yx);

    std::size_t ii=0;
    for (std::size_t i=0 ; i<Ax.size1() ; i++)
      for (std::size_t j=0 ; j<Ax.size2() ; j++)
      {        
        Ax_(i,j) += epsilon;
        AD_Ax_BsT<T> ad(&Ax_, &Bs, &Yx_, &AxInd, &YxInd);
        ad.f();
        std::size_t jj=0;
        for (std::size_t iii=0 ; iii<Yx.size1() ; iii++)
          for (std::size_t jjj=0 ; jjj<Yx.size2() ; jjj++)
            Jacobian(ii, jj++) = (Yx_(iii,jjj) - Yx(iii,jjj))/epsilon;
        ii++;
        Ax_(i,j) = Ax(i,j);
        if (jj != nOut)
          throw ExcKal(std::string("AutoDiff: : Incorrect size in gradNum"));
      }
    if (ii != nIn)
      throw ExcKal(std::string("AutoDiff: : Incorrect size in gradNum"));
  }

  const matrix<T>& Ax;
  const SparseMatrix<T>& Bs;
  matrix<T>& Yx;
  const matrix<std::size_t>& AxInd;
  const matrix<std::size_t>& YxInd;
};

template <class T>
class AD_Ax_AxT : public AutoDiff<T>
{ // Y(x) = (A(x) + A(x)^T)/2 
public:
  AD_Ax_AxT(const matrix<T>* Axp, matrix<T>* Yxp, 
    const matrix<std::size_t>* AxIndp, const matrix<std::size_t>* YxIndp) 
    : Ax(*Axp), Yx(*Yxp), AxInd(*AxIndp), YxInd(*YxIndp) 
  {
    // Check dimensions
    std::size_t m = Ax.size1();
    std::size_t n = Ax.size2();
    std::string strFunc("AD_Ax_AxT");
    checkDim(strFunc, "Ax", Ax, m, n);
    checkDim(strFunc, "Ax", Ax, n, m);
    checkDim(strFunc, "Yx", Yx, m, n);
    checkDim(strFunc, "AxInd", AxInd, m, n);
    checkDim(strFunc, "YxInd", YxInd, m, n);
  }
  virtual ~AD_Ax_AxT() {}

  virtual void count(std::size_t& rMax, std::size_t& dMax) const 
  {
    std::size_t m = Ax.size1();
    std::size_t n = Ax.size2();
    rMax += m*n; 
    dMax += 2*m*n; 
  }
  virtual void f() 
  {
    std::size_t m = Ax.size1();
    std::size_t n = Ax.size2();

    for (std::size_t i=0 ; i<m ; i++) 
      for (std::size_t j=0 ; j<n ; j++)
      {
        T avg = (Ax(i,j)+Ax(j,i))/2;
        Yx(i,j) = avg;
      }
  }
  virtual void fAD(vector<T>& d, vector<std::size_t>& col, vector<std::size_t>& row, vector<std::size_t>& rInd, std::size_t& nr, std::size_t& nd) 
  {
    std::size_t m = Ax.size1();
    std::size_t n = Ax.size2();

    for (std::size_t i=0 ; i<m ; i++) 
      for (std::size_t j=0 ; j<n ; j++)
      {
        rInd(nr) = YxInd(i,j);
        row(nr++) = nd;
        T avg = (Ax(i,j)+Ax(j,i))/2;
        Yx(i,j) = avg;
        d(nd) = 0.5;                              // d Yx(i,j) / d Ax(i,j) 
        col(nd++) = AxInd(i,j);
        d(nd) = 0.5;                              // d Yx(i,j) / d Ax(j,i) 
        col(nd++) = AxInd(j,i);
      }
  }
protected:
  virtual void countSize(std::size_t& nIn, std::size_t& nOut)
  {
    nIn = Ax.size1() * Ax.size2();
    nOut = Yx.size1() * Yx.size2();
  }

  virtual AutoDiff<T>* createLocalObj(std::size_t& nc)
  { // Deallocated in freeLocalObj
    matrix<std::size_t>* AxInd_= new matrix<std::size_t>(Ax.size1(), Ax.size2()); 
    matrix<std::size_t>* YxInd_= new matrix<std::size_t>(Yx.size1(), Yx.size2()); 
    nc = 0;
    for (std::size_t i=0 ; i<Ax.size1() ; i++)
      for (std::size_t j=0 ; j<Ax.size2() ; j++)
        (*AxInd_)(i,j) = nc++;
    for (std::size_t i=0 ; i<Yx.size1() ; i++)
      for (std::size_t j=0 ; j<Yx.size2() ; j++)
        (*YxInd_)(i,j) = nc++;
    std::size_t nIn;
    std::size_t nOut;
    countSize(nIn, nOut);
    if (nc != nIn+nOut)
      throw ExcKal(std::string("AutoDiff: : Incorrect size in createLocalObj"));

    return new AD_Ax_AxT(&Ax, &Yx, AxInd_, YxInd_);
  }

  virtual void freeLocalObj()
  { // Free objects created by createLocalObj
    delete &AxInd;
    delete &YxInd;
  }

  virtual void testADgrad(const vector<T>& dOdx, vector<T>& grad) // Compute gradient with respect to input arguments
  { 
    std::size_t k=0;
    for (std::size_t i=0 ; i<Ax.size1() ; i++)
      for (std::size_t j=0 ; j<Ax.size2() ; j++)
        grad(k++) = dOdx(AxInd(i,j));
    if (k != grad.size())
      throw ExcKal(std::string("AutoDiff: : Incorrect size in grad"));
  }

  virtual void gradNum(const T epsilon, matrix<T>& Jacobian)
  {
    std::size_t nIn;
    std::size_t nOut;
    countSize(nIn, nOut);
    Jacobian.resize(nIn, nOut);
    f();
    matrix<T> Ax_(Ax);
    matrix<T> Yx_(Yx);

    std::size_t ii=0;
    for (std::size_t i=0 ; i<Ax.size1() ; i++)
      for (std::size_t j=0 ; j<Ax.size2() ; j++)
      {        
        Ax_(i,j) += epsilon;
        AD_Ax_AxT<T> ad(&Ax_, &Yx_, &AxInd, &YxInd);
        ad.f();
        std::size_t jj=0;
        for (std::size_t iii=0 ; iii<Yx.size1() ; iii++)
          for (std::size_t jjj=0 ; jjj<Yx.size2() ; jjj++)
            Jacobian(ii, jj++) = (Yx_(iii,jjj) - Yx(iii,jjj))/epsilon;
        ii++;
        Ax_(i,j) = Ax(i,j);
        if (jj != nOut)
          throw ExcKal(std::string("AutoDiff: : Incorrect size in gradNum"));
      }
    if (ii != nIn)
      throw ExcKal(std::string("AutoDiff: : Incorrect size in gradNum"));
  }

  const matrix<T>& Ax;
  matrix<T>& Yx;
  const matrix<std::size_t>& AxInd;
  const matrix<std::size_t>& YxInd;
};

template <class T>
class AD_Ax_BsT_Cx : public AutoDiff<T>
{ // Y(x) = A(x) * transpose(Bsparse) + C(x) - note that Bs is not the transpose
public:
  AD_Ax_BsT_Cx(const matrix<T>* Axp, const SparseMatrix<T>* Bsp, const matrix<T>* Cxp, matrix<T>* Yxp, 
    const matrix<std::size_t>* AxIndp, const matrix<std::size_t>* CxIndp, const matrix<std::size_t>* YxIndp) 
    : Ax(*Axp), Bs(*Bsp), Cx(*Cxp), Yx(*Yxp), AxInd(*AxIndp), CxInd(*CxIndp), YxInd(*YxIndp) 
  {
    // Check dimensions
    std::size_t m = Ax.size1();
    std::size_t K = Ax.size2();
    std::size_t n = Bs.size1();
    std::string strFunc("AD_Ax_BsT");
    checkDim(strFunc, "Ax", Ax, m, K);
    checkDim(strFunc, "Bs", Bs, n, K);
    checkDim(strFunc, "Cx", Cx, m, n);
    checkDim(strFunc, "Yx", Yx, m, n);
    checkDim(strFunc, "AxInd", AxInd, m, K);
    checkDim(strFunc, "CxInd", CxInd, m, n);
    checkDim(strFunc, "YxInd", YxInd, m, n);
  }
  virtual ~AD_Ax_BsT_Cx() {}

  virtual void count(std::size_t& rMax, std::size_t& dMax) const 
  {
    std::size_t m = Ax.size1();
    std::size_t K = Ax.size2();
    std::size_t n = Bs.size1();
    const vector<T>& Bv = Bs.val();
    rMax += m*n; 
    dMax += m*n + m*Bv.size(); 
  }
  virtual void f() 
  {
    std::size_t m = Ax.size1();
    std::size_t K = Ax.size2();
    std::size_t n = Bs.size1();

    const vector<std::size_t>& Br = Bs.row();
    const vector<std::size_t>& Bc = Bs.col();
    const vector<T>& Bv = Bs.val();

    for (std::size_t i=0 ; i<m ; i++)
      for (std::size_t j=0 ; j<n ; j++)
      {
        Yx(i,j) = Cx(i,j);
        T dot = 0;
        for (std::size_t k=Br(j) ; k<Br(j+1) ; k++)
          dot += Ax(i, Bc(k))*Bv(k);
        Yx(i,j) += dot;
      }

    for (std::size_t i=0 ; i<m ; i++) // Utilize that matrix (in application) is symmetric - improves numerical stability even though mathematically equivalent
      for (std::size_t j=i+1 ; j<n ; j++)
      {
        T avg = (Yx(i,j)+Yx(j,i))/2;
        Yx(i,j) = avg;
        Yx(j,i) = avg;
      }
  }
  virtual void fAD(vector<T>& d, vector<std::size_t>& col, vector<std::size_t>& row, vector<std::size_t>& rInd, std::size_t& nr, std::size_t& nd) 
  {
    std::size_t m = Ax.size1();
    std::size_t K = Ax.size2();
    std::size_t n = Bs.size1();

    const vector<std::size_t>& Br = Bs.row();
    const vector<std::size_t>& Bc = Bs.col();
    const vector<T>& Bv = Bs.val();

    for (std::size_t i=0 ; i<m ; i++)
      for (std::size_t j=0 ; j<n ; j++)
      {
        rInd(nr) = YxInd(i,j);
        row(nr++) = nd;
        Yx(i,j) = Cx(i,j);
        d(nd) = 1;                              // d Yx(i,j) / d Cx(i,j) 
        col(nd++) = CxInd(i,j);
        T dot = 0;
        for (std::size_t k=Br(j) ; k<Br(j+1) ; k++)
        {
          dot += Ax(i, Bc(k))*Bv(k);
          d(nd) = Bv(k);                        // d Yx(i,j) / d Ax(i, Bc(k)) 
          col(nd++) = AxInd(i, Bc(k));
        }
        Yx(i,j) += dot;
      }

    for (std::size_t i=0 ; i<m ; i++) // Utilize that matrix (in application) is symmetric - improves numerical stability even though mathematically equivalent
      for (std::size_t j=i+1 ; j<n ; j++)
      {
        T avg = (Yx(i,j)+Yx(j,i))/2;
        Yx(i,j) = avg;
        Yx(j,i) = avg;
      }
  }
protected:
  virtual void countSize(std::size_t& nIn, std::size_t& nOut)
  {
    nIn = Ax.size1() * Ax.size2() + Cx.size1() * Cx.size2();
    nOut = Yx.size1() * Yx.size2();
  }

  virtual AutoDiff<T>* createLocalObj(std::size_t& nc)
  { // Deallocated in freeLocalObj
    matrix<std::size_t>* AxInd_= new matrix<std::size_t>(Ax.size1(), Ax.size2()); 
    matrix<std::size_t>* CxInd_= new matrix<std::size_t>(Cx.size1(), Cx.size2()); 
    matrix<std::size_t>* YxInd_= new matrix<std::size_t>(Yx.size1(), Yx.size2()); 
    nc = 0;
    for (std::size_t i=0 ; i<Ax.size1() ; i++)
      for (std::size_t j=0 ; j<Ax.size2() ; j++)
        (*AxInd_)(i,j) = nc++;
    for (std::size_t i=0 ; i<Cx.size1() ; i++)
      for (std::size_t j=0 ; j<Cx.size2() ; j++)
        (*CxInd_)(i,j) = nc++;
    for (std::size_t i=0 ; i<Yx.size1() ; i++)
      for (std::size_t j=0 ; j<Yx.size2() ; j++)
        (*YxInd_)(i,j) = nc++;
    std::size_t nIn;
    std::size_t nOut;
    countSize(nIn, nOut);
    if (nc != nIn+nOut)
      throw ExcKal(std::string("AutoDiff: : Incorrect size in createLocalObj"));

    return new AD_Ax_BsT_Cx(&Ax, &Bs, &Cx, &Yx, AxInd_, CxInd_, YxInd_);
  }

  virtual void freeLocalObj()
  { // Free objects created by createLocalObj
    delete &AxInd;
    delete &CxInd;
    delete &YxInd;
  }

  virtual void testADgrad(const vector<T>& dOdx, vector<T>& grad) // Compute gradient with respect to input arguments
  { 
    std::size_t k=0;
    for (std::size_t i=0 ; i<Ax.size1() ; i++)
      for (std::size_t j=0 ; j<Ax.size2() ; j++)
        grad(k++) = dOdx(AxInd(i,j));
    for (std::size_t i=0 ; i<Cx.size1() ; i++)
      for (std::size_t j=0 ; j<Cx.size2() ; j++)
        grad(k++) = dOdx(CxInd(i,j));
    if (k != grad.size())
      throw ExcKal(std::string("AutoDiff: : Incorrect size in grad"));
  }

  virtual void gradNum(const T epsilon, matrix<T>& Jacobian)
  {
    std::size_t nIn;
    std::size_t nOut;
    countSize(nIn, nOut);
    Jacobian.resize(nIn, nOut);
    f();
    matrix<T> Ax_(Ax);
    matrix<T> Cx_(Cx);
    matrix<T> Yx_(Yx);

    std::size_t ii=0;
    for (std::size_t i=0 ; i<Ax.size1() ; i++)
      for (std::size_t j=0 ; j<Ax.size2() ; j++)
      {        
        Ax_(i,j) += epsilon;
        AD_Ax_BsT_Cx<T> ad(&Ax_, &Bs, &Cx_, &Yx_, &AxInd, &CxInd, &YxInd);
        ad.f();
        std::size_t jj=0;
        for (std::size_t iii=0 ; iii<Yx.size1() ; iii++)
          for (std::size_t jjj=0 ; jjj<Yx.size2() ; jjj++)
            Jacobian(ii, jj++) = (Yx_(iii,jjj) - Yx(iii,jjj))/epsilon;
        ii++;
        Ax_(i,j) = Ax(i,j);
        if (jj != nOut)
          throw ExcKal(std::string("AutoDiff: : Incorrect size in gradNum"));
      }
    for (std::size_t i=0 ; i<Cx.size1() ; i++)
      for (std::size_t j=0 ; j<Cx.size2() ; j++)
      {        
        Cx_(i,j) += epsilon;
        AD_Ax_BsT_Cx<T> ad(&Ax_, &Bs, &Cx_, &Yx_, &AxInd, &CxInd, &YxInd);
        ad.f();
        std::size_t jj=0;
        for (std::size_t iii=0 ; iii<Yx.size1() ; iii++)
          for (std::size_t jjj=0 ; jjj<Yx.size2() ; jjj++)
            Jacobian(ii, jj++) = (Yx_(iii,jjj) - Yx(iii,jjj))/epsilon;
        ii++;
        Cx_(i,j) = Cx(i,j);
        if (jj != nOut)
          throw ExcKal(std::string("AutoDiff: : Incorrect size in gradNum"));
      }
    if (ii != nIn)
      throw ExcKal(std::string("AutoDiff: : Incorrect size in gradNum"));
  }

  const matrix<T>& Ax;
  const SparseMatrix<T>& Bs;
  const matrix<T>& Cx;
  matrix<T>& Yx;
  const matrix<std::size_t>& AxInd;
  const matrix<std::size_t>& CxInd;
  const matrix<std::size_t>& YxInd;
};

template <class T>
class AD_adx_Bx_cdx : public AutoDiff<T>
{ // y(x) = diag(a(x)) * B(x) * diag(c(x))
public:
  AD_adx_Bx_cdx(const vector<T>* adxp, const matrix<T>* Bxp, const vector<T>* cdxp, matrix<T>* Yxp, 
    const vector<std::size_t>* adxIndp, const matrix<std::size_t>* BxIndp, const vector<std::size_t>* cdxIndp, const matrix<std::size_t>* YxIndp) 
    : adx(*adxp), Bx(*Bxp), cdx(*cdxp), Yx(*Yxp), adxInd(*adxIndp), BxInd(*BxIndp), cdxInd(*cdxIndp), YxInd(*YxIndp) 
  {
    // Check dimensions
    std::size_t m = Bx.size1();
    std::size_t n = Bx.size2();
    std::string strFunc("AD_adx_Bx_cdx");
    checkDim(strFunc, "adx", adx, m);
    checkDim(strFunc, "Bx", Bx, m, n);
    checkDim(strFunc, "cdx", cdx, n);
    checkDim(strFunc, "Yx", Yx, m, n);
    checkDim(strFunc, "adxInd", adxInd, m);
    checkDim(strFunc, "BxInd", BxInd, m, n);
    checkDim(strFunc, "cdxInd", cdxInd, n);
    checkDim(strFunc, "YxInd", YxInd, m, n);
  }
  virtual ~AD_adx_Bx_cdx() {}

  virtual void count(std::size_t& rMax, std::size_t& dMax) const 
  {
    std::size_t m = Bx.size1();
    std::size_t n = Bx.size2();
    rMax += m*n;
    dMax += 3*m*n;
  }
  virtual void f() 
  {
    std::size_t m = Bx.size1();
    std::size_t n = Bx.size2();

    for (std::size_t i=0 ; i<m ; i++)
      for (std::size_t j=0 ; j<n ; j++)
        Yx(i,j) = adx(i)*Bx(i,j)*cdx(j);
  }
  virtual void fAD(vector<T>& d, vector<std::size_t>& col, vector<std::size_t>& row, vector<std::size_t>& rInd, std::size_t& nr, std::size_t& nd) 
  {
    std::size_t m = Bx.size1();
    std::size_t n = Bx.size2();

    for (std::size_t i=0 ; i<m ; i++)
      for (std::size_t j=0 ; j<n ; j++)
      {
        rInd(nr) = YxInd(i,j);
        row(nr++) = nd;
        Yx(i,j) = adx(i)*Bx(i,j)*cdx(j);
        d(nd) = Bx(i,j)*cdx(j);           // d Yx(i,j) / d adx(i) 
        col(nd++) = adxInd(i);
        d(nd) = adx(i)*cdx(j);            // d Yx(i,j) / d Bx(i,j) 
        col(nd++) = BxInd(i,j);
        d(nd) = Bx(i,j)*adx(i);           // d Yx(i,j) / d cdx(j) 
        col(nd++) = cdxInd(j);
      }
  }
protected:
  virtual void countSize(std::size_t& nIn, std::size_t& nOut)
  {
    nIn = adx.size() + Bx.size1() * Bx.size2() + cdx.size();
    nOut = Yx.size1() * Yx.size2();
  }

  virtual AutoDiff<T>* createLocalObj(std::size_t& nc)
  { // Deallocated in freeLocalObj
    vector<std::size_t>* adxInd_= new vector<std::size_t>(adx.size()); 
    matrix<std::size_t>* BxInd_= new matrix<std::size_t>(Bx.size1(), Bx.size2()); 
    vector<std::size_t>* cdxInd_= new vector<std::size_t>(cdx.size()); 
    matrix<std::size_t>* YxInd_= new matrix<std::size_t>(Yx.size1(), Yx.size2()); 
    nc = 0;
    for (std::size_t i=0 ; i<adx.size() ; i++)
      (*adxInd_)(i) = nc++;
    for (std::size_t i=0 ; i<Bx.size1() ; i++)
      for (std::size_t j=0 ; j<Bx.size2() ; j++)
        (*BxInd_)(i,j) = nc++;
    for (std::size_t i=0 ; i<cdx.size() ; i++)
      (*cdxInd_)(i) = nc++;
    for (std::size_t i=0 ; i<Yx.size1() ; i++)
      for (std::size_t j=0 ; j<Yx.size2() ; j++)
        (*YxInd_)(i,j) = nc++;
    std::size_t nIn;
    std::size_t nOut;
    countSize(nIn, nOut);
    if (nc != nIn+nOut)
      throw ExcKal(std::string("AutoDiff: : Incorrect size in createLocalObj"));


    return new AD_adx_Bx_cdx(&adx, &Bx, &cdx, &Yx, adxInd_, BxInd_, cdxInd_, YxInd_);
  }

  virtual void freeLocalObj()
  { // Free objects created by createLocalObj
    delete &adxInd;
    delete &BxInd;
    delete &cdxInd;
    delete &YxInd;
  }

  virtual void testADgrad(const vector<T>& dOdx, vector<T>& grad) // Compute gradient with respect to input arguments
  { 
    std::size_t k=0;
    for (std::size_t i=0 ; i<adx.size() ; i++)
      grad(k++) = dOdx(adxInd(i));
    for (std::size_t i=0 ; i<Bx.size1() ; i++)
      for (std::size_t j=0 ; j<Bx.size2() ; j++)
        grad(k++) = dOdx(BxInd(i,j));
    for (std::size_t i=0 ; i<cdx.size() ; i++)
      grad(k++) = dOdx(cdxInd(i));
  }

  virtual void gradNum(const T epsilon, matrix<T>& Jacobian)
  {
    std::size_t nIn;
    std::size_t nOut;
    countSize(nIn, nOut);
    Jacobian.resize(nIn, nOut);
    f();
    vector<T> adx_(adx);
    matrix<T> Bx_(Bx);
    vector<T> cdx_(cdx);
    matrix<T> Yx_(Yx);

    std::size_t ii=0;
    for (std::size_t i=0 ; i<adx.size() ; i++)
    {        
      adx_(i) += epsilon;
      AD_adx_Bx_cdx<T> ad(&adx_, &Bx_, &cdx_, &Yx_, &adxInd, &BxInd, &cdxInd, &YxInd);
      ad.f();
      std::size_t jj=0;
      for (std::size_t iii=0 ; iii<Yx.size1() ; iii++)
        for (std::size_t jjj=0 ; jjj<Yx.size2() ; jjj++)
          Jacobian(ii, jj++) = (Yx_(iii,jjj) - Yx(iii,jjj))/epsilon;
      ii++;
      adx_(i) = adx(i);
      if (jj != nOut)
        throw ExcKal(std::string("AutoDiff: : Incorrect size in gradNum"));
    }
    for (std::size_t i=0 ; i<Bx.size1() ; i++)
      for (std::size_t j=0 ; j<Bx.size2() ; j++)
      {        
        Bx_(i,j) += epsilon;
        AD_adx_Bx_cdx<T> ad(&adx_, &Bx_, &cdx_, &Yx_, &adxInd, &BxInd, &cdxInd, &YxInd);
        ad.f();
        std::size_t jj=0;
        for (std::size_t iii=0 ; iii<Yx.size1() ; iii++)
          for (std::size_t jjj=0 ; jjj<Yx.size2() ; jjj++)
            Jacobian(ii, jj++) = (Yx_(iii,jjj) - Yx(iii,jjj))/epsilon;
        ii++;
        Bx_(i,j) = Bx(i,j);
        if (jj != nOut)
          throw ExcKal(std::string("AutoDiff: : Incorrect size in gradNum"));
      }
    for (std::size_t i=0 ; i<cdx.size() ; i++)
    {        
      cdx_(i) += epsilon;
      AD_adx_Bx_cdx<T> ad(&adx_, &Bx_, &cdx_, &Yx_, &adxInd, &BxInd, &cdxInd, &YxInd);
      ad.f();
      std::size_t jj=0;
      for (std::size_t iii=0 ; iii<Yx.size1() ; iii++)
        for (std::size_t jjj=0 ; jjj<Yx.size2() ; jjj++)
          Jacobian(ii, jj++) = (Yx_(iii,jjj) - Yx(iii,jjj))/epsilon;
      ii++;
      cdx_(i) = cdx(i);
    if (jj != nOut)
      throw ExcKal(std::string("AutoDiff: : Incorrect size in gradNum"));
    }
    if (ii != nIn)
      throw ExcKal(std::string("AutoDiff: : Incorrect size in gradNum"));
  }

  const vector<T>& adx;
  const matrix<T>& Bx;
  const vector<T>& cdx;
  matrix<T>& Yx;
  const vector<std::size_t>& adxInd;
  const matrix<std::size_t>& BxInd;
  const vector<std::size_t>& cdxInd;
  const matrix<std::size_t>& YxInd;
};

template <class T>
class AD_abx_Cs_cx : public AutoDiff<T>
{ // y(x) = a - b(x) - Cs * d(x) - Special implementation for Kalman filter
public:
  AD_abx_Cs_cx(const vector<T>* ap, const vector<T>* bxp, const SparseMatrix<T>* Csp, const std::size_t _od, const vector<T>* dxp, vector<T>* yxp,
               const vector<std::size_t>* bxIndp, const vector<std::size_t>* dxIndp, const vector<std::size_t>* yxIndp) 
    : a(*ap), bx(*bxp), Cs(*Csp), od(_od), dx(*dxp), yx(*yxp), bxInd(*bxIndp), dxInd(*dxIndp), yxInd(*yxIndp)   
  {
    // Check dimensions
    const std::size_t m = a.size();
    const std::size_t n = Cs.size2();
    const std::string strFunc("AD_abx_Cs_cx");
    checkDim(strFunc, "a", a, m);
    checkDim(strFunc, "bx", bx, m);
    checkDim(strFunc, "Cs", Cs, m, n);
    checkDim(strFunc, "dx", dx, od+n);
    checkDim(strFunc, "yx", yx, m);
    checkDim(strFunc, "bxInd", bxInd, m);
    checkDim(strFunc, "dxInd", dxInd, od+n);
    checkDim(strFunc, "yxInd", yxInd, m);
  }
  virtual ~AD_abx_Cs_cx() {}

  virtual void count(std::size_t& rMax, std::size_t& dMax) const 
  {
    const std::size_t m = a.size();
    const std::size_t n = Cs.size2();

    const vector<T>& Cv = Cs.val();

    rMax += m; 
    dMax += m + Cv.size();
  }
  virtual void f() 
  {
    const std::size_t m = a.size();
    const std::size_t n = Cs.size2();

    const vector<std::size_t>& Cr = Cs.row();
    const vector<std::size_t>& Cc = Cs.col();
    const vector<T>& Cv = Cs.val();

    for (std::size_t i=0 ; i<m ; i++)
    {
      yx(i) = a(i) - bx(i);
      T dot = 0;
      for (std::size_t k=Cr(i) ; k<Cr(i+1) ; k++)
        dot += Cv(k)*dx(od+Cc(k));
      yx(i) -= dot;
    }
  }
  virtual void fAD(vector<T>& d, vector<std::size_t>& col, vector<std::size_t>& row, vector<std::size_t>& rInd, std::size_t& nr, std::size_t& nd) 
  {
    const std::size_t m = a.size();
    const std::size_t n = Cs.size2();

    const vector<std::size_t>& Cr = Cs.row();
    const vector<std::size_t>& Cc = Cs.col();
    const vector<T>& Cv = Cs.val();

    for (std::size_t i=0 ; i<m ; i++)
    {
      rInd(nr) = yxInd(i);
      row(nr++) = nd;
      yx(i) = a(i) - bx(i);
      d(nd) = -1;                               // d yx(i) / d bx(i) 
      col(nd++) = bxInd(i);
      T dot = 0;
      for (std::size_t k=Cr(i) ; k<Cr(i+1) ; k++)
      {
        dot += Cv(k)*dx(od+Cc(k));
        d(nd) = -Cv(k);                        // d yx(i) / d dx(od+Cc(k)) 
        col(nd++) = dxInd(od+Cc(k));
      }
      yx(i) -= dot;
    }
  }
protected:
  virtual void countSize(std::size_t& nIn, std::size_t& nOut)
  {
    nIn = bx.size() + dx.size();
    nOut = yx.size();
  }

  virtual AutoDiff<T>* createLocalObj(std::size_t& nc)
  { // Deallocated in freeLocalObj
    vector<std::size_t>* bxInd_= new vector<std::size_t>(bx.size()); 
    vector<std::size_t>* dxInd_= new vector<std::size_t>(dx.size()); 
    vector<std::size_t>* yxInd_= new vector<std::size_t>(yx.size()); 
    nc = 0;
    for (std::size_t i=0 ; i<bx.size() ; i++)
      (*bxInd_)(i) = nc++;
    for (std::size_t i=0 ; i<dx.size() ; i++)
      (*dxInd_)(i) = nc++;
    for (std::size_t i=0 ; i<yx.size() ; i++)
      (*yxInd_)(i) = nc++;
    std::size_t nIn;
    std::size_t nOut;
    countSize(nIn, nOut);
    if (nc != nIn+nOut)
      throw ExcKal(std::string("AutoDiff: : Incorrect size in createLocalObj"));

    return new AD_abx_Cs_cx(&a, &bx, &Cs, od, &dx, &yx, bxInd_, dxInd_, yxInd_);
  }

  virtual void freeLocalObj()
  { // Free objects created by createLocalObj
    delete &bxInd;
    delete &dxInd;
    delete &yxInd;
  }

  virtual void testADgrad(const vector<T>& dOdx, vector<T>& grad) // Compute gradient with respect to input arguments
  { 
    std::size_t k=0;
    for (std::size_t i=0 ; i<bx.size() ; i++)
      grad(k++) = dOdx(bxInd(i));
    for (std::size_t i=0 ; i<dx.size() ; i++)
      grad(k++) = dOdx(dxInd(i));
  }

  virtual void gradNum(const T epsilon, matrix<T>& Jacobian)
  {
    std::size_t nIn;
    std::size_t nOut;
    countSize(nIn, nOut);
    Jacobian.resize(nIn, nOut);
    vector<T> bx_(bx);
    vector<T> dx_(dx);
    vector<T> yx_(yx);

    std::size_t ii=0;
    for (std::size_t i=0 ; i<bx.size() ; i++)
    {        
      bx_(i) += epsilon;
      AD_abx_Cs_cx<T> ad(&a, &bx_, &Cs, od, &dx_, &yx_, &bxInd, &dxInd, &yxInd);
      ad.f();
      std::size_t jj=0;
      for (std::size_t iii=0 ; iii<yx.size() ; iii++)
        Jacobian(ii, jj++) = (yx_(iii) - yx(iii))/epsilon;
      ii++;
      bx_(i) = bx(i);
    if (jj != nOut)
      throw ExcKal(std::string("AutoDiff: : Incorrect size in gradNum"));
    }
    for (std::size_t i=0 ; i<dx.size() ; i++)
    {        
      dx_(i) += epsilon;
      AD_abx_Cs_cx<T> ad(&a, &bx_, &Cs, od, &dx_, &yx_, &bxInd, &dxInd, &yxInd);
      ad.f();
      std::size_t jj=0;
      for (std::size_t iii=0 ; iii<yx.size() ; iii++)
        Jacobian(ii, jj++) = (yx_(iii) - yx(iii))/epsilon;
      ii++;
      dx_(i) = dx(i);
    if (jj != nOut)
      throw ExcKal(std::string("AutoDiff: : Incorrect size in gradNum"));
    }
    if (ii != nIn)
      throw ExcKal(std::string("AutoDiff: : Incorrect size in gradNum"));
  }


  const vector<T>& a;
  const vector<T>& bx;
  const SparseMatrix<T>& Cs;
  const std::size_t od; // offset in d-vector
  const vector<T>& dx;
  vector<T>& yx;
  const vector<std::size_t>& bxInd;
  const vector<std::size_t>& dxInd;
  const vector<std::size_t>& yxInd;
};

template <class T>
class AD_AxBs_Cx : public AutoDiff<T>
{ // Y(x) = [A(x) Bs] * C(x)
public:
  AD_AxBs_Cx(const matrix<T>* Axp, const SparseMatrix<T>* Bsp, const matrix<T>* Cxp, matrix<T>* Yxp, 
    const matrix<std::size_t>* AxIndp, const matrix<std::size_t>* CxIndp, const matrix<std::size_t>* YxIndp) 
    : Ax(*Axp), Bs(*Bsp), Cx(*Cxp), Yx(*Yxp), AxInd(*AxIndp), CxInd(*CxIndp), YxInd(*YxIndp) 
  {
    // Check dimensions
    const std::size_t m = Ax.size1();
    const std::size_t K1 = Ax.size2();
    const std::size_t K2 = Bs.size2();
    const std::size_t n = Cx.size2();
    const std::string strFunc("AD_AxBs_Cx");
    checkDim(strFunc, "Ax", Ax, m, K1);
    checkDim(strFunc, "Bs", Bs, m, K2);
    checkDim(strFunc, "Cx", Cx, K1+K2, n);
    checkDim(strFunc, "Yx", Yx, m, n);
    checkDim(strFunc, "AxInd", AxInd, m, K1);
    checkDim(strFunc, "CxInd", CxInd, K1+K2, n);
    checkDim(strFunc, "YxInd", YxInd, m, n);
  }
  virtual ~AD_AxBs_Cx() {}

  virtual void count(std::size_t& rMax, std::size_t& dMax) const 
  {
    const std::size_t m = Ax.size1();
    const std::size_t K1 = Ax.size2();
    const std::size_t n = Cx.size2();

    const vector<T>& Bv = Bs.val();

    rMax += m*n; 
    dMax += 2*m*n*K1 + n*Bv.size();
  }
  virtual void f() 
  {
    const std::size_t m = Ax.size1();
    const std::size_t K1 = Ax.size2();
    const std::size_t K2 = Bs.size2();
    const std::size_t n = Cx.size2();

    const vector<std::size_t>& Br = Bs.row();
    const vector<std::size_t>& Bc = Bs.col();
    const vector<T>& Bv = Bs.val();

    for (std::size_t i=0 ; i<m ; i++)
      for (std::size_t j=0 ; j<n ; j++)
      {
        T dot = 0;
        for (std::size_t k=0 ; k<K1 ; k++)
          dot += Ax(i,k)*Cx(k,j);
        Yx(i,j) = dot;
        dot = 0;
        for (std::size_t k=Br(i) ; k<Br(i+1) ; k++)
          dot += Bv(k)*Cx(K1+Bc(k), j);
        Yx(i,j) += dot;
      }
  }
  virtual void fAD(vector<T>& d, vector<std::size_t>& col, vector<std::size_t>& row, vector<std::size_t>& rInd, std::size_t& nr, std::size_t& nd) 
  {
    const std::size_t m = Ax.size1();
    const std::size_t K1 = Ax.size2();
    const std::size_t K2 = Bs.size2();
    const std::size_t n = Cx.size2();

    const vector<std::size_t>& Br = Bs.row();
    const vector<std::size_t>& Bc = Bs.col();
    const vector<T>& Bv = Bs.val();

    for (std::size_t i=0 ; i<m ; i++)
      for (std::size_t j=0 ; j<n ; j++)
      {
        rInd(nr) = YxInd(i,j);
        row(nr++) = nd;
        T dot = 0;
        for (std::size_t k=0 ; k<K1 ; k++)
        {
          dot += Ax(i,k)*Cx(k,j);
          d(nd) = Cx(k,j);                      // d Yx(i,j) / d Ax(i,k) 
          col(nd++) = AxInd(i,k);
          d(nd) = Ax(i,k);                      // d Yx(i,j) / d Cx(k,j) 
          col(nd++) = CxInd(k,j);
        }
        Yx(i,j) = dot;
        dot = 0;
        for (std::size_t k=Br(i) ; k<Br(i+1) ; k++)
        {
          dot += Bv(k)*Cx(K1+Bc(k), j);
          d(nd) = Bv(k);                        // d Yx(i,j) / d Cx(nxs(t)+Bc(k), j) 
          col(nd++) = CxInd(K1+Bc(k), j);
        }
        Yx(i,j) += dot;
      }
  }
protected:
  virtual void countSize(std::size_t& nIn, std::size_t& nOut)
  {
    nIn = Ax.size1() * Ax.size2() + Cx.size1() * Cx.size2();
    nOut = Yx.size1() * Yx.size2();
  }

  virtual AutoDiff<T>* createLocalObj(std::size_t& nc)
  { // Deallocated in freeLocalObj
    matrix<std::size_t>* AxInd_= new matrix<std::size_t>(Ax.size1(), Ax.size2()); 
    matrix<std::size_t>* CxInd_= new matrix<std::size_t>(Cx.size1(), Cx.size2()); 
    matrix<std::size_t>* YxInd_= new matrix<std::size_t>(Yx.size1(), Yx.size2()); 
    nc = 0;
    for (std::size_t i=0 ; i<Ax.size1() ; i++)
      for (std::size_t j=0 ; j<Ax.size2() ; j++)
        (*AxInd_)(i,j) = nc++;
    for (std::size_t i=0 ; i<Cx.size1() ; i++)
      for (std::size_t j=0 ; j<Cx.size2() ; j++)
        (*CxInd_)(i,j) = nc++;
    for (std::size_t i=0 ; i<Yx.size1() ; i++)
      for (std::size_t j=0 ; j<Yx.size2() ; j++)
        (*YxInd_)(i,j) = nc++;
    std::size_t nIn;
    std::size_t nOut;
    countSize(nIn, nOut);
    if (nc != nIn+nOut)
      throw ExcKal(std::string("AutoDiff: : Incorrect size in createLocalObj"));

    return new AD_AxBs_Cx(&Ax, &Bs, &Cx, &Yx, AxInd_, CxInd_, YxInd_);
  }

  virtual void freeLocalObj()
  { // Free objects created by createLocalObj
    delete &AxInd;
    delete &CxInd;
    delete &YxInd;
  }

  virtual void testADgrad(const vector<T>& dOdx, vector<T>& grad) // Compute gradient with respect to input arguments
  { 
    std::size_t k=0;
    for (std::size_t i=0 ; i<Ax.size1() ; i++)
      for (std::size_t j=0 ; j<Ax.size2() ; j++)
        grad(k++) = dOdx(AxInd(i,j));
    for (std::size_t i=0 ; i<Cx.size1() ; i++)
      for (std::size_t j=0 ; j<Cx.size2() ; j++)
        grad(k++) = dOdx(CxInd(i,j));
    if (k != grad.size())
      throw ExcKal(std::string("AutoDiff: : Incorrect size in grad"));
  }

  virtual void gradNum(const T epsilon, matrix<T>& Jacobian)
  {
    std::size_t nIn;
    std::size_t nOut;
    countSize(nIn, nOut);
    Jacobian.resize(nIn, nOut);
    f();
    matrix<T> Ax_(Ax);
    matrix<T> Cx_(Cx);
    matrix<T> Yx_(Yx);

    std::size_t ii=0;
    for (std::size_t i=0 ; i<Ax.size1() ; i++)
      for (std::size_t j=0 ; j<Ax.size2() ; j++)
      {        
        Ax_(i,j) += epsilon;
        AD_AxBs_Cx<T> ad(&Ax_, &Bs, &Cx_, &Yx_, &AxInd, &CxInd, &YxInd);
        ad.f();
        std::size_t jj=0;
        for (std::size_t iii=0 ; iii<Yx.size1() ; iii++)
          for (std::size_t jjj=0 ; jjj<Yx.size2() ; jjj++)
            Jacobian(ii, jj++) = (Yx_(iii,jjj) - Yx(iii,jjj))/epsilon;
        ii++;
        Ax_(i,j) = Ax(i,j);
        if (jj != nOut)
          throw ExcKal(std::string("AutoDiff: : Incorrect size in gradNum"));
      }
    for (std::size_t i=0 ; i<Cx.size1() ; i++)
      for (std::size_t j=0 ; j<Cx.size2() ; j++)
      {        
        Cx_(i,j) += epsilon;
        AD_AxBs_Cx<T> ad(&Ax_, &Bs, &Cx_, &Yx_, &AxInd, &CxInd, &YxInd);
        ad.f();
        std::size_t jj=0;
        for (std::size_t iii=0 ; iii<Yx.size1() ; iii++)
          for (std::size_t jjj=0 ; jjj<Yx.size2() ; jjj++)
            Jacobian(ii, jj++) = (Yx_(iii,jjj) - Yx(iii,jjj))/epsilon;
        ii++;
        Cx_(i,j) = Cx(i,j);
        if (jj != nOut)
          throw ExcKal(std::string("AutoDiff: : Incorrect size in gradNum"));
      }
    if (ii != nIn)
      throw ExcKal(std::string("AutoDiff: : Incorrect size in gradNum"));
  }

  const matrix<T>& Ax;
  const SparseMatrix<T>& Bs;
  const matrix<T>& Cx;
  matrix<T>& Yx;
  const matrix<std::size_t>& AxInd;
  const matrix<std::size_t>& CxInd;
  const matrix<std::size_t>& YxInd;
};

template <class T>
class AD_Ax_BxCsT : public AutoDiff<T>
{ // Y(x) = A(x) * [B(x) Cs]^T
public:
  AD_Ax_BxCsT(const matrix<T>* Axp, const matrix<T>* Bxp, const SparseMatrix<T>* Csp, matrix<T>* Yxp, 
    const matrix<std::size_t>* AxIndp, const matrix<std::size_t>* BxIndp, const matrix<std::size_t>* YxIndp) 
    : Ax(*Axp), Bx(*Bxp), Cs(*Csp), Yx(*Yxp), AxInd(*AxIndp), BxInd(*BxIndp), YxInd(*YxIndp) 
  {
    // Check dimensions
    const std::size_t m = Ax.size1();
    const std::size_t K1 = Bx.size2();
    const std::size_t K2 = Cs.size2();
    const std::size_t n = Bx.size1();
    const std::string strFunc("AD_Ax_BxCsT");
    checkDim(strFunc, "Ax", Ax, m, K1+K2);
    checkDim(strFunc, "Bx", Bx, n, K1);
    checkDim(strFunc, "Cs", Cs, n, K2);
    checkDim(strFunc, "Yx", Yx, m, n);
    checkDim(strFunc, "AxInd", AxInd, m, K1+K2);
    checkDim(strFunc, "BxInd", BxInd, n, K1);
    checkDim(strFunc, "YxInd", YxInd, m, n);
  }
  virtual ~AD_Ax_BxCsT() {}

  virtual void count(std::size_t& rMax, std::size_t& dMax) const 
  {
    const std::size_t m = Ax.size1();
    const std::size_t K1 = Bx.size2();
    const std::size_t n = Bx.size1();

    const vector<T>& Cv = Cs.val();

    rMax += m*n; 
    dMax += 2*m*n*K1 + m*Cv.size();
  }
  virtual void f() 
  {
    const std::size_t m = Ax.size1();
    const std::size_t K1 = Bx.size2();
    const std::size_t K2 = Cs.size2();
    const std::size_t n = Bx.size1();

    const vector<std::size_t>& Cr = Cs.row();
    const vector<std::size_t>& Cc = Cs.col();
    const vector<T>& Cv = Cs.val();

    for (std::size_t i=0 ; i<m ; i++)
      for (std::size_t j=0 ; j<n ; j++)
      {
        T dot = 0;
        for (std::size_t k=0 ; k<K1 ; k++)
          dot += Ax(i,k)*Bx(j,k);
        Yx(i,j) = dot;
        dot = 0;
        for (std::size_t k=Cr(j) ; k<Cr(j+1) ; k++)
          dot += Ax(i, K1+Cc(k))*Cv(k);
        Yx(i,j) += dot;
      }
  }
  virtual void fAD(vector<T>& d, vector<std::size_t>& col, vector<std::size_t>& row, vector<std::size_t>& rInd, std::size_t& nr, std::size_t& nd) 
  {
    const std::size_t m = Ax.size1();
    const std::size_t K1 = Bx.size2();
    const std::size_t K2 = Cs.size2();
    const std::size_t n = Bx.size1();

    const vector<std::size_t>& Cr = Cs.row();
    const vector<std::size_t>& Cc = Cs.col();
    const vector<T>& Cv = Cs.val();

    for (std::size_t i=0 ; i<m ; i++)
      for (std::size_t j=0 ; j<n ; j++)
      {
        rInd(nr) = YxInd(i,j);
        row(nr++) = nd;
        T dot = 0;
        for (std::size_t k=0 ; k<K1 ; k++)
        {
          dot += Ax(i,k)*Bx(j,k);
          d(nd) = Ax(i,k);                      // d Yx(i,j) / d Bx(j,k) 
          col(nd++) = BxInd(j,k);
          d(nd) = Bx(j,k);                // d Yx(i,j) / d Ax(i,k) 
          col(nd++) = AxInd(i,k);
        }
        Yx(i,j) = dot;
        dot = 0;
        for (std::size_t k=Cr(j) ; k<Cr(j+1) ; k++)
        {
          dot += Ax(i, K1+Cc(k))*Cv(k);
          d(nd) = Cv(k);                       // d Yx(i,j) / d Ax(i, K1+Cc(k)) 
          col(nd++) = AxInd(i, K1+Cc(k));
        }
        Yx(i,j) += dot;
      }
  }
protected:
  virtual void countSize(std::size_t& nIn, std::size_t& nOut)
  {
    nIn = Ax.size1() * Ax.size2() + Bx.size1() * Bx.size2();
    nOut = Yx.size1() * Yx.size2();
  }

  virtual AutoDiff<T>* createLocalObj(std::size_t& nc)
  { // Deallocated in freeLocalObj
    matrix<std::size_t>* AxInd_= new matrix<std::size_t>(Ax.size1(), Ax.size2()); 
    matrix<std::size_t>* BxInd_= new matrix<std::size_t>(Bx.size1(), Bx.size2()); 
    matrix<std::size_t>* YxInd_= new matrix<std::size_t>(Yx.size1(), Yx.size2()); 
    nc = 0;
    for (std::size_t i=0 ; i<Ax.size1() ; i++)
      for (std::size_t j=0 ; j<Ax.size2() ; j++)
        (*AxInd_)(i,j) = nc++;
    for (std::size_t i=0 ; i<Bx.size1() ; i++)
      for (std::size_t j=0 ; j<Bx.size2() ; j++)
        (*BxInd_)(i,j) = nc++;
    for (std::size_t i=0 ; i<Yx.size1() ; i++)
      for (std::size_t j=0 ; j<Yx.size2() ; j++)
        (*YxInd_)(i,j) = nc++;
    std::size_t nIn;
    std::size_t nOut;
    countSize(nIn, nOut);
    if (nc != nIn+nOut)
      throw ExcKal(std::string("AutoDiff: : Incorrect size in createLocalObj"));

    return new AD_Ax_BxCsT(&Ax, &Bx, &Cs, &Yx, AxInd_, BxInd_, YxInd_);
  }

  virtual void freeLocalObj()
  { // Free objects created by createLocalObj
    delete &AxInd;
    delete &BxInd;
    delete &YxInd;
  }

  virtual void testADgrad(const vector<T>& dOdx, vector<T>& grad) // Compute gradient with respect to input arguments
  { 
    std::size_t k=0;
    for (std::size_t i=0 ; i<Ax.size1() ; i++)
      for (std::size_t j=0 ; j<Ax.size2() ; j++)
        grad(k++) = dOdx(AxInd(i,j));
    for (std::size_t i=0 ; i<Bx.size1() ; i++)
      for (std::size_t j=0 ; j<Bx.size2() ; j++)
        grad(k++) = dOdx(BxInd(i,j));
    if (k != grad.size())
      throw ExcKal(std::string("AutoDiff: : Incorrect size in grad"));
  }

  virtual void gradNum(const T epsilon, matrix<T>& Jacobian)
  {
    std::size_t nIn;
    std::size_t nOut;
    countSize(nIn, nOut);
    Jacobian.resize(nIn, nOut);
    f();
    matrix<T> Ax_(Ax);
    matrix<T> Bx_(Bx);
    matrix<T> Yx_(Yx);

    std::size_t ii=0;
    for (std::size_t i=0 ; i<Ax.size1() ; i++)
      for (std::size_t j=0 ; j<Ax.size2() ; j++)
      {        
        Ax_(i,j) += epsilon;
        AD_Ax_BxCsT<T> ad(&Ax_, &Bx_, &Cs, &Yx_, &AxInd, &BxInd, &YxInd);
        ad.f();
        std::size_t jj=0;
        for (std::size_t iii=0 ; iii<Yx.size1() ; iii++)
          for (std::size_t jjj=0 ; jjj<Yx.size2() ; jjj++)
            Jacobian(ii, jj++) = (Yx_(iii,jjj) - Yx(iii,jjj))/epsilon;
        ii++;
        Ax_(i,j) = Ax(i,j);
        if (jj != nOut)
          throw ExcKal(std::string("AutoDiff: : Incorrect size in gradNum"));
      }
    for (std::size_t i=0 ; i<Bx.size1() ; i++)
      for (std::size_t j=0 ; j<Bx.size2() ; j++)
      {        
        Bx_(i,j) += epsilon;
        AD_Ax_BxCsT<T> ad(&Ax_, &Bx_, &Cs, &Yx_, &AxInd, &BxInd, &YxInd);
        ad.f();
        std::size_t jj=0;
        for (std::size_t iii=0 ; iii<Yx.size1() ; iii++)
          for (std::size_t jjj=0 ; jjj<Yx.size2() ; jjj++)
            Jacobian(ii, jj++) = (Yx_(iii,jjj) - Yx(iii,jjj))/epsilon;
        ii++;
        Bx_(i,j) = Bx(i,j);
        if (jj != nOut)
          throw ExcKal(std::string("AutoDiff: : Incorrect size in gradNum"));
      }
    if (ii != nIn)
      throw ExcKal(std::string("AutoDiff: : Incorrect size in gradNum"));
  }

  const matrix<T>& Ax;
  const matrix<T>& Bx;
  const SparseMatrix<T>& Cs;
  matrix<T>& Yx;
  const matrix<std::size_t>& AxInd;
  const matrix<std::size_t>& BxInd;
  const matrix<std::size_t>& YxInd;
};

template <class T>
class AD_cholesky : public AutoDiff<T>
{ // L(x) = cholesky(A(x))
public:
  AD_cholesky(const matrix<T>* Axp, matrix<T>* Lxp, const matrix<std::size_t>* AxIndp, const matrix<std::size_t>* LxIndp) 
    : Ax(*Axp), Lx(*Lxp), AxInd(*AxIndp), LxInd(*LxIndp) 
  {
    // Check dimensions
    const std::size_t m = Ax.size1();
    const std::size_t n = Ax.size2();
    const std::string strFunc("AD_cholesky");
    if (m != n)
      throw ExcKal(std::string("AutoDiff:") + strFunc + std::string(" has to be a square matrix "));
    checkDim(strFunc, "Ax", Ax, m, n);
    checkDim(strFunc, "Lx", Lx, m, n);
    checkDim(strFunc, "AxInd", AxInd, m, n);
    checkDim(strFunc, "LxInd", LxInd, m, n);
  }
  virtual ~AD_cholesky() {}

  virtual void count(std::size_t& rMax, std::size_t& dMax) const 
  {
    const std::size_t m = Ax.size1();

    for (std::size_t j=0 ; j<m ; j++)
    {  
      rMax += m-j; 
      dMax += j+1 + 2*(m-(j+1)) + 2*(m-(j+1))*j; 
    }
  }
  virtual void f() 
  {
    const std::size_t m = Ax.size1();

    for (std::size_t j=0 ; j<m ; j++)
      for (std::size_t i=0 ; i<j ; i++)
			  Lx(i,j) = 0.0;

    for (std::size_t j=0 ; j<m ; j++)                // form column j of Lx
    {
      T dot = 0;

      for (std::size_t i=0 ; i<j ; i++)
        dot += Lx(j,i)*Lx(j,i);

      /*if (Ax(j,j) - dot <= 0.0) 
      {
	      //std::cerr << "Fixing Cholesky " << Ax(j,j) - dot << std::endl;
        Lx(j,j) = sqrt(dot - Ax(j,j) + 1);
      }
      else*/
        Lx(j,j) = sqrt( Ax(j,j) - dot );

      for (std::size_t i=j+1 ; i<m ; i++)
      {
        T dot = 0;
        for (std::size_t k=0 ; k<j ; k++)
	        dot += Lx(i,k)*Lx(j,k);
        Lx(i,j) = (Ax(j,i) - dot)/Lx(j,j);
      }
    }
  }
  virtual void fAD(vector<T>& d, vector<std::size_t>& col, vector<std::size_t>& row, vector<std::size_t>& rInd, std::size_t& nr, std::size_t& nd) 
  {
    const std::size_t m = Ax.size1();

    for (std::size_t j=0 ; j<m ; j++)
      for (std::size_t i=0 ; i<j ; i++)
			  Lx(i,j) = 0.0;

    for (std::size_t j=0 ; j<m ; j++)                // form column j of Lx
    {
      T dot = 0;

      for (std::size_t i=0 ; i<j ; i++)
        dot += Lx(j,i)*Lx(j,i);

      rInd(nr) = LxInd(j,j);
      row(nr++) = nd;
      /*if (Ax(j,j) - dot <= 0.0) 
      {
	      //std::cerr << "Fixing Cholesky " << Ax(j,j) - dot << std::endl;
        Lx(j,j) = sqrt(dot - Ax(j,j) + 1);
        d(nd) = -0.5/Lx(j,j);                    // d Lx(j,j) / d Ax(j,j) 
        col(nd++) = AxInd(j,j);
        for (std::size_t i=0 ; i<j ; i++)
        {
          d(nd) = Lx(j,i)/Lx(j,j);                // d Lx(j,j) / d Lx(j,i) 
          col(nd++) = LxInd(j,i);
        }
      }
      else*/
      {
        Lx(j,j) = sqrt( Ax(j,j) - dot );
        d(nd) = 0.5/Lx(j,j);                     // d Lx(j,j) / d Ax(j,j) 
        col(nd++) = AxInd(j,j);
        for (std::size_t i=0 ; i<j ; i++)
        {
          d(nd) = -Lx(j,i)/Lx(j,j);               // d Lx(j,j) / d Lx(j,i) 
          col(nd++) = LxInd(j,i);
        }
      }

      for (std::size_t i=j+1 ; i<m ; i++)
      {
        rInd(nr) = LxInd(i,j);
        row(nr++) = nd;
        T dot = 0;
        for (std::size_t k=0 ; k<j ; k++)
	        dot += Lx(i,k)*Lx(j,k);
        Lx(i,j) = (Ax(j,i) - dot)/Lx(j,j);

        d(nd) = 1/Lx(j,j);                       // d Lx(i,j) / d Ax(j,i) 
        col(nd++) = AxInd(j,i);
        d(nd) = -Lx(i,j)/Lx(j,j);                 // d Lx(i,j) / d Lx(j,j) 
        col(nd++) = LxInd(j,j);
        for (std::size_t k=0 ; k<j ; k++)
        {
          d(nd) = -Lx(j,k)/Lx(j,j);               // d Lx(i,j) / d Lx(i,k) 
          col(nd++) = LxInd(i,k);
          d(nd) = -Lx(i,k)/Lx(j,j);               // d Lx(i,j) / d Lx(j,k) 
          col(nd++) = LxInd(j,k);
        }
      }
    }
  }
protected:
  virtual void countSize(std::size_t& nIn, std::size_t& nOut)
  {
    nIn = Ax.size1() * Ax.size2();
    nOut = Lx.size1() * Lx.size2();
  }

  virtual AutoDiff<T>* createLocalObj(std::size_t& nc)
  { // Deallocated in freeLocalObj
    matrix<std::size_t>* AxInd_= new matrix<std::size_t>(Ax.size1(), Ax.size2()); 
    matrix<std::size_t>* LxInd_= new matrix<std::size_t>(Lx.size1(), Lx.size2()); 
    nc = 0;
    for (std::size_t i=0 ; i<Ax.size1() ; i++)
      for (std::size_t j=0 ; j<Ax.size2() ; j++)
        (*AxInd_)(i,j) = nc++;
    for (std::size_t i=0 ; i<Lx.size1() ; i++)
      for (std::size_t j=0 ; j<Lx.size2() ; j++)
        (*LxInd_)(i,j) = nc++;
    std::size_t nIn;
    std::size_t nOut;
    countSize(nIn, nOut);
    if (nc != nIn+nOut)
      throw ExcKal(std::string("AutoDiff: : Incorrect size in createLocalObj"));

    return new AD_cholesky(&Ax, &Lx, AxInd_, LxInd_);
  }

  virtual void freeLocalObj()
  { // Free objects created by createLocalObj
    delete &AxInd;
    delete &LxInd;
  }

  virtual void testADgrad(const vector<T>& dOdx, vector<T>& grad) // Compute gradient with respect to input arguments
  { 
    std::size_t k=0;
    for (std::size_t i=0 ; i<Ax.size1() ; i++)
      for (std::size_t j=0 ; j<Ax.size2() ; j++)
        grad(k++) = dOdx(AxInd(i,j));
    if (k != grad.size())
      throw ExcKal(std::string("AutoDiff: : Incorrect size in grad"));
  }

  virtual void gradNum(const T epsilon, matrix<T>& Jacobian)
  {
    std::size_t nIn;
    std::size_t nOut;
    countSize(nIn, nOut);
    Jacobian.resize(nIn, nOut);
    f();
    matrix<T> Ax_(Ax);
    matrix<T> Lx_(Lx);

    std::size_t ii=0;
    for (std::size_t i=0 ; i<Ax.size1() ; i++)
      for (std::size_t j=0 ; j<Ax.size2() ; j++)
      {        
        Ax_(i,j) += epsilon;
        AD_cholesky<T> ad(&Ax_, &Lx_, &AxInd, &LxInd);
        ad.f();
        std::size_t jj=0;
        for (std::size_t iii=0 ; iii<Lx.size1() ; iii++)
          for (std::size_t jjj=0 ; jjj<Lx.size2() ; jjj++)
            Jacobian(ii, jj++) = (Lx_(iii,jjj) - Lx(iii,jjj))/epsilon;
        ii++;
        Ax_(i,j) = Ax(i,j);
        if (jj != nOut)
          throw ExcKal(std::string("AutoDiff: : Incorrect size in gradNum"));
      }
    if (ii != nIn)
      throw ExcKal(std::string("AutoDiff: : Incorrect size in gradNum"));
  }

  const matrix<T>& Ax;
  matrix<T>& Lx;
  const matrix<std::size_t>& AxInd;
  const matrix<std::size_t>& LxInd;
};

template <class T>
class AD_eqSolveTriang1_LxLxTyx_ax : public AutoDiff<T>
{ // Solve triangular equation system L(x) L(x)^T y(x) = a(x) - step 1 

public:
  AD_eqSolveTriang1_LxLxTyx_ax(const vector<T>* axp, const matrix<T>* Lxp, vector<T>* yxp, const vector<std::size_t>* axIndp, const matrix<std::size_t>* LxIndp, const vector<std::size_t>* yxIndp) 
    : ax(*axp), Lx(*Lxp), yx(*yxp), axInd(*axIndp), LxInd(*LxIndp), yxInd(*yxIndp) 
  {
    // Check dimensions
    const std::size_t m = Lx.size1();
    const std::string strFunc("AD_eqSolveTriang1_LxLxTyx_ax");
    checkDim(strFunc, "ax", ax, m);
    checkDim(strFunc, "Lx", Lx, m, m);
    checkDim(strFunc, "yx", yx, m);
    checkDim(strFunc, "axInd", axInd, m);
    checkDim(strFunc, "LxInd", LxInd, m, m);
    checkDim(strFunc, "yxInd", yxInd, m);
  }
  virtual ~AD_eqSolveTriang1_LxLxTyx_ax() {}

  virtual void count(std::size_t& rMax, std::size_t& dMax) const 
  {
    const std::size_t m = Lx.size1();
    for (std::size_t i=0 ; i<m ; i++)
    {  
      rMax += 1; 
      dMax += 2 + 2*i; 
    } 
  }
  virtual void f() 
  {
    const std::size_t m = Lx.size1();

    for (std::size_t i=0 ; i<m ; i++)
	  {
		  T ytmp=0;
		  for (std::size_t j=0; j<i; j++)
			  ytmp += Lx(i,j)*yx(j);
		  yx(i) = (ax(i) - ytmp)/Lx(i,i);
    }
  }
  virtual void fAD(vector<T>& d, vector<std::size_t>& col, vector<std::size_t>& row, vector<std::size_t>& rInd, std::size_t& nr, std::size_t& nd) 
  {
    const std::size_t m = Lx.size1();

    for (std::size_t i=0 ; i<m ; i++)
	  {
      rInd(nr) = yxInd(i);
      row(nr++) = nd;
		  T ytmp=0;
		  for (std::size_t j=0; j<i; j++)
			  ytmp += Lx(i,j)*yx(j);
		  yx(i) = (ax(i) - ytmp)/Lx(i,i);

      d(nd) = 1/Lx(i,i);                         // d yx(i) / d ax(i) 
      col(nd++) = axInd(i);
      d(nd) = -yx(i)/Lx(i,i);          // d yx(i) / d Lx(i,i) 
      col(nd++) = LxInd(i,i);
		  for (std::size_t j=0; j<i; j++)
		  {
        d(nd) = -yx(j)/Lx(i,i);        // d yx(i) / d Lx(i,j)
        col(nd++) = LxInd(i,j);
        d(nd) = -Lx(i,j)/Lx(i,i);                 // d yx(i) / d yx(j)
        col(nd++) = yxInd(j);
		  }
    }
  }
protected:
  virtual void countSize(std::size_t& nIn, std::size_t& nOut)
  {
    nIn = ax.size() + Lx.size1() * Lx.size2();
    nOut = yx.size();
  }

  virtual AutoDiff<T>* createLocalObj(std::size_t& nc)
  { // Deallocated in freeLocalObj
    vector<std::size_t>* axInd_= new vector<std::size_t>(ax.size()); 
    matrix<std::size_t>* LxInd_= new matrix<std::size_t>(Lx.size1(), Lx.size2()); 
    vector<std::size_t>* yxInd_= new vector<std::size_t>(yx.size()); 
    nc = 0;
    for (std::size_t i=0 ; i<ax.size() ; i++)
      (*axInd_)(i) = nc++;
    for (std::size_t i=0 ; i<Lx.size1() ; i++)
      for (std::size_t j=0 ; j<Lx.size2() ; j++)
        (*LxInd_)(i,j) = nc++;
    for (std::size_t i=0 ; i<yx.size() ; i++)
      (*yxInd_)(i) = nc++;
    std::size_t nIn;
    std::size_t nOut;
    countSize(nIn, nOut);
    if (nc != nIn+nOut)
      throw ExcKal(std::string("AutoDiff: : Incorrect size in createLocalObj"));

    return new AD_eqSolveTriang1_LxLxTyx_ax(&ax, &Lx, &yx, axInd_, LxInd_, yxInd_);
  }

  virtual void freeLocalObj()
  { // Free objects created by createLocalObj
    delete &axInd;
    delete &LxInd;
    delete &yxInd;
  }

  virtual void testADgrad(const vector<T>& dOdx, vector<T>& grad) // Compute gradient with respect to input arguments
  { 
    std::size_t k=0;
    for (std::size_t i=0 ; i<ax.size() ; i++)
      grad(k++) = dOdx(axInd(i));
    for (std::size_t i=0 ; i<Lx.size1() ; i++)
      for (std::size_t j=0 ; j<Lx.size2() ; j++)
        grad(k++) = dOdx(LxInd(i,j));
  }

  virtual void gradNum(const T epsilon, matrix<T>& Jacobian)
  {
    std::size_t nIn;
    std::size_t nOut;
    countSize(nIn, nOut);
    Jacobian.resize(nIn, nOut);
    vector<T> ax_(ax);
    matrix<T> Lx_(Lx);
    vector<T> yx_(yx);

    std::size_t ii=0;
    for (std::size_t i=0 ; i<ax.size() ; i++)
    {        
      ax_(i) += epsilon;
      AD_eqSolveTriang1_LxLxTyx_ax<T> ad(&ax_, &Lx_, &yx_, &axInd, &LxInd, &yxInd);
      ad.f();
      std::size_t jj=0;
      for (std::size_t iii=0 ; iii<yx.size() ; iii++)
        Jacobian(ii, jj++) = (yx_(iii) - yx(iii))/epsilon;
      ii++;
      ax_(i) = ax(i);
      if (jj != nOut)
        throw ExcKal(std::string("AutoDiff: : Incorrect size in gradNum"));
    }
    for (std::size_t i=0 ; i<Lx.size1() ; i++)
      for (std::size_t j=0 ; j<Lx.size2() ; j++)
      {        
        Lx_(i,j) += epsilon;
        AD_eqSolveTriang1_LxLxTyx_ax<T> ad(&ax_, &Lx_, &yx_, &axInd, &LxInd, &yxInd);
        ad.f();
        std::size_t jj=0;
        for (std::size_t iii=0 ; iii<yx.size() ; iii++)
          Jacobian(ii, jj++) = (yx_(iii) - yx(iii))/epsilon;
        ii++;
        Lx_(i,j) = Lx(i,j);
        if (jj != nOut)
          throw ExcKal(std::string("AutoDiff: : Incorrect size in gradNum"));
      }
    if (ii != nIn)
      throw ExcKal(std::string("AutoDiff: : Incorrect size in gradNum"));
  }

  const vector<T>& ax;
  const matrix<T>& Lx;
  vector<T>& yx;
  const vector<std::size_t>& axInd;
  const matrix<std::size_t>& LxInd;
  const vector<std::size_t>& yxInd;
};

template <class T>
class AD_eqSolveTriang2_LxLxTyx_ax : public AutoDiff<T>
{ // Solve triangular equation system L(x) L(x)^T y(x) = a(x) - step 2 

public:
  AD_eqSolveTriang2_LxLxTyx_ax(const vector<T>* axp, const matrix<T>* Lxp, vector<T>* yxp, const vector<std::size_t>* axIndp, const matrix<std::size_t>* LxIndp, const vector<std::size_t>* yxIndp) 
    : ax(*axp), Lx(*Lxp), yx(*yxp), axInd(*axIndp), LxInd(*LxIndp), yxInd(*yxIndp) 
  {
    // Check dimensions
    const std::size_t m = Lx.size1();
    const std::string strFunc("AD_eqSolveTriang2_LxLxTyx_ax");
    checkDim(strFunc, "ax", ax, m);
    checkDim(strFunc, "Lx", Lx, m, m);
    checkDim(strFunc, "yx", yx, m);
    checkDim(strFunc, "axInd", axInd, m);
    checkDim(strFunc, "LxInd", LxInd, m, m);
    checkDim(strFunc, "yxInd", yxInd, m);
  }
  virtual ~AD_eqSolveTriang2_LxLxTyx_ax() {}

  virtual void count(std::size_t& rMax, std::size_t& dMax) const 
  {
    const std::size_t m = Lx.size1();
    for (std::size_t i=0 ; i<m ; i++)
    {  
      rMax += 1; 
      dMax += 2 + 2*(m-(i+1));
    } 
  }
  virtual void f() 
  {
    const std::size_t m = Lx.size1();

    for (std::size_t i=m ; i-->0 ; )
	  {
		  T tmp=0;
		  for (std::size_t j=i+1; j<m; j++)
			  tmp += Lx(j,i)*yx(j);
		  yx(i) = (ax(i) - tmp)/Lx(i,i);
    }
  }
  virtual void fAD(vector<T>& d, vector<std::size_t>& col, vector<std::size_t>& row, vector<std::size_t>& rInd, std::size_t& nr, std::size_t& nd) 
  {
    const std::size_t m = Lx.size1();

    for (std::size_t i=m ; i-->0 ; )
	  {
      rInd(nr) = yxInd(i);
      row(nr++) = nd;
		  T tmp=0;
		  for (std::size_t j=i+1; j<m; j++)
			  tmp += Lx(j,i)*yx(j);
		  yx(i) = (ax(i) - tmp)/Lx(i,i);

      d(nd) = 1/Lx(i,i);                         // d yx(i) / d ax(i) 
      col(nd++) = axInd(i);
      d(nd) = -yx(i)/Lx(i,i);           // d yx(i) / d Lx(i,i) 
      col(nd++) = LxInd(i,i);        
      for (std::size_t j=i+1; j<m; j++)
		  {
        d(nd) = -yx(j)/Lx(i,i);         // d yx(i) / d Lx(j,i)
        col(nd++) = LxInd(j,i);
        d(nd) = -Lx(j,i)/Lx(i,i);                 // d yx(i) / d yx(j)
        col(nd++) = yxInd(j);
		  }
    }
  }
protected:
  virtual void countSize(std::size_t& nIn, std::size_t& nOut)
  {
    nIn = ax.size() + Lx.size1() * Lx.size2();
    nOut = yx.size();
  }

  virtual AutoDiff<T>* createLocalObj(std::size_t& nc)
  { // Deallocated in freeLocalObj
    vector<std::size_t>* axInd_= new vector<std::size_t>(ax.size()); 
    matrix<std::size_t>* LxInd_= new matrix<std::size_t>(Lx.size1(), Lx.size2()); 
    vector<std::size_t>* yxInd_= new vector<std::size_t>(yx.size()); 
    nc = 0;
    for (std::size_t i=0 ; i<ax.size() ; i++)
      (*axInd_)(i) = nc++;
    for (std::size_t i=0 ; i<Lx.size1() ; i++)
      for (std::size_t j=0 ; j<Lx.size2() ; j++)
        (*LxInd_)(i,j) = nc++;
    for (std::size_t i=0 ; i<yx.size() ; i++)
      (*yxInd_)(i) = nc++;
    std::size_t nIn;
    std::size_t nOut;
    countSize(nIn, nOut);
    if (nc != nIn+nOut)
      throw ExcKal(std::string("AutoDiff: : Incorrect size in createLocalObj"));

    return new AD_eqSolveTriang2_LxLxTyx_ax(&ax, &Lx, &yx, axInd_, LxInd_, yxInd_);
  }

  virtual void freeLocalObj()
  { // Free objects created by createLocalObj
    delete &axInd;
    delete &LxInd;
    delete &yxInd;
  }

  virtual void testADgrad(const vector<T>& dOdx, vector<T>& grad) // Compute gradient with respect to input arguments
  { 
    std::size_t k=0;
    for (std::size_t i=0 ; i<ax.size() ; i++)
      grad(k++) = dOdx(axInd(i));
    for (std::size_t i=0 ; i<Lx.size1() ; i++)
      for (std::size_t j=0 ; j<Lx.size2() ; j++)
        grad(k++) = dOdx(LxInd(i,j));
  }

  virtual void gradNum(const T epsilon, matrix<T>& Jacobian)
  {
    std::size_t nIn;
    std::size_t nOut;
    countSize(nIn, nOut);
    Jacobian.resize(nIn, nOut);
    vector<T> ax_(ax);
    matrix<T> Lx_(Lx);
    vector<T> yx_(yx);

    std::size_t ii=0;
    for (std::size_t i=0 ; i<ax.size() ; i++)
    {        
      ax_(i) += epsilon;
      AD_eqSolveTriang2_LxLxTyx_ax<T> ad(&ax_, &Lx_, &yx_, &axInd, &LxInd, &yxInd);
      ad.f();
      std::size_t jj=0;
      for (std::size_t iii=0 ; iii<yx.size() ; iii++)
        Jacobian(ii, jj++) = (yx_(iii) - yx(iii))/epsilon;
      ii++;
      ax_(i) = ax(i);
      if (jj != nOut)
        throw ExcKal(std::string("AutoDiff: : Incorrect size in gradNum"));
    }
    for (std::size_t i=0 ; i<Lx.size1() ; i++)
      for (std::size_t j=0 ; j<Lx.size2() ; j++)
      {        
        Lx_(i,j) += epsilon;
        AD_eqSolveTriang2_LxLxTyx_ax<T> ad(&ax_, &Lx_, &yx_, &axInd, &LxInd, &yxInd);
        ad.f();
        std::size_t jj=0;
        for (std::size_t iii=0 ; iii<yx.size() ; iii++)
          Jacobian(ii, jj++) = (yx_(iii) - yx(iii))/epsilon;
        ii++;
        Lx_(i,j) = Lx(i,j);
        if (jj != nOut)
          throw ExcKal(std::string("AutoDiff: : Incorrect size in gradNum"));
      }
    if (ii != nIn)
      throw ExcKal(std::string("AutoDiff: : Incorrect size in gradNum"));
  }

  const vector<T>& ax;
  const matrix<T>& Lx;
  vector<T>& yx;
  const vector<std::size_t>& axInd;
  const matrix<std::size_t>& LxInd;
  const vector<std::size_t>& yxInd;
};


template <class T>
class AD_eqSolveTriang1_YxLxLxT_Ax : public AutoDiff<T>
{ // Solve triangular equation system Y(x) L(x) L(x)^T = A(x) - step 1 (special implementation to skip transpose)

public:
  AD_eqSolveTriang1_YxLxLxT_Ax(const matrix<T>* Axp, const matrix<T>* Lxp, matrix<T>* Yxp, const matrix<std::size_t>* AxIndp, const matrix<std::size_t>* LxIndp, const matrix<std::size_t>* YxIndp) 
    : Ax(*Axp), Lx(*Lxp), Yx(*Yxp), AxInd(*AxIndp), LxInd(*LxIndp), YxInd(*YxIndp) 
  {
    // Check dimensions
    const std::size_t m = Lx.size1();
    const std::size_t n = Lx.size2();
    const std::size_t K = Ax.size2();
    const std::string strFunc("AD_eqSolveTriang1_YxLxLxT_Ax");
    if (m != n)
      throw ExcKal(std::string("AutoDiff:") + strFunc + std::string(" has to be a square matrix "));
    checkDim(strFunc, "Ax", Ax, m, K);
    checkDim(strFunc, "Lx", Lx, m, n);
    checkDim(strFunc, "Yx", Yx, K, m);
    checkDim(strFunc, "AxInd", AxInd, m, K);
    checkDim(strFunc, "LxInd", LxInd, m, n);
    checkDim(strFunc, "YxInd", YxInd, K, m);
  }
  virtual ~AD_eqSolveTriang1_YxLxLxT_Ax() {}

  virtual void count(std::size_t& rMax, std::size_t& dMax) const 
  {
    const std::size_t m = Lx.size1();
    const std::size_t K = Ax.size2();

    for (std::size_t i=0 ; i<m ; i++)
    {  
      rMax += K; 
      dMax += K*2 + K*2*i; 
    } 
  }
  virtual void f() 
  {
    const std::size_t m = Lx.size1();
    const std::size_t K = Ax.size2();
    
    for (std::size_t k=0 ; k<K ; k++)
      for (std::size_t i=0 ; i<m ; i++)
	    {
		    T ytmp=0;
		    for (std::size_t j=0; j<i; j++)
			    ytmp += Lx(i,j)*Yx(k,j);
		    Yx(k,i) = (Ax(i,k) - ytmp)/Lx(i,i);
	    }
  }
  virtual void fAD(vector<T>& d, vector<std::size_t>& col, vector<std::size_t>& row, vector<std::size_t>& rInd, std::size_t& nr, std::size_t& nd) 
  {
    const std::size_t m = Lx.size1();
    const std::size_t K = Ax.size2();
    
    for (std::size_t k=0 ; k<K ; k++)
      for (std::size_t i=0 ; i<m ; i++)
	    {
        rInd(nr) = YxInd(k,i);
        row(nr++) = nd;
		    T ytmp=0;
		    for (std::size_t j=0; j<i; j++)
			    ytmp += Lx(i,j)*Yx(k,j);
		    Yx(k,i) = (Ax(i,k) - ytmp)/Lx(i,i);

        d(nd) = 1/Lx(i,i);                       // d Yx(k,i) / d Ax(i,k) 
        col(nd++) = AxInd(i,k);
        d(nd) = -Yx(k,i)/Lx(i,i);                // d Yx(k,i) / d Lx(i,i) 
        col(nd++) = LxInd(i,i);
		    for (std::size_t j=0; j<i; j++)
		    {
          d(nd) = -Yx(k,j)/Lx(i,i);              // d Yx(k,i) / d Lx(i,j)
          col(nd++) = LxInd(i,j);
          d(nd) = -Lx(i,j)/Lx(i,i);               // d Yx(k,i) / d Yx(k,j)
          col(nd++) = YxInd(k,j);
		    }
	    }
  }
protected:
  virtual void countSize(std::size_t& nIn, std::size_t& nOut)
  {
    nIn = Ax.size1() * Ax.size2() + Lx.size1() * Lx.size2();
    nOut = Yx.size1() * Yx.size2();
  }

  virtual AutoDiff<T>* createLocalObj(std::size_t& nc)
  { // Deallocated in freeLocalObj
    matrix<std::size_t>* AxInd_= new matrix<std::size_t>(Ax.size1(), Ax.size2()); 
    matrix<std::size_t>* LxInd_= new matrix<std::size_t>(Lx.size1(), Lx.size2()); 
    matrix<std::size_t>* YxInd_= new matrix<std::size_t>(Yx.size1(), Yx.size2()); 
    nc = 0;
    for (std::size_t i=0 ; i<Ax.size1() ; i++)
      for (std::size_t j=0 ; j<Ax.size2() ; j++)
        (*AxInd_)(i,j) = nc++;
    for (std::size_t i=0 ; i<Lx.size1() ; i++)
      for (std::size_t j=0 ; j<Lx.size2() ; j++)
        (*LxInd_)(i,j) = nc++;
    for (std::size_t i=0 ; i<Yx.size1() ; i++)
      for (std::size_t j=0 ; j<Yx.size2() ; j++)
        (*YxInd_)(i,j) = nc++;
    std::size_t nIn;
    std::size_t nOut;
    countSize(nIn, nOut);
    if (nc != nIn+nOut)
      throw ExcKal(std::string("AutoDiff: : Incorrect size in createLocalObj"));

    return new AD_eqSolveTriang1_YxLxLxT_Ax(&Ax, &Lx, &Yx, AxInd_, LxInd_, YxInd_);
  }

  virtual void freeLocalObj()
  { // Free objects created by createLocalObj
    delete &AxInd;
    delete &LxInd;
    delete &YxInd;
  }

  virtual void testADgrad(const vector<T>& dOdx, vector<T>& grad) // Compute gradient with respect to input arguments
  { 
    std::size_t k=0;
    for (std::size_t i=0 ; i<Ax.size1() ; i++)
      for (std::size_t j=0 ; j<Ax.size2() ; j++)
        grad(k++) = dOdx(AxInd(i,j));
    for (std::size_t i=0 ; i<Lx.size1() ; i++)
      for (std::size_t j=0 ; j<Lx.size2() ; j++)
        grad(k++) = dOdx(LxInd(i,j));
    if (k != grad.size())
      throw ExcKal(std::string("AutoDiff: : Incorrect size in grad"));
  }

  virtual void gradNum(const T epsilon, matrix<T>& Jacobian)
  {
    std::size_t nIn;
    std::size_t nOut;
    countSize(nIn, nOut);
    Jacobian.resize(nIn, nOut);
    f();
    matrix<T> Ax_(Ax);
    matrix<T> Lx_(Lx);
    matrix<T> Yx_(Yx);

    std::size_t ii=0;
    for (std::size_t i=0 ; i<Ax.size1() ; i++)
      for (std::size_t j=0 ; j<Ax.size2() ; j++)
      {        
        Ax_(i,j) += epsilon;
        AD_eqSolveTriang1_YxLxLxT_Ax<T> ad(&Ax_, &Lx_, &Yx_, &AxInd, &LxInd, &YxInd);
        ad.f();
        std::size_t jj=0;
        for (std::size_t iii=0 ; iii<Yx.size1() ; iii++)
          for (std::size_t jjj=0 ; jjj<Yx.size2() ; jjj++)
            Jacobian(ii, jj++) = (Yx_(iii,jjj) - Yx(iii,jjj))/epsilon;
        ii++;
        Ax_(i,j) = Ax(i,j);
        if (jj != nOut)
          throw ExcKal(std::string("AutoDiff: : Incorrect size in gradNum"));
      }
    for (std::size_t i=0 ; i<Lx.size1() ; i++)
      for (std::size_t j=0 ; j<Lx.size2() ; j++)
      {        
        Lx_(i,j) += epsilon;
        AD_eqSolveTriang1_YxLxLxT_Ax<T> ad(&Ax_, &Lx_, &Yx_, &AxInd, &LxInd, &YxInd);
        ad.f();
        std::size_t jj=0;
        for (std::size_t iii=0 ; iii<Yx.size1() ; iii++)
          for (std::size_t jjj=0 ; jjj<Yx.size2() ; jjj++)
            Jacobian(ii, jj++) = (Yx_(iii,jjj) - Yx(iii,jjj))/epsilon;
        ii++;
        Lx_(i,j) = Lx(i,j);
        if (jj != nOut)
          throw ExcKal(std::string("AutoDiff: : Incorrect size in gradNum"));
      }
    if (ii != nIn)
      throw ExcKal(std::string("AutoDiff: : Incorrect size in gradNum"));
  }

  const matrix<T>& Ax;
  const matrix<T>& Lx;
  matrix<T>& Yx;
  const matrix<std::size_t>& AxInd;
  const matrix<std::size_t>& LxInd;
  const matrix<std::size_t>& YxInd;
};

template <class T>
class AD_eqSolveTriang2_YxLxLxT_Ax : public AutoDiff<T>
{ // Solve triangular equation system Y(x) L(x) L(x)^T = A(x) - step 2 (special implementation to skip transpose)

public:
  AD_eqSolveTriang2_YxLxLxT_Ax(const matrix<T>* Axp, const matrix<T>* Lxp, matrix<T>* Yxp, const matrix<std::size_t>* AxIndp, const matrix<std::size_t>* LxIndp, const matrix<std::size_t>* YxIndp) 
    : Ax(*Axp), Lx(*Lxp), Yx(*Yxp), AxInd(*AxIndp), LxInd(*LxIndp), YxInd(*YxIndp) 
  {
    // Check dimensions
    const std::size_t m = Lx.size1();
    const std::size_t n = Lx.size2();
    const std::size_t K = Ax.size1();
    const std::string strFunc("AD_eqSolveTriang1_YxLxLxT_Ax");
    if (m != n)
      throw ExcKal(std::string("AutoDiff:") + strFunc + std::string(" has to be a square matrix "));
    checkDim(strFunc, "Ax", Ax, K, m);
    checkDim(strFunc, "Lx", Lx, m, n);
    checkDim(strFunc, "Yx", Yx, K, m);
    checkDim(strFunc, "AxInd", AxInd, K, m);
    checkDim(strFunc, "LxInd", LxInd, m, n);
    checkDim(strFunc, "YxInd", YxInd, K, m);
  }
  virtual ~AD_eqSolveTriang2_YxLxLxT_Ax() {}

  virtual void count(std::size_t& rMax, std::size_t& dMax) const 
  {
    const std::size_t m = Lx.size1();
    const std::size_t K = Ax.size1();

    for (std::size_t i=0 ; i<m ; i++)
    {  
      rMax += K; 
      dMax += K*2 + K*2*(m-(i+1)); 
    } 
  }
  virtual void f() 
  {
    const std::size_t m = Lx.size1();
    const std::size_t K = Ax.size1();
    
    for (std::size_t k=0 ; k<K ; k++)
      for (std::size_t i=m ; i-->0 ; )
	    {
		    T tmp=0;
		    for (std::size_t j=i+1 ; j<m ; j++)
			    tmp += Lx(j,i)*Yx(k,j);
		    Yx(k,i) = (Ax(k,i) - tmp)/Lx(i,i);
	    }
  }
  virtual void fAD(vector<T>& d, vector<std::size_t>& col, vector<std::size_t>& row, vector<std::size_t>& rInd, std::size_t& nr, std::size_t& nd) 
  {
    const std::size_t m = Lx.size1();
    const std::size_t K = Ax.size1();
    
    for (std::size_t k=0 ; k<K ; k++)
      for (std::size_t i=m ; i-->0 ; )
	    {
        rInd(nr) = YxInd(k,i);
        row(nr++) = nd;
		    T tmp=0;
		    for (std::size_t j=i+1 ; j<m ; j++)
			    tmp += Lx(j,i)*Yx(k,j);
		    Yx(k,i) = (Ax(k,i) - tmp)/Lx(i,i);

        d(nd) = 1/Lx(i,i);                       // d Yx(k,i) / d Ax(k,i) 
        col(nd++) = AxInd(k,i);
        d(nd) = -Yx(k,i)/Lx(i,i);                 // d Yx(k,i) / d Lx(i,i) 
        col(nd++) = LxInd(i,i);        
        for (std::size_t j=i+1 ; j<m ; j++)
		    {
          d(nd) = -Yx(k,j)/Lx(i,i);               // d Yx(k,i) / d Lx(j,i)
          col(nd++) = LxInd(j,i);
          d(nd) = -Lx(j,i)/Lx(i,i);               // d Yx(k,i) / d Yx(k,j)
          col(nd++) = YxInd(k,j);
		    }
	    }
  }
protected:
  virtual void countSize(std::size_t& nIn, std::size_t& nOut)
  {
    nIn = Ax.size1() * Ax.size2() + Lx.size1() * Lx.size2();
    nOut = Yx.size1() * Yx.size2();
  }

  virtual AutoDiff<T>* createLocalObj(std::size_t& nc)
  { // Deallocated in freeLocalObj
    matrix<std::size_t>* AxInd_= new matrix<std::size_t>(Ax.size1(), Ax.size2()); 
    matrix<std::size_t>* LxInd_= new matrix<std::size_t>(Lx.size1(), Lx.size2()); 
    matrix<std::size_t>* YxInd_= new matrix<std::size_t>(Yx.size1(), Yx.size2()); 
    nc = 0;
    for (std::size_t i=0 ; i<Ax.size1() ; i++)
      for (std::size_t j=0 ; j<Ax.size2() ; j++)
        (*AxInd_)(i,j) = nc++;
    for (std::size_t i=0 ; i<Lx.size1() ; i++)
      for (std::size_t j=0 ; j<Lx.size2() ; j++)
        (*LxInd_)(i,j) = nc++;
    for (std::size_t i=0 ; i<Yx.size1() ; i++)
      for (std::size_t j=0 ; j<Yx.size2() ; j++)
        (*YxInd_)(i,j) = nc++;
    std::size_t nIn;
    std::size_t nOut;
    countSize(nIn, nOut);
    if (nc != nIn+nOut)
      throw ExcKal(std::string("AutoDiff: : Incorrect size in createLocalObj"));

    return new AD_eqSolveTriang2_YxLxLxT_Ax(&Ax, &Lx, &Yx, AxInd_, LxInd_, YxInd_);
  }

  virtual void freeLocalObj()
  { // Free objects created by createLocalObj
    delete &AxInd;
    delete &LxInd;
    delete &YxInd;
  }

  virtual void testADgrad(const vector<T>& dOdx, vector<T>& grad) // Compute gradient with respect to input arguments
  { 
    std::size_t k=0;
    for (std::size_t i=0 ; i<Ax.size1() ; i++)
      for (std::size_t j=0 ; j<Ax.size2() ; j++)
        grad(k++) = dOdx(AxInd(i,j));
    for (std::size_t i=0 ; i<Lx.size1() ; i++)
      for (std::size_t j=0 ; j<Lx.size2() ; j++)
        grad(k++) = dOdx(LxInd(i,j));
    if (k != grad.size())
      throw ExcKal(std::string("AutoDiff: : Incorrect size in grad"));
  }

  virtual void gradNum(const T epsilon, matrix<T>& Jacobian)
  {
    std::size_t nIn;
    std::size_t nOut;
    countSize(nIn, nOut);
    Jacobian.resize(nIn, nOut);
    f();
    matrix<T> Ax_(Ax);
    matrix<T> Lx_(Lx);
    matrix<T> Yx_(Yx);

    std::size_t ii=0;
    for (std::size_t i=0 ; i<Ax.size1() ; i++)
      for (std::size_t j=0 ; j<Ax.size2() ; j++)
      {        
        Ax_(i,j) += epsilon;
        AD_eqSolveTriang2_YxLxLxT_Ax<T> ad(&Ax_, &Lx_, &Yx_, &AxInd, &LxInd, &YxInd);
        ad.f();
        std::size_t jj=0;
        for (std::size_t iii=0 ; iii<Yx.size1() ; iii++)
          for (std::size_t jjj=0 ; jjj<Yx.size2() ; jjj++)
            Jacobian(ii, jj++) = (Yx_(iii,jjj) - Yx(iii,jjj))/epsilon;
        ii++;
        Ax_(i,j) = Ax(i,j);
        if (jj != nOut)
          throw ExcKal(std::string("AutoDiff: : Incorrect size in gradNum"));
      }
    for (std::size_t i=0 ; i<Lx.size1() ; i++)
      for (std::size_t j=0 ; j<Lx.size2() ; j++)
      {        
        Lx_(i,j) += epsilon;
        AD_eqSolveTriang2_YxLxLxT_Ax<T> ad(&Ax_, &Lx_, &Yx_, &AxInd, &LxInd, &YxInd);
        ad.f();
        std::size_t jj=0;
        for (std::size_t iii=0 ; iii<Yx.size1() ; iii++)
          for (std::size_t jjj=0 ; jjj<Yx.size2() ; jjj++)
            Jacobian(ii, jj++) = (Yx_(iii,jjj) - Yx(iii,jjj))/epsilon;
        ii++;
        Lx_(i,j) = Lx(i,j);
        if (jj != nOut)
          throw ExcKal(std::string("AutoDiff: : Incorrect size in gradNum"));
      }
    if (ii != nIn)
      throw ExcKal(std::string("AutoDiff: : Incorrect size in gradNum"));
  }

  const matrix<T>& Ax;
  const matrix<T>& Lx;
  matrix<T>& Yx;
  const matrix<std::size_t>& AxInd;
  const matrix<std::size_t>& LxInd;
  const matrix<std::size_t>& YxInd;
};

template <class T>
class AD_ax_Bx_cx : public AutoDiff<T>
{ // y(x) = a(x) + B(x) * c(x) 
public:
  AD_ax_Bx_cx(const vector<T>* axp, const matrix<T>* Bxp, const vector<T>* cxp, vector<T>* yxp,
               const vector<std::size_t>* axIndp, const matrix<std::size_t>* BxIndp, const vector<std::size_t>* cxIndp, const vector<std::size_t>* yxIndp) 
    : ax(*axp), Bx(*Bxp), cx(*cxp), yx(*yxp), axInd(*axIndp), BxInd(*BxIndp), cxInd(*cxIndp), yxInd(*yxIndp)   
  {
    // Check dimensions
    const std::size_t m = ax.size();
    const std::size_t n = Bx.size2();
    const std::string strFunc("AD_ax_Bx_cx");
    checkDim(strFunc, "ax", ax, m);
    checkDim(strFunc, "Bx", Bx, m, n);
    checkDim(strFunc, "cx", cx, n);
    checkDim(strFunc, "yx", yx, m);
    checkDim(strFunc, "axInd", axInd, m);
    checkDim(strFunc, "BxInd", BxInd, m, n);
    checkDim(strFunc, "cxInd", cxInd, n);
    checkDim(strFunc, "yxInd", yxInd, m);
  }
  virtual ~AD_ax_Bx_cx() {}

  virtual void count(std::size_t& rMax, std::size_t& dMax) const 
  {
    const std::size_t m = ax.size();
    const std::size_t n = Bx.size2();

    rMax += m; 
    dMax += m + 2*m*n;
  }
  virtual void f() 
  {
    const std::size_t m = ax.size();
    const std::size_t n = Bx.size2();

    for (std::size_t i=0 ; i<m ; i++)
    {
      yx(i) = ax(i);
      T dot = 0;
      for (std::size_t j=0 ; j<n ; j++)
        dot += Bx(i,j)*cx(j);
      yx(i) += dot;
    }
  }
  virtual void fAD(vector<T>& d, vector<std::size_t>& col, vector<std::size_t>& row, vector<std::size_t>& rInd, std::size_t& nr, std::size_t& nd) 
  {
    const std::size_t m = ax.size();
    const std::size_t n = Bx.size2();

    for (std::size_t i=0 ; i<m ; i++)
    {
      rInd(nr) = yxInd(i);
      row(nr++) = nd;
      yx(i) = ax(i);
      d(nd) = 1;                                 // d yx(i) / d ax(i) 
      col(nd++) = axInd(i);
      T dot = 0;
      for (std::size_t j=0 ; j<n ; j++)
      {
        dot += Bx(i,j)*cx(j);
        d(nd) = cx(j);                           // d yx(i) / d Bx(i,j) 
        col(nd++) = BxInd(i,j);
        d(nd) = Bx(i,j);                         // d yx(i) / d cx(j) 
        col(nd++) = cxInd(j);
      }
      yx(i) += dot;
    }
  }
protected:
  virtual void countSize(std::size_t& nIn, std::size_t& nOut)
  {
    nIn = ax.size() + Bx.size1() * Bx.size2() + cx.size();
    nOut = yx.size();
  }

  virtual AutoDiff<T>* createLocalObj(std::size_t& nc)
  { // Deallocated in freeLocalObj
    vector<std::size_t>* axInd_= new vector<std::size_t>(ax.size()); 
    matrix<std::size_t>* BxInd_= new matrix<std::size_t>(Bx.size1(), Bx.size2()); 
    vector<std::size_t>* cxInd_= new vector<std::size_t>(cx.size()); 
    vector<std::size_t>* yxInd_= new vector<std::size_t>(yx.size()); 
    nc = 0;
    for (std::size_t i=0 ; i<ax.size() ; i++)
      (*axInd_)(i) = nc++;
    for (std::size_t i=0 ; i<Bx.size1() ; i++)
      for (std::size_t j=0 ; j<Bx.size2() ; j++)
        (*BxInd_)(i,j) = nc++;
    for (std::size_t i=0 ; i<cx.size() ; i++)
      (*cxInd_)(i) = nc++;
    for (std::size_t i=0 ; i<yx.size() ; i++)
      (*yxInd_)(i) = nc++;
    std::size_t nIn;
    std::size_t nOut;
    countSize(nIn, nOut);
    if (nc != nIn+nOut)
      throw ExcKal(std::string("AutoDiff: : Incorrect size in createLocalObj"));

    return new AD_ax_Bx_cx(&ax, &Bx, &cx, &yx, axInd_, BxInd_, cxInd_, yxInd_);
  }

  virtual void freeLocalObj()
  { // Free objects created by createLocalObj
    delete &axInd;
    delete &BxInd;
    delete &cxInd;
    delete &yxInd;
  }

  virtual void testADgrad(const vector<T>& dOdx, vector<T>& grad) // Compute gradient with respect to input arguments
  { 
    std::size_t k=0;
    for (std::size_t i=0 ; i<ax.size() ; i++)
      grad(k++) = dOdx(axInd(i));
    for (std::size_t i=0 ; i<Bx.size1() ; i++)
      for (std::size_t j=0 ; j<Bx.size2() ; j++)
        grad(k++) = dOdx(BxInd(i,j));
    for (std::size_t i=0 ; i<cx.size() ; i++)
      grad(k++) = dOdx(cxInd(i));
  }

  virtual void gradNum(const T epsilon, matrix<T>& Jacobian)
  {
    std::size_t nIn;
    std::size_t nOut;
    countSize(nIn, nOut);
    Jacobian.resize(nIn, nOut);
    vector<T> ax_(ax);
    matrix<T> Bx_(Bx);
    vector<T> cx_(cx);
    vector<T> yx_(yx);

    std::size_t ii=0;
    for (std::size_t i=0 ; i<ax.size() ; i++)
    {        
      ax_(i) += epsilon;
      AD_ax_Bx_cx<T> ad(&ax_, &Bx_, &cx_, &yx_, &axInd, &BxInd, &cxInd, &yxInd);
      ad.f();
      std::size_t jj=0;
      for (std::size_t iii=0 ; iii<yx.size() ; iii++)
        Jacobian(ii, jj++) = (yx_(iii) - yx(iii))/epsilon;
      ii++;
      ax_(i) = ax(i);
      if (jj != nOut)
        throw ExcKal(std::string("AutoDiff: : Incorrect size in gradNum"));
    }
    for (std::size_t i=0 ; i<Bx.size1() ; i++)
      for (std::size_t j=0 ; j<Bx.size2() ; j++)
      {        
        Bx_(i,j) += epsilon;
        AD_ax_Bx_cx<T> ad(&ax_, &Bx_, &cx_, &yx_, &axInd, &BxInd, &cxInd, &yxInd);
        ad.f();
        std::size_t jj=0;
        for (std::size_t iii=0 ; iii<yx.size() ; iii++)
          Jacobian(ii, jj++) = (yx_(iii) - yx(iii))/epsilon;
        ii++;
        Bx_(i,j) = Bx(i,j);
        if (jj != nOut)
          throw ExcKal(std::string("AutoDiff: : Incorrect size in gradNum"));
      }
    for (std::size_t i=0 ; i<cx.size() ; i++)
    {        
      cx_(i) += epsilon;
      AD_ax_Bx_cx<T> ad(&ax_, &Bx_, &cx_, &yx_, &axInd, &BxInd, &cxInd, &yxInd);
      ad.f();
      std::size_t jj=0;
      for (std::size_t iii=0 ; iii<yx.size() ; iii++)
        Jacobian(ii, jj++) = (yx_(iii) - yx(iii))/epsilon;
      ii++;
      cx_(i) = cx(i);
      if (jj != nOut)
        throw ExcKal(std::string("AutoDiff: : Incorrect size in gradNum"));
    }
    if (ii != nIn)
      throw ExcKal(std::string("AutoDiff: : Incorrect size in gradNum"));
  }

  const vector<T>& ax;
  const matrix<T>& Bx;
  const vector<T>& cx;
  vector<T>& yx;
  const vector<std::size_t>& axInd;
  const matrix<std::size_t>& BxInd;
  const vector<std::size_t>& cxInd;
  const vector<std::size_t>& yxInd;
};

template <class T>
class AD_Ax_Bx_Cx : public AutoDiff<T>
{ // Y(x) = A(x) - B(x) * C(x) 
public:
  AD_Ax_Bx_Cx(const matrix<T>* Axp, const matrix<T>* Bxp, const matrix<T>* Cxp, matrix<T>* Yxp,
               const matrix<std::size_t>* AxIndp, const matrix<std::size_t>* BxIndp, const matrix<std::size_t>* CxIndp, const matrix<std::size_t>* YxIndp) 
    : Ax(*Axp), Bx(*Bxp), Cx(*Cxp), Yx(*Yxp), AxInd(*AxIndp), BxInd(*BxIndp), CxInd(*CxIndp), YxInd(*YxIndp)   
  {
    // Check dimensions
    const std::size_t m = Ax.size1();
    const std::size_t n = Ax.size2();
    const std::size_t K = Bx.size2();
    const std::string strFunc("AD_Ax_Bx_Cx");
    checkDim(strFunc, "Ax", Ax, m, n);
    checkDim(strFunc, "Bx", Bx, m, K);
    checkDim(strFunc, "Cx", Cx, K, n);
    checkDim(strFunc, "Yx", Yx, m, m);
    checkDim(strFunc, "AxInd", AxInd, m, n);
    checkDim(strFunc, "BxInd", BxInd, m, K);
    checkDim(strFunc, "CxInd", CxInd, K, n);
    checkDim(strFunc, "YxInd", YxInd, m, n);
  }
  virtual ~AD_Ax_Bx_Cx() {}

  virtual void count(std::size_t& rMax, std::size_t& dMax) const 
  {
    const std::size_t m = Ax.size1();
    const std::size_t n = Ax.size2();
    const std::size_t K = Bx.size2();

    rMax += m*n; 
    dMax += m*n + 2*m*n*K;
  }
  virtual void f() 
  {
    const std::size_t m = Ax.size1();
    const std::size_t n = Ax.size2();
    const std::size_t K = Bx.size2();

    for (std::size_t i=0 ; i<m ; i++)
      for (std::size_t j=0 ; j<n ; j++)
      {
        Yx(i,j) = Ax(i,j);
        T dot = 0;
        for (std::size_t k=0 ; k<K ; k++)
          dot += Bx(i,k)*Cx(k,j);
        Yx(i,j) -= dot;
      }
  }
  virtual void fAD(vector<T>& d, vector<std::size_t>& col, vector<std::size_t>& row, vector<std::size_t>& rInd, std::size_t& nr, std::size_t& nd) 
  {
    const std::size_t m = Ax.size1();
    const std::size_t n = Ax.size2();
    const std::size_t K = Bx.size2();

    for (std::size_t i=0 ; i<m ; i++)
      for (std::size_t j=0 ; j<n ; j++)
      {
        rInd(nr) = YxInd(i,j);
        row(nr++) = nd;
        Yx(i,j) = Ax(i,j);
        d(nd) = 1;                              // d Yx(i,j) / d Ax(i,j) 
        col(nd++) = AxInd(i,j);
        T dot = 0;
        for (std::size_t k=0 ; k<K ; k++)
        {
          dot += Bx(i,k)*Cx(k,j);
          d(nd) = -Cx(k,j);                     // d Yx(i,j) / d Bx(i,k) 
          col(nd++) = BxInd(i,k);
          d(nd) = -Bx(i,k);                      // d Yx(i,j) / d Cx(k,j) 
          col(nd++) = CxInd(k,j);
        }
        Yx(i,j) -= dot;
      }
  }
protected:
  virtual void countSize(std::size_t& nIn, std::size_t& nOut)
  {
    nIn = Ax.size1() * Ax.size2() + Bx.size1() * Bx.size2() + Cx.size1() * Cx.size2();
    nOut = Yx.size1() * Yx.size2();
  }

  virtual AutoDiff<T>* createLocalObj(std::size_t& nc)
  { // Deallocated in freeLocalObj
    matrix<std::size_t>* AxInd_= new matrix<std::size_t>(Ax.size1(), Ax.size2()); 
    matrix<std::size_t>* BxInd_= new matrix<std::size_t>(Bx.size1(), Bx.size2()); 
    matrix<std::size_t>* CxInd_= new matrix<std::size_t>(Cx.size1(), Cx.size2()); 
    matrix<std::size_t>* YxInd_= new matrix<std::size_t>(Yx.size1(), Yx.size2()); 
    nc = 0;
    for (std::size_t i=0 ; i<Ax.size1() ; i++)
      for (std::size_t j=0 ; j<Ax.size2() ; j++)
        (*AxInd_)(i,j) = nc++;
    for (std::size_t i=0 ; i<Bx.size1() ; i++)
      for (std::size_t j=0 ; j<Bx.size2() ; j++)
        (*BxInd_)(i,j) = nc++;
    for (std::size_t i=0 ; i<Cx.size1() ; i++)
      for (std::size_t j=0 ; j<Cx.size2() ; j++)
        (*CxInd_)(i,j) = nc++;
    for (std::size_t i=0 ; i<Yx.size1() ; i++)
      for (std::size_t j=0 ; j<Yx.size2() ; j++)
        (*YxInd_)(i,j) = nc++;
    std::size_t nIn;
    std::size_t nOut;
    countSize(nIn, nOut);
    if (nc != nIn+nOut)
      throw ExcKal(std::string("AutoDiff: : Incorrect size in createLocalObj"));

    return new AD_Ax_Bx_Cx(&Ax, &Bx, &Cx, &Yx, AxInd_, BxInd_, CxInd_, YxInd_);
  }

  virtual void freeLocalObj()
  { // Free objects created by createLocalObj
    delete &AxInd;
    delete &BxInd;
    delete &CxInd;
    delete &YxInd;
  }

  virtual void testADgrad(const vector<T>& dOdx, vector<T>& grad) // Compute gradient with respect to input arguments
  { 
    std::size_t k=0;
    for (std::size_t i=0 ; i<Ax.size1() ; i++)
      for (std::size_t j=0 ; j<Ax.size2() ; j++)
        grad(k++) = dOdx(AxInd(i,j));
    for (std::size_t i=0 ; i<Bx.size1() ; i++)
      for (std::size_t j=0 ; j<Bx.size2() ; j++)
        grad(k++) = dOdx(BxInd(i,j));
    for (std::size_t i=0 ; i<Cx.size1() ; i++)
      for (std::size_t j=0 ; j<Cx.size2() ; j++)
        grad(k++) = dOdx(CxInd(i,j));
    if (k != grad.size())
      throw ExcKal(std::string("AutoDiff: : Incorrect size in grad"));
  }

  virtual void gradNum(const T epsilon, matrix<T>& Jacobian)
  {
    std::size_t nIn;
    std::size_t nOut;
    countSize(nIn, nOut);
    Jacobian.resize(nIn, nOut);
    f();
    matrix<T> Ax_(Ax);
    matrix<T> Bx_(Bx);
    matrix<T> Cx_(Cx);
    matrix<T> Yx_(Yx);

    std::size_t ii=0;
    for (std::size_t i=0 ; i<Ax.size1() ; i++)
      for (std::size_t j=0 ; j<Ax.size2() ; j++)
      {        
        Ax_(i,j) += epsilon;
        AD_Ax_Bx_Cx<T> ad(&Ax_, &Bx_, &Cx_, &Yx_, &AxInd, &BxInd, &CxInd, &YxInd);
        ad.f();
        std::size_t jj=0;
        for (std::size_t iii=0 ; iii<Yx.size1() ; iii++)
          for (std::size_t jjj=0 ; jjj<Yx.size2() ; jjj++)
            Jacobian(ii, jj++) = (Yx_(iii,jjj) - Yx(iii,jjj))/epsilon;
        ii++;
        Ax_(i,j) = Ax(i,j);
        if (jj != nOut)
          throw ExcKal(std::string("AutoDiff: : Incorrect size in gradNum"));
      }
    for (std::size_t i=0 ; i<Bx.size1() ; i++)
      for (std::size_t j=0 ; j<Bx.size2() ; j++)
      {        
        Bx_(i,j) += epsilon;
        AD_Ax_Bx_Cx<T> ad(&Ax_, &Bx_, &Cx_, &Yx_, &AxInd, &BxInd, &CxInd, &YxInd);
        ad.f();
        std::size_t jj=0;
        for (std::size_t iii=0 ; iii<Yx.size1() ; iii++)
          for (std::size_t jjj=0 ; jjj<Yx.size2() ; jjj++)
            Jacobian(ii, jj++) = (Yx_(iii,jjj) - Yx(iii,jjj))/epsilon;
        ii++;
        Bx_(i,j) = Bx(i,j);
        if (jj != nOut)
          throw ExcKal(std::string("AutoDiff: : Incorrect size in gradNum"));
      }
    for (std::size_t i=0 ; i<Cx.size1() ; i++)
      for (std::size_t j=0 ; j<Cx.size2() ; j++)
      {        
        Cx_(i,j) += epsilon;
        AD_Ax_Bx_Cx<T> ad(&Ax_, &Bx_, &Cx_, &Yx_, &AxInd, &BxInd, &CxInd, &YxInd);
        ad.f();
        std::size_t jj=0;
        for (std::size_t iii=0 ; iii<Yx.size1() ; iii++)
          for (std::size_t jjj=0 ; jjj<Yx.size2() ; jjj++)
            Jacobian(ii, jj++) = (Yx_(iii,jjj) - Yx(iii,jjj))/epsilon;
        ii++;
        Cx_(i,j) = Cx(i,j);
        if (jj != nOut)
          throw ExcKal(std::string("AutoDiff: : Incorrect size in gradNum"));
      }
    if (ii != nIn)
      throw ExcKal(std::string("AutoDiff: : Incorrect size in gradNum"));
  }

  const matrix<T>& Ax;
  const matrix<T>& Bx;
  const matrix<T>& Cx;
  matrix<T>& Yx;
  const matrix<std::size_t>& AxInd;
  const matrix<std::size_t>& BxInd;
  const matrix<std::size_t>& CxInd;
  const matrix<std::size_t>& YxInd;
};


template <class T>
class AD_InstrPricing : public AutoDiff<T>
{ // Y(x) = Instrument pricing
public:
  AD_InstrPricing(const vector<std::size_t>* atp, const vector<std::size_t>* oTp, const vector<std::size_t>* oIp, const vector<T>* tcp, const vector<std::size_t>* datesStepRelp, 
                  const matrix<T>* intEfp, const vector<T>* xpp, matrix<T>* op, vector<T>* eoxp, vector<T>* denp, vector<T>* gp, matrix<T>* Jacobianp, 
                  const matrix<std::size_t>* intEfIndp, const vector<std::size_t>* xpIndp, matrix<std::size_t>* oIndp, vector<std::size_t>* eoxIndp, vector<std::size_t>* denIndp, vector<std::size_t>* gIndp, matrix<std::size_t>* JacobianIndp) 
    : at(*atp), oT(*oTp), oI(*oIp), tc(*tcp), datesStepRel(*datesStepRelp),
    intEf(*intEfp), xp(*xpp), o(*op), eox(*eoxp), den(*denp), g(*gp), Jacobian(*Jacobianp), intEfInd(*intEfIndp), xpInd(*xpIndp), oInd(*oIndp), eoxInd(*eoxIndp), denInd(*denIndp), gInd(*gIndp), JacobianInd(*JacobianIndp)
  {
    // Check dimensions
    const std::size_t nT = intEf.size1()-1;
    const std::size_t np = intEf.size2();
    const std::size_t noT = oT.size();
    const std::size_t nz = at.size();
    const std::size_t nxc = datesStepRel.size();
    const std::size_t nxs = np + nxc;
    const std::size_t nx = xp.size();

    std::string strFunc("AD_InstrPricing");
    checkDim(strFunc, "at", at, nz);
    checkDim(strFunc, "oT", oT, noT);
    checkDim(strFunc, "oI", oI, nz+1);
    checkDim(strFunc, "tc", tc, noT);
    checkDim(strFunc, "datesStepRel", datesStepRel, nxc);

    checkDim(strFunc, "intEf", intEf, nT+1, np);
    checkDim(strFunc, "xp", xp, nx);
    checkDim(strFunc, "o", o, noT, nxs);
    checkDim(strFunc, "eox", eox, noT);
    checkDim(strFunc, "den", den, nz);
    checkDim(strFunc, "g", g, nz);
    checkDim(strFunc, "Jacobian", Jacobian, nz, nxs);

    checkDim(strFunc, "intEfInd", intEfInd, nT+1, np);
    checkDim(strFunc, "xpInd", xpInd, nx);
    checkDim(strFunc, "oInd", oInd, noT, nxs);
    checkDim(strFunc, "eoxInd", eoxInd, noT);
    checkDim(strFunc, "denInd", denInd, nz);
    checkDim(strFunc, "gInd", gInd, nz);
    checkDim(strFunc, "JacobianInd", JacobianInd, nz, nxs);
  }
  virtual ~AD_InstrPricing() {}

  virtual void count(std::size_t& rMax, std::size_t& dMax) const 
  {
    const std::size_t np = intEf.size2();
    const std::size_t noT = oT.size();
    const std::size_t nz = at.size();
    const std::size_t nxc = datesStepRel.size();
    const std::size_t nxs = np + nxc;

    rMax += noT*np; dMax += noT*np;  // o
    for (std::size_t i=0 ; i<nz ; i++)
    {
      if (at(i) == 20) // IBOR
      {
        rMax += 1; dMax += 2*nxs;  // eox
        rMax += 1; dMax += 1;  // g
        rMax += nxs; dMax += 2*nxs;  // Jacobian
      }
      else if (at(i) == 24) // OIS
      {
        for (std::size_t k=oI(i) ; k<oI(i+1) ; k++)
        {  rMax += 1; dMax += 2*nxs; } // eox
        rMax += 1; dMax += oI(i+1)-oI(i)-1;  // den
        rMax += 1; dMax += 3;  // g
        for (std::size_t j=0 ; j<nxs ; j++)
        {  rMax += 1; dMax += 5 + 2*(oI(i+1) - oI(i)-1); } // Jacobian
      }
      else
        throw ExcKal("Hardcoded instrument types did not work - check!");
    }
  }
  virtual void f() 
  {
    const std::size_t nT = intEf.size1()-1;
    const std::size_t np = intEf.size2();
    const std::size_t noT = oT.size();
    const std::size_t nz = at.size();
    const std::size_t nxc = datesStepRel.size();
    const std::size_t nxs = np + nxc;

    //pretty_print(std::cout, tc, "tc", 5, 12);

    /*T dt = 1.0/365;
    matrix<T> intEs(nT+1, nxc, 0.0);
    for (std::size_t j=0 ; j<nxc ; j++)
      for (std::size_t i=datesStepRel(j)+1 ; i<nT+1 ; i++)
        //intEs(i,j) = intEs(i-1,j) + dt;
        intEs(i,j) = (i-datesStepRel(j))*dt;*/

    /*for (std::size_t i=0 ; i<100 ; i++)
    {
      std::cout << i << " ";
      for (std::size_t j=0 ; j<nxc ; j++)
        std::cout << intEs(i,j) << " ";
      std::cout << std::endl;
    }*/

    for (std::size_t i=0 ; i<noT ; i++)
    {
      for (std::size_t j=0 ; j<np ; j++)
        o(i,j) = -intEf(oT(i), j);
      for (std::size_t j=0 ; j<nxc ; j++)
        // o(i,np+j) = -intEs(oT(i), j);
        if (oT(i) >= datesStepRel(j))
          o(i,np+j) = -((oT(i)-datesStepRel(j))/365.0);
        else
          o(i,np+j) = 0;
    }    

    //pretty_print(std::cout, o, "o", 5, 12);


    for (std::size_t i=0 ; i<nz ; i++)
    {
      if (at(i) == 20) // IBOR
      {
        std::size_t k=oI(i);

        T sum = 0;
        for (std::size_t j=0 ; j<nxs ; j++)
          sum += o(k,j)*xp(j);
        eox(k) = exp(-sum);

        g(i) = (eox(k)-1)/tc(k);

        for (std::size_t j=0 ; j<nxs ; j++)
          Jacobian(i,j) = -eox(k)*o(k,j)/tc(k);
      }
      else if (at(i) == 24) // OIS
      {
        for (std::size_t k=oI(i) ; k<oI(i+1) ; k++)
        {
          T sum = 0;
          for (std::size_t j=0 ; j<nxs ; j++)
            sum += o(k,j)*xp(j);
          eox(k) = exp(sum);
        }

        T sum = 0;
        for (std::size_t k=oI(i)+1 ; k<oI(i+1) ; k++)
          sum += tc(k)*eox(k);
        den(i) = sum;

        g(i) = (eox(oI(i)) - eox(oI(i+1)-1))/den(i);

        for (std::size_t j=0 ; j<nxs ; j++)
        {
          T sum2 = 0;
          for (std::size_t k=oI(i)+1 ; k<oI(i+1) ; k++)
            sum2 += tc(k)*eox(k)*o(k,j);
          //Jacobian(i,j) = (eox(oI(i))*o(oI(i),j)-eox(oI(i+1)-1)*o(oI(i+1)-1,j))/den(i) - (eox(oI(i))-eox(oI(i+1)-1))/(den(i)*den(i))*sum2;
          T J1 = (eox(oI(i))*o(oI(i),j)-eox(oI(i+1)-1)*o(oI(i+1)-1,j))/den(i);
          T J2 = (eox(oI(i))-eox(oI(i+1)-1))/(den(i)*den(i))*sum2;
          Jacobian(i,j) = J1 - J2;
        }
      }
      else
        throw ExcKal("Hardcoded instrument types did not work - check!");
    }
  }
  virtual void fAD(vector<T>& d, vector<std::size_t>& col, vector<std::size_t>& row, vector<std::size_t>& rInd, std::size_t& nr, std::size_t& nd) 
  {
    const std::size_t nT = intEf.size1()-1;
    const std::size_t np = intEf.size2();
    const std::size_t noT = oT.size();
    const std::size_t nz = at.size();
    const std::size_t nxc = datesStepRel.size();
    const std::size_t nxs = np + nxc;

    //pretty_print(std::cout, tc, "tc", 5, 12);

    /*T dt = 1.0/365;
    matrix<T> intEs(nT+1, nxc, 0.0);
    for (std::size_t j=0 ; j<nxc ; j++)
      for (std::size_t i=datesStepRel(j)+1 ; i<nT+1 ; i++)
        //intEs(i,j) = intEs(i-1,j) + dt;
        intEs(i,j) = (i-datesStepRel(j))*dt;*/

    /*for (std::size_t i=0 ; i<100 ; i++)
    {
      std::cout << i << " ";
      for (std::size_t j=0 ; j<nxc ; j++)
        std::cout << intEs(i,j) << " ";
      std::cout << std::endl;
    }*/

    for (std::size_t i=0 ; i<noT ; i++)
    {
      for (std::size_t j=0 ; j<np ; j++)
      {
        rInd(nr) = oInd(i,j);
        row(nr++) = nd;
        o(i,j) = -intEf(oT(i), j);
        d(nd) = -1;                             // d o(i,j) / d intEf(oT(i), j) 
        col(nd++) = intEfInd(oT(i), j);
      }
      for (std::size_t j=0 ; j<nxc ; j++)
        // o(i,np+j) = -intEs(oT(i), j);
        if (oT(i) >= datesStepRel(j))
          o(i,np+j) = -((oT(i)-datesStepRel(j))/365.0);
        else
          o(i,np+j) = 0;

    }    

    //pretty_print(std::cout, o, "o", 5, 12);


    for (std::size_t i=0 ; i<nz ; i++)
    {
      if (at(i) == 20) // IBOR
      {
        std::size_t k=oI(i);

        rInd(nr) = eoxInd(k);
        row(nr++) = nd;
        T sum = 0;
        for (std::size_t j=0 ; j<nxs ; j++)
          sum += o(k,j)*xp(j);
        eox(k) = exp(-sum);
        for (std::size_t j=0 ; j<nxs ; j++)
        {
          d(nd) = -eox(k)*xp(j);                // d eox(k) / d o(k,j) 
          col(nd++) = oInd(k,j);
          d(nd) = -eox(k)*o(k,j);               // d eox(k) / d xp(j) 
          col(nd++) = xpInd(j);        
        }

        rInd(nr) = gInd(i);
        row(nr++) = nd;
        g(i) = (eox(k)-1)/tc(k);
        d(nd) = 1/tc(k);                        // d g(i) / d eox(k) 
        col(nd++) = eoxInd(k);

        for (std::size_t j=0 ; j<nxs ; j++)
        {
          rInd(nr) = JacobianInd(i,j);
          row(nr++) = nd;
          Jacobian(i,j) = -eox(k)*o(k,j)/tc(k);
          d(nd) = -o(k,j)/tc(k);                // d Jacobian(i,j) / d eox(k) 
          col(nd++) = eoxInd(k);
          d(nd) = -eox(k)/tc(k);                // d Jacobian(i,j) / d o(k,j) 
          col(nd++) = oInd(k,j);
        }
      }
      else if (at(i) == 24) // OIS
      {
        for (std::size_t k=oI(i) ; k<oI(i+1) ; k++)
        {
          rInd(nr) = eoxInd(k);
          row(nr++) = nd;
          T sum = 0;
          for (std::size_t j=0 ; j<nxs ; j++)
            sum += o(k,j)*xp(j);
          eox(k) = exp(sum);
          for (std::size_t j=0 ; j<nxs ; j++)
          {
            d(nd) = eox(k)*xp(j);               // d eox(k) / d o(k,j) 
            col(nd++) = oInd(k,j);
            d(nd) = eox(k)*o(k,j);              // d eox(k) / d xp(j) 
            col(nd++) = xpInd(j);        
          }
        }

        rInd(nr) = denInd(i);
        row(nr++) = nd;
        T sum = 0;
        for (std::size_t k=oI(i)+1 ; k<oI(i+1) ; k++)
        {
          sum += tc(k)*eox(k);
          d(nd) = tc(k);                        // d den(i) / d eox(k) 
          col(nd++) = eoxInd(k);
        }
        den(i) = sum;

        rInd(nr) = gInd(i);
        row(nr++) = nd;
        g(i) = (eox(oI(i)) - eox(oI(i+1)-1))/den(i);
        d(nd) = 1/den(i);                       // d g(i) / d eox(oI(i)) 
        col(nd++) = eoxInd(oI(i));
        d(nd) =  - 1/den(i);                    // d g(i) / d eox(oI(i+1)-1) 
        col(nd++) = eoxInd(oI(i+1)-1);
        d(nd) = -g(i)/den(i);                   // d g(i) / d den(i) 
        col(nd++) = denInd(i);

        for (std::size_t j=0 ; j<nxs ; j++)
        {
          rInd(nr) = JacobianInd(i,j);
          row(nr++) = nd;
          T sum2 = 0;
          for (std::size_t k=oI(i)+1 ; k<oI(i+1) ; k++)
            sum2 += tc(k)*eox(k)*o(k,j);
          //Jacobian(i,j) = (eox(oI(i))*o(oI(i),j)-eox(oI(i+1)-1)*o(oI(i+1)-1,j))/den(i) - (eox(oI(i))-eox(oI(i+1)-1))/(den(i)*den(i))*sum2;
          T J1 = (eox(oI(i))*o(oI(i),j)-eox(oI(i+1)-1)*o(oI(i+1)-1,j))/den(i);
          T J2 = (eox(oI(i))-eox(oI(i+1)-1))/(den(i)*den(i))*sum2;
          Jacobian(i,j) = J1 - J2;
          d(nd) = o(oI(i),j)/den(i) - 1/(den(i)*den(i))*sum2; // d Jacobian(i,j) / d eox(oI(i)) 
          col(nd++) = eoxInd(oI(i));
          d(nd) = eox(oI(i))/den(i);            // d Jacobian(i,j) / d o(oI(i),j) 
          col(nd++) = oInd(oI(i),j);
          d(nd) = -o(oI(i+1)-1,j)/den(i) + 1/(den(i)*den(i))*sum2; // d Jacobian(i,j) / d eox(oI(i+1)-1) 
          col(nd++) = eoxInd(oI(i+1)-1);
          d(nd) = -eox(oI(i+1)-1)/den(i);       // d Jacobian(i,j) / d o(oI(i+1)-1,j) 
          col(nd++) = oInd(oI(i+1)-1,j);
          d(nd) = -J1/den(i) + 2* J2/den(i);    // d Jacobian(i,j) / d den(i) 
          col(nd++) = denInd(i);
          T J2a = (eox(oI(i))-eox(oI(i+1)-1))/(den(i)*den(i));
          for (std::size_t k=oI(i)+1 ; k<oI(i+1) ; k++)
          {
            d(nd) = -J2a*tc(k)*o(k,j);           // d Jacobian(i,j) / d eox(k) 
            col(nd++) = eoxInd(k);
            d(nd) = -J2a*tc(k)*eox(k);           // d Jacobian(i,j) / d o(k,j) 
            col(nd++) = oInd(k,j);          
          }

        }
      }
      else
        throw ExcKal("Hardcoded instrument types did not work - check!");
    }
  }
protected:
  virtual void countSize(std::size_t& nIn, std::size_t& nOut)
  {
    nIn = intEf.size1() * intEf.size2() + xp.size();
    nOut = o.size1() * o.size2() + eox.size() + den.size() + g.size() + Jacobian.size1() * Jacobian.size2();
  }

  virtual AutoDiff<T>* createLocalObj(std::size_t& nc)
  { // Deallocated in freeLocalObj
    matrix<std::size_t>* intEfInd_= new matrix<std::size_t>(intEf.size1(), intEf.size2()); 
    vector<std::size_t>* xpInd_= new vector<std::size_t>(xp.size()); 
    matrix<std::size_t>* oInd_= new matrix<std::size_t>(o.size1(), o.size2()); 
    vector<std::size_t>* eoxInd_= new vector<std::size_t>(eox.size()); 
    vector<std::size_t>* denInd_= new vector<std::size_t>(den.size()); 
    vector<std::size_t>* gInd_= new vector<std::size_t>(g.size()); 
    matrix<std::size_t>* JacobianInd_= new matrix<std::size_t>(Jacobian.size1(), Jacobian.size2()); 
    nc = 0;
    for (std::size_t i=0 ; i<intEf.size1() ; i++)
      for (std::size_t j=0 ; j<intEf.size2() ; j++)
        (*intEfInd_)(i,j) = nc++;
    for (std::size_t i=0 ; i<xp.size() ; i++)
      (*xpInd_)(i) = nc++;
    for (std::size_t i=0 ; i<o.size1() ; i++)
      for (std::size_t j=0 ; j<o.size2() ; j++)
        (*oInd_)(i,j) = nc++;
    for (std::size_t i=0 ; i<eox.size() ; i++)
      (*eoxInd_)(i) = nc++;
    for (std::size_t i=0 ; i<den.size() ; i++)
      (*denInd_)(i) = nc++;
    for (std::size_t i=0 ; i<g.size() ; i++)
      (*gInd_)(i) = nc++;
    for (std::size_t i=0 ; i<Jacobian.size1() ; i++)
      for (std::size_t j=0 ; j<Jacobian.size2() ; j++)
        (*JacobianInd_)(i,j) = nc++;
    std::size_t nIn;
    std::size_t nOut;
    countSize(nIn, nOut);
    if (nc != nIn+nOut)
      throw ExcKal(std::string("AutoDiff: : Incorrect size in createLocalObj"));

    return new AD_InstrPricing(&at, &oT, &oI, &tc, &datesStepRel,
                               &intEf, &xp, &o, &eox, &den, &g, &Jacobian, 
                               intEfInd_, xpInd_, oInd_, eoxInd_, denInd_, gInd_, JacobianInd_);
  }

  virtual void freeLocalObj()
  { // Free objects created by createLocalObj
    delete &intEfInd;
    delete &xpInd;
    delete &oInd;
    delete &eoxInd;
    delete &denInd;
    delete &gInd;
    delete &JacobianInd;
  }

  virtual void testADgrad(const vector<T>& dOdx, vector<T>& grad) // Compute gradient with respect to input arguments
  { 
    std::size_t k=0;
    for (std::size_t i=0 ; i<intEf.size1() ; i++)
      for (std::size_t j=0 ; j<intEf.size2() ; j++)
        grad(k++) = dOdx(intEfInd(i,j));
    for (std::size_t i=0 ; i<xp.size() ; i++)
      grad(k++) = dOdx(xpInd(i));
  }

  virtual void gradNum(const T epsilon, matrix<T>& JacobianAll)
  {
    std::size_t nIn;
    std::size_t nOut;
    countSize(nIn, nOut);
    JacobianAll.resize(nIn, nOut);
    f();

    matrix<T> intEf_(intEf);
    vector<T> xp_(xp);
    matrix<T> o_(o);
    vector<T> eox_(eox);
    vector<T> den_(den);
    vector<T> g_(g);
    matrix<T> Jacobian_(Jacobian);


    std::size_t ii=0;
    for (std::size_t i=0 ; i<intEf.size1() ; i++)
      for (std::size_t j=0 ; j<intEf.size2() ; j++)
      {        
        intEf_(i,j) += epsilon;
        AD_InstrPricing<T> ad(&at, &oT, &oI, &tc, &datesStepRel,
                              &intEf_, &xp_, &o_, &eox_, &den_, &g_, &Jacobian_, 
                              &intEfInd, &xpInd, &oInd, &eoxInd, &denInd, &gInd, &JacobianInd);
        ad.f();
        std::size_t jj=0;
        for (std::size_t iii=0 ; iii<o.size1() ; iii++)
          for (std::size_t jjj=0 ; jjj<o.size2() ; jjj++)
            JacobianAll(ii, jj++) = (o_(iii,jjj) - o(iii,jjj))/epsilon;

        for (std::size_t iii=0 ; iii<eox.size() ; iii++)
            JacobianAll(ii, jj++) = (eox_(iii) - eox(iii))/epsilon;

        for (std::size_t iii=0 ; iii<den.size() ; iii++)
            JacobianAll(ii, jj++) = (den_(iii) - den(iii))/epsilon;

        for (std::size_t iii=0 ; iii<g.size() ; iii++)
            JacobianAll(ii, jj++) = (g_(iii) - g(iii))/epsilon;

        for (std::size_t iii=0 ; iii<Jacobian.size1() ; iii++)
          for (std::size_t jjj=0 ; jjj<Jacobian.size2() ; jjj++)
            JacobianAll(ii, jj++) = (Jacobian_(iii,jjj) - Jacobian(iii,jjj))/epsilon;
        ii++;
        intEf_(i,j) = intEf(i,j);
        if (jj != nOut)
          throw ExcKal(std::string("AutoDiff: : Incorrect size in gradNum"));
      }
    for (std::size_t i=0 ; i<xp.size() ; i++)
    {        
      xp_(i) += epsilon;
      AD_InstrPricing<T> ad(&at, &oT, &oI, &tc, &datesStepRel,
                            &intEf_, &xp_, &o_, &eox_, &den_, &g_, &Jacobian_, 
                            &intEfInd, &xpInd, &oInd, &eoxInd, &denInd, &gInd, &JacobianInd);
      ad.f();
      std::size_t jj=0;
      for (std::size_t iii=0 ; iii<o.size1() ; iii++)
        for (std::size_t jjj=0 ; jjj<o.size2() ; jjj++)
          JacobianAll(ii, jj++) = (o_(iii,jjj) - o(iii,jjj))/epsilon;

      for (std::size_t iii=0 ; iii<eox.size() ; iii++)
          JacobianAll(ii, jj++) = (eox_(iii) - eox(iii))/epsilon;

      for (std::size_t iii=0 ; iii<den.size() ; iii++)
          JacobianAll(ii, jj++) = (den_(iii) - den(iii))/epsilon;

      for (std::size_t iii=0 ; iii<g.size() ; iii++)
          JacobianAll(ii, jj++) = (g_(iii) - g(iii))/epsilon;

      for (std::size_t iii=0 ; iii<Jacobian.size1() ; iii++)
        for (std::size_t jjj=0 ; jjj<Jacobian.size2() ; jjj++)
          JacobianAll(ii, jj++) = (Jacobian_(iii,jjj) - Jacobian(iii,jjj))/epsilon;
      ii++;
      xp_(i) = xp(i);
    if (jj != nOut)
      throw ExcKal(std::string("AutoDiff: : Incorrect size in gradNum"));
    }
    if (ii != nIn)
      throw ExcKal(std::string("AutoDiff: : Incorrect size in gradNum"));
  }

  const vector<std::size_t>& at;
  const vector<std::size_t>& oT;
  const vector<std::size_t>& oI;
  const vector<T>& tc;
  const vector<std::size_t>& datesStepRel;

  const matrix<T>& intEf;
  const vector<T>& xp;
  matrix<T>& o;
  vector<T>& eox;
  vector<T>& den;
  vector<T>& g;
  matrix<T>& Jacobian;
  const matrix<std::size_t>& intEfInd;
  const vector<std::size_t>& xpInd;
  matrix<std::size_t>& oInd;
  vector<std::size_t>& eoxInd;
  vector<std::size_t>& denInd;
  vector<std::size_t>& gInd;
  matrix<std::size_t>& JacobianInd;
};

template <class T>
class AD_LogLikelihood : public AutoDiff<T>
{ // sum_t (log |L(t)| + vepsilon(t)^T vSinvepsilon(t))
public:
  AD_LogLikelihood(const vector<matrix<T>*>* vLp, const vector<vector<T>*>* vepsilonp, const vector<vector<T>*>* vSinvepsilonp, T* logLikelihoodp, matrix<T>* fAllp,
  const vector<matrix<std::size_t>*>* vLIndp, const vector<vector<std::size_t>*>* vepsilonIndp, const vector<vector<std::size_t>*>* vSinvepsilonIndp, const std::size_t* logLikelihoodIndp) 
    : vL(*vLp), vepsilon(*vepsilonp), vSinvepsilon(*vSinvepsilonp), logLikelihood(*logLikelihoodp), fAll(*fAllp), vLInd(*vLIndp), vepsilonInd(*vepsilonIndp), vSinvepsilonInd(*vSinvepsilonIndp), logLikelihoodInd(*logLikelihoodIndp)
  {
    // Check dimensions
    const std::size_t nt = vL.size();

    const std::string strFunc("AD_LogLikelihood");
    checkDim(strFunc, "vL", vL, nt);
    checkDim(strFunc, "vepsilon", vepsilon, nt);
    checkDim(strFunc, "vSinvepsilon", vSinvepsilon, nt);
    checkDim(strFunc, "fAll", fAll, nt, 1);

    for (std::size_t t=0 ; t<nt ; t++)
    {
      const matrix<T>& L = *vL(t);
      const matrix<std::size_t>& LInd = *vLInd(t);
      const vector<T>& epsilon = *vepsilon(t);
      const vector<std::size_t>& epsilonInd = *vepsilonInd(t);
      const vector<T>& Sinvepsilon = *vSinvepsilon(t);
      const vector<std::size_t>& SinvepsilonInd = *vSinvepsilonInd(t);

      const std::size_t m = L.size1();

      checkDim(strFunc, "L", L, m, m);
      checkDim(strFunc, "epsilon", epsilon, m);
      checkDim(strFunc, "Sinvepsilon", Sinvepsilon, m);
      checkDim(strFunc, "LInd", LInd, m, m);
      checkDim(strFunc, "epsilonInd", epsilonInd, m);
      checkDim(strFunc, "SinvepsilonInd", SinvepsilonInd, m);    
    }
  }
  virtual ~AD_LogLikelihood() {}

  virtual void count(std::size_t& rMax, std::size_t& dMax) const 
  {
    rMax += 1;
    const std::size_t nt = vL.size();

    for (std::size_t t=0 ; t<nt ; t++)
    {
      const matrix<T>& L = *vL(t);
      const std::size_t m = L.size1();
      dMax += 3*m;
    }
  }
  virtual void f() 
  {
    const std::size_t nt = vL.size();

    T f = 0;

    for (std::size_t t=0 ; t<nt ; t++)
    {
      const matrix<T>& L = *vL(t);
      const vector<T>& epsilon = *vepsilon(t);
      const vector<T>& Sinvepsilon = *vSinvepsilon(t);

      const std::size_t m = L.size1();

      T ft = -0.5*m*log(2*M_PI);

      T dot = 0;
      for (std::size_t i=0 ; i<m ; i++)
        dot += log(L(i,i));
      ft -= dot;

      dot = 0;
      for (std::size_t i=0 ; i<m ; i++)
        dot += 0.5*epsilon(i)*Sinvepsilon(i);
      ft -= dot;

      // std::cout << "f = " << f << std::endl;
      fAll(t,0) = ft;
      f += ft;
    }
    logLikelihood = f;
  }
  virtual void fAD(vector<T>& d, vector<std::size_t>& col, vector<std::size_t>& row, vector<std::size_t>& rInd, std::size_t& nr, std::size_t& nd) 
  {
    const std::size_t nt = vL.size();

    T f = 0;

    rInd(nr) = logLikelihoodInd;
    row(nr++) = nd;
    for (std::size_t t=0 ; t<nt ; t++)
    {
      const matrix<T>& L = *vL(t);
      const matrix<std::size_t>& LInd = *vLInd(t);
      const vector<T>& epsilon = *vepsilon(t);
      const vector<std::size_t>& epsilonInd = *vepsilonInd(t);
      const vector<T>& Sinvepsilon = *vSinvepsilon(t);
      const vector<std::size_t>& SinvepsilonInd = *vSinvepsilonInd(t);

      const std::size_t m = L.size1();

      T ft = -0.5*m*log(2*M_PI);

      T dot = 0;
      for (std::size_t i=0 ; i<m ; i++)
      {
        dot += log(L(i,i));
        d(nd) = -1/L(i,i);                        // d f / d L(i,i) 
        col(nd++) = LInd(i,i);
      }
      ft -= dot;

      dot = 0;
      for (std::size_t i=0 ; i<m ; i++)
      {
        dot += 0.5*epsilon(i)*Sinvepsilon(i);
        d(nd) = -0.5*Sinvepsilon(i);              // d f / d epsilon(i) 
        col(nd++) = epsilonInd(i);
        d(nd) = -0.5*epsilon(i);                  // d f / d Sinvepsilon(i) 
        col(nd++) = SinvepsilonInd(i);
      }
      ft -= dot;
      // std::cout << "f = " << f << std::endl;
      fAll(t,0) = ft;
      f += ft;
    }
    logLikelihood = f;
  }
protected:
  virtual void countSize(std::size_t& nIn, std::size_t& nOut)
  {
    const std::size_t nt = vL.size();
    nIn = 0;
    for (std::size_t t=0 ; t<nt ; t++)
    {
      const matrix<T>& L = *vL(t);
      const vector<T>& epsilon = *vepsilon(t);
      const vector<T>& Sinvepsilon = *vSinvepsilon(t);
      nIn += L.size1()*L.size2() + epsilon.size() + Sinvepsilon.size();
    }
    nOut = 1;
  }

  virtual AutoDiff<T>* createLocalObj(std::size_t& nc)
  { // Deallocated in freeLocalObj
    const std::size_t nt = vL.size();

    nc = 0;
    vector<matrix<std::size_t>*>* vLInd_ = new vector<matrix<std::size_t>*>(nt);
    vector<vector<std::size_t>*>* vepsilonInd_ = new vector<vector<std::size_t>*>(nt);
    vector<vector<std::size_t>*>* vSinvepsilonInd_ = new vector<vector<std::size_t>*>(nt);
    for (std::size_t t=0 ; t<nt ; t++)
    {
      const matrix<T>& L = *vL(t);
      const vector<T>& epsilon = *vepsilon(t);
      const vector<T>& Sinvepsilon = *vSinvepsilon(t);

      (*vLInd_)(t) = new matrix<std::size_t>(L.size1(), L.size2()); 
      (*vepsilonInd_)(t) = new vector<std::size_t>(epsilon.size()); 
      (*vSinvepsilonInd_)(t) = new vector<std::size_t>(Sinvepsilon.size()); 
      matrix<std::size_t>& LInd = *(*vLInd_)(t);
      vector<std::size_t>& epsilonInd = *(*vepsilonInd_)(t);
      vector<std::size_t>& SinvepsilonInd = *(*vSinvepsilonInd_)(t);
      for (std::size_t i=0 ; i<L.size1() ; i++)
        for (std::size_t j=0 ; j<L.size2() ; j++)
          LInd(i,j) = nc++;
      for (std::size_t i=0 ; i<epsilon.size() ; i++)
        epsilonInd(i) = nc++;
      for (std::size_t i=0 ; i<Sinvepsilon.size() ; i++)
        SinvepsilonInd(i) = nc++;
    }
    std::size_t* logLikelihoodInd_ = new std::size_t;
    *logLikelihoodInd_ = nc++;
    std::size_t nIn;
    std::size_t nOut;
    countSize(nIn, nOut);
    if (nc != nIn+nOut)
      throw ExcKal(std::string("AutoDiff: : Incorrect size in createLocalObj"));

    return new AD_LogLikelihood<T>(&vL, &vepsilon, &vSinvepsilon, &logLikelihood, &fAll, 
                                   vLInd_, vepsilonInd_, vSinvepsilonInd_, logLikelihoodInd_);

  }

  virtual void freeLocalObj()
  { // Free objects created by createLocalObj
    const std::size_t nt = vL.size();
    vector<matrix<std::size_t>*>* vLInd_ = new vector<matrix<std::size_t>*>(nt);
    vector<vector<std::size_t>*>* vepsilonInd_ = new vector<vector<std::size_t>*>(nt);
    vector<vector<std::size_t>*>* vSinvepsilonInd_ = new vector<vector<std::size_t>*>(nt);
    for (std::size_t t=0 ; t<nt ; t++)
    {
      delete vLInd(t); 
      delete vepsilonInd(t); 
      delete vSinvepsilonInd(t); 
    }
    delete &vLInd; 
    delete &vepsilonInd; 
    delete &vSinvepsilonInd;
    delete &logLikelihoodInd;
  }

  virtual void testADgrad(const vector<T>& dOdx, vector<T>& grad) // Compute gradient with respect to input arguments
  { 
    std::size_t k=0;

    const std::size_t nt = vL.size();

    vector<matrix<std::size_t>*>* vLInd_ = new vector<matrix<std::size_t>*>(nt);
    vector<vector<std::size_t>*>* vepsilonInd_ = new vector<vector<std::size_t>*>(nt);
    vector<vector<std::size_t>*>* vSinvepsilonInd_ = new vector<vector<std::size_t>*>(nt);
    for (std::size_t t=0 ; t<nt ; t++)
    {
      const matrix<T>& L = *vL(t);
      const vector<T>& epsilon = *vepsilon(t);
      const vector<T>& Sinvepsilon = *vSinvepsilon(t);

      matrix<std::size_t>& LInd = *(vLInd(t));
      vector<std::size_t>& epsilonInd = *(vepsilonInd(t));
      vector<std::size_t>& SinvepsilonInd = *(vSinvepsilonInd(t));

      for (std::size_t i=0 ; i<L.size1() ; i++)
        for (std::size_t j=0 ; j<L.size2() ; j++)
          grad(k++) = dOdx(LInd(i,j));
      for (std::size_t i=0 ; i<epsilon.size() ; i++)
        grad(k++) = dOdx(epsilonInd(i));
      for (std::size_t i=0 ; i<Sinvepsilon.size() ; i++)
        grad(k++) = dOdx(SinvepsilonInd(i));
    }

  }

  virtual void gradNum(const T epsilonNum, matrix<T>& Jacobian)
  {
    std::size_t nIn;
    std::size_t nOut;
    countSize(nIn, nOut);
    Jacobian.resize(nIn, nOut);
    f();
    const std::size_t nt = vL.size();

    vector<matrix<T>*> vL_(nt);
    vector<vector<T>*> vepsilon_(nt);
    vector<vector<T>*> vSinvepsilon_(nt);
    for (std::size_t t=0 ; t<nt ; t++)
    {
      const matrix<T>& L = *vL(t);
      const vector<T>& epsilon = *vepsilon(t);
      const vector<T>& Sinvepsilon = *vSinvepsilon(t);

      vL_(t) = new matrix<T>(L); 
      vepsilon_(t) = new vector<T>(epsilon); 
      vSinvepsilon_(t) = new vector<T>(Sinvepsilon); 
    }

    std::size_t ii=0;
    for (std::size_t t=0 ; t<nt ; t++)
    {
      const matrix<T>& L = *vL(t);
      const vector<T>& epsilon = *vepsilon(t);
      const vector<T>& Sinvepsilon = *vSinvepsilon(t);

      matrix<T>& L_ = *vL_(t);
      vector<T>& epsilon_ = *vepsilon_(t);
      vector<T>& Sinvepsilon_ = *vSinvepsilon_(t);
      T logLikelihood_;

      for (std::size_t i=0 ; i<L.size1() ; i++)
        for (std::size_t j=0 ; j<L.size2() ; j++)
        {
          L_(i,j) += epsilonNum;
          AD_LogLikelihood<T> ad(&vL_, &vepsilon_, &vSinvepsilon_, &logLikelihood_, &fAll, 
                                 &vLInd, &vepsilonInd, &vSinvepsilonInd, &logLikelihoodInd);
          ad.f();
          std::size_t jj=0;
          Jacobian(ii, jj++) = (logLikelihood_ - logLikelihood)/epsilonNum;
          ii++;
          L_(i,j) = L(i,j);
          if (jj != nOut)
            throw ExcKal(std::string("AutoDiff: : Incorrect size in gradNum"));        
        }
      for (std::size_t i=0 ; i<epsilon.size() ; i++)
      {        
        epsilon_(i) += epsilonNum;
        AD_LogLikelihood<T> ad(&vL_, &vepsilon_, &vSinvepsilon_, &logLikelihood_, &fAll, 
                               &vLInd, &vepsilonInd, &vSinvepsilonInd, &logLikelihoodInd);
        ad.f();
        std::size_t jj=0;
          Jacobian(ii, jj++) = (logLikelihood_ - logLikelihood)/epsilonNum;
        ii++;
        epsilon_(i) = epsilon(i);
        if (jj != nOut)
          throw ExcKal(std::string("AutoDiff: : Incorrect size in gradNum"));
      }
      for (std::size_t i=0 ; i<Sinvepsilon.size() ; i++)
      {        
        Sinvepsilon_(i) += epsilonNum;
        AD_LogLikelihood<T> ad(&vL_, &vepsilon_, &vSinvepsilon_, &logLikelihood_, &fAll, 
                               &vLInd, &vepsilonInd, &vSinvepsilonInd, &logLikelihoodInd);
        ad.f();
        std::size_t jj=0;
          Jacobian(ii, jj++) = (logLikelihood_ - logLikelihood)/epsilonNum;
        ii++;
        Sinvepsilon_(i) = Sinvepsilon(i);
        if (jj != nOut)
          throw ExcKal(std::string("AutoDiff: : Incorrect size in gradNum"));
      }
    }
    if (ii != nIn)
      throw ExcKal(std::string("AutoDiff: : Incorrect size in gradNum"));

    for (std::size_t t=0 ; t<nt ; t++)
    {
      delete vL_(t);
      delete vepsilon_(t); 
      delete vSinvepsilon_(t); 
    }
  }

  const vector<matrix<T>*>& vL;
  const vector<vector<T>*>& vepsilon;
  const vector<vector<T>*>& vSinvepsilon;
  T& logLikelihood;
  matrix<T>& fAll;
  const vector<matrix<std::size_t>*>& vLInd;
  const vector<vector<std::size_t>*>& vepsilonInd;
  const vector<vector<std::size_t>*>& vSinvepsilonInd;
  const std::size_t& logLikelihoodInd;
};


template<class T>
void gradKalmanLogLikelihood(const vector<vector<T>*>& vz, const vector<SparseMatrix<T>*>& vA, const vector<SparseMatrix<T>*>& vB, const vector<SparseMatrix<T>*>& vD, const vector<SparseMatrix<T>*>& vG, const vector<SparseMatrix<T>*>& vIz,
                             const vector<T>& thetaF, const vector<T>& Sigmav, const vector<T>& Sigmaw, const vector<T>& ax, const vector<T>& Sigma0,
                             const vector<vector<std::size_t>*>& vat, const vector<vector<std::size_t>*>& voT, const vector<vector<std::size_t>*>& voI, const vector<vector<T>*>& vtc, const matrix<T>& Ef, const vector<vector<std::size_t>*>& vdatesStepRel,
											       T& f, T& obj, matrix<T>& gradObjEf, vector<T>& gradObjthetaF, vector<T>& gradObjSigmav, vector<T>& gradObjSigmaw, vector<T>& gradObjax, vector<T>& gradObjSigma0, matrix<T>& retM)
{
  std::size_t nt = vz.size(); // Number of historical dates
  std::size_t np;             // Number of principal components (held constant)
  std::size_t nthetaF;        // Number of thetaF-variables
  std::size_t nSigmav;        // Number of Sigmav-variables
  std::size_t nSigmaw;        // Number of Sigmaw-variables
  std::size_t nax;            // Number of ax-variables
  std::size_t nSigma0;        // Number of Sigma0-variables
  vector<std::size_t> nx(nt); // Number of x-variables
  vector<std::size_t> nxs(nt);// Number of systematic x-variables (np + nxc(t))
  vector<std::size_t> nxu(nt);// Number of unsystematic x-variables (number of instruments)
  vector<std::size_t> nxc(nt);// Number of xc-variables (central bank steps/cuts)
  vector<std::size_t> nz(nt); // Number of z-variables
  vector<std::size_t> noT(nt);// Number of components to compute instrument prices
  std::size_t nT;             // Number of time elemements in eigenvectors
  
  nT = Ef.size1();
  np = Ef.size2();
  nthetaF = thetaF.size();
  nSigmav = Sigmav.size();
  nSigmaw = Sigmaw.size();
  nax = ax.size();
  nSigma0 = Sigma0.size();
  for (std::size_t t=0 ; t<nt ; t++)
  {
    nx(t) = vA(t)->size1();
    nxc(t) = vdatesStepRel(t)->size();
    nz(t) = vz(t)->size();
    nxu(t) = vIz(t)->size2();
    noT(t) = voT(t)->size();
    nxs(t) = np + nxc(t);
    if (vA(t)->size2() != nthetaF)
      throw ExcKal("Incorrect number of columns in A");
    if (vB(t)->size1() != nthetaF)
      throw ExcKal("Incorrect number of rows in B");
    if (t==0)
    {
      if (vB(t)->size2() != nSigma0)
        throw ExcKal("Incorrect number of columns in B");
    }
    else
      if (vB(t)->size2() != nx(t))
        throw ExcKal("Incorrect number of columns in B");
    if (vD(t)->size1() != nthetaF)
      throw ExcKal("Incorrect number of rows in D");
    if (vD(t)->size2() != nSigmaw)
      throw ExcKal("Incorrect number of columns in D");
    if (vG(t)->size1() != nz(t))
      throw ExcKal("Incorrect number of rows in G");
    if (vG(t)->size2() != nSigmav)
      throw ExcKal("Incorrect number of columns in G");
    if (vIz(t)->size1() != nz(t))
      throw ExcKal("Incorrect number of rows in Iz");
    if (vIz(t)->size2() != nxu(t))
      throw ExcKal("Incorrect number of columns in Iz");
    if (vat(t)->size() != nz(t))
      throw ExcKal("Incorrect size of at");
    if (voI(t)->size() != nz(t)+1)
      throw ExcKal("Incorrect size of oI");
    if (voT(t)->size() != noT(t))
      throw ExcKal("Incorrect size of oT");
    if (vtc(t)->size() != noT(t))
      throw ExcKal("Incorrect size of tc");
  }

  matrix<T> intEf(nT+1, np);

  // Predictor step
  vector<vector<T>*> vxp(nt);
  vector<matrix<T>*> vBPu(nt);
  vector<matrix<T>*> vBPuBT(nt);
  vector<matrix<T>*> vPp1(nt);
  vector<matrix<T>*> vDSigmaw(nt);
  vector<matrix<T>*> vPp(nt);

  // Instrument pricing
  vector<matrix<T>*> vo(nt);
  vector<vector<T>*> veox(nt);
  vector<vector<T>*> vden(nt);
  vector<vector<T>*> vg(nt); // Instrument price 
  vector<matrix<T>*> vJacobian(nt); // Jacobian of instrument prices

  // Update step
  vector<vector<T>*> vepsilon(nt);
  vector<matrix<T>*> vHP(nt);
  vector<matrix<T>*> vS1(nt);
  vector<matrix<T>*> vGSigmav(nt);
  vector<matrix<T>*> vS2(nt);
  vector<matrix<T>*> vS(nt);
  vector<matrix<T>*> vL(nt);
  vector<matrix<T>*> vK1(nt);
  vector<matrix<T>*> vK(nt);
  vector<vector<T>*> vxu(nt);
  vector<matrix<T>*> vPu(nt);
  vector<vector<T>*> vSinvepsilon1(nt);
  vector<vector<T>*> vSinvepsilon(nt);

  // Index variables for automatic differentiation

  // Initial variables
  matrix<std::size_t> EfInd(nT, np);
  vector<std::size_t> thetaFInd(nthetaF);
  vector<std::size_t> SigmavInd(nSigmav);
  vector<std::size_t> SigmawInd(nSigmaw);
  vector<std::size_t> axInd(nax);
  vector<std::size_t> Sigma0Ind(nSigma0);

  matrix<std::size_t> intEfInd(nT+1, np);

  // Predictor step
  vector<vector<std::size_t>*> vxpInd(nt);
  vector<matrix<std::size_t>*> vBPuInd(nt);
  vector<matrix<std::size_t>*> vBPuBTInd(nt);
  vector<matrix<std::size_t>*> vPp1Ind(nt);
  vector<matrix<std::size_t>*> vDSigmawInd(nt);
  vector<matrix<std::size_t>*> vPpInd(nt);

  // Instrument pricing
  vector<matrix<std::size_t>*> voInd(nt);
  vector<vector<std::size_t>*> veoxInd(nt);
  vector<vector<std::size_t>*> vdenInd(nt);
  vector<vector<std::size_t>*> vgInd(nt); // Instrument price 
  vector<matrix<std::size_t>*> vJacobianInd(nt); // Jacobian of instrument prices

  // Update step
  vector<vector<std::size_t>*> vepsilonInd(nt);
  vector<matrix<std::size_t>*> vHPInd(nt);
  vector<matrix<std::size_t>*> vS1Ind(nt);
  vector<matrix<std::size_t>*> vGSigmavInd(nt);
  vector<matrix<std::size_t>*> vS2Ind(nt);
  vector<matrix<std::size_t>*> vSInd(nt);
  vector<matrix<std::size_t>*> vLInd(nt);
  vector<matrix<std::size_t>*> vK1Ind(nt);
  vector<matrix<std::size_t>*> vKInd(nt);
  vector<vector<std::size_t>*> vxuInd(nt);
  vector<matrix<std::size_t>*> vPuInd(nt);
  vector<vector<std::size_t>*> vSinvepsilon1Ind(nt);
  vector<vector<std::size_t>*> vSinvepsilonInd(nt);

  // Log-likelihood value

  std::size_t fInd;

  // Automatic differentiation variables

  std::size_t nc=0;

  // Create indexing for automatic differentiation
  // Initial variables

  for (std::size_t i=0 ; i<nT ; i++)
    for (std::size_t j=0 ; j<np ; j++)
      EfInd(i,j) = nc++;
  for (std::size_t i=0 ; i<nthetaF ; i++)
    thetaFInd(i) = nc++;
  for (std::size_t i=0 ; i<nSigmav ; i++)
    SigmavInd(i) = nc++;
  for (std::size_t i=0 ; i<nSigmaw ; i++)
    SigmawInd(i) = nc++;
  for (std::size_t i=0 ; i<nax ; i++)
    axInd(i) = nc++;
  for (std::size_t i=0 ; i<nSigma0 ; i++)
    Sigma0Ind(i) = nc++;

  for (std::size_t i=0 ; i<nT+1 ; i++)
    for (std::size_t j=0 ; j<np ; j++)
      intEfInd(i,j) = nc++;

  for (std::size_t t=0 ; t<nt ; t++)
  {
    // Predictor step
    vxp(t) = new vector<T>(nx(t));
    if (t==0)
      vBPu(t) = new matrix<T>(nthetaF, nSigma0);
    else
      vBPu(t) = new matrix<T>(nthetaF, nx(t));
    vBPuBT(t) = new matrix<T>(nthetaF, nthetaF);
    vPp1(t) = new matrix<T>(nx(t), nx(t));
    vDSigmaw(t) = new matrix<T>(nthetaF, nSigmaw);
    vPp(t) = new matrix<T>(nx(t), nx(t));

    // Instrument pricing
    vo(t) = new matrix<T>(noT(t), nxs(t));
    veox(t) = new vector<T>(noT(t));
    vden(t) = new vector<T>(nz(t));
    vg(t) = new vector<T>(nz(t));
    vJacobian(t) = new matrix<T>(nz(t), nxs(t));

    // Update step
    vepsilon(t) = new vector<T>(nz(t));
    vHP(t) = new matrix<T>(nz(t), nx(t));
    vS1(t) = new matrix<T>(nz(t), nz(t));
    vGSigmav(t) = new matrix<T>(nz(t), nSigmav);
    vS2(t) = new matrix<T>(nz(t), nz(t));
    vS(t) = new matrix<T>(nz(t), nz(t));
    vL(t) = new matrix<T>(nz(t), nz(t));
    vK1(t) = new matrix<T>(nx(t), nz(t));
    vK(t) = new matrix<T>(nx(t), nz(t));
    vxu(t) = new vector<T>(nx(t));
    vPu(t) = new matrix<T>(nx(t), nx(t));
    vSinvepsilon1(t) = new vector<T>(nz(t));
    vSinvepsilon(t) = new vector<T>(nz(t));

    // Memory allocation of index variables for automatic differentiation

    // Predictor step
    vxpInd(t) = new vector<std::size_t>(nx(t));
    if (t==0)
      vBPuInd(t) = new matrix<std::size_t>(nthetaF, nSigma0);
    else
      vBPuInd(t) = new matrix<std::size_t>(nthetaF, nx(t));
    vBPuBTInd(t) = new matrix<std::size_t>(nthetaF, nthetaF);
    vPp1Ind(t) = new matrix<std::size_t>(nx(t), nx(t));
    vDSigmawInd(t) = new matrix<std::size_t>(nthetaF, nSigmaw);
    vPpInd(t) = new matrix<std::size_t>(nx(t), nx(t));

    // Instrument pricing
    voInd(t) = new matrix<std::size_t>(noT(t), nxs(t));
    veoxInd(t) = new vector<std::size_t>(noT(t));
    vdenInd(t) = new vector<std::size_t>(nz(t));
    vgInd(t) = new vector<std::size_t>(nz(t));
    vJacobianInd(t) = new matrix<std::size_t>(nz(t), nxs(t));

    // Update step
    vepsilonInd(t) = new vector<std::size_t>(nz(t));
    vHPInd(t) = new matrix<std::size_t>(nz(t), nx(t));
    vS1Ind(t) = new matrix<std::size_t>(nz(t), nz(t));
    vGSigmavInd(t) = new matrix<std::size_t>(nz(t), nSigmav);
    vS2Ind(t) = new matrix<std::size_t>(nz(t), nz(t));
    vSInd(t) = new matrix<std::size_t>(nz(t), nz(t));
    vLInd(t) = new matrix<std::size_t>(nz(t), nz(t));
    vK1Ind(t) = new matrix<std::size_t>(nx(t), nz(t));
    vKInd(t) = new matrix<std::size_t>(nx(t), nz(t));
    vxuInd(t) = new vector<std::size_t>(nx(t));
    vPuInd(t) = new matrix<std::size_t>(nx(t), nx(t));
    vSinvepsilon1Ind(t) = new vector<std::size_t>(nz(t));
    vSinvepsilonInd(t) = new vector<std::size_t>(nz(t));

    // Predictor step
    vector<std::size_t>& xpInd = *vxpInd(t);
    matrix<std::size_t>& BPuInd = *vBPuInd(t);
    matrix<std::size_t>& BPuBTInd = *vBPuBTInd(t);
    matrix<std::size_t>& Pp1Ind = *vPp1Ind(t);
    matrix<std::size_t>& DSigmawInd = *vDSigmawInd(t);
    matrix<std::size_t>& PpInd = *vPpInd(t);

    // Instrument pricing
    matrix<std::size_t>& oInd = *voInd(t);
    vector<std::size_t>& eoxInd = *veoxInd(t);
    vector<std::size_t>& denInd = *vdenInd(t);
    vector<std::size_t>& gInd = *vgInd(t);
    matrix<std::size_t>& JacobianInd = *vJacobianInd(t);

    // Update step
    vector<std::size_t>& epsilonInd = *vepsilonInd(t);
    matrix<std::size_t>& HPInd = *vHPInd(t);
    matrix<std::size_t>& S1Ind = *vS1Ind(t);
    matrix<std::size_t>& GSigmavInd = *vGSigmavInd(t);
    matrix<std::size_t>& S2Ind = *vS2Ind(t);
    matrix<std::size_t>& SInd = *vSInd(t);
    matrix<std::size_t>& LInd = *vLInd(t);
    matrix<std::size_t>& K1Ind = *vK1Ind(t);
    matrix<std::size_t>& KInd = *vKInd(t);
    vector<std::size_t>& xuInd = *vxuInd(t);
    matrix<std::size_t>& PuInd = *vPuInd(t);
    vector<std::size_t>& Sinvepsilon1Ind = *vSinvepsilon1Ind(t);
    vector<std::size_t>& SinvepsilonInd = *vSinvepsilonInd(t);

    // Setup indexing for automatic differentiation

    // Predictor step
    for (std::size_t i=0 ; i<nx(t) ; i++)
      xpInd(i) = nc++;
    if (t==0)
      for (std::size_t i=0 ; i<nthetaF ; i++)
        for (std::size_t j=0 ; j<nSigma0 ; j++)
          BPuInd(i,j) = nc++;
    else
      for (std::size_t i=0 ; i<nthetaF ; i++)
        for (std::size_t j=0 ; j<nx(t) ; j++)
          BPuInd(i,j) = nc++;
    for (std::size_t i=0 ; i<nthetaF ; i++)
      for (std::size_t j=0 ; j<nthetaF ; j++)
        BPuBTInd(i,j) = nc++;
    for (std::size_t i=0 ; i<nx(t) ; i++)
      for (std::size_t j=0 ; j<nx(t) ; j++)
        Pp1Ind(i,j) = nc++;
    for (std::size_t i=0 ; i<nthetaF ; i++)
      for (std::size_t j=0 ; j<nSigmaw ; j++)
        DSigmawInd(i,j) = nc++;
    for (std::size_t i=0 ; i<nx(t) ; i++)
      for (std::size_t j=0 ; j<nx(t) ; j++)
        PpInd(i,j) = nc++;

    // Instrument pricing
    for (std::size_t i=0 ; i<noT(t) ; i++)
      for (std::size_t j=0 ; j<nxs(t) ; j++)
        oInd(i,j) = nc++;
    for (std::size_t i=0 ; i<noT(t) ; i++)
      eoxInd(i) = nc++;
    for (std::size_t i=0 ; i<nz(t) ; i++)
      denInd(i) = nc++;
    for (std::size_t i=0 ; i<nz(t) ; i++)
      gInd(i) = nc++;
    for (std::size_t i=0 ; i<nz(t) ; i++)
      for (std::size_t j=0 ; j<nxs(t) ; j++)
        JacobianInd(i,j) = nc++;

    // Update step
    for (std::size_t i=0 ; i<nz(t) ; i++)
      epsilonInd(i) = nc++;
    for (std::size_t i=0 ; i<nz(t) ; i++)
      for (std::size_t j=0 ; j<nx(t) ; j++)
        HPInd(i,j) = nc++;
    for (std::size_t i=0 ; i<nz(t) ; i++)
      for (std::size_t j=0 ; j<nz(t) ; j++)
        S1Ind(i,j) = nc++;
    for (std::size_t i=0 ; i<nz(t) ; i++)
      for (std::size_t j=0 ; j<nSigmav ; j++)
        GSigmavInd(i,j) = nc++;
    for (std::size_t i=0 ; i<nz(t) ; i++)
      for (std::size_t j=0 ; j<nz(t) ; j++)
        S2Ind(i,j) = nc++;
    for (std::size_t i=0 ; i<nz(t) ; i++)
      for (std::size_t j=0 ; j<nz(t) ; j++)
        SInd(i,j) = nc++;
    for (std::size_t i=0 ; i<nz(t) ; i++)
      for (std::size_t j=0 ; j<nz(t) ; j++)
        LInd(i,j) = nc++;
    for (std::size_t i=0 ; i<nx(t) ; i++)
      for (std::size_t j=0 ; j<nz(t) ; j++)
        K1Ind(i,j) = nc++;
    for (std::size_t i=0 ; i<nx(t) ; i++)
      for (std::size_t j=0 ; j<nz(t) ; j++)
        KInd(i,j) = nc++;
    for (std::size_t i=0 ; i<nx(t) ; i++)
      xuInd(i) = nc++;
    for (std::size_t i=0 ; i<nx(t) ; i++)
      for (std::size_t j=0 ; j<nx(t) ; j++)
        PuInd(i,j) = nc++;
    for (std::size_t i=0 ; i<nz(t) ; i++)
      Sinvepsilon1Ind(i) = nc++;
    for (std::size_t i=0 ; i<nz(t) ; i++)
      SinvepsilonInd(i) = nc++;
  }

  fInd = nc++;

  std::size_t adMax = 1+20*nt+1;
  vector<AutoDiff<T>*> ad(adMax, NULL);

  std::size_t rMax = 0;
  std::size_t dMax = 0;

  std::size_t na = 0;

  ad(na++) = new AD_a_cumsumBx<T>(365, &Ef, &intEf, &EfInd, &intEfInd); // intEf
  //ad(na-1)->testAD(); // Test if Jacobian is correct

  for (std::size_t t=0 ; t<nt ; t++)
  {
    const SparseMatrix<T>& B = *vB(t);
    const vector<T>& Bv = B.val();
    const SparseMatrix<T>& D = *vD(t);
    const vector<T>& Dv = D.val();

    if (t==0)
    {
      ad(na++) = new AD_adx_Bs_cx<T>(&thetaF, vB(t), &ax, vxp(t), &thetaFInd, &axInd, vxpInd(t));  // 0 xp = diag(thetaF) * B * ax
      ad(na++) = new AD_As_bdx<T>(vB(t), &Sigma0, vBPu(t), &Sigma0Ind, vBPuInd(t));  // 1 BPu = B * diag(Sigma0)
    }
    else
    {
      ad(na++) = new AD_adx_Bs_cx<T>(&thetaF, vB(t), vxu(t-1), vxp(t), &thetaFInd, vxuInd(t-1), vxpInd(t));  // 0 xp = diag(thetaF) * B * vxu(t-1)
      ad(na++) = new AD_As_Bx<T>(vB(t), vPu(t-1), vBPu(t), vPuInd(t-1), vBPuInd(t));  // 1 BPu = B * Pu(t-1)
    }
    ad(na++) = new AD_Ax_BsT<T>(vBPu(t), vB(t), vBPuBT(t), vBPuInd(t), vBPuBTInd(t));  // 2 BPuBT = BPu * B^T
    ad(na++) = new AD_adx_Bx_cdx<T>(&thetaF, vBPuBT(t), &thetaF, vPp1(t), &thetaFInd, vBPuBTInd(t), &thetaFInd, vPp1Ind(t));  // 3 Pp1 = diag(thetaF) * BPuBT * diag(thetaF)
    ad(na++) = new AD_As_bdx<T>(vD(t), &Sigmaw, vDSigmaw(t), &SigmawInd, vDSigmawInd(t));  // 4 DSigmaw = D * diag(sigmaw)
    ad(na++) = new AD_Ax_BsT_Cx<T>(vDSigmaw(t), vD(t), vPp1(t), vPp(t), vDSigmawInd(t), vPp1Ind(t), vPpInd(t));  // 5 Pp = Pp1 + DSigmaw * D^T

    ad(na++) = new AD_InstrPricing<T>(vat(t), voT(t), voI(t), vtc(t), vdatesStepRel(t), &intEf, vxp(t), vo(t), veox(t), vden(t), vg(t), vJacobian(t), &intEfInd, vxpInd(t), voInd(t), veoxInd(t), vdenInd(t), vgInd(t), vJacobianInd(t)); // 6 Instrument pricing g, Jacobian


    const vector<std::size_t>& at = *vat(t);
    const vector<std::size_t>& oI = *voI(t);


    const SparseMatrix<T>& Iz = *vIz(t);
    const vector<T>& Izv = Iz.val();

    ad(na++) = new AD_abx_Cs_cx<T>(vz(t), vg(t), vIz(t), nxs(t), vxp(t), vepsilon(t), vgInd(t), vxpInd(t), vepsilonInd(t)); // 7 epsilon = z - H * xp - g 
    ad(na++) = new AD_AxBs_Cx<T>(vJacobian(t), vIz(t), vPp(t), vHP(t), vJacobianInd(t), vPpInd(t), vHPInd(t));  // 8 HP = H Pp

    ad(na++) = new AD_Ax_BxCsT<T>(vHP(t), vJacobian(t), vIz(t), vS1(t), vHPInd(t), vJacobianInd(t), vS1Ind(t));  // 9 S1 = HP * H^T

    const SparseMatrix<T>& G = *vG(t);
    const vector<T>& Gv = G.val();

    ad(na++) = new AD_As_bdx<T>(vG(t), &Sigmav, vGSigmav(t), &SigmavInd, vGSigmavInd(t));  // 10 GSigmav = G * Sigmav
    ad(na++) = new AD_Ax_BsT_Cx<T>(vGSigmav(t), vG(t), vS1(t), vS2(t), vGSigmavInd(t), vS1Ind(t), vS2Ind(t));  // 11 S2 = S1 + GSigmav * G^T 
    ad(na++) = new AD_Ax_AxT<T>(vS2(t), vS(t), vS2Ind(t), vSInd(t));  // 11 S = (S2+S2^T)/2
    ad(na++) = new AD_cholesky<T>(vS(t), vL(t), vSInd(t), vLInd(t)); // 12 L such that L L^T = S
    ad(na++) = new AD_eqSolveTriang1_YxLxLxT_Ax<T>(vHP(t), vL(t), vK1(t), vHPInd(t), vLInd(t), vK1Ind(t)); // 13 Solve triangular equation system K L L^T = P H - step 1
    ad(na++) = new AD_eqSolveTriang2_YxLxLxT_Ax<T>(vK1(t), vL(t), vK(t), vK1Ind(t), vLInd(t), vKInd(t)); // 14 Solve triangular equation system K L L^T = P H - step 2

    ad(na++) = new AD_ax_Bx_cx<T>(vxp(t), vK(t), vepsilon(t), vxu(t), vxpInd(t), vKInd(t), vepsilonInd(t), vxuInd(t)); // 15 xu = xp + K * epsilon
    ad(na++) = new AD_Ax_Bx_Cx<T>(vPp(t), vK(t), vHP(t), vPu(t), vPpInd(t), vKInd(t), vHPInd(t), vPuInd(t)); // 16 Pu = Pp - K * HP
    ad(na++) = new AD_eqSolveTriang1_LxLxTyx_ax<T>(vepsilon(t), vL(t), vSinvepsilon1(t), vepsilonInd(t), vLInd(t), vSinvepsilon1Ind(t)); // 17 Solve triangular equation system L L^T y = epsilon - step 1
    ad(na++) = new AD_eqSolveTriang2_LxLxTyx_ax<T>(vSinvepsilon1(t), vL(t), vSinvepsilon(t), vSinvepsilon1Ind(t), vLInd(t), vSinvepsilonInd(t)); // 18 Solve triangular equation system L L^T y = epsilon - step 2

  }
  matrix<T> fAll(nt,1);
  ad(na++) = new AD_LogLikelihood<T>(&vL, &vepsilon, &vSinvepsilon, &f, &fAll, &vLInd, &vepsilonInd, &vSinvepsilonInd, &fInd); 

  if (na != ad.size())
  {
    std::cerr << "Incorrect size of ad" << std::endl;
    return;
  }

  for (std::size_t i=0 ; i<ad.size() ; i++)
    if (ad(i) == NULL)
      std::cerr << "Incorrect size of ad" << std::endl;
    else
      ad(i)->count(rMax, dMax);

                      // Remove these later, only here to create equality while debugging
  dMax -= 2450; // xp, BPu
  dMax -= 12323500; // DSigmaw
  dMax -= 3799737; // GSigmav

  //std::cout << rMax << " " << dMax << std::endl;

  std::size_t debugM = 50000000; // Ensures sufficient amount of memory while debugging, remove
  rMax += debugM;
  dMax += debugM;

  // Initialize computations for automatic differentiation

  vector<T> d(dMax, 0.0);
  vector<std::size_t> col(dMax);
  vector<std::size_t> row(rMax+1);
  vector<std::size_t> rInd(rMax);
  std::size_t nr=0;
  std::size_t nd=0;

  /*for (std::size_t i=0 ; i<100 ; i++)
  {
    std::cout << i << " ";
    for (std::size_t j=0 ; j<np ; j++)
      std::cout << Ef(i,j) << " ";
    std::cout << std::endl;
  }*/

  //AD_a_cumsumBx<T> ad(1.0/365, &Ef, &intEf, &EfInd, &intEfInd);
  //ad.fAD(d, col, row, rInd, nr, nd);

  na = 0;
  ad(na++)->fAD(d, col, row, rInd, nr, nd);

  for (std::size_t t=0 ; t<nt ; t++)
  {
    if (t==0)
    {
      ad(na++)->fAD(d, col, row, rInd, nr, nd);  // xp
      ad(na++)->fAD(d, col, row, rInd, nr, nd);  // BPu

    }
    else
    {
      ad(na++)->fAD(d, col, row, rInd, nr, nd);  // xp
      ad(na++)->fAD(d, col, row, rInd, nr, nd);  // BPu
    }

    //pretty_print(std::cout, *vxp(t), "xp", 5, 12);

    ad(na++)->fAD(d, col, row, rInd, nr, nd);  // BPuBT


    // Compute Pp1 = A thetaF B Pu BT thetaF AT
    // Assume that A is identity
    if (nx(t) != nthetaF)
      throw ExcKal("Assumes that A is identity matrix");

    ad(na++)->fAD(d, col, row, rInd, nr, nd);  // Pp1

    ad(na++)->fAD(d, col, row, rInd, nr, nd);  // DSigmaw

    //pretty_print(std::cout, *vDSigmaw(t), "DSigmaw", 5, 12);

    ad(na++)->fAD(d, col, row, rInd, nr, nd);  // Pp

    //pretty_print(std::cout, *vPp(t), "Pp", 5, 12);



    ad(na++)->fAD(d, col, row, rInd, nr, nd);  // Instrument pricing

    // Instrument pricing
    //pretty_print(std::cout, *vtc(t), "tc", 5, 12);

    //pretty_print(std::cout, *vo(t), "o", 5, 12);


    //pretty_print(std::cout, *vg(t), "g", 5, 12);
    //pretty_print(std::cout, *vJacobian(t), "G", 5, 12);

    ad(na++)->fAD(d, col, row, rInd, nr, nd);  // epsilon

    //pretty_print(std::cout, *vepsilon(t), "epsilon", 5, 12);

    ad(na++)->fAD(d, col, row, rInd, nr, nd);  // HP

    //pretty_print(std::cout, *vHP(t), "HP", 5, 12);

    ad(na++)->fAD(d, col, row, rInd, nr, nd);  // S1 = H P HT


    //pretty_print(std::cout, *vS1(t), "HPHT", 5, 12);

    ad(na++)->fAD(d, col, row, rInd, nr, nd);  // GSigmav
    ad(na++)->fAD(d, col, row, rInd, nr, nd);  // S2
    ad(na++)->fAD(d, col, row, rInd, nr, nd);  // S

    //pretty_print(std::cout, *vS(t), "S", 5, 12);

    // Cholesky factorization

    ad(na++)->fAD(d, col, row, rInd, nr, nd);  // L

    //if (t==0)
    //  pretty_print(std::cout, *vL(t), "L", 5, 12);

    ad(na++)->fAD(d, col, row, rInd, nr, nd);  // Solve triangular equation system K L L^T = P H - step 1
    ad(na++)->fAD(d, col, row, rInd, nr, nd);  // Solve triangular equation system K L L^T = P H - step 2

    //pretty_print(std::cout, *vK(t), "K", 5, 12);

    ad(na++)->fAD(d, col, row, rInd, nr, nd);  // xu = xp + K * epsilon

    //pretty_print(std::cout, *vxu(t), "xu", 5, 12);

    ad(na++)->fAD(d, col, row, rInd, nr, nd);  // Pu = Pp - K * HP
    //pretty_print(std::cout, *vPu(t), "Pu", 5, 12);
    //pretty_print(std::cout, *vL(t), "L", 5, 12);

    
    // Solve triangular equation system L L^T y = epsilon - step 1

 
    ad(na++)->fAD(d, col, row, rInd, nr, nd);  // Solve triangular equation system L L^T y = epsilon - step 1
    ad(na++)->fAD(d, col, row, rInd, nr, nd);  // Solve triangular equation system L L^T y = epsilon - step 2

  }
  ad(na++)->fAD(d, col, row, rInd, nr, nd);  // Compute log likelihood value

  row(nr) = nd; // Final column

  retM = fAll;

  /*retM.resize(vxp(0)->size(), 1);
  for (std::size_t i=0 ; i<retM.size1() ; i++)
    retM(i,0) = (*vxp(0))(i);*/

  //retM = *vo(0);
  //retM = *vJacobian(0);
  /*retM.resize(vepsilon(0)->size(), 1);
  for (std::size_t i=0 ; i<retM.size1() ; i++)
    retM(i,0) = (*vepsilon(0))(i);*/

  //retM = *vPu(0);
  /*retM.resize(vL.size(), vL(0)->size1());
  for (std::size_t i=0 ; i<retM.size1() ; i++)
    for (std::size_t j=0 ; j<retM.size2() ; j++)
      retM(i,j) = (*vL(i))(j,j);*/

  //if (nr+debugM != rMax)
  if (nr+debugM > rMax)
    std::cout << "AD: Incorrect number of rows:" << nr << "!=" << rMax-debugM << " diff= " << rMax-debugM-nr << std::endl;
  //if (nd+debugM != dMax)
  if (nd+debugM > dMax)
    std::cout << "AD: Incorrect number of data:" << nd << "!=" << dMax-debugM << " diff= " << dMax-debugM-nd << std::endl;


  // Test derivatives one at a time

  /*for (std::size_t i=1 ; i<ad.size()-1 ; i++)
  {
    std::size_t t = (i-1)/19;
    std::size_t j = i-1 - t*19;
    //if (j==0)
    //if (j==1)
    //if (j==2)
    //if (j==3)
    //if (j==4)
    //if (j==5)
    //if (j==6)
    //if (j==7)
    //if (j==8)
    //if (j==9 && (t%50==0 || t<100))
    //if (j==10 && (t%50==0 || t<100)) 
    //if (j==11 && (t%50==0 || t<100)) 
    //if (j==12 && (t%50==0 || t<100))
    //if (j==13 && (t%50==0 || t<100))
    //if (j==14 && (t%50==0 || t<100))
    //if (j==15 && (t%50==0 || t<100))
    //if (j==16 && (t%50==0 || t<100))
    //if (j==17 && (t%50==0 || t<100))
    if (j==18 && (t%50==0 || t<100))
    {
      T norm = ad(i)->testAD(); // Test if Jacobian is correct
      if (norm <= 1E-5)
        ;
      else
        std::cout << t << " " << norm << std::endl;    
    }
  }*/

//  T norm = ad(ad.size()-1)->testAD(); // Test if Jacobian is correct
//  std::cout << "likelihood norm " << norm << std::endl;    

  // Run automatic differentiation

  vector<T> dOdx(nc, 0.0);

  // Debug AD
  int testCase = 0;
//std::cout << "test = " << testCase << " ";
  if (testCase == 0) // This computes the gradient of the log-likelihood function
  {
    dOdx(fInd) = 1;
    obj = f;
  }
  else if (testCase == 1)
  {
    dOdx(intEfInd(nT, 0)) = 1;
    obj = intEf(nT, 0);
  }
  else if (testCase == 2)
  {
    dOdx((*vxpInd(0))(0)) = 1;
    obj = (*vxp(0))(0);
  }
  else if (testCase == 3)
  {
    dOdx((*vdenInd(0))(1)) = 1;
    obj = (*vden(0))(1);
  }
  else if (testCase == 4)
  {
    dOdx((*veoxInd(0))(2)) = 1;
    obj = (*veox(0))(2);
  }
  else if (testCase == 5)
  {
    dOdx((*vgInd(0))(2)) = 1;
    obj = (*vg(0))(2);
  }
  else if (testCase == 6)
  {
    dOdx((*vJacobianInd(0))(1,0)) = 1;
    obj = (*vJacobian(0))(1,0);
  }
  else if (testCase == 7)
  {
    dOdx((*vepsilonInd(0))(0)) = 1;
    obj = (*vepsilon(0))(0);
  }
  else if (testCase == 8)
  {
    dOdx((*vSInd(0))(1,1)) = 1;
    obj = (*vS(0))(1,1);
  }
  else if (testCase == 9)
  {
    dOdx((*vLInd(0))(2,1)) = 1;
    obj = (*vL(0))(2,1);
  }
  else if (testCase == 10)
  {
    dOdx((*vxuInd(0))(2)) = 1;
    obj = (*vxu(0))(2);
  }
  else if (testCase == 11)
  {
    dOdx((*vPuInd(0))(2,1)) = 1;
    obj = (*vPu(0))(2,1);
  }
  else if (testCase == 12)
  {
    dOdx((*vSinvepsilonInd(0))(2)) = 1;
    obj = (*vSinvepsilon(0))(2);
  }


  for (std::size_t i=nr ; i-- > 0 ; )
  {
    const T v1 = dOdx(rInd(i));
    if (v1 != 0)
    {
//std::cerr << i << " " << rInd(i) << " : " << row(i) << " ->" << row(i+1) << " mult " << v1 << std::endl; 
      for (std::size_t j=row(i) ; j<row(i+1) ; j++)
      {  
        T& v2 = dOdx(col(j));
//std::cerr << xn(rInd(i)) << " " << rInd(i) << " : "  << xn(col(j)) << " = " << v2 << " + " << v1*d(j) << " = " << v2 + v1*d(j) << std::endl; 
        if (v2==0)
          v2 = v1*d(j);
        else
          v2 += v1*d(j);
      }
    }
  }

  gradObjEf.resize(nT, np);
  for (std::size_t i=0 ; i<nT ; i++)
    for (std::size_t j=0 ; j<np ; j++)
      gradObjEf(i,j) = dOdx(EfInd(i,j));

  gradObjthetaF.resize(nthetaF);
  for (std::size_t i=0 ; i<nthetaF ; i++)
    gradObjthetaF(i) = dOdx(thetaFInd(i));

  gradObjSigmav.resize(nSigmav);
  for (std::size_t i=0 ; i<nSigmav ; i++)
    gradObjSigmav(i) = dOdx(SigmavInd(i));

  gradObjSigmaw.resize(nSigmaw);
  for (std::size_t i=0 ; i<nSigmaw ; i++)
    gradObjSigmaw(i) = dOdx(SigmawInd(i));

  gradObjax.resize(nax);
  for (std::size_t i=0 ; i<nax ; i++)
    gradObjax(i) = dOdx(axInd(i));

  gradObjSigma0.resize(nSigma0);
  for (std::size_t i=0 ; i<nSigma0 ; i++)
    gradObjSigma0(i) = dOdx(Sigma0Ind(i));


  for (std::size_t i=0 ; i<ad.size() ; i++)
    if (ad(i) == NULL)
      std::cerr << "Incorrect size of ad" << std::endl;
    else
      delete ad(i);

  for (std::size_t t=0 ; t<nt ; t++)
  {
    // Predictor step
    delete vxp(t);
    delete vBPu(t);
    delete vBPuBT(t);
    delete vPp1(t);
    delete vDSigmaw(t);
    delete vPp(t);

    // Instrument pricing
    delete vo(t);
    delete veox(t);
    delete vden(t);
    delete vg(t);
    delete vJacobian(t);

    // Update step
    delete vepsilon(t);
    delete vHP(t);
    delete vS1(t);
    delete vGSigmav(t);
    delete vS2(t);
    delete vS(t);
    delete vL(t);
    delete vK1(t);
    delete vK(t);
    delete vxu(t);
    delete vPu(t);
    delete vSinvepsilon1(t);
    delete vSinvepsilon(t);

    // Memory allocation of index variables for automatic differentiation

    // Predictor step
    delete vxpInd(t);
    delete vBPuInd(t);
    delete vBPuBTInd(t);
    delete vPp1Ind(t);
    delete vDSigmawInd(t);
    delete vPpInd(t);

    // Instrument pricing
    delete voInd(t);
    delete veoxInd(t);
    delete vdenInd(t);
    delete vgInd(t);
    delete vJacobianInd(t);

    // Update step
    delete vepsilonInd(t);
    delete vHPInd(t);
    delete vS1Ind(t);
    delete vGSigmavInd(t);
    delete vS2Ind(t);
    delete vSInd(t);
    delete vLInd(t);
    delete vK1Ind(t);
    delete vKInd(t);
    delete vxuInd(t);
    delete vPuInd(t);
    delete vSinvepsilon1Ind(t);
    delete vSinvepsilonInd(t);
  }

}


#endif // KALMANLOGLIKELIHOODAD_H




