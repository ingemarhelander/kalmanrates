//#include <windows.h>
#include <iomanip>
#include "mex.hpp"
#include "mexAdapter.hpp"
#include<string>
#include<memory>

#include "boost/numeric/ublas/matrix.hpp"
#include "boost/numeric/ublas/io.hpp"
#include "boost/lexical_cast.hpp"
//#include "XiMath/boostHelper.h"
//#include "XiMath/sparseMatrix.h"
#include "XiFinPortfolio/kalmanLogLikelihoodAd.h"
//#include "mexConvertTypes.hpp"


using namespace boost::numeric::ublas;
using namespace XiMath;
using namespace XiFin;
using namespace matlab::mex;
using namespace matlab::data;

class MexFunction : public Function {
private:
  std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr;
public:
  /* Constructor for the class. */
  MexFunction()
  {
    matlabPtr = getEngine();
  }
  
  /* Helper function to print output string on MATLAB command prompt. */
  void displayOnMATLAB(std::ostringstream stream)
  {
    ArrayFactory factory;
    matlabPtr->feval(u"fprintf", 0, std::vector<Array>
            ({ factory.createScalar(stream.str())}));
  }
  
  /* Helper function to generate an error message from given string,
   * and display it over MATLAB command prompt.
   */
  void displayError(std::string errorMessage)
  {
    ArrayFactory factory;
    matlabPtr->feval(u"error", 0, std::vector<Array>({
      factory.createScalar(errorMessage) }));
  }
  
  /* Helper function to information about an empty field in the structure. */
  void emptyFieldInformation(std::string fieldName, size_t index)
  {
    std::ostringstream stream;
    stream<<"Field: "<<std::string(fieldName)<<" of the element at index: "
          <<index+1<<" is empty."<<std::endl;
    displayOnMATLAB(std::move(stream));
  }
  
  /* Helper function to information about an invalid field in the structure. */
  void invalidFieldInformation(std::string fieldName, size_t index)
  {
    std::ostringstream stream;
    stream<<"Field: "<<std::string(fieldName)<<" of the element at index: "
          <<index+1<<" contains wrong value."<<std::endl;
    displayOnMATLAB(std::move(stream));
  }
  
  template<class T>
  void createVector(const TypedArray<T>& inArray, vector<T>& v)
  {
    v.resize(inArray.getNumberOfElements());
    std::size_t j=0;
    for (auto& elem : inArray) 
      v(j++) = elem;
  }

  template<class T>
  void createMatrix(const TypedArray<T>& inArray, matrix<T>& mat)
  {
    ArrayDimensions dim = inArray.getDimensions();
    std::size_t m = dim[0];
    std::size_t n = dim[1];

    mat.resize(m,n);

    std::size_t k=0;
    for (auto& elem : inArray) 
    {
      std::size_t j = k/m;
      std::size_t i = k-m*j;
      mat(i,j) = elem;
      k++;
    }
  }

  template<class T>
  void createVecVector(CellArray& imageCell, vector<vector<T>*>& vv)
  {
    vv.resize(imageCell.getNumberOfElements());
    for (std::size_t i=0 ; i<imageCell.getNumberOfElements() ; i++)
    {
      ArrayDimensions dim = imageCell[i].getDimensions();
      if (dim[1] != 1)
        displayError("Should be a column vector");
      vv(i) = new vector<T>(dim[0]);
      vector<T>& v = *vv(i);

      const TypedArray<double> inArray = imageCell[i];
      std::size_t j=0;
      for (auto& elem : inArray) 
        v(j++) = elem;

      // std::cout << dim[0] << "x" << dim[1] << std::endl;
    }
  }

  template<class T>
  void createVecMatrix(CellArray& imageCell, vector<matrix<T>*>& vm)
  {
    vm.resize(imageCell.getNumberOfElements());
    for (std::size_t i=0 ; i<imageCell.getNumberOfElements() ; i++)
    {
      ArrayDimensions dim = imageCell[i].getDimensions();
      std::size_t m = dim[0];
      std::size_t n = dim[1];

      vm(i) = new matrix<T>(m, n);
      matrix<T>& v = *vm(i);

      const TypedArray<double> inArray = imageCell[i];
      std::size_t k=0;
      for (auto& elem : inArray) 
      {
        std::size_t i = k/m;
        std::size_t j = k-m*i;
        v(i,j) = elem;
        k++;
      }

      // std::cout << dim[0] << "x" << dim[1] << std::endl;
    }
  }

  template<class T>
  void createVecSparseMatrix(CellArray& imageCell, vector<SparseMatrix<T>*>& vm)
  {
    // vr = row, vc = column, vv = values
    // length of vr[i], vc[i] and vv[i] are the same 
    vm.resize(imageCell.getNumberOfElements());
    for (std::size_t i=0 ; i<imageCell.getNumberOfElements() ; i++)
    {
      ArrayDimensions dim = imageCell[i].getDimensions();
      SparseArray<T> spArrayInp = std::move(imageCell[i]);
      std::size_t nz = spArrayInp.getNumberOfNonZeroElements();

      vector<std::size_t> r(nz);
      vector<std::size_t> c(nz);
      vector<T> v(nz);
      auto idx = spArrayInp.getIndex(spArrayInp.begin());
      std::size_t itr=0;
      for (TypedIterator<T> inpSpArrayPtr = spArrayInp.begin(); inpSpArrayPtr != spArrayInp.end(); ++inpSpArrayPtr)
      {
          idx = spArrayInp.getIndex(inpSpArrayPtr);
          r(itr) = idx.first;
          c(itr) = idx.second;
          v(itr++) = *inpSpArrayPtr;
//std::cout << r(itr-1) << " " << c(itr-1) << " " << v(itr-1) << " " << std::endl;
      }

      vm(i) = new SparseMatrix<T>;
      sparseMatlab2sparse(r, c, v, dim[0], dim[1], *vm(i));
      /*matrix<T> B;
      sparse2matrix(*vm(i), B);
      for (std::size_t i=0 ; i<B.size1() ; i++)
      {
        for (std::size_t j=0 ; j<B.size2() ; j++)
          std::cout << std::setw(10) << B(i,j) << " ";
        std::cout << std::endl;
      }*/
      // std::cout << dim[0] << "x" << dim[1] << std::endl;
    }
  }

  /* This is the gateway routine for the MEX-file. */
  void 
  operator()(ArgumentList outputs, ArgumentList inputs) { 
    
    // z{t} A{t} B{t} D{t} G{t} Iz{t} thetaF Sigmav Sigmaw ax Sigma0 
    // at{t} oT{t} oInd{t} tc{t} E nEigs

    if (inputs.size() != 17)
      displayError("Seventeen inputs required.");
    if (outputs.size() != 1 && outputs.size() != 9)
      displayError("Too many outputs specified.");

    if (inputs[0].getType() != ArrayType::CELL)
      displayError("Input 1 must be a cell.");
    CellArray imageCell = std::move(inputs[0]);
    vector<vector<double>*> z;
    createVecVector<double>(imageCell, z);

    if (inputs[1].getType() != ArrayType::CELL)
      displayError("Input 2 must be a cell.");
    imageCell = std::move(inputs[1]);
    vector<SparseMatrix<double>*> A;
    createVecSparseMatrix<double>(imageCell, A);

    if (inputs[2].getType() != ArrayType::CELL)
      displayError("Input 3 must be a cell.");
    imageCell = std::move(inputs[2]);
    vector<SparseMatrix<double>*> B;
    createVecSparseMatrix<double>(imageCell, B);

    if (inputs[3].getType() != ArrayType::CELL)
      displayError("Input 4 must be a cell.");
    imageCell = std::move(inputs[3]);
    vector<SparseMatrix<double>*> D;
    createVecSparseMatrix<double>(imageCell, D);

    if (inputs[4].getType() != ArrayType::CELL)
      displayError("Input 5 must be a cell.");
    imageCell = std::move(inputs[4]);
    vector<SparseMatrix<double>*> G;
    createVecSparseMatrix<double>(imageCell, G);

    if (inputs[5].getType() != ArrayType::CELL)
      displayError("Input 6 must be a cell.");
    imageCell = std::move(inputs[5]);
    vector<SparseMatrix<double>*> Iz;
    createVecSparseMatrix<double>(imageCell, Iz);

    if (inputs[6].getType() != ArrayType::DOUBLE)
      displayError("Input 7 must be a double.");
    TypedArray<double> imageDouble = std::move(inputs[6]);
    vector<double> thetaF;
    createVector(imageDouble, thetaF);

    /*for (std::size_t i=0 ; i<thetaF.size() ; i++)
      std::cout << thetaF(i) << " ";
    std::cout << std::endl;*/

    if (inputs[7].getType() != ArrayType::DOUBLE)
      displayError("Input 8 must be a double.");
    imageDouble = std::move(inputs[7]);
    vector<double> Sigmav;
    createVector<double>(imageDouble, Sigmav);

    if (inputs[8].getType() != ArrayType::DOUBLE)
      displayError("Input 9 must be a double.");
    imageDouble = std::move(inputs[8]);
    vector<double> Sigmaw;
    createVector<double>(imageDouble, Sigmaw);

    if (inputs[9].getType() != ArrayType::DOUBLE)
      displayError("Input 10 must be a double.");
    imageDouble = std::move(inputs[9]);
    vector<double> ax;
    createVector<double>(imageDouble, ax);
  
    if (inputs[10].getType() != ArrayType::DOUBLE)
      displayError("Input 11 must be a double.");
    imageDouble = std::move(inputs[10]);
    matrix<double> Sigma0;
    createMatrix<double>(imageDouble, Sigma0);

    if (inputs[11].getType() != ArrayType::CELL)
      displayError("Input 12 must be a cell.");
    imageCell = std::move(inputs[11]);
    vector<vector<std::size_t>*> at;
    createVecVector<std::size_t>(imageCell, at);

    if (inputs[12].getType() != ArrayType::CELL)
      displayError("Input 13 must be a cell.");
    imageCell = std::move(inputs[12]);
    vector<vector<std::size_t>*> oT;
    createVecVector<std::size_t>(imageCell, oT);
    for (std::size_t t=0 ; t<oT.size() ; t++) // Change from Matlabs 1-based indexing to 0-based indexing
    {
      vector<std::size_t>& tmp = *oT(t);
      for (std::size_t i=0 ; i<tmp.size() ; i++)
        tmp(i)--;
    }


    if (inputs[13].getType() != ArrayType::CELL)
      displayError("Input 14 must be a cell.");
    imageCell = std::move(inputs[13]);
    vector<vector<std::size_t>*> oInd;
    createVecVector<std::size_t>(imageCell, oInd);
    for (std::size_t t=0 ; t<oInd.size() ; t++) // Change from Matlabs 1-based indexing to 0-based indexing
    {
      vector<std::size_t>& tmp = *oInd(t);
      for (std::size_t i=0 ; i<tmp.size() ; i++)
        tmp(i)--;
    }


    if (inputs[14].getType() != ArrayType::CELL)
      displayError("Input 15 must be a cell.");
    imageCell = std::move(inputs[14]);
    vector<vector<double>*> tc;
    createVecVector<double>(imageCell, tc);

    if (inputs[15].getType() != ArrayType::DOUBLE)
      displayError("Input 16 must be a double.");
    imageDouble = std::move(inputs[15]);
    matrix<double> Ef;
    createMatrix<double>(imageDouble, Ef);
  
    if (inputs[16].getType() != ArrayType::CELL)
      displayError("Input 17 must be a cell.");
    imageCell = std::move(inputs[16]);
    vector<vector<std::size_t>*> datesStepRel;
    createVecVector<std::size_t>(imageCell, datesStepRel);

    double f, obj;
    matrix<double> gradObjEf;
    vector<double> gradObjthetaF;
    vector<double> gradObjSigmav;
    vector<double> gradObjSigmaw;
    vector<double> gradObjax;
    matrix<double> gradObjSigma0;
    matrix<double> retM;
    try {
    /*gradKalmanLogLikelihood(z, A, B, D, G, Iz,
                            thetaF, Sigmav, Sigmaw, ax, Sigma0,
                            at, oT, oInd, tc, Ef, datesStepRel,
											      f, obj, gradObjEf, gradObjthetaF, gradObjSigmav, gradObjSigmaw, gradObjax, gradObjSigma0, retM);*/
    bool calcGrad;
    if (outputs.size() == 1)
      calcGrad = false;
    else
      calcGrad = true;
    KalmanData<double> kd(z, A, B, D, G, Iz,
                  thetaF, Sigmav, Sigmaw, ax, Sigma0,
                  at, oT, oInd, tc, Ef, datesStepRel,
									calcGrad);
    gradKalmanLogLikelihood(kd);


    f = kd.f;
    obj = kd.obj;
    gradObjEf = kd.gradObjEf;
    gradObjthetaF = kd.gradObjthetaF;
    gradObjSigmav = kd.gradObjSigmav;
    gradObjSigmaw = kd.gradObjSigmaw;
    gradObjax = kd.gradObjax;
    gradObjSigma0 = kd.gradObjSigma0;
    retM = kd.retM;
    }
    catch (ExcKal& e)
    {
      e.print();
    }

    for (std::size_t t=0 ; t<z.size() ; t++)
    {
      delete z(t);
      delete A(t);
      delete B(t);
      delete D(t);
      delete G(t);
      delete Iz(t);
      delete at(t);
      delete oT(t);
      delete oInd(t);
      delete tc(t);
      delete datesStepRel(t);
    }

    ArrayFactory factory;

    outputs[0] = factory.createScalar<double>(f);

    if (outputs.size() == 9)
    {
      outputs[1] = factory.createScalar<double>(obj);

      std::size_t m = gradObjEf.size1();
      std::size_t n = gradObjEf.size2();
      matlab::data::TypedArray<double> outputMatrix1 = factory.createArray<double>({m, n});
      for (std::size_t i=0 ; i<m ; i++)
        for (std::size_t j=0 ; j<n; j++)
          outputMatrix1[i][j] = gradObjEf(i,j);
      outputs[2] = outputMatrix1;

      m = gradObjthetaF.size();
      n = 1;
      matlab::data::TypedArray<double> outputMatrix2 = factory.createArray<double>({m, n});
      for (std::size_t i=0 ; i<m ; i++)
        for (std::size_t j=0 ; j<n; j++)
          outputMatrix2[i][j] = gradObjthetaF(i);
      outputs[3] = outputMatrix2;

      m = gradObjSigmav.size();
      n = 1;
      matlab::data::TypedArray<double> outputMatrix3 = factory.createArray<double>({m, n});
      for (std::size_t i=0 ; i<m ; i++)
        for (std::size_t j=0 ; j<n; j++)
          outputMatrix3[i][j] = gradObjSigmav(i);
      outputs[4] = outputMatrix3;

      m = gradObjSigmaw.size();
      n = 1;
      matlab::data::TypedArray<double> outputMatrix4 = factory.createArray<double>({m, n});
      for (std::size_t i=0 ; i<m ; i++)
        for (std::size_t j=0 ; j<n; j++)
          outputMatrix4[i][j] = gradObjSigmaw(i);
      outputs[5] = outputMatrix4;

      m = gradObjax.size();
      n = 1;
      matlab::data::TypedArray<double> outputMatrix5 = factory.createArray<double>({m, n});
      for (std::size_t i=0 ; i<m ; i++)
        for (std::size_t j=0 ; j<n; j++)
          outputMatrix5[i][j] = gradObjax(i);
      outputs[6] = outputMatrix5;

      m = gradObjSigma0.size1();
      n = gradObjSigma0.size2();
      matlab::data::TypedArray<double> outputMatrix6 = factory.createArray<double>({m, n});
      for (std::size_t i=0 ; i<m ; i++)
        for (std::size_t j=0 ; j<n; j++)
          outputMatrix6[i][j] = gradObjSigma0(i,j);
      outputs[7] = outputMatrix6;

      m = retM.size1();
      n = retM.size2();
      matlab::data::TypedArray<double> outputMatrix7 = factory.createArray<double>({m, n});
      for (std::size_t i=0 ; i<m ; i++)
        for (std::size_t j=0 ; j<n; j++)
          outputMatrix7[i][j] = retM(i,j);
      outputs[8] = outputMatrix7;
    }    
  }
  
};

