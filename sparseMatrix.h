#ifndef SPARSEMATRIX_H
#define SPARSEMATRIX_H

#include <iomanip>

#include "boost/numeric/ublas/matrix.hpp"
#include "boost/numeric/ublas/vector.hpp"

using namespace boost::numeric::ublas;

namespace XiFin {


template<class T>
class SparseMatrix 
{
public:
	SparseMatrix() : m_(0), n_(0) { }
	SparseMatrix(std::size_t _m, std::size_t _n) : m_(_m), n_(_n), row_(_m+1) { }
  virtual ~SparseMatrix() { }
  
  SparseMatrix<T>& operator=(const SparseMatrix<T>& A)
  {
    m_ = A.size1();
    n_ = A.size2();
    row_ = A.row();
    col_ = A.col();
    val_ = A.val();

    return *this;
  }

  void resize(std::size_t _m, std::size_t _n) 
  { 
    m_ = _m;
    n_ = _n;
    row_.resize(_m+1);
  }
  inline const std::size_t& m() const { return m_; }
  inline const std::size_t& n() const { return n_; }
  inline const std::size_t& size1() const { return m_; }
  inline const std::size_t& size2() const { return n_; }
  inline vector<std::size_t>& row()  { return row_; }
  inline vector<std::size_t>& col()  { return col_; }
  inline vector<T>& val()             { return val_; }

  inline const vector<std::size_t>& row() const { return row_; }
  inline const vector<std::size_t>& col() const { return col_; }
  inline const vector<T>& val()            const { return val_; }

  void printWhichElement(const T* p)
  {
    for (std::size_t i=0 ; i<m_ ; i++)
      for (std::size_t j=row_(i) ; j<row_(i+1) ; j++)
        if (&val_(j) == p) // Slow implementation
          std::cout << "(" << i << "," << col_(j) << ")";
  }
  void checkValid()
  {
    if (row_(m_) != col_.size())
    {
      std::cout << "Inorrect last row index" << std::endl;
      assert(0);
    }
    for (std::size_t i=0 ; i<m_ ; i++)
    {
      if (row_(i) >= col_.size())
      {
        std::cout << "Inorrect row index" << std::endl;
        assert(0);
      }
      vector<std::size_t> vec(n_, 0);
      for (std::size_t k=row_(i) ; k<row_(i+1) ; k++)
      {
        std::size_t j = col_(k);
        if (vec(j) == 1)
        {
          std::cout << "In row " << i << " column " << j << " has been set several times - not possible in sparse matrix" << std::endl;
          assert(0);        
        }
        vec(j) = 1;
      }
    }
  }

  void initLastInd(vector<std::size_t>& ind) const
  { // Used for efficient loop over all non-zero columns in reverse order
    ind.resize(m_);
    for (std::size_t i=0 ; i<m_ ; i++)
      if (row_(i) == row_(i+1))
        ind(i) = std::numeric_limits<std::size_t>::max();
      else
        ind(i) = row_(i+1)-1;
  }

  void findPrevInd(const std::size_t& i, vector<std::size_t>& ind) const
  { // Used for efficient loop over all non-zero columns in reverse order
    ind.resize(m_);
    if (ind(i) == std::numeric_limits<std::size_t>::max())
      return;
    if (ind(i) == row_(i))
      ind(i) = std::numeric_limits<std::size_t>::max();
    else
      ind(i)--;
  }

  const std::size_t& col(const std::size_t& ind) const { return col_(ind); }

protected:
  std::size_t m_;
  std::size_t n_;
  vector<std::size_t> row_; // This stores the first position of each row
  vector<std::size_t> col_; 
  vector<T> val_; 
};

template<class T>
void sparse2matrix(const SparseMatrix<T>& As, matrix<T>& A)
{
  const vector<std::size_t>& row = As.row();
  const vector<std::size_t>& col = As.col();
  const vector<T>& val = As.val();
  A = zero_matrix<T>(As.m(),As.n());
  for (std::size_t i=0 ; i<row.size()-1 ; i++)
    for (std::size_t j=row(i) ; j<row(i+1) ; j++)
      A(i,col(j)) = val(j);
}

template<class T>
void sparseMatlab2sparse(const vector<std::size_t>& r, const vector<std::size_t>& c, const vector<T>& v, const std::size_t m, const std::size_t n, SparseMatrix<T>& Bs)
{ // Matlab stores all row indices + by column - works as long as sorted by either rows or columns
  std::size_t nz = r.size();
  if (c.size() != nz || v.size() != nz)
    std::cerr << "sparseMatlab2sparse: All vectors should have the same size" << std::endl;

  Bs.resize(m, n);
  vector<std::size_t>& Brow = Bs.row();
  vector<std::size_t>& Bcol = Bs.col();
  vector<T>& Bval = Bs.val();

  Brow.resize(Bs.m()+1);
  Bcol.resize(nz);
  Bval.resize(nz);

  vector<std::size_t> nElementsInRow(m, 0);
  for (std::size_t i=0 ; i<nz ; i++)
    nElementsInRow(r(i))++;

  Brow(0) = 0;
  for (std::size_t i=0 ; i<Bs.m() ; i++)
    Brow(i+1) = Brow(i) + nElementsInRow(i);

  //nElementsInRow = zero_vector<std::size_t>(m);
  for (std::size_t i=0 ; i<m ; i++) // zero_vector not working...
    nElementsInRow(i) = 0;

  for (std::size_t i=0 ; i<nz ; i++)
  {
    std::size_t ind_ = Brow(r(i)) + nElementsInRow(r(i));
    Bcol(ind_) = c(i);
    Bval(ind_) = v(i);
    nElementsInRow(r(i))++;
  }
}


template<class T>
void transposeSparseMatrix(const SparseMatrix<T>& As, SparseMatrix<T>& Bs)
{
  const vector<std::size_t>& Arow = As.row();
  const vector<std::size_t>& Acol = As.col();
  const vector<T>& Aval = As.val();

  Bs.resize(As.n(), As.m());
  vector<std::size_t>& Brow = Bs.row();
  vector<std::size_t>& Bcol = Bs.col();
  vector<T>& Bval = Bs.val();

  Brow.resize(Bs.m()+1);
  std::size_t nz = Aval.size();
  Bcol.resize(nz);
  Bval.resize(nz);
  //ind.resize(nz);

  vector<std::size_t> nElementsInCol(As.n(), 0);
  for (std::size_t i=0 ; i<nz ; i++)
    nElementsInCol(Acol(i))++;

  Brow(0) = 0;
  for (std::size_t i=0 ; i<Bs.m() ; i++)
    Brow(i+1) = Brow(i) + nElementsInCol(i);

  //nElementsInCol = zero_vector<std::size_t>(As.n());
  for (std::size_t i=0 ; i<As.n() ; i++) // zero_vector not working...
    nElementsInCol(i) = 0;

  for (std::size_t i=0 ; i<Arow.size()-1 ; i++)
    for (std::size_t j=Arow(i) ; j<Arow(i+1) ; j++)
    {
      std::size_t ind_ = Brow(Acol(j)) + nElementsInCol(Acol(j));
      Bcol(ind_) = i;
      Bval(ind_) = Aval(j);
      nElementsInCol(Acol(j))++;
    }
}

template<class T>
void eliminateZeros(SparseMatrix<T>& As)
{
  vector<std::size_t>& Arow = As.row();
  vector<std::size_t>& Acol = As.col();
  vector<T>& Aval = As.val();

  std::size_t nz=0;
  vector<std::size_t> row(As.size1()+1);
  row(0) = 0;
  for (std::size_t i=0 ; i<As.size1() ; i++)
  {
    for (std::size_t j=Arow(i) ; j<Arow(i+1) ; j++)
      if (Aval(j) == 0.0)
        ; // Keep NaN, etc.
      else if (j != nz)
      {
        Aval(nz) = Aval(j);
        Acol(nz++) = Acol(j);
      }
    row(i+1) = nz;
  }
  Arow = row;
  Acol.resize(nz, true);
  Aval.resize(nz, true);
}


template<class T>
void transposeSparseMatrixInd(const SparseMatrix<T>& As, SparseMatrix<T>& Bs, vector<std::size_t>& ind)
{ // This function only builds up the structure to be able to transpose matrix As
  // for (std::size_t i=0 ; i<As.val().size() ; i++)
  //   Bs.val()(ind(i)) = As.val()(i);
  const vector<std::size_t>& Arow = As.row();
  const vector<std::size_t>& Acol = As.col();
  const vector<T>& Aval = As.val();

  Bs.resize(As.n(), As.m());
  vector<std::size_t>& Brow = Bs.row();
  vector<std::size_t>& Bcol = Bs.col();
  vector<T>& Bval = Bs.val();

  Brow.resize(Bs.m()+1);
  std::size_t nz = Aval.size();
  Bcol.resize(nz);
  Bval.resize(nz);
  ind.resize(nz);

  vector<std::size_t> nElementsInCol(As.n(), 0);
  for (std::size_t i=0 ; i<nz ; i++)
    nElementsInCol(Acol(i))++;

  Brow(0) = 0;
  for (std::size_t i=0 ; i<Bs.m() ; i++)
    Brow(i+1) = Brow(i) + nElementsInCol(i);

  //nElementsInCol = zero_vector<std::size_t>(As.n());
  for (std::size_t i=0 ; i<As.n() ; i++) // zero_vector not working...
    nElementsInCol(i) = 0;

  for (std::size_t i=0 ; i<Arow.size()-1 ; i++)
    for (std::size_t j=Arow(i) ; j<Arow(i+1) ; j++)
    {
      std::size_t ind_ = Brow(Acol(j)) + nElementsInCol(Acol(j));
      Bcol(ind_) = i;
      Bval(ind_) = Aval(j);
      ind(j) = ind_;
      nElementsInCol(Acol(j))++;
    }
}



} // namespace XiFin

#endif // SPARSEMATRIX_H
