#include "Matrix.h"
#include "MatrixException.h"
#include <string>

#define DIMENSION_INVALID "Invalid matrix dimensions.\n"
#define DIV_BY_ZERO "Division by zero.\n"
#define RANGE_INVALID "Index out of range.\n"
#define STREAM_INVALID "Error loading from input stream.\n"

/*
 * A helper function to (Matrix &m) +, += operators.
 */
bool MatricesAreInSameDimensions (const Matrix &m_1,
                                          const Matrix &m_2) {
  return (m_1.GetRows() == m_2.GetRows() &&
      m_1.GetCols() == m_2.GetCols());
}

/*
 * A helper function to (Matrix &m) *, *= operators. 
 */
bool MatricesCanBeMultiplied (const Matrix &lhs, const Matrix &rhs) {
  return lhs.GetCols() == rhs.GetRows();
}

/*
 * See header.
 */
Matrix & Matrix::Vectorize() noexcept{
  this->rows_ = size_;
  this->cols_ = 1;
  return *this;
}

/*
 * See header.
 */
void Matrix::Print() const noexcept{
 std::cout << *this;
}

/*
 * See header.
 */
void Matrix::Init() noexcept{
  for (int i=0; i<this->size_; ++i){
    this->mat_[i] = DEFAULT_VAL;
  }
}

/*
 * See header.
 */
void Matrix::CopyMat(const Matrix &m, const float add_to_each,
                                const float multiply_by) noexcept{
  for (int i=0; i<this->size_; ++i){
    this->mat_[i] = add_to_each + (multiply_by * m.mat_[i]);
  }
}

/*
 *  A helper function to =, != operators. 
 */
bool ArraysCompare(float const *arr_1, float const *arr_2, int size) noexcept{
  for (int i=0; i<size;++i){
    if (arr_1[i] != arr_2[i]) {
      return false;
    }
  }
  return true;
}

/*
 *  A helper function to *, *= (Matrix &m) operators.
 */
float MultCellValue(const Matrix &lhs,
                    const Matrix &rhs, int row, int col) noexcept{
  float res = 0;
  for(int k=0; k < lhs.GetCols(); ++k){
    res += lhs(row, k)*rhs(k,col);
  }
  return res;
}

/*
 * See header.
 */
bool Matrix::operator==(const Matrix &rhs) const noexcept{
  return !(*this != rhs);
}

/*
 * See header.
 */
bool Matrix::operator!=(const Matrix &rhs) const noexcept{
  return !(MatricesAreInSameDimensions(*this, rhs)) ||
      !ArraysCompare(this->mat_, rhs.mat_, this->size_);
}

/*
 * See header.
 */
Matrix &Matrix::operator=(const Matrix &rhs) noexcept(false){
  if (this != &rhs){
    delete[] mat_;
    this->rows_ = rhs.rows_;
    this->cols_ = rhs.cols_;
    this->size_ = rhs.size_;
    this->mat_ = new float[rhs.size_];
    CopyMat(rhs,NO_ADD,NO_MULTIPLICATION);
  }
  return *this;
}

/*
 * See header.
 */
Matrix Matrix::operator*(const Matrix &rhs) const noexcept(false) {
  if (!MatricesCanBeMultiplied(*this, rhs)){
    throw MatrixException(DIMENSION_INVALID);
  }
  Matrix m = Matrix(this->rows_,rhs.cols_);
  for (int i=0; i < this->rows_; ++i){
    for (int j=0; j < rhs.cols_; ++j){
      m(i,j) = MultCellValue(*this, rhs, i, j);
    }
  }
  return m;
}

/*
 * See header.
 */
Matrix &Matrix::operator*=(const Matrix &rhs) noexcept(false) {
  if (!MatricesCanBeMultiplied(*this, rhs)){
    throw MatrixException(DIMENSION_INVALID);
  }
  Matrix mult = ((*this)*rhs);
  *this = mult;
  return *this;
}

/*
 * See header.
 */
Matrix Matrix::operator*(float scalar) const noexcept{
  return Matrix(*this, NO_ADD, scalar);
}

/*
 * See header.
 */
Matrix operator*(const float scalar, const Matrix &m) noexcept{
  return m*scalar;
}

/*
 * See header.
 */
Matrix &Matrix::operator*=(const float scalar) noexcept {
  this->CopyMat(*this, NO_ADD, scalar);
  return *this;
}

/*
 * See header.
 */
Matrix Matrix::operator/(const float scalar) const noexcept(false){
  if (scalar == 0){
    throw MatrixException(DIV_BY_ZERO);
  }
  return Matrix(*this, NO_ADD,  ((float) 1)/scalar);
}

/*
 * See header.
 */
Matrix &Matrix::operator/=(const float scalar) noexcept(false){
  if (scalar == 0){
    throw MatrixException(DIV_BY_ZERO);
  }
  return (*this)*=(((float) 1)/scalar);
}

/*
 * See header.
 */
Matrix Matrix::operator+(const Matrix &rhs) const noexcept(false){
  if (!MatricesAreInSameDimensions(*this, rhs)){
    throw MatrixException(DIMENSION_INVALID);
  }
  Matrix m = Matrix(*this);
  return (m+=rhs);
}

/*
 * See header.
 */
Matrix &Matrix::operator+=(const Matrix &rhs) noexcept(false){
  if (!MatricesAreInSameDimensions(*this, rhs)){
    throw MatrixException(DIMENSION_INVALID);
  }
  for (int i=0; i<this->size_; ++i){
    this->mat_[i] += rhs.mat_[i];
  }
  return *this;
}

/*
 * See header.
 */
Matrix &Matrix::operator+=(const float scalar) noexcept{
  this->CopyMat(*this, scalar, NO_MULTIPLICATION);
  return *this;
}

/*
 * See header.
 */
float &Matrix::operator()(int row, int col) noexcept(false){
  if (row < 0 || col < 0 || row >= rows_ || col >= cols_) {
    throw MatrixException(RANGE_INVALID);
  }
  return this->mat_[row*cols_+col];
}

/*
 * See header.
 */
float Matrix::operator()(int row, int col) const noexcept(false){
  if (row < 0 || col < 0 || row >= rows_ || col >= cols_) {
    throw MatrixException(RANGE_INVALID);
  }
  return this->mat_[row*cols_+col];
}

/*
 * See header.
 */
float &Matrix::operator[](int idx) noexcept(false){
  if (idx < 0 || idx >= size_) {
    throw MatrixException(RANGE_INVALID);
  }
  return this->mat_[idx];
}

/*
 * See header.
 */
  float Matrix::operator[](int idx) const noexcept(false){
     if (idx < 0 || idx >= size_) {
    throw MatrixException(RANGE_INVALID);
  }
  return this->mat_[idx];
  }

/*
 * See header.
 */
std::istream &operator>>(std::istream &is, Matrix &matrix) noexcept(false){
  for(int i=0; i<matrix.size_; ++i) {
    if (is.bad() || is.fail()) {
      throw MatrixException(STREAM_INVALID);
    }
    is >> matrix[i];
  }
  return is;
}

/*
 * See header.
 */
std::ostream &operator<<(std::ostream &os, const Matrix &matrix) noexcept{
    for (int i=0; i<matrix.rows_; ++i){
      for (int j=0; j<matrix.cols_; ++j){
        os << matrix.mat_[i*matrix.cols_ + j];
        if (j < matrix.cols_-1){
          os << SEP_CELL;
        }
      }

      if (i<matrix.rows_-1){
        os << "\n";
      }
    }
  return os;
}

