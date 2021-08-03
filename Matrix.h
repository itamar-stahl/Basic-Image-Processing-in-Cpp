#ifndef EX5__MATRIX_H_
#define EX5__MATRIX_H_
#include <iostream>

#define DEFAULT_VAL 0
#define DEFAULT_COLROW_SIZE 1
#define SEP_CELL " "
#define NO_ADD 0
#define NO_MULTIPLICATION 1

class Matrix {
 private:
  //****************** members *********************//
  int rows_;
  int cols_;
  int size_;
  float *mat_;

  //**************** privat constructors ******************//

  /*
   * Private copy ctor. See CopyMat doc.
   */
  Matrix(const Matrix &m, const float add_to_each, const float multiply_by) :
      rows_(m.rows_), cols_(m.cols_), size_(m.size_),
      mat_(new float [m.size_]) {CopyMat(m, add_to_each, multiply_by);}
  //****************** private methods ********************///

  /*
   * initialize all cells to zeros.
   */
  void Init() noexcept;

  /**
   * Copy the table of the given matrix m. For each cell: add: add_to_each,
   * multiply by: multiply_by
   * @param m
   * @param add_to_each
   * @param multiply_by
   */
  void CopyMat(const Matrix &m, const float add_to_each,
                              const float multiply_by) noexcept;

 public:

  //**************** public constructors ******************//

  /**
   * Return new matrix in size of rows*cols; initialize all cells to zeros.
   * @param rows
   * @param cols
   */
  Matrix(int rows, int cols) : rows_(rows),
      cols_(cols), size_(cols*rows), mat_(new float [cols*rows]) {Init();}

  /**
 * Default ctor. return new matrix in size 1*1; initialize all cells to zeros.
 * @param rows
 * @param cols
 */
  Matrix() : Matrix (DEFAULT_COLROW_SIZE, DEFAULT_COLROW_SIZE) {}

  /**
* Copy ctor.
* @param rows
* @param cols
*/
  Matrix(const Matrix &m) : Matrix(m, NO_ADD, NO_MULTIPLICATION) {}

  //****************** public methods ********************///

  /**
   * Transforms a matrix into a column vector. Supports function calling
   */
  Matrix &Vectorize() noexcept;

  /**
   * Prints matrix elements, no return value (void).
   * Prints space after each element (not including the last element in the row)
   * Prints a new line after each row (not including the last row).
   */
  void Print() const noexcept;


  //****************** operators *******************//

  /**
   * Matrix assignment
   * @param rhs
   * @return
   */
  Matrix &operator=(const Matrix &rhs) noexcept(false);

  /**
   * Matrix multiplication
   * @param rhs
   * @return
   */
  Matrix operator*(const Matrix &rhs) const noexcept(false);

  /**
   * multiply by a scalar (from the right)
   * @param scalar
   * @return
   */
  Matrix operator*(float const scalar) const noexcept;

  /**
   * multiply by a scalar (from the left)
   * @param scalar
   * @param m - matrix
   * @return
   */
  friend Matrix operator*(const float scalar, const Matrix &m) noexcept;

  /**
   * Matrix multiplied by rhs matrix
   * @param rhs
   * @return
   */
  Matrix & operator*=(const Matrix &rhs) noexcept(false);

  /**
   * multiply by a scalar from the right
   * @param scalar
   * @return
   */
  Matrix & operator*=(const float scalar) noexcept;

  /**
   * Division matrix by a scalar - from the right
   * @param scalar
   * @return
   */
  Matrix operator/(const float scalar) const noexcept(false);

  /**
   *Division Matrix itself by a scalar
   * @param scalar
   * @return
   */
  Matrix &operator/=(const float scalar) noexcept(false);

  /**
   * Matrix Additon
   * @param rhs
   * @return
   */
  Matrix operator+(const Matrix &rhs) const noexcept(false);

  /**
   * Matrix addition on itself - rhs matrix
   * @param rhs
   * @return
   */
  Matrix &operator+=(const Matrix &rhs) noexcept(false);

  /**
   * Add scalar to Matrix
   * @param scalar
   * @return
   */
  Matrix &operator+=(const float scalar) noexcept;

  /**
   * Non-const access by row, col indexes
   * @param row
   * @param col
   * @return
   */
  float &operator()(int row, int col) noexcept(false);

  /**
   * Const access by row, col indexes
   * @param row
   * @param col
   * @return
   */
  float operator()(int row, int col) const noexcept(false);

  /**
   * Non-const access by a single index
   * @param row
   * @param col
   * @return
   */
  float &operator[](int idx) noexcept(false);

   /**
   * Const access by a single index
   * @param row
   * @param col
   * @return
   */
  float operator[](int idx) const noexcept(false);

  /**
   * Equality
   * @param rhs
   * @return
   */
  bool operator==(const Matrix &rhs) const noexcept;

  /**
   * Non-Equality
   * @param rhs
   * @return
   */
  bool operator!=(const Matrix &rhs) const noexcept;
  /**
   * Input stream
   * @param is
   * @param matrix
   * @return
   */
  friend std::istream &operator>>(std::istream &is,
                                  Matrix &matrix) noexcept(false);

  /**
   * Output Stream
   * @param os
   * @param matrix
   * @return
   */
  friend std::ostream &operator<<(std::ostream &os,
                                  const Matrix &matrix) noexcept;

  //****************** getters & setters *******************//
  int GetRows() const { return this->rows_; }
  int GetCols() const { return this->cols_; }

  //************************************************//
  //**************** destructor ******************//
  //************************************************//
  ~Matrix() {
    delete[] mat_;
  }
};


#endif //EX5__MATRIX_H_
