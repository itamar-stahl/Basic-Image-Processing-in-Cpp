#include "Filters.h"
#include "Matrix.h"
#include <cmath>

#define MAX_VALUE 255
#define MIN_VALUE 0
#define BLUR_MATRIX_VALS {1, 2, 1, 2, 4, 2, 1, 2, 1}
#define CONVO_MATRIX_DIM {3, 3}
#define CONVO_MATRIX_CENTER {1, 1}
#define CONVO_MATRIX_SIZE 9
#define BLUR_CONST (1.0/16)
#define SOBEL_CONST (1.0/8)
#define SOBEL_MATRIX_X_VALS {1, 0, -1, 2, 0, -2, 1, 0, -1}
#define SOBEL_MATRIX_Y_VALS {1, 2, 1, 0, 0, 0, -1, -2, -1}

/**
 * Performs quantization on the input image by the given number of levels.
 * Returns new matrix which is the result of running the operator on the image.
 * @param image
 * @param levels
 * @return
 */

/*
 * A helper function which is used after convolution.
 */
inline void PutCellInRange(float &val){
  if (val > MAX_VALUE) {val = MAX_VALUE;}
  if (val < MIN_VALUE) {val = MIN_VALUE;}
}

/*
 * A helper function which is used after convolution. 
 */
void PutMatrixInRange(Matrix &m, int m_size){
  for (int i=0; i<m_size; ++i){
    PutCellInRange(m[i]);
  }
}

/*
 * A helper function which is used during Quantization. 
 */
inline int FlooredMean(int a, int b){
  return (int) floor((a+b)/2);
}

/*
 * A helper function which is used during Quantization.
 * Creates a table with conversion value for any value in the range 0-255.
 */
void InstallTable(int levels, float table[MAX_VALUE+1]){
  int base_interval = (MAX_VALUE+1)/levels;
  int last_break = 0;
  int curr_interval_average = FlooredMean(last_break,
                                          last_break+base_interval-1);
  for (int i=0;i <= MAX_VALUE; ++i){
      if (i==last_break+base_interval){
        last_break += base_interval;
        curr_interval_average = FlooredMean(last_break,
                                            last_break+base_interval-1);
      }
      table[i] = curr_interval_average;
  }
}

/*
 * A helper function which is used during convolution.
 */
Matrix GetConvolutionMatrix(float constant, float const mat_vals[]){
  int dimensions[] = CONVO_MATRIX_DIM;
  Matrix convolution_mat(dimensions[0], dimensions[1]);
  for (int i=0; i < CONVO_MATRIX_SIZE; ++i){
    convolution_mat[i] = (constant)*mat_vals[i];
  }
  return convolution_mat;
}

/*
 * A helper function which is used during convolution as part of the Sobel filter.
 */
inline Matrix GetSobelXMatrix(){
  float matrix_vals[CONVO_MATRIX_SIZE] = SOBEL_MATRIX_X_VALS;
  return GetConvolutionMatrix(SOBEL_CONST, matrix_vals);
}

/*
 * A helper function which is used during convolution as part of the Sobel filter.
 */
inline Matrix GetSobelYMatrix(){
  float matrix_vals[CONVO_MATRIX_SIZE] = SOBEL_MATRIX_Y_VALS;
  return GetConvolutionMatrix(SOBEL_CONST, matrix_vals);
}

/*
 * A helper function which is used during convolution as part of the blur filter.
 *
 */
inline Matrix GetBlurMatrix(){
  float matrix_vals[CONVO_MATRIX_SIZE] = BLUR_MATRIX_VALS;
  return GetConvolutionMatrix(BLUR_CONST, matrix_vals);
}

/*
 * A helper function which is used during convolution.
 * Allows accessing to a matrix's values, return 0 when indexes are out of
 * the matrix's range.
 */
float PaddedMatrixAccessor(const Matrix &m, int row, int col){
  if (row < 0 || col < 0 || row >= m.GetRows() || col >= m.GetCols()){
    return 0;
  }
  return m(row, col);
}

/*
 * A helper function which is used during convolution.
 */
float CellConvolution(const Matrix &image, const Matrix &convo_mat, int center[],
                      int row, int col){

  float result = 0.0;
  for (int i =0; i < convo_mat.GetRows(); ++i){
    for(int j=0; j < convo_mat.GetCols(); ++j){
      result += convo_mat(i,j)*PaddedMatrixAccessor(image,
                                                    row + (i-center[0]),
                                                    col + (j-center[1]));
    }
  }
  return rintf(result);
}

/*
 * A helper function to Sobel and blur filters.
 * Performs a convolution of two matrices and puts the result into &output.
 */
void Convolution(Matrix &output, const Matrix &image, const Matrix &convo_mat){
  int center[] = CONVO_MATRIX_CENTER;
  for (int i=0; i < image.GetRows(); ++i){
    for(int j=0; j <image.GetCols(); ++j){
      output(i, j) = CellConvolution(image, convo_mat, center, i, j);
    }
  }
}

/**
 * Performs quantization on the input image by the given number of levels.
 * Returns new matrix which is the result of running the operator on the image.
 * @param image
 * @param levels
 * @return
 */
Matrix Quantization(const Matrix& image,int levels){
float conversion_table[MAX_VALUE+1];
InstallTable(levels, conversion_table);
int image_size = image.GetCols()*image.GetRows();
Matrix quantized_img(image);
for (int i=0; i<image_size; ++i){
  quantized_img[i] = conversion_table[(int) quantized_img[i]];
}
  return quantized_img;
}

/**
 * Performs Gaussian blurring on the input image.
 * Returns new matrix which is the result of running the operator on the image.
 * @param image
 * @return
 */
Matrix Blur(const Matrix& image){
  Matrix blurred_img(image.GetRows(), image.GetCols());
  Matrix blur_m = GetBlurMatrix();
  Convolution(blurred_img, image, blur_m);
  PutMatrixInRange(blurred_img,
                   blurred_img.GetRows()*blurred_img.GetCols());
  return blurred_img;
}

/**
 * Performs Sobel edge detection on the input image.
 * Returns new matrix which is the result of running the operator on the image.
 * @param image
 * @return
 */
Matrix Sobel(const Matrix& image){
  int rows = image.GetRows(); int  cols = image.GetCols();
  Matrix m_x(rows, cols);
  Matrix m_y(rows, cols);
  Matrix convo_x = GetSobelXMatrix();
  Matrix convo_y = GetSobelYMatrix();
  Convolution(m_x, image, convo_x);
  Convolution(m_y, image, convo_y);
  m_x += m_y;
  PutMatrixInRange(m_x, rows*cols);
  return m_x;
}


