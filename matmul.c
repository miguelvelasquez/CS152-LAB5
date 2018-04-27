// See LICENSE for license details.

#include "dataset.h"
#include "util.h"
#include <stddef.h>

void matmul(const size_t coreid, const size_t ncores, const size_t lda,  const data_t A[], const data_t B[], data_t C[])
{
  size_t i, j, k;

  for (i = 0; i < lda; i++) {
    for (j = coreid; j < lda; j += ncores) {
      data_t sum = 0;
      for (k = 0; k < lda; k++)
        sum += A[j*lda + k] * B[k*lda + i];
      C[i + j*lda] = sum;
    }
  }
}

void matmul_opt(const size_t coreid, const size_t ncores, const size_t lda,  const data_t A[], const data_t B[], data_t C[])
{
	size_t i, j, k;

  for (i = 0; i < lda; i++) {
    for (j = coreid * 4; j < lda; j = j + 8) {
      data_t sum = 0;
      for (k = 0; k < lda; k++)
        sum += A[j*lda + k] * B[k*lda + i];
        sum += A[(j+1)*lda + k] * B[k*lda + i];
        sum += A[(j+2)*lda + k] * B[k*lda + i];
        sum += A[(j+3)*lda + k] * B[k*lda + i];
      C[i + j*lda] = sum;
      C[i + (j+1)*lda] = sum;
      C[i + (j+2)*lda] = sum;
      C[i + (j+3)*lda] = sum;
    }
  }
}
