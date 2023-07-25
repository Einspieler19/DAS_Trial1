
#include <cmath>
#include <hls_math.h>
#include "My_MatInverse_Top.hpp"

// Baseline alg.
void My_MatInverse_Baseline(hls::stream<MInv_T>& matrixInStrm,
		hls::stream<MInv_TOut>& matrixOutStrm){

	//	MInv_T temp=0;
	MInv_TOut Sum[ColsA];
	MInv_T A[RowsA][ColsA];
	MInv_TOut L[RowsA][ColsB];


	// 写入流
	loop_strmAIn:
	for (int r = 0; r < RowsA; r++) {
		for (int c = 0; c < ColsA; c++) {
			matrixInStrm.read(A[r][c]);
		}
	}



	for (int c = 0; c < ColsA; ++c) {
		for (int r = c; r < ColsB; ++r) {
			Sum[c] = 0;
			k_loop:
			for (int k = 0; k < c; ++k) {
				Sum[c] += L[r][k] * myconj(L[c][k]);
			}
			if (r == c) {
				mysqrt_real((A[r][c] - Sum[c]), L[r][c]);
			} else {
				L[r][c] = (A[r][c] - Sum[c]) / L[c][c];
			}
		}
	}

	// 写入流
	loop_strmOut:
	for (int r = 0; r < RowsA; r++) {
		for (int c = 0; c < ColsB; c++) {
			matrixOutStrm.write(L[r][c]);
		}
	}

}




void Vitis_chol_Basic(hls::stream<MInv_T>& matrixInStrm,
		hls::stream<MInv_TOut>& matrixOutStrm){

	MInv_T  L_internal[RowsColsA][RowsColsA];
	MInv_T  retrieved_L;
	//    MInv_T sum[RowsColsA];


	PROD_T prod;
	ACCUM_T sum[RowsColsA];
	ACCUM_T A_cast_to_sum;    // A with the same dimensions as sum.
	ACCUM_T prod_cast_to_sum; // prod with the same dimensions as sum.

	ADD_T A_minus_sum;
	DIAG_T new_L_diag;         // sqrt(A_minus_sum)
	OFF_DIAG_T new_L_off_diag; // sum/L
	OFF_DIAG_T L_cast_to_new_L_off_diag;

	MInv_TOut new_L;


	//	MInv_T temp=0;
	MInv_TOut Sum[ColsA];
	MInv_T A[RowsA][ColsA];
	MInv_TOut L[RowsA][ColsB];

	// 写入流
	loop_strmAIn:
	for (int r = 0; r < RowsA; r++) {
		for (int c = 0; c < ColsA; c++) {
			matrixInStrm.read(A[r][c]);
		}
	}

	for (int r = 0; r < RowsA; r++) {
		for (int c = 0; c < ColsA; c++) {
			if (r==c){
				L[r][c]=1;
			}
			else L[r][c]=0;
		}
	}

	col_loop:
	for (int j = 0; j < RowsColsA; j++) {
		sum[j] = 0;

		// Calculate the diagonal value for this column
		diag_loop:
		for (int k = 0; k < RowsColsA; k++) {
			if (k <= (j - 1)) {
				retrieved_L = L[j][k];
				sum[j] += hls::x_conj(retrieved_L) * retrieved_L;
			}
		}
		A_cast_to_sum = A[j][j];
		A_minus_sum = A_cast_to_sum - sum[j];

		mysqrt_real(A_minus_sum, L[j][j]);

		// Calculate the off diagonal values for this column
		off_diag_loop:
		for (int i = 0; i < RowsColsA; i++) {
			if (i > j) {
				sum[j] = A[i][j];
				sum_loop:
				for (int k = 0; k < RowsColsA; k++) {
#pragma HLS PIPELINE II=1
					if (k <= (j - 1)) {
						prod = -L[i][k] * hls::x_conj(L[j][k]);
						prod_cast_to_sum = prod;
						sum[j] += prod_cast_to_sum;
					}
				}
				new_L_off_diag = sum[j];

				L_cast_to_new_L_off_diag = L[j][j];

				// Diagonal is always real, avoid complex division
				new_L_off_diag = new_L_off_diag / hls::x_real(L_cast_to_new_L_off_diag);

				// Round to target format using method specifed by traits defined types.
				new_L = new_L_off_diag;

				L[i][j] = new_L;
				L_internal[i][j] = new_L;
			} else if (i < j) {
				L[i][j] = 0;
			}
		}
	}

	// 写入流
	loop_strmOut:
	for (int r = 0; r < RowsA; r++) {
		for (int c = 0; c < ColsB; c++) {
			matrixOutStrm.write(L[r][c]);
		}
	}


}









void Vitis_chol_Alt1(hls::stream<MInv_T>& matrixInStrm,
		hls::stream<MInv_TOut>& matrixOutStrm){

	MInv_T  L_internal[(RowsColsA * RowsColsA - RowsColsA) / 2];
	MInv_T  retrieved_L;
	//    MInv_T sum[RowsColsA];


	MInv_T square_sum;
	MInv_T product_sum;
	PROD_T prod;
	ACCUM_T sum[RowsColsA];
	ACCUM_T A_cast_to_sum;    // A with the same dimensions as sum.
	ACCUM_T prod_cast_to_sum; // prod with the same dimensions as sum.

	ADD_T A_minus_sum;
	MInv_T new_L_diag;         // sqrt(A_minus_sum)
	OFF_DIAG_T new_L_off_diag; // sum/L
	OFF_DIAG_T L_cast_to_new_L_off_diag;

	MInv_TOut new_L;

	MInv_T0 diag_internal[RowsColsA];
	MInv_T0 L_diag_recip;
	MInv_T0 new_L_diag_recip;
	MInv_T A_minus_sum_cast_diag;
	MInv_T prod_cast_to_off_diag;


	//	MInv_T temp=0;
	MInv_TOut Sum[ColsA];
	MInv_T A[RowsA][ColsA];
	MInv_TOut L[RowsA][ColsB];

	// 写入流
	loop_strmAIn:
	for (int r = 0; r < RowsA; r++) {
		for (int c = 0; c < ColsA; c++) {
			matrixInStrm.read(A[r][c]);
		}
	}
	row_loop:
	for (int i = 0; i < RowsColsA; i++) {
		// Index generation for optimized/packed L_internal memory
		int i_sub1 = i - 1;
		int i_off = ((i_sub1 * i_sub1 - i_sub1) / 2) + i_sub1;

		// Off diagonal calculation
		square_sum = 0;
		col_loop:
		for (int j = 0; j < i; j++) {
			//#pragma HLS loop_tripcount max = 1 + RowsColsA / 2
			// Index generation
			int j_sub1 = j - 1;
			int j_off = ((j_sub1 * j_sub1 - j_sub1) / 2) + j_sub1;
			// Prime the off-diagonal sum with target elements A value.
			product_sum = A[i][j];

			sum_loop:
			for (int k = 0; k < j; k++) {
				//#pragma HLS loop_tripcount max = 1 + RowsColsA / 2
#pragma HLS PIPELINE II = 1
				prod = -L_internal[i_off + k] * hls::x_conj(L_internal[j_off + k]);
				prod_cast_to_sum = prod;
				product_sum += prod_cast_to_sum;
			}
			prod_cast_to_off_diag = product_sum;
			// Fetch diagonal value
			L_diag_recip = diag_internal[j];
			// Diagonal is stored in its reciprocal form so only need to multiply the product sum
			mymult_complex_real(prod_cast_to_off_diag, L_diag_recip, new_L_off_diag);


			// Round to target format using method specifed by traits defined types.
			new_L = new_L_off_diag;
			// Build sum for use in diagonal calculation for this row.
			square_sum += hls::x_conj(new_L) * new_L;
			// Store result
			L_internal[i_off + j] = new_L;

			L[i][j] = new_L; // store in lower triangle

			L[j][i] = 0;     // Zero upper
		}
		// Diagonal calculation
		A_cast_to_sum = A[i][i];
		A_minus_sum = A_cast_to_sum - square_sum;
		mysqrt_real(A_minus_sum, new_L_diag);

		// Round to target format using method specifed by traits defined types.
		new_L = new_L_diag;
		L[i][i] = new_L;

		// Generate the reciprocal of the diagonal for internal use to aviod the latency of a divide in every
		// off-diagonal calculation
		A_minus_sum_cast_diag = A_minus_sum;
		new_L_diag_recip = hls::sqrt(A_minus_sum_cast_diag.real());
		// Store diagonal value
		diag_internal[i] = 1/new_L_diag_recip;
	}

	// 写入流
	loop_strmOut:
	for (int r = 0; r < RowsA; r++) {
		for (int c = 0; c < ColsB; c++) {
			matrixOutStrm.write(L[r][c]);
		}


	}
}











void Vitis_chol_Alt2(hls::stream<MInv_T>& matrixInStrm,
		hls::stream<MInv_TOut>& matrixOutStrm){

	MInv_T  L_internal[RowsColsA][RowsColsA];
	MInv_T  retrieved_L;
	//    MInv_T sum[RowsColsA];


	MInv_T square_sum;
	MInv_T product_sum;
	PROD_T prod;
	ACCUM_T sum[RowsColsA];
	ACCUM_T A_cast_to_sum;    // A with the same dimensions as sum.
	ACCUM_T prod_cast_to_sum; // prod with the same dimensions as sum.

	ADD_T A_minus_sum;
	MInv_T new_L_diag;         // sqrt(A_minus_sum)
	OFF_DIAG_T new_L_off_diag; // sum/L
	OFF_DIAG_T L_cast_to_new_L_off_diag;

	MInv_T square_sum_array[RowsColsA];
	MInv_TOut new_L;

	MInv_T0 diag_internal[RowsColsA];
	MInv_T0 L_diag_recip;
	MInv_T prod_column_top;
	MInv_T product_sum_array[RowsColsA];
	MInv_T0 new_L_diag_recip;
	MInv_T A_minus_sum_cast_diag;
	MInv_T prod_cast_to_off_diag;

	int ARCH2_ZERO_LOOP = 1;

	//	MInv_T temp=0;
	MInv_TOut Sum[ColsA];
	MInv_T A[RowsA][ColsA];
	MInv_TOut L[RowsA][ColsB];

	// 写入流
	loop_strmAIn:
	for (int r = 0; r < RowsA; r++) {
		for (int c = 0; c < ColsA; c++) {
			matrixInStrm.read(A[r][c]);
		}
	}


	col_loop:
	for (int j = 0; j < RowsColsA; j++) {
		// Diagonal calculation
		A_cast_to_sum = A[j][j];
		if (j == 0) {
			A_minus_sum = A_cast_to_sum;
		} else {
			A_minus_sum = A_cast_to_sum - square_sum_array[j];
		}
		mysqrt_real(A_minus_sum, new_L_diag);
		// Round to target format using method specifed by traits defined types.
		new_L = new_L_diag;
		// Generate the reciprocal of the diagonal for internal use to aviod the latency of a divide in every
		// off-diagonal calculation
		//        A_minus_sum_cast_diag = A_minus_sum;
		//        cholesky_rsqrt(hls::x_real(A_minus_sum_cast_diag), new_L_diag_recip);
		//        // Store diagonal value
		L[j][j] = new_L;
		new_L_diag_recip = 1/new_L.real();
		sum_loop:
		for (int k = 0; k <= j; k++) {
			// Define average trip count for reporting, loop reduces in length for every iteration of col_loop
#pragma HLS loop_tripcount max = (1 + (RowsColsA/2))
			// Same value used in all calcs
			// o Implement -1* here
			prod_column_top = -hls::x_conj(L_internal[j][k]);

			// NOTE: Using a fixed loop length combined with a "if" to implement reducing loop length
			// o Ensures the inner loop can achieve the maximum II (1)
			// o May introduce a small overhead resolving the "if" statement but HLS struggled to schedule when the variable
			//   loop bound expression was used.
			// o Will report inaccurate trip count as it will reduce by one with the col_loop
			// o Variable loop bound code: row_loop: for(int i = j+1; i < RowsColsA; i++) {
			row_loop:
			for (int i = 0; i < RowsColsA; i++) {
				// IMPORTANT: row_loop must not merge with sum_loop as the merged loop becomes variable length and HLS will struggle
				// with scheduling
				//#pragma HLS LOOP_FLATTEN off
#pragma HLS PIPELINE II = 5
				//#pragma HLS UNROLL FACTOR = UNROLL_FACTOR

				if (i > j) {
					prod = L_internal[i][k] * prod_column_top;
					prod_cast_to_sum = prod;

					if (k == 0) {
						// Prime first sum
						A_cast_to_sum = A[i][j];
						product_sum = A_cast_to_sum;
					} else {
						product_sum = product_sum_array[i];
					}

					if (k < j) {
						// Accumulate row sum of columns
						product_sum_array[i] = product_sum + prod_cast_to_sum;
					} else {
						// Final calculation for off diagonal value
						prod_cast_to_off_diag = product_sum;
						// Diagonal is stored in its reciprocal form so only need to multiply the product sum

						mymult_complex_real(prod_cast_to_off_diag, new_L_diag_recip, new_L_off_diag);

						// Round to target format using method specifed by traits defined types.
						new_L = new_L_off_diag;
						// Build sum for use in diagonal calculation for this row.
						square_sum_array[i] += hls::x_conj(new_L) * new_L;
						// Store result
						L_internal[i][j] = new_L;
						// NOTE: Use the upper/lower triangle zeroing in the subsequent loop so the double memory access
						// does not
						// become a bottleneck
						// o Results in a further increase of DSP resources due to the higher II of this loop.
						// o Retaining the zeroing operation here can give this a loop a max II of 2 and HLS will
						// resource share.
						L[i][j] = new_L;                                  // Store in lower triangle
						if (!ARCH2_ZERO_LOOP) L[j][i] = 0; // Zero upper
					}
				}
			}
		}
	}
	// Zero upper/lower triangle
	// o Use separate loop to ensure main calcuation can achieve an II of 1
	// o As noted above this may increase the DSP resources.
	// o Required when unrolling the inner loop due to array dimension access
	if (ARCH2_ZERO_LOOP) {
		zero_rows_loop:
		for (int i = 0; i < RowsColsA - 1; i++) {
			zero_cols_loop:
			for (int j = i + 1; j < RowsColsA; j++) {
				// Define average trip count for reporting, loop reduces in length for every iteration of zero_rows_loop
#pragma HLS loop_tripcount max = (1 + (RowsColsA/2))
				//#pragma HLS PIPELINE
				L[i][j] = 0; // Zero upper
			}
		}
	}
	// 写入流
	loop_strmOut:
	for (int r = 0; r < RowsA; r++) {
		for (int c = 0; c < ColsB; c++) {
			matrixOutStrm.write(L[r][c]);
		}

	}
}




void My_LDL(hls::stream<MInv_T>& matrixInStrm,
		hls::stream<MInv_TOut>& matrixOutStrm){

	MInv_T  L_internal[RowsColsA][RowsColsA];
	MInv_T	D[RowsColsA];

	MInv_T product_sum;
	PROD_T prod;
	ACCUM_T sum[RowsColsA];
	MInv_T0 A_cast_to_sum;    // A with the same dimensions as sum.
	MInv_T0 A_minus_sum;


	ACCUM_T prod_cast_to_sum; // prod with the same dimensions as sum.


	MInv_T new_L_diag;         // sqrt(A_minus_sum)
	MInv_T0 new_L_d;
	OFF_DIAG_T new_L_off_diag; // sum/L
	OFF_DIAG_T L_cast_to_new_L_off_diag;

	MInv_T0 square_sum_array[RowsColsA];
	MInv_T0 square_sum;
	MInv_TOut new_L;

	MInv_T0 diag_internal[RowsColsA];
	MInv_T0 L_diag_recip;
	MInv_T prod_column_top;
	MInv_T product_sum_array[RowsColsA];
	MInv_T0 new_L_diag_recip;
	MInv_T A_minus_sum_cast_diag;
	MInv_T prod_cast_to_off_diag;

	int ARCH2_ZERO_LOOP = 1;

	//	MInv_T temp=0;
	MInv_TOut Sum[ColsA];
	MInv_T0 Diag[ColsA];
	MInv_T A[RowsA][ColsA];
	MInv_TOut L[RowsA][ColsB];







	// 写入流
	loop_strmAIn:
	for (int i = 0; i < RowsA; i++) {
		for (int j = 0; j < ColsA; j++) {
			matrixInStrm.read(A[i][j]);
			L_internal[i][j]=0;
		}
		square_sum_array[i]=0;
		product_sum_array[i]=0;
	}



	col_loop:
	for (int j = 0; j < RowsColsA; j++) {
		// 取到A对角元
		A_cast_to_sum = A[j][j].real();

		// 完成减法
		if (j == 0) {
			A_minus_sum = A_cast_to_sum;
		} else {
			A_minus_sum = A_cast_to_sum - square_sum_array[j];
		}
		// 顺便赋L对角：考虑放在哪里更好
		L[j][j] = 1;

		new_L_d = A_minus_sum;

		//存对角元素
		Diag[j] = new_L_d;

		// 存逆元素
		new_L_diag_recip = 1/new_L_d;


		sum_loop:
		for (int k = 0; k <= j; k++) {
#pragma HLS loop_tripcount max = (1 + (RowsColsA/2))
//#pragma HLS DEPENDENCE dependent=false type=intra variable=product_sum_array
//#pragma HLS DEPENDENCE dependent=false type=intra variable=square_sum_array

			//预存 i * d
//			mymult_complex_real(-hls::x_conj(L_internal[j][k]),Diag[k],prod_column_top);
			prod_column_top = -hls::x_conj(L_internal[j][k]);
			row_loop:
			for (int i = 0; i < RowsColsA; i++) {
				// 不要和sum loop合并，调度困难？
				#pragma HLS LOOP_FLATTEN off

				#pragma HLS PIPELINE II = 1
//				#pragma HLS UNROLL FACTOR = 1

				if (i > j) {
					prod = L_internal[i][k] * prod_column_top;
					prod_cast_to_sum = prod;

					// 遍历起始: 找到A元
					if (k == 0) {
						// Prime first sum

						product_sum = A[i][j];
					}

					// 遍历中途: 累加
					else {
						product_sum = product_sum_array[i];
					}

					if (k < j) {
						// 遍历中途: 累加
						product_sum_array[i] = product_sum + prod_cast_to_sum;

					}
					// 遍历结束: 累加
					else {

						// Final calculation for off diagonal value
						prod_cast_to_off_diag = product_sum;

						// Diagonal is stored in its reciprocal form so only need to multiply the product sum

						mymult_complex_real(prod_cast_to_off_diag, new_L_diag_recip, new_L_off_diag);

						// 定点数round
						new_L = new_L_off_diag;


						// 存储 i*d*i
//						mymult_conj(new_L, square_sum);
						square_sum = square_sum* Diag[k];
//
						square_sum_array[i] += square_sum;
						// off元素最终赋值
						L_internal[i][j] = new_L;
						// 另起置0循环
						// o 提高II，DSP要求多
						// o 置0留在这里的话最多II=2

						L[i][j] = new_L;
						if (!ARCH2_ZERO_LOOP) L[j][i] = 0; // Zero upper
					}
				}
			}
		}
	}

	// Zero upper/lower triangle
	// o 单拎达到 II=1
	// o DSP 增多
	// o Required when unrolling the inner loop due to array dimension access
	if (ARCH2_ZERO_LOOP) {
		zero_rows_loop:
		for (int i = 0; i < RowsColsA - 1; i++) {
			zero_cols_loop:
			for (int j = i + 1; j < RowsColsA; j++) {
				// 1+N/2 次循环
#pragma HLS loop_tripcount max = (1 + (RowsColsA/2))
				//#pragma HLS PIPELINE
				L[i][j] = 0; // Zero upper
			}
		}
	}

	// 写入流
	loop_strmOut:
	for (int r = 0; r < RowsA; r++) {
		for (int c = 0; c < ColsB; c++) {
			matrixOutStrm.write(L[r][c]);
		}

	}
}



