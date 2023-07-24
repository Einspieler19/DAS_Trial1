
#include <cstdio>
#include <fstream>
#include "utils.hpp"
#include "My_MatInverse_Top.hpp"



int main(){
	
	int pass_fail=0;
	
	MInv_T A[RowsA][ColsA];
 	MInv_T B[RowsB][ColsB];
 	MInv_TOut DUT_C[RowsA][ColsB];
	MInv_T REF_C[RowsA][ColsB];
		
		
	MInv_T_RD  Ar[NUM_Tests][RowsA][ColsA];
	MInv_T_RD  Br[NUM_Tests][RowsB][ColsB];
	MInv_T_RD  REF_Cr[NUM_Tests][RowsA][ColsB];

	
	// 文件size, 读取的参数
	int A_size = RowsA * ColsA * NUM_Tests;
	int B_size = RowsB * ColsB * NUM_Tests;
	int C_size = RowsA * ColsB * NUM_Tests;
		
	MInv_T_RD* A_ptr = reinterpret_cast<MInv_T_RD*>(Ar);
	MInv_T_RD* B_ptr = reinterpret_cast<MInv_T_RD*>(Br);
	MInv_T_RD* C_ptr = reinterpret_cast<MInv_T_RD*>(REF_Cr);
	
	
	std::string file_A = "Matrix_A.txt";
	std::string file_B = "Matrix_B.txt";
	std::string file_C = "MInv_gold.txt";
		
	readTxt(file_A, A_ptr, A_size);
//	readTxt(file_B, B_ptr, B_size);
	readTxt(file_C, C_ptr, C_size);

		
	hls::stream<MInv_T> matrixInStrm;
	hls::stream<MInv_T> matrixDUTOutStrm;
//	hls::stream<MInv_T> matrixREF_CStrm;
//	hls::stream<MInv_TOut> matrixDUT_CStrm;
		
	for (int i = 0; i < NUM_Tests; i++) {

		// 写入流
	for (int r = 0; r < RowsA; r++) {
		for (int c = 0; c < ColsA; c++) {
			A[r][c] = MInv_T(Ar[i][r][c]);
			matrixInStrm.write(A[r][c]);
		}
	}


		
	My_MatInverse_Top(matrixInStrm,matrixDUTOutStrm);

	
		// 读出流
	for (int r = 0; r < RowsA; r++) {
		for (int c = 0; c < ColsB; c++) {
		matrixDUTOutStrm.read(DUT_C[r][c]);
		}	
	}
	
	
	std::cout<<"DUT:"<<std::endl;	
	
	for (int r = 0; r < RowsA; r++) {
		for (int c = 0; c < ColsB; c++) {
		std::cout<<DUT_C[r][c]<<"     ";
		}	
		std::cout<<std::endl;	
	}
	
		
	
	std::cout<<"REF:"<<std::endl;	
	
	for (int r = 0; r < RowsA; r++) {
		for (int c = 0; c < ColsB; c++) {
		std::cout<<REF_Cr[i][r][c]<<"     ";
		}	
		std::cout<<std::endl;	
	}
	
			
	MInv_T_RD0 Error_Sum = 0;
	
	std::cout<<"error:"<<std::endl;	
	
	for (int r = 0; r < RowsA; r++) {
		for (int c = 0; c < ColsB; c++) {
		std::cout<<MInv_T_RD(DUT_C[r][c])- REF_Cr[i][r][c]<<" ";
		Error_Sum+=csquf(MInv_T_RD(DUT_C[r][c])- REF_Cr[i][r][c]);
		}
		std::cout<<std::endl;
	}
	
	std::cout<<Error_Sum<<std::endl;
	if (Error_Sum > 1)
	pass_fail=1;
	
	}


	
	
	if (pass_fail == 1) {
	std::cout << "TB:Fail" << std::endl;
} else {
	std::cout << "TB:Pass" << std::endl;
}
std::cout << "" << std::endl;
return (pass_fail);
	

}
