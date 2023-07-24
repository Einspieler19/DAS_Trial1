
#include "My_MatInverse_Top.hpp"

// top-level of HLS design


void My_MatInverse_Top(hls::stream<MInv_T>& matrixInStrm,
					hls::stream<MInv_TOut>& matrixOutStrm){
						
//	My_MatInverse_Baseline(matrixInStrm,matrixOutStrm);
//	Vitis_chol_Basic(matrixInStrm,matrixOutStrm);
//	Vitis_chol_Alt1(matrixInStrm,matrixOutStrm);
	Vitis_chol_Alt2(matrixInStrm,matrixOutStrm);
//	My_LDL(matrixInStrm,matrixOutStrm);

}
