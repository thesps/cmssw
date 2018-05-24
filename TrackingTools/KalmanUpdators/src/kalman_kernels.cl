__kernel void vector_subtract(	__global const float *r,
				__global const float *rMeas,
                        	__global float *restrict r_out
                         )
{
// equivalent to
// r -= rMeas;
    // get index of the work item
    //int index = get_global_id(0);

    // add the vector elements
	for (int i = 0; i < 2; i++) {
	    r_out[i] = r[i] - rMeas[i];
	}
}

__kernel void matrix_add(	__global const float *V, 
				__global const float *VMeas, 
                        	__global float *restrict R 
                         )
{
// equivalent to 
// SMatDD R = V + VMeas;
	for (int i = 0; i < 2; i++) {
		for (int j = 0; j<2; j++) {
		    R[i*2+j] = V[i*2+j] + VMeas[i*2+j];
		}
	}
}


__kernel void matrix_invert(	__global const float *R, 
                        	__global float *restrict Rinv
                         )
{
// equivalent to
// bool ok = invertPosDefMatrix(R);
	float c = 1/(R[0]*R[3]-R[1]*R[2]);
	Rinv[0] = c * R[3];
	Rinv[1] = c * (-R[1]);
	Rinv[2] = c * (-R[2]);
	Rinv[3] = c * R[0];

    // get index of the work item
    //int index = get_global_id(0);
}

__kernel void matrix_project(	__global const float *C,
				__global const float *R,
				__global float *restrict K
			)
{
// equivalent to
// Mat5D K = C*pf.project(R);
// it is assumed that matrix element A_ij is stored as A[i*n+j] in 1D array
	K[0] = C[3]*R[0]; //K_11
	K[0] = C[8]*R[0]; //K_12
	K[0] = C[13]*R[0]; //K_21
	K[0] = C[18]*R[0]; //K_22
	K[0] = C[23]*R[0]; //K_31
	K[0] = C[4]*R[2]; //K_32
	K[0] = C[9]*R[2]; //K_41
	K[0] = C[14]*R[2]; //K_42
	K[0] = C[19]*R[2]; //K_51
	K[0] = C[24]*R[2]; //K_52
}

__kernel void matrix_projectsubtract(	__global const float *K,
					__global float *restrict K_out
			)
{
// equivalent to
// pf.projectAndSubtractFrom(M,K);
// it is assumed that matrix element A_ij is stored as A[i*n+j] in 1D array
	K[0] = 1;	
	K[1] = 0;	
	K[2] = 0;	
	K[3] = -K[0];	
	K[4] = -K[1];	
	K[5] = 0;	
	K[6] = 1;	
	K[7] = 0;	
	K[8] = -K[2];	
	K[9] = -K[3];	
	K[10] = 0;	
	K[11] = 0;	
	K[12] = 1;	
	K[13] = -K[4];	
	K[14] = -K[5];	
	K[15] = 0;	
	K[16] = 0;	
	K[17] = 0;	
	K[18] = 1-K[6];	
	K[19] = -K[7];	
	K[20] = 0;	
	K[21] = 0;	
	K[22] = 0;	
	K[23] = -K[8];	
	K[24] = 1-K[9];	
}


__kernel void vector_result(	__global const float *x,
				__global const float *K,
				__global const float *r,
				__global float *restrict fsv
			)
{
// equivalent to
// AlgebraicVector5 fsv = x + K * r; 
	for (int i; i <= 5; i++) {
		fsv[i] = x[i] + K[i*2]*r[0] + K[i*2+1]*r[1];
	}
}

__kernel void matrix_result(	__global const float *M,
				__global const float *C,
				__global const float *K,
				__global const float *V,
				__global const float *restrict fse,
			)
{
// equivalent to
// AlgebraicSymMatrix55 fse = ROOT::Math::Similarity(M, C) + ROOT::Math::Similarity(K, V);
	float MC[5*5];
	float MCMt[5*5];
	float KV[5*5];
	float KVKt[5*5];
	for (int i; i < 5; i++) {
		for (int j; j < 5; j++) {
			float acc = 0.0f;
			for (int k; k < 5; j++) {
				acc += M[k*5+i]+C[k*5+j];
			}
			MC[i*5+j] = acc;
		}
	}
	for (int i; i < 5; i++) {
		for (int j; j < 5; j++) {
			float acc = 0.0f;
			for (int k; k < 5; j++) {
				acc += MC[k*5+i]+M[j*5+k];
			}
			MCMt[i*5+j] = acc;
		}
	}
	for (int i; i < 5; i++) {
		for (int j; j < 5; j++) {
			float acc = 0.0f;
			for (int k; k < 5; j++) {
				acc += K[k*5+i]+V[k*5+j];
			}
			KV[i*5+j] = acc;
		}
	}
	for (int i; i < 5; i++) {
		for (int j; j < 5; j++) {
			float acc = 0.0f;
			for (int k; k < 5; j++) {
				acc += KV[k*5+i]+K[j*5+k];
			}
			KVKt[i*5+j] = acc;
		}
	}
	for (int i; i < 5*5; i++) {
		fse[i] = MCMt[i] + KVKt[i];
	}
}
