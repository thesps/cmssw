__kernel void vector_subtract(	__global const float *r, 
				__global const float *rMeas, 
                        	__global float *r_out
                         )
{
    // get index of the work item
    //int index = get_global_id(0);

    // add the vector elements
	for (int i=0; i<2; i++) {
	    r_out[i] = r[i] + rMeas[i];
	}
}

__kernel void matrix_add(	__global const float *V, 
				__global const float *VMeas, 
                        	__global float *R 
                         )
{
    // get index of the work item
    //int index = get_global_id(0);

    // add the vector elements
	for (int i=0; i<2; i++) {
		for (int j=0; j<2; j++) {
		    R[i*2+j] = V[i*2+j] + VMeas[i*2+j];
		}
	}
}


__kernel void matrix_invert(	__global const float *R, 
                        	__global float *Rinv
                         )
{
	float c = 1/(R[0]*R[3]-R[1]*R[2]);
	Rinv[0] = c * R[3];
	Rinv[1] = c * (-R[1]);
	Rinv[2] = c * (-R[2]);
	Rinv[3] = c * R[0];

    // get index of the work item
    //int index = get_global_id(0);
}
