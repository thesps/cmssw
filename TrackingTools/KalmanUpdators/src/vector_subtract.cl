__kernel void vector_subtract(	__global const float *r, 
				__global const float *rMeas, 
                        	__global float *r_out, 
                         )
{
    // get index of the work item
    //int index = get_global_id(0);

    // add the vector elements
	for (int i=0; i<2; i++) {
	    r_out[i] = r[i] + rMeas[i];
	}
}

