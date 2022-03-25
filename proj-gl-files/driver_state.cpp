#include "driver_state.h"
#include <cstring>

using namespace std;

float triangleArea(vec2 a, vec2 b, vec2 c);
void newVertex(driver_state& state, data_geometry* tri, const data_geometry* x, const data_geometry* y, int plane, bool pos, float* newd);

driver_state::driver_state()
{
}

driver_state::~driver_state()
{
    delete [] image_color;
    delete [] image_depth;
}

// This function should allocate and initialize the arrays that store color and
// depth.  This is not done during the constructor since the width and height
// are not known when this class is constructed.
void initialize_render(driver_state& state, int width, int height)
{
    state.image_width=width;
    state.image_height=height;
    int size = width * height;
    //cout<<"TODO: allocate and initialize state.image_color and state.image_depth."<<endl;
	state.image_color = new pixel[size];
	state.image_depth = new float[size];

	for (int i = 0; i < (size); i++) {
		state.image_color[i] = make_pixel(0, 0, 0);
        state.image_depth[i] = 2;
    }
}

// This function will be called to render the data that has been stored in this class.
// Valid values of type are:
//   render_type::triangle - Each group of three vertices corresponds to a triangle.
//   render_type::indexed -  Each group of three indices in index_data corresponds
//                           to a triangle.  These numbers are indices into vertex_data.
//   render_type::fan -      The vertices are to be interpreted as a triangle fan.
//   render_type::strip -    The vertices are to be interpreted as a triangle strip.
void render(driver_state& state, render_type type)
{
    //cout<<"TODO: implement rendering."<<endl;
    
    switch (type) {           
        case(render_type::triangle): {

            data_vertex* vertex_seg[state.num_vertices];
            for (int i = 0, j = 0; i < state.num_vertices* state.floats_per_vertex; i += state.floats_per_vertex, j++) {
                vertex_seg[j] = new data_vertex[MAX_FLOATS_PER_VERTEX];
                vertex_seg[j]->data = &state.vertex_data[i];
            }


            data_geometry* vertex_geo[state.num_vertices];
            for (int i = 0; i < state.num_vertices; i++) {
                vertex_geo[i] = new data_geometry();
                vertex_geo[i]->data = vertex_seg[i]->data;
                state.vertex_shader(*vertex_seg[i], *vertex_geo[i], state.uniform_data);
            }

            for (int i = 0; i < state.num_vertices / 3; i++) {
                const data_geometry** pair_vertices = const_cast<const data_geometry**>(vertex_geo + (3 * i));
                clip_triangle(state, *pair_vertices[0], *pair_vertices[1], *pair_vertices[2], 0);		
            }

            break;
        }
        case(render_type::indexed):{
		data_vertex* vertex_seg[state.num_vertices];
		for (int i = 0, j = 0; i < state.num_vertices* state.floats_per_vertex; i += state.floats_per_vertex, j++) {
			vertex_seg[j] = new data_vertex[MAX_FLOATS_PER_VERTEX];
			vertex_seg[j]->data = &state.vertex_data[i];
		}

		data_geometry* vertex_geo[state.num_triangles * 3];
		for (int i = 0; i < state.num_triangles * 3; i++) {
			vertex_geo[i] = new data_geometry();
			vertex_geo[i]->data = vertex_seg[state.index_data[i]]->data;
			state.vertex_shader(*vertex_seg[state.index_data[i]], *vertex_geo[i], state.uniform_data);
		}

		for (int i = 0; i < state.num_triangles; i++) {
			const data_geometry** pair_vertices = const_cast<const data_geometry**>(vertex_geo + (3 * i));
			clip_triangle(state, *pair_vertices[0], *pair_vertices[1], *pair_vertices[2], 0);
		}
		break;
		}
        case(render_type::fan):{
		data_vertex* vertex_seg[state.num_vertices];
		for (int i = 0, j = 0; i < state.num_vertices* state.floats_per_vertex; i += state.floats_per_vertex, j++) {
			vertex_seg[j] = new data_vertex[MAX_FLOATS_PER_VERTEX];
			vertex_seg[j]->data = &state.vertex_data[i];
		}

		data_geometry* vertex_geo[state.num_vertices];
		for (int i = 0; i < state.num_vertices; i++) {
			vertex_geo[i] = new data_geometry();
			vertex_geo[i]->data = vertex_seg[i]->data;
			state.vertex_shader(*vertex_seg[i], *vertex_geo[i], state.uniform_data);
		}

		for (int i = 1, j = 2; i < state.num_vertices - 1; i++, j++) {
			data_geometry* triangleGrp[3];
			triangleGrp[0] = vertex_geo[0];
			triangleGrp[1] = vertex_geo[i];
			triangleGrp[2] = vertex_geo[j];

			const data_geometry** pair_vertices = const_cast<const data_geometry**>(triangleGrp);
			clip_triangle(state, *pair_vertices[0], *pair_vertices[1], *pair_vertices[2], 0);
		}
		break;
		}
        case(render_type::strip):{
		data_vertex* vertex_seg[state.num_vertices];
		for (int i = 0, j = 0; i < state.num_vertices* state.floats_per_vertex; i += state.floats_per_vertex, j++) {
			vertex_seg[j] = new data_vertex[MAX_FLOATS_PER_VERTEX];
			vertex_seg[j]->data = &state.vertex_data[i];
		}

		data_geometry* vertex_geo[state.num_vertices];
		for (int i = 0; i < state.num_vertices; i++) {
			vertex_geo[i] = new data_geometry();
			vertex_geo[i]->data = vertex_seg[i]->data;
			state.vertex_shader(*vertex_seg[i], *vertex_geo[i], state.uniform_data);
		}

		for (int i = 0, j = 1, k = 2; i < state.num_vertices - 2; i++, j++, k++) {
			data_geometry* triangleGrp[3];
			triangleGrp[0] = vertex_geo[i];
			triangleGrp[1] = vertex_geo[j];
			triangleGrp[2] = vertex_geo[k];

			const data_geometry** pair_vertices = const_cast<const data_geometry**>(triangleGrp);
			clip_triangle(state,  *pair_vertices[0], *pair_vertices[1], *pair_vertices[2], 0);
		}
		break;
		}
        default:;
	}

}



// This function clips a triangle (defined by the three vertices in the "in" array).
// It will be called recursively, once for each clipping face (face=0, 1, ..., 5) to
// clip against each of the clipping faces in turn.  When face=6, clip_triangle should
// simply pass the call on to rasterize_triangle.
void clip_triangle(driver_state& state, const data_geometry& v0,const data_geometry& v1, const data_geometry& v2,int face)
{
    bool pos = true;
    bool a = false, b = false, c = false;
	int plane = 0;

    if(face==6)
    {
        rasterize_triangle(state, v0, v1, v2);
        return;
    }
    //cout<<"TODO: implement clipping. (The current code passes the triangle through without clipping them.)"<<endl;
	switch (face)
	{
	case 0:
		pos = true;
		plane = 0;
		a = ( v0.gl_Position[plane] <=  v0.gl_Position[3]) ? true : false;
		b = ( v1.gl_Position[plane] <=  v1.gl_Position[3]) ? true : false;
		c = ( v2.gl_Position[plane] <=  v2.gl_Position[3]) ? true : false;
		break;
	case 1:
		pos = false;
		plane = 0;
		a = ( v0.gl_Position[plane] >= - v0.gl_Position[3]) ? true : false;
		b = ( v1.gl_Position[plane] >= - v1.gl_Position[3]) ? true : false;
		c = ( v2.gl_Position[plane] >= - v2.gl_Position[3]) ? true : false;
		break;
	case 2:
		pos = true;
		plane = 1;
		a = ( v0.gl_Position[plane] <=  v0.gl_Position[3]) ? true : false;
		b = ( v1.gl_Position[plane] <=  v1.gl_Position[3]) ? true : false;
		c = ( v2.gl_Position[plane] <=  v2.gl_Position[3]) ? true : false;
		break;
	case 3:
		pos = false;
		plane = 1;
		a = ( v0.gl_Position[plane] >= - v0.gl_Position[3]) ? true : false;
		b = ( v1.gl_Position[plane] >= - v1.gl_Position[3]) ? true : false;
		c = ( v2.gl_Position[plane] >= - v2.gl_Position[3]) ? true : false;
		break;
	case 4:
		pos = true;
		plane = 2;
		a = ( v0.gl_Position[plane] <=  v0.gl_Position[3]) ? true : false;
		b = ( v1.gl_Position[plane] <=  v1.gl_Position[3]) ? true : false;
		c = ( v2.gl_Position[plane] <=  v2.gl_Position[3]) ? true : false;
		break;
	case 5:
		pos = false;
		plane = 2;
		a = ( v0.gl_Position[plane] >= - v0.gl_Position[3]) ? true : false;
		b = ( v1.gl_Position[plane] >= - v1.gl_Position[3]) ? true : false;
		c = ( v2.gl_Position[plane] >= - v2.gl_Position[3]) ? true : false;
		break;
	default:
		break;
	}

	data_geometry* t1[3];
	data_geometry* t2[3];
	for (int i = 0; i < 3; i++) {
		t1[i] = new data_geometry();
		t2[i] = new data_geometry();
	}
	float* data1 = new float[MAX_FLOATS_PER_VERTEX];
	float* data2 = new float[MAX_FLOATS_PER_VERTEX];

	if (a && b && c) {
	//All in
		clip_triangle(state,  v0,v1,v2, face + 1);
	}
	else if (a && b && !c) {
	//AB in

		t1[0]->gl_Position =  v0.gl_Position; 
		t1[0]->data = v0.data;				

		t1[1]->gl_Position =  v1.gl_Position;
		t1[1]->data = v1.data;				

		newVertex(state, t1[2], &v0, &v2, plane, pos, data1); 

		t2[0]->gl_Position =  v1.gl_Position; 
		t2[0]->data = v1.data;				

		newVertex(state, t2[1], &v1, &v2, plane, pos, data2); 

		t2[2]->gl_Position = t1[2]->gl_Position;  
		t2[2]->data = t1[2]->data;			


		clip_triangle(state, *t1[0], *t1[1], *t1[2], face + 1);
		clip_triangle(state, *t2[0], *t2[1], *t2[2], face + 1);
	}
	else if (a && !b && c) {
	//AC in

		t1[0]->gl_Position =  v2.gl_Position;	
		t1[0]->data = v2.data;				

		t1[1]->gl_Position =  v0.gl_Position; 
		t1[1]->data = v0.data;				

		newVertex(state, t1[2], &v2, &v1, plane, pos, data1); 

		t2[0]->gl_Position =  v0.gl_Position; 
		t2[0]->data = v0.data;				

		newVertex(state, t2[1], &v0, &v1, plane, pos, data2); 

		t2[2]->gl_Position = t1[2]->gl_Position;  
		t2[2]->data = t1[2]->data;				


		clip_triangle(state, *t1[0], *t1[1], *t1[2], face + 1);
		clip_triangle(state, *t2[0], *t2[1], *t2[2], face + 1);
	}
	else if (!a && b && c) {
	//BC in

		t1[0]->gl_Position =  v1.gl_Position;	
		t1[0]->data = v1.data;				

		t1[1]->gl_Position =  v2.gl_Position; 
		t1[1]->data = v2.data;				

		newVertex(state, t1[2], &v1, &v0, plane, pos, data1); 

		t2[0]->gl_Position =  v2.gl_Position; 
		t2[0]->data = v2.data;				

		newVertex(state, t2[1], &v2, &v0, plane, pos, data2); 

		t2[2]->gl_Position = t1[2]->gl_Position;	
		t2[2]->data = t1[2]->data;				


		clip_triangle(state, *t1[0],*t1[1],*t1[2], face + 1);
		clip_triangle(state, *t2[0],*t2[1],*t2[2], face + 1);
	}
	else if (a && !b && !c) {
    //A in
		t1[0]->gl_Position = v0.gl_Position; 
		t1[0]->data = v0.data;				
		newVertex(state, t1[1], &v0, &v1, plane, pos, data1); 
		newVertex(state, t1[2], &v0, &v2, plane, pos, data2); 


		clip_triangle(state, *t1[0], *t1[1], *t1[2], face + 1);
	}
	else if (!a && b && !c) {
	//B in

		t1[0]->gl_Position = v1.gl_Position; 
		t1[0]->data = v1.data;				
		newVertex(state, t1[1], &v1, &v2, plane, pos, data1);
		newVertex(state, t1[2], &v1, &v0, plane, pos, data2); 


		clip_triangle(state, *t1[0], *t1[1], *t1[2], face + 1);
	}
	else if (!a && !b && c) {
	// C in

		t1[0]->gl_Position = v2.gl_Position; 
		t1[0]->data = v2.data;				

		newVertex(state, t1[1], &v2, &v0, plane, pos, data1); 
		newVertex(state, t1[2], &v2, &v1, plane, pos, data2); 

		clip_triangle(state, *t1[0], *t1[1], *t1[2], face + 1);
	}
}

// Rasterize the triangle defined by the three vertices in the "in" array.  This
// function is responsible for rasterization, interpolation of data to
// fragments, calling the fragment shader, and z-buffering.
void rasterize_triangle(driver_state& state, const data_geometry& v0,
    const data_geometry& v1, const data_geometry& v2)
{
    //cout<<"TODO: implement rasterization"<<endl;
    int width = state.image_width;
	int height = state.image_height;


    vec2 p_ind[3];
    p_ind[0][0] = ((width / 2) * v0.gl_Position[0] / v0.gl_Position[3]) + ((width / 2) - (0.5));
	p_ind[0][1] = ((height / 2) * v0.gl_Position[1] / v0.gl_Position[3]) + ((height / 2) - (0.5)); 
    p_ind[1][0] = ((width / 2) * v1.gl_Position[0] / v1.gl_Position[3]) + ((width / 2) - (0.5));
	p_ind[1][1] = ((height / 2) * v1.gl_Position[1] / v1.gl_Position[3]) + ((height / 2) - (0.5)); 
    p_ind[2][0] = ((width / 2) * v2.gl_Position[0] / v2.gl_Position[3]) + ((width / 2) - (0.5));
	p_ind[2][1] = ((height / 2) * v2.gl_Position[1] / v2.gl_Position[3]) + ((height / 2) - (0.5)); 

	int min_x = min(min(p_ind[0][0], p_ind[1][0]), p_ind[2][0]);
	int min_y = min(min(p_ind[0][1], p_ind[1][1]), p_ind[2][1]);
	int max_x = max(max(p_ind[0][0], p_ind[1][0]), p_ind[2][0]);
	int max_y = max(max(p_ind[0][1], p_ind[1][1]), p_ind[2][1]);

	float ABC_area = triangleArea(p_ind[0], p_ind[1], p_ind[2]);

	for (int i = min_x; i <= max_x; i++) {
		for (int j = min_y; j <= max_y; j++) {

			float alpha = 0;
            float beta = 0;
            float gamma = 0;

			vec2 curr_p(i, j);
			
			alpha = triangleArea(curr_p, p_ind[1], p_ind[2]) / ABC_area;
			beta = triangleArea(curr_p, p_ind[2], p_ind[0]) / ABC_area;
			gamma = triangleArea(curr_p, p_ind[0], p_ind[1]) / ABC_area;

			if ((alpha >= 0) && (beta >= 0) && (gamma >= 0)) {

				float z = (alpha * v0.gl_Position[2] / v0.gl_Position[3]) +(beta * v1.gl_Position[2] / v1.gl_Position[3]) + (gamma * v2.gl_Position[2] / v2.gl_Position[3]);

				if (z < state.image_depth[(width * j) + i]) {

					state.image_depth[(width * j) + i] = z;

					data_output* lastc = new data_output();
					data_fragment* cdata = new data_fragment();
					float* inter_cdata = new float[MAX_FLOATS_PER_VERTEX];
					for (int k = 0; k < state.floats_per_vertex; k++) {

						switch (state.interp_rules[k]) {
						case(interp_type::noperspective): { 

							inter_cdata[k] = (alpha * v0.data[k] + beta * v1.data[k] + gamma * v2.data[k]);
							break;
						}
						case(interp_type::smooth): {
							float alpha1 = 0;
							float beta1 = 0;
							float gamma1 = 0;
							float s = 0;

							s = (alpha / v0.gl_Position[3]) + (beta / v1.gl_Position[3]) + (gamma / v2.gl_Position[3]);

							alpha1 = alpha / (v0.gl_Position[3] * s);
							beta1 = beta / (v1.gl_Position[3] * s);
							gamma1 = gamma / (v2.gl_Position[3] * s);

							inter_cdata[k] = (alpha1 * v0.data[k] + beta1 * v1.data[k] + gamma1 * v2.data[k]);
							break;
						}
						case(interp_type::flat): {
							inter_cdata[k] = v0.data[k];
							break;							
						}
						default:;
						}
					}
					cdata->data = inter_cdata;

					const data_fragment* c_cdata = const_cast<const data_fragment*>(cdata);
					state.fragment_shader(*c_cdata, *lastc, state.uniform_data);

					state.image_color[(width * j) + i] = make_pixel(lastc->output_color[0] * 255, lastc->output_color[1] * 255, lastc->output_color[2] * 255);

				}
			}

		}
	}
 
}

float triangleArea(vec2 a, vec2 b, vec2 c) {
    float A = ((b[0] * c[1]) - (c[0] * b[1]));
    float B =((a[0] * c[1]) - (c[0] * a[1]));
    float C = ((a[0] * b[1]) - (b[0] * a[1]));
	return 0.5 * (A - B + C);
}

void newVertex(driver_state& state, data_geometry* tri, const data_geometry* x, const data_geometry* y, int plane, bool pos, float* newd) {

	float alpha_s = 0, alpha_np = 0;

	if (pos)
		alpha_s = (y->gl_Position[3] - y->gl_Position[plane]) / (x->gl_Position[plane] - x->gl_Position[3] + y->gl_Position[3] - y->gl_Position[plane]);
	else
		alpha_s = (-y->gl_Position[3] - y->gl_Position[plane]) / (x->gl_Position[plane] + x->gl_Position[3] - y->gl_Position[3] - y->gl_Position[plane]);

	tri->gl_Position = alpha_s * x->gl_Position + (1 - alpha_s) * y->gl_Position;

	alpha_np = alpha_s * x->gl_Position[3] / (alpha_s * x->gl_Position[3] + (1 - alpha_s) * y->gl_Position[3]);

	for (int i = 0; i < state.floats_per_vertex; i++) {
		switch (state.interp_rules[i]) {
		case(interp_type::noperspective): {
			newd[i] = alpha_np * x->data[i] + (1 - alpha_np) * y->data[i];
			break;
		}
		case(interp_type::smooth): {
			newd[i] = alpha_s * x->data[i] + (1 - alpha_s) * y->data[i];

			break;
		}
		case(interp_type::flat): {
			newd[i] = x->data[i];
			break;
		}
		default:;
		}

	}
	tri->data = newd;

}