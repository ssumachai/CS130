#include "driver_state.h"
#include <cstring>
#include <limits>

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
    state.image_color=0;
    state.image_depth=0;

    int allocate_size = width * height;
    //std::cout<<"TODO: allocate and initialize state.image_color and state.image_depth."<<std::endl;

    state.image_color = new pixel[allocate_size];
    state.image_depth = new float[allocate_size];

    for(int i = 0; i < allocate_size; ++i){
        state.image_color[i] = make_pixel(0, 0, 0);
        state.image_depth[i] = std::numeric_limits<float>::max();
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
    //std::cout<<"TODO: implement rendering."<<std::endl;

    int triangles;
    int index = 0;


    data_geometry* triArray = new data_geometry[3];
    data_vertex vertexStuff;

    switch(type){
        case render_type::triangle:
            triangles = state.num_vertices / 3;

            for(int i = 0; i < triangles; ++i){
                for(int j = 0; j < 3; ++j){
                    triArray[j].data = state.vertex_data + index;
                    index += state.floats_per_vertex;
                }
                for(int k = 0; k < 3; ++k){
                    vertexStuff.data = triArray[k].data;
                    state.vertex_shader(vertexStuff, triArray[k], state.uniform_data);
                }
                rasterize_triangle(state, triArray[0], triArray[1], triArray[2]);
            }
            break;
        case render_type::indexed:
        case render_type::fan:
        case render_type::strip:
            break;
        default:
            break;
    }
    delete[] triArray;
}

// This function clips a triangle (defined by the three vertices in the "in" array).
// It will be called recursively, once for each clipping face (face=0, 1, ..., 5) to
// clip against each of the clipping faces in turn.  When face=6, clip_triangle should
// simply pass the call on to rasterize_triangle.
void clip_triangle(driver_state& state, const data_geometry& v0,
    const data_geometry& v1, const data_geometry& v2,int face)
{
    if(face==6)
    {
        rasterize_triangle(state, v0, v1, v2);
        return;
    }
    std::cout<<"TODO: implement clipping. (The current code passes the triangle through without clipping them.)"<<std::endl;
    clip_triangle(state, v0, v1, v2,face+1);
}

// Rasterize the triangle defined by the three vertices in the "in" array.  This
// function is responsible for rasterization, interpolation of data to
// fragments, calling the fragment shader, and z-buffering.
void rasterize_triangle(driver_state& state, const data_geometry& v0,
    const data_geometry& v1, const data_geometry& v2)
{
    //std::cout<<"TODO: implement rasterization"<<std::endl;

    int x[3];
    int y[3];
    int z[3];

    int min_x, max_x, min_y, max_y = 0;

    int width_over_two = state.image_width/2.0f;
    int height_over_two = state.image_height/2.0f;

    x[0] = width_over_two * (v0.gl_Position[0]/v0.gl_Position[3]) + (width_over_two - 0.5);
    y[0] = height_over_two * (v0.gl_Position[1]/v0.gl_Position[3]) + (height_over_two - 0.5);
    z[0] = width_over_two * (v0.gl_Position[2]/v0.gl_Position[3]) + (width_over_two - 0.5);
    x[1] = width_over_two * (v1.gl_Position[0]/v1.gl_Position[3]) + (width_over_two - 0.5);
    y[1] = height_over_two * (v1.gl_Position[1]/v1.gl_Position[3]) + (height_over_two - 0.5);
    z[1] = width_over_two * (v1.gl_Position[2]/v1.gl_Position[3]) + (width_over_two - 0.5);
    x[2] = width_over_two * (v2.gl_Position[0]/v2.gl_Position[3]) + (width_over_two - 0.5);
    y[2] = height_over_two * (v2.gl_Position[1]/v2.gl_Position[3]) + (height_over_two - 0.5);
    z[2] = width_over_two * (v2.gl_Position[2]/v2.gl_Position[3]) + (width_over_two - 0.5);


    float* data_stuff = new float[MAX_FLOATS_PER_VERTEX];
    data_fragment fragment_data{data_stuff};
    data_output fragment_output;

    min_x = std::min(std::min(x[0], x[1]), x[2]);
    max_x = std::max(std::max(x[0], x[1]), x[2]);
    min_y = std::min(std::min(y[0], y[1]), y[2]);
    max_y = std::max(std::max(y[0], y[1]), y[2]);

    if(min_x < 0){min_x = 0;}
    if(max_x > state.image_width){max_x = state.image_width - 1;}
    if(min_y < 0){min_y = 0;}
    if(max_y > state.image_height){max_y = state.image_height - 1;}

    //Given from lab pdf, where area = 0.5f((b_x*c_y - c_x*b_y) + (c_x*a_y - a_x*c_y) + (a_x*b_y - b_x*a_y))
    float triArea = (0.5f * ((x[1]*y[2] - x[2]*y[1]) + (x[2]*y[0] - x[0]*y[2]) + (x[0]*y[1] - x[1]*y[0])));


    for(int i = min_x; i < max_x + 1; ++i){
        for(int j = min_y; j < max_y + 1; ++j){
            float alpha_prime = (0.5f * ((x[1]*y[2] - x[2]*y[1]) + (y[1] - y[2])*i + (x[2] - x[1])*j)) / triArea;
            float beta_prime =  (0.5f * ((x[2]*y[0] - x[0]*y[2]) + (y[2] - y[0])*i + (x[0] - x[2])*j)) / triArea;
            float gamma_prime = (0.5f * ((x[0]*y[1] - x[1]*y[0]) + (y[0] - y[1])*i + (x[1] - x[0])*j)) / triArea;

            if(alpha_prime >= 0  && beta_prime >= 0 && gamma_prime >= 0){
                    float alpha_true = alpha_prime;
                    float beta_true = beta_prime;
                    float gamma_true = gamma_prime;
                    float z_buff = (alpha_true * z[0]) + (beta_true * z[1]) + (gamma_true * z[2]);

                    if(z_buff < state.image_depth[i + j * state.image_width]){
                        state.image_depth[i + j * state.image_width] = z_buff;

                        for(int k = 0; k < state.floats_per_vertex; k++){
                            float temp_k = 0;
                            
                            switch(state.interp_rules[k]){
                                case interp_type::flat:
                                    fragment_data.data[k] = v0.data[k];
                                    break;
                                case interp_type::smooth:
                                    temp_k = ((alpha_prime / v0.gl_Position[3]) + (beta_prime / v1.gl_Position[3]) + (gamma_prime / v2.gl_Position[3]));
                                    alpha_prime = alpha_true / (temp_k * v0.gl_Position[3]);
                                    beta_prime = beta_true / (temp_k * v1.gl_Position[3]);
                                    gamma_prime = gamma_true / (temp_k * v2.gl_Position[3]);
                                    break;
                                case interp_type::noperspective:
                                    fragment_data.data[k] = ((alpha_prime * v0.data[k]) + (beta_prime * v1.data[k]) + (gamma_prime * v2.data[k]));
                                    break;
                                default:
                                    break;
                            }
                        }

                        state.fragment_shader(fragment_data, fragment_output, state.uniform_data);
                        int color_r = fragment_output.output_color[0] * 255;
                        int color_g = fragment_output.output_color[1] * 255;
                        int color_b = fragment_output.output_color[2] * 255;
                        state.image_color[i + j * state.image_width] = make_pixel(color_r, color_g, color_b);

                    }
                }
            }
        }

        delete[] data_stuff;
}

    

