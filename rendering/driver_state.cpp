#include "driver_state.h"
#include <cstring>

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

    for(int i = 0; i < allocate_size; ++i){
        state.image_color[i] = make_pixel(0, 0, 0);
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

    int width_over_two = state.image_width/2.0f;
    int height_over_two = state.image_height/2.0f;

    x[0] = width_over_two * v0.gl_Position[0] + (width_over_two - 0.5);
    y[0] = height_over_two * v0.gl_Position[1] + (height_over_two - 0.5);
    x[1] = width_over_two * v1.gl_Position[0] + (width_over_two - 0.5);
    y[1] = height_over_two * v1.gl_Position[1] + (height_over_two - 0.5);
    x[2] = width_over_two * v2.gl_Position[0] + (width_over_two - 0.5);
    y[2] = height_over_two * v2.gl_Position[1] + (height_over_two - 0.5);

    //Given from lab pdf, where area = 0.5f((b_x*c_y - c_x*b_y) + (c_x*a_y - a_x*c_y) + (a_x*b_y - b_x*a_y))
    float triArea = (0.5f * ((x[1]*y[2] - x[2]*y[1]) + (x[2]*y[0] - x[0]*y[2]) + (x[0]*y[1] - x[1]*y[0])));


    for(int i = 0; i < state.image_width - 1; ++i){
        for(int j = 0; j < state.image_height - 1; ++j){
            float alpha = (0.5f * ((x[1]*y[2] - x[2]*y[1]) + (y[1] - y[2])*i + (x[2] - x[1])*j)) / triArea;
            float beta =  (0.5f * ((x[2]*y[0] - x[0]*y[2]) + (y[2] - y[0])*i + (x[0] - x[2])*j)) / triArea;
            float gamma = (0.5f * ((x[0]*y[1] - x[1]*y[0]) + (y[0] - y[1])*i + (x[1] - x[0])*j)) / triArea;

            if(alpha >= 0  && beta >= 0 && gamma >= 0){
                    state.image_color[i + j * state.image_width] = make_pixel(255, 255, 255);
                }
            }
        }
}

    

