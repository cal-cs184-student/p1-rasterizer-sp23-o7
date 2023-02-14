#include "rasterizer.h"

using namespace std;



  float inside_line(float line_x0, float line_y0, float line_x1, float line_y1, float pt_x, float pt_y){
    float d_x = line_x1 - line_x0;
    float d_y = line_y1 - line_y0;
    return (-(pt_x - line_x0)*d_y + (pt_y - line_y0)*d_x);
  }

  bool comp(int a, int b)
  {
    return (a < b);
  }

float line_dist(float x, float y, float xp, float yp, float xq, float yq){
  return -(x-xp)*(yq-yp)+(y-yp)*(xq-xp);
}

namespace CGL {

  RasterizerImp::RasterizerImp(PixelSampleMethod psm, LevelSampleMethod lsm,
    size_t width, size_t height,
    unsigned int sample_rate) {
    this->psm = psm;
    this->lsm = lsm;
    this->width = width;
    this->height = height;
    this->sample_rate = sample_rate;

    sample_buffer.resize(width * height * sample_rate, Color::White);
  }

  // Used by rasterize_point and rasterize_line
  void RasterizerImp::fill_pixel(size_t x, size_t y, Color c) {
    // TODO: Task 2: You might need to this function to fix points and lines (such as the black rectangle border in test4.svg)
    // NOTE: You are not required to implement proper supersampling for points and lines
    // It is sufficient to use the same color for all supersamples of a pixel for points and lines (not triangles)

    for(int z = 0; z < sample_rate; z++){
      sample_buffer[z*(width*height)+(y * width + x)] = c; //num row * row size + col num
    }
   
    
  }

  // Rasterize a point: simple example to help you start familiarizing
  // yourself with the starter code.
  //
  void RasterizerImp::rasterize_point(float x, float y, Color color) {
    // fill in the nearest pixel
    int sx = (int)floor(x);
    int sy = (int)floor(y);

    // check bounds
    if (sx < 0 || sx >= width) return;
    if (sy < 0 || sy >= height) return;

    fill_pixel(sx, sy, color);
    return;
  }

  // Rasterize a line.
  void RasterizerImp::rasterize_line(float x0, float y0,
    float x1, float y1,
    Color color) {
    if (x0 > x1) {
      swap(x0, x1); swap(y0, y1);
    }

    float pt[] = { x0,y0 };
    float m = (y1 - y0) / (x1 - x0);
    float dpt[] = { 1,m };
    int steep = abs(m) > 1;
    if (steep) {
      dpt[0] = x1 == x0 ? 0 : 1 / abs(m);
      dpt[1] = x1 == x0 ? (y1 - y0) / abs(y1 - y0) : m / abs(m);
    }

    while (floor(pt[0]) <= floor(x1) && abs(pt[1] - y0) <= abs(y1 - y0)) {
      rasterize_point(pt[0], pt[1], color);
      pt[0] += dpt[0]; pt[1] += dpt[1];
    }
  }

  // Rasterize a triangle.
  void RasterizerImp::rasterize_triangle(float x0, float y0,
    float x1, float y1,
    float x2, float y2,
    Color color) {
    // TODO: Task 1: Implement basic triangle rasterization here, no supersampling
    // rasterize the lines 
    rasterize_line(x0, y0, x1, y1, color);
    rasterize_line(x1, y1, x2, y2, color);
    rasterize_line(x2, y2, x0, y0, color);

    float sample_sqrt = 1;    
    if(sample_rate != 1){
      sample_sqrt = sqrt(sample_rate);
    }
    

    float y_lower_bound = floor(std::min({y0, y1, y2}, comp));
    float y_upper_bound = floor(std::max({y0, y1, y2}, comp));
    float x_lower_bound = floor(std::min({x0, x1, x2}, comp));
    float x_upper_bound = floor(std::max({x0, x1, x2}, comp));

    for (float x = (x_lower_bound ); x < x_upper_bound; x = x + 1.0){
      for (float y = (y_lower_bound); y < y_upper_bound; y = y + 1.0 ){
        int z = 0;
        for(int j = 0; j < sample_sqrt; j++){
          for(int i = 0; i < sample_sqrt; i++){
            float new_x = x + (i+0.5)/sample_sqrt;
            float new_y = y + (j+0.5)/sample_sqrt;
            float in_one = inside_line(x0, y0, x1, y1, new_x, new_y);
            float in_two = inside_line(x1, y1, x2, y2, new_x, new_y);
            float in_three = inside_line(x2, y2, x0, y0, new_x, new_y);
            
            if((in_one >= 0.0) && (in_two >= 0.0) && (in_three >= 0.0)){
              //rasterize_point(x, y, color);
              sample_buffer[z*(width*height)+(y * width + x)] = color;
            }
            z++;
          }
        }
        
      }
    }
    

    // TODO: Task 2: Update to implement super-sampled rasterization


  }




  void RasterizerImp::rasterize_interpolated_color_triangle(float x0, float y0, Color c0,
    float x1, float y1, Color c1,
    float x2, float y2, Color c2)
  {
  //   // TODO: Task 4: Rasterize the triangle, calculating barycentric coordinates and using them to interpolate vertex colors across the triangle
  //   // Hint: You can reuse code from rasterize_triangle
  //   // V = aV1 + bV2 + cV3

    float sample_sqrt = 1;    
    if(sample_rate != 1){
      sample_sqrt = sqrt(sample_rate);
    }
    
    float y_lower_bound = floor(std::min({y0, y1, y2}, comp));
    float y_upper_bound = floor(std::max({y0, y1, y2}, comp));
    float x_lower_bound = floor(std::min({x0, x1, x2}, comp));
    float x_upper_bound = floor(std::max({x0, x1, x2}, comp));

    for (float x = (x_lower_bound); x <= x_upper_bound; x = x + 1.0){
      for (float y = (y_lower_bound); y <= y_upper_bound; y = y + 1.0 ){
        int z = 0;
        for(int j = 0; j < sample_sqrt; j++){
          for(int i = 0; i < sample_sqrt; i++){
            float new_x = x + (i+0.5)/sample_sqrt;
            float new_y = y + (j+0.5)/sample_sqrt;
        
            float in_one = inside_line(x0, y0, x1, y1, new_x, new_y);
            float in_two = inside_line(x1, y1, x2, y2, new_x, new_y);
            float in_three = inside_line(x2, y2, x0, y0, new_x, new_y);
                
            if((in_one >= 0.0) && (in_two >= 0.0) && (in_three >= 0.0)){
              float alpha = line_dist(new_x, new_y, x1,y1,x2,y2)/line_dist(x0, y0, x1,y1,x2,y2);
              float beta = line_dist(new_x, new_y, x0,y0,x2,y2)/line_dist(x1, y1, x0,y0,x2,y2);
              float gamma = line_dist(new_x, new_y, x1,y1,x0,y0)/line_dist(x2, y2, x1,y1,x0,y0);
              Color pt_color = alpha*c0 + beta*c1 + gamma*c2;
              sample_buffer[z*(width*height)+(y * width + x)] = pt_color;
            }
            z++;
          }
        }
        
      }
    }

  }


  void RasterizerImp::rasterize_textured_triangle(float x0, float y0, float u0, float v0,
    float x1, float y1, float u1, float v1,
    float x2, float y2, float u2, float v2,
    Texture& tex)
  {
    // TODO: Task 5: Fill in the SampleParams struct and pass it to the tex.sample function.
    // TODO: Task 6: Set the correct barycentric differentials in the SampleParams struct.
    // Hint: You can reuse code from rasterize_triangle/rasterize_interpolated_color_triangle

    SampleParams
    tex.sample()


    float sample_sqrt = 1;    
    if(sample_rate != 1){
      sample_sqrt = sqrt(sample_rate);
    }
    
    float y_lower_bound = floor(std::min({y0, y1, y2}, comp));
    float y_upper_bound = floor(std::max({y0, y1, y2}, comp));
    float x_lower_bound = floor(std::min({x0, x1, x2}, comp));
    float x_upper_bound = floor(std::max({x0, x1, x2}, comp));

    for (float x = (x_lower_bound); x <= x_upper_bound; x = x + 1.0){
      for (float y = (y_lower_bound); y <= y_upper_bound; y = y + 1.0 ){
        int z = 0;
        for(int j = 0; j < sample_sqrt; j++){
          for(int i = 0; i < sample_sqrt; i++){
            float new_x = x + (i+0.5)/sample_sqrt;
            float new_y = y + (j+0.5)/sample_sqrt;
        
            float in_one = inside_line(x0, y0, x1, y1, new_x, new_y);
            float in_two = inside_line(x1, y1, x2, y2, new_x, new_y);
            float in_three = inside_line(x2, y2, x0, y0, new_x, new_y);
                
            if((in_one >= 0.0) && (in_two >= 0.0) && (in_three >= 0.0)){
              float alpha = line_dist(new_x, new_y, x1,y1,x2,y2)/line_dist(x0, y0, x1,y1,x2,y2);
              float beta = line_dist(new_x, new_y, x0,y0,x2,y2)/line_dist(x1, y1, x0,y0,x2,y2);
              float gamma = line_dist(new_x, new_y, x1,y1,x0,y0)/line_dist(x2, y2, x1,y1,x0,y0);
              uv->x = new_x;
              uv->y = new_y;
              c = tex.sample(samp_para)

              get_texel
              Color pt_color = alpha*c0 + beta*c1 + gamma*c2;
              sample_buffer[z*(width*height)+(y * width + x)] = pt_color;
            }
            z++;
          }
        }
        
      }
    }



  }

  void RasterizerImp::set_sample_rate(unsigned int rate) {
    // TODO: Task 2: You may want to update this function for supersampling support

    this->sample_rate = rate;


    this->sample_buffer.resize(sample_rate * width * height, Color::White);
  }



  void RasterizerImp::set_framebuffer_target(unsigned char* rgb_framebuffer,
    size_t width, size_t height)
  {
    // TODO: Task 2: You may want to update this function for supersampling support

    this->width = width; 
    this->height = height;
    this->rgb_framebuffer_target = rgb_framebuffer;


    this->sample_buffer.resize((width * height), Color::White);
  }


  void RasterizerImp::clear_buffers() {
    std::fill(rgb_framebuffer_target, rgb_framebuffer_target + 3 * width * height, 255);
    std::fill(sample_buffer.begin(), sample_buffer.end(), Color::White);
  }


  // This function is called at the end of rasterizing all elements of the
  // SVG file.  If you use a supersample buffer to rasterize SVG elements
  // for antialising, you could use this call to fill the target framebuffer
  // pixels from the supersample buffer data.
  //
  void RasterizerImp::resolve_to_framebuffer() {
    // TODO: Task 2: You will likely want to update this function for supersampling support

    //where we want to sum down?
    float sample_sqrt = 1;    
    if(sample_rate != 1){
      sample_sqrt = sqrt(sample_rate);
    } 

    for (int x = 0; x < width; x++) {
      for (int y = 0; y < height; y++) {
        Color col = (0,0,0);
        for(int z = 0; z < sample_rate; z++){
          col = col + sample_buffer[z*(width*height) + (y * width + x)];
        }

        for (int k = 0; k < 3; ++k) {
          this->rgb_framebuffer_target[3 * (y * width + x) + k] = ((&col.r)[k] * float(1)/float(sample_rate)) * 255;
        }

      }
    }

  }

  Rasterizer::~Rasterizer() { }


}// CGL
