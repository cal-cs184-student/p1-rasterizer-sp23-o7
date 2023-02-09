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


    sample_buffer[y * width + x] = c;
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

    float y_lower_bound = floor(std::min({y0, y1, y2}, comp));
    float y_upper_bound = floor(std::max({y0, y1, y2}, comp));
    float x_lower_bound = floor(std::min({x0, x1, x2}, comp));
    float x_upper_bound = floor(std::max({x0, x1, x2}, comp));
   
    // for (float y = (y_lower_bound + 0.5); y < y_upper_bound; y = y + 1.0){
    //   for (float x = (x_lower_bound +0.5); x < x_upper_bound; x = x + 1.0 ){
    //     float in_one = inside_line(x0, y0, x1, y1, x, y);
    //     float in_two = inside_line(x1, y1, x2, y2, x, y);
    //     float in_three = inside_line(x2, y2, x0, y0, x, y);
        
    //     if((in_one >= 0.0) && (in_two >= 0.0) && (in_three >= 0.0)){
    //       rasterize_point(x, y, color);
    //     }
    //   }
    // }

    // TODO: Task 2: Update to implement super-sampled rasterization
    //std::cout<<sample_rate<<flush;
    // set i to be for x axis,j for y axis 

    set_sample_rate(sample_rate);
    int al = sizeof(sample_buffer); ///sizeof(sample_buffer[0]); //length calculation
    //cout << "The length of the array is: " <<al;
    float sample_sqrt = 1;    
    if(sample_rate != 1){
      sample_sqrt = sqrt(sample_rate);
    }

    for(float j = 1; j <= sample_sqrt; j++ ){
      for(float i = 1; i <= sample_sqrt; i++){
        for (float y = (y_lower_bound + (j/(sample_sqrt +1))); y < y_upper_bound; y = y + 1.0){
          for (float x = (x_lower_bound +  (i/(sample_sqrt +1))); x < x_upper_bound; x = x + 1.0 ){
            float in_one = inside_line(x0, y0, x1, y1, x, y);
            float in_two = inside_line(x1, y1, x2, y2, x, y);
            float in_three = inside_line(x2, y2, x0, y0, x, y);
            
            if((in_one >= 0.0) && (in_two >= 0.0) && (in_three >= 0.0)){
              rasterize_point(x, y, color);
            }
            
          }
        }
      }
    }
    
  }




  void RasterizerImp::rasterize_interpolated_color_triangle(float x0, float y0, Color c0,
    float x1, float y1, Color c1,
    float x2, float y2, Color c2)
  {
    // TODO: Task 4: Rasterize the triangle, calculating barycentric coordinates and using them to interpolate vertex colors across the triangle
    // Hint: You can reuse code from rasterize_triangle



  }


  void RasterizerImp::rasterize_textured_triangle(float x0, float y0, float u0, float v0,
    float x1, float y1, float u1, float v1,
    float x2, float y2, float u2, float v2,
    Texture& tex)
  {
    // TODO: Task 5: Fill in the SampleParams struct and pass it to the tex.sample function.
    // TODO: Task 6: Set the correct barycentric differentials in the SampleParams struct.
    // Hint: You can reuse code from rasterize_triangle/rasterize_interpolated_color_triangle




  }

  void RasterizerImp::set_sample_rate(unsigned int rate) {
    // TODO: Task 2: You may want to update this function for supersampling support

    this->sample_rate = rate;


    // this->sample_buffer.resize(width * height, Color::White);
    //this->sample_buffer.resize((sqrt(rate) * width) * (sqrt(rate) * height), Color::White);
    this->sample_buffer.resize(rate * (width* height), Color::White);
  }



  void RasterizerImp::set_framebuffer_target(unsigned char* rgb_framebuffer,
    size_t width, size_t height)
  {
    // TODO: Task 2: You may want to update this function for supersampling support

    this->width = width; 
    this->height = height;
    this->rgb_framebuffer_target = rgb_framebuffer;


    this->sample_buffer.resize(width * height, Color::White);
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

    for (int x = 0; x < width; ++x) {
      for (int y = 0; y < height; ++y) {
        Color col = sample_buffer[y * width + x];

        for (int k = 0; k < 3; ++k) {
          this->rgb_framebuffer_target[3 * (y * width + x) + k] = (&col.r)[k] * 255;
        }
      }
    }

  }

  Rasterizer::~Rasterizer() { }


}// CGL
