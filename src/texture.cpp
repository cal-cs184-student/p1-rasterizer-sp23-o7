#include "texture.h"
#include "CGL/color.h"

#include <cmath>
#include <algorithm>

bool newcomp(int a, int b)
  {
    return (a < b);
  }



namespace CGL {

  Color Texture:: lerp(float x, Color v0, Color v1){
    // Color diff = (v1[0]-v0[0], v1[1]-v0[1],v1[2]-v0[2]);
    // diff = (x * diff[0], diff[1] * x, diff[2]*x);
    // Color sum = (v0[0]+diff[0], v0[1] + diff[1], v0[2]+diff[2]);
    Color sum = (1-x)*v0 + x*v1;
    return sum;
  }

  Color Texture::sample(const SampleParams& sp) {
    // TODO: Task 6: Fill this in.

    float lvl = get_level(sp);
    if(sp.lsm == 0){
      int level = 0;
      if (sp.psm == 0){
        return sample_nearest(sp.p_uv, level);
      }else{
        return sample_bilinear(sp.p_uv, level);
      }
    }else if(sp.lsm ==1){
      if (sp.psm == 0){
        return sample_nearest(sp.p_uv, round(lvl));
      }else{
        return sample_bilinear(sp.p_uv, round(lvl));
      }
    }else if (sp.lsm ==2){
      int lower = std::floor(lvl);
      int higher = lower + 1;
      float diff = lvl - (lower);
      
      if (sp.psm == 0){
        Color lowerC = sample_nearest(sp.p_uv, lower);
        Color higherC = sample_nearest(sp.p_uv, higher);
        //Color temp = (1.0-diff)*lowerC + diff*higherC;
        return (1.0-diff)*lowerC + diff*higherC;
      }else{
        Color lowerC = sample_bilinear(sp.p_uv, lower);
        Color higherC = sample_bilinear(sp.p_uv, higher);
        //Color temp = (1.0-diff)*lowerC + diff*higherC;
        return (1.0-diff)*lowerC + diff*higherC;
      }
    }
    //std::cout<< level;
    // return magenta for invalid level
    return Color(0, 0, 1);
  }

  float Texture::get_level(const SampleParams& sp) {
    // TODO: Task 6: Fill this in.
    float first_sqrt = float(sqrt(float(pow(sp.p_dx_uv.x,2)) + float(pow(sp.p_dx_uv.y,2))));
    float second_sqrt = float(sqrt(float(pow(sp.p_dy_uv.x,2)) + float(pow(sp.p_dy_uv.y,2))));
    float max = log2(std::max({first_sqrt, second_sqrt}, newcomp));

    if (max < 0){
      return 0.;
    }else if (max > mipmap.size()-1){
      return float(mipmap.size()-1);
    }
    //cout<<max;
    return max;
  }

  Color MipLevel::get_texel(int tx, int ty) {
    return Color(&texels[tx * 3 + ty * width * 3]);
  }

  Color Texture::sample_nearest(Vector2D uv, int level) {
    // TODO: Task 5: Fill this in.
    auto& mip = mipmap[level];
    //std::cout<< mip;
    float tx = (mip.width-1) * uv.x;
    float ty = (mip.height-1)  * uv.y;
    // return magenta for invalid level
    return mip.get_texel(int(round(tx)), int(round(ty))); //if again, clamp with width and height
  }

  Color Texture::sample_bilinear(Vector2D uv, int level) {
    // TODO: Task 5: Fill this in.
    auto& mip = mipmap[level];

    float tx = (mip.width-1) * uv.x;
    float ty = (mip.height-1)  * uv.y;

    float top_right_x = std::ceil(tx);
    float top_right_y = std::ceil(ty);

    Color top_right = mip.get_texel(top_right_x, top_right_y);

    float top_left_x = std::floor(tx);
    float top_left_y = std::ceil(ty);

    Color top_left = mip.get_texel(top_left_x, top_left_y);

    float bottom_right_x = std::ceil(tx);
    float bottom_right_y = std::floor(ty);

    Color bottom_right = mip.get_texel(bottom_right_x, bottom_right_y);

    float bottom_left_x = std::floor(tx);
    float bottom_left_y = std::floor(ty);

    Color bottom_left = mip.get_texel(bottom_left_x, bottom_left_y);

    float s = tx - bottom_left_x;
    float t = ty - bottom_left_y;

    Color u0_x = lerp(s ,bottom_left, bottom_right);
    //Color u0_y = lerp(s ,bottom_left_y, bottom_right_y);


    Color u1_x = lerp(s ,top_left, top_right);
    //Color u1_y = lerp(s ,top_left_y, top_right_y);

    Color f_x = lerp(t ,u1_x, u0_x);
    //Color f_y = lerp(t ,u1_y, u0_y);

    // return magenta for invalid level
    return f_x;
  }



  /****************************************************************************/

  // Helpers

  inline void uint8_to_float(float dst[3], unsigned char* src) {
    uint8_t* src_uint8 = (uint8_t*)src;
    dst[0] = src_uint8[0] / 255.f;
    dst[1] = src_uint8[1] / 255.f;
    dst[2] = src_uint8[2] / 255.f;
  }

  inline void float_to_uint8(unsigned char* dst, float src[3]) {
    uint8_t* dst_uint8 = (uint8_t*)dst;
    dst_uint8[0] = (uint8_t)(255.f * max(0.0f, min(1.0f, src[0])));
    dst_uint8[1] = (uint8_t)(255.f * max(0.0f, min(1.0f, src[1])));
    dst_uint8[2] = (uint8_t)(255.f * max(0.0f, min(1.0f, src[2])));
  }

  void Texture::generate_mips(int startLevel) {

    // make sure there's a valid texture
    if (startLevel >= mipmap.size()) {
      std::cerr << "Invalid start level";
    }

    // allocate sublevels
    int baseWidth = mipmap[startLevel].width;
    int baseHeight = mipmap[startLevel].height;
    int numSubLevels = (int)(log2f((float)max(baseWidth, baseHeight)));

    numSubLevels = min(numSubLevels, kMaxMipLevels - startLevel - 1);
    mipmap.resize(startLevel + numSubLevels + 1);

    int width = baseWidth;
    int height = baseHeight;
    for (int i = 1; i <= numSubLevels; i++) {

      MipLevel& level = mipmap[startLevel + i];

      // handle odd size texture by rounding down
      width = max(1, width / 2);
      //assert (width > 0);
      height = max(1, height / 2);
      //assert (height > 0);

      level.width = width;
      level.height = height;
      level.texels = vector<unsigned char>(3 * width * height);
    }

    // create mips
    int subLevels = numSubLevels - (startLevel + 1);
    for (int mipLevel = startLevel + 1; mipLevel < startLevel + subLevels + 1;
      mipLevel++) {

      MipLevel& prevLevel = mipmap[mipLevel - 1];
      MipLevel& currLevel = mipmap[mipLevel];

      int prevLevelPitch = prevLevel.width * 3; // 32 bit RGB
      int currLevelPitch = currLevel.width * 3; // 32 bit RGB

      unsigned char* prevLevelMem;
      unsigned char* currLevelMem;

      currLevelMem = (unsigned char*)&currLevel.texels[0];
      prevLevelMem = (unsigned char*)&prevLevel.texels[0];

      float wDecimal, wNorm, wWeight[3];
      int wSupport;
      float hDecimal, hNorm, hWeight[3];
      int hSupport;

      float result[3];
      float input[3];

      // conditional differentiates no rounding case from round down case
      if (prevLevel.width & 1) {
        wSupport = 3;
        wDecimal = 1.0f / (float)currLevel.width;
      }
      else {
        wSupport = 2;
        wDecimal = 0.0f;
      }

      // conditional differentiates no rounding case from round down case
      if (prevLevel.height & 1) {
        hSupport = 3;
        hDecimal = 1.0f / (float)currLevel.height;
      }
      else {
        hSupport = 2;
        hDecimal = 0.0f;
      }

      wNorm = 1.0f / (2.0f + wDecimal);
      hNorm = 1.0f / (2.0f + hDecimal);

      // case 1: reduction only in horizontal size (vertical size is 1)
      if (currLevel.height == prevLevel.height) {
        //assert (currLevel.height == 1);

        for (int i = 0; i < currLevel.width; i++) {
          wWeight[0] = wNorm * (1.0f - wDecimal * i);
          wWeight[1] = wNorm * 1.0f;
          wWeight[2] = wNorm * wDecimal * (i + 1);

          result[0] = result[1] = result[2] = 0.0f;

          for (int ii = 0; ii < wSupport; ii++) {
            uint8_to_float(input, prevLevelMem + 3 * (2 * i + ii));
            result[0] += wWeight[ii] * input[0];
            result[1] += wWeight[ii] * input[1];
            result[2] += wWeight[ii] * input[2];
          }

          // convert back to format of the texture
          float_to_uint8(currLevelMem + (3 * i), result);
        }

        // case 2: reduction only in vertical size (horizontal size is 1)
      }
      else if (currLevel.width == prevLevel.width) {
        //assert (currLevel.width == 1);

        for (int j = 0; j < currLevel.height; j++) {
          hWeight[0] = hNorm * (1.0f - hDecimal * j);
          hWeight[1] = hNorm;
          hWeight[2] = hNorm * hDecimal * (j + 1);

          result[0] = result[1] = result[2] = 0.0f;
          for (int jj = 0; jj < hSupport; jj++) {
            uint8_to_float(input, prevLevelMem + prevLevelPitch * (2 * j + jj));
            result[0] += hWeight[jj] * input[0];
            result[1] += hWeight[jj] * input[1];
            result[2] += hWeight[jj] * input[2];
          }

          // convert back to format of the texture
          float_to_uint8(currLevelMem + (currLevelPitch * j), result);
        }

        // case 3: reduction in both horizontal and vertical size
      }
      else {

        for (int j = 0; j < currLevel.height; j++) {
          hWeight[0] = hNorm * (1.0f - hDecimal * j);
          hWeight[1] = hNorm;
          hWeight[2] = hNorm * hDecimal * (j + 1);

          for (int i = 0; i < currLevel.width; i++) {
            wWeight[0] = wNorm * (1.0f - wDecimal * i);
            wWeight[1] = wNorm * 1.0f;
            wWeight[2] = wNorm * wDecimal * (i + 1);

            result[0] = result[1] = result[2] = 0.0f;

            // convolve source image with a trapezoidal filter.
            // in the case of no rounding this is just a box filter of width 2.
            // in the general case, the support region is 3x3.
            for (int jj = 0; jj < hSupport; jj++)
              for (int ii = 0; ii < wSupport; ii++) {
                float weight = hWeight[jj] * wWeight[ii];
                uint8_to_float(input, prevLevelMem +
                  prevLevelPitch * (2 * j + jj) +
                  3 * (2 * i + ii));
                result[0] += weight * input[0];
                result[1] += weight * input[1];
                result[2] += weight * input[2];
              }

            // convert back to format of the texture
            float_to_uint8(currLevelMem + currLevelPitch * j + 3 * i, result);
          }
        }
      }
    }
  }

}
