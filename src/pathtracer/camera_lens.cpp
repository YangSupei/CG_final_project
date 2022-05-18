#include "camera.h"

#include <iostream>
#include <sstream>
#include <fstream>

#include "CGL/misc.h"
#include "CGL/vector2D.h"
#include "CGL/vector3D.h"

using std::cout;
using std::endl;
using std::max;
using std::min;
using std::ifstream;
using std::ofstream;

namespace CGL {

using Collada::CameraInfo;

Ray Camera::generate_ray_for_thin_lens(double x, double y, double rndR, double rndTheta) const {

  // TODO Assignment 7: Part 4
  // compute position and direction of ray from the input sensor sample coordinate.
  // Note: use rndR and rndTheta to uniformly sample a unit disk.

  double width = tan((hFov*PI/180)*0.5)*2, height = tan((vFov*PI/180)*0.5)*2;
  Vector3D d,redDir = Vector3D(x*width-width*0.5,y*height-height*0.5,-1);

  Vector3D pLens = Vector3D(lensRadius * sqrt(rndR) * cos(rndTheta),lensRadius * sqrt(rndR) * sin(rndTheta),0);
  Vector3D pFocus = redDir / (-1) * (-focalDistance);

  d = pFocus - pLens;
  d.normalize();

  Ray ret = Ray(c2w * pLens + pos, c2w * d);
  ret.max_t = fClip;
  ret.min_t = nClip;
  
  return ret;

}


} // namespace CGL
