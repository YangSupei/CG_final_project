#include "bbox.h"

#include "GL/glew.h"

#include <algorithm>
#include <iostream>

namespace CGL {

bool BBox::intersect(const Ray& r, double& t0, double& t1) const {

  // TODO (Part 2.2):
  // Implement ray - bounding box intersection test
  // If the ray intersected the bouding box within the range given by
  // t0, t1, update t0 and t1 with the new intersection times.

  double tx0,tx1,ty0,ty1,tz0,tz1,t00,t11;
  tx0 = std::min((min[0] - r.o[0]) / r.d[0],(max[0] - r.o[0]) / r.d[0]);
  tx1 = std::max((min[0] - r.o[0]) / r.d[0],(max[0] - r.o[0]) / r.d[0]);
  ty0 = std::min((min[1] - r.o[1]) / r.d[1],(max[1] - r.o[1]) / r.d[1]);
  ty1 = std::max((min[1] - r.o[1]) / r.d[1],(max[1] - r.o[1]) / r.d[1]);
  tz0 = std::min((min[2] - r.o[2]) / r.d[2],(max[2] - r.o[2]) / r.d[2]);
  tz1 = std::max((min[2] - r.o[2]) / r.d[2],(max[2] - r.o[2]) / r.d[2]);
  t00 = std::max(tx0,ty0);
  t00 = std::max(t00,tz0);
  t11 = std::min(tx1,ty1);
  t11 = std::min(t11,tz1);
  if(t00 <= t11 && std::max(t00,r.min_t) <= std::min(t11,r.max_t)){
    t0 = t00;
    t1 = t11;
    return true;
  }
  else{
    return false;
  };
  


}

void BBox::draw(Color c, float alpha) const {

  glColor4f(c.r, c.g, c.b, alpha);

  // top
  glBegin(GL_LINE_STRIP);
  glVertex3d(max.x, max.y, max.z);
  glVertex3d(max.x, max.y, min.z);
  glVertex3d(min.x, max.y, min.z);
  glVertex3d(min.x, max.y, max.z);
  glVertex3d(max.x, max.y, max.z);
  glEnd();

  // bottom
  glBegin(GL_LINE_STRIP);
  glVertex3d(min.x, min.y, min.z);
  glVertex3d(min.x, min.y, max.z);
  glVertex3d(max.x, min.y, max.z);
  glVertex3d(max.x, min.y, min.z);
  glVertex3d(min.x, min.y, min.z);
  glEnd();

  // side
  glBegin(GL_LINES);
  glVertex3d(max.x, max.y, max.z);
  glVertex3d(max.x, min.y, max.z);
  glVertex3d(max.x, max.y, min.z);
  glVertex3d(max.x, min.y, min.z);
  glVertex3d(min.x, max.y, min.z);
  glVertex3d(min.x, min.y, min.z);
  glVertex3d(min.x, max.y, max.z);
  glVertex3d(min.x, min.y, max.z);
  glEnd();

}

std::ostream& operator<<(std::ostream& os, const BBox& b) {
  return os << "BBOX(" << b.min << ", " << b.max << ")";
}

} // namespace CGL
