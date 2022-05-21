#include "sphere.h"

#include <cmath>

#include "pathtracer/bsdf.h"
#include "util/sphere_drawing.h"

namespace CGL {
namespace SceneObjects {

bool Sphere::test(const Ray &r, double &t1, double &t2) const {

  // TODO (Part 1.4):
  // Implement ray - sphere intersection test.
  // Return true if there are intersections and writing the
  // smaller of the two intersection times in t1 and the larger in t2.
  double a,b,c,t,delta; 
  a = dot(r.d,r.d);
  b = 2 * dot(r.o - this->o,r.d);
  c = dot(r.o - this->o,r.o - this->o) - this->r2;
  delta = b*b-4*a*c;
  
  if(delta < 0) return false;
  else{
    t1 = (-b - sqrt(delta)) / (2 * a);
    t2 = (-b + sqrt(delta)) / (2 * a);
    if(t1 >= r.min_t && t1 <= r.max_t) return true;
  }
  

}

bool Sphere::has_intersection(const Ray &r) const {

  // TODO (Part 1.4):
  // Implement ray - sphere intersection.
  // Note that you might want to use the the Sphere::test helper here.
  double t1,t2;
  if(test(r,t1,t2)){
    if(r.min_t < t1 && r.max_t > t1){
      r.max_t = t1;
      return true;
    }
    else if(r.min_t < t2 && r.max_t > t2){
      r.max_t  = t2;
      return true;
    }
  }
  else return false;

}

bool Sphere::intersect(const Ray &r, Intersection *i) const {

  // TODO (Part 1.4):
  // Implement ray - sphere intersection.
  // Note again that you might want to use the the Sphere::test helper here.
  // When an intersection takes place, the Intersection data should be updated
  // correspondingly.
  double t1,t2;
  if(test(r,t1,t2)){
    Vector3D intersect_point = r.o + t1 * r.d;
    Vector3D intersect_normal = (intersect_point - o);
    intersect_normal.normalize();
    r.max_t = t1;
    i->t = t1;
    i->n = intersect_normal;
    i->bsdf = get_bsdf();
    i->primitive = this;
    return true;
  }
  else return false;
}

bool Sphere::intersect(double radius,const Vector3D &center,const Ray &r,double &t1,double &t2) {
  Vector3D o_to_center = r.o - center;

  double a = dot(r.d, r.d);
  double b = 2 * dot(o_to_center, r.d);
  double c = dot(o_to_center, o_to_center) - radius * radius;

  double delta = b * b - 4 * a * c;
  if (delta < 0) return false;

  t1 = (-b - sqrt(delta)) / (2 * a);
  t2 = (-b + sqrt(delta)) / (2 * a);
  
  return true;
}

void Sphere::draw(const Color &c, float alpha) const {
  Misc::draw_sphere_opengl(o, r, c);
}

void Sphere::drawOutline(const Color &c, float alpha) const {
  // Misc::draw_sphere_opengl(o, r, c);
}

} // namespace SceneObjects
} // namespace CGL
