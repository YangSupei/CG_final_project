#include "camera.h"

#include <iostream>
#include <sstream>
#include <fstream>
#include <random>
#include <chrono>
#include "pathtracer/sampler.h"
#include "bsdf.h"
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

void Camera::init_lens(){
  this->elts.clear();
  double radius[11] = {29.475,84.83,19.275,40.77,12.75,0,-14.495,40.77,-20.385,437.065,-39.73};
  double axpos[11] = {3.76,0.12,4.025,3.275,5.705,4.5,1.18,6.065,0.19,3.22,0};
  double N[11] = {1.67,1,1.67,1.699,1,0,1.603,1.658,1,1.717,1};
  double aperture[11] = {25.2,25.2,23,23,18,17.1,17,20,20,20,20};
  for(int i = 0; i < 11 ; i++){
    LensElement e;
    e.radius = radius[i];
    e.center = axpos[i];
    e.ior = N[i];
    e.aperture = aperture[i];
    elts.push_back(e);
  }

}

Vector3D Camera::back_lens_sample() const {

  // Part 1 Task 2: Implement this. Should return a point randomly sampled
  // on the back element of the lens (the element closest to the sensor)
  double aperture = elts.front().aperture;
  double radius = elts.front().radius;
  double center = elts.front().center;
  

  double theta = 2.0*M_PI*random_uniform();
  double r = aperture*.5*sqrt(random_uniform());
  double x = r*cos(theta);
  double y = r*sin(theta);
  double better_z = center - (radius > 0 ? 1 : -1) * sqrt(radius * radius - aperture * aperture * 0.25);
  return Vector3D(x,y,better_z);
}

bool Sphere_intersect(double radius,const Vector3D &center,const Ray &r,double &t1,double &t2) {
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

bool LensElement::pass_through(Ray &r, double &prev_ior) const {
  if (radius == 0) {
    double t = (center - r.o.z) / r.d.z;
    if (t < 0) return false;
    Vector3D p_intersect = r.at_time(t);
    // double actual_aperture = aperture_override != 0 ? aperture_override : aperture;
    double actual_aperture = aperture;

    if (4 * (p_intersect.x * p_intersect.x + p_intersect.y * p_intersect.y) > actual_aperture * actual_aperture) {
      return false;
    }
    return true;
  }

  // intersect with the surface (sphere) defined by (center, radius)
  double t1, t2;
  double t;
  Vector3D p_center(0, 0, center);

  if (!Sphere_intersect(radius, p_center, r, t1, t2)) {
    return false;
  }

  if (r.d.z * radius < 0) {
    // check t2 (the further-away one first)
    if (t2 >= 0) {
      t = t2;
    } else if (t1 >= 0) {
      t = t1;
    } else return false;
  } else {
    if (t1 >= 0) {
      t = t1;
    } else if (t2 >= 0) {
      t = t2;
    } else return false;
  }

  Vector3D p_intersect = r.at_time(t);

  // distance to the z-axis is greater than (aperture / 2)
  if (4 * (p_intersect.x * p_intersect.x + p_intersect.y * p_intersect.y) > aperture * aperture) {
    return false;
  }

  // refract using the Snell's law
  Vector3D normal = p_intersect - p_center;
  normal.normalize();

  // regularize normal to point in the opposite direction of the incoming ray
  if (dot(normal, r.d) > 0) {
    normal = -normal;
  }
  // assuming r.d is normalized
  double cos_theta_i = fabs(dot(r.d, normal));

  double real_ior = 1;
  double dis = 1;
  double B1, B2, B3, C1, C2, C3;
  double w2 = r.waveLength * r.waveLength;
  B1 = 1.03961212;
  B2 = 0.231792344;
  B3 = 1.01046945;
  C1 = 0.00600069867;
  C2 = 0.0200179144;
  C3 = 103.560653;
  real_ior = sqrt(1 + ((B1 * w2) / (w2 - C1)) + ((B2 * w2) / (w2 - C2)) + ((B3 * w2) / (w2 - C3)));
  
  double k = prev_ior / real_ior;
  double sin_theta_o_2 = k * k * (1 - cos_theta_i * cos_theta_i);

  // total internal reflection
  if (sin_theta_o_2 > 1) {
    return false;
  }

  Vector3D old_direction = r.d;
  r.d = k * r.d + (k * cos_theta_i - sqrt(1 - sin_theta_o_2)) * normal;
  r.d.normalize();
  
  r.o = p_intersect;

  prev_ior = real_ior;

  return true;

}

bool Camera::trace(Ray &r) const{

  double current_ior = 1; // air
  r.d.normalize();

  for (int i = 0; i < elts.size(); i++) {
    if (!elts[i].pass_through(r, current_ior)) return false;
    current_ior = elts[i].ior;
  }
  return true;
}


bool Camera::generate_ray_real_len(Ray &sampled_ray,double x, double y,int color) const {


    // Vector3D sample = lens.back_lens_sample();
    
  LensElement closest = elts[0];

  UniformGridSampler2D sampler;
  Vector2D sample = sampler.get_sample();
  double r = sample.x;
  double theta = sample.y * 2.0 * PI;
  double sensor_depth = 46.2;
  Vector3D random_position = Vector3D(r * closest.radius * cos(theta),r * closest.radius * sin(theta),closest.center - closest.radius);

  double h_fov_radian = hFov * PI / 180;
  double v_fov_radian = vFov * PI / 180;
  double w_half = tan(h_fov_radian / 2) * sensor_depth;
  double h_half = tan(v_fov_radian / 2) * sensor_depth;
  Vector3D sensor_position(w_half - x * w_half * 2,h_half - y * h_half * 2,sensor_depth);

  Vector3D direction = random_position - sensor_position;
  direction.normalize();

  sampled_ray.o = sensor_position;
  sampled_ray.d = direction;

  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  std::default_random_engine gen(seed);
  std::normal_distribution<double> R(600,25), G(550,25), B(450,15);
  switch (color){
    case 0:
      sampled_ray.waveLength = R(gen);
      break;
    case 1:
      sampled_ray.waveLength = G(gen);
      break;
    case 2:
      sampled_ray.waveLength = B(gen);
      break;
  }
  sampled_ray.color = color;
  // (cos theta)^4 for an unbiased estimate of the radiance
  // coeff = direction.z;
  // coeff *= coeff;
  // trace the sensor ray through the compound lens
  // cout <<"before" << x << " " << y << "" << sampled_ray.d << endl;
  if (trace(sampled_ray)) {
    return false;
  }
  Vector3D sampled_dir = c2w * sampled_ray.d;
  sampled_dir.normalize(); 
  sampled_ray.o = pos + c2w * sampled_ray.o;
  sampled_ray.d = sampled_dir;
  sampled_ray.min_t = nClip;
  sampled_ray.max_t = fClip;
  // cout <<"after:"<< x << " " << y << "" << sampled_ray.d << endl;

  return true;
    
 
}

//-----------------------------------------------------------------------------------------------------

// bool LensElement::pass_through(Ray &r, double &prev_ior) const {
//   // Part 1 Task 1: Implement this. It takes r and passes it through this lens element.
//   Vector3D hit_p;
//   if (intersect(r, &hit_p)){

//     if (refract(r, hit_p, prev_ior)){
//       if(radius != 0.){
//    // prev_ior = this->ior;
//        prev_ior = ior;}

//       return true;
//     }
//   }
//   return false;
// }

// bool LensElement::intersect(const Ray &r, Vector3D *hit_p) const {
//   // Part 1 Task 1: Implement this. It intersects r with this spherical lens elemnent 
//   // (or aperture diaphragm). You'll want to reuse some sphere intersect code.

//   //if aperature element: plane intersection else sphere element
//   Vector3D cen = Vector3D(0,0,center);
//   double t = 0.0;
//   if (radius == 0){
//     double t = (center-r.o[2])/r.d[2]; //intersects plane
//     return true;
//   }
//   //sphere intersect
//   double a = dot(r.d,r.d);
//   double b = dot(2*(r.o-cen), r.d);
//   double c = dot((r.o-cen),(r.o-cen))-radius*radius;
//   double d = b*b-4*a*c;
//   if (d<0)return false;
//   double tp = (-b+sqrt(d))/(2*a);
//   double tm = (-b-sqrt(d))/(2*a);
//   double t1 = min(tp,tm);
//   double t2 = max(tp,tm);
    
//   if (radius*r.d.z>0){
//     t = t1; //debugging journey
//   }else{
//     t = t2;  
//   }
//   *hit_p = r.o+r.d*t;
//   if ((hit_p->x*hit_p->x + hit_p->y*hit_p->y) > (aperture*aperture*.25)) { // missed element, debugging journey squared aperature
//       return false;
//   }
//   return true;
// }

// bool LensElement::refract(Ray& r, const Vector3D& hit_p, const double& prev_ior) const {
//   // Part 1 Task 1: Implement this. It refracts the Ray r with this lens element or 
//   // does nothing at the aperture element.
//   // You'll want to consult your refract function from the previous assignment.
//   if (radius == 0){
//     return true;
//   }

//   Matrix3x3 o2w;
//   Vector3D cen = Vector3D(0,0,center);
//   make_coord_space(o2w,(hit_p-cen)/radius);
//   Matrix3x3 w2o = o2w.T();
//   Vector3D wo = (w2o*r.d);

//   float sign = 1.0f;
//   double ni, no;
//   if (r.d.z < 0) { //forward ray needs to flip normal
//     no = prev_ior; 
//     ni = ior;
//     sign = -1.0f; 
//   }else{//backward ray
//     ni = prev_ior; 
//     no = ior;
//   }

//   Vector3D n = Vector3D(0,0,sign);
//   double p = no/ni; //no/ni
//   double c = dot(-n, wo); //- vs positive
//   Vector3D v = p*wo + (p*c - sqrt(1-p*p*(1-c*c)))*n;
//   Ray g = Ray(hit_p, o2w*v);//-wo, parentheses d= v debugging journey v to world coordinates
  
//   if (no*sin_theta(wo) >= ni){
//     return false;
//   }
//   r = g;
//   return true;


// }

// bool Camera::trace(Ray &r) const {
//   // Part 1 Task 1: Implement this. It traces a ray from the sensor out into the world.
//   double prev_ior = 1.0;
//   bool passed = true;
//   //bool passed = elts[0].pass_through(r,start_ior);  
//   for (int i=0; i < elts.size(); i++){ //elts.size()  //i=1
//     if (passed){
//       //cout << "Ray Origin: x "<< r.o.x << " y " << r.o.y << " z " << r.o.z <<endl;
//       passed = elts[i].pass_through(r, prev_ior); //elts[i-1].ior
//     }else{
//       return false;
//     }
//   }
//   return passed;
// }

// Ray Camera::generate_ray_real_len_2(double x, double y, int& rays_tried, double& cos4_term) const {

//   Ray r = Ray(Vector3D(),Vector3D());


//     // Part 1 Task 2: Implement this. It generates a ray from sensor pixel (x,y)
//     // pointing toward the back element of the lens (use back_lens_sample) and traces
//     // it through the Lens (using your "trace" function)
//     double film_d = sqrt(24*24+36*36);
//     double screen_d = sqrt(screenW*screenW + screenH*screenH);
//     double film_w = film_d * screenW / screen_d;
//     double film_h = film_d * screenH / screen_d;
//     double sensor_depth = 46.2;
//     Vector3D sensor_point(-(x-0.5)*film_w, -(y-0.5)*film_h, sensor_depth);

//     r = Ray(sensor_point, back_lens_sample()-sensor_point);
//     r.d.normalize();
//     cos4_term = pow(r.d.z,4);     
//     bool b = trace(r);      
//     while (!b) {
//        r.d = back_lens_sample()-sensor_point;
//        r.d.normalize();
//        cos4_term = pow(r.d.z,4);
//        r.o = sensor_point; //debugging error
//        b = trace(r);

//        rays_tried+=1;
//       if (rays_tried > 20){
//         r.o = Vector3D(); //debugging error
//         r.d = Vector3D(0,0,1);
//         cos4_term = 0;
//         break;
//       } 
//     }
//     //cout<<cos4_term<<endl;
//      //r.o = sensor_point;

//     /***** end of your code ******/

//     // This code converts the ray you traced through the lens into world coordinates.
//     double scale = 1;
//     r.o = pos + c2w * r.o * scale;
//     r.d = (c2w * r.d).unit();
//      //cout<<r.o.x<<"world y " << r.o.y << " z "<< r.o.z<<endl;
//     //cout<<"direction z "<<r.d.z<<endl;


//   //rays_tried = 1;
//   r.min_t = nClip; r.max_t = fClip;
//   return r;
// }


} // namespace CGL
