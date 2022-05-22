#include "light.h"

#include <iostream>

#include "pathtracer/sampler.h"

namespace CGL { namespace SceneObjects {

// Directional Light //

DirectionalLight::DirectionalLight(const Vector3D rad,
                                   const Vector3D lightDir)
    : radiance(rad) {
  dirToLight = -lightDir.unit();
}

double DirectionalLight::sample_L(const Vector3D p, Vector3D* wi,
                                    double* distToLight, double* pdf,int color, double waveLength) const {
  *wi = dirToLight;
  *distToLight = INF_D;
  *pdf = 1.0;
  return 0;
}

// Infinite Hemisphere Light //

InfiniteHemisphereLight::InfiniteHemisphereLight(const Vector3D rad)
    : radiance(rad) {
  sampleToWorld[0] = Vector3D(1,  0,  0);
  sampleToWorld[1] = Vector3D(0,  0, -1);
  sampleToWorld[2] = Vector3D(0,  1,  0);
}

double InfiniteHemisphereLight::sample_L(const Vector3D p, Vector3D* wi,
                                           double* distToLight,
                                           double* pdf,int color, double waveLength) const {
  Vector3D dir = sampler.get_sample();
  *wi = sampleToWorld* dir;
  *distToLight = INF_D;
  *pdf = 1.0 / (2.0 * PI);
  return 0;
}

// Point Light //

PointLight::PointLight(const Vector3D rad, const Vector3D pos) : 
  radiance(rad), position(pos) { }

double PointLight::sample_L(const Vector3D p, Vector3D* wi,
                             double* distToLight,
                             double* pdf,int color, double waveLength) const {
  Vector3D d = position - p;
  *wi = d.unit();
  *distToLight = d.norm();
  *pdf = 1.0;
  return 0;
}


// Spot Light //

SpotLight::SpotLight(const Vector3D rad, const Vector3D pos,
                     const Vector3D dir, double angle) {

}

double SpotLight::sample_L(const Vector3D p, Vector3D* wi,
                             double* distToLight, double* pdf,int color, double waveLength) const {
  return 0;
}


// Area Light //





AreaLight::AreaLight(const Vector3D rad, 
                     const Vector3D pos,   const Vector3D dir, 
                     const Vector3D dim_x, const Vector3D dim_y)
  : radiance(rad), position(pos), direction(dir),
    dim_x(dim_x), dim_y(dim_y), area(dim_x.norm() * dim_y.norm()) { }

double AreaLight::sample_L(const Vector3D p, Vector3D* wi, 
                             double* distToLight, double* pdf,int color, double waveLength) const {

  Vector2D sample = sampler.get_sample() - Vector2D(0.5f, 0.5f);
  Vector3D d = position + sample.x * dim_x + sample.y * dim_y - p;
  double cosTheta = dot(d, direction);
  double sqDist = d.norm2();
  double dist = sqrt(sqDist);
  *wi = d / dist;
  *distToLight = dist;
  *pdf = sqDist / (area * fabs(cosTheta));
  
  double energy;
  double c = 2.998 * pow(10,8), h = 6.626 * pow(0.1,34), k = 1.3806505 * pow(0.1,23);
  double T = 7500;
  energy = 10;
  // energy = (2 * h * c * c) / pow(waveLength,5) / (exp(h*c/(waveLength*k*T*pow(0.1,6)))-1) * pow(10,30);
  // printf("%lf\n",energy);
  double radiance;
  double temp = 2700;
  temp /= 100;
  // Set Temperature = Temperature \ 100
  if(color == 0){
    double red;
    if(temp <= 66) red = 255;
    else{
      red = temp - 60;
      red = 329.698727446 * pow(red,-0.1332047592);
      if(red < 0) red = 0;
      else red = 255;
    }
    radiance = red;
  }
  else if(color == 1){
    double green;
    if(temp <= 66){
      green = temp;
      green = 99.4708025861 * log(green) - 161.1195681661;
      if(green < 0) green = 0;
      if(green > 255) green = 255;
    }
    else{
      green = temp - 60;
      green = 288.1221695283 *  pow(green,-0.0755148492);
      if(green < 0) green = 0;
      if(green > 255) green = 255;
    }
    radiance = green;
  }
  else{
    double blue;
    if(temp >= 66) blue = 255;
    else{
      if(temp <= 19) blue = 0;
      else{
        blue = temp - 10;
        blue = 138.5177312231 * log(blue) - 305.0447927307;
        if(blue < 0) blue = 0;
        if(blue > 255) blue = 255;
      }
    }
    radiance = blue;
  }
  radiance /= 255.0;
  radiance *= 10;
  return cosTheta < 0 ? radiance : 0;
};


// Sphere Light //

SphereLight::SphereLight(const Vector3D rad, const SphereObject* sphere) {

}

double SphereLight::sample_L(const Vector3D p, Vector3D* wi, 
                               double* distToLight, double* pdf,int color, double waveLength) const {

  return 0;
}

// Mesh Light

MeshLight::MeshLight(const Vector3D rad, const Mesh* mesh) {

}

double MeshLight::sample_L(const Vector3D p, Vector3D* wi, 
                             double* distToLight, double* pdf,int color, double waveLength) const {
  return 0;
}

} // namespace SceneObjects
} // namespace CGL
