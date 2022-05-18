#include "bsdf.h"

#include <algorithm>
#include <iostream>
#include <utility>

#include "application/visual_debugger.h"

using std::max;
using std::min;
using std::swap;

namespace CGL {

// Mirror BSDF //

Vector3D MirrorBSDF::f(const Vector3D wo, const Vector3D wi) {
  return Vector3D();
}

Vector3D MirrorBSDF::sample_f(const Vector3D wo, Vector3D* wi, double* pdf) {

  // TODO Assignment 7: Part 1
  // Implement MirrorBSDF
  BSDF::reflect(wo,wi);
  *pdf = 1;

  return reflectance / abs_cos_theta(*wi);
}

void MirrorBSDF::render_debugger_node()
{
  if (ImGui::TreeNode(this, "Mirror BSDF"))
  {
    DragDouble3("Reflectance", &reflectance[0], 0.005);
    ImGui::TreePop();
  }
}

// Microfacet BSDF //

double MicrofacetBSDF::G(const Vector3D wo, const Vector3D wi) {
  return 1.0 / (1.0 + Lambda(wi) + Lambda(wo));
}

double MicrofacetBSDF::D(const Vector3D h) {
  // TODO Assignment 7: Part 2
  // Compute Beckmann normal distribution function (NDF) here.
  // You will need the roughness alpha.

  double theta = acos(h.z);
  double over = exp(-pow(tan(theta),2)/pow(alpha,2));
  double below = PI * pow(alpha,2) * pow(cos(theta),4);

  return over / below;
}

Vector3D MicrofacetBSDF::F(const Vector3D wi) {
  // TODO Assignment 7: Part 2
  // Compute Fresnel term for reflection on dielectric-conductor interface.
  // You will need both eta and etaK, both of which are Vector3D.
  Vector3D Rs,Rp;
  Rs = ((eta * eta + k * k) - 2 * eta * wi.z + wi.z * wi.z) /
        ((eta * eta + k * k) + 2 * eta * wi.z + wi.z * wi.z);
  Rp = ((eta * eta + k * k) * wi.z * wi.z - 2 * eta * wi.z + 1) /
        ((eta * eta + k * k) * wi.z * wi.z + 2 * eta * wi.z + 1);

  return (Rs + Rp) / 2;

}

Vector3D MicrofacetBSDF::f(const Vector3D wo, const Vector3D wi) {
  // TODO Assignment 7: Part 2
  // Implement microfacet model here.
  if( wo.z < 0 || wi.z < 0) return Vector3D();
  Vector3D h = (wo + wi);
  h.normalize();
  return (F(wi) * G(wo,wi) * D(h)) / (4 * wo.z * wi.z);
}

Vector3D MicrofacetBSDF::sample_f(const Vector3D wo, Vector3D* wi, double* pdf) {
  // TODO Assignment 7: Part 2
  // *Importance* sample Beckmann normal distribution function (NDF) here.
  // Note: You should fill in the sampled direction *wi and the corresponding *pdf,
  //       and return the sampled BRDF value.

  // *wi = cosineHemisphereSampler.get_sample(pdf);
  // return MicrofacetBSDF::f(wo, *wi);
  
  
  double r1 = sampler.get_sample().x, r2 = sampler.get_sample().y; 
  double theta_h = atan(sqrt(-alpha * alpha * log(1-r1)));
  double phi_h = 2 * PI * r2;
  Vector3D h = Vector3D(sin(theta_h) * cos(phi_h), sin(theta_h) * sin(phi_h), cos(theta_h));
  *wi = 2 * dot(h, wo) * h - wo;
  if(wi->z <= 0){
    *pdf = 0;
    return Vector3D();
  }
  double p_theta, p_phi, p_h;
  p_theta = ( 2 * sin(theta_h) / ( pow(alpha,2) * pow(cos(theta_h),3) ) ) * exp( -pow(tan(theta_h),2) / pow(alpha,2) );
  p_phi = 0.5 / PI;
  *pdf = p_theta * p_phi / sin(theta_h) / 4 / dot(*wi,h);

  return f(wo,*wi);
}

void MicrofacetBSDF::render_debugger_node()
{
  if (ImGui::TreeNode(this, "Micofacet BSDF"))
  {
    DragDouble3("eta", &eta[0], 0.005);
    DragDouble3("K", &k[0], 0.005);
    DragDouble("alpha", &alpha, 0.005);
    ImGui::TreePop();
  }
}

// Refraction BSDF //

Vector3D RefractionBSDF::f(const Vector3D wo, const Vector3D wi) {
  

  return Vector3D();
}

Vector3D RefractionBSDF::sample_f(const Vector3D wo, Vector3D* wi, double* pdf) {
  // TODO Assignment 7: Part 1
  // Implement RefractionBSDF
  if(!BSDF::refract(wo,wi,ior))
    return Vector3D();
  double eta;
  *pdf = 1;
  wo.z > 0 ? eta = 1 / ior : eta = ior; 
  return transmittance / abs_cos_theta(*wi) / (eta*eta);
}

void RefractionBSDF::render_debugger_node()
{
  if (ImGui::TreeNode(this, "Refraction BSDF"))
  {
    DragDouble3("Transmittance", &transmittance[0], 0.005);
    DragDouble("ior", &ior, 0.005);
    ImGui::TreePop();
  }
}

// Glass BSDF //

Vector3D GlassBSDF::f(const Vector3D wo, const Vector3D wi) {
  return Vector3D();
}

Vector3D GlassBSDF::sample_f(const Vector3D wo, Vector3D* wi, double* pdf) {

  // TODO Assignment 7: Part 1
  // Compute Fresnel coefficient and either reflect or refract based on it.

  // compute Fresnel coefficient and use it as the probability of reflection
  // - Fundamentals of Computer Graphics page 305

  double eta;
  if(wo.z > 0){
    eta = 1 / ior;
  }
  else{
    eta = ior;
  }

  if(!BSDF::refract(wo,wi,ior)){
    BSDF::reflect(wo,wi);
    *pdf = 1;
    return reflectance / abs_cos_theta(*wi);
  }
  else{
    double r0,r,cosTheta;
    cosTheta = (wo.z > 0) ? wo.z : -wo.z;
    r0 = pow((1-ior)/(1+ior),2);
    r = r0 + (1 - r0) * pow(1-cosTheta,5);
    if(coin_flip(r)){
      BSDF::reflect(wo,wi);
      *pdf = r;
      return r * reflectance / abs_cos_theta(*wi);
    }
    else{
      BSDF::refract(wo,wi,ior);
      *pdf = 1 - r;
      return (1 - r) * transmittance / abs_cos_theta(*wi) / pow(eta,2);
    }
  }


}

void GlassBSDF::render_debugger_node()
{
  if (ImGui::TreeNode(this, "Refraction BSDF"))
  {
    DragDouble3("Reflectance", &reflectance[0], 0.005);
    DragDouble3("Transmittance", &transmittance[0], 0.005);
    DragDouble("ior", &ior, 0.005);
    ImGui::TreePop();
  }
}

void BSDF::reflect(const Vector3D wo, Vector3D* wi) {

  // TODO Assignment 7: Part 1
  // Implement reflection of wo about normal (0,0,1) and store result in wi.
  wi->x = -wo.x;
  wi->y = -wo.y;
  wi->z = wo.z;
  return;

}

bool BSDF::refract(const Vector3D wo, Vector3D* wi, double ior) {

  // TODO Assignment 7: Part 1
  // Use Snell's Law to refract wo surface and store result ray in wi.
  // Return false if refraction does not occur due to total internal reflection
  // and true otherwise. When dot(wo,n) is positive, then wo corresponds to a
  // ray entering the surface through vacuum.
  double eta,delta;
  wo.z > 0 ? eta = 1 / ior : eta = ior;
  wi->x = -1 * eta * wo.x;
  wi->y = -1 * eta * wo.y;
  delta = 1 - eta * eta * (1 - wo.z*wo.z);
  if(delta < 0) return false;
  wi->z = sqrt(delta);
  if(wo.z > 0) wi->z = 0 - wi->z;
  return true;

}

} // namespace CGL
