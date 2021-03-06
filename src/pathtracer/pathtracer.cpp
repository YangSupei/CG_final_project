#include "pathtracer.h"

#include "scene/light.h"
#include "scene/sphere.h"
#include "scene/triangle.h"


using namespace CGL::SceneObjects;

namespace CGL {

PathTracer::PathTracer() {
  gridSampler = new UniformGridSampler2D();
  hemisphereSampler = new UniformHemisphereSampler3D();

  tm_gamma = 2.2f;
  tm_level = 1.0f;
  tm_key = 0.18;
  tm_wht = 5.0f;
}

PathTracer::~PathTracer() {
  delete gridSampler;
  delete hemisphereSampler;
}

void PathTracer::set_frame_size(size_t width, size_t height) {
  sampleBuffer.resize(width, height);
  sampleCountBuffer.resize(width * height);
}

void PathTracer::clear() {
  bvh = NULL;
  scene = NULL;
  camera = NULL;
  sampleBuffer.clear();
  sampleCountBuffer.clear();
  sampleBuffer.resize(0, 0);
  sampleCountBuffer.resize(0, 0);
}

void PathTracer::write_to_framebuffer(ImageBuffer &framebuffer, size_t x0,
                                      size_t y0, size_t x1, size_t y1) {
  sampleBuffer.toColor(framebuffer, x0, y0, x1, y1);
}

Vector3D
PathTracer::estimate_direct_lighting_hemisphere(const Ray &r,
                                                const Intersection &isect) {
  // Estimate the lighting from this intersection coming directly from a light.
  // For this function, sample uniformly in a hemisphere.

  // Note: When comparing Cornel Box (CBxxx.dae) results to importance sampling, you may find the "glow" around the light source is gone.
  // This is totally fine: the area lights in importance sampling has directionality, however in hemisphere sampling we don't model this behaviour.

  // make a coordinate system for a hit point
  // with N aligned with the Z direction.
  Matrix3x3 o2w;
  make_coord_space(o2w, isect.n);
  Matrix3x3 w2o = o2w.T();

  // w_out points towards the source of the ray (e.g.,
  // toward the camera if this is a primary ray)
  const Vector3D hit_p = r.o + r.d * isect.t;
  const Vector3D w_out = w2o * (-r.d);

  // This is the same number of total samples as
  // estimate_direct_lighting_importance (outside of delta lights). We keep the
  // same number of samples for clarity of comparison.
  int num_samples = scene->lights.size() * ns_area_light;
  Vector3D L_out(0,0,0);

  // TODO (Part 3): Write your sampling loop here
  // TODO BEFORE YOU BEGIN
  // UPDATE `est_radiance_global_illumination` to return direct lighting instead of normal shading 
  for(int i = 0; i < num_samples; i++){
    Vector3D w_in = hemisphereSampler->get_sample();
    Vector3D bsdfValue = isect.bsdf->f(w_out,w_in);
    double cosValue = dot(w_in,Vector3D(0,0,1));
    w_in = o2w * w_in; // change to world space to compute L value
    Ray w_in_ray(hit_p,w_in,1);
    Intersection w_in_intersection;
    w_in_ray.min_t = EPS_F;
    if(bvh->intersect(w_in_ray,&w_in_intersection)){
      L_out += (bsdfValue * w_in_intersection.bsdf->get_emission() * cosValue) * (2 * PI);
    }
  }
  L_out /= (double) num_samples;
  return L_out;

  // return Vector3D(1.0);

}

Vector3D
PathTracer::estimate_direct_lighting_importance(const Ray &r,
                                                const Intersection &isect) {
  // Estimate the lighting from this intersection coming directly from a light.
  // To implement importance sampling, sample only from lights, not uniformly in
  // a hemisphere.

  // make a coordinate system for a hit point
  // with N aligned with the Z direction.
  Matrix3x3 o2w;
  make_coord_space(o2w, isect.n);
  Matrix3x3 w2o = o2w.T();

  // w_out points towards the source of the ray (e.g.,
  // toward the camera if this is a primary ray)
  const Vector3D hit_p = r.o + r.d * isect.t;
  const Vector3D w_out = w2o * (-r.d);
  Vector3D L_out;

  for(auto light = scene->lights.begin(); light != scene->lights.end(); light++){
    int num_samples;
    Vector3D L_one_light;
    (*light)->is_delta_light() ? num_samples = 1 : num_samples = ns_area_light;
    for(int i = 0; i < num_samples; i++){
      Vector3D w_in,L_value;
      double dist,pdf; 
      L_value = (*light)->sample_L(hit_p,&w_in,&dist,&pdf);
      Vector3D bsdfValue = isect.bsdf->f(w_out,w_in);
      Ray w_in_ray(hit_p,w_in,1);
      w_in = w2o * w_in;
      double cosValue = dot(w_in,Vector3D(0,0,1));
      Intersection w_in_intersection;
      w_in_ray.min_t = EPS_D;
      w_in_ray.max_t = dist - EPS_D;
      if(!bvh->intersect(w_in_ray,&w_in_intersection) && cosValue >= 0){
        L_one_light += (bsdfValue * L_value * cosValue) / pdf;
      }
    }
    L_out += L_one_light / (double) num_samples;
  }
  return L_out;



  // return Vector3D(1.0);

}

Vector3D PathTracer::zero_bounce_radiance(const Ray &r,
                                          const Intersection &isect) {
  // TODO: Part 3, Task 2
  // Returns the light that results from no bounces of light


  return isect.bsdf->get_emission();


}

Vector3D PathTracer::one_bounce_radiance(const Ray &r,
                                         const Intersection &isect) {
  // TODO: Part 3, Task 3
  // Returns either the direct illumination by hemisphere or importance sampling
  // depending on `direct_hemisphere_sample`

  if(direct_hemisphere_sample) return estimate_direct_lighting_hemisphere(r,isect);
  else return estimate_direct_lighting_importance(r,isect);

}

Vector3D PathTracer::at_least_one_bounce_radiance(const Ray &r,
                                                  const Intersection &isect) {

  Matrix3x3 o2w;
  make_coord_space(o2w, isect.n);
  Matrix3x3 w2o = o2w.T();

  Vector3D hit_p = r.o + r.d * isect.t;
  Vector3D w_out = w2o * (-r.d);

  Vector3D L_out(0, 0, 0);

    if (!isect.bsdf->is_delta()) {
        L_out += one_bounce_radiance(r, isect);
    }
    
    if (r.depth < max_ray_depth && coin_flip(0.7)) {
        Vector3D w_in;
        double pdf;
        Vector3D f = isect.bsdf->sample_f(w_out, &w_in, &pdf);            
        Ray sample_ray = Ray(hit_p, (o2w * w_in).unit(), INF_D, r.depth - 1);
        sample_ray.min_t = EPS_D;
        Intersection sample_i;
        if (bvh->intersect(sample_ray, &sample_i)) {
            Vector3D l = at_least_one_bounce_radiance(sample_ray, sample_i);
            if (isect.bsdf->is_delta()) {
                l += zero_bounce_radiance(sample_ray, sample_i);
            }
            Vector3D f_l = f * l;
            L_out += f_l * abs_cos_theta(w_in) / pdf / 0.7;
        }
    }

  return L_out;

}

Vector3D PathTracer::est_radiance_global_illumination(const Ray &r) {
  Intersection isect;
  Vector3D L_out;

  // You will extend this in assignment 3-2.
  // If no intersection occurs, we simply return black.
  // This changes if you implement hemispherical lighting for extra credit.

  // The following line of code returns a debug color depending
  // on whether ray intersection with triangles or spheres has
  // been implemented.
  //
  // REMOVE THIS LINE when you are ready to begin Part 3.
  
  if (!bvh->intersect(r, &isect))
    return envLight ? envLight->sample_dir(r) : L_out;


  // L_out = (isect.t == INF_D) ? debug_shading(r.d) : normal_shading(isect.n);

  // TODO (Part 3): Return the direct illumination.
  // L_out = zero_bounce_radiance(r,isect);
  L_out = zero_bounce_radiance(r,isect) + at_least_one_bounce_radiance(r,isect);
  // TODO (Part 4): Accumulate the "direct" and "indirect"
  // parts of global illumination into L_out rather than just direct

  return L_out;
}

void PathTracer::raytrace_pixel(size_t x, size_t y) {
  // TODO (Part 1.2):
  // Make a loop that generates num_samples camera rays and traces them
  // through the scene. Return the average Vector3D.
  // You should call est_radiance_global_illumination in this function.

  // TODO (Part 5):
  // Modify your implementation to include adaptive sampling.
  // Use the command line parameters "samplesPerBatch" and "maxTolerance"

  int num_samples = ns_aa;          // total samples to evaluate
  Vector2D origin = Vector2D(x, y); // bottom left corner of the pixel
  Vector3D sum(0,0,0);
  int cnt = 0;
  double s1 = 0, s2 = 0, sigma2 = 0, miu = 0;
  for(int i = 0; i < num_samples; i++){
    cnt++;
    if(cnt % samplesPerBatch == 0){
      double I;
      miu = s1 / (double) cnt; 
      sigma2 = (1 / ((double)cnt - 1) ) * (s2 - s1 * s1 / (double)cnt);
      I = 1.96 * sqrt(sigma2 / (double)cnt);
      if(I <= maxTolerance * miu) break;
    }
    Vector2D samplePosition = origin + gridSampler->get_sample();
    Ray sampledRay;
    Vector3D temp;
    samplePosition.x /= sampleBuffer.w;
    samplePosition.y /= sampleBuffer.h;
    sampledRay = camera->generate_ray(samplePosition.x,samplePosition.y);
    sampledRay.depth = 0;
    temp = PathTracer::est_radiance_global_illumination(sampledRay);
    sum += temp;

    s1 += temp.illum();
    s2 += temp.illum() * temp.illum();
  }
  sum /= (double) cnt;
  sampleBuffer.update_pixel(sum, x, y);
  sampleCountBuffer[x + y * sampleBuffer.w] = cnt;
}

void PathTracer::autofocus(Vector2D loc) {
  Ray r = camera->generate_ray(loc.x / sampleBuffer.w, loc.y / sampleBuffer.h);
  Intersection isect;

  bvh->intersect(r, &isect);

  camera->focalDistance = isect.t;
}

} // namespace CGL
