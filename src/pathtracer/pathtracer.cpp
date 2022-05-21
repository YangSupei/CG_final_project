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

//to be converted to wavelength dependent
double
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
  double L_out = 0;

  for(int i = 0; i < num_samples; i++){
    Vector3D w_in = hemisphereSampler->get_sample();
    double bsdfValue = isect.bsdf->f(w_out,w_in,r.waveLength,r.color);
    double cosValue = dot(w_in,Vector3D(0,0,1));
    w_in = o2w * w_in; // change to world space to compute L value
    Ray w_in_ray(hit_p,w_in,1);
    Intersection w_in_intersection;
    w_in_ray.min_t = EPS_D;
    if(bvh->intersect(w_in_ray,&w_in_intersection)){
      L_out += (bsdfValue * w_in_intersection.bsdf->get_emission()[r.color] * cosValue) * (2 * PI);
    }
  }
  L_out /= (double) num_samples;
  return L_out;

  // return Vector3D(1.0);

}

double
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
  double L_out = 0;

  for(auto light = scene->lights.begin(); light != scene->lights.end(); light++){
    int num_samples;
    double L_one_light = 0;
    (*light)->is_delta_light() ? num_samples = 1 : num_samples = ns_area_light;
    for(int i = 0; i < num_samples; i++){
      Vector3D w_in;
      double dist,pdf,L_value; 
      L_value = (*light)->sample_L(hit_p,&w_in,&dist,&pdf,r.waveLength);
      double bsdfValue = isect.bsdf->f(w_out,w_in,r.waveLength,r.color);
      Ray w_in_ray(hit_p,w_in,1);
      w_in_ray.waveLength = r.waveLength;
      w_in_ray.color = r.color;
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
  // printf("%lf\n",L_out);
  return L_out;



  // return Vector3D(1.0);

}

double PathTracer::zero_bounce_radiance(const Ray &r,
                                          const Intersection &isect) {

  return isect.bsdf->get_emission()[r.color];

}

double PathTracer::one_bounce_radiance(const Ray &r,
                                         const Intersection &isect) {
  // TODO: Part 3, Task 3
  // Returns either the direct illumination by hemisphere or importance sampling
  // depending on `direct_hemisphere_sample`

  // if(direct_hemisphere_sample) return estimate_direct_lighting_hemisphere(r,isect);
  // else return estimate_direct_lighting_importance(r,isect);

  return estimate_direct_lighting_importance(r,isect);

}

double PathTracer::at_least_one_bounce_radiance(const Ray &r,
                                                  const Intersection &isect) {

  Matrix3x3 o2w;
  make_coord_space(o2w, isect.n);
  Matrix3x3 w2o = o2w.T();

  Vector3D hit_p = r.o + r.d * isect.t;
  Vector3D w_out = w2o * (-r.d);

  double L_out = 0;

    if (!isect.bsdf->is_delta()) {
        L_out += one_bounce_radiance(r, isect);
    }
    
    if (r.depth < max_ray_depth && coin_flip(0.7)) {
        Vector3D w_in;
        double pdf;
        double f = isect.bsdf->sample_f(w_out, &w_in, &pdf,r.waveLength,r.color);    
        Ray sample_ray = Ray(hit_p, (o2w * w_in).unit(), INF_D, r.depth + 1);
        sample_ray.min_t = EPS_D;
        sample_ray.waveLength = r.waveLength;
        sample_ray.color = r.color;
        Intersection sample_i;
        if (bvh->intersect(sample_ray, &sample_i)) {
            double l = at_least_one_bounce_radiance(sample_ray, sample_i);
            if (isect.bsdf->is_delta()) {
                l += zero_bounce_radiance(sample_ray, sample_i);
            }
            double f_l = f * l;
            L_out += f_l * abs_cos_theta(w_in) / pdf / 0.7;
        }
    }
  return L_out;

}

double PathTracer::est_radiance_global_illumination(const Ray &r) {
  Intersection isect;
  double L_out;

  
  if (!bvh->intersect(r, &isect))
    // return envLight ? envLight->sample_dir(r) : L_out;
    return 0;

  L_out = zero_bounce_radiance(r,isect) + at_least_one_bounce_radiance(r,isect);

  return L_out;
}


Vector3D intensity2RGB(Vector3D intensity){
  double x,y,z;
  x = intensity[0] * 1.2;
  y = intensity[1] * 1.0;
  z = intensity[2] * 1.75;

  x = x / (x + y + z);
  y = y / (x + y + z);
  z = z / (x + y + z);

  Matrix3x3 tran(0.4815, -0.1587, -0.0828,-0.0912,  0.2524,  0.0157, 0.0009, -0.0025 ,  0.1786);
  Vector3D xyz = Vector3D(x,y,z);
  Vector3D rgb = tran * xyz;
  // cout << intensity << endl;
  return rgb;
}


void PathTracer::raytrace_pixel(size_t x, size_t y) {


  int num_samples = ns_aa;          // total samples to evaluate
  Vector2D origin = Vector2D(x, y); // bottom left corner of the pixel
  Vector3D sum(0,0,0);
  int cnt = 0;
  double s1 = 0, s2 = 0, sigma2 = 0, miu = 0;
  double intensity[3] = {0,0,0};
  for(int i = 0; i < num_samples; i++){
    cnt++;
    for(int color = 0; color < 3; color++){
      Vector2D samplePosition = origin + gridSampler->get_sample();
      Ray sampledRay;
      Ray sampledRay2;
      Vector3D temp;
      samplePosition.x /= sampleBuffer.w;
      samplePosition.y /= sampleBuffer.h;
      sampledRay = camera->generate_ray(samplePosition.x,samplePosition.y,color);
      // while(!camera->generate_ray_real_len(sampledRay2,x,y,color));
      double cos_term;
      int tried;
      // while(!camera->generate_ray_real_len(samplePosition.x,samplePosition.y))
      // cout << sampledRay.d << " " << sampledRay2.d << endl;
      sampledRay.depth = 0;
      sampledRay.color = color;
      intensity[color] += PathTracer::est_radiance_global_illumination(sampledRay);
    }
    // Vector3D intensity2RGB;
    // printf("#%.2lf %.2lf %.2lf#\n",intensity[0],intensity[1],intensity[2]);
    // intensity2RGB.x = intensity[0] / (intensity[0] + intensity[1] + intensity[2]);
    // intensity2RGB.y = intensity[1] / (intensity[0] + intensity[1] + intensity[2]);
    // intensity2RGB.z = intensity[2] / (intensity[0] + intensity[1] + intensity[2]);
    // printf("#%.2lf %.2lf %.2lf#\n",intensity2RGB[0],intensity2RGB[1],intensity2RGB[2]);
    // sum += intensity2RGB;
  }

  intensity[0] /= cnt;
  intensity[1] /= cnt;
  intensity[2] /= cnt;
  // printf("%lf %lf %lf \n",intensity[0],intensity[1],intensity[2]);

  // sum = intensity2RGB(Vector3D(intensity[0],intensity[1],intensity[2]));
  // sum /= (double) cnt;
  sum = Vector3D(intensity[0],intensity[1],intensity[2]);
  sampleBuffer.update_pixel(sum, x, y);
  sampleCountBuffer[x + y * sampleBuffer.w] = cnt;
}

void PathTracer::autofocus(Vector2D loc) {
  Ray r = camera->generate_ray(loc.x / sampleBuffer.w, loc.y / sampleBuffer.h,0); //todo
  Intersection isect;

  bvh->intersect(r, &isect);

  camera->focalDistance = isect.t;
}

} // namespace CGL
