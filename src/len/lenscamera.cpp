#include "lenscamera.h"
#include "scene/sphere.h"
#include "util/image.h"
#include <random>
#include <chrono>
using namespace std;

namespace CGL {


/****** Helpers ******/
  

// Extract the R, G, or B channel out of an RGBA color stored in a single 32bit integer
static uint32_t red_channel(uint32_t color) {
    return (255 & (color >> (0)));
}

static uint32_t green_channel(uint32_t color) {
    return (255 & (color >> (8)));
}

static uint32_t blue_channel(uint32_t color) {
    return (255 & (color >> (16)));
}

// Convert from millimeters to meters
static const double scale = .001;




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


/****** LensElement functions ******/


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
bool LensElement::intersect(const Ray &r, Vector3D *hit_p) const {
  // Part 1 Task 1: Implement this. It intersects r with this spherical lens elemnent 
  // (or aperture diaphragm). You'll want to reuse some sphere intersect code.


  return true;
  
}
bool LensElement::refract(Ray& r, const Vector3D& hit_p, const double& prev_ior) const {
  // Part 1 Task 1: Implement this. It refracts the Ray r with this lens element or 
  // does nothing at the aperture element.
  // You'll want to consult your refract function from the previous assignment.



  return true;

}






/****** Lens functions ******/



void Lens::parse_lens_file(std::string filename) {

  ifstream infile(filename);
  cout << filename << "d";
  string line;
  double z_coord = 0;
  double z_ap;
  vector<LensElement> backwards;
  elts.clear();
  bool first = true;
  while (getline(infile, line)) {
    if (first) {
      cout << "[Lens] Loading lens file " << line << endl;
      first = false;
    }
    if (line[0] == '#')
      continue;
    stringstream ss(line);
    LensElement lens;
    double offset;
    ss >> lens.radius >> offset >> lens.ior >> lens.aperture;
    lens.center = z_coord;
    if (!lens.radius) {
      z_ap = z_coord;
    }
    z_coord += offset;
    backwards.push_back(lens);
  }
  for (int i = backwards.size() - 1; i >= 0; --i) {
    LensElement l = backwards[i];
    l.center = (l.center - z_ap) + l.radius;
    if (i) l.ior = backwards[i-1].ior;
    else l.ior = 1;
    if (!l.ior) l.ior = 1;
    elts.push_back(l);
    if (!l.radius)
      ap_i = elts.size()-1;
    // cout << "Lens element edge first " << (l.center - l.radius) << " " 
    //   << l.radius << " " << l.center << " " << l.ior << " " << l.aperture << endl;
  }
  double c = elts.front().center, r = elts.front().radius, a = elts.front().aperture * .5;
  back_elt = c - (r>0?1:-1) * sqrt(r*r-a*a);
  ap_radius = ap_original = elts[ap_i].aperture;

  // Get infinity and close focus depths, also get focal length.
  set_focus_params();
  // Focus at infinity to start.
  sensor_depth = infinity_focus;
       
}


void Lens::set_focus_params() {

  // Part 1 Task 2: Implement this. 
  // After this function is called, the three variables
  // infinity_focus, near_focus, and focal_length
  // should be set correctly.



  cout << "[Lens] Infinity focus depth is " << infinity_focus << endl;
  cout << "[Lens] Close focus depth is " << near_focus << endl;
  cout << "[Lens] True focal length is " << focal_length << endl;
}




bool Lens::trace(Ray &r, std::vector<Vector3D> *trace) const {
  // Part 1 Task 1: Implement this. It traces a ray from the sensor out into the world.

  double current_ior = 1; // air
  r.d.normalize();

  for (int i = 0; i < elts.size(); i++) {
    if (!elts[i].pass_through(r, current_ior)) return false;
    current_ior = elts[i].ior;
    if (trace) trace->push_back(r.o);
  }
  return true;
}

bool Lens::trace_backwards(Ray &r, std::vector<Vector3D> *trace) const {
  // Part 1 Task 1: Implement this. It traces a ray from the world backwards through 
  // the lens towards the sensor.
  double current_ior = 1; // air
  r.d.normalize();

  for (int i = elts.size() - 1; i >= 0; i--) {
    // printf("%lf => %lf (%lf %lf)\n", current_ior, elem.ior, elem.center, elem.radius);
    if (elts[i].pass_through(r, current_ior)) return false;
    current_ior = elts[i].ior;
    if (trace) trace->push_back(r.o);
  }
  
  return true;
}

float Lens::focus_depth(float d) const {

  // Part 1 Task 2: Implement this. Should find the conjugate of a ray
  // starting from the sensor at depth d.

  return 0;
}

Vector3D Lens::back_lens_sample() const {

  // Part 1 Task 2: Implement this. Should return a point randomly sampled
  // on the back element of the lens (the element closest to the sensor)

  return Vector3D();

}



/****** LensCamera functions ******/


LensCamera::LensCamera(): pt(NULL) {
  // string path = string(__FILE__).substr(0,string(__FILE__).find_last_of('/')+1) + "../lenses/";
  string path = "/home/yang/Desktop/CG_final_project/lenses/";
  static const vector<string> lens_files = {"dgauss.50mm.dat", "wide.22mm.dat", "telephoto.250mm.dat", "fisheye.10mm.dat"};
  for (string lens_file : lens_files)
    lenses.emplace_back(path + lens_file);

  mount_lens(0);
}


// bool LensCamera::generate_ray(Ray &sampled_ray,double x, double y,int color) const {

//   if (lens_ind >= 0) {

//     // sample a point on the rearest lens
//     Lens lens = lenses[lens_ind];
//     Vector3D sample = lens.back_lens_sample();

//     // compute the chef ray
//     double h_fov_radian = hFov * PI / 180;
//     double v_fov_radian = vFov * PI / 180;

//     double w_half = tan(h_fov_radian / 2) * lens.sensor_depth;
//     double h_half = tan(v_fov_radian / 2) * lens.sensor_depth;

//     // generate a ray from x, y on the sensor plane (at sensor_depth)
//     Vector3D sensor_position(w_half - x * w_half * 2,h_half - y * h_half * 2,lens.sensor_depth);

//     Vector3D direction = sample - sensor_position;
//     direction.normalize();

//     Ray sampled_ray(sensor_position, direction);

//     unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
//     std::default_random_engine gen(seed);
//     std::normal_distribution<double> R(600,25), G(550,25), B(450,15);
//     switch (color){
//       case 0:
//         sampled_ray.waveLength = R(gen);
//         break;
//       case 1:
//         sampled_ray.waveLength = G(gen);
//         break;
//       case 2:
//         sampled_ray.waveLength = B(gen);
//         break;
//     }
//     sampled_ray.color = color;
//     // (cos theta)^4 for an unbiased estimate of the radiance
//     // coeff = direction.z;
//     // coeff *= coeff;


//     // trace the sensor ray through the compound lens
//     if (!lens.trace(sampled_ray, NULL)) {
//       return false;
//     }

//     Vector3D sampled_dir = c2w * sampled_ray.d;
//     sampled_dir.normalize(); 


//     sampled_ray.o = pos + c2w * sampled_ray.o;
//     sampled_ray.d = sampled_dir;
//     sampled_ray.min_t = nClip;
//     sampled_ray.max_t = fClip;


//     return true;
    
//   } 
//   else {

//     // Generate ray for a pinhole camera. Same as in the previous assignment.
//     x = 2*(x-.5); y = 2*(y-.5);
//     sampled_ray = Ray(pos,(c2w*Vector3D(x*tan(radians(hFov)*.5),y*tan(radians(vFov)*.5),-1)).unit());
//     return true;

//   }
//   // return true;
 
// }

Ray LensCamera::generate_ray(double x, double y, int& rays_tried, double& cos4_term) const {

  Ray r = Ray(Vector3D(),Vector3D());

    

    double film_d = sqrt(24*24+36*36);
    double screen_d = sqrt(screenW*screenW + screenH*screenH);
    double film_w = film_d * screenW / screen_d;
    double film_h = film_d * screenH / screen_d;
    Vector3D sensor_point(-(x-0.5)*film_w, -(y-0.5)*film_h, lenses[0].sensor_depth);

    r = Ray(sensor_point, lenses[0].back_lens_sample()-sensor_point);
    r.d.normalize();
    cos4_term = pow(r.d.z,4);     
    bool b = (lenses[0]).trace(r);      
    while (!b) {
       r.d = lenses[0].back_lens_sample()-sensor_point;
       r.d.normalize();
       cos4_term = pow(r.d.z,4);
       r.o = sensor_point; //debugging error
       b = (lenses[0]).trace(r);

       rays_tried+=1;
      if (rays_tried > 20){
        r.o = Vector3D(); //debugging error
        r.d = Vector3D(0,0,1);
        cos4_term = 0;
        break;
      } 
    }
    //cout<<cos4_term<<endl;
     //r.o = sensor_point;

    /***** end of your code ******/

    // This code converts the ray you traced through the lens into world coordinates.
    r.o = pos + c2w * r.o * scale;
    r.d = (c2w * r.d).unit();
     //cout<<r.o.x<<"world y " << r.o.y << " z "<< r.o.z<<endl;
    //cout<<"direction z "<<r.d.z<<endl;

  //rays_tried = 1;
  r.min_t = nClip; r.max_t = fClip;
  return r;
}

void LensCamera::move_sensor(float delta) {
  if (lens_ind < 0) return;
  curr_lens().sensor_depth += delta;
  cout << "[LensCamera] Sensor plane moved to " << curr_lens().sensor_depth
       << ", focus now at " << lenses[lens_ind].focus_depth(lenses[lens_ind].sensor_depth) << endl;
}

void LensCamera::stop_down(float ratio) {
  float ap = curr_lens().ap_radius * ratio;
  if (ap > curr_lens().ap_original) ap = curr_lens().ap_original;
  curr_lens().ap_radius = ap;
  cout << "[LensCamera] Aperture is now " << curr_lens().ap_radius << "mm" << endl;
}

void LensCamera::mount_lens(int i) {
  lens_ind = i;
  if (i >= 0) {
    cout << "[LensCamera] Switched to lens #" << (i+1) 
         << " with focal length " << curr_lens().focal_length << "mm" << endl;
  } else {
    cout << "[LensCamera] Switched to pinhole camera" << endl;
  }
}



// A dummy function to demonstrate how to work with the image buffer.
// Calculates the average value of the green color channel in the image.
// You'll have to remember your 2D array indexing in order to take differences
// of neighboring pixels in a more sophisticated metric function.
static double mean_green(const ImageBuffer& ib) {
  double sum = 0;
  for (int i = 0; i < ib.w * ib.h; ++i) {
      sum += green_channel(ib.data[i]);
  }
  double mean = sum / (ib.w * ib.h);
  
  return mean;
}

double LensCamera::focus_metric(const ImageBuffer& ib) const {

  // Part 2 Task 1: Implement this. Design a metric to judge how "in-focus"
  // the image patch stored in the provided ImageBuffer is.

  return mean_green(ib); //  A meaningless standin
}


void LensCamera::autofocus() {


  // Part 2 Task 2: Implement this. Design a global search using your 
  // focus metric to set the sensor to be at the depth where the 
  // render cell is most "in focus". Provided code shows how to 
  // move the sensor, request a render of the cell, and evaluate the focus metric.

  // This call ensures that your pathtracer is rendering at high enough quality.
  // Increase samples per pixel to 16 and samples per light to 16.
  // pt->bump_settings();

  // Example code. Nothing to do with your actual implementation except to 
  // demonstrate functionality.
  // ImageBuffer ib;
  // curr_lens().sensor_depth += 1;
  // pt->raytrace_cell(ib);
  // cout << "[LensCamera] The mean green is " << focus_metric(ib) << endl;


  
}





void LensCamera::dump_settings(string filename) {
  ofstream file(filename);
  file << hFov << " " << vFov << " " << ar << " " << nClip << " " << fClip << endl;
  for (int i = 0; i < 3; ++i)
    file << pos[i] << " ";
  for (int i = 0; i < 3; ++i)
    file << targetPos[i] << " ";
  file << endl;
  file << phi << " " << theta << " " << r << " " << minR << " " << maxR << endl;
  for (int i = 0; i < 9; ++i)
    file << c2w(i/3, i%3) << " ";
  file << endl;
  file << screenW << " " << screenH << " " << screenDist << endl;

  file << lens_ind << endl;
  for (Lens &lens : lenses) {
    file << lens.sensor_depth << " ";
  }
  file << endl;

  cout << "[LensCamera] Dumped settings to " << filename << endl;
}

void LensCamera::load_settings(string filename) {
  ifstream file(filename);

  file >> hFov >> vFov >> ar >> nClip >> fClip;
  for (int i = 0; i < 3; ++i)
    file >> pos[i];
  for (int i = 0; i < 3; ++i)
    file >> targetPos[i];
  file >> phi >> theta >> r >> minR >> maxR;
  for (int i = 0; i < 9; ++i)
    file >> c2w(i/3, i%3);
  file >> screenW >> screenH >> screenDist;

  file >> lens_ind;
  for (Lens &lens : lenses) {
    file >> lens.sensor_depth;
  }

  cout << "[LensCamera] Loaded settings from " << filename << endl;
}


} // namespace CGL

