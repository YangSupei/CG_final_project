#include "bvh.h"

#include "CGL/CGL.h"
#include "triangle.h"

#include <iostream>
#include <stack>

using namespace std;

namespace CGL {
namespace SceneObjects {

BVHAccel::BVHAccel(const std::vector<Primitive *> &_primitives,
                   size_t max_leaf_size) {

  primitives = std::vector<Primitive *>(_primitives);
  root = construct_bvh(primitives.begin(), primitives.end(), max_leaf_size);
}

BVHAccel::~BVHAccel() {
  if (root)
    delete root;
  primitives.clear();
}

BBox BVHAccel::get_bbox() const { return root->bb; }

void BVHAccel::draw(BVHNode *node, const Color &c, float alpha) const {
  if (node->isLeaf()) {
    for (auto p = node->start; p != node->end; p++) {
      (*p)->draw(c, alpha);
    }
  } else {
    draw(node->l, c, alpha);
    draw(node->r, c, alpha);
  }
}

void BVHAccel::drawOutline(BVHNode *node, const Color &c, float alpha) const {
  if (node->isLeaf()) {
    for (auto p = node->start; p != node->end; p++) {
      (*p)->drawOutline(c, alpha);
    }
  } else {
    drawOutline(node->l, c, alpha);
    drawOutline(node->r, c, alpha);
  }
}

BVHNode *BVHAccel::construct_bvh(std::vector<Primitive *>::iterator start,
                                 std::vector<Primitive *>::iterator end,
                                 size_t max_leaf_size) {

  // TODO (Part 2.1):
  // Construct a BVH from the given vector of primitives and maximum leaf
  // size configuration. The starter code build a BVH aggregate with a
  // single leaf node (which is also the root) that encloses all the
  // primitives.


  BBox bbox;
  Vector3D range,centroid(0,0,0);
  std::vector<Primitive *> left,right;
  

  int chosen_idx = 0, prim_num = 0;
  for (auto p = start; p != end; p++) {
    BBox bb = (*p)->get_bbox();
    bbox.expand(bb);
    centroid += bb.centroid();
    prim_num++;
  }
  BVHNode* root = new BVHNode(bbox);
  root->start = start;
  root->end = end;

  centroid /= (double) prim_num;
  range = bbox.extent;
  if(range[1] > range[0] && range[1] > range[2]) chosen_idx = 1;
  else if(range[2] > range[0] && range[2] > range[1]) chosen_idx = 2;

  // printf("%d ",prim_num);

  if(prim_num <= max_leaf_size){
    root->start = start;
    root->end = end;
    return root;
  }
  else{
    std::vector<Primitive*>::iterator mid = start;
    for(std::vector<Primitive*>::iterator it = start; it != end; it++){
      if((*it)->get_bbox().centroid()[chosen_idx] <= centroid[chosen_idx]){
        left.push_back(*it);
      }
      else right.push_back(*it);
    }
    if(left.size() == 0 || right.size() == 0){ //consider as leaf node
      root->start = start;
      root->end = end;
      return root;
    }
    else{
      //reorder
      int cnt = 0;
      for(auto it = start; it != end; it++){
        if(cnt < left.size()){
          *it = left[cnt];
          cnt++;
          mid++;
        }
        else{
          *it = right[cnt-left.size()];
          cnt++;
        }
      }


      root->l = construct_bvh(start,mid,max_leaf_size);
      root->r = construct_bvh(mid,end,max_leaf_size);
      return root;
    }
  }

}

bool BVHAccel::has_intersection(const Ray &ray, BVHNode *node) const {
  // TODO (Part 2.3):
  // Fill in the intersect function.
  // Take note that this function has a short-circuit that the
  // Intersection version cannot, since it returns as soon as it finds
  // a hit, it doesn't actually have to find the closest hit.
  if(node == NULL || !node->bb.intersect(ray,ray.min_t,ray.max_t)) return false;

  if(node->isLeaf()){
    bool ret = false;
    for(auto prim_it = node->start; prim_it != node->end; prim_it++){
      total_isects++;
      if((*prim_it)->has_intersection(ray)){
        ret = true;
        break;
      }
    }
    return ret;
  }
  
  return BVHAccel::has_intersection(ray,node->l) || BVHAccel::has_intersection(ray,node->r);
  

}

bool BVHAccel::intersect(const Ray &ray, Intersection *i, BVHNode *node) const {
  // TODO (Part 2.3):
  // Fill in the intersect function.
  double t0,t1;
  if(node == NULL || !node->bb.intersect(ray,t0,t1)) return false;
  bool ret = false;
  if(node->isLeaf()){
    for(auto prim_it = node->start; prim_it != node->end; prim_it++){
      total_isects++;
      if((*prim_it)->intersect(ray,i)){
        ret = true;
        
      }
    }
    return ret;
  }
  bool lFlag, rFlag;
  lFlag = BVHAccel::intersect(ray,i,node->l);
  rFlag = BVHAccel::intersect(ray,i,node->r);

  return lFlag || rFlag;

  
}

} // namespace SceneObjects
} // namespace CGL
