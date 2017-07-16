
#include "ksp_plugin/interface.hpp"

#include "glog/logging.h"

namespace principia {
namespace interface {

using ksp_plugin::Planetarium;

Planetarium* principia__PlanetariumCreate(XYZ const x_in_camera,
                                          XYZ const y_in_camera,
                                          XYZ const z_in_camera,
                                          XYZ const origin_in_camera,
                                          double const focal) {
  //LOG(WARNING) << "X: " << x_in_camera.x << " " << x_in_camera.y << " "
  //             << x_in_camera.z;
  //LOG(WARNING) << "Y: " << y_in_camera.x << " " << y_in_camera.y << " "
  //             << y_in_camera.z;
  //LOG(WARNING) << "Z: " << z_in_camera.x << " " << z_in_camera.y << " "
  //             << z_in_camera.z;
  //LOG(WARNING) << "O: " << origin_in_camera.x << " " << origin_in_camera.y
  //             << " " << origin_in_camera.z;
  //LOG(WARNING) << focal;
  return nullptr;
}

}  // namespace interface
}  // namespace principia
