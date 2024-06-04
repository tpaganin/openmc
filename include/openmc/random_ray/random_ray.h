#ifndef OPENMC_RANDOM_RAY_H
#define OPENMC_RANDOM_RAY_H

#include "openmc/memory.h"
#include "openmc/particle.h"
#include "openmc/random_ray/flat_source_domain.h"
#include "openmc/source.h"

namespace openmc {

/*
 * The RandomRay class encompasses data and methods for transporting random rays
 * through the model. It is a small extension of the Particle class.
 */

// TODO: Inherit from GeometryState instead of Particle
class RandomRay : public Particle {
public:
  //----------------------------------------------------------------------------
  // Constructors
  RandomRay();
  RandomRay(uint64_t ray_id, FlatSourceDomain* domain);

  //----------------------------------------------------------------------------
  // Methods
  void event_advance_ray();
  void attenuate_flux(double distance, bool is_active);
  void initialize_ray(uint64_t ray_id, FlatSourceDomain* domain);
  uint64_t transport_history_based_single_ray();
  //void event_advance_ray_first_collided();
  //void attenuate_flux_first_collided(double distance, bool is_active);
  //void initialize_ray_first_collided(uint64_t ray_id, FlatSourceDomain* domain);
  //uint64_t transport_history_based_single_ray_first_collided();

  //----------------------------------------------------------------------------
  // Static data members
  static double distance_inactive_;      // Inactive (dead zone) ray length
  static double distance_active_;        // Active ray length
  static unique_ptr<Source> ray_source_; // Starting source for ray sampling

  //----------------------------------------------------------------------------
  // Public data members
  vector<float> angular_flux_;

private:
  //----------------------------------------------------------------------------
  // Private data members
  vector<float> delta_psi_;
  int negroups_;
  FlatSourceDomain* domain_ {nullptr}; // pointer to domain that has flat source
                                       // data needed for ray transport
  double distance_travelled_ {0};
  bool is_active_ {false};
  bool is_alive_ {true};
}; // class RandomRay

} // namespace openmc

#endif // OPENMC_RANDOM_RAY_H
