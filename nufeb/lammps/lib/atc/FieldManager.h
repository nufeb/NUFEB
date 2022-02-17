#ifndef FIELD_MANAGER_H
#define FIELD_MANAGER_H

#include <string>

#include "ATC_TypeDefs.h"
#include "ATC_Error.h"
#include "PerAtomQuantity.h"

namespace ATC {
  class ATC_Method;
  enum CanonicalName {ATOMIC_TWICE_FLUCTUATING_KINETIC_ENERGY,
                      ATOMIC_TWICE_KINETIC_ENERGY,
                      ATOMIC_FLUCTUATING_VELOCITY,
                      ATOMIC_CHARGE_VELOCITY,
                      ATOMIC_SPECIES_VELOCITY,
                      PROLONGED_VELOCITY};
  typedef PerAtomQuantity<double> PAQ;
  /**
   *  @class  FieldManager 
   *  @brief  Manager for constructing fields from atomic data
   */
  class FieldManager{

  public:
    FieldManager(ATC_Method * atc);
    virtual ~FieldManager(void){};
    /** this function returns a (density) field derived from atomic data */
    DENS_MAN * nodal_atomic_field(FieldName fieldName, 
                                  std::string name = "default") {
      switch (fieldName) {
      case CHARGE_DENSITY:        return charge_density(name);
      case MASS_DENSITY:          return mass_density(name);
      case SPECIES_CONCENTRATION: return species_concentration(name);
      case NUMBER_DENSITY:        return number_density(name);
      case MOMENTUM:              return momentum(name);
      case VELOCITY:              return velocity(name);
      case PROJECTED_VELOCITY:    return projected_velocity(name);
      case DISPLACEMENT:          return displacement(name);
      case REFERENCE_POTENTIAL_ENERGY:   return reference_potential_energy(name);
      case POTENTIAL_ENERGY:      return potential_energy(name);
      case THERMAL_ENERGY:        return thermal_energy(name);
      case KINETIC_ENERGY:        return kinetic_energy(name);
      case TEMPERATURE:           return temperature(name);
      case KINETIC_TEMPERATURE:   return kinetic_temperature(name);
      case CHARGE_FLUX:           return charge_flux(name);
      case SPECIES_FLUX:          return species_flux(name);
      case INTERNAL_ENERGY:       return internal_energy(name);
      case ENERGY:                return energy(name);
      default: throw ATC_Error("FieldManager:: unknown field"); return NULL;
      }
    }
    CanonicalName string_to_canonical_name(std::string name){
       if      (name == "AtomicTwiceFluctuatingKineticEnergy") 
         return ATOMIC_TWICE_FLUCTUATING_KINETIC_ENERGY;
       else if (name == "AtomicTwiceKineticEnergy") 
         return ATOMIC_TWICE_KINETIC_ENERGY;
       else if (name == "AtomicTwiceKineticEnergy") 
         return ATOMIC_TWICE_KINETIC_ENERGY;
       else if (name == "AtomicFluctuatingVelocity") 
         return ATOMIC_FLUCTUATING_VELOCITY;
       else if (name == "AtomicChargeVelocity") // ionic current
         return ATOMIC_CHARGE_VELOCITY;
       else if (name == "AtomicSpeciesVelocity") // per species momentum 
         return ATOMIC_SPECIES_VELOCITY;
       else if (name == field_to_prolongation_name(VELOCITY))
         return PROLONGED_VELOCITY;
       else
         throw ATC_Error("unknown canonical name "+name);
    }
    PAQ * per_atom_quantity(std::string name) {
      switch (string_to_canonical_name(name)) {
      case ATOMIC_TWICE_FLUCTUATING_KINETIC_ENERGY: 
        return atomic_twice_fluctuating_kinetic_energy();
      case ATOMIC_TWICE_KINETIC_ENERGY: 
        return atomic_twice_kinetic_energy();
      case ATOMIC_FLUCTUATING_VELOCITY: 
        return atomic_fluctuating_velocity();
      case ATOMIC_CHARGE_VELOCITY:
        return atomic_charge_velocity();
      case ATOMIC_SPECIES_VELOCITY:
        return atomic_species_velocity();
      case PROLONGED_VELOCITY:
        return prolonged_field(VELOCITY);
      default: 
        throw ATC_Error("FieldManager:: unknown PAQ"); return NULL;
      }
    }
    /** this function returns a restriction of atomic data */
    DENS_MAN * restricted_atom_quantity(FieldName field, std::string name = "default", PAQ * atomi = NULL);
  protected:
    ATC_Method * atc_;
    InterscaleManager & interscaleManager_;
    // nodal atomic fields
    DENS_MAN * charge_density(std::string name);
    DENS_MAN * mass_density(std::string name);
    DENS_MAN * species_concentration(std::string name);
    DENS_MAN * number_density(std::string name);
    DENS_MAN * momentum(std::string name);
    DENS_MAN * velocity(std::string name);
    DENS_MAN * projected_velocity(std::string name);
    DENS_MAN * displacement(std::string name);
    DENS_MAN * reference_potential_energy(std::string name);
    DENS_MAN * potential_energy(std::string name);
    DENS_MAN * thermal_energy(std::string name);
    DENS_MAN * kinetic_energy(std::string name);
    DENS_MAN * temperature(std::string name);
    DENS_MAN * kinetic_temperature(std::string name);
    DENS_MAN * charge_flux(std::string name);
    DENS_MAN * species_flux(std::string name);
    DENS_MAN * internal_energy(std::string name);
    DENS_MAN * energy(std::string name);

    // non intrinsic per atom quantities (intrinsic are handled elsewhere)
    PAQ * atomic_twice_kinetic_energy();
    PAQ * atomic_twice_fluctuating_kinetic_energy();
    PAQ * atomic_fluctuating_velocity();
    PAQ * atomic_charge_velocity();
    PAQ * atomic_species_velocity();
    PAQ * atomic_species_vector();

    // internal functions
    DENS_MAN * projected_atom_quantity(FieldName field,std::string name, PAQ * atomic,  DIAG_MAN * normalization = NULL);
    DENS_MAN * scaled_projected_atom_quantity(FieldName field,std::string name, PAQ * atomic, double scale, DIAG_MAN * normalization = NULL);
    DENS_MAN * referenced_projected_atom_quantity(FieldName field, std::string name, PAQ * atomic, DENS_MAN * reference, DIAG_MAN * normalization = NULL);
    DENS_MAN * inferred_atom_quantity(FieldName /* field */, std::string /* name */, PAQ * /* atomic */){return NULL;};
    PAQ * prolonged_field(FieldName field);
  private:
    FieldManager(void);
  };
}
#endif
