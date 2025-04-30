#include "namelist.h"

#namelist target,struct
          STRING material = "Fe";
          double length_cm = 30.0;
          double radius_cm = 5.0;
          long include_attenuation = 1;
#end

#namelist primary_shield,struct
          STRING material = "Concrete";
          double thickness_cm = 60.0;
          double distance_cm = 400.0;
          double angle_deg = 90.0;
#end

#namelist secondary_shield,struct
          STRING material = NULL;
          double thickness_cm = 0.0;
#end

#namelist run,struct
          double beam_energy_GeV = 1.0;
          double theta_min_deg = 0;
          double theta_max_deg = 90;
          double n_theta = 91;
          STRING output = NULL;
#end

