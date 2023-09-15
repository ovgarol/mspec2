#include "fabm_driver.h"

!! zooplankton for mixotrophic adaptation and optimality: a component of mspec

module hereon_zoo
   use fabm_types

   implicit none

   private

   type, extends(type_base_model), public :: type_hereon_zoo
      ! Variable identifiers
      type (type_state_variable_id) :: id_m, id_q
      type (type_state_variable_id) :: id_exctarget, id_morttarget, id_grztarget, id_grztargetq

      ! Model parameters
      real(rk) :: m0,  rmax,  iv,  eff, rmn, rmd
      logical :: ia
   contains
      procedure :: initialize
      procedure :: do
      !procedure :: do_ppdd
   end type

contains

   subroutine initialize(self, configunit)
      class (type_hereon_zoo), intent(inout), target :: self
      integer,                        intent(in)            :: configunit

      real(rk), parameter :: d_per_s = 1.0_rk/86400.0_rk

      ! Store parameter values in our own derived type
      ! NB: all rates must be provided in values per day and are converted here to values per second.
      call self%get_parameter(self%m0,   'm0',   'mmol m-3',  'background concentration',      default=0.0225_rk)
      call self%get_parameter(self%rmax, 'rmax', 'd-1',       'maximum specific growth rate', default=1.5_rk, scale_factor=d_per_s)
      call self%get_parameter(self%iv,   'iv',   'm3 mmol-1', 'Ivlev grazing constant',        default=1.1_rk)
      call self%get_parameter(self%eff,  'eff',  'non-dim',   'grazing efficiency',        default=0.91_rk)
      call self%get_parameter(self%rmn,  'rmn',  'd-1',       'excretion rate',                default=0.01_rk, scale_factor=d_per_s)
      call self%get_parameter(self%rmd,  'rmd',  'd-1',       'mortality',                     default=0.02_rk, scale_factor=d_per_s)
      call self%get_parameter(self%ia, 'ia', '-', 'instantaneous acclimation', default=.false.)

      ! Register state variables
      call self%register_state_variable(self%id_m, 'c', 'mmol m-3', 'concentration', 0.0_rk, minimum=0.0_rk, vertical_movement=-0.0_rk)
      call self%register_state_variable(self%id_q, 'q', 'non-dim', 'nutrient quota', 0.5_rk, minimum=0.0_rk, maximum=1.0_rk)

      ! Register contribution of state to global aggregate variables.
      call self%add_to_aggregate_variable(standard_variables%total_nitrogen, self%id_m)

      ! Register dependencies on external state variables.
      call self%register_state_dependency(self%id_grztarget,  'grazing_target',   'mmol m-3', 'prey source')
      call self%register_state_dependency(self%id_grztargetq, 'grazing_target_nutrient',   'mmol m-3', 'prey source nutrient content')
      call self%register_state_dependency(self%id_exctarget,  'excretion_target', 'mmol m-3', 'sink for excreted matter')
      call self%register_state_dependency(self%id_morttarget, 'mortality_target', 'mmol m-3', 'sink for dead matter')
   end subroutine initialize

   subroutine do(self, _ARGUMENTS_DO_)
      class (type_hereon_zoo), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      real(rk), parameter :: d_per_s = 1.0_rk/86400.0_rk

      real(rk) :: m, q, p, qp
      real(rk) :: graz, nutuptake
      real(rk) :: ch, cpn, gphot

      ! Enter spatial loops (if any)
      _LOOP_BEGIN_

         ! Retrieve current (local) state variable values.
         _GET_(self%id_m,m)         ! mixotroph
         _GET_(self%id_q,q)         ! nutritional state
         _GET_(self%id_grztarget,p) ! prey
         _GET_(self%id_grztargetq,qp) ! prey nutrients

         ! adaptation and growth rate calculation
         !!!   AQUI ME QUEDE
        
         ch = max(1.0_rk - exp(-(self%iv*p)**2),0.0_rk)
         cpn = ch*qp

         gphot = self%rmax!*(1.0_rk-exp(-q**2)) !here is multiplied a temperature dependence term
             
         graz = ch*gphot
         nutuptake = cpn*gphot

         !!!
         ! Set temporal derivatives
         _ADD_SOURCE_(self%id_m, graz*self%eff*m - self%rmd*m)
         _SET_ODE_(self%id_q, nutuptake*self%eff-q*graz*self%eff)
         
         _ADD_SOURCE_(self%id_grztarget, -graz*m)
         _ADD_SOURCE_(self%id_morttarget, self%rmd*m+(1.0_rk-self%eff)*graz*m)
         _ADD_SOURCE_(self%id_exctarget, self%rmn*m*q)

      ! Leave spatial loops (if any)
      _LOOP_END_
   end subroutine do

end module hereon_zoo

!-----------------------------------------------------------------------
! Copyright OG-O - GNU Public License - www.gnu.org
!-----------------------------------------------------------------------