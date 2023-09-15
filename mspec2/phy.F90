#include "fabm_driver.h"

!! phytoplankton for mixotrophic adaptation and optimality: a component of mspec

module hereon_phy
   use fabm_types

   implicit none

   private

   type, extends(type_base_model), public :: type_hereon_phy
      ! Variable identifiers
      type (type_state_variable_id) :: id_m, id_q
      type (type_state_variable_id) :: id_exctarget, id_morttarget, id_upttarget
      type (type_dependency_id)          :: id_par
      type (type_surface_dependency_id)  :: id_I_0
      type (type_diagnostic_variable_id) :: id_dPAR

      ! Model parameters
      real(rk) :: m0, tau, rmax, kc, alpha, vmax, clight,  nstar, rr, rmn, rmd, i_min
      logical :: ia
   contains
      procedure :: initialize
      procedure :: do
      !procedure :: do_ppdd
   end type

contains

   subroutine initialize(self, configunit)
      class (type_hereon_phy), intent(inout), target :: self
      integer,                        intent(in)            :: configunit

      real(rk), parameter :: d_per_s = 1.0_rk/86400.0_rk

      ! Store parameter values in our own derived type
      ! NB: all rates must be provided in values per day and are converted here to values per second.
      call self%get_parameter(self%m0,   'm0',   'mmol m-3',  'background concentration',      default=0.0225_rk)
      call self%get_parameter(self%tau,  'tau',  'non-dim',   'trade-off coefficient',        default=0.7_rk)
      call self%get_parameter(self%rmax, 'rmax', 'd-1',       'maximum specific growth rate', default=1.5_rk, scale_factor=d_per_s)
      call self%get_parameter(self%kc,   'kc',   'm2 mmol-1', 'specific light extinction',  default=0.03_rk)
      call self%get_parameter(self%alpha, 'alpha', 'mmol m-3',  'half-saturation nutrient concentration',   default=.6_rk)
      call self%get_parameter(self%vmax, 'vmax', 'd-1',  'maximum nutrient uptake rate',   default=1.8_rk, scale_factor=d_per_s)
      call self%get_parameter(self%clight,   'clight',   'PAR units', 'Light growth parameter',        default=0.04_rk)
      call self%get_parameter(self%nstar,  'nstar',  'non-dim',   'independence',        default=5.0_rk)
      call self%get_parameter(self%rr,  'rr',  'd-1',   'nutrient uptake respiration costs',        default=0.1_rk, scale_factor=d_per_s)
      call self%get_parameter(self%rmn,  'rmn',  'd-1',       'excretion rate',                default=0.01_rk, scale_factor=d_per_s)
      call self%get_parameter(self%rmd,  'rmd',  'd-1',       'mortality',                     default=0.02_rk, scale_factor=d_per_s)
      call self%get_parameter(self%i_min, 'i_min', 'W m-2',     'minimum light intensity in euphotic zone', default=25.0_rk)
      call self%get_parameter(self%ia, 'ia', '-', 'instantaneous acclimation', default=.false.)

      ! Register state variables
      call self%register_state_variable(self%id_m, 'c', 'mmol m-3', 'concentration', 0.0_rk, minimum=0.0_rk, vertical_movement=-0.0_rk)
      call self%register_state_variable(self%id_q, 'q', 'non-dim', 'nutrient quota', 0.5_rk, minimum=0.0_rk, maximum=1.0_rk)

      ! Register contribution of state to global aggregate variables.
      call self%add_to_aggregate_variable(standard_variables%total_nitrogen, self%id_m)

      ! Register diagnostic variables
      call self%register_diagnostic_variable(self%id_dPAR, 'PAR', 'W m-2',        'photosynthetically active radiation')

      ! Register environmental dependencies
      call self%register_dependency(self%id_par, standard_variables%downwelling_photosynthetic_radiative_flux)
      call self%register_dependency(self%id_I_0, standard_variables%surface_downwelling_photosynthetic_radiative_flux)

      ! Contribute to light attentuation
      call self%add_to_aggregate_variable(standard_variables%attenuation_coefficient_of_photosynthetic_radiative_flux, self%id_m, scale_factor=self%kc)
      call self%add_to_aggregate_variable(standard_variables%attenuation_coefficient_of_photosynthetic_radiative_flux, self%m0 * self%kc)

      ! Register dependencies on external state variables.
      call self%register_state_dependency(self%id_upttarget,  'uptake_target',    'mmol m-3', 'nutrient source')
      call self%register_state_dependency(self%id_exctarget,  'excretion_target', 'mmol m-3', 'sink for excreted matter')
      call self%register_state_dependency(self%id_morttarget, 'mortality_target', 'mmol m-3', 'sink for dead matter')
   end subroutine initialize

   subroutine do(self, _ARGUMENTS_DO_)
      class (type_hereon_phy), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      real(rk), parameter :: d_per_s = 1.0_rk/86400.0_rk

      real(rk) :: m, q, n, par, I_0
      real(rk) :: mu, nutuptake, graz
      real(rk) :: ci, cn, ca, gphot, ind

      ! Enter spatial loops (if any)
      _LOOP_BEGIN_

         ! Retrieve current (local) state variable values.
         _GET_(self%id_m,m)         ! mass
         _GET_(self%id_q,q)         ! nutritional state
         _GET_(self%id_upttarget,n) ! nutrients

         ! Retrieve current environmental conditions.
         _GET_(self%id_par,par)          ! local photosynthetically active radiation
         _GET_SURFACE_(self%id_I_0,I_0)  ! surface photosynthetically active radiation

         ! adaptation and growth rate calculation
         !!!   AQUI ME QUEDE
        
         !f = 0.1_rk+0.8_rk*f

         ci = 1.0_rk - exp(-(self%clight*par)**2)
         cn = n/(n+self%alpha)

         gphot = self%rmax !here is multiplied a temperature dependence term
             
         ind = self%nstar*(1.0_rk+q)
         ca = gfunc(ci,q,ind)
          
         mu = ca*gphot  - self%rr*self%vmax*(1.0_rk-q)*cn
         nutuptake = self%vmax*(1.0_rk-q)*cn 

         !mu = max(0.0_rk,mu)

         !!!
         ! Set temporal derivatives
         _ADD_SOURCE_(self%id_m, mu*m - self%rmd*m)
         _SET_ODE_(self%id_q, nutuptake-q*mu)
         _SET_DIAGNOSTIC_(self%id_dPAR,par)

         _ADD_SOURCE_(self%id_upttarget,-nutuptake*m)
         _ADD_SOURCE_(self%id_morttarget, self%rmd*m)
         _ADD_SOURCE_(self%id_exctarget, self%rmn*m*q)

      ! Leave spatial loops (if any)
      _LOOP_END_
   end subroutine do

   ! colimitation gfunction
   elemental real(rk) function gfunc(x,y,n)
      real(rk), intent(in) :: x,y,n
      real(rk) :: r
      gfunc = x * y
      !r = (x+0.001_rk)/(y+0.001_rk)
      !r = r*(1.0_rk-r**n)/(1.0_rk-r**(n+1.0_rk))
      !gfunc = y*r*(1.0_rk+y*x/n+log(4.0_rk**(-1.0_rk/n) + 1.0_rk/(2.0_rk*n)))
   end function gfunc

   ! derivative of the colimitation gfunction
   elemental real(rk) function devg(x,y,n)
      real(rk), intent(in) :: x,y,n
      devg = (gfunc(x,y+0.01_rk,n)-gfunc(x,y-0.01_rk,n))/(0.02_rk)
   end function devg

end module hereon_phy

!-----------------------------------------------------------------------
! Copyright OG-O - GNU Public License - www.gnu.org
!-----------------------------------------------------------------------