module hereon_model_library

   use fabm_types, only: type_base_model_factory, type_base_model
   
   use hereon_mspec
   use hereon_phy
   use hereon_zoo
   ! Add use statements for new models here

   implicit none

   private

   type, extends(type_base_model_factory) :: type_factory
   contains
      procedure :: create
   end type

   type (type_factory), save, target, public :: hereon_model_factory

contains

   subroutine create(self, name, model)

      class (type_factory), intent(in) :: self
      character(*),         intent(in) :: name
      class (type_base_model), pointer :: model

      select case (name)
         ! Add case statements for new models here
         case ('mspec'); allocate(type_hereon_mspec::model)
         case ('phy'); allocate(type_hereon_phy::model)
         case ('zoo'); allocate(type_hereon_zoo::model)
         case default
            call self%type_base_model_factory%create(name, model)
      end select

   end subroutine create

end module