module real_subtraction_strata
  implicit none
  private

  integer,parameter,public :: real_stratum_regular=1
  integer,parameter,public :: real_stratum_migration=2
  integer,parameter,public :: n_real_strata=2

  public :: classify_real_subtraction_stratum
  public :: split_real_subtraction_weight
  public :: migration_pt_distance

contains

  pure integer function classify_real_subtraction_stratum(real_pass,alpha_active,mapped_pass) result(stratum)
    logical,intent(in) :: real_pass
    logical,intent(in) :: alpha_active(:),mapped_pass(:)

    stratum=real_stratum_regular
    if (size(alpha_active).ne.size(mapped_pass)) return
    if (any(alpha_active .and. (mapped_pass .neqv. real_pass))) stratum=real_stratum_migration
  end function classify_real_subtraction_stratum

  pure subroutine split_real_subtraction_weight(full_weight,stratum,regular_weight,migration_weight)
    real(kind=8),intent(in) :: full_weight
    integer,intent(in) :: stratum
    real(kind=8),intent(out) :: regular_weight,migration_weight

    regular_weight=0d0
    migration_weight=0d0
    if (stratum.eq.real_stratum_regular) then
       regular_weight=full_weight
    elseif (stratum.eq.real_stratum_migration) then
       migration_weight=full_weight
    endif
  end subroutine split_real_subtraction_weight

  pure real(kind=8) function migration_pt_distance(real_margin,mapped_margins,real_pass,&
       alpha_active,mapped_pass) result(distance)
    real(kind=8),intent(in) :: real_margin,mapped_margins(:)
    logical,intent(in) :: real_pass,alpha_active(:),mapped_pass(:)
    integer :: idip
    logical :: real_pt_pass,mapped_pt_pass

    distance=-1d0
    if (size(mapped_margins).ne.size(alpha_active) .or. &
         size(mapped_pass).ne.size(alpha_active)) return
    if (abs(real_margin).gt.0.25d0*huge(1d0)) return
    real_pt_pass=real_margin.gt.0d0
    do idip=1,size(alpha_active)
       if (.not.alpha_active(idip)) cycle
       if (mapped_pass(idip) .eqv. real_pass) cycle
       if (abs(mapped_margins(idip)).gt.0.25d0*huge(1d0)) cycle
       mapped_pt_pass=mapped_margins(idip).gt.0d0
       if (mapped_pt_pass .eqv. real_pt_pass) cycle
       if (distance.lt.0d0) then
          distance=min(abs(real_margin),abs(mapped_margins(idip)))
       else
          distance=min(distance,abs(real_margin),abs(mapped_margins(idip)))
       endif
    enddo
  end function migration_pt_distance

end module real_subtraction_strata
