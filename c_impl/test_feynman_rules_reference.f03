! test_feynman_rules_reference.f03
!
! Reference driver: calls every routine in the FeynmanRules module using the
! same inputs as the C++ test suite, then writes the results to the binary
! file "fortran_reference.dat".
!
! Each record in the file has the layout (no Fortran record markers –
! opened with access='stream'):
!
!   bytes  0-63  : tag string, null-padded to 64 bytes
!   bytes 64-67  : int32  n   (number of double-precision values that follow)
!   bytes 68-... : n * 8 bytes of real(kind=8) values
!                  complex(kind=8) values are stored as pairs (re, im)
!
! Build:
!   gfortran -O2 test_feynman_rules_reference.f03 feynmanrules.f03 \
!            -o fortran_ref
! Run:
!   ./fortran_ref
! Then run the C++ test binary to compare:
!   ./test_feynman_rules

program reference_driver
  use FeynmanRules
  implicit none

  ! ── reference momenta (same values as C++ pGeneric, pMassless, etc.) ──────
  real(kind=8), parameter :: pGeneric(0:3)   = [ 10d0,  3d0,  4d0,  5d0 ]
  real(kind=8), parameter :: pMassless(0:3)  = [ 10d0,  6d0,  0d0,  8d0 ]
  real(kind=8), parameter :: pAlongZ(0:3)    = [  5d0,  0d0,  0d0,  5d0 ]
  real(kind=8), parameter :: pAlongNZ(0:3)   = [  5d0,  0d0,  0d0, -5d0 ]
  real(kind=8), parameter :: pIncoming(0:3)  = [ -dsqrt(50d0),  3d0,  4d0,  5d0 ]
  real(kind=8), parameter :: pMass2(0:3)     = [  7d0,  1d0,  2d0,  3d0 ]

  ! ── working arrays ──────────────────────────────────────────────────────────
  complex(kind=8) :: wf1(6), wf2(6), wf3(6)
  complex(kind=8) :: wfq(6), wfaq(6)
  complex(kind=8) :: wfT(6)
  complex(kind=8) :: wfOut(6)
  complex(kind=8) :: wfs(1)
  complex(kind=8) :: wfsIn(1)
  real(kind=8)    :: wfR(4), wfR2(4) 
  real(kind=8)    :: wfTr(6)
  real(kind=8)    :: coupl(2)
  real(kind=8)    :: couplFlag(2)

  integer, parameter :: U = 10   ! file unit

  ! Many FR subroutines take a dimension(4) dummy wf while these locals are
  ! declared (6) (sized for tensor-gluon compatibility); components 5-6 are
  ! never written by such calls. Zero them up front so write_cmplx dumps a
  ! deterministic pad instead of leaking uninitialized stack memory, matching
  ! the C++ test's zero-initialized arrays.
  wf1 = (0d0,0d0); wf2 = (0d0,0d0); wf3 = (0d0,0d0)
  wfq = (0d0,0d0); wfaq = (0d0,0d0)
  wfT = (0d0,0d0); wfOut = (0d0,0d0)

  open(unit=U, file='fortran_reference.dat', &
       access='stream', form='unformatted', status='replace')

  ! ══════════════════════════════════════════════════════════════════════════
  ! ext_gluon_cmplx  – four momenta × two helicities
  ! ══════════════════════════════════════════════════════════════════════════

  call ext_gluon_cmplx(pMassless,  1, 0, wf1)
  call write_cmplx(U, 'ext_gluon_cmplx_massless_1', wf1, 6)

  call ext_gluon_cmplx(pMassless, -1, 0, wf1)
  call write_cmplx(U, 'ext_gluon_cmplx_massless_-1', wf1, 6)

  call ext_gluon_cmplx(pGeneric,   1, 0, wf1)
  call write_cmplx(U, 'ext_gluon_cmplx_generic_1', wf1, 6)

  call ext_gluon_cmplx(pGeneric,  -1, 0, wf1)
  call write_cmplx(U, 'ext_gluon_cmplx_generic_-1', wf1, 6)

  call ext_gluon_cmplx(pAlongZ,    1, 0, wf1)
  call write_cmplx(U, 'ext_gluon_cmplx_alongZ_1', wf1, 6)

  call ext_gluon_cmplx(pAlongZ,   -1, 0, wf1)
  call write_cmplx(U, 'ext_gluon_cmplx_alongZ_-1', wf1, 6)

  call ext_gluon_cmplx(pAlongNZ,   1, 0, wf1)
  call write_cmplx(U, 'ext_gluon_cmplx_alongNZ_1', wf1, 6)

  call ext_gluon_cmplx(pAlongNZ,  -1, 0, wf1)
  call write_cmplx(U, 'ext_gluon_cmplx_alongNZ_-1', wf1, 6)

  ! ══════════════════════════════════════════════════════════════════════════
  ! ext_gluon_real
  ! ══════════════════════════════════════════════════════════════════════════

  call ext_gluon_real(pMassless,  1, 0, wfR)
  call write_real(U, 'ext_gluon_real_1', wfR, 4)

  call ext_gluon_real(pMassless, -1, 0, wfR)
  call write_real(U, 'ext_gluon_real_-1', wfR, 4)

  ! ══════════════════════════════════════════════════════════════════════════
  ! ext_gluon_mass  – nhel = +1, -1, 0
  ! ══════════════════════════════════════════════════════════════════════════

  call ext_gluon_mass(pGeneric,  1,  1, wf1, 5d0)
  call write_cmplx(U, 'ext_gluon_mass_1', wf1, 6)

  call ext_gluon_mass(pGeneric, -1, -1, wf1, 5d0)
  call write_cmplx(U, 'ext_gluon_mass_m1', wf1, 6)

  call ext_gluon_mass(pGeneric,  0,  0, wf1, 5d0)
  call write_cmplx(U, 'ext_gluon_mass_0', wf1, 6)

  ! ══════════════════════════════════════════════════════════════════════════
  ! ext_quark  – outgoing (pGeneric) and incoming (pIncoming),
  !              massless (0.0) and massive (1.5)
  ! ══════════════════════════════════════════════════════════════════════════

  call ext_quark(pGeneric,   1, 0, wf1, 0d0)
  call write_cmplx(U, 'ext_quark_1_0', wf1, 6)

  call ext_quark(pGeneric,  -1, 0, wf1, 0d0)
  call write_cmplx(U, 'ext_quark_-1_0', wf1, 6)

  call ext_quark(pGeneric,   1, 0, wf1, 1.5d0)
  call write_cmplx(U, 'ext_quark_1_15', wf1, 6)

  call ext_quark(pGeneric,  -1, 0, wf1, 1.5d0)
  call write_cmplx(U, 'ext_quark_-1_15', wf1, 6)

  call ext_quark(pIncoming,  1, 0, wf1, 0d0)
  call write_cmplx(U, 'ext_quark_inc_1_0', wf1, 6)

  call ext_quark(pIncoming, -1, 0, wf1, 0d0)
  call write_cmplx(U, 'ext_quark_inc_-1_0', wf1, 6)

  call ext_quark(pIncoming,  1, 0, wf1, 1.5d0)
  call write_cmplx(U, 'ext_quark_inc_1_15', wf1, 6)

  call ext_quark(pIncoming, -1, 0, wf1, 1.5d0)
  call write_cmplx(U, 'ext_quark_inc_-1_15', wf1, 6)

  ! ══════════════════════════════════════════════════════════════════════════
  ! ext_antiquark
  ! ══════════════════════════════════════════════════════════════════════════

  call ext_antiquark(pGeneric,   1, 0, wf1, 0d0)
  call write_cmplx(U, 'ext_antiquark_1_0', wf1, 6)

  call ext_antiquark(pGeneric,  -1, 0, wf1, 0d0)
  call write_cmplx(U, 'ext_antiquark_-1_0', wf1, 6)

  call ext_antiquark(pGeneric,   1, 0, wf1, 1.5d0)
  call write_cmplx(U, 'ext_antiquark_1_15', wf1, 6)

  call ext_antiquark(pGeneric,  -1, 0, wf1, 1.5d0)
  call write_cmplx(U, 'ext_antiquark_-1_15', wf1, 6)

  ! ══════════════════════════════════════════════════════════════════════════
  ! Vertex functions
  ! Replicate exactly the C++ test inputs:
  !   wf1 = ext_gluon_cmplx(pGeneric,  +1)
  !   wf2 = ext_gluon_cmplx(pMassless, -1)
  !   wf3 = ext_gluon_cmplx(pMass2,    +1)   (for tensor tests)
  !   wfq  = ext_quark(pGeneric, +1, fmass=0)
  !   wfaq = ext_antiquark(pMass2, -1, fmass=0)
  ! ══════════════════════════════════════════════════════════════════════════

  call ext_gluon_cmplx(pGeneric,  +1, 0, wf1)
  call ext_gluon_cmplx(pMassless, -1, 0, wf2)
  call ext_gluon_cmplx(pMass2,    +1, 0, wf3)
  call ext_quark(pGeneric,    +1, 0, wfq,  0d0)
  call ext_antiquark(pMass2,  -1, 0, wfaq, 0d0)
  coupl = [ 0.7d0, 0.3d0 ]

  ! ThreeGluon(wf1, pGeneric, wf2, pMassless, wfOut)
  call ThreeGluon(wf1, pGeneric, wf2, pMassless, wfOut)
  call write_cmplx(U, 'ThreeGluon', wfOut, 6)

  ! FourGluon(wf1, wf2, wf1, wfOut)
  call FourGluon(wf1, wf2, wf1, wfOut)
  call write_cmplx(U, 'FourGluon', wfOut, 6)

  ! TwoGluontoTensor(wf1, wf2, wfT)
  call TwoGluontoTensor(wf1, wf2, wfT)
  call write_cmplx(U, 'TwoGluontoTensor', wfT, 6)

  ! TensorGluontoGluon(wfT, wf3, wfOut)
  call TensorGluontoGluon(wfT, wf3, wfOut)
  call write_cmplx(U, 'TensorGluontoGluon', wfOut, 6)

  ! GluonTensortoGluon(wf3, wfT, wfOut)
  call GluonTensortoGluon(wf3, wfT, wfOut)
  call write_cmplx(U, 'GluonTensortoGluon', wfOut, 6)

  ! QuarkGluontoQuark(wfq, wf2, wfOut)
  call QuarkGluontoQuark(wfq, wf2, wfOut)
  call write_cmplx(U, 'QuarkGluontoQuark', wfOut, 6)

  ! GluonQuarktoQuark(wf2, wfq, wfOut)
  call GluonQuarktoQuark(wf2, wfq, wfOut)
  call write_cmplx(U, 'GluonQuarktoQuark', wfOut, 6)

  ! AquarkGluontoAquark(wfaq, wf2, wfOut)
  call AquarkGluontoAquark(wfaq, wf2, wfOut)
  call write_cmplx(U, 'AQuarkGluontoAQuark', wfOut, 6)

  ! GluonAquarktoAquark(wf2, wfaq, wfOut)
  call GluonAquarktoAquark(wf2, wfaq, wfOut)
  call write_cmplx(U, 'GluonAQuarktoAQuark', wfOut, 6)

  ! QuarkAquarktoGluon(wfq, wfaq, wfOut, coupl)
  call QuarkAquarktoGluon(wfq, wfaq, wfOut, coupl)
  call write_cmplx(U, 'QuarkAQuarktoGluon', wfOut, 6)

  ! AquarkQuarktoGluon(wfaq, wfq, wfOut)
  call AquarkQuarktoGluon(wfaq, wfq, wfOut)
  call write_cmplx(U, 'AQuarkQuarktoGluon', wfOut, 6)

  ! GluonGluontoScalar(wf1, wf2, wfs, coupl)
  call GluonGluontoScalar(wf1, wf2, wfs, coupl)
  call write_cmplx(U, 'GluonGluontoScalar', wfs, 1)

  ! ScalarGluontoGluon(wfs, wf2, wfOut, coupl)
  call ScalarGluontoGluon(wfs, wf2, wfOut, coupl)
  call write_cmplx(U, 'ScalarGluontoGluon', wfOut, 6)

  ! ══════════════════════════════════════════════════════════════════════════
  ! Lepton-pair vertices
  !   Same fermion inputs as the quark/anti-quark gluon vertices:
  !     wfq  = ext_quark(pGeneric, +1, fmass=0)
  !     wfaq = ext_antiquark(pMass2, -1, fmass=0)
  !   coupl = [0.7, 0.3]  (both L and R entries are exercised)
  ! ══════════════════════════════════════════════════════════════════════════

  ! LeptonAleptontoGluon(wfq, wfaq, wfOut, coupl)
  wfOut = (0d0, 0d0)
  call LeptonAleptontoGluon(wfq, wfaq, wfOut, coupl)
  call write_cmplx(U, 'LeptonALeptontoGluon', wfOut, 6)

  ! AleptonLeptontoGluon(wfaq, wfq, wfOut, coupl)
  wfOut = (0d0, 0d0)
  call AleptonLeptontoGluon(wfaq, wfq, wfOut, coupl)
  call write_cmplx(U, 'ALeptonLeptontoGluon', wfOut, 6)

  ! ══════════════════════════════════════════════════════════════════════════
  ! Scalar vertices
  !   wfs holds the GluonGluontoScalar output from above; build a fresh,
  !   well-defined scalar input wfsIn = (1,0) via ext_scalar for the
  !   scalar-fed vertices.  ScalarScalartoScalar is exercised twice to cover
  !   both branches of the coupl(2) flag.
  ! ══════════════════════════════════════════════════════════════════════════

  call ext_scalar(pGeneric, 0, wfsIn)

  ! QuarkScalartoQuark(wfq, wfsIn, wfOut, coupl)
  wfOut = (0d0, 0d0)
  call QuarkScalartoQuark(wfq, wfsIn, wfOut, coupl)
  call write_cmplx(U, 'QuarkScalartoQuark', wfOut, 6)

  ! GluonScalartoGluon(wf1, wfsIn, wfOut, coupl)
  wfOut = (0d0, 0d0)
  call GluonScalartoGluon(wf1, wfsIn, wfOut, coupl)
  call write_cmplx(U, 'GluonScalartoGluon', wfOut, 6)

  ! ScalarScalartoScalar(wfsIn, wfsIn, wfs, coupl)  – default branch
  call ScalarScalartoScalar(wfsIn, wfsIn, wfs, coupl)
  call write_cmplx(U, 'ScalarScalartoScalar', wfs, 1)

  ! ScalarScalartoScalar with coupl(2) = -10  – extra factor i branch
  couplFlag = [ 0.7d0, -10d0 ]
  call ScalarScalartoScalar(wfsIn, wfsIn, wfs, couplFlag)
  call write_cmplx(U, 'ScalarScalartoScalar_flag', wfs, 1)

  ! ══════════════════════════════════════════════════════════════════════════
  ! Propagators
  ! ══════════════════════════════════════════════════════════════════════════

  ! GluonPropagator: start from ext_gluon_cmplx(pGeneric,+1) then apply
  call ext_gluon_cmplx(pGeneric, +1, 0, wf1)
  call GluonPropagator(wf1, pGeneric)
  call write_cmplx(U, 'GluonPropagator', wf1, 6)

  ! GluonPropagator_real: fixed input {1, 0.5, -0.3, 0.8}
  wfR = [ 1d0, 0.5d0, -0.3d0, 0.8d0 ]
  call GluonPropagator_real(wfR, pGeneric)
  call write_real(U, 'GluonPropagator_real', wfR, 4)

  ! GluonPropagator_mass: fixed input wfg = {1, 0.5+0.3i, -0.2+0.1i, 0.4}
  wf1 = [ (1d0,0d0), (0.5d0,0.3d0), (-0.2d0,0.1d0), (0.4d0,0d0), (0d0,0d0), (0d0,0d0) ]
  call GluonPropagator_mass(wf1, pGeneric, 91.2d0, 2.5d0)
  call write_cmplx(U, 'GluonPropagator_mass', wf1, 6)

  ! QuarkPropagator: start from ext_quark(pGeneric,+1,fmass=0)
  call ext_quark(pGeneric, +1, 0, wf1, 0d0)
  call QuarkPropagator(wf1, pMass2, 1.5d0, 0.1d0)
  call write_cmplx(U, 'QuarkPropagator', wf1, 6)

  ! AquarkPropagator: start from ext_antiquark(pGeneric,+1,fmass=0)
  call ext_antiquark(pGeneric, +1, 0, wf1, 0d0)
  call AquarkPropagator(wf1, pMass2, 1.5d0, 0.1d0)
  call write_cmplx(U, 'AQuarkPropagator', wf1, 6)

  ! ScalarPropagator: wfs = (1,0)
  wfs(1) = (1d0, 0d0)
  call ScalarPropagator(wfs, pGeneric, 125d0, 0.004d0)
  call write_cmplx(U, 'ScalarPropagator', wfs, 1)

  close(U)
  write(*,*) 'fortran_reference.dat written successfully.'

contains

  ! ── write a complex(kind=8) array ─────────────────────────────────────────
  ! Layout: 64-byte name | int32 count=2*n | 2*n doubles (re0,im0,re1,im1,...)
  subroutine write_cmplx(unit, tag, arr, n)
    integer,           intent(in) :: unit, n
    character(len=*),  intent(in) :: tag
    complex(kind=8),   intent(in) :: arr(n)
    character(len=64)              :: name64
    integer(kind=4)                :: cnt
    real(kind=8)                   :: buf(2*n)
    integer                        :: i
    name64 = ' '
    name64(1:len_trim(tag)) = tag
    cnt = int(2*n, kind=4)
    do i = 1, n
      buf(2*i-1) = real(arr(i),  kind=8)
      buf(2*i)   = aimag(arr(i))
    end do
    write(unit) name64, cnt, buf
  end subroutine write_cmplx

  ! ── write a real(kind=8) array ────────────────────────────────────────────
  ! Layout: 64-byte name | int32 count=n | n doubles
  subroutine write_real(unit, tag, arr, n)
    integer,          intent(in) :: unit, n
    character(len=*), intent(in) :: tag
    real(kind=8),     intent(in) :: arr(n)
    character(len=64)             :: name64
    integer(kind=4)               :: cnt
    name64 = ' '
    name64(1:len_trim(tag)) = tag
    cnt = int(n, kind=4)
    write(unit) name64, cnt, arr(1:n)
  end subroutine write_real

end program reference_driver
