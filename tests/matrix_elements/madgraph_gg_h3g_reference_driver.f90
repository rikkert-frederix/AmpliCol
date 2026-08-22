! Reference-value driver for a MadGraph standalone `g g > h g g g` output.
!
! Compile this from SubProcesses/P1_gg_hggg, linking the generated matrix and
! model libraries.  The momentum point and alpha_s are the first event of the
! 100k-event sample documented in tests/cross_sections/heft_gg_h3g.
program madgraph_gg_h3g_reference_driver
  implicit none

  real(kind=8) :: p(0:3,6), value, summed, averaged
  integer :: hel(6), ic(6)
  integer :: h1, h2, h4, h5, h6
  real(kind=8), external :: matrix

  call setpara('../../Cards/param_card.dat')

  p(:,1)=[373.75312326566302d0, 0d0, 0d0, 373.75312326566302d0]
  p(:,2)=[403.47072902643276d0, 0d0, 0d0,-403.47072902643276d0]
  p(:,3)=[216.54829708118396d0,-42.341877309465147d0,&
       147.54305366710500d0, 87.785976723259580d0]
  p(:,4)=[200.03287849382329d0,-34.923666820043678d0,&
       -20.350832369065799d0,-195.90644092590080d0]
  p(:,5)=[220.58980288703404d0, 76.031960247160669d0,&
       -203.34116228757387d0,-39.132772432492416d0]
  p(:,6)=[140.05287383005444d0, 1.2335838823481566d0,&
       76.148940989534665d0,117.53563087436389d0]

  ic=1
  summed=0d0
  do h1=-1,1,2
     do h2=-1,1,2
        do h4=-1,1,2
           do h5=-1,1,2
              do h6=-1,1,2
                 hel=[h1,h2,0,h4,h5,h6]
                 value=matrix(p,hel,ic)
                 summed=summed+value
                 write (*,'(5(i3,1x),es26.17)') h1,h2,h4,h5,h6,value
              enddo
           enddo
        enddo
     enddo
  enddo

  call smatrix(p,averaged)
  write (*,'(a,es26.17)') 'MG_SUMMED_MATRIX2=',summed
  write (*,'(a,es26.17)') 'MG_AVERAGED_MATRIX2=',averaged
  write (*,'(a,es26.17)') 'MG_AVERAGED_TIMES_1536=',1536d0*averaged
end program madgraph_gg_h3g_reference_driver
