module seleccio_de_potencial
      Character (Len=80) :: selec_gs='Ba_plus_gs'
      Real      (Kind=8) :: r_cutoff_gs=1.05d0
      Real      (Kind=8) :: umax_gs=6.836177d3
      Character (Len=80) :: selec_pi='Ba_plus_2D1'
      Real      (Kind=8) :: r_cutoff_pi=1.05d0
      Real      (Kind=8) :: umax_pi=5.216077d4
      Character (Len=80) :: selec_sigma='Ba_plus_2D0'
      Real      (Kind=8) :: r_cutoff_sigma=1.05d0
      Real      (Kind=8) :: umax_sigma=9.674836d4
      Character (Len=80) :: selec_delta='Ba_plus_2D2'
      Real      (Kind=8) :: r_cutoff_delta=1.05d0
      Real      (Kind=8) :: umax_delta=6.260015d4
end module seleccio_de_potencial

Module Modificacio_De_Select_Pot

Integer*4 Npg_Ba_plus_gs_fixC4, Npg_Ba_plus_pi_fixC4, Npg_Ba_plus_sigma_fixC4

Parameter(Npg_Ba_plus_gs_fixC4=13)
Parameter(Npg_Ba_plus_pi_fixC4=13)
Parameter(Npg_Ba_plus_sigma_fixC4=13)

Real*16 Pg_Ba_plus_gs_fixC4
Real*16 Pg_Ba_plus_pi_fixC4
Real*16 Pg_Ba_plus_sigma_fixC4

Integer*4 k0_Ba_plus_gs_fixC4
Integer*4 k0_Ba_plus_pi_fixC4
Integer*4 k0_Ba_plus_sigma_fixC4

Real*8 rcutoff_Ba_plus_gs_fixC4, fcutoff_Ba_plus_gs_fixC4,r0_Ba_plus_gs_fixC4,aa_Ba_plus_gs_fixC4,  &
       bb_Ba_plus_gs_fixC4,cc_Ba_plus_gs_fixC4
Real*8 rcutoff_Ba_plus_pi_fixC4, fcutoff_Ba_plus_pi_fixC4,r0_Ba_plus_pi_fixC4,aa_Ba_plus_pi_fixC4,  &
       bb_Ba_plus_pi_fixC4,cc_Ba_plus_pi_fixC4
Real*8 rcutoff_Ba_plus_sigma_fixC4, fcutoff_Ba_plus_sigma_fixC4,r0_Ba_plus_sigma_fixC4,aa_Ba_plus_sigma_fixC4, &
       bb_Ba_plus_sigma_fixC4,cc_Ba_plus_sigma_fixC4
Common/Param_LJ_V_Ba_plus_gs_fixC4/rcutoff_Ba_plus_gs_fixC4,fcutoff_Ba_plus_gs_fixC4,               &
       r0_Ba_plus_gs_fixC4,aa_Ba_plus_gs_fixC4,bb_Ba_plus_gs_fixC4,cc_Ba_plus_gs_fixC4,             &
       Pg_Ba_plus_gs_fixC4(Npg_Ba_plus_gs_fixC4),k0_Ba_plus_gs_fixC4
Common/Param_LJ_V_Ba_plus_pi_fixC4/rcutoff_Ba_plus_pi_fixC4,fcutoff_Ba_plus_pi_fixC4,               &
       r0_Ba_plus_pi_fixC4,aa_Ba_plus_pi_fixC4,bb_Ba_plus_pi_fixC4,cc_Ba_plus_pi_fixC4,             &
       Pg_Ba_plus_pi_fixC4(Npg_Ba_plus_pi_fixC4),k0_Ba_plus_pi_fixC4
Common/Param_LJ_V_Ba_plus_sigma_fixC4/rcutoff_Ba_plus_sigma_fixC4,fcutoff_Ba_plus_sigma_fixC4,      &
       r0_Ba_plus_sigma_fixC4,aa_Ba_plus_sigma_fixC4,bb_Ba_plus_sigma_fixC4,cc_Ba_plus_sigma_fixC4, &
       Pg_Ba_plus_sigma_fixC4(Npg_Ba_plus_sigma_fixC4),k0_Ba_plus_sigma_fixC4

  Real (Kind=8) :: Quita_C4_Ba_plus_gs_fix_C4=0.0d0
  Real (Kind=8) :: Quita_C4_Ba_plus_pi_fix_C4=0.0d0
  Real (Kind=8) :: Quita_C4_Ba_plus_sigma_fix_C4=0.0d0

  Save

End Module Modificacio_De_Select_Pot

Double Precision Function  V_gs(r)
  use seleccio_de_potencial
      Implicit none
      Real      (Kind=8) :: Select_pot
      Real      (Kind=8) :: r
      V_gs=Select_pot(selec_gs,r,r_cutoff_gs,umax_gs)
  End Function V_gs
  Double Precision Function V_pi(r)
  use seleccio_de_potencial
      Real      (Kind=8) :: Select_pot
      Real      (Kind=8) :: r
      V_pi=Select_pot(selec_pi,r,r_cutoff_pi,umax_pi)
  End Function V_pi
  Double Precision Function V_sigma(r)
  use seleccio_de_potencial
      Implicit none
      Real      (Kind=8) :: Select_pot
      Real      (Kind=8) :: r
      V_sigma=Select_pot(selec_sigma,r,r_cutoff_sigma,umax_sigma)
  End Function V_sigma
  Double Precision Function V_delta(r)
  use seleccio_de_potencial
      Implicit none
      Real      (Kind=8) :: Select_pot
      Real      (Kind=8) :: r
      V_delta=Select_pot(selec_delta,r,r_cutoff_delta,umax_delta)
  End Function V_delta
  Double Precision Function V_Sigma_Pi(r)
  use seleccio_de_potencial
      Implicit none
      Real      (Kind=8) :: r, V_Sigma, V_Pi
      V_Sigma_Pi=(V_Sigma(r) - V_Pi(r))/r**2
  End Function V_Sigma_Pi

Double Precision Function Select_Pot(selec,r,r_cutoff,umax)

Implicit None

Logical       :: Lcontrol=.false.        ! Per verificar el funcionament correcte dels nous potencials
Character (Len=80) :: selec
Real (Kind=8) :: r, r_cutoff, umax
Real (Kind=8) :: V_null=0.d0             ! selec='null'    Potencial de interacci� nul
Real (Kind=8) :: V_Ar_Ar                 ! selec='Ar_Ar'          Potencial Xe-He Tang & Toennies J.Chem. Phys. 118 (2003)  4976
Real (Kind=8) :: V_Xe_Xe                 ! selec='Xe_Xe'          Potencial Xe-He Tang & Toennies J.Chem. Phys. 118 (2003)  4976
Real (Kind=8) :: drV_Ar_Ar               ! selec='dr_Ar_Ar'       Potencial Xe-He Tang & Toennies J.Chem. Phys. 118 (2003)  4976
Real (Kind=8) :: drV_Xe_Xe               ! selec='dr_Xe_Xe'       Potencial Xe-He Tang & Toennies J.Chem. Phys. 118 (2003)  4976
Real (Kind=8) :: V_Ar_He                 ! selec='Ar_He'          Potencial Xe-He Tang & Toennies J.Chem. Phys. 118 (2003) 4976.
Real (Kind=8) :: V_Sr_He                 ! selec='Sr_He'          Potencial Sr-He Lovallo et al. Jcp 120, (2004) 246, Version SR
Real (Kind=8) :: V_LJ_OT                 ! selec='LJ_OT'   Potencial de Lenard-Jones capat a la Orsay-Trento
Real (Kind=8) :: V_Aziz_He               ! selec='Aziz_He' Potencial de Aziz pel He
Real (Kind=8) :: V_alka                  ! selec='LI', 'NA', 'K', 'RB' o 'CS' Potencials de Patil
Real (Kind=8) :: V_Na_Sigma              ! selec='Na_Sigma'          Pascale potential for Na(3Sigma)-He
Real (Kind=8) :: V_Na_Pi                 ! selec='Na_Pi'             Pascale potential for Na(3Pi)-He
Real (Kind=8) :: V_Kerkines_Li_gs        ! selec='Kerkines_Li_gs'    Kerkines-Mavidris Potential for Li-He gs
Real (Kind=8) :: V_Kerkines_Li_Pi        ! selec='Kerkines_Li_Pi'    Kerkines-Mavidris Potential for Li-He 2Pi
Real (Kind=8) :: V_Kerkines_Li_Sigma     ! selec='Kerkines_Li_Sigma' Kerkines-Mavidris Potential for Li-He 2Sigma
Real (Kind=8) :: V_He2s_fci              ! selec='He2S_fci', Potencial He(2S) del Eloranta (fci) (5-11-2015), Nist: He(2S): 329179.7623 cm-1
Real (Kind=8) :: V_Eloranta_He2_plus_gs  ! selec='He_plus', Potencial He2+ de Eloranta (6-6-2016)
Real (Kind=8) :: V_Koutselos_Cs_plus_gs  ! selec='Cs_plus_Koutselos'
Real (Kind=8) :: V_Koutselos_K_plus_gs   ! selec='K_plus_Koutselos'
Real (Kind=8) :: V_Koutselos_Rb_plus_gs  ! selec='Rb_plus_Koutselos'
Real (Kind=8) :: V_Koutselos_Na_plus_gs  ! selec='Na_plus_Koutselos'
Real (Kind=8) :: V_Fausto_Rb_plus_gs     ! selec='Rb_plus_Fausto'    New: computed by Fausto
Real (Kind=8) :: V_Ca                    ! selec='Ca'          Mauschick and Meyer 1989
Real (Kind=8) :: V_Ba_gs                 ! selec='Ba_gs'          Lovallo potential
Real (Kind=8) :: V_Ba_plus_old_gs        ! selec='Ba_plus_old_gs'
Real (Kind=8) :: V_Ba_plus_gs_fixC4      ! selec='Ba_plus_gs_fix_C4'
Real (Kind=8) :: V_Ba_plus_gs            ! selec='Ba_plus_gs'
Real (Kind=8) :: V_Ba_plus_pi            ! selec='Ba_plus_pi'
Real (Kind=8) :: V_Ba_plus_pi_fixC4      ! selec='Ba_plus_pi_fix_C4'
Real (Kind=8) :: V_Ba_plus_sigma         ! selec='Ba_plus_sigma'
Real (Kind=8) :: V_Ba_plus_sigma_fixC4   ! selec='Ba_plus_sigma_fix_C4'
Real (Kind=8) :: V_Ba_plus_2D0           ! selec='Ba_plus_2D0'
Real (Kind=8) :: V_Ba_plus_2D1           ! selec='Ba_plus_2D1'
Real (Kind=8) :: V_Ba_plus_2D2           ! selec='Ba_plus_2D2'
Real (Kind=8) :: V_CH3I_0                ! selec='CH3I_0'
Real (Kind=8) :: V_CH3I_1                ! selec='CH3I_1'
Real (Kind=8) :: V_CH3I_2                ! selec='CH3I_2'
Real (Kind=8) :: V_CH3I_3                ! selec='CH3I_3'
Real (Kind=8) :: V_CH3I_4                ! selec='CH3I_4'
Real (Kind=8) :: V_Ag_gs                 ! selec='Ag_gs'
Real (Kind=8) :: V_Ag_Pi                 ! selec='Ag_Pi'
Real (Kind=8) :: V_Ag_Sig                ! selec='Ag_Sig'
Real (Kind=8) :: V_K_4p_Sigma            ! selec='K_4p_Sigma'     Potencial de Pascale
Real (Kind=8) :: V_K_4p_Pi               ! selec='K_4p_Pi'            "     "     "       
Real (Kind=8) :: V_K_5s                  ! selec='K_5s'               "     "     "       
Real (Kind=8) :: V_Cs_7s                 ! selec='Cs_7s'              "     "     "
Real (Kind=8) :: V_Rb_6s                 ! selec='Rb_6s'              "     "     "
Real (Kind=8) :: V_Rb_6p_Sigma           ! selec='Rb_6p_sigma'        "     "     "
Real (Kind=8) :: V_Rb_6p_Pi              ! selec='Rb_6p_pi'           "     "     "
Real (Kind=8) :: V_Cs_sigma              ! selec='Cs_sigma'           "     "     "
Real (Kind=8) :: V_Cs_pi                 ! selec='Cs_pi'              "     "     "
Real (Kind=8) :: V_Fausto_Cs_gs          ! selec='Cs_gs_Fausto'       "     "  Fausto Cargnoni
Real (Kind=8) :: V_Fausto_Cs_6p0         ! selec='Cs_6p_Sigma_Fausto' "     "  Fausto Cargnoni
Real (Kind=8) :: V_Fausto_Cs_6p1         ! selec='Cs_6p_Pi_Fausto'    "     "  Fausto Cargnoni
Real (Kind=8) :: V_Fausto_Cs_7s          ! selec='Cs_7s_Fausto'       "     "  Fausto Cargnoni
Real (Kind=8) :: V_Rb_5p_sigma           ! selec='Rb_sigma'           "     "  Pascale
Real (Kind=8) :: V_Rb_5p_pi              ! selec='Rb_pi'              "     "      "
Real (Kind=8) :: V_Fausto_Rb_gs          ! selec='Rb_gs_Fausto'       "     "  Fausto Cargnoni
Real (Kind=8) :: V_Fausto_Rb_5p0         ! selec='Rb_5p_Sigma_Fausto' "     "  Fausto Cargnoni
Real (Kind=8) :: V_Fausto_Rb_5p1         ! selec='Rb_5p_Pi_Fausto'    "     "  Fausto Cargnoni
Real (Kind=8) :: V_Fausto_Rb_6s          ! selec='Rb_6s_Fausto'       "     "  Fausto Cargnoni
Real (Kind=8) :: V_grafeno                   ! selec='Grafeno'
Real (Kind=8) :: V_grafeno_sin_dispersion    ! selec='Grafeno_sin_dispersion'
Real (Kind=8) :: V_TiO2                  ! selec='Tio2'
Real (Kind=8) :: V_Au                    ! selec='Au'
Real (Kind=8) :: V_Au_TiO2               ! selec='Au_TiO2'                     potencial Au-TiO2( Pilar )
Real (Kind=8) :: V_dAu_TiO2              ! selec='dAu_TiO2'       Derivada del potencial Au-TiO2( Pilar )
Real (Kind=8) :: V_pH2_He                ! selec='pH2_He'         Potencial pH2-He ajutadode Gianturco
Real (Kind=8) :: V_Xe_He                 ! selec='Xe_He'          Potencial Xe-He Tang & Toennies J.Chem. Phys. 118 (2003)  4976
Character (Len=3)  :: elem
!
! Aqui comen�a la seleccio del potencial
!
If(Trim(selec).Eq.'Aziz_He')Then
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = V_Aziz_He(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_Aziz_He")')
   Endif
   Return
ElseIf(Trim(selec).Eq.'null')Then
   If(r.lt.r_cutoff)Then
     Select_pot = V_null
   Else
     Select_pot = V_null
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_null")')
   Endif
   Return
ElseIf(Trim(selec).Eq.'He2S_fci')Then
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = V_He2s_fci(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_He2s_fci")')
   Endif
   Return
ElseIf(Trim(selec).Eq.'He_plus')Then
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = V_Eloranta_He2_plus_gs(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_Eloranta_He2_plus_gs")')
   Endif
   Return
ElseIf(Trim(selec).Eq.'LJ_OT')Then
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = V_LJ_OT(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_LJ_OT")')
   Endif
   Return
ElseIf(Trim(selec).Eq.'LI'.Or.   &
       Trim(selec).Eq.'NA'.Or.   &
       Trim(selec).Eq.'K' .Or.   &
       Trim(selec).Eq.'RB'.Or.   &
       Trim(selec).Eq.'CS')Then
       If(r.lt.r_cutoff)Then
       Select_pot = umax
   Else
     If(Len_Trim(selec).Eq.1)elem=(Trim(selec)//'  ')
     If(Len_Trim(selec).Eq.2)elem=(Trim(selec)//' ')
     Select_pot = V_alka(r,elem)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_alka amb el parametre:",A3)')elem
   Endif
   Return
ElseIf(Trim(selec).Eq.'Cs_plus_Koutselos')Then
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = V_Koutselos_Cs_plus_gs(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_Koutselos_Cs_plus_gs")')
   Endif
   Return
ElseIf(Trim(selec).Eq.'Na_plus_Koutselos')Then
!
!  r_cutoff=1.0d0
!  umax=9.012825d4
!
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = V_Koutselos_Na_plus_gs(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_Koutselos_Na_plus_gs")')
   Endif
   Return
ElseIf(Trim(selec).Eq.'Rb_plus_Koutselos')Then
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = V_Koutselos_Rb_plus_gs(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_Koutselos_Rb_plus_gs")')
   Endif
   Return
ElseIf(Trim(selec).Eq.'Rb_plus_Fausto')Then
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = V_Fausto_Rb_plus_gs(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_Fausto_Rb_plus_gs")')
   Endif
   Return
ElseIf(Trim(selec).Eq.'K_plus_Koutselos')Then
!
!  r_cutoff=1.0d0, 2.0d0
!  umax=1.8628633412d5, 3.2762532479d3 
!
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = V_Koutselos_K_plus_gs(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_Koutselos_K_plus_gs")')
   Endif
   Return
ElseIf(Trim(selec).Eq.'Rb_plus_Fausto')Then
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = V_Fausto_Rb_plus_gs(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_Fausto_Rb_plus_gs")')
   Endif
   Return
ElseIf(Trim(selec).Eq.'Ba_gs')Then
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = V_Ba_gs(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_Ba_gs")')
   Endif
   Return
ElseIf(Trim(selec).Eq.'Ba_plus_old_gs')Then
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = V_Ba_plus_old_gs(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_Ba_plus_old_gs")')
   Endif
   Return
ElseIf(Trim(selec).Eq.'Ba_plus_gs')Then
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = V_Ba_plus_gs(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_Ba_plus_gs")')
   Endif
   Return
ElseIf(Trim(selec).Eq.'Ba_plus_gs_fix_C4')Then
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = V_Ba_plus_gs_fixC4(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_Ba_plus_gs_fixC4")')
   Endif
   Return
ElseIf(Trim(selec).Eq.'Ba_plus_pi')Then
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = V_Ba_plus_pi(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_Ba_plus_pi")')
   Endif
   Return
ElseIf(Trim(selec).Eq.'Ba_plus_pi_fix_C4')Then
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = V_Ba_plus_pi_fixC4(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_Ba_plus_pi_fixC4")')
   Endif
   Return
ElseIf(Trim(selec).Eq.'Ba_plus_sigma')Then
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = V_Ba_plus_sigma(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_Ba_plus_sigma")')
   Endif
   Return
ElseIf(Trim(selec).Eq.'Ba_plus_sigma_fix_C4')Then
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = V_Ba_plus_sigma_fixC4(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_Ba_plus_sigma_fixC4")')
   Endif
   Return
ElseIf(Trim(selec).Eq.'Ba_plus_2D0')Then
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = V_Ba_plus_2D0(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_Ba_plus_2D0")')
   Endif
   Return
ElseIf(Trim(selec).Eq.'Ba_plus_2D1')Then
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = V_Ba_plus_2D1(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_Ba_plus_2D1")')
   Endif
   Return
ElseIf(Trim(selec).Eq.'Ba_plus_2D2')Then
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = V_Ba_plus_2D2(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_Ba_plus_2D2")')
   Endif
   Return
ElseIf(Trim(selec).Eq.'Ag_gs')Then
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = V_Ag_gs(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_Ag_gs")')
   Endif
   Return
ElseIf(Trim(selec).Eq.'Ag_Pi')Then
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = V_Ag_Pi(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_Ag_Pi")')
   Endif
   Return
ElseIf(Trim(selec).Eq.'Ag_Sig')Then
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = V_Ag_Sig(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_Ag_Sig")')
   Endif
   Return
ElseIf(Trim(selec).Eq.'Cs_7s')Then
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = V_Cs_7s(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_Cs_7s")')
   Endif
   Return
ElseIf(Trim(selec).Eq.'Rb_6s')Then
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = V_Rb_6s(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_Rb_6s")')
   Endif
   Return
ElseIf(Trim(selec).Eq.'Cs_sigma')Then
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = V_Cs_sigma(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_Cs_sigma")')
   Endif
   Return
ElseIf(Trim(selec).Eq.'Cs_pi')Then
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = V_Cs_pi(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_Cs_pi")')
   Endif
   Return
ElseIf(Trim(selec).Eq.'Cs_gs_Fausto')Then
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = V_Fausto_Cs_gs(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_Fausto_Cs_gs")')
   Endif
   Return
ElseIf(Trim(selec).Eq.'Cs_6p_Sigma_Fausto')Then
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = V_Fausto_Cs_6p0(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_Fausto_Cs_6p0")')
   Endif
   Return
ElseIf(Trim(selec).Eq.'Cs_6p_Pi_Fausto')Then
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = V_Fausto_Cs_6p1(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_Fausto_Cs_6p1")')
   Endif
   Return
ElseIf(Trim(selec).Eq.'Cs_7s_Fausto')Then
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = V_Fausto_Cs_7s(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_Fausto_Cs_7s")')
   Endif
   Return
ElseIf(Trim(selec).Eq.'Rb_sigma')Then
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = V_Rb_5p_sigma(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_Rb_5p_sigma")')
   Endif
   Return
ElseIf(Trim(selec).Eq.'Rb_pi')Then
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = V_Rb_5p_pi(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_Rb_5p_pi")')
   Endif
   Return
ElseIf(Trim(selec).Eq.'Rb_6p_sigma')Then
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = V_Rb_6p_Sigma(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_Rb_6p_sigma")')
   Endif
   Return
ElseIf(Trim(selec).Eq.'Rb_6p_pi')Then
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = V_Rb_6p_Pi(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_Rb_6p_Pi")')
   Endif
   Return
ElseIf(Trim(selec).Eq.'Rb_gs_Fausto')Then
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = V_Fausto_Rb_gs(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_Fausto_Rb_gs")')
   Endif
   Return
ElseIf(Trim(selec).Eq.'Rb_5p_Sigma_Fausto')Then
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = V_Fausto_Rb_5p0(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_Fausto_Rb_5p0")')
   Endif
   Return
ElseIf(Trim(selec).Eq.'Rb_5p_Pi_Fausto')Then
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = V_Fausto_Rb_5p1(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_Fausto_Rb_5p1")')
   Endif
   Return
ElseIf(Trim(selec).Eq.'Rb_6s_Fausto')Then
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = V_Fausto_Rb_6s(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_Fausto_Rb_6s")')
   Endif
   Return
ElseIf(Trim(selec).Eq.'Grafeno')Then
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = V_grafeno(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_grafeno")')
   Endif
   Return
ElseIf(Trim(selec).Eq.'Grafeno_sin_dispersion')Then
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = V_grafeno_sin_dispersion(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_grafeno_sin_dispersion")')
   Endif
   Return
ElseIf(Trim(selec).Eq.'TiO2')Then
!
!  r_cutoff=2.5d0
!  umax=1.887964E+04
!
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = V_TiO2(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_Tio2")')
   Endif
   Return
ElseIf(Trim(selec).Eq.'Au')Then
!
!  r_cutoff=1.d0
!  umax=6.645232E+04
!
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = V_Au(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_TAu")')
   Endif
   Return
ElseIf(Trim(selec).Eq.'Au_TiO2')Then
!
!  r_cutoff=2.d0
!  umax=21450.4867067095d0
!
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = V_Au_TiO2(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_Au_TiO2")')
   Endif
   Return
ElseIf(Trim(selec).Eq.'dAu_TiO2')Then
!
!  r_cutoff=2.d0
!  umax=-6.2842194d4
!
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = V_dAu_TiO2(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_dAu_TiO2")')
   Endif
   Return
ElseIf(Trim(selec).Eq.'pH2_He')Then
!
!  r_cutoff=
!  umax=-6.2842194d4
!
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = V_pH2_He(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_pH2_He")')
   Endif
ElseIf(Trim(selec).Eq.'Xe_He')Then
!
!  r_cutoff= 2.0d0
!  umax= 21450.4867067095d0
!
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = V_Xe_He(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_Xe_He")')
   Endif
ElseIf(Trim(selec).Eq.'Ar_He')Then
!
!  r_cutoff= 1.9d0
!  umax= 1.226915474d4
!
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = V_Ar_He(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_Ar_He")')
   Endif
ElseIf(Trim(selec).Eq.'CH3I_0')Then
!
!  r_cutoff= 3.6d0
!  umax    = 3.1648655845d3
!
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = V_CH3I_0(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_CH3I_0(r)")')
   Endif
ElseIf(Trim(selec).Eq.'CH3I_1')Then
!
!  r_cutoff= 3.6d0
!  umax    = -5.1115736616d3
!
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = V_CH3I_1(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_CH3I_1(r)")')
   Endif
ElseIf(Trim(selec).Eq.'CH3I_2')Then
!
!  r_cutoff= 3.6d0
!  umax    = 5.5581502684d3
!
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = V_CH3I_2(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_CH3I_2(r)")')
   Endif
ElseIf(Trim(selec).Eq.'CH3I_3')Then
!
!  r_cutoff= 3.6d0
!  umax    = -2.7398540818d3
!
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = V_CH3I_3(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_CH3I_3(r)")')
   Endif
ElseIf(Trim(selec).Eq.'CH3I_4')Then
!
!  r_cutoff= 3.6d0
!  umax    = 2.4025334047d3
!
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = V_CH3I_4(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_CH3I_4(r)")')
   Endif
ElseIf(Trim(selec).Eq.'K_4p_Sigma')Then
!
!  r_cutoff= 2.2d0
!  umax    = 5.494421d3
!
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = V_K_4p_Sigma(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_K_4p_Sigma(r)")')
   Endif
ElseIf(Trim(selec).Eq.'K_4p_Pi')Then
!
!  r_cutoff= 1.9d0
!  umax    = 5.348304d3
!
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = V_K_4p_Pi(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_K_4p_Pi(r)")')
   Endif
ElseIf(Trim(selec).Eq.'K_5s')Then
!
!  r_cutoff= 2.d0
!  umax    = 2.729808d3
!
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = V_K_5s(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_K_5s(r)")')
   Endif
ElseIf(Trim(selec).Eq.'Ca')Then
!
!  r_cutoff= 2.5d0
!  umax    = 7.7661179647d3
!
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = V_Ca(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_Ca(r)")')
   Endif
ElseIf(Trim(selec).Eq.'Ar_Ar')Then
!
!  r_cutoff= 2.5d0
!  umax    = 7.7661179647d3
!
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = V_Ar_Ar(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_Ar_Ar(r)")')
   Endif
ElseIf(Trim(selec).Eq.'Xe_Xe')Then
!
!  r_cutoff= 2.5d0
!  umax    = 7.7661179647d3
!
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = V_Xe_Xe(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_Xe_Xe(r)")')
   Endif
ElseIf(Trim(selec).Eq.'dr_Xe_Xe')Then
!
!  r_cutoff= 2.5d0
!  umax    = 7.7661179647d3
!
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = drV_Xe_Xe(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a drV_Xe_Xe(r)")')
   Endif
ElseIf(Trim(selec).Eq.'dr_Ar_Ar')Then
!
!  r_cutoff= 2.5d0
!  umax    = 7.7661179647d3
!
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = drV_Ar_Ar(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a drV_Ar_Ar(r)")')
   Endif
ElseIf(Trim(selec).Eq.'Sr_He')Then
!
!  r_cutoff= 2.5d0
!  umax    = 6.423137d03
!
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = V_Sr_He(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_Sr_He(r)")')
   Endif
ElseIf(Trim(selec).Eq.'Kerkines_Li_gs')Then
!
!  r_cutoff= 1.d0
!  umax    = 8.962011d3
!
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = V_Kerkines_Li_gs(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_Kerkines_Li_gs(r)")')
   Endif
ElseIf(Trim(selec).Eq.'Kerkines_Li_Pi')Then
!
!  r_cutoff= 1.d0
!  umax    = 3.453693d4
!
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = V_Kerkines_Li_Pi(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_Kerkines_Li_Pi(r)")')
   Endif
ElseIf(Trim(selec).Eq.'Kerkines_Li_Sigma')Then
!
!  r_cutoff= 1.d0
!  umax    = 1.501675d5
!
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = V_Kerkines_Li_Sigma(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_Kerkines_Li_Sigma(r)")')
   Endif
ElseIf(Trim(selec).Eq.'Na_Sigma')Then
!
!  r_cutoff= 1.7d0
!  umax    = 1.0654129292E+04
!
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = V_Na_Sigma(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_Na_Sigma(r)")')
   Endif
ElseIf(Trim(selec).Eq.'Na_Pi')Then
!
!  r_cutoff= 1.65d0
!  umax    = 1.012716E+04
!
   If(r.lt.r_cutoff)Then
     Select_pot = umax
   Else
     Select_pot = V_Na_Pi(r)
   Endif
   If(Lcontrol)Then
     Write(6,'("selec....:",A80)')selec
     Write(6,'("Hem trucat a V_Na_Pi(r)")')
   Endif
Else
   Write(6,'(A,": aquest potencial no est� definit, de part de: Selec_Pot")')Trim(Selec)
   Stop 'From V_impur: 0001'
Endif
End Function Select_Pot
!
!  Potencial de Aziz pel He
!
double precision function V_Aziz_He(rr)
implicit none
real      (kind=8) :: Eps    = 10.948d0
real      (kind=8) :: A      = 1.8443101d5
real      (kind=8) :: alpha  = 10.43329537d0
real      (kind=8) :: beta   = -2.27965105d0
real      (kind=8) :: D      = 1.4826d0
real      (kind=8) :: C6     = 1.36745214d0
real      (kind=8) :: C8     = 0.42123807d0
real      (kind=8) :: C10    = 0.17473318d0
real      (kind=8) :: rm     = 2.963d0
real (kind=8) :: rr,r,ff
 ff=1.d0
 r=rr/rm
 if(r.le.D)ff=dexp(-(D/r-1.d0)**2)
 if(r==0)then
   V_Aziz_He= A*Eps
 else
   V_Aziz_He=(A*dexp( - alpha*r + beta*r**2) - ff*(C6/r**6 + C8/r**8 + C10/r**10 ))*Eps
 endif
end function

!
!  Potencial de Lenard-Jones capat a la Orsay-Trento
!
Block Data Inicio_LJ_V_LJ_OT
Implicit Real*8(A-H,O-Z)
Parameter(Npg=7)
!Real*16 Pg
Real*8 Pg
Integer*4 k0
Common/Param_LJ_V_LJ_OT/rcutoff,fcutoff,r0,aa,bb,cc,Pg(Npg),k0
!Data Pg/ -1.139924227407542Q4, 0., 0., 0., 0., 0., 3.178638072971337Q6/
Data Pg/ -1.139924227407542D4, 0., 0., 0., 0., 0., 3.178638072971337D6/
Data rcutoff/ 2.190323254008688d0/
Data fcutoff /1.574659520268255d2/
Data r0/ 0.d0/
Data aa/ 0.d0/
Data bb/ 0.d0/
Data cc/ 0.d0/
Data k0/   5/
End
Double Precision Function V_LJ_OT(x)
Parameter(Npg=7)
Implicit Real*8(A-H,O-Z)
!Real*16 Pg,Sto,xq
Real*8 Pg,Sto,xq
Integer*4 k0
Common/Param_LJ_V_LJ_OT/rcutoff,fcutoff,r0,aa,bb,cc,Pg(Npg),k0
If(x.le.rcutoff)Then
  V_LJ_OT=fcutoff
  Return
Endif
If(x.le.r0)Then
  V_LJ_OT=Abs(aa*x**2+bb*x+cc)
  Return
Endif
Sto=0.0q0
xq=1.0q0/x
Do k=1,Npg
  Sto=Sto+pg(k)*xq**(k+k0)
EndDo
V_LJ_OT=Sto
Return
End
!
!  Potencial de la Plata pel g.s.
!
double precision function V_Ag_gs(r)
implicit none
real      (kind=8) :: Ags    = 78741.3877
real      (kind=8) :: alphags= 0.65299
real      (kind=8) :: betags = 0.36367
real      (kind=8) :: Dgs    = 12.093
real      (kind=8) :: C6gs   = 137042.07d0
real      (kind=8) :: C8gs   = 189386.217d0
real      (kind=8) :: C10gs  = 111358462.d0
real      (kind=8) :: C12gs  = 1.17670377d10
real (kind=8) :: r,ff
 ff=1.d0
 if(r.le.Dgs)ff=dexp(-(Dgs/r-1.d0)**2)
 if(r==0)then
   V_Ag_gs= Ags
 else
   V_Ag_gs=Ags*dexp(-alphags*r-betags*r**2) - ff*(C6gs/r**6 + C8gs/r**8 + C10gs/r**10 + C12gs/r**12 )
 endif
end function

double precision function V_Ag_Pi(r)
implicit none
real      (kind=8) :: Api    = 0.7095092d6
real      (kind=8) :: alphapi= 1.9928
real      (kind=8) :: betapi = 0.45500
real      (kind=8) :: Dpi    = 4.0562
real      (kind=8) :: C6pi   = 0.21661134d6
real      (kind=8) :: C8pi   = 123454.7952
real      (kind=8) :: C10pi  = 0.307155024d6
real      (kind=8) :: C12pi  = 0.740464032d-10
real (kind=8) :: r,ff
 ff=1.d0
 if(r.le.Dpi)ff=dexp(-(Dpi/r-1.d0)**2)
 if(r==0)then
   V_Ag_Pi= Api
 else
   V_Ag_Pi=Api*dexp(-alphapi*r-betapi*r**2) - ff*(C6pi/r**6 + C8pi/r**8 + C10pi/r**10 + C12pi/r**12 )
 endif
end function

double precision function V_Ag_Sig(r)
implicit none
real      (kind=8) :: Asi    = 40005.2
real      (kind=8) :: alphasi= 0.41529
real      (kind=8) :: betasi = 0.14888
real      (kind=8) :: Dsi    = 25.730
real      (kind=8) :: C6si   = 1.97443d-9
real      (kind=8) :: C8si   = 3.37393d-6
real      (kind=8) :: C10si  = 0.00788363
real      (kind=8) :: C12si  = 7.52624d12
real (kind=8) :: r,ff
 ff=1.d0
 if(r.le.Dsi)ff=dexp(-(Dsi/r-1.d0)**2)
 if(r==0)then
   V_Ag_Sig= Asi
 else
   V_Ag_Sig=Asi*dexp(-alphasi*r-betasi*r**2) - ff*(C6si/r**6 + C8si/r**8 + C10si/r**10 + C12si/r**12 )
 endif
end function
double precision function V_alka(dist,elem)

implicit none

real (kind=8),    parameter :: evk     = 11604.448    ! 1eV   = 11604.448  K
real (kind=8),    parameter :: bohr    = 0.529177249  ! 1Bohr = 0.529 \AA
real (kind=8),    parameter :: hartree = 27.211608    ! 1H    = 27.2  eV
real (kind=8),    parameter :: factor  = evk*hartree

real      (kind=8), intent(IN) :: dist
character (len=3) , intent(IN) :: elem

real (kind=8) :: abar,u,v,roa
real (kind=8) :: a ,b ,c6 ,c8 ,agran
real (kind=8) :: av(5),bv(5),c6v(5),c8v(5),agranv(5)
real (kind=8) :: r2v(5),r4v(5),r6v(5)
real (kind=8) :: d,r2,r4,r6
real (kind=8) :: s
real (kind=8) :: vimp
real (kind=8) :: pi
real (kind=8) :: f(2)

character*3 celem(5)

integer (kind=4) ::  l,l7
integer (kind=4) ::  n,ielem

data celem /  'LI ' ,  'NA ' ,  'K  ' , 'RB '  , 'CS ' /
data av    / 1.588d0, 1.627d0, 1.771d0, 1.805d0, 1.869d0/
data bv    / 0.744d0, 0.744d0, 0.744d0, 0.744d0, 0.744d0/
data c6v   /  22.5d0,  24.7d0,  38.9d0,  44.6d0,  51.2d0/
data c8v   /  1.06d3,  1.29d3,  2.66d3,  3.18d3,  4.34d3/
data r2v   /  2.37d0,  2.37d0,  2.37d0,  2.37d0,  2.37d0/
data r4v   /  7.78d0,  7.78d0,  7.78d0,  7.78d0,  7.78d0/
data r6v   /  50.5d0,  50.5d0,  50.5d0,  50.5d0,  50.5d0/
data agranv/0.0519d0,0.0450d0,0.0259d0,0.0226d0,0.0173d0/


do ielem=1,7
  if(ielem.eq.7) STOP 'From V_impur: This impurity does not exist'
  if(elem.ne.celem(ielem)) cycle
  a     = av(ielem)
  b     = bv(ielem)
  c6    = c6v(ielem)
  c8    = c8v(ielem)
  r2    = r2v(ielem)
  r4    = r4v(ielem)
  r6    = r6v(ielem)
  agran = agranv(ielem)
  exit
end do

d=dist/bohr



pi    = 4.0d0*datan(1.0d0)
abar  = (a+b)/4.d0
u     = a-1.d0
v     = -0.5d0*a**2*u
roa   = agran*d**(2.d0*u)*(1.d0+v/d)**2*dexp(-2.d0*d/a)

do l=1,2
  s  = 0.d0
  l7 = 2*l+7
  do n=0,l7
    s = s+(d/abar)**n/fac(n)
  end do
  f(l)=1.d0-s*dexp(-d/abar)
end do
vimp = 4.d0*pi/3.d0*roa*(r2+37.d0/90.d0*r4/a**2 +            &
       28.d0/525.d0*r6/a**4)-c6*f(1)/d**6-c8*f(2)/d**8

V_alka=vimp*factor

return

contains

!....................................................................
!. Function fac 
!....................................................................

! Small version of factorial

double precision function fac(n)
implicit none
integer :: n,i
fac=1.d0
if (n.eq.0.or.n.eq.1) return
do i=2,n
  fac=fac*i
end do
return
end function fac
end function V_alka


! Everything in Angstrom Kelvin 

function V_Ba_gs(x)
!..........................................................!
! Ba-He potential. Fit with Mathematica in:
! /u/carraca/dmateo/barium/lovallo_Ba/Fit_potential.nb
!..........................................................!
implicit real*8(a-h,r-z)
    if(x<=3.7d0) then
        V_Ba_gs = 1297.1d0
    else
        V_Ba_gs = -6.146426949008034d14/x**18 + 1.6797251118695925d14/x**16 - 1.6513058293510484d13/x**14 + &
                7.015334514979032d11/x**12 - 1.1066880774828484d10/x**10 + 4.900282758445466d7/x**8 - &
                450294.6878721156d0/x**6
    endif
return
end
!------------------------------------------------------
!------------------------------------------------------
!------------------------------------------------------
!
! Potencial 2D Sigma del Ba+
!
Block Data Inicio_LJ_V_Ba_plus_2D0
Implicit Real*8(A-H,O-Z)
Parameter(Npg=13)
Real*16 Pg
Integer*4 k0
Common/Param_LJ_V_Ba_plus_2D0/rcutoff,fcutoff,r0,aa,bb,cc,Pg(Npg),k0
Data Pg/                           &
 -1087334.00470567619271105807407299Q0,  84828727.5548886783972206269121162Q0, -2706837566.34944393727530187448557Q0,  &
  47482174965.3049172156448636792347Q0, -519244658310.373320149026960630547Q0,  3776794752968.39098091679878465104Q0,  &
 -18931184609496.2743530917730707405Q0,  66467512849641.5867602022345173376Q0, -163349259737503.439380947298152502Q0,  &
  275471652871895.368050105366543138Q0, -303799513108731.163572027673779711Q0,  197342887711895.715873097734435815Q0,  &
 -57247315878299.6181404889603871296Q0/ 
Data rcutoff/ 2.0000000000000000D+00/
Data fcutoff/ 2.3097950000000001D+04/
Data r0/ 2.0000000000000000D+00/
Data aa/ 2.4148340053884458D+04/
Data bb/-1.4167500343195186D+05/
Data cc/ 2.0985459566015401D+05/
Data k0/   3/
End
Double Precision Function V_Ba_plus_2D0(x)
Parameter(Npg=13)
Implicit Real*8(A-H,O-Z)
Real*16 Pg,Sto,xq,Aux
Integer*4 k0
Common/Param_LJ_V_Ba_plus_2D0/rcutoff,fcutoff,r0,aa,bb,cc,Pg(Npg),k0
If(x.le.rcutoff)Then
  V_Ba_plus_2D0=fcutoff
  Return
Endif
If(x.le.r0)Then
  V_Ba_plus_2D0=Abs(aa*x**2+bb*x+cc)
  Return
Endif
Sto=0.0q0
xq=1.0q0/x
Aux=xq**(k0+1)
Do k=1,Npg
  Sto=Sto+pg(k)*Aux
  Aux=Aux*xq
EndDo
V_Ba_plus_2D0=Sto
Return
End
!------------------------------------------------------
!------------------------------------------------------
!------------------------------------------------------
Block Data Inicio_LJ_V_Ba_plus_2D1
Implicit Real*8(A-H,O-Z)
Parameter(Npg=13)
Real*16 Pg
Integer*4 k0
Common/Param_LJ_V_Ba_plus_2D1/rcutoff,fcutoff,r0,aa,bb,cc,Pg(Npg),k0
Data Pg/                           &
 -17806.3555611251407400764850972719Q0,  1823098.97448459689302902398930546Q0, -85830730.9944556694819916079148077Q0,  &
  1725318458.71278993648593521565983Q0, -19436563702.5972378924691700268134Q0,  135815448941.441168085385272965158Q0,  &
 -614470726947.401636604960682655298Q0,  1813794100230.25415271043202962245Q0, -3408845825677.80602502264479414645Q0,  &
  3767306350864.97856677416953497814Q0, -1870543050026.55759258430089080476Q0, -263830499030.803962487458419711050Q0,  &
  490468253633.822352576865392543355Q0/ 
Data rcutoff/ 2.0000000000000000D+00/
Data fcutoff/ 8.9667530000000006D+03/
Data r0/ 2.0000000000000000D+00/
Data aa/ 1.6993055762720283D+04/
Data bb/-9.4455519278364838D+04/
Data cc/ 1.2990556820533771D+05/
Data k0/   3/
End
Double Precision Function V_Ba_plus_2D1(x)
Parameter(Npg=13)
Implicit Real*8(A-H,O-Z)
Real*16 Pg,Sto,xq
Integer*4 k0
Common/Param_LJ_V_Ba_plus_2D1/rcutoff,fcutoff,r0,aa,bb,cc,Pg(Npg),k0
If(x.le.rcutoff)Then
  V_Ba_plus_2D1=fcutoff
  Return
Endif
If(x.le.r0)Then
  V_Ba_plus_2D1=Abs(aa*x**2+bb*x+cc)
  Return
Endif
Sto=0.0q0
xq=1.0q0/x
Do k=1,Npg
  Sto=Sto+pg(k)*xq**(k+k0)
EndDo
V_Ba_plus_2D1=Sto
Return
End

Block Data Inicio_LJ_V_Ba_plus_2D2
Implicit Real*8(A-H,O-Z)
Parameter(Npg=13)
Real*16 Pg
Integer*4 k0
Common/Param_LJ_V_Ba_plus_2D2/rcutoff,fcutoff,r0,aa,bb,cc,Pg(Npg),k0
Data Pg/                           &
  3438.08194066326485411722455335026Q0,  0.00000000000000000000000000000000Q0, -22271935.4558369081774443426658153Q0,  &
  509815799.148361882435180461889140Q0, -5103275255.84719311793286237036308Q0,  24347543963.7285486850837178167031Q0,  &
 -22098022739.0755793944143793174099Q0, -371505231363.903388585217865272476Q0,  2188591221460.46490291662922516704Q0,  &
 -6008459393013.55312031551791210851Q0,  9242835510174.84636285804976040549Q0, -7681295358440.52888750989831902749Q0,  &
  2696390604822.62910359945307774040Q0/ 
Data rcutoff/ 2.0000000000000000D+00/
Data fcutoff/ 1.1687330000000000D+04/
Data r0/ 2.0000000000000000D+00/
Data aa/ 1.9539148260089081D+04/
Data bb/-1.0968122164089163D+05/
Data cc/ 1.5289317408495379D+05/
Data k0/   3/
End
Double Precision Function V_Ba_plus_2D2(x)
Parameter(Npg=13)
Implicit Real*8(A-H,O-Z)
Real*16 Pg,Sto,xq
Integer*4 k0
Common/Param_LJ_V_Ba_plus_2D2/rcutoff,fcutoff,r0,aa,bb,cc,Pg(Npg),k0
If(x.le.rcutoff)Then
  V_Ba_plus_2D2=fcutoff
  Return
Endif
If(x.le.r0)Then
  V_Ba_plus_2D2=Abs(aa*x**2+bb*x+cc)
  Return
Endif
Sto=0.0q0
xq=1.0q0/x
Do k=1,Npg
  Sto=Sto+pg(k)*xq**(k+k0)
EndDo
V_Ba_plus_2D2=Sto
Return
End
!
! Nou ajust del Ba+-He amb mes punts a curtes distancies(4-06-15)
!
Block Data Inicio_LJ_V_Ba_plus_gs
Implicit Real*8(A-H,O-Z)
Parameter(Npg=13)
Real*16 Pg
Integer*4 k0
Common/Param_LJ_V_Ba_plus_gs/rcutoff,fcutoff,r0,aa,bb,cc,Pg(Npg),k0
Data Pg/                           &
 -1233025.18431181221429176595319949Q0,  84965045.0899243014427808833458580Q0, -2329830901.26373970562255233678155Q0,  &
  34487035482.8661381911727426332743Q0, -312315418217.722824241631231266832Q0,  1837177830634.86715823984693486481Q0,  &
 -7222749564186.92692954840041234856Q0,  19113448483679.6067074196915946358Q0, -33537885793985.8415962556220185604Q0,  &
  37188922413954.1348677597526953193Q0, -23134511849820.2719590286530591524Q0,  5477820267454.80225927552742559374Q0,  &
  657017642818.351074515977126633909Q0/ 
Data rcutoff/ 1.0000000000000000D+00/
Data fcutoff/ 2.7557869999999999D+04/
Data r0/ 2.0000000000000000D+00/
Data aa/-5.1918208117033464D+03/
Data bb/ 0.0000000000000000D+00/
Data cc/ 3.2749691986506852D+04/
Data k0/   3/
End
Double Precision Function V_Ba_plus_gs(x)
Parameter(Npg=13)
Implicit Real*8(A-H,O-Z)
Real*16 Pg,Sto,xq,Aux
Integer*4 k0
Common/Param_LJ_V_Ba_plus_gs/rcutoff,fcutoff,r0,aa,bb,cc,Pg(Npg),k0
If(x.le.rcutoff)Then
  V_Ba_plus_gs=fcutoff
  Return
Endif
If(x.le.r0)Then
  V_Ba_plus_gs=Abs(aa*x**2+bb*x+cc)
  Return
Endif
Sto=0.0q0
xq=1.0q0/x
Aux=xq**(k0+1)
Do k=1,Npg
  Sto=Sto+pg(k)*Aux
  Aux=Aux*xq
EndDo
V_Ba_plus_gs=Sto
Return
End

Block Data Inicio_LJ_V_Ba_plus_gs_fixC4
Implicit Real*8(A-H,O-Z)
Parameter(Npg=13)
Real*16 Pg
Integer*4 k0
Common/Param_LJ_V_Ba_plus_gs_fixC4/rcutoff,fcutoff,r0,aa,bb,cc,Pg(Npg),k0
Data Pg/                           &
 -17123.9045474000000000000000000000Q0,  1608302.51173741722899670327302374Q0, -79502305.7654667384046894884587665Q0,  &
  1903942512.84206481338934071885394Q0, -30134081435.0207804471353184171407Q0,  340887729307.024694401976841468399Q0,  &
 -2749060057173.10100779144500705256Q0,  15421639364129.7422511135464734090Q0, -58747373560398.8645333932507612952Q0,  &
  147892600745628.790351767793542405Q0, -234781354096552.082032333269632683Q0,  212676835684307.885943839636875957Q0,  &
 -83791034270656.6190293274527463146Q0/ 
Data rcutoff/ 1.0000000000000000D+00/
Data fcutoff/ 6.9240280000000002D+03/
Data r0/ 2.3999999999999999D+00/
Data aa/-8.5702212115060604D+02/
Data bb/ 0.0000000000000000D+00/
Data cc/ 7.7810505607043342D+03/
Data k0/   3/
End
Double Precision Function V_Ba_plus_gs_fixC4(x)
Parameter(Npg=13)
Implicit Real*8(A-H,O-Z)
Real*16 Pg,Sto,xq
Integer*4 k0
Common/Param_LJ_V_Ba_plus_gs_fixC4/rcutoff,fcutoff,r0,aa,bb,cc,Pg(Npg),k0
If(x.le.rcutoff)Then
  V_Ba_plus_gs_fixC4=fcutoff
  Return
Endif
If(x.le.r0)Then
  V_Ba_plus_gs_fixC4=Abs(aa*x**2+bb*x+cc)
  Return
Endif
Sto=0.0q0
xq=1.0q0/x
Do k=1,Npg
  Sto=Sto+pg(k)*xq**(k+k0)
EndDo
V_Ba_plus_gs_fixC4=Sto
Return
End
      Double Precision Function V_Ba_plus_old_gs(x)
      Implicit Real*8(A-H,O-Z)
      V_Ba_plus_old_gs=Func_Ba_plus_gs(x)
      Return
      End

      Block Data Inicio_LJ_Ba_plus_gs
      Implicit Real*8(A-H,O-Z)
      Parameter(Npg=16)
      Common/Param_LJ_Ba_plus_gs/r0,aa,bb,cc,Pg(Npg)
      Data Pg/                                                           &
      -5.75372193971343D+06, 7.15787621605274D+08,-4.96564765588082D+10, &
       1.98124418444111D+12,-4.71450937592845D+13, 6.52366641357938D+14, &
      -3.61296282269969D+15,-3.34341761554747D+16, 6.27946406049051D+17, &
      -1.25469318152625D+18,-2.47420248203689D+19,-1.56040616253124D+20, &
       4.94566594588086D+21,-2.19511834000164D+22,-6.01807193986940D+22, &
       4.58579988466807D+23                                              &
      /
      Data r0/ 3.000000000000000D+00/
      Data aa/ 1.406590885854125D+06/
      Data bb/-8.567381411545813D+06/
      Data cc/ 1.304908921806432D+07/
      End
      Double Precision Function Func_Ba_plus_gs(x)
      Parameter(Npg=16)
      Implicit Real*8(A-H,O-Z)
      Common/Param_LJ_Ba_plus_gs/r0,aa,bb,cc,Pg(Npg)
      If(x.le.r0)Then
        Func_Ba_plus_gs=aa*x**2+bb*x+cc
        Return
      Endif
      Sto=0.0d0
      Do k=1,Npg
        Sto=Sto+pg(k)/x**(2*(k-1)+6)
      EndDo
      Func_Ba_plus_gs=Sto
      Return
      End
Block Data Inicio_LJ_V_Ba_plus_pi
Implicit Real*8(A-H,O-Z)
Parameter(Npg=13)
Real*16 Pg
Integer*4 k0
Common/Param_LJ_V_Ba_plus_pi/rcutoff,fcutoff,r0,aa,bb,cc,Pg(Npg),k0
Data Pg/                           &
 -9826.50263754921899133897618912485Q0,  0.00000000000000000000000000000000Q0,  7905587.64166108193848904251702271Q0,  &
 -486504109.501964824574972292520630Q0,  10201523709.2410213874782756339362Q0, -115537603753.397220352368797609468Q0,  &
  811129915505.566723750323994437416Q0, -3741881843963.65413253722118668059Q0,  11556765251738.7569992644472524436Q0,  &
 -23666525048784.0598736494386530051Q0,  30832393059680.8435265222087894345Q0, -23137225457502.0957616538305554638Q0,  &
  7615700353700.85986113395012178657Q0/ 
Data rcutoff/ 2.0000000000000000D+00/
Data fcutoff/ 7.953992E+03/
Data r0/ 2.1250000000000000D+00/
Data aa/ 1.1631685135869398D+04/
Data bb/-6.7226852121137432D+04/
Data cc/ 9.5880955855918975D+04/
Data k0/   3/
End
Double Precision Function V_Ba_plus_pi(x)
Parameter(Npg=13)
Implicit Real*8(A-H,O-Z)
Real*16 Pg,Sto,xq
Integer*4 k0
Common/Param_LJ_V_Ba_plus_pi/rcutoff,fcutoff,r0,aa,bb,cc,Pg(Npg),k0
If(x.le.rcutoff)Then
  V_Ba_plus_pi=fcutoff
  Return
Endif
If(x.le.r0)Then
  V_Ba_plus_pi=Abs(aa*x**2+bb*x+cc)
  Return
Endif
Sto=0.0q0
xq=1.0q0/x
Do k=1,Npg
  Sto=Sto+pg(k)*xq**(k+k0)
EndDo
V_Ba_plus_pi=Sto
Return
End
!
! Comentari sobre la funcio
!
Block Data Inicio_LJ_V_Ba_plus_pi_fixC4
Implicit Real*8(A-H,O-Z)
Parameter(Npg=13)
Real*16 Pg
Integer*4 k0
Common/Param_LJ_V_Ba_plus_pi_fixC4/rcutoff,fcutoff,r0,aa,bb,cc,Pg(Npg),k0
Data Pg/                           &
 -17123.9045474000000000000000000000Q0,  1480585.38629277157051414278594678Q0, -66524410.6117297064088341969433940Q0,  &
  1248284301.42647870704024354840855Q0, -13069642195.0685028205930753960471Q0,  83381807453.7017144683188590909548Q0,  &
 -330934385581.979371931451499683451Q0,  772746327472.086796494677428176487Q0, -781686189760.313352024174457834920Q0,  &
 -739084882072.026661092648940188249Q0,  3162914410401.90196849006392816034Q0, -3570205775497.43736747130271625283Q0,  &
  1460004701781.63332377962082886275Q0/ 
Data rcutoff/ 2.0000000000000000D+00/
Data fcutoff/ 1.0377790000000001D+04/
Data r0/ 2.0000000000000000D+00/
Data aa/ 1.9282248948487551D+04/
Data bb/-1.0746305336734858D+05/
Data cc/ 1.4817490133941581D+05/
Data k0/   3/
End
Double Precision Function V_Ba_plus_pi_fixC4(x)
Parameter(Npg=13)
Implicit Real*8(A-H,O-Z)
Real*16 Pg,Sto,xq,Aux
Integer*4 k0
Common/Param_LJ_V_Ba_plus_pi_fixC4/rcutoff,fcutoff,r0,aa,bb,cc,Pg(Npg),k0
If(x.le.rcutoff)Then
  V_Ba_plus_pi_fixC4=fcutoff
  Return
Endif
If(x.le.r0)Then
  V_Ba_plus_pi_fixC4=Abs(aa*x**2+bb*x+cc)
  Return
Endif
Sto=0.0q0
xq=1.0q0/x
Aux=xq**(k0+1)
Do k=1,Npg
  Sto=Sto+pg(k)*Aux
  Aux=Aux*xq
EndDo
V_Ba_plus_pi_fixC4=Sto
Return
End
!
! Adjust of Ba+ 2P sigma (from Fausto calculations)
!
Block Data Inicio_LJ_V_Ba_plus_sigma
Implicit Real*8(A-H,O-Z)
Parameter(Npg=13)
Real*16 Pg
Integer*4 k0
Common/Param_LJ_V_Ba_plus_sigma/rcutoff,fcutoff,r0,aa,bb,cc,Pg(Npg),k0
Data Pg/                           &
 -1946800.56568338068309116827100120Q0,  160853137.584288414167976754320470Q0, -5443802098.75045627541231247577964Q0,  &
  101026984958.806322366443362695122Q0, -1156732968598.98543549109440671187Q0,  8665994393409.78419887735371608207Q0,  &
 -43862571678823.8062435100392360317Q0,  152374645410796.393876923253474426Q0, -363469581548844.057191969295126660Q0,  &
  584580735231448.559942865666177553Q0, -605026417604047.637639655063714044Q0,  363284554181176.601081977921602535Q0,  &
 -95975010650140.9085581326839809273Q0/ 
Data rcutoff/ 1.5000000000000000D+00/
Data fcutoff/ 2.2626640000000000D+07/
Data r0/ 1.0000000000000000D+00/
Data aa/ 1.5500784610027036D+04/
Data bb/-9.5332029461583734D+04/
Data cc/ 1.5198462687793490D+05/
Data k0/   3/
End
Double Precision Function V_Ba_plus_sigma(x)
Parameter(Npg=13)
Implicit Real*8(A-H,O-Z)
Real*16 Pg,Sto,xq,Aux
Integer*4 k0
Common/Param_LJ_V_Ba_plus_sigma/rcutoff,fcutoff,r0,aa,bb,cc,Pg(Npg),k0
If(x.le.rcutoff)Then
  V_Ba_plus_sigma=fcutoff
  Return
Endif
If(x.le.r0)Then
  V_Ba_plus_sigma=Abs(aa*x**2+bb*x+cc)
  Return
Endif
Sto=0.0q0
xq=1.0q0/x
Aux=xq**(k0+1)
Do k=1,Npg
  Sto=Sto+pg(k)*Aux
  Aux=Aux*xq
EndDo
V_Ba_plus_sigma=Sto
Return
End
!
! Comentari sobre la funcio
!
Block Data Inicio_LJ_V_Ba_plus_sigma_fixC4
Implicit Real*8(A-H,O-Z)
Parameter(Npg=13)
Real*16 Pg
Integer*4 k0
Common/Param_LJ_V_Ba_plus_sigma_fixC4/rcutoff,fcutoff,r0,aa,bb,cc,Pg(Npg),k0
Data Pg/                           &
 -17123.9045474000000000000000000000Q0,  17447561.0137612249712543003457326Q0, -1721357727.60473744645852747379753Q0,  &
  71971618157.3904453097848751064117Q0, -1705493726202.32479039775818481073Q0,  25631441812788.7446995165917324296Q0,  &
 -257596356972754.910546506340985159Q0,  1771762993803105.76996486810206154Q0, -8363689636005407.15424774634164701Q0,  &
  26640527998718191.1824700127470625Q0, -54721138485295325.7984174767437262Q0,  65456519239493414.9665950639016363Q0,  &
 -34651280867863688.9220625829424204Q0/ 
Data rcutoff/ 1.0000000000000000D+00/
Data fcutoff/ 2.0907770000000000D+04/
Data r0/ 3.5000000000000000D+00/
Data aa/ 1.4415808830195720D+03/
Data bb/-1.3729666852098217D+04/
Data cc/ 3.3158457362188419D+04/
Data k0/   3/
End
Double Precision Function V_Ba_plus_sigma_fixC4(x)
Parameter(Npg=13)
Implicit Real*8(A-H,O-Z)
Real*16 Pg,Sto,xq,Aux
Integer*4 k0
Common/Param_LJ_V_Ba_plus_sigma_fixC4/rcutoff,fcutoff,r0,aa,bb,cc,Pg(Npg),k0
If(x.le.rcutoff)Then
  V_Ba_plus_sigma_fixC4=fcutoff
  Return
Endif
If(x.le.r0)Then
  V_Ba_plus_sigma_fixC4=Abs(aa*x**2+bb*x+cc)
  Return
Endif
Sto=0.0q0
xq=1.0q0/x
Aux=xq**(k0+1)
Do k=1,Npg
  Sto=Sto+pg(k)*Aux
  Aux=Aux*xq
EndDo
V_Ba_plus_sigma_fixC4=Sto
Return
End
  Block Data Inicio_AP_Fit_APG_Cs_7s0                                                                  
      Implicit Real*8(A-H,O-Z)
      Parameter(Npg=12)
      Common/Param_AP_Fit_APG_Cs_7s0/r0,aa,bb,cc,Pg(Npg),Ps(Npg)
      Data Pg/								&
     -6.99807816089155D+07, 1.23439533336424D+04, 7.41947367149584D+06, &
      5.14036819939642D+07,-2.54666356840017D+13, 4.53268712277415D+09, &
      4.66829432137196D+12, 1.23144121825973D+14, 3.89994535976293D+12, &
      2.25670962111522D+16,-2.36027013243917D-04,-5.02039698946873D-04  &
     /
      Data Ps/								&
      1.29148499223708D+00, 7.42273770077802D-01, 1.83956945157573D+00, &
      3.01474078952469D+00, 6.17493897637321D+00, 4.97136805912161D+00, &
      6.59152034486467D+00, 8.41614866661145D+00, 7.36125612301946D+00, &
      1.34956539961155D+01, 1.18675333603722D+00, 1.59183920812119D+00  &
     /
      Data r0/ 2.100000000000000D+00/
      Data aa/ 2.420229505061327D+07/
      Data bb/-1.036963443170083D+08/
      Data cc/ 1.110675394653385D+08/
      End
      Double Precision Function V_Cs_7s(xx)
      Parameter(Npg=12)
      Implicit Real*8(A-H,O-Z)
      Common/Param_AP_Fit_APG_Cs_7s0/r0,aa,bb,cc,Pg(Npg),Ps(Npg)
      data x0/2.5d0/
      If(xx.lt.x0)then
        x=x0
      Else
         x=xx
      Endif   
      If(x.le.r0)Then
        V_Cs_7s=aa*x**2+bb*x+cc
        Return
      Endif
      Sto=pg(1)*Dexp(-Ps(1)*x)/x
      Do k=2,Npg
        Sto=Sto+pg(k)*x**(k-1)*Dexp(-Ps(k)*x)
      EndDo
      V_Cs_7s=Sto
      Return
      End

      Double Precision Function V_Koutselos_Cs_plus_gs(xx)
!
!      Potencial He-Cs+ de Koutselos et al, JCP 93, 7125 (1990)
!
!                     xx: Entrada en Angs.
! V_Koutselos_Cs_plus_gs: Salida en K
!
      Implicit none
      Real  (Kind=8)  :: aux,sto,x2,x4,x6,x8,CC4,CC6,CC8
      Real  (Kind=8)  :: xx,x,h_r,Rs,V_ex
!
!     Ctes para convertir de Angs. a bohr y de atomicas a K
!
      Real  (Kind=8)  :: a_bohr=0.529177d0, Hartree=27.2114d0, eV_K=11604.448d0
!
! Parametros adimensionales para los potenciales iones alkalinos-gas noble
!
      Real  (Kind=8)  :: A=146.98d0, a_a=1.5024d0, B=70.198d0, b_b=1.4041d0
!
! Los parametors estan en unidades atomicas
!

! Empiezan los parametros para el He-Cs+
      Real  (Kind=8)  :: v_0=0.5092d0, rho=0.9241d0, R_m=6.38d0
      Real  (Kind=8)  :: C4=0.0d0, C6=13.69d0, C8=236.4d0
      Real  (Kind=8)  :: a_d=1.3831d0, a_q=2.4434d0, a_o=10.614d0
! Acaban los parametros para el Cs+

      CC4 = C4 + a_d/2.0d0; CC6 = C6 + a_q/2.0d0; CC8 =C8 + a_o/2.0d0
      x=xx/a_bohr  ! Transformem els Angs en unitats at�miques de llongitut
      Rs=x/rho     
      V_ex=v_0*(A*Dexp(-a_a*Rs) - B*Dexp(-b_b*Rs))
      Aux=R_m*1.28d0
      If(x.Ge.Aux)Then
        h_r=1.0d0
      Else
        h_r=0.0d0
        If(x.Ne.0.0d0)h_r=Dexp(-(Aux/x - 1.d0)**2)
      Endif
      Sto=V_ex
      If(x.Ne.0.0d0)Then
        x2=x*x; x4=x2*x2; x6=x4*x2; x8=x6*x2
        Sto = Sto + (-CC4/x4 -CC6/x6 -CC8/x8)*h_r
      Endif
      V_Koutselos_Cs_plus_gs=Sto*Hartree*eV_k
      End Function V_Koutselos_Cs_plus_gs
      Double Precision Function V_Koutselos_Na_plus_gs(xx)
!
!                     xx: Entrada en Angs.
! V_Koutselos_Na_plus_gs: Salida en K
!
      Implicit none
      Real  (Kind=8)  :: aux,sto,x2,x4,x6,x8,CC4,CC6,CC8
      Real  (Kind=8)  :: xx,x,h_r,Rs,V_ex
!
!     Ctes para convertir de Angs. a bohr y de atomicas a K
!
      Real  (Kind=8)  :: a_bohr=0.529177d0, Hartree=27.2114d0, eV_K=11604.448d0
!
! Parametros adimensionales para los potenciales iones alkalinos-gas noble
!
      Real  (Kind=8)  :: A=146.98d0, a_a=1.5024d0, B=70.198d0, b_b=1.4041d0
!
! Los parametors estan en unidades atomicas
!

! Empiezan los parametros para el He-Na+
      Real  (Kind=8)  :: v_0=0.2063d0, rho=0.7640d0, R_m=4.55d0
      Real  (Kind=8)  :: C4=0.0d0, C6=1.424d0, C8=13.17d0
! Acaban los parametros para el Rb+
!
! Parametres propis del He
      Real  (Kind=8)  :: a_d=1.3831d0, a_q=2.4434d0, a_o=10.614d0

      CC4 = C4 + a_d/2.0d0; CC6 = C6 + a_q/2.0d0; CC8 =C8 + a_o/2.0d0
      x=xx/a_bohr  ! Transformem els Angs en unitats at�miques de llongitut
      Rs=x/rho     
      V_ex=v_0*(A*Dexp(-a_a*Rs) - B*Dexp(-b_b*Rs))
      Aux=R_m*1.28d0
      If(x.Ge.Aux)Then
        h_r=1.0d0
      Else
        h_r=0.0d0
        If(x.Ne.0.0d0)h_r=Dexp(-(Aux/x - 1.d0)**2)
      Endif
      Sto=V_ex
      If(x.Ne.0.0d0)Then
        x2=x*x; x4=x2*x2; x6=x4*x2; x8=x6*x2
        Sto = Sto + (-CC4/x4 -CC6/x6 -CC8/x8)*h_r
      Endif
      V_Koutselos_Na_plus_gs=Sto*Hartree*eV_k
      End Function V_Koutselos_Na_plus_gs
      Double Precision Function V_Koutselos_Rb_plus_gs(xx)
!
!                     xx: Entrada en Angs.
! V_Koutselos_Rb_plus_gs: Salida en K
!
      Implicit none
      Real  (Kind=8)  :: aux,sto,x2,x4,x6,x8,CC4,CC6,CC8
      Real  (Kind=8)  :: xx,x,h_r,Rs,V_ex
!
!     Ctes para convertir de Angs. a bohr y de atomicas a K
!
      Real  (Kind=8)  :: a_bohr=0.529177d0, Hartree=27.2114d0, eV_K=11604.448d0
!
! Parametros adimensionales para los potenciales iones alkalinos-gas noble
!
      Real  (Kind=8)  :: A=146.98d0, a_a=1.5024d0, B=70.198d0, b_b=1.4041d0
!
! Los parametors estan en unidades atomicas
!

! Empiezan los parametros para el He-Rb+
      Real  (Kind=8)  :: v_0=0.3474d0, rho=0.8868d0, R_m=5.85d0
      Real  (Kind=8)  :: C4=0.0d0, C6=8.8d0, C8=126.2d0
      Real  (Kind=8)  :: a_d=1.3831d0, a_q=2.4434d0, a_o=10.614d0
! Acaban los parametros para el Rb+

      CC4 = C4 + a_d/2.0d0; CC6 = C6 + a_q/2.0d0; CC8 =C8 + a_o/2.0d0
      x=xx/a_bohr  ! Transformem els Angs en unitats at�miques de llongitut
      Rs=x/rho     
      V_ex=v_0*(A*Dexp(-a_a*Rs) - B*Dexp(-b_b*Rs))
      Aux=R_m*1.28d0
      If(x.Ge.Aux)Then
        h_r=1.0d0
      Else
        h_r=0.0d0
        If(x.Ne.0.0d0)h_r=Dexp(-(Aux/x - 1.d0)**2)
      Endif
      Sto=V_ex
      If(x.Ne.0.0d0)Then
        x2=x*x; x4=x2*x2; x6=x4*x2; x8=x6*x2
        Sto = Sto + (-CC4/x4 -CC6/x6 -CC8/x8)*h_r
      Endif
      V_Koutselos_Rb_plus_gs=Sto*Hartree*eV_k
      End Function V_Koutselos_Rb_plus_gs
      Double Precision Function V_Koutselos_K_plus_gs(xx)
!
!                     xx: Entrada en Angs.
! V_Koutselos_Rb_plus_gs: Salida en K
!
      Implicit none
      Real  (Kind=8)  :: aux,sto,x2,x4,x6,x8,CC4,CC6,CC8
      Real  (Kind=8)  :: xx,x,h_r,Rs,V_ex
!
!     Ctes para convertir de Angs. a bohr y de atomicas a K
!
      Real  (Kind=8)  :: a_bohr=0.529177d0, Hartree=27.2114d0, eV_K=11604.448d0
!
! Parametros adimensionales para los potenciales iones alkalinos-gas noble
!
      Real  (Kind=8)  :: A=146.98d0, a_a=1.5024d0, B=70.198d0, b_b=1.4041d0
!
! Los parametors estan en unidades atomicas
!

! Empiezan los parametros para el He-Rb+
      Real  (Kind=8)  :: v_0=0.2690d0, rho=0.8589d0, R_m=5.49d0
      Real  (Kind=8)  :: C4=0.0d0, C6=5.83d0, C8=72.05d0
      Real  (Kind=8)  :: a_d=1.3831d0, a_q=2.4434d0, a_o=10.614d0
! Acaban los parametros para el Rb+

      CC4 = C4 + a_d/2.0d0; CC6 = C6 + a_q/2.0d0; CC8 =C8 + a_o/2.0d0
      x=xx/a_bohr  ! Transformem els Angs en unitats at�miques de llongitut
      Rs=x/rho     
      V_ex=v_0*(A*Dexp(-a_a*Rs) - B*Dexp(-b_b*Rs))
      Aux=R_m*1.28d0
      If(x.Ge.Aux)Then
        h_r=1.0d0
      Else
        h_r=0.0d0
        If(x.Ne.0.0d0)h_r=Dexp(-(Aux/x - 1.d0)**2)
      Endif
      Sto=V_ex
      If(x.Ne.0.0d0)Then
        x2=x*x; x4=x2*x2; x6=x4*x2; x8=x6*x2
        Sto = Sto + (-CC4/x4 -CC6/x6 -CC8/x8)*h_r
      Endif
      V_Koutselos_K_plus_gs=Sto*Hartree*eV_k
      End Function V_Koutselos_K_plus_gs
Block Data Inicio_LJ_Rb_6s
Implicit Real*8(A-H,O-Z)
Parameter(Npg=11)
Real*16 Pg
Common/Param_LJ_Rb_6s/rcutoff,fcutoff,r0,aa,bb,cc,Pg(Npg)
Data Pg/                           &
 -926758671.131109233668080149984260    ,  17751542326.9239124078481901178123    ,  2294371208218.10069268028221768684    ,  &
 -127580931840930.107534925325246956    ,  3010274903906904.53291694447191940    , -40628265880993902.4497466283119501    ,  &
  341766354107911552.266136623510917    , -1829660300635868428.84732815469319    ,  6079663336799721325.17593352865994    ,  &
 -11458268191609309405.5213761979227    ,  9376679057615973867.05134498742729    /
Data rcutoff/ 4.5000000000000000D+00/
Data fcutoff/ 8.2666990000000005D+03/
Data r0/ 4.5000000000000000D+00/
Data aa/ 8.8941971711569586D+03/
Data bb/-8.8209867027930217D+04/
Data cc/ 2.1983426802106350D+05/
End
Double Precision Function V_Rb_6s(x)
Parameter(Npg=11)
Implicit Real*8(A-H,O-Z)
Real*16 Pg,Sto
Common/Param_LJ_Rb_6s/rcutoff,fcutoff,r0,aa,bb,cc,Pg(Npg)
If(x.le.rcutoff)Then
  V_Rb_6s=fcutoff
  Return
Endif
If(x.le.r0)Then
  V_Rb_6s=Abs(aa*x**2+bb*x+cc)
  Return
Endif
Sto=0.0q0
Do k=1,Npg
  Sto=Sto+pg(k)/x**(k+5)
EndDo
V_Rb_6s=Sto
Return
End
Block Data Inicio_LJ_V_Cs_sigma
Implicit Real*8(A-H,O-Z)
Parameter(Npg=11)
Real*16 Pg
Common/Param_LJ_V_Cs_sigma/rcutoff,fcutoff,r0,aa,bb,cc,Pg(Npg)
Data Pg/                           &
  218972446.498381935851370116102201    , -20328460725.5889253542688904081392    ,  792064635008.609454121016549101098    ,  &
 -16988413450311.4536905549403585190    ,  220781047872801.503342285464821515    , -1818628789150376.23619182940382982    ,  &
  9689518147992622.13883260672871893    , -33234418435959134.0581023528704222    ,  70598735214769031.0197085111104892    ,  &
 -83973545896479856.9435608557993711    ,  42292199911291753.0362255745407528    /
Data rcutoff/ 3.9688300000000001D+00/
Data fcutoff/ 1.4295530000000001D+03/
Data r0/ 3.7999999999999998D+00/
Data aa/-5.2840987940870477D+01/
Data bb/ 0.0000000000000000D+00/
Data cc/ 2.2618840928664013D+03/
End
Double Precision Function V_Cs_sigma(x)
Parameter(Npg=11)
Implicit Real*8(A-H,O-Z)
Real*16 Pg,Sto
Common/Param_LJ_V_Cs_sigma/rcutoff,fcutoff,r0,aa,bb,cc,Pg(Npg)
If(x.le.rcutoff)Then
  V_Cs_sigma=fcutoff
  Return
Endif
If(x.le.r0)Then
  V_Cs_sigma=Abs(aa*x**2+bb*x+cc)
  Return
Endif
Sto=0.0q0
Do k=1,Npg
  Sto=Sto+pg(k)/x**(k+5)
EndDo
V_Cs_sigma=Sto
Return
End
Block Data Inicio_LJ_Cs_pi
Implicit Real*8(A-H,O-Z)
Parameter(Npg=11)
Real*16 Pg
Common/Param_LJ_Cs_pi/rcutoff,fcutoff,r0,aa,bb,cc,Pg(Npg)
Data Pg/                           &
 -271926.174454643323364441958407885    , -18233848.3121256520794825424693032    ,  953937957.545172931044438821385762    ,  &
 -27223560403.3593039055504509753318    ,  364166058669.852409015532291696960    , -2714770540249.03637313691657921121    ,  &
  12344772907874.9689365059244245234    , -35232833044851.8149649455556894484    ,  61848915415723.4235994491545114415    ,  &
 -61192851763498.3879914475300334623    ,  26164252463999.8193908959197343657    /
Data rcutoff/ 2.3813000000000000D+00/
Data fcutoff/ 2.5438229999999999D+03/
Data r0/ 2.3799999999999999D+00/
Data aa/-4.6845761002283029D+02/
Data bb/ 0.0000000000000000D+00/
Data cc/ 5.2002536629712176D+03/
End
Double Precision Function V_Cs_pi(x)
Parameter(Npg=11)
Implicit Real*8(A-H,O-Z)
Real*16 Pg,Sto
Common/Param_LJ_Cs_pi/rcutoff,fcutoff,r0,aa,bb,cc,Pg(Npg)
If(x.le.rcutoff)Then
  V_Cs_pi=fcutoff
  Return
Endif
If(x.le.r0)Then
  V_Cs_pi=Abs(aa*x**2+bb*x+cc)
  Return
Endif
Sto=0.0q0
Do k=1,Npg
  Sto=Sto+pg(k)/x**(k+5)
EndDo
V_Cs_pi=Sto
Return
End
!
! Function for He-Cs gs(6s) potential, from Fausto Cargnoni
!
Block Data Inicio_LJ_V_Fausto_Cs_gs
Implicit Real*8(A-H,O-Z)
Parameter(Npg=14)
Real*16 Pg
Integer*4 k0
Common/Param_LJ_V_Fausto_Cs_gs/rcutoff,fcutoff,r0,aa,bb,cc,Pg(Npg),k0
Data Pg/                           &
 -692107812.206136401684998300995763Q0,  53612999836.2373067596697398126653Q0, -1870657670551.24152231789523420133Q0,  &
  38854872489697.0224741123416213730Q0, -535664068036490.962168104702289226Q0,  5175823532151103.79719017070590004Q0,  &
 -36081508593491482.7030173645408564Q0,  183930078114942817.654215375874979Q0, -686481003395905469.369448400855558Q0,  &
  1855128368497788766.28473667982076Q0, -3532497504702993458.11254607147916Q0,  4493344993463040582.19605074943107Q0,  &
 -3425510460926691233.35152150484476Q0,  1183127973060179869.81759271465610Q0/ 
Data rcutoff/ 2.7500000000000000D+00/
Data fcutoff/ 2.1223260000000000D+03/
Data r0/ 2.5000000000000000D+00/
Data aa/ 2.1884620798613876D+04/
Data bb/-1.2336001555678876D+05/
Data cc/ 1.7538255947080589D+05/
Data k0/   5/
End
Double Precision Function V_Fausto_Cs_gs(x)
Parameter(Npg=14)
Implicit Real*8(A-H,O-Z)
Real*16 Pg,Sto,xq,Aux
Integer*4 k0
Common/Param_LJ_V_Fausto_Cs_gs/rcutoff,fcutoff,r0,aa,bb,cc,Pg(Npg),k0
If(x.le.rcutoff)Then
  V_Fausto_Cs_gs=fcutoff
  Return
Endif
If(x.le.r0)Then
  V_Fausto_Cs_gs=Abs(aa*x**2+bb*x+cc)
  Return
Endif
Sto=0.0q0
xq=1.0q0/x
Aux=xq**(k0+1)
Do k=1,Npg
  Sto=Sto+pg(k)*Aux
  Aux=Aux*xq
EndDo
V_Fausto_Cs_gs=Sto
Return
End
!
! Function for He-Cs(6p_Sigma), Fausto Cargnoni: more points (21-04-2015)
!
Block Data Inicio_LJ_V_Fausto_Cs_6p0
Implicit Real*8(A-H,O-Z)
Parameter(Npg=14)
Real*16 Pg
Integer*4 k0
Common/Param_LJ_V_Fausto_Cs_6p0/rcutoff,fcutoff,r0,aa,bb,cc,Pg(Npg),k0
Data Pg/                           &
  983012682.924529743263329958061894Q0, -84137127560.7181859398054466255159Q0,  3147363940329.22685865151335515968Q0,  &
 -68037224947910.0508055503708771544Q0,  948500179815737.587517335900477079Q0, -9034983578856463.08414720204540198Q0,  &
  60838778286068358.5669137120545770Q0, -294962476058276873.288394571504344Q0,  1035102455473919644.13005948263137Q0,  &
 -2608241792508316463.19876570809066Q0,  4603440277322316076.27999884747413Q0, -5404901585348550051.78209124683525Q0,  &
  3792737249528805624.44606122782020Q0, -1203728882114200188.77827500640065Q0/ 
Data rcutoff/ 2.6000000000000001D+00/
Data fcutoff/ 2.3170690000000000D+03/
Data r0/ 2.5000000000000000D+00/
Data aa/ 8.0635334129850098D+03/
Data bb/-4.6993574624164961D+04/
Data cc/ 7.0021714221092683D+04/
Data k0/   5/
End
Double Precision Function V_Fausto_Cs_6p0(x)
Parameter(Npg=14)
Implicit Real*8(A-H,O-Z)
Real*16 Pg,Sto,xq,Aux
Integer*4 k0
Common/Param_LJ_V_Fausto_Cs_6p0/rcutoff,fcutoff,r0,aa,bb,cc,Pg(Npg),k0
If(x.le.rcutoff)Then
  V_Fausto_Cs_6p0=fcutoff
  Return
Endif
If(x.le.r0)Then
  V_Fausto_Cs_6p0=Abs(aa*x**2+bb*x+cc)
  Return
Endif
Sto=0.0q0
xq=1.0q0/x
Aux=xq**(k0+1)
Do k=1,Npg
  Sto=Sto+pg(k)*Aux
  Aux=Aux*xq
EndDo
V_Fausto_Cs_6p0=Sto
Return
End
!
! Function for He-Cs(7s), Fausto Cargnoni: more points (21-04-2015)
!
Block Data Inicio_LJ_V_Fausto_Cs_7s
Implicit Real*8(A-H,O-Z)
Parameter(Npg=14)
Real*16 Pg
Integer*4 k0
Common/Param_LJ_V_Fausto_Cs_7s/rcutoff,fcutoff,r0,aa,bb,cc,Pg(Npg),k0
Data Pg/                           &
  3825265036.06471340558006213331397Q0, -106391699401.629756164352751452134Q0, -2415651814425.17303291100897069416Q0,  &
  160947890567832.110284939811381405Q0, -3494992921526184.29345177011789542Q0,  43853249263146779.6148549547055530Q0,  &
 -362634101182696589.346055675040595Q0,  2082630471929281529.33806179721184Q0, -8475927125141034661.39181490017403Q0,  &
  24430137804313288176.7891778391113Q0, -48844600104781726185.3454398788720Q0,  64493166634665917457.9924699139579Q0,  &
 -50603788390020951995.5383689991595Q0,  17874696930230973941.5742344069639Q0/ 
Data rcutoff/ 2.8999999999999999D+00/
Data fcutoff/ 1.3262239999999999D+03/
Data r0/ 2.7000000000000002D+00/
Data aa/ 2.7945811089763697D+03/
Data bb/-1.9344908344361000D+04/
Data cc/ 3.3924031337637469D+04/
Data k0/   5/
End
Double Precision Function V_Fausto_Cs_7s(x)
Parameter(Npg=14)
Implicit Real*8(A-H,O-Z)
Real*16 Pg,Sto,xq,Aux
Integer*4 k0
Common/Param_LJ_V_Fausto_Cs_7s/rcutoff,fcutoff,r0,aa,bb,cc,Pg(Npg),k0
If(x.le.rcutoff)Then
  V_Fausto_Cs_7s=fcutoff
  Return
Endif
If(x.le.r0)Then
  V_Fausto_Cs_7s=Abs(aa*x**2+bb*x+cc)
  Return
Endif
Sto=0.0q0
xq=1.0q0/x
Aux=xq**(k0+1)
Do k=1,Npg
  Sto=Sto+pg(k)*Aux
  Aux=Aux*xq
EndDo
V_Fausto_Cs_7s=Sto
Return
End
!
! Function for He-Cs(6p_Pi), Fausto Cargnoni: more points (21-04-2015)
!
Block Data Inicio_LJ_V_Fausto_Cs_6p1
Implicit Real*8(A-H,O-Z)
Parameter(Npg=14)
Real*16 Pg
Integer*4 k0
Common/Param_LJ_V_Fausto_Cs_6p1/rcutoff,fcutoff,r0,aa,bb,cc,Pg(Npg),k0
Data Pg/                           &
 -28916389.8646917185914187350216402Q0,  1698212759.47458891331191572250302Q0, -46334669587.0885737998848348571058Q0,  &
  757625158504.461142433946867564510Q0, -8235739328760.26112437184337510191Q0,  62619687246020.4804186731767136627Q0,  &
 -342217765797573.666436871785982512Q0,  1361896256113801.28031727108845724Q0, -3953071532136512.29016059435905051Q0,  &
  8282683882411200.31626276523433664Q0, -12202657183519187.5983436477674604Q0,  11995345991890722.2461560807048455Q0,  &
 -7064895065092973.46664535404252966Q0,  1885957227405596.86900447464232992Q0/ 
Data rcutoff/ 2.3999999999999999D+00/
Data fcutoff/ 2.6430850000000000D+03/
Data r0/ 2.1000000000000001D+00/
Data aa/ 3.6295506204203572D+04/
Data bb/-1.8721605892180325D+05/
Data cc/ 2.4278672572368380D+05/
Data k0/   5/
End
Double Precision Function V_Fausto_Cs_6p1(x)
Parameter(Npg=14)
Implicit Real*8(A-H,O-Z)
Real*16 Pg,Sto,xq,Aux
Integer*4 k0
Common/Param_LJ_V_Fausto_Cs_6p1/rcutoff,fcutoff,r0,aa,bb,cc,Pg(Npg),k0
If(x.le.rcutoff)Then
  V_Fausto_Cs_6p1=fcutoff
  Return
Endif
If(x.le.r0)Then
  V_Fausto_Cs_6p1=Abs(aa*x**2+bb*x+cc)
  Return
Endif
Sto=0.0q0
xq=1.0q0/x
Aux=xq**(k0+1)
Do k=1,Npg
  Sto=Sto+pg(k)*Aux
  Aux=Aux*xq
EndDo
V_Fausto_Cs_6p1=Sto
Return
End
!
! Function for He-Rb gs(5s), from Fausto Cargnoni
!
Block Data Inicio_LJ_V_Fausto_Rb_gs
Implicit Real*8(A-H,O-Z)
Parameter(Npg=14)
Real*16 Pg
Integer*4 k0
Common/Param_LJ_V_Fausto_Rb_gs/rcutoff,fcutoff,r0,aa,bb,cc,Pg(Npg),k0
Data Pg/                           &
  19100771.2956833483968680974522071Q0, -1172294215.25983631779261524084372Q0,  33435877600.2318788619551266974640Q0,  &
 -626736755694.868731078365328255351Q0,  8889903709250.47980290206377580454Q0, -98809797408403.8989098633620172560Q0,  &
  838535230285386.179156473882388379Q0, -5250871933035878.90801737880175866Q0,  23702922123542559.7263712283775712Q0,  &
 -75653235518269118.4852893687161765Q0,  166080556623113322.510632192031701Q0, -238324245940614178.157128370536680Q0,  &
  201208460378793600.224928931873042Q0, -75780945879534690.8963660066902431Q0/ 
Data rcutoff/ 2.7500000000000000D+00/
Data fcutoff/ 1.7084310000000000D+03/
Data r0/ 2.5000000000000000D+00/
Data aa/ 2.5322093281985394D+03/
Data bb/-1.6991936870778129D+04/
Data cc/ 2.9338073941076589D+04/
Data k0/   5/
End
Double Precision Function V_Fausto_Rb_gs(x)
Parameter(Npg=14)
Implicit Real*8(A-H,O-Z)
Real*16 Pg,Sto,xq,Aux
Integer*4 k0
Common/Param_LJ_V_Fausto_Rb_gs/rcutoff,fcutoff,r0,aa,bb,cc,Pg(Npg),k0
If(x.le.rcutoff)Then
  V_Fausto_Rb_gs=fcutoff
  Return
Endif
If(x.le.r0)Then
  V_Fausto_Rb_gs=Abs(aa*x**2+bb*x+cc)
  Return
Endif
Sto=0.0q0
xq=1.0q0/x
Aux=xq**(k0+1)
Do k=1,Npg
  Sto=Sto+pg(k)*Aux
  Aux=Aux*xq
EndDo
V_Fausto_Rb_gs=Sto
Return
End
!
! Function for He-Rb(5p_Sigma), Fausto Cargnoni: more points (21-04-2015)
!
Block Data Inicio_LJ_V_Fausto_Rb_5p0
Implicit Real*8(A-H,O-Z)
Parameter(Npg=14)
Real*16 Pg
Integer*4 k0
Common/Param_LJ_V_Fausto_Rb_5p0/rcutoff,fcutoff,r0,aa,bb,cc,Pg(Npg),k0
Data Pg/                           &
  755939876.954161676057804724548838Q0, -61892314758.9221075767972745682921Q0,  2217888483872.67461384962575591948Q0,  &
 -45968257647250.0207867967385818922Q0,  614701763499008.868032004191219577Q0, -5618036039629895.48401381180406123Q0,  &
  36301923962954331.4203297556619919Q0, -168897353040772297.040299997362214Q0,  568755905574416789.719347690523984Q0,  &
 -1375110600761493082.00782733743110Q0,  2328499818364179948.01741052654689Q0, -2622699174021865946.26117767851036Q0,  &
  1765450061925162442.37646295715979Q0, -537479750788251198.126286046944468Q0/ 
Data rcutoff/ 2.3999999999999999D+00/
Data fcutoff/ 3.3863560000000002D+03/
Data r0/ 2.3999999999999999D+00/
Data aa/ 5.1234896077967460D+03/
Data bb/-2.8817394513216812D+04/
Data cc/ 4.3036802950585261D+04/
Data k0/   5/
End
Double Precision Function V_Fausto_Rb_5p0(x)
Parameter(Npg=14)
Implicit Real*8(A-H,O-Z)
Real*16 Pg,Sto,xq,Aux
Integer*4 k0
Common/Param_LJ_V_Fausto_Rb_5p0/rcutoff,fcutoff,r0,aa,bb,cc,Pg(Npg),k0
If(x.le.rcutoff)Then
  V_Fausto_Rb_5p0=fcutoff
  Return
Endif
If(x.le.r0)Then
  V_Fausto_Rb_5p0=Abs(aa*x**2+bb*x+cc)
  Return
Endif
Sto=0.0q0
xq=1.0q0/x
Aux=xq**(k0+1)
Do k=1,Npg
  Sto=Sto+pg(k)*Aux
  Aux=Aux*xq
EndDo
V_Fausto_Rb_5p0=Sto
Return
End
!
! Function for He-Rb(6p_Pi), Fausto Cargnoni: more points (21-04-2015)
!
Block Data Inicio_LJ_V_Fausto_Rb_5p1
Implicit Real*8(A-H,O-Z)
Parameter(Npg=14)
Real*16 Pg
Integer*4 k0
Common/Param_LJ_V_Fausto_Rb_5p1/rcutoff,fcutoff,r0,aa,bb,cc,Pg(Npg),k0
Data Pg/                           &
  21065388.4588612130967750242379954Q0, -1356769725.00046943557636990593514Q0,  38031190615.3371756235401161001838Q0,  &
 -628459250285.123546375190971037063Q0,  6849776778387.97809323521638065669Q0, -52146078017554.2500446984422353406Q0,  &
  285792287469017.676532214877262703Q0, -1143063591357831.95828165866077991Q0,  3339589997615048.00774407054849081Q0,  &
 -7045930377343015.98657569250566483Q0,  10445469637201842.5898884645426793Q0, -10315641071334019.0899871943233809Q0,  &
  6090035673546989.48603755367194306Q0, -1625247766723261.80485107267140556Q0/ 
Data rcutoff/ 2.2000000000000002D+00/
Data fcutoff/ 2.7312600000000002D+03/
Data r0/ 1.0000000000000000D+00/
Data aa/ 3.2376196740691164D+04/
Data bb/-1.5866800201408935D+05/
Data cc/ 1.9510007184216869D+05/
Data k0/   5/
End
Double Precision Function V_Fausto_Rb_5p1(x)
Parameter(Npg=14)
Implicit Real*8(A-H,O-Z)
Real*16 Pg,Sto,xq,Aux
Integer*4 k0
Common/Param_LJ_V_Fausto_Rb_5p1/rcutoff,fcutoff,r0,aa,bb,cc,Pg(Npg),k0
If(x.le.rcutoff)Then
  V_Fausto_Rb_5p1=fcutoff
  Return
Endif
If(x.le.r0)Then
  V_Fausto_Rb_5p1=Abs(aa*x**2+bb*x+cc)
  Return
Endif
Sto=0.0q0
xq=1.0q0/x
Aux=xq**(k0+1)
Do k=1,Npg
  Sto=Sto+pg(k)*Aux
  Aux=Aux*xq
EndDo
V_Fausto_Rb_5p1=Sto
Return
End
!
! Potencial del Rb+ calculat per en Fasuto et al. amb long range del tipus 1/r**4
!
Block Data Inicio_LJ_V_Fausto_Rb_plus_gs
Implicit Real*8(A-H,O-Z)
Parameter(Npg=13)
Real*16 Pg
Integer*4 k0
Common/Param_LJ_V_Fausto_Rb_plus_gs/rcutoff,fcutoff,r0,aa,bb,cc,Pg(Npg),k0
Data Pg/                           &
 -406622.375170394176501817185213738Q0,  21332931.6604765915055670610776133Q0, -515840595.362288581156823966609460Q0,  &
  7286963557.49977496017331947024218Q0, -67050629086.2414609732865947739718Q0,  423545364741.341210425847019599572Q0,  &
 -1884804263958.89912314686654731090Q0,  5961212677076.96608868548002263051Q0, -13327238149555.0874553559671429585Q0,  &
  20600169905046.3280051663330937625Q0, -20966862545953.0513329124868195299Q0,  12657369730747.1509804505111593368Q0,  &
 -3437256517661.33873687152088791033Q0/ 
Data rcutoff/ 2.0000000000000000D+00/
Data fcutoff/ 6.9130180000000000D+03/
Data r0/ 2.1000000000000001D+00/
Data aa/ 2.2004608525779306D+04/
Data bb/-1.1197802949237564D+05/
Data cc/ 1.4285064271139764D+05/
Data k0/   3/
End
Double Precision Function V_Fausto_Rb_plus_gs(x)
Parameter(Npg=13)
Implicit Real*8(A-H,O-Z)
Real*16 Pg,Sto,xq
Integer*4 k0
Common/Param_LJ_V_Fausto_Rb_plus_gs/rcutoff,fcutoff,r0,aa,bb,cc,Pg(Npg),k0
If(x.le.rcutoff)Then
  V_Fausto_Rb_plus_gs=fcutoff
  Return
Endif
If(x.le.r0)Then
  V_Fausto_Rb_plus_gs=Abs(aa*x**2+bb*x+cc)
  Return
Endif
Sto=0.0q0
xq=1.0q0/x
Do k=1,Npg
  Sto=Sto+pg(k)*xq**(k+k0)
EndDo
V_Fausto_Rb_plus_gs=Sto
Return
End
!
! Function for He-Rb(6s), Fausto Cargnoni: more points (21-04-2015)
!
Block Data Inicio_LJ_V_Fausto_Rb_6s
Implicit Real*8(A-H,O-Z)
Parameter(Npg=14)
Real*16 Pg
Integer*4 k0
Common/Param_LJ_V_Fausto_Rb_6s/rcutoff,fcutoff,r0,aa,bb,cc,Pg(Npg),k0
Data Pg/                           &
 -4780502820.90050189325416969698484Q0,  401314863378.456544916850650103099Q0, -14753972047325.5303952564485971199Q0,  &
  314957742848247.648994876344500017Q0, -4375394053042189.01319139625006427Q0,  42025746888366128.5201928564392822Q0,  &
 -288764135760081719.682880070083058Q0,  1443952790541320699.20726573664655Q0, -5274837452555578768.65696176481809Q0,  &
  13946175032217141852.1624775870371Q0, -26003985198509507719.3486720078825Q0,  32447777154760679012.5074313031215Q0,  &
 -24326312738770487365.6660811096881Q0,  8287435576981232747.15789107114440Q0/ 
Data rcutoff/ 3.0000000000000000D+00/
Data fcutoff/ 2.5755680000000002D+03/
Data r0/ 2.8999999999999999D+00/
Data aa/-8.1144161152598926D+01/
Data bb/ 2.1803190565823604D+02/
Data cc/ 2.6582835372153691D+03/
Data k0/   5/
End
Double Precision Function V_Fausto_Rb_6s(x)
Parameter(Npg=14)
Implicit Real*8(A-H,O-Z)
Real*16 Pg,Sto,xq,Aux
Integer*4 k0
Common/Param_LJ_V_Fausto_Rb_6s/rcutoff,fcutoff,r0,aa,bb,cc,Pg(Npg),k0
If(x.le.rcutoff)Then
  V_Fausto_Rb_6s=fcutoff
  Return
Endif
If(x.le.r0)Then
  V_Fausto_Rb_6s=Abs(aa*x**2+bb*x+cc)
  Return
Endif
Sto=0.0q0
xq=1.0q0/x
Aux=xq**(k0+1)
Do k=1,Npg
  Sto=Sto+pg(k)*Aux
  Aux=Aux*xq
EndDo
V_Fausto_Rb_6s=Sto
Return
End

Double Precision Function V_grafeno(z)
!
! Potencial He - Superficie de Grafeno, de Pilar de Lara
!
Implicit Real*8(A-H,O-Z)
!
! Parametres del potencial
!
!  z_cutoff=1.4, umax=3.282645d4
!
!
! A     = 554811 meV
! alpha = 3.21193 Ang.-1
! C4    = 1531.83 meV. Angs4
! C6    = 22259 meV Angs6
!
Data C_meV_to_K/11.604448d0/
Data A/554811.d0/, alpha/3.21193d0/, C4/1531.83d0/, C6/22259.d0/

V_average=A*exp(-alpha*z)-C4/(z**4)-C6/(z**6)

V_grafeno=V_average*C_meV_to_K
Return
End
Double precision Function V_grafeno_sin_dispersion(z)
!
! Potencial He - Superficie de Grafeno, de Pilar de Lara
!
Implicit Real*8(A-H,O-Z)
!
! Parametres del potencial
!
!  z_cutoff=1.4, umax=5.107576d4
!
! alpha = 1.54605 Angs.-1
! D     = 0.493191 meV
! Ze    = 4.34871 Angs.
!
Data C_meV_to_K/11.604448d0/
Data D/0.493191d0/, alpha/1.54605d0/, Ze/4.34871d0/

V_Morse=D*(1.0d0-dexp(-alpha*(z-Ze)))**2-D

V_grafeno_sin_dispersion=V_Morse*C_meV_to_K
Return
End
Double Precision Function V_TiO2(Z)
!
!  Ultim potencial He-TiO2, calculat per la Pilar
!
!     Z es la distancia
!     al Ti(5f) y está en AA, VZ en cm-1

IMPLICIT REAL*8(A-H,O-Z)
Data cmm1toK/1.4387770d0/
Data de/95.9583d0/,alpha/1.47089d0/,ze/4.22534d0/
Data z0/4.77765d0/,c3/8.85163d0/, aa/1.06161d0/
vz= de*(1-dexp(-alpha*(z-ze)))**2-de 
vzd= - 0.5d0 * (1.d0 + dtanh ((z-z0)/aa)*c3*(c3/z)**3)
V_TiO2=(vz+vzd)*cmm1toK
return
end
Double Precision Function V_Au(rr)
Implicit None
!Implicit Real*8(A-H,R-Z)
!
! Interaccion Oro-He: entrada en Ang, salida en K
!

Real  (Kind=8)  :: rr,a1,a2,a3,a4,a5,a6,a7,a8,a9,r,Re,De


Data a1/1.6417d0/,a2/-0.6164d0/,a3/0.3770d0/,a4/-0.0549d0/
Data a5/0.0018d0/,a6/0.0046d0/,a7/-0.0010d0/,a8/0.0001d0/
Data a9/2.7d-6/,Re/4.124d0/,De/20.301143d0/
! rms = 0.0034
! De = 14.110 ! in cm^-1

r = rr - Re
!V_Au = -De*(1.d0 + a1*r + a2*r**2 + a3*r**3 + a4*r**4 + a5*r**5 + a6*r**6 + a7*r**7 + a8*r**8 + a9*r**9)*dexp(-a1*r)

V_Au = -De*(1.d0 + r*(a1 + r*(a2 + r*(a3 + r*(a4 + r*(a5 + r*(a6 + r*(a7 + r*(a8 + r*a9)))))))))*dexp(-a1*r)

End Function V_Au
      Double Precision Function V_Au_TiO2(z)
      Implicit None
      Real (Kind=8)  :: C6=-4.18518d6, C8=6.54459D7, C10=-3.73112D8, C12=9.44942D8, C14=-8.97265D8
      Real (kind=8)  :: z, Aux, Aux2, Aux6, CmeV_to_K=11.604448d0
      If(z.Gt.2.)Then
        Aux=1.d0/z
        Aux2=Aux**2
        Aux6=Aux2**3
        V_Au_TiO2 = (((((C14*Aux2+C12)*Aux2+C10)*Aux2+C8)*Aux2+C6)*Aux6)*CmeV_to_K
      Else
        V_Au_TiO2 = 1821.44287109375d0*CmeV_to_K
      Endif
      Return
      End

      Double Precision Function V_dAu_TiO2(z)
      Implicit None
      Real (Kind=8)  :: C6=-4.18518d6, C8=6.54459D7, C10=-3.73112D8, C12=9.44942D8, C14=-8.97265D8
      Real (kind=8)  :: z, Aux, Aux2, Aux7, CmeV_to_K=11.604448d0
      If(z.Gt.2.)Then
        Aux=1.d0/z
        Aux2=Aux**2
        Aux7=Aux2**3*Aux
        V_dAu_TiO2 = (((((-14.0d0*C14*Aux2-12.0d0*C12)*Aux2-10.0d0*C10)*Aux2-8.0d0*C8)*Aux2-6.0d0*C6)*Aux7)*CmeV_to_K
      Else
        V_dAu_TiO2 = -5415.35400390625d0*CmeV_to_K
      Endif
      Return
      End
!
! Potencial pH2-He de Gianturco ajustado (nouevo promedio 29/5/2015)
!
Block Data Inicio_LJ_V_ph2_He
Implicit Real*8(A-H,O-Z)
Parameter(Npg=20)
Real*16 Pg
Integer*4 k0
Common/Param_LJ_V_ph2_He/rcutoff,fcutoff,r0,aa,bb,cc,Pg(Npg),k0
Data Pg/                           &
 -1163439.97300733181599835248769270Q0,  47292172.4090617880619017432408668Q0, -841086014.436697107304093367838305Q0,  &
  8465900672.40128294507703823716660Q0, -53828298395.5314532464705938241342Q0,  227557134285.265792622701141972689Q0,  &
 -650430073601.165531760983033096775Q0,  1218814073699.63463676672914897444Q0, -1216214794159.68959709558137823136Q0,  &
 -567844642080.136519610799404824244Q0,  4746308391191.44071750176385675985Q0, -10143616872740.8356624374938660844Q0,  &
  14089430994274.8504311602782253396Q0, -14374323621667.8460913194491051522Q0,  11078725941332.5159330085651217405Q0,  &
 -6402732017618.46278859102735729547Q0,  2688141511623.61334745699379906995Q0, -771949591132.651962769099634129229Q0,  &
  135129364439.500429837521587722877Q0, -10838377535.0308515558156365736038Q0/ 
Data rcutoff/ 1.0000000000000000D+00/
Data fcutoff/ 6.5305870000000003D+04/
Data r0/ 5.0000000000000000D-01/
Data aa/-3.6541132591800960D+06/
Data bb/ 0.0000000000000000D+00/
Data cc/ 8.5553041453631725D+06/
Data k0/   5/
End
Double Precision Function V_ph2_He(x)
Parameter(Npg=20)
Implicit Real*8(A-H,O-Z)
Real*16 Pg,Sto,xq,Aux
Integer*4 k0
Common/Param_LJ_V_ph2_He/rcutoff,fcutoff,r0,aa,bb,cc,Pg(Npg),k0
If(x.le.rcutoff)Then
  V_ph2_He=fcutoff
  Return
Endif
If(x.le.r0)Then
  V_ph2_He=Abs(aa*x**2+bb*x+cc)
  Return
Endif
Sto=0.0q0
xq=1.0q0/x
Aux=xq**(k0+1)
Do k=1,Npg
  Sto=Sto+pg(k)*Aux
  Aux=Aux*xq
EndDo
V_ph2_He=Sto
Return
End
!
! Function for He-Rb(6p_Sigma), From Pascale calculation (06-05-2015)
!
Block Data Inicio_LJ_V_Pascale_Rb_6p0
Implicit Real*8(A-H,O-Z)
Parameter(Npg=18)
Real*16 Pg
Integer*4 k0
Common/Param_LJ_V_Pascale_Rb_6p0/rcutoff,fcutoff,r0,aa,bb,cc,Pg(Npg),k0
Data Pg/                           &
 -27238209654.0970551997286039358492Q0,  3155212616112.62367579232721486859Q0, -159215352636228.495876363247390109Q0,  &
  4657247444544819.83119004845414435Q0, -89015024933389427.4025630902664066Q0,  1190266214999619563.79458870093428Q0,  &
 -11626907744717667292.0430483506713Q0,  85342874325414564884.2727262312123Q0, -479206558620624734149.448129236565Q0,  &
  2079235750526732316320.70782350523Q0, -6995096564286351837337.16158701126Q0,  18191735476673693147864.7261298775Q0,  &
 -36198334606432096939869.4997435395Q0,  54039097780171191560513.9501513440Q0, -58516086287882353120407.4077308311Q0,  &
  43341928857121805849873.8181966873Q0, -19613221630184108600706.5386003175Q0,  4082648274688574570356.65959774196Q0/ 
Data rcutoff/ 2.3500000000000001D+00/
Data fcutoff/ 2.2097469999999998D+03/
Data r0/-3.2000000000000002D+00/
Data aa/ 1.2723431741929555D+03/
Data bb/ 0.0000000000000000D+00/
Data cc/-6.9728651248261986D+03/
Data k0/   5/
End
Double Precision Function V_Pascale_Rb_6p0(x)
Parameter(Npg=18)
Implicit Real*8(A-H,O-Z)
Real*16 Pg,Sto,xq,Aux
Integer*4 k0
Common/Param_LJ_V_Pascale_Rb_6p0/rcutoff,fcutoff,r0,aa,bb,cc,Pg(Npg),k0
If(x.le.rcutoff)Then
  V_Pascale_Rb_6p0=fcutoff
  Return
Endif
If(x.le.r0)Then
  V_Pascale_Rb_6p0=Abs(aa*x**2+bb*x+cc)
  Return
Endif
Sto=0.0q0
xq=1.0q0/x
Aux=xq**(k0+1)
Do k=1,Npg
  Sto=Sto+pg(k)*Aux
  Aux=Aux*xq
EndDo
V_Pascale_Rb_6p0=Sto
Return
End
!
! Function for He-Rb(6p_Pi), From Pascale calculation (06-05-2015)
!
Block Data Inicio_LJ_V_Pascale_Rb_6p1
Implicit Real*8(A-H,O-Z)
Parameter(Npg=18)
Real*16 Pg
Integer*4 k0
Common/Param_LJ_V_Pascale_Rb_6p1/rcutoff,fcutoff,r0,aa,bb,cc,Pg(Npg),k0
Data Pg/                           &
  473687428.390644638435766046481170Q0, -64347148419.7845939200849296048072Q0,  3637808075975.30386305837579355438Q0,  &
 -117799828212905.378494464152006355Q0,  2485273579964282.73699052550747524Q0, -36637637078170660.1214956082180008Q0,  &
  393964578545592969.739990894637550Q0, -3176400344113005754.50295649136964Q0,  19540934127522939346.8829848452431Q0,  &
 -92640862627641172344.0073533069276Q0,  339664614593611552754.420694144245Q0, -960564218001611305059.888936433560Q0,  &
  2074959812887202322685.47117302676Q0, -3359492547653334985410.92994317689Q0,  3944818157081899283205.45238007755Q0,  &
 -3171181660093470107341.15793633147Q0,  1560672268303756719160.88518302906Q0, -354536702973118696499.223764295971Q0/ 
Data rcutoff/ 2.1499999999999999D+00/
Data fcutoff/ 3.0895619999999999D+03/
Data r0/-3.2000000000000002D+00/
Data aa/ 1.5092579075989399D+03/
Data bb/ 0.0000000000000000D+00/
Data cc/-6.9138253346225456D+03/
Data k0/   5/
End
Double Precision Function V_Pascale_Rb_6p1(x)
Parameter(Npg=18)
Implicit Real*8(A-H,O-Z)
Real*16 Pg,Sto,xq,Aux
Integer*4 k0
Common/Param_LJ_V_Pascale_Rb_6p1/rcutoff,fcutoff,r0,aa,bb,cc,Pg(Npg),k0
If(x.le.rcutoff)Then
  V_Pascale_Rb_6p1=fcutoff
  Return
Endif
If(x.le.r0)Then
  V_Pascale_Rb_6p1=Abs(aa*x**2+bb*x+cc)
  Return
Endif
Sto=0.0q0
xq=1.0q0/x
Aux=xq**(k0+1)
Do k=1,Npg
  Sto=Sto+pg(k)*Aux
  Aux=Aux*xq
EndDo
V_Pascale_Rb_6p1=Sto
Return
End
!  Funcio que calcula el potencial Ar-He, entrada \AA, sortida K (Tang & Toennies JCP118 4976 (3003))
!
!  r_cutoff = 1.9d0; umax = 1.5101972249d4
!
!
Double precision function V_Ar_He(r)
Implicit none
real      (kind=8) :: AA     = 124.3d0
real      (kind=8) :: BB     = 2.153d0
real      (kind=8) :: C6     = 9.538d0
real      (kind=8) :: C8     = 167.5d0
real      (kind=8) :: C10    = 3701.d0 
real      (kind=8) :: rbohr  = 0.5292d0
real      (kind=8) :: Hartree_to_K  = 27.21d0*11604.d0
real (kind=8) :: r, rr, Stoe, Stor, F6, F8, F10, Sto, Aux
real (kind=8) :: c1or2, c1or6
real (kind=8) :: r_cutoff=1.9d0, f_cutoff=1.5101972249d4
real (kind=8) :: Fact(0:10)
Logical  :: Lfirst = .true.
Integer (Kind=4) :: i
Save     :: Lfirst, Fact
if(r.Le.r_cutoff)Then
  V_Ar_He = f_cutoff
  Return
Endif        
If(Lfirst)Then
!  Write(6,'("#Control per veure si nomes hi pasem una vegada")')
  Fact(0) = 1.d0; Fact(1) = 1.d0
  Do i=2, 10
    Fact(i) = Fact(i-1)*i
  EndDo
  Lfirst = .false.
Endif
rr = r/rbohr
c1or2 = 1.d0/(rr*rr)
c1or6 = c1or2*c1or2*c1or2
Stor = bb*rr
Stoe = Exp(-Stor)
Sto  = 1.d0
Aux = Stor 
Do i=1, 6
  Sto = Sto + Aux/Fact(i)
  Aux = Aux*Stor
Enddo
F6 = 1.d0 - Stoe*Sto
Do i=7, 8
  Sto = Sto + Aux/Fact(i)
  Aux = Aux*Stor
Enddo
F8 = 1.d0 - Stoe*Sto
Do i=9, 10
  Sto = Sto + Aux/Fact(i)
  Aux = Aux*Stor
Enddo
F10 = 1.d0 - Stoe*Sto
Aux = AA*Stoe - (C6*F6 + (C8*F8 + C10*F10*c1or2)*c1or2)*c1or6
V_Ar_He = Aux*Hartree_to_K
Return
End
!
!  Funcio que calcula el potencial Xe-He, entrada \AA, sortida K (Tang & Toennies J. Chem. Phys. 118 (2003) 4976
!
!
Double precision function V_Xe_He(r)
Implicit none
real      (kind=8) :: AA     = 95.9d0
real      (kind=8) :: BB     = 1.853d0
real      (kind=8) :: C6     = 19.54d0
real      (kind=8) :: C8     = 525.d0
real      (kind=8) :: C10    = 16670.d0 
real      (kind=8) :: rbohr  = 0.5292d0
real      (kind=8) :: Hartree_to_K  = 27.21d0*11604.d0
real (kind=8) :: r, rr, Stoe, Stor, F6, F8, F10, Sto, Aux
real (kind=8) :: c1or2, c1or6
real (kind=8) :: r_cutoff=2.2d0, f_cutoff=1.1678125552d4
real (kind=8) :: Fact(0:10)
Logical  :: Lfirst = .true.
Integer (Kind=4) :: i
Save     :: Lfirst, Fact
If(r.Le.r_cutoff)Then
  V_Xe_He = f_cutoff
  Return
EndIf  
If(Lfirst)Then
!  Write(6,'("#Control per veure si nomes hi pasem una vegada")')
  Fact(0) = 1.d0; Fact(1) = 1.d0
  Do i=2, 10
    Fact(i) = Fact(i-1)*i
  EndDo
  Lfirst = .false.
Endif
rr = r/rbohr
c1or2 = 1.d0/(rr*rr)
c1or6 = c1or2*c1or2*c1or2
Stor = bb*rr
Stoe = Exp(-Stor)
Sto  = 1.d0
Aux = Stor 
Do i=1, 6
  Sto = Sto + Aux/Fact(i)
  Aux = Aux*Stor
Enddo
F6 = 1.d0 - Stoe*Sto
Do i=7, 8
  Sto = Sto + Aux/Fact(i)
  Aux = Aux*Stor
Enddo
F8 = 1.d0 - Stoe*Sto
Do i=9, 10
  Sto = Sto + Aux/Fact(i)
  Aux = Aux*Stor
Enddo
F10 = 1.d0 - Stoe*Sto
Aux = AA*Stoe - (C6*F6 + (C8*F8 + C10*F10*c1or2)*c1or2)*c1or6
V_Xe_He = Aux*Hartree_to_K
Return
End
!
! Ajust potencial 5p_Sigma de Pascale (29-10-2015) (unitats correctes)!!
!
Block Data Inicio_LJ_V_Rb_5p_sigma
Implicit Real*8(A-H,O-Z)
Parameter(Npg=11)
Real*16 Pg
Integer*4 k0
Common/Param_LJ_V_Rb_5p_sigma/rcutoff,fcutoff,r0,aa,bb,cc,Pg(Npg),k0
Data Pg/                           &
 -16369756.1679657643374850977670248Q0,  358314586.951384914502183247424336Q0,  25272036141.5453163057909074234291Q0,  &
 -1211017031542.79911725659637667590Q0,  21309419454172.0588581060167965039Q0, -197508142095030.863285799109515175Q0,  &
  1082206330864640.68613317062793854Q0, -3643715772086083.49000422254676492Q0,  7432352868423233.80679404088747201Q0,  &
 -8450047598780025.77211999338779863Q0,  4117216576621054.92026930080355380Q0/ 
Data rcutoff/ 3.1750600000000002D+00/
Data fcutoff/ 2.7347030000000000D+03/
Data r0/ 2.5000000000000000D+00/
Data aa/-1.4702350246680987D+02/
Data bb/ 0.0000000000000000D+00/
Data cc/ 4.3142207107069107D+03/
Data k0/   5/
End
Double Precision Function V_Rb_5p_sigma(x)
Parameter(Npg=11)
Implicit Real*8(A-H,O-Z)
Real*16 Pg,Sto,xq,Aux
Integer*4 k0
Common/Param_LJ_V_Rb_5p_sigma/rcutoff,fcutoff,r0,aa,bb,cc,Pg(Npg),k0
If(x.le.rcutoff)Then
  V_Rb_5p_sigma=fcutoff
  Return
Endif
If(x.le.r0)Then
  V_Rb_5p_sigma=Abs(aa*x**2+bb*x+cc)
  Return
Endif
Sto=0.0q0
xq=1.0q0/x
Aux=xq**(k0+1)
Do k=1,Npg
  Sto=Sto+pg(k)*Aux
  Aux=Aux*xq
EndDo
V_Rb_5p_sigma=Sto
Return
End
!
! Ajust potencial 5_Pi de Pascale (29-10-2015) (unitats correctes)!!
!
Block Data Inicio_LJ_V_Rb_5p_pi
Implicit Real*8(A-H,O-Z)
Parameter(Npg=11)
Real*16 Pg
Integer*4 k0
Common/Param_LJ_V_Rb_5p_pi/rcutoff,fcutoff,r0,aa,bb,cc,Pg(Npg),k0
Data Pg/                           &
 -8741753.92488395686001702665966285Q0,  486176414.868784278332186303239805Q0, -11591463068.8021123878989034988742Q0,  &
  148126647810.351946447172854810821Q0, -1164177194524.65854784396663864098Q0,  5992071282192.61641098759405074505Q0,  &
 -20640521486966.5879282525969265615Q0,  47219174940119.0277918081596800239Q0, -68878071746805.1901021018758721574Q0,  &
  57983399307116.7535300990213004201Q0, -21431560362273.3380432160483734331Q0/ 
Data rcutoff/ 2.1167099999999999D+00/
Data fcutoff/ 4.8593370000000004D+03/
Data r0/ 2.0000000000000000D+00/
Data aa/-3.0976037301688293D+03/
Data bb/ 0.0000000000000000D+00/
Data cc/ 1.8738030191588245D+04/
Data k0/   5/
End
Double Precision Function V_Rb_5p_pi(x)
Parameter(Npg=11)
Implicit Real*8(A-H,O-Z)
Real*16 Pg,Sto,xq,Aux
Integer*4 k0
Common/Param_LJ_V_Rb_5p_pi/rcutoff,fcutoff,r0,aa,bb,cc,Pg(Npg),k0
If(x.le.rcutoff)Then
  V_Rb_5p_pi=fcutoff
  Return
Endif
If(x.le.r0)Then
  V_Rb_5p_pi=Abs(aa*x**2+bb*x+cc)
  Return
Endif
Sto=0.0q0
xq=1.0q0/x
Aux=xq**(k0+1)
Do k=1,Npg
  Sto=Sto+pg(k)*Aux
  Aux=Aux*xq
EndDo
V_Rb_5p_pi=Sto
Return
End
!
! Function for He-Rb(6p_Sigma), From Pascale calculation (30-09-2015)unitats corre
!
Block Data Inicio_LJ_V_Rb_6p_Sigma
Implicit Real*8(A-H,O-Z)
Parameter(Npg=18)
Real*16 Pg
Integer*4 k0
Common/Param_LJ_V_Rb_6p_Sigma/rcutoff,fcutoff,r0,aa,bb,cc,Pg(Npg),k0
Data Pg/                           &
 -39855258349.8562779111031444973651Q0,  4616917836286.01917946198638628870Q0, -233048345693666.247944477233277838Q0,  &
  6821255093289797.55038593246614875Q0, -130498399702449880.755284684685938Q0,  1747068821575881055.36281494557793Q0,  &
 -17090582433391375124.6575563934023Q0,  125655099724694531869.495505442861Q0, -706879748387175611849.223558587160Q0,  &
  3073458276584230945887.87133405050Q0, -10363697041788686021395.9401520884Q0,  27020880156532784314594.6844200979Q0,  &
 -53918859334910322896215.0051685128Q0,  80747814635804792925159.0352864619Q0, -87748853076398042247783.7300727085Q0,  &
  65257436010976109146648.5575394867Q0, -29668229502584876388311.8179117254Q0,  6209323876594826013383.91461241976Q0/ 
Data rcutoff/ 2.3500000000000001D+00/
Data fcutoff/ 3.0013080000000000D+03/
Data r0/-3.2000000000000002D+00/
Data aa/ 1.8414125504060087D+03/
Data bb/ 0.0000000000000000D+00/
Data cc/-1.0102382215546504D+04/
Data k0/   5/
End
Double Precision Function V_Rb_6p_Sigma(x)
Parameter(Npg=18)
Implicit Real*8(A-H,O-Z)
Real*16 Pg,Sto,xq,Aux
Integer*4 k0
Common/Param_LJ_V_Rb_6p_Sigma/rcutoff,fcutoff,r0,aa,bb,cc,Pg(Npg),k0
If(x.le.rcutoff)Then
  V_Rb_6p_Sigma=fcutoff
  Return
Endif
If(x.le.r0)Then
  V_Rb_6p_Sigma=Abs(aa*x**2+bb*x+cc)
  Return
Endif
Sto=0.0q0
xq=1.0q0/x
Aux=xq**(k0+1)
Do k=1,Npg
  Sto=Sto+pg(k)*Aux
  Aux=Aux*xq
EndDo
V_Rb_6p_Sigma=Sto
Return
End
!
! Function for He-Rb(6p_Pi), From Pascale calculation (30-09-2015) unitast correct
!
Block Data Inicio_LJ_V_Rb_6p_Pi
Implicit Real*8(A-H,O-Z)
Parameter(Npg=18)
Real*16 Pg
Integer*4 k0
Common/Param_LJ_V_Rb_6p_Pi/rcutoff,fcutoff,r0,aa,bb,cc,Pg(Npg),k0
Data Pg/                           &
  840075983.737711643987850807092476Q0, -109565252185.238302267920541824536Q0,  6040041697493.86281523568261964119Q0,  &
 -192096523030023.441530962612113581Q0,  3997548804666069.64682575835922540Q0, -58293811547018697.3031323722230219Q0,  &
  621277997758149591.165075250840335Q0, -4971929235445194262.19280905259923Q0,  30392879922429262840.1034288721564Q0,  &
 -143297441127511031342.841171428469Q0,  522870652282050925787.809971272961Q0, -1472390377555816831804.88187901110Q0,  &
  3168561475038540774092.35836032836Q0, -5112734757976406776681.55293208016Q0,  5985212425925539459439.17268289123Q0,  &
 -4798131130638083481859.59074371543Q0,  2355427891789631538140.18011927531Q0, -533852506199621938645.045228250685Q0/ 
Data rcutoff/ 2.1499999999999999D+00/
Data fcutoff/ 4.4447330000000002D+03/
Data r0/-3.2000000000000002D+00/
Data aa/ 2.1849245258600372D+03/
Data bb/ 0.0000000000000000D+00/
Data cc/-1.0021627430299464D+04/
Data k0/   5/
End
Double Precision Function V_Rb_6p_Pi(x)
Parameter(Npg=18)
Implicit Real*8(A-H,O-Z)
Real*16 Pg,Sto,xq,Aux
Integer*4 k0
Common/Param_LJ_V_Rb_6p_Pi/rcutoff,fcutoff,r0,aa,bb,cc,Pg(Npg),k0
If(x.le.rcutoff)Then
  V_Rb_6p_Pi=fcutoff
  Return
Endif
If(x.le.r0)Then
  V_Rb_6p_Pi=Abs(aa*x**2+bb*x+cc)
  Return
Endif
Sto=0.0q0
xq=1.0q0/x
Aux=xq**(k0+1)
Do k=1,Npg
  Sto=Sto+pg(k)*Aux
  Aux=Aux*xq
EndDo
V_Rb_6p_Pi=Sto
Return
End
!
! Potencial He(2S) del Eloranta (fci) (5-11-2015), Nist: He(2S): 329179.7623 cm-1
!
Block Data Inicio_LJ_V_He2s_fci
Implicit Real*8(A-H,O-Z)
Parameter(Npg=11)
Real*16 Pg
Integer*4 k0
Common/Param_LJ_V_He2s_fci/rcutoff,fcutoff,r0,aa,bb,cc,Pg(Npg),k0
Data Pg/                           &
 -84641112.1038843174868366388591186Q0,  5501813411.89974475880756271682097Q0, -158511366392.752661259845175521651Q0,  &
  2656659334094.79191838532512036583Q0, -28623286583317.0595767811072611642Q0,  206395570574320.545052095736748454Q0,  &
 -1005349880685171.95006561484544998Q0,  3261504889875669.31289981948240792Q0, -6745427649591000.99921588607933893Q0,  &
  8040867895208215.67859457559937689Q0, -4202197034021318.61309730280594230Q0/ 
Data rcutoff/ 3.1750200000000000D+00/
Data fcutoff/ 5.8000649999999996D+02/
Data r0/ 0.0000000000000000D+00/
Data aa/ 4.3188872751600064D+16/
Data bb/ 0.0000000000000000D+00/
Data cc/-1.3088461815357608D+17/
Data k0/   5/
End
Double Precision Function V_He2s_fci(x)
Parameter(Npg=11)
Implicit Real*8(A-H,O-Z)
Real*16 Pg,Sto,xq,Aux
Integer*4 k0
Common/Param_LJ_V_He2s_fci/rcutoff,fcutoff,r0,aa,bb,cc,Pg(Npg),k0
If(x.le.rcutoff)Then
  V_He2s_fci=fcutoff
  Return
Endif
If(x.le.r0)Then
  V_He2s_fci=Abs(aa*x**2+bb*x+cc)
  Return
Endif
Sto=0.0q0
xq=1.0q0/x
Aux=xq**(k0+1)
Do k=1,Npg
  Sto=Sto+pg(k)*Aux
  Aux=Aux*xq
EndDo
V_He2s_fci=Sto
Return
End
!
! Function for He-He+ potential, from Jussi Eloranta (6/6/2016)
!
Block Data Inicio_LJ_V_Eloranta_He2_plus_gs
Implicit Real*8(A-H,O-Z)
Parameter(Npg=17)
Real*16 Pg
Integer*4 k0
Common/Param_LJ_V_Eloranta_He2_plus_gs/rcutoff,fcutoff,r0,aa,bb,cc,Pg(Npg),k0
Data Pg/                           &
 -17123.9045474000000000000000000000Q0,  588108.342377894438906335100185956Q0, -12708070.5256971354232632580211002Q0,  &
  46593490.4538337953748571219487504Q0,  903614394.043264050274603032772751Q0, -11755047059.0947744782716618472676Q0,  &
  66544986310.3971244712772558577499Q0, -230023812497.692645654143106555337Q0,  540723165541.193441108405864043649Q0,  &
 -908192814328.760315407843038697665Q0,  1115091165095.81869466120443671210Q0, -1006768340907.62729321655030853465Q0,  &
  662532608654.360730944885641242397Q0, -309620499229.935417801050724076383Q0,  97461010879.4325864938071757129741Q0,  &
 -18542204834.6544171733101740365183Q0,  1611683836.74445477113707981227279Q0/ 
Data rcutoff/ 5.2917700000000001D-01/
Data fcutoff/ 1.1804110000000001D+04/
Data r0/ 7.3999999999999999D-01/
Data aa/-3.9258895086799792D+04/
Data bb/ 0.0000000000000000D+00/
Data cc/ 2.2797707492869155D+04/
Data k0/   3/
End
Double Precision Function V_Eloranta_He2_plus_gs(x)
Parameter(Npg=17)
Implicit Real*8(A-H,O-Z)
Real*16 Pg,Sto,xq,Aux
Integer*4 k0
Common/Param_LJ_V_Eloranta_He2_plus_gs/rcutoff,fcutoff,r0,aa,bb,cc,Pg(Npg),k0
If(x.le.rcutoff)Then
  V_Eloranta_He2_plus_gs=fcutoff
  Return
Endif
If(x.le.r0)Then
  V_Eloranta_He2_plus_gs=Abs(aa*x**2+bb*x+cc)
  Return
Endif
Sto=0.0q0
xq=1.0q0/x
Aux=xq**(k0+1)
Do k=1,Npg
  Sto=Sto+pg(k)*Aux
  Aux=Aux*xq
EndDo
V_Eloranta_He2_plus_gs=Sto
Return
End
!
!  Potencial de la molecula CH3I con lambda=0
! 
double precision function V_CH3I_0(r)
implicit none
Real (kind=8) :: Radpot, r, r_au, r_cutoff=3.6d0,  umax=3.1648655845d3
!
!     Ctes para convertir de Angs. a bohr y de atomicas a K
!
      Real  (Kind=8)  :: a_bohr=0.529177d0, Hartree=27.2114d0, eV_K=11604.448d0
If(r.Le.r_cutoff)Then
  V_CH3I_0 = Umax
Else      
  r_au     = r/a_bohr      
  V_CH3I_0 = Radpot(0,r_au)*Hartree*eV_K
EndIf
Return
end function
!
!  Potencial de la molecula CH3I con lambda=1
! 
double precision function V_CH3I_1(r)
implicit none
Real (kind=8) :: Radpot, r, r_au, r_cutoff=3.6, umax=-5.1115736616d3
!
!     Ctes para convertir de Angs. a bohr y de atomicas a K
!
      Real  (Kind=8)  :: a_bohr=0.529177d0, Hartree=27.2114d0, eV_K=11604.448d0
If(r.Le.r_cutoff)Then
  V_CH3I_1 = Umax
Else      
  r_au     = r/a_bohr      
  V_CH3I_1 = Radpot(1,r_au)*Hartree*eV_K
EndIf
Return
end function
!
!  Potencial de la molecula CH3I con lambda=2
! 
double precision function V_CH3I_2(r)
implicit none
Real (kind=8) :: Radpot, r, r_au, r_cutoff=3.6d0, umax=5.5581502684d3
!
!     Ctes para convertir de Angs. a bohr y de atomicas a K
!
      Real  (Kind=8)  :: a_bohr=0.529177d0, Hartree=27.2114d0, eV_K=11604.448d0
If(r.Le.r_cutoff)Then
  V_CH3I_2 = Umax
Else      
  r_au     = r/a_bohr      
  V_CH3I_2 = Radpot(2,r_au)*Hartree*eV_K
EndIf
Return
end function
!
!  Potencial de la molecula CH3I con lambda=3
! 
double precision function V_CH3I_3(r)
implicit none
Real (kind=8) :: Radpot, r, r_au, r_cutoff=3.6d0, umax=-2.7398540818d3
!
!     Ctes para convertir de Angs. a bohr y de atomicas a K
!
      Real  (Kind=8)  :: a_bohr=0.529177d0, Hartree=27.2114d0, eV_K=11604.448d0
If(r.Le.r_cutoff)Then
  V_CH3I_3 = Umax
Else      
  r_au     = r/a_bohr      
  V_CH3I_3 = Radpot(3,r_au)*Hartree*eV_K
EndIf
Return
end function
!
!  Potencial de la molecula CH3I con lambda=4
! 
double precision function V_CH3I_4(r)
implicit none
Real (kind=8) :: Radpot, r, r_au, r_cutoff=3.6, umax=2.4025334047d3
!
!     Ctes para convertir de Angs. a bohr y de atomicas a K
!
      Real  (Kind=8)  :: a_bohr=0.529177d0, Hartree=27.2114d0, eV_K=11604.448d0
If(r.Le.r_cutoff)Then
  V_CH3I_4 = Umax
Else      
  r_au     = r/a_bohr      
  V_CH3I_4 = Radpot(4,r_au)*Hartree*eV_K
EndIf
Return
end function
      function radpot(l,r)
      implicit none
!
!     the He-CH3I interaction potential iin grd state, only m=0
!


      integer l
      real*8  radpot,r

      real*8  d,b,c


      select case (l) 

!  fort.50 2:3
      case (0)
      d               = 0.000185445   !  +/- 2.003e-06    (1.08%)
      b               = 0.782913      !  +/- 0.005996     (0.7659%)
      c               = 9.5243        !  +/- 0.01072      (0.1125%)

!  fort.50 2:4
      case (1)
      d               = -4.51657e-05  !  +/- 5.099e-07    (1.129%)
      b               = 0.837569      !  +/- 0.004355     (0.52%)
      c               = 10.3772       !  +/- 0.007319     (0.07053%)

!  fort.50 2:5
      case (2)
      d               = 2.53688e-05   !  +/- 2.933e-07    (1.156%)
      b               = 0.849398      !  +/- 0.002787     (0.3281%)
      c               = 10.6988       !  +/- 0.004548     (0.04251%)

!  fort.50 2:6
      case (3)
      d               = -6.51825e-06  !  +/- 4.938e-08    (0.7576%)
      b               = 0.779254      !  +/- 0.0061       (0.7828%)
      c               = 11.454        !  +/- 0.007429     (0.06486%)

!  fort.50 2:7
      case (4)
      d               = 1.5032e-06    !  +/- 2.002e-08    (1.332%)
      b               = 0.787553      !  +/- 0.0054       (0.6856%)
      c               = 12.236        !  +/- 0.01018      (0.08318%)

      case default
      stop 'in radpot: l out of range'
      end select

!     Morse function
      radpot=d*(1.d0-dexp(-b*(r-c)))**2-d

      return 
      end
!
! Function for He-K(4p_Sigma), From Pascale data potential
!
Block Data Inicio_LJ_V_K_4p_Sigma
Implicit Real*8(A-H,O-Z)
Parameter(Npg=18)
Real*16 Pg
Integer*4 k0
Common/Param_LJ_V_K_4p_Sigma/rcutoff,fcutoff,r0,aa,bb,cc,Pg(Npg),k0
Data Pg/                           &
  123270229.922828519094376988757267Q0, -12085480784.0001990705905344017368Q0,  493376520440.846742309241896876082Q0,  &
 -11130615786861.2195078000310223820Q0,  155081464596787.947785296266388633Q0, -1423219373370643.41282606175556080Q0,  &
  8982444221608931.65773550693292715Q0, -40017473420474543.0863116482077444Q0,  127156098740259262.973889637297298Q0,  &
 -286289630234854071.709781745781639Q0,  443029672222162513.944773691873259Q0, -435778229290691161.879891464441776Q0,  &
  216297197302294980.110743898796384Q0, -1493759170257458.23144984879814468Q0,  3521195510768782.35743497581930372Q0,  &
 -105695005528937470.204215869914495Q0,  101576258377617829.001630054279239Q0, -30021724906855482.6250169475416685Q0/ 
Data rcutoff/ 2.2000000000000002D+00/
Data fcutoff/ 5.4944210000000003D+03/
Data r0/-2.0000000000000000D+00/
Data aa/-1.1846279611402621D+13/
Data bb/ 2.7476162593191215D+13/
Data cc/-1.5810332386853908D+13/
Data k0/   5/
End
Double Precision Function V_K_4p_Sigma(x)
Parameter(Npg=18)
Implicit Real*8(A-H,O-Z)
Real*16 Pg,Sto,xq,Aux
Integer*4 k0
Common/Param_LJ_V_K_4p_Sigma/rcutoff,fcutoff,r0,aa,bb,cc,Pg(Npg),k0
If(x.le.rcutoff)Then
  V_K_4p_Sigma=fcutoff
  Return
Endif
If(x.le.r0)Then
  V_K_4p_Sigma=Abs(aa*x**2+bb*x+cc)
  Return
Endif
Sto=0.0q0
xq=1.0q0/x
Aux=xq**(k0+1)
Do k=1,Npg
  Sto=Sto+pg(k)*Aux
  Aux=Aux*xq
EndDo
V_K_4p_Sigma=Sto
Return
End
!
! Function for He-K(4p_Pi), From Pascale data potential
!
Block Data Inicio_LJ_V_K_4p_Pi
Implicit Real*8(A-H,O-Z)
Parameter(Npg=18)
Real*16 Pg
Integer*4 k0
Common/Param_LJ_V_K_4p_Pi/rcutoff,fcutoff,r0,aa,bb,cc,Pg(Npg),k0
Data Pg/                           &
  7463823.30565045302747409754900998Q0, -718023677.708100048067350320776803Q0,  28527955672.4973076800039235370599Q0,  &
 -649780615925.820487084894875824681Q0,  9494248151960.38025074812833414839Q0, -95423287484468.9225481965073472623Q0,  &
  688611363436857.181033088432509010Q0, -3659192698303041.57349366577052330Q0,  14504208321803535.5304968565259890Q0,  &
 -43002303720081040.5938560703278609Q0,  94613503021063793.2518086792502445Q0, -151026197668471021.305996587780281Q0,  &
  166277977852662387.539989377379199Q0, -110888302469029095.889154928172435Q0,  22770818008615233.2321789074210083Q0,  &
  27284389361289459.6643521776459617Q0, -23355843002431104.0076310950131271Q0,  5879948233732488.71454983137084160Q0/ 
Data rcutoff/ 1.8999999999999999D+00/
Data fcutoff/ 5.3483040000000001D+03/
Data r0/-2.0000000000000000D+00/
Data aa/ 9.6643483830893945D+11/
Data bb/-2.2413727796818149D+12/
Data cc/ 1.2896495964123530D+12/
Data k0/   5/
End
Double Precision Function V_K_4p_Pi(x)
Parameter(Npg=18)
Implicit Real*8(A-H,O-Z)
Real*16 Pg,Sto,xq,Aux
Integer*4 k0
Common/Param_LJ_V_K_4p_Pi/rcutoff,fcutoff,r0,aa,bb,cc,Pg(Npg),k0
If(x.le.rcutoff)Then
  V_K_4p_Pi=fcutoff
  Return
Endif
If(x.le.r0)Then
  V_K_4p_Pi=Abs(aa*x**2+bb*x+cc)
  Return
Endif
Sto=0.0q0
xq=1.0q0/x
Aux=xq**(k0+1)
Do k=1,Npg
  Sto=Sto+pg(k)*Aux
  Aux=Aux*xq
EndDo
V_K_4p_Pi=Sto
Return
End
!!
!! Function for He-K(5s), From Pascale data potential (old Fit by Marti)
!!
!Block Data Inicio_LJ_V_K_5s
!Implicit Real*8(A-H,O-Z)
!Parameter(Npg=20)
!Real*16 Pg
!Integer*4 k0
!Common/Param_LJ_V_K_5s/rcutoff,fcutoff,r0,aa,bb,cc,Pg(Npg),k0
!Data Pg/                           &
!  22013374032.6756992466989907471231Q0, -2607473486972.96115845928542359422Q0,  137278559465513.504318508804350585Q0,  &
! -4284085369698061.73976060584229719Q0,  89194755667933302.9197442829881481Q0, -1322401314191165810.56518310532255Q0,  &
!  14549200062430442236.7345712007885Q0, -122132199482850453201.933729576988Q0,  797237288979113873142.546850234427Q0,  &
! -4098091530865115675119.54131140732Q0,  16712660424005787663159.1277264440Q0, -54228388406852069443691.0811371931Q0,  &
!  139769338235053275092644.042876545Q0, -284264086314955716327432.688486366Q0,  450393287366663046416997.448177291Q0,  &
! -544150565669437497681315.135493803Q0,  484098752094243106107299.506791578Q0, -298854633397265157385967.162411014Q0,  &
!  114324971752633317584036.132040739Q0, -20408914798270344038969.9560381946Q0/ 
!Data rcutoff/ 2.0000000000000000D+00/
!Data fcutoff/ 2.7298080000000000D+03/
!Data r0/-2.0000000000000000D+00/
!Data aa/-1.1118034495930548D+05/
!Data bb/ 4.7153152909318462D+05/
!Data cc/-4.9833086495060596D+05/
!Data k0/   5/
!End
!Double Precision Function V_K_5s(x)
!Parameter(Npg=20)
!Implicit Real*8(A-H,O-Z)
!Real*16 Pg,Sto,xq,Aux
!Integer*4 k0
!Common/Param_LJ_V_K_5s/rcutoff,fcutoff,r0,aa,bb,cc,Pg(Npg),k0
!If(x.le.rcutoff)Then
!  V_K_5s=fcutoff
!  Return
!Endif
!If(x.le.r0)Then
!  V_K_5s=Abs(aa*x**2+bb*x+cc)
!  Return
!Endif
!Sto=0.0q0
!xq=1.0q0/x
!Aux=xq**(k0+1)
!Do k=1,Npg
!  Sto=Sto+pg(k)*Aux
!  Aux=Aux*xq
!EndDo
!V_K_5s=Sto
!Return
!End
!
! Function for He-K(5s), From Pascale data potential (new fit from Maxime Martinez)
!
Block Data Inicio_LJ_V_K_5s
Implicit Real*8(A-H,O-Z)
Parameter(Npg=20)
Real*16 Pg
Integer*4 k0
Common/Param_LJ_V_K_5s/rcutoff,fcutoff,r0,aa,bb,cc,Pg(Npg),k0
Data Pg/                           &
 -1168397771.45542931556701660156250Q0,  118762977694.007186889648437500000Q0, -5067632180013.47558593750000000000Q0,  &
  119820387630240.453125000000000000Q0, -1770902027637011.75000000000000000Q0,  17565225586361198.0000000000000000Q0,  &
 -122078985717169440.000000000000000Q0,  608590024952670464.000000000000000Q0, -2190107282135850496.00000000000000Q0,  &
  5613121256762224640.00000000000000Q0, -9769037125165985792.00000000000000Q0,  9996223050685026304.00000000000000Q0,  &
 -2343730974240809984.00000000000000Q0, -7106519104385017856.00000000000000Q0,  5621562826636969984.00000000000000Q0,  &
  4921614491484340224.00000000000000Q0, -6352243461664580608.00000000000000Q0, -4076025095303220224.00000000000000Q0,  &
  8675010493626412032.00000000000000Q0, -3523847995886618112.00000000000000Q0/
Data rcutoff/ 1.980000000000000D+00/
Data fcutoff/ 2.946118600000D+03/
Data r0/-2.0000000000000000D+00/
Data aa/-1.1118034495930548D+05/
Data bb/ 4.7153152909318462D+05/
Data cc/-4.9833086495060596D+05/
Data k0/   5/
End
Double Precision Function V_K_5s(x)
Parameter(Npg=20)
Implicit Real*8(A-H,O-Z)
Real*16 Pg,Sto,xq,Aux
Integer*4 k0
Common/Param_LJ_V_K_5s/rcutoff,fcutoff,r0,aa,bb,cc,Pg(Npg),k0
If(x.le.rcutoff)Then
  V_K_5s=fcutoff
  Return
Endif
If(x.le.r0)Then
  V_K_5s=Abs(aa*x**2+bb*x+cc)
  Return
Endif
Sto=0.0q0
xq=1.0q0/x
Aux=xq**(k0+1)
Do k=1,Npg
  Sto=Sto+pg(k)*Aux
  Aux=Aux*xq
EndDo
V_K_5s=Sto
Return
End
!....................
!.. Funcion for Ca ..
!....................
double precision function V_ca(rr)

!...ground state potential for Ca-He. Mauschick and Meyer 1989

implicit none

real (kind=8),    parameter :: bohr    = 0.529177249  ! 1Bohr = 0.529 \AA
integer (kind=4), parameter :: mlr=5
integer (kind=4), parameter :: mp =6
real    (kind=8)            :: clr(mlr),p(mp)
real    (kind=8)            :: r,rr
integer (kind=4)            :: k,m

data p/ 9.397738d-6,1.233746d0,0.06802583d0,    &
        0.4441812d0,1.641748d0,4.73051d0/

data clr/36.404,2.0834d+3,1.2405d+5,7.68460D+6,4.95275D+8/ 

interface  
function gam(l,x)
   integer (kind=4) :: l
   real    (kind=8) :: x
   real    (kind=8) :: gam
end function gam
end interface
r=rr/bohr
V_ca = p(1)*(1+p(3)*r**p(4))*exp(-p(2)*(r-11)) ! Repulsion

do k=1,mlr
  m    = 4+2*k
  V_ca = V_ca-clr(k)/r**m*gam(m,p(5)*(r-p(6))) ! Attraction
end do
V_ca = V_ca/3.1667d-6
return

end

!...................
!.. Function gam ...
!...................

double precision function gam(l,x)
implicit none

integer (kind=4) :: i,l
real    (kind=8) :: x,xx


gam = 0.0d0
if(x.lt.0.d0) return
gam = 1.0d0
xx  = 1.0d0

do i=1,l
  xx  = xx*x/i
  gam = gam+xx
end do
gam = 1.0d0-gam*exp(-x)

return
end 
!
!  Funcio que calcula el potencial Ar-Ar, entrada \AA, sortida K
!
!  r_cutoff = 2.0d0; umax = 8.0808095081d3
!
!
Double precision function V_Ar_Ar(r)
Implicit none
real      (kind=8) :: AA     = 748.3d0
real      (kind=8) :: BB     = 2.031d0
real      (kind=8) :: C6     = 64.3d0
real      (kind=8) :: C8     = 1623.d0
real      (kind=8) :: C10    = 49060.d0 
real      (kind=8) :: rbohr  = 0.5292d0
real      (kind=8) :: Hartree_to_K  = 27.21d0*11604.d0
real (kind=8) :: r, rr, Stoe, Stor, F6, F8, F10, Sto, Aux
real (kind=8) :: c1or2, c1or6
real (kind=8) :: r_cutoff=2.5d0, f_cutoff=1.224418d4
real (kind=8) :: Fact(0:10)
Logical  :: Lfirst = .true.
Integer (Kind=4) :: i
Save     :: Lfirst, Fact
if(r.Le.r_cutoff)Then
  V_Ar_Ar = f_cutoff
  Return
Endif        
If(Lfirst)Then
!  Write(6,'("#Control per veure si nomes hi pasem una vegada")')
  Fact(0) = 1.d0; Fact(1) = 1.d0
  Do i=2, 10
    Fact(i) = Fact(i-1)*i
  EndDo
  Lfirst = .false.
Endif
rr = r/rbohr
c1or2 = 1.d0/(rr*rr)
c1or6 = c1or2*c1or2*c1or2
Stor = bb*rr
Stoe = Exp(-Stor)
Sto  = 1.d0
Aux = Stor 
Do i=1, 6
  Sto = Sto + Aux/Fact(i)
  Aux = Aux*Stor
Enddo
F6 = 1.d0 - Stoe*Sto
Do i=7, 8
  Sto = Sto + Aux/Fact(i)
  Aux = Aux*Stor
Enddo
F8 = 1.d0 - Stoe*Sto
Do i=9, 10
  Sto = Sto + Aux/Fact(i)
  Aux = Aux*Stor
Enddo
F10 = 1.d0 - Stoe*Sto
Aux = AA*Stoe - (C6*F6 + (C8*F8 + C10*F10*c1or2)*c1or2)*c1or6
V_Ar_Ar = Aux*Hartree_to_K
Return
End
!
!  Funcio que calcula el potencial Xe-Xe, entrada \AA, sortida K
!
!
Double precision function V_Xe_Xe(r)
Implicit none
real      (kind=8) :: AA     = 951.8d0
real      (kind=8) :: BB     = 1.681d0
real      (kind=8) :: C6     = 285.9d0
real      (kind=8) :: C8     = 12810.d0
real      (kind=8) :: C10    = 619800.d0 
real      (kind=8) :: rbohr  = 0.5292d0
real      (kind=8) :: Hartree_to_K  = 27.21d0*11604.d0
real (kind=8) :: r, rr, Stoe, Stor, F6, F8, F10, Sto, Aux
real (kind=8) :: c1or2, c1or6
real (kind=8) :: r_cutoff=2.5d0, f_cutoff=8.858633d4
real (kind=8) :: Fact(0:10)
Logical  :: Lfirst = .true.
Integer (Kind=4) :: i
Save     :: Lfirst, Fact
If(r.Le.r_cutoff)Then
  V_Xe_Xe = f_cutoff
  Return
EndIf  
If(Lfirst)Then
!  Write(6,'("#Control per veure si nomes hi pasem una vegada")')
  Fact(0) = 1.d0; Fact(1) = 1.d0
  Do i=2, 10
    Fact(i) = Fact(i-1)*i
  EndDo
  Lfirst = .false.
Endif
rr = r/rbohr
c1or2 = 1.d0/(rr*rr)
c1or6 = c1or2*c1or2*c1or2
Stor = bb*rr
Stoe = Exp(-Stor)
Sto  = 1.d0
Aux = Stor 
Do i=1, 6
  Sto = Sto + Aux/Fact(i)
  Aux = Aux*Stor
Enddo
F6 = 1.d0 - Stoe*Sto
Do i=7, 8
  Sto = Sto + Aux/Fact(i)
  Aux = Aux*Stor
Enddo
F8 = 1.d0 - Stoe*Sto
Do i=9, 10
  Sto = Sto + Aux/Fact(i)
  Aux = Aux*Stor
Enddo
F10 = 1.d0 - Stoe*Sto
Aux = AA*Stoe - (C6*F6 + (C8*F8 + C10*F10*c1or2)*c1or2)*c1or6
V_Xe_Xe = Aux*Hartree_to_K
Return
End
!
!  Funcio que calcula la derivada del potencial Ar-Ar, entrada \AA, sortida K
!
!  r_cutoff = 2.0d0; umax = 8.0808095081d3
!
!
Double precision function drV_Ar_Ar(r)
Implicit none
real      (kind=8) :: AA     = 748.3d0
real      (kind=8) :: BB     = 2.031d0
real      (kind=8) :: C6     = 64.3d0
real      (kind=8) :: C8     = 1623.d0
real      (kind=8) :: C10    = 49060.d0 
real      (kind=8) :: rbohr  = 0.5292d0
real      (kind=8) :: Hartree_to_K  = 27.21d0*11604.d0
real (kind=8) :: r, rr, Stoe, Stor, F6, F8, F10, Sto, Aux
real (kind=8) :: c1or, c1or2, c1or7, b2, b6, b8, b10
real (kind=8) :: r_cutoff=2.5d0, f_cutoff=-5.266478d4
real (kind=8) :: Fact(0:10)
Logical  :: Lfirst = .true.
Integer (Kind=4) :: i
Save     :: Lfirst, Fact
If(r.Le.r_cutoff)Then
  drV_Ar_Ar = f_cutoff
  Return
EndIf  
If(Lfirst)Then
!  Write(6,'("#Control per veure si nomes hi pasem una vegada")')
  Fact(0) = 1.d0; Fact(1) = 1.d0
  Do i=2, 10
    Fact(i) = Fact(i-1)*i
  EndDo
  Lfirst = .false.
Endif
b2=bb*bb
b6=b2*b2*b2
b8=b6*b2
b10=b8*b2
rr = r/rbohr
c1or = 1.d0/rr
c1or2 = c1or*c1or
c1or7 = c1or2*c1or2*c1or2*c1or
Stor = bb*rr
Stoe = Exp(-Stor)
Sto  = 1.d0
Aux = Stor 
Do i=1, 6
  Sto = Sto + Aux/Fact(i)
  Aux = Aux*Stor
Enddo
F6 = 1.d0 - Stoe*Sto
Do i=7, 8
  Sto = Sto + Aux/Fact(i)
  Aux = Aux*Stor
Enddo
F8 = 1.d0 - Stoe*Sto
Do i=9, 10
  Sto = Sto + Aux/Fact(i)
  Aux = Aux*Stor
Enddo
F10 = 1.d0 - Stoe*Sto
Aux = -bb*AA*Stoe + (6.d0*C6*F6 + (8.d0*C8*F8 + 10.d0*C10*F10*c1or2)*c1or2)*c1or7      &
    -  bb*Stoe*(b6*C6/Fact(6) + b8*C8/Fact(8) + b10*C10/fact(10))
drV_Ar_Ar = Aux*Hartree_to_K/rbohr
Return
End
!
!  Funcio que calcula la derivada del potencial Ar-Ar, entrada \AA, sortida K
!
!
Double precision function drV_Xe_Xe(r)
Implicit none
real      (kind=8) :: AA     = 951.8d0
real      (kind=8) :: BB     = 1.681d0
real      (kind=8) :: C6     = 285.9d0
real      (kind=8) :: C8     = 12810.d0
real      (kind=8) :: C10    = 619800.d0 
real      (kind=8) :: rbohr  = 0.5292d0
real      (kind=8) :: Hartree_to_K  = 27.21d0*11604.d0
real (kind=8) :: r, rr, Stoe, Stor, F6, F8, F10, Sto, Aux
real (kind=8) :: c1or, c1or2, c1or7, b2, b6, b8, b10
real (kind=8) :: r_cutoff=2.5d0, f_cutoff=-3.017394d5
real (kind=8) :: Fact(0:10)
Logical  :: Lfirst = .true.
Integer (Kind=4) :: i
Save     :: Lfirst, Fact
If(r.Le.r_cutoff)Then
  drV_Xe_Xe = f_cutoff
  Return
EndIf  
If(Lfirst)Then
!  Write(6,'("#Control per veure si nomes hi pasem una vegada")')
  Fact(0) = 1.d0; Fact(1) = 1.d0
  Do i=2, 10
    Fact(i) = Fact(i-1)*i
  EndDo
  Lfirst = .false.
Endif
b2=bb*bb
b6=b2*b2*b2
b8=b6*b2
b10=b8*b2
rr = r/rbohr
c1or = 1.d0/rr
c1or2 = c1or*c1or
c1or7 = c1or2*c1or2*c1or2*c1or
Stor = bb*rr
Stoe = Exp(-Stor)
Sto  = 1.d0
Aux = Stor 
Do i=1, 6
  Sto = Sto + Aux/Fact(i)
  Aux = Aux*Stor
Enddo
F6 = 1.d0 - Stoe*Sto
Do i=7, 8
  Sto = Sto + Aux/Fact(i)
  Aux = Aux*Stor
Enddo
F8 = 1.d0 - Stoe*Sto
Do i=9, 10
  Sto = Sto + Aux/Fact(i)
  Aux = Aux*Stor
Enddo
F10 = 1.d0 - Stoe*Sto
Aux = -bb*AA*Stoe + (6.d0*C6*F6 + (8.d0*C8*F8 + 10.d0*C10*F10*c1or2)*c1or2)*c1or7      &
    -  bb*Stoe*(b6*C6/Fact(6) + b8*C8/Fact(8) + b10*C10/fact(10))
drV_Xe_Xe = Aux*Hartree_to_K/rbohr
Return
End
!
! Function for He-Sr(5s2), From Lovallo data potential
!
Block Data Inicio_LJ_V_Sr_He
Implicit Real*8(A-H,O-Z)
Parameter(N=16)
Integer (Kind=4) :: k0
Common/Param_Spline_LJ_V_Sr_He/ak,a12,r0,rf,aa,cc,X(N),A(N),B(N),C(N),D(N),k0
Data X/                                                                                                   &
  2.00000000000000D0,  3.00000000000000D0,  4.00000000000000D0,  5.00000000000000D0,  5.30000000000000D0, &
  6.30000000000000D0,  6.40000000000000D0,  6.50000000000000D0,  6.60000000000000D0,  6.70000000000000D0, &
  6.80000000000000D0,  7.50000000000000D0,  8.00000000000000D0,  9.00000000000000D0,  10.0000000000000D0, &
  12.0000000000000D0/ 
Data A/                                                                                                   &
  19682.2373956429D0,  135538.884184861D0, -28755.0411506975D0,  82682.1648908597D0, -1028.21398783985D0, &
  12723.4933807024D0, -3297.73585567146D0,  5511.51880624365D0, -1397.11955900432D0,  555.849036582276D0, &
  824.670054265175D0, -628.064864640832D0, -6.25389576350048D0, -82.9390390800452D0, -21.4578500114020D0, &
  0.00000000000000D0/ 
Data B/                                                                                                   &
  12722.4786274321D0, -103134.168161786D0,  20086.2758398830D0, -46776.0477850513D0,  607.185542514463D0, &
 -5941.24653774376D0,  1568.70466680649D0, -2497.10517715432D0,  643.184988867479D0, -231.278561395176D0, &
 -349.876069196455D0,  231.217898365948D0, -1.96121496305145D0,  23.6004994757968D0,  5.15614275520381D0, &
  0.00000000000000D0/ 
Data C/                                                                                                   &
 -12360.7671325954D0,  26258.1151304771D0, -4546.99586994003D0,  8825.46885504684D0, -114.763848267464D0, &
  924.669815265587D0, -248.760060445390D0,  376.749146317813D0, -99.0523939885207D0,  31.4645538118754D0, &
  48.9053637826517D0, -28.5738318923353D0, 0.573557273789606D0, -2.26663321941575D0,-0.422197547356456D0, &
  0.00000000000000D0/ 
Data D/                                                                                                   &
  2060.12785543257D0, -2230.85906268660D0,  336.233520681496D0, -555.264127650962D0,  7.01465872100645D0, &
 -47.9818314130174D0,  13.1343079469293D0, -18.9430872716965D0,  5.08729355185571D0, -1.40608693572618D0, &
 -2.26102860096032D0,  1.18249120681688D0,-3.198334177165370D0, 7.320889871743369D0, 1.172770964879045D0, &
  0.00000000000000D0/ 
Data ak/-4.6914394784477353D+05/
Data a12/ 1.2247267783468750D+10/
Data rf/ 8.0000000000000000D+00/
Data r0/ 2.0000000000000000D+00/
Data aa/-2.9997639094396982D+03/
Data cc/ 2.4164204601344718D+04/
Data k0/   6/
End
Double Precision Function V_Sr_He(xx)
Parameter(N=16)
Implicit Real*8(A-H,O-Z)
Integer*4 k0
Common/Param_Spline_LJ_V_Sr_He/ak,a12,r0,rf,aa,cc,X(N),A(N),B(N),C(N),D(N),k0
If(xx.le.r0)Then
  V_Sr_He=aa*xx**2+cc
  Return
Endif
If(xx.gt.rf)Then
  V_Sr_He=ak/xx**k0 + a12/xx**12
  Return
Endif
Call Findi(X,xx,N,ix)
V_Sr_He= a(ix) + (b(ix) + xx*(c(ix) + xx*d(ix)))*xx
Return
End
double precision function V_Kerkines_Li_gs(r)
implicit none
real      (kind=8) :: A      = 1.6725d4
real      (kind=8) :: alpha  = 0.42685d0
real      (kind=8) :: beta   = 0.19706d0
real      (kind=8) :: D      = 13.988d0
real      (kind=8) :: C6     = 1.5992d5
!real      (kind=8) :: C8     = 0.42123807d0
real      (kind=8) :: C12    = 3.2761d10
real (kind=8) :: r,ff
 If(r==0)then
   V_Kerkines_Li_gs = A
   Return
 EndIf  
 ff=1.d0
 If(r.le.D)ff=dexp(-(D/r-1.d0)**2)
 V_Kerkines_Li_gs = A*dexp( -alpha*r -beta*r**2) - ff*(C6/r**6 + C12/r**12)
 Return
end function
double precision function V_Kerkines_Li_Pi(r)
implicit none
real      (kind=8) :: A      = 6.2669d5
real      (kind=8) :: alpha  = 0.d0
real      (kind=8) :: beta   = 2.8981d0
real      (kind=8) :: D      = 4.5734d0
real      (kind=8) :: C6     = 1.9521d5
!real      (kind=8) :: C8     = 0.42123807d0
real      (kind=8) :: C10    = 3.7221d6
!real      (kind=8) :: C12    = 1.7527d12
real (kind=8) :: r,ff
 If(r==0)then
   V_Kerkines_Li_Pi = A
   Return
 EndIf  
 ff=1.d0
 If(r.le.D)ff=dexp(-(D/r-1.d0)**2)
 V_Kerkines_Li_Pi = A*dexp( -Beta*r**2) - ff*(C6/r**6 + C10/r**10)
 Return
end function
double precision function V_Kerkines_Li_Sigma(r)
implicit none
real      (kind=8) :: A      = 8.6936d5
real      (kind=8) :: alpha  = 1.6978d0
real      (kind=8) :: beta   = 0.d0
real      (kind=8) :: D      = 21.456d0
!real      (kind=8) :: C6     = 1.5992d5
!real      (kind=8) :: C8     = 0.42123807d0
real      (kind=8) :: C10    = 7.0088d3
real      (kind=8) :: C12    = 1.7527d12
real (kind=8) :: r,ff
 If(r==0)then
   V_Kerkines_Li_Sigma = A
   Return
 EndIf  
 ff=1.d0
 If(r.le.D)ff=dexp(-(D/r-1.d0)**2)
 V_Kerkines_Li_Sigma = A*dexp( -alpha*r) - ff*(C10/r**10 + C12/r**12)
 Return
end function
!
! Function for He-Na(3p0), From Pascale data potential
!
Block Data Inicio_LJ_V_Na_Sigma
Implicit Real*8(A-H,O-Z)
Parameter(N=34)
Integer (Kind=4) :: k0
Common/Param_Spline_LJ_V_Na_Sigma/ak,a12,r0,rf,aa,cc,X(N),A(N),B(N),C(N),D(N),k0
Data X/                                                                                                                      &
  1.058354000000000D+00,   1.322942500000000D+00,   1.587531000000000D+00,   1.852119500000000D+00,   2.116708000000000D+00, &
  2.381296500000000D+00,   2.645885000000000D+00,   2.910473500000000D+00,   3.175062000000000D+00,   3.439650500000000D+00, &
  3.704239000000000D+00,   3.968827500000000D+00,   4.233416000000000D+00,   4.762593000000000D+00,   5.291770000000000D+00, &
  5.820947000000000D+00,   6.350124000000000D+00,   6.879301000000000D+00,   7.408478000000000D+00,   7.937655000000000D+00, &
  8.466832000000000D+00,   8.996009000000001D+00,   9.525186000000000D+00,   1.005436300000000D+01,   1.058354000000000D+01, &
  1.190648250000000D+01,   1.322942500000000D+01,   1.455236750000000D+01,   1.587531000000000D+01,   1.719825250000000D+01, &
  1.852119500000000D+01,   2.116708000000000D+01,   2.381296500000000D+01,   2.645885000000000D+01  /
Data A/                                                                                                                      &
 -3.890915094145514D+05,   2.842048955360034D+06,  -9.387352154073080D+04,   4.860945168944493D+05,  -3.682190378428085D+04, &
  3.674455587895416D+04,  -1.721438688442292D+04,   2.567156632954004D+03,   1.303065248252341D+04,   1.997174092036358D+04, &
  2.697079033034076D+04,   2.768983007896994D+04,   2.813607007945593D+04,   2.551325162844106D+04,   2.110508160272032D+04, &
  1.557964097612617D+04,   1.047575045279671D+04,   6.640191937636850D+03,   4.126079405769719D+03,   2.436647851185728D+03, &
  1.417993363266505D+03,   8.060063945920606D+02,   3.925260172221136D+02,   2.248211095194046D+02,   3.369649265259418D+01, &
 -2.516352812110249D+01,  -1.056065295577870D+01,  -1.446393496628078D+01,   2.398632234516816D+00,  -7.540154642234392D+00, &
  3.404083314497583D-02,  -1.285427052711144D+00,   3.352404111107370D-01,   0.000000000000000D+00  /
Data B/                                                                                                                      &
  2.038171376789891D+06,  -5.288996201788737D+06,   2.590954138724962D+05,  -6.803170351112438D+05,   6.080987612837456D+04, &
 -3.187046795717086D+04,   2.931013070457417D+04,   8.920104630764952D+03,  -9.664700404409080D+02,  -7.020359909054157D+03, &
 -1.268877067572756D+04,  -1.323228617139665D+04,  -1.354851306746780D+04,  -1.189637622666603D+04,  -9.397305012844907D+03, &
 -6.549603104584538D+03,  -4.138363328796164D+03,  -2.465711478632961D+03,  -1.447643044657535D+03,  -8.091302012860393D+02, &
 -4.481966828510839D+02,  -2.441105257537164D+02,  -1.138830286623377D+02,  -6.384358577489398D+01,  -9.667586950263340D+00, &
  5.162994820724863D+00,   1.851541337601345D+00,   2.656210820522679D+00,  -5.303538262396706D-01,   1.203331653157357D+00, &
 -2.351069832904046D-02,   1.634968462906406D-01,  -4.067774479726644D-02,   0.000000000000000D+00  /
Data C/                                                                                                                      &
 -2.228002766666602D+06,   3.310535437812148D+06,  -1.842571776777128D+05,   3.229522377966883D+05,  -2.717962320616065D+04, &
  1.174049618491284D+04,  -1.138242815297223D+04,  -4.376686278390921D+03,  -1.262865297507440D+03,   4.971652255393781D+02, &
  2.027414965546275D+03,   2.164361077633301D+03,   2.239058885755905D+03,   1.892160286442092D+03,   1.419904082219095D+03, &
  9.306878244073738D+02,   5.509724399852153D+02,   3.078297354919530D+02,   1.704103581278810D+02,   8.996936599966477D+01, &
  4.734025531990225D+01,   2.465395552880960D+01,   1.098204433553960D+01,   6.005155905368167D+00,   8.862638499093484D-01, &
 -3.593249938980653D-01,  -1.090152878357445D-01,  -1.643100350939714D-01,   3.641452693110717D-02,  -6.439138223885225D-02, &
  1.848531086740251D-03,  -6.986298498619706D-03,   1.587795113322655D-03,   0.000000000000000D+00  /
Data D/                                                                                                                      &
  7.017194519875840D+05,  -6.937904824902311D+05,   4.000991059182749D+04,  -5.127459084457212D+04,   3.863208197281073D+03, &
 -1.584821113792885D+03,   1.328246327143395D+03,   5.258864047709351D+02,   1.989824492068674D+02,   2.841952689379911D+01, &
 -1.092830475537368D+02,  -1.207848586122947D+02,  -1.266664609339661D+02,  -1.023870680826480D+02,  -7.263922667211827D+01, &
 -4.462452632572438D+01,  -2.469234944108620D+01,  -1.291098267658196D+01,  -6.728004321060446D+00,  -3.349970878848767D+00, &
 -1.671692207449495D+00,  -8.310861197539162D-01,  -3.526383780415761D-01,  -1.876390827050994D-01,  -2.641728434466360D-02, &
  8.454163311614218D-03,   2.147270254097337D-03,   3.413839588565351D-03,  -8.007670768599447D-04,   1.153034317980596D-03, &
 -3.911182099860281D-05,   1.000166367348203D-04,  -2.000332734696399D-05,   0.000000000000000D+00  /
Data ak/-6.7791488520017383D+05/
Data a12/-1.2075208543818171D+12/
Data rf/ 1.5000000000000000D+01/
Data r0/ 1.3000000000000000D+00/
Data aa/-7.5737767294701451D+04/
Data cc/ 1.6488106749051623D+05/
Data k0/   6/
End
Double Precision Function V_Na_Sigma(xx)
Parameter(N=34)
Implicit Real*8(A-H,O-Z)
Integer*4 k0
Common/Param_Spline_LJ_V_Na_Sigma/ak,a12,r0,rf,aa,cc,X(N),A(N),B(N),C(N),D(N),k0
If(xx.le.r0)Then
  V_Na_Sigma=aa*xx**2+cc
  Return
Endif
If(xx.gt.rf)Then
  V_Na_Sigma=ak/xx**k0 + a12/xx**12
  Return
Endif
Call Findi(X,xx,N,ix)
V_Na_Sigma= a(ix) + (b(ix) + xx*(c(ix) + xx*d(ix)))*xx
Return
End
!
! Function for He-Na(3p1), From Pascale data potential
!
Block Data Inicio_LJ_V_Na_Pi
Implicit Real*8(A-H,O-Z)
Parameter(N=34)
Integer (Kind=4) :: k0
Common/Param_Spline_LJ_V_Na_Pi/ak,a12,r0,rf,aa,cc,X(N),A(N),B(N),C(N),D(N),k0
Data X/                                                                                                                      &
  1.058354000000000D+00,   1.322942500000000D+00,   1.587531000000000D+00,   1.852119500000000D+00,   2.116708000000000D+00, &
  2.381296500000000D+00,   2.645885000000000D+00,   2.910473500000000D+00,   3.175062000000000D+00,   3.439650500000000D+00, &
  3.704239000000000D+00,   3.968827500000000D+00,   4.233416000000000D+00,   4.762593000000000D+00,   5.291770000000000D+00, &
  5.820947000000000D+00,   6.350124000000000D+00,   6.879301000000000D+00,   7.408478000000000D+00,   7.937655000000000D+00, &
  8.466832000000000D+00,   8.996009000000001D+00,   9.525186000000000D+00,   1.005436300000000D+01,   1.058354000000000D+01, &
  1.190648250000000D+01,   1.322942500000000D+01,   1.455236750000000D+01,   1.587531000000000D+01,   1.719825250000000D+01, &
  1.852119500000000D+01,   2.116708000000000D+01,   2.381296500000000D+01,   2.645885000000000D+01  /
Data A/                                                                                                                      &
 -3.829612265209504D+05,   2.858488533719470D+06,  -8.850505962450517D+04,   4.985161221791919D+05,  -9.374435653208049D+03, &
  7.506304846840330D+04,   1.048141319844631D+04,   2.389023908616203D+04,   1.761808333336876D+04,   1.864857334491393D+04, &
  1.600855163427982D+04,   1.423083304231262D+04,   1.009103672139864D+04,   5.509370743020555D+03,   2.782969659106633D+03, &
  1.139314871927326D+03,   3.791845333915003D+02,   5.516886143493618D+01,  -3.906518534900924D+01,  -6.883528596251043D+01, &
 -6.820116280662592D+01,  -4.342108704388850D+01,  -3.307391379617849D+01,  -2.971075725994607D+01,  -1.441692751124337D+01, &
 -5.637482470438435D+00,  -4.645791803909031D-01,  -6.208047560804336D+00,   3.766387287683818D+00,  -2.791993659373083D+00, &
  1.963023548370577D-01,  -6.934189615083726D-01,   1.987803926698458D-01,   0.000000000000000D+00  /
Data B/                                                                                                                      &
  2.042259675306722D+06,  -5.308286006551154D+06,   2.607269954197847D+05,  -6.901098946464834D+05,   2.972093393128923D+04, &
 -7.665492155954697D+04,  -3.429930371392179D+03,  -1.725121328743561D+04,  -1.132488263361537D+04,  -1.222365709765917D+04, &
 -1.008554853287642D+04,  -8.741787491132578D+03,  -5.808130850079128D+03,  -2.922098401970659D+03,  -1.376452379611076D+03, &
 -5.293454807614694D+02,  -1.702357979257522D+02,  -2.893553851426660D+01,   9.223735273390762D+00,   2.047520713006675D+01, &
  2.025052221041162D+01,   1.198683246885931D+01,   8.727943902995086D+00,   7.724452223940203D+00,   3.389278052905926D+00, &
  1.177177617499013D+00,   4.132691499388059D-03,   1.188160317330051D+00,  -6.967354450496733D-01,   4.472844395967870D-01, &
 -3.674955727168822D-02,   8.935021412031723D-02,  -2.305065920791614D-02,   0.000000000000000D+00  /
Data C/                                                                                                                      &
 -2.238951111328112D+06,   3.317258385197909D+06,  -1.907127989685347D+05,   3.226649231850672D+05,  -1.740599286843306D+04, &
  2.726541007985003D+04,  -4.096367147592592D+02,   4.339172341926340D+03,   2.472647954742643D+03,   2.733946148859692D+03, &
  2.156740259935474D+03,   1.818161412216886D+03,   1.125185186621894D+03,   5.192059546975057D+02,   2.271210714997131D+02, &
  8.159373739933278D+01,   2.504213575569288D+01,   4.502220521356484D+00,  -6.485370009925553D-01,  -2.066017586464712D+00, &
 -2.039480515733283D+00,  -1.120885420780364D+00,  -7.787515699700785D-01,  -6.789449806262401D-01,  -2.693302230844084D-01, &
 -8.354063867886183D-02,   5.128742341071029D-03,  -7.623448778836754D-02,   4.249678501069743D-02,  -2.402275728852516D-02, &
  2.111301386872158D-03,  -3.846051795144703D-03,   8.741026807147035D-04,   0.000000000000000D+00  /
Data D/                                                                                                                      &
  7.051676821832495D+05,  -6.947947743674628D+05,   4.177271820778849D+04,  -5.062191798309119D+04,   2.931511218212643D+03, &
 -3.321581435065812D+03,   1.649705496966011D+02,  -3.789053358692815D+02,  -1.829485968424848D+02,  -2.082707619982226D+02, &
 -1.563297390310924D+02,  -1.278932257519952D+02,  -7.332923412010020D+01,  -3.091681775193069D+01,  -1.251812927291563D+01, &
 -4.184586403208442D+00,  -1.216056358564700D+00,  -2.208043102177145D-01,   1.094628337470038D-02,   7.047186473225737D-02, &
  6.942711856937026D-02,   3.538998443752289D-02,   2.341703011025180D-02,   2.010813196179497D-02,   7.207131680797889D-03, &
  2.005775069513611D-03,  -2.283704311394874D-04,   1.635317845686456D-03,  -8.576785703219029D-04,   4.315908346020876D-04, &
 -3.875424258105890D-05,   5.506051098937314D-05,  -1.101210219787461D-05,   0.000000000000000D+00  /
Data ak/-2.9611722908468673D+05/
Data a12/ 8.8120472293789087D+11/
Data rf/ 1.5000000000000000D+01/
Data r0/ 1.3000000000000000D+00/
Data aa/-7.8389640568190327D+04/
Data cc/ 1.6988086355011899D+05/
Data k0/   6/
End
Double Precision Function V_Na_Pi(xx)
Parameter(N=34)
Implicit Real*8(A-H,O-Z)
Integer*4 k0
Common/Param_Spline_LJ_V_Na_Pi/ak,a12,r0,rf,aa,cc,X(N),A(N),B(N),C(N),D(N),k0
If(xx.le.r0)Then
  V_Na_Pi=aa*xx**2+cc
  Return
Endif
If(xx.gt.rf)Then
  V_Na_Pi=ak/xx**k0 + a12/xx**12
  Return
Endif
Call Findi(X,xx,N,ix)
V_Na_Pi= a(ix) + (b(ix) + xx*(c(ix) + xx*d(ix)))*xx
Return
End
