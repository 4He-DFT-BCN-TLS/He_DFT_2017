 &input
  Lbackflow=.true.
  Min_rho=1.d-5
  Max_rho=3.5d-2
  Lvortex=.true.
  Vortex_Axis='Z'
  nv=1
  precie=1.d-5
  ximp = 0.d0
  yimp = 0.d0
  zimp = 26.55d0
  xc=0.0d0
  yc=0.0d0
  zc=0.0d0
  gwf=1.
  afermi=1.5d0
  rimpur=5.0d0
  limp=.true.
  lexternal=.false.
  Selec_gs='NA', r_cutoff_gs=2.0d0,  umax_gs=6.7131669574d3
  lsolid=.false.
  core4='OTC'
  nq=1500
  leepot='NO '
  mimpur=22.98977d0
  nthread=4
  title   = '4He_1000 + Na'
  n4      = 1000
  mode    = 0
  filedenin  = 'den.inp'
  fileimpin  = 'denx.inp'
  filedenout = 'den.out'
  fileimpout = 'denx.out'
  fileout = 'OT-FFT-128x128x128'
  nsfiles=1
  nx = 192
  ny = 192
  nz = 192
  xmax = 35.d0
  ymax = 35.d0
  zmax = 35.d0
  fftwplan="FFTW_PATIENT"
  npd=13
  niter=6000000
  icon=13
  pener=200
  pchem=50
  ppot=1
  pdenpar=500
  denmin=-1
  psimin=-5
  iron=1
  ironx=1
  pafl=0.01
  pafl=0.05
  pafl=0.01
  lpaflv=.false.
&End
100  0.01
200  0.05
400  0.10
1000 0.20
