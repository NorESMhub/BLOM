import xarray as xr

def init_parameters(units):
  global a11
  global a12
  global a13
  global a14
  global a15
  global a16
  global b11
  global b12
  global b13
  global a21
  global a22
  global a23
  global a24
  global a25
  global a26
  global b21
  global b22
  global b23


  match units:
    case 'MKS':
      a11 =  9.9985372432159340e+02
      a12 =  1.0380621928183473e+01
      a13 =  1.7073577195684715e+00
      a14 = -3.6570490496333680e-02
      a15 = -7.3677944503527477e-03
      a16 = -3.5529175999643348e-03
      b11 =  1.7083494994335439e-06
      b12 =  7.1567921402953455e-09
      b13 =  1.2821026080049485e-09
      a21 =  1.
      a22 =  1.0316374535350838e-02
      a23 =  8.9521792365142522e-04
      a24 = -2.8438341552142710e-05
      a25 = -1.1887778959461776e-05
      a26 = -4.0163964812921489e-06
      b21 =  1.1995545126831476e-09
      b22 =  5.5234008384648383e-12
      b23 =  8.4310335919950873e-13

    case 'CGS':
      a11 =  9.9985372432159340e-01
      a12 =  1.0380621928183473e-02
      a13 =  1.7073577195684715e-03
      a14 = -3.6570490496333680e-05
      a15 = -7.3677944503527477e-06
      a16 = -3.5529175999643348e-06
      b11 =  1.7083494994335439e-10
      b12 =  7.1567921402953455e-13
      b13 =  1.2821026080049485e-13
      a21 =  1.0
      a22 =  1.0316374535350838e-02
      a23 =  8.9521792365142522e-04
      a24 = -2.8438341552142710e-05
      a25 = -1.1887778959461776e-05
      a26 = -4.0163964812921489e-06
      b21 =  1.1995545126831476e-10
      b22 =  5.5234008384648383e-13
      b23 =  8.4310335919950873e-14

    case _:
      print('Unit system:   ' +units+'   not known to eos - exit now')
      exit()

def rho(p,T,S,units='CGS'):
  '''
   pure real(r8) function rho(p, th, s)
   ! ---------------------------------------------------------------------------
   ! In situ density [g cm-3].
   ! ---------------------------------------------------------------------------

      real(r8), intent(in) :: &
         p,    & ! Pressure [g cm-1 s-2].
         th,   & ! Potental temperature [deg C].
         s       ! Salinity [g kg-1].

      rho =  ( a11 + (a12 + a14*th + a15*s)*th + (a13 + a16*s)*s &
             + (b11 + b12*th + b13*s)*p) &
            /( a21 + (a22 + a24*th + a25*s)*th + (a23 + a26*s)*s &
             + (b21 + b22*th + b23*s)*p)

   end function rho
  '''
  init_parameters(units)
  density =  ( a11 + (a12 + a14*T + a15*S)*T + (a13 + a16*S)*S + (b11 + b12*T + b13*S)*p) \
             /( a21 + (a22 + a24*T + a25*S)*T + (a23 + a26*S)*S + (b21 + b22*T + b23*S)*p)

  return density


def p_alpha(p1,p2,T,S):
  '''
  pure real(r8) function p_alpha(p1, p2, th, s)
  ! ---------------------------------------------------------------------------
  ! Integral of specific volume with respect to pressure [cm2 s-2].
  ! ---------------------------------------------------------------------------
   
      real(r8), intent(in) :: &
         p1,   & ! Lower pressure bound [g cm-1 s-2].
         p2,   & ! Upper pressure bound [g cm-1 s-2].
         th,   & ! Potential temperature [deg C].
         s       ! Salinity [g kg-1].

      real(r8), parameter :: &
         r1_3 = 1._r8/3._r8, &
         r1_5 = 1._r8/5._r8, &
         r1_7 = 1._r8/7._r8, &
         r1_9 = 1._r8/9._r8

      real(r8) :: a1, a2, b1, b2, pm, r, q, qq

      a1 = a11 + (a12 + a14*th + a15*s)*th + (a13 + a16*s)*s
      a2 = a21 + (a22 + a24*th + a25*s)*th + (a23 + a26*s)*s
      b1 = b11 + b12*th + b13*s
      b2 = b21 + b22*th + b23*s

      ! The analytic solution of the integral is
      !    p_alpha = ( b2*(p2 - p1)
      !              + (a2 - a1*b2/b1)*log((a1 + b1*p2)/(a1 + b1*p1)))/b1
      ! A truncated series expansion of the integral is used that provides
      ! better computational efficiency and accuarcy for most relevant
      ! parameters.

      pm = .5_r8*(p2 + p1)
      r = .5_r8*(p2 - p1)/(a1 + b1*pm)
      q = b1*r
      qq = q*q

      p_alpha = 2._r8*r*( a2 + b2*pm &
                        + (a2 - a1*b2/b1)*qq*( r1_3 &
                                             + qq*( r1_5 &
                                                  + qq*( r1_7 &
                                                       + qq*  r1_9))))

  end function p_alpha
  '''
  r1_3 = 1./3.
  r1_5 = 1./5.
  r1_7 = 1./7.
  r1_9 = 1./9.

  a1 = a11 + (a12 + a14*T + a15*S)*T + (a13 + a16*S)*S
  a2 = a21 + (a22 + a24*T + a25*S)*T + (a23 + a26*S)*S
  b1 = b11 + b12*T + b13*S
  b2 = b21 + b22*T + b23*S

  pm = 0.5*(p2 + p1)
  r  = 0.5*(p2 - p1)/(a1 + b1*pm)
  q  = b1*r
  qq = q*q
  
  palpha = 2.*r*( a2 + b2*pm \
                    + (a2 - a1*b2/b1)*qq*( r1_3 \
                                             + qq*( r1_5 \
                                                  + qq*( r1_7 \
                                                       + qq*  r1_9))))
  return palpha


