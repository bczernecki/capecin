#' Compute CAPE/CIN according to George H. Bryan formula
#'
#' TODO
#' To be completed..
#'
#' @param pinc = 2.0 # Pressure increment (Pa)
#' @param source Source parcel: 1 = surface, 2 = most unstable (max theta-e), 3 = mixed-layer (specify ml_depth)
#' @param ml_depth  depth (m) of mixed layer (for source = 3); default 200
#' @param adiabat = 4   #! Formulation of moist adiabat: 1 = pseudoadiabatic, liquid only, 2 = reversible, liquid only, 3 = pseudoadiabatic, with ice,  4 = reversible, with ice
#' @export
#'
#' @examples
#' p_in = c(1000, 855, 700, 500, 300, 100, 10)
#' t_in = c(25, 10, 0, -15, -30, -50, -92)
#' td_in = c(20, 5, -5, -30, -55, -80, -99)
#'
#'

source("R/getthe.R")


cape_cin = function(p_in, t_in, td_in,
                    pinc = 100, source = 1, ml_depth = 200, adiabat = 4) {

nk = length(p_in)
g     = 9.81
p00   = 100000.0
cp    = 1005.7
rd    = 287.04
rv    = 461.5
xlv   = 2501000.0
xls   = 2836017.0
t0    = 273.15
cpv   = 1875.0
cpl   = 4190.0
cpi   = 2118.636
lv1   = xlv + (cpl - cpv)*t0
lv2   = cpl - cpv
ls1   = xls + (cpi - cpv)*t0
ls2   = cpi - cpv
rp00  = 1.0/p00
eps   = rd/rv
reps  = rv/rd
rddcp = rd/cp
cpdrd = cp/rd
cpdg  = cp/g

converge = 0.0002
debug_level =   1000


#!---- convert p,t,td to mks units; get pi,q,th,thv ----!
#!print *, "nk = ", nk
#cat(paste("nk = ", nk))
#do k=1,nk
#cat(paste("p_in(k) = ", p_in[i]))
p = 100.0*p_in
t = 273.15 + t_in
td = 273.15 + td_in
pi = (p*rp00)^rddcp
q = getqvs(p,td)
th = t/pi
thv = th*(1.0 + reps*q)/(1.0 + q)

# !print *, "p = ", p
# !print *, "t = ", t

# obliczenia wysokosci:
#!---- get height using the hydrostatic equation ----!
z = 0.0
for (k in 2:nk) {
  dz = -cpdg*0.5*(thv[k] + thv[k - 1])*(pi[k] - pi[k - 1])
  z[k] = z[k - 1] + dz
}


  #!---- find source parcel ----!
  if (source == 1) {
    #! use surface parcel
  kmax = 1

  } else if (source == 2) {
    #! use most unstable parcel (max theta-e)
    if (p[1] < 50000.0) {
      #! first report is above 500 mb ... just use the first level reported
      kmax = 1
      maxthe = getthe(p[1],t[1],td[1],q[1])
    } else {
        # ! find max thetae below 500 mb
      maxthe = 0.0
      for (k in 1:nk) {
        if (p[k] >= 50000.0) {
          the = getthe(p[k], t[k], td[k], q[k])
          if (the > maxthe) {
            maxthe = the
            kmax = k
          }
        } # koniec ifa
      } # koniec for'a
    } # koniec else

  } else if ( source == 3 ) {
    #! use mixed layer
    if ((z[2] - z[1]) > ml_depth ) {
      #! the second level is above the mixed-layer depth:  just use the
      #! lowest level
      avgth = th[1]
      avgqv = q[1]
      kmax = 1
    } else if (z[nk] < ml_depth) {
      #! the top-most level is within the mixed layer:  just use the
      #! upper-most level
      avgth = th[nk]
      avgqv = q[nk]
      kmax = nk
    } else {
      #! calculate the mixed-layer properties:
      avgth = 0.0
      avgqv = 0.0
      k = 2
      if (debug_level > 100) print(paste('  ml_depth = ', ml_depth))
      if (debug_level > 100) print(paste('  k,z,th,q:'))
      if (debug_level > 100) print(paste(z[1], th[1], q[1]))

      while ((z[k] <= ml_depth) & (k <= nk)) {
        if (debug_level >= 100) print(paste(k, z[k], th[k], q[k]))
        avgth = avgth + 0.5*(z[k] - z[k - 1])*(th[k] + th[k - 1])
        avgqv = avgqv + 0.5*(z[k] - z[k - 1])*(q[k] + q[k - 1])
        k = k + 1
      }

      th2 = th[k - 1] + (th[k] - th[k - 1])*(ml_depth - z[k - 1])/(z[k] - z[k - 1])
      qv2 =  q[k - 1] + (q[k] - q[k - 1])*(ml_depth - z[k - 1])/(z[k] - z[k - 1])

      if (debug_level > 100) print(paste("L127\t", ml_depth,th2,qv2))

      avgth = avgth + 0.5*(ml_depth - z[k - 1])*(th2 + th[k - 1])
      avgqv = avgqv + 0.5*(ml_depth - z[k - 1])*(qv2 + q[k - 1])

      if (debug_level > 100) print(paste("L132", k, z[k], th[k], q[k]))

      avgth = avgth/ml_depth
      avgqv = avgqv/ml_depth

      kmax = 1
    }

    } else {
        cat(paste("source not defined properly"))
      stop()
    } # koniec dlugiego ifa do liczenia source (parcel: 1,2,3,4)



#!---- define parcel properties at initial location ----!
narea = 0.0
if ( any(source %in% 1:2 )) {
  k    = kmax
  th2  = th[kmax]
  pi2  = pi[kmax]
  p2   = p[kmax]
  t2   = t[kmax]
  thv2 = thv[kmax]
  qv2  = q[kmax]
  b2   = 0.0
} else if ( source == 3 ) {
  k    = kmax
  th2  = avgth
  qv2  = avgqv
  thv2 = th2*(1.0 + reps*qv2)/(1.0 + qv2)
  pi2  = pi[kmax]
  p2   = p[kmax]
  t2   = th2*pi2
  b2   = g*(thv2 - thv[kmax])/thv[kmax]
}


ql2 = 0.0
qi2 = 0.0
qt  = qv2

cape = 0.0
cin  = 0.0
lfc  = 0.0

doit = TRUE
cloud = FALSE
  if (any(adiabat %in% 1:2)) {
    ice = FALSE
  } else {
    ice = TRUE
  }

the = getthe(p2,t2,t2,qv2)
# if(debug_level.ge.100) print *,'  the = ',the













#!---- begin ascent of parcel ----!
if (debug_level > 100) {
  print('  Start loop:')
  print(paste('  p2,th2,qv2 = ',p2,th2,qv2))
}

while (doit & (k < nk)) {

  k = k + 1
  b1 = b2

  dp = p[k - 1] - p[k]

  if ( dp < pinc ) {
    nloop = 1
  } else {
    # TODO:
    nloop = 1 + trunc( dp/pinc )
    dp = dp/nloop
  }

  for (n in 1:nloop) {

    p1 = p2
    t1 = t2
    pi1 = pi2
    th1 = th2
    qv1 = qv2
    ql1 = ql2
    qi1 = qi2
    thv1 = thv2

    p2 = p2 - dp
    pi2 = (p2*rp00)^rddcp

    thlast = th1
    i = 0
    not_converged = TRUE

    while (not_converged) {
      i = i + 1
      t2 = thlast*pi2
      if (ice) {
        # TODO:
        fliq = max(c(min(c((t2 - 233.15)/(273.15 - 233.15), 1.0)), 0.0))
        fice = 1.0 - fliq
      } else {
        fliq = 1.0
        fice = 0.0
      }
      qv2 = min( c(qt, fliq*getqvs(p2,t2) + fice*getqvi(p2, t2) ))
      qi2 = max( c(fice*(qt - qv2), 0.0))
      ql2 = max( c(qt - qv2 - qi2, 0.0))

      tbar  = 0.5*(t1 + t2)
      qvbar = 0.5*(qv1 + qv2)
      qlbar = 0.5*(ql1 + ql2)
      qibar = 0.5*(qi1 + qi2)

      lhv = lv1 - lv2*tbar
      lhs = ls1 - ls2*tbar
      lhf = lhs - lhv

      rm = rd + rv*qvbar
      cpm = cp + cpv*qvbar + cpl*qlbar + cpi*qibar
      th2 = th1*exp(lhv*(ql2 - ql1)/(cpm*tbar) + lhs*(qi2 - qi1)/(cpm*tbar) + (rm/cpm - rd/cp)*log(p2/p1))

      if (i > 90) print(paste(i, th2, thlast, th2 - thlast))
      if (i > 100) {
        print('  Error:  lack of convergence')
        print('  ... stopping iteration ')
        stop()
      }

      if ( abs(th2 - thlast) > converge ) {
        thlast = thlast + 0.3*(th2 - thlast)
      } else {
        not_converged = FALSE
      }

    } # koniec petli while (?)

#! Latest pressure increment is complete.  Calculate some important stuff:

  if ( ql2 >= 1.0e-10 ) cloud = TRUE

  if (any(adiabat %in% c(1, 3))) {
#! pseudoadiabat
    qt  = qv2
    ql2 = 0.0
    qi2 = 0.0
  } else if (adiabat <= 0 | adiabat >= 5) {
    print('  Undefined adiabat')
    stop()
  }

  } # koniec petli po nloop?

  thv2 = th2*(1.0 + reps*qv2)/(1.0 + qv2 + ql2 + qi2)
  b2 = g*(thv2 - thv[k])/thv[k]
  dz = -cpdg*0.5*(thv[k] + thv[k - 1])*(pi[k] - pi[k - 1])
  the = getthe(p2, t2, t2, qv2)

#! Get contributions to CAPE and CIN:

  if ( (b2 >= 0.0) & (b1 < 0.0) ) {
#! first trip into positive area
    ps = p[k - 1] + (p[k] - p[k - 1])*(0.0 - b1)/(b2 - b1)
    frac = b2/(b2 - b1)
    parea =  0.5*b2*dz*frac
    narea = narea - 0.5*b1*dz*(1.0 - frac)

    if (debug_level > 200) {
      print(paste('      b1,b2 = ', b1,b2))
      print(paste('      p1,ps,p2 = ', p[k - 1], ps, p[k]))
      print(paste('      frac = ', frac))
      print(paste('      parea = ', parea))
      print(paste('      narea = ', narea))
    }

    cin  = cin  + narea
    narea = 0.0

    } else if ( (b2 < 0.0) & (b1 > 0.0) ) {
#! first trip into neg area
      ps = p[k - 1] + (p[k] - p[k - 1])*(0.0 - b1)/(b2 - b1)
      frac = b1/(b1 - b2)
      parea =  0.5*b1*dz*frac
      narea = -0.5*b2*dz*(1.0 - frac)

      if (debug_level >= 200) {
        print(paste('      b1,b2 = ', b1, b2))
        print(paste('      p1,ps,p2 = ', p[k - 1], ps, p[k]))
        print(paste('      frac = ', frac))
        print(paste('      parea = ', parea))
        print(paste('      narea = ', narea))
      }
    } else if ( b2 < 0.0 ) {
#! still collecting negative buoyancy
      parea =  0.0
      narea = narea-0.5*dz*(b1 + b2)
    } else {
# ! still collecting positive buoyancy
      parea =  0.5*dz*(b1 + b2)
      narea =  0.0
    }

print(paste( "cape & parea = ", cape, parea))
cape = cape + max(c(0.0, parea))

if (debug_level >= 200) {
  print(paste(p2,b1,b2,cape,cin,cloud))
}

if( (p[k] <= 10000.0) & (b2 < 0.0) ) {
  #! stop if b < 0 and p < 100 mb
  doit = FALSE
}

} # koniec glownej petli


# enddo
#
# !---- All done ----!
#
#   return
# end subroutine getcape
return(list(cape, cin))
} # tu dopisac koniec funkcji
