function genplt (args)

                           n = 1
expid   = subwrd(args,n) ; n = n + 1
EXPORT  = subwrd(args,n) ; n = n + 1
GC      = subwrd(args,n) ; n = n + 1
alias   = subwrd(args,n) ; n = n + 1
season  = subwrd(args,n) ; n = n + 1
output  = subwrd(args,n) ; n = n + 1
level   = subwrd(args,n) ; n = n + 1
nmod    = subwrd(args,n) ; n = n + 1
nobs    = subwrd(args,n) ; n = n + 1
expfile = subwrd(args,n) ; n = n + 1
anafile = subwrd(args,n) ; n = n + 1
anal    = subwrd(args,n) ; n = n + 1
obsname = subwrd(args,n) ; n = n + 1
debug   = subwrd(args,n) ; n = n + 1
expdsc  = subwrd(args,n) ; n = n + 1
stat    = subwrd(args,n) ; n = n + 1
 
blak   = 0

* Get Dates
* ---------
'run getenv "BEGDATEO"'
         begdateo = result
'run getenv "ENDDATEO"'
         enddateo = result

'run getenv "BEGDATE"'
         begdate  = result
'run getenv "ENDDATE"'
         enddate  = result
if( begdate = "NULL" )
   'set dfile 'expfile
   'set t    '1
   'getinfo date'
         begdate = result
endif
if( enddate = "NULL" )
   'set dfile 'expfile
   'getinfo tdim'
            tdim = result
   'set t  'tdim
   'getinfo date'
         enddate = result
endif

'run getenv "CLIMATE"'
         climate = result
if( begdate = begdateo & enddate = enddateo ) 
         climate = 'Actual'
endif

* Check for Contour Level Type
* ----------------------------
'run getenv "LEVTYPE"'
           CLEVS  = result

* Get Plotting Values from Resource File
* --------------------------------------
'run getenv "GEOSUTIL"'
         geosutil = result
PLOTRC = geosutil'/plots/grads_util/plot.rc'

    PRFX = ''
if( stat = 'STD' )
    PRFX = 'STD_'
endif
if( stat = 'RMS' )
    PRFX = 'RMS_'
endif
if( stat = 'BIAS' )
    PRFX = 'BIAS_'
endif

say 'LEVTYPE: 'CLEVS
                        'getresource 'PLOTRC' 'PRFX''EXPORT'_'GC'_'level'_CBSCALE'
if( result = 'NULL' ) ; 'getresource 'PLOTRC' 'PRFX''EXPORT'_'GC'_CBSCALE' ; endif
                                                                  cbscale = result

                        'getresource 'PLOTRC' 'PRFX''EXPORT'_'GC'_'level'_FACTOR'
if( result = 'NULL' ) ; 'getresource 'PLOTRC' 'PRFX''EXPORT'_'GC'_FACTOR' ; endif
                                                                  fact = result

                        'getresource 'PLOTRC' 'PRFX''EXPORT'_'GC'_FIXED_PLOT_FACTOR'
                                                                  fixpltfact = result
                        'getresource 'PLOTRC' 'PRFX''EXPORT'_'GC'_FIXED_PLOT_CINT'
                                                                  fixpltcint = result

                        'getresource 'PLOTRC' 'PRFX''EXPORT'_'GC'_'level'_TITLE'
if( result = 'NULL' ) ; 'getresource 'PLOTRC' 'PRFX''EXPORT'_'GC'_TITLE' ; endif
                                                                  title = result

                        'getresource 'PLOTRC' 'PRFX''EXPORT'_'GC'_'level'_CCOLS'
if( result = 'NULL' ) ; 'getresource 'PLOTRC' 'PRFX''EXPORT'_'GC'_CCOLS' ; endif
                                                                  ccols = result

                        'getresource 'PLOTRC' 'PRFX''EXPORT'_'GC'_'level'_'CLEVS
if( result = 'NULL' ) ; 'getresource 'PLOTRC' 'PRFX''EXPORT'_'GC'_'level'_CLEVS' ; endif
if( result = 'NULL' ) ; 'getresource 'PLOTRC' 'PRFX''EXPORT'_'GC'_'CLEVS ; endif
if( result = 'NULL' ) ; 'getresource 'PLOTRC' 'PRFX''EXPORT'_'GC'_CLEVS' ; endif
                                                                  clevs = result

                        'getresource 'PLOTRC' 'mname'_'GC'_'level'_DPCT'
if( result = 'NULL' ) ; 'getresource 'PLOTRC' 'mname'_'GC'_DPCT' ; endif
                                                           dpct = result

                        'getresource 'PLOTRC' 'PRFX''EXPORT'_'GC'_REGRID'
                                                                  method = result
                        'getresource 'PLOTRC' 'PRFX''EXPORT'_'GC'_MASK'
                                                                  mask   = result

* Remove possible BLANKS from mask
* --------------------------------
DUMMY = ''
length = strlen(result)
i = 1
while( i<=length )
  bit = substr(result,i,1)
  if( bit != ' ' )
      if( DUMMY = '' )
          DUMMY = bit
      else
          DUMMY = DUMMY''bit
      endif
  endif
i = i+1
endwhile
mask = DUMMY


if( title   = 'NULL' )
    title   =  EXPORT':'GC
else
    title   =  EXPORT':'GC'  'title
endif
if( cbscale = 'NULL' ) ; cbscale =  0.8         ; endif
if( clab    = 'NULL' ) ; clab    =  on          ; endif
if( dpct    = 'NULL' ) ; dpct    =  0.1         ; endif

                     'd 'dpct
                         dpct    = subwrd(result,4)
say ''
say 'title: 'title
say ' fact: ' fact
say 'ccols: 'ccols
say 'clevs: 'clevs
say ' mask: ' mask
say ''

* Remove possible BLANKS from FACTOR
* ----------------------------------
DESC = ''
length = getlength(fact)
i = 1
while( i<=length )
  bit = substr(fact,i,1)
  if( bit != ' ' )
      if( DESC = '' )
          DESC = bit
      else
          DESC = DESC''bit
      endif
  endif
i = i+1
endwhile
fact = DESC

* Set Default Plotting Variables
********************************

if( fact    = 'NULL' ) ; fact    = 1         ; endif

* Plot Mean Field
* ---------------
'c'
'set display color white'
'set vpage off'
'set parea off'
'set grid  off'
'set mproj scaled'
'set frame on'
'set xlopts 1 3 .11'
'set ylopts 1 3 .11'
'rgbset'
'set rgb 84 204 204 204'
'set rgb 85 137 137 137'

'set dfile 'expfile
'setlons'
'setlats'
if( level = 0 )
   'set z 1'
else
   'set lev 'level
endif
'set t 1'
'q dims'
say 'EXP DIMS Environment: 'result

* Get Dimension of Environment
* ----------------------------
'getinfo lonmin'
         lonbeg = result
'getinfo lonmax'
         lonend = result
'getinfo latmin'
         latbeg = result
'getinfo latmax'
         latend = result

say 'Environment Dimension: 'lonbeg' 'lonend' 'latbeg' 'latend

'define qmod  = mod'season'*'fact
'define qmod  = regrid2( qmod,0.25,0.25,bs_p1,'lonbeg','latbeg')'

if( mask != NULL )
   'setmask'
   'define lwmask = regrid2( lwmask,0.25,0.25,bs_p1,'lonbeg','latbeg')'
   if( mask = 'LAND' )
       say 'define qmod = maskout( qmod, 0.5-lwmask )'
           'define qmod = maskout( qmod, 0.5-lwmask )'
   endif
   if( mask = 'OCEAN' )
       say 'define qmod = maskout( qmod, lwmask-0.5 )'
           'define qmod = maskout( qmod, lwmask-0.5 )'
   endif
endif

'define maskm = 1 + qmod-qmod'



* Determine DLAT & DLON of Analysis
* ---------------------------------
'set dfile 'anafile
'set z 1'
'set t 1'
'getinfo dlat'
         dlat = result
'getinfo dlon'
         dlon = result
'set gxout shaded'

say 'Analysis DLAT: 'dlat
say 'Analysis DLON: 'dlon

'set lon 'lonbeg' 'lonend
'set lat 'latbeg' 'latend
'define qobs  = obs'season'*'fact
'define qobs  = regrid2( qobs,0.25,0.25,bs_p1,'lonbeg','latbeg')'
'define qobs  = maskout(qobs,maskm)'

   'set gxout stat'
   'd qmod'
   qminmax = sublin(result,8)
   qmin    = subwrd(qminmax,4)
   qmax    = subwrd(qminmax,5)
   say 'QMin Value: 'qmin
   say 'QMax Value: 'qmax
   'set gxout shaded'
   'd abs('qmin')'
           qmin = subwrd(result,4)
   'd abs('qmax')'
           qmax = subwrd(result,4)
   if( qmin > qmax ) ; qmax = qmin ; endif

   say 'Absolute QMAX: 'qmax

m = 0
if( ccols = NULL )
   if( qmax > 0 )
      'd log10('qmax')'
       m = subwrd(result,4)
   else
       m = 0
   endif
   say '    Log Factor: 'm
   if( m<0 ) ; m = m-2 ; endif
   'getint 'm
            m = result
   if( m>0 )
       if( m<=2 )
           m = 0
       else
           m = m-2
       endif
   endif
   say 'Field Scaling Factor: 'm
   if( m<0 )
       j = -1 * m
      'define qmod = qmod * 1e'j
      'define qobs = qobs * 1e'j
   else
      'define qmod = qmod / 1e'm
      'define qobs = qobs / 1e'm
   endif
endif

'set dfile 'expfile
'set lon 'lonbeg' 'lonend
'set lat 'latbeg' 'latend


* Make Plot: Top Panel
* --------------------
'set vpage 0 8.5 0.0 11'
'set parea 1.5 7.0 7.70 10.50'
'set grads off'

   'set gxout stat'
   'd qmod'
   qmodminmax = sublin(result,8)
   qmodmin    = subwrd(qmodminmax,4)
   qmodmax    = subwrd(qmodminmax,5)
   say 'QMOD_Max Value: 'qmodmax
   say 'QMOD_Min Value: 'qmodmin
   'set gxout shaded'
   'd abs('qmodmin')'
          aqmodmin = subwrd(result,4)
   'd abs('qmodmax')'
          aqmodmax = subwrd(result,4)
   if( aqmodmin > aqmodmax ) ; aqmodmax = aqmodmin ; endif
   say 'Absolute QMOD_MAX: ' aqmodmax

   'set gxout shaded'
if( ccols != NULL )
   'set clevs 'clevs
   'set ccols 'ccols
   'set clab  'clab
else
   'shades 'qmod' 0'
    cint = result*2
    if( clevs != NULL ) ; 'set clevs 'clevs ; endif
endif
   'd qmod'

* Make Plot: Middle Panel
* -----------------------
'set parea off'
'set vpage 0 8.5 0.0 11'
'set parea 1.5 7.0 4.30 7.10'
'set grads off'

   'set gxout stat'
   'd qobs'
   qobsminmax = sublin(result,8)
   qobsmin    = subwrd(qobsminmax,4)
   qobsmax    = subwrd(qobsminmax,5)
   say 'QOBS_Max Value: 'qobsmax
   say 'QOBS_Min Value: 'qobsmin
   'set gxout shaded'
   'd abs('qobsmin')'
          aqobsmin = subwrd(result,4)
   'd abs('qobsmax')'
          aqobsmax = subwrd(result,4)
   if( aqobsmin > aqobsmax ) ; aqobsmax = aqobsmin ; endif
   say 'Absolute QOBS_MAX: ' aqobsmax

if( ccols != NULL )
   'set gxout shaded'
   'set clevs 'clevs
   'set ccols 'ccols
else
   'shades 'qmod' 0'
    if( clevs != NULL ) ; 'set clevs 'clevs ; endif
endif
   'd qobs'

   'cbarn -vert -snum 0.8 -ymid 6.4 -scaley 0.9 '

* Make Plot: Bottom Panel
* -----------------------
'set parea off'
'set vpage 0 8.5 0.0 11'
'set parea 1.5 7.0 0.90 3.70'
'set grads off'
'getinfo lon'
         lon = result
'set dfile 'expfile
if( level = 0 )
   'set z 1'
else
   'set lev 'level
endif
'set t 1'
'q dims'

'stats maskout(qmod,abs(qobs))'
 avgmod = subwrd(result,1)
 stdmod = subwrd(result,2)

'stats maskout(qobs,abs(qobs))'
 avgobs = subwrd(result,1)
 stdobs = subwrd(result,2)

'stats maskout(qmod-qobs,abs(qobs))'
 avgdif = subwrd(result,1)
 stddif = subwrd(result,2)

  dqmax = stddif/3

'set gxout shaded'

     say 'Absolute DQMAX: 'dqmax'  qmax: 'qmax
     dqrel = dqmax / qmax  * 100 * 100
    'getint 'dqrel
             dqrel = result/100

     say 'QMAX: 'qmax'  DQMAX: 'dqmax
     say 'Relative PCT Difference: 'dqrel' (100*DQMAX/QMAX)'
     say ' Minimum PCT for  Plots: 'dpct

     if( dqrel < dpct )
         dqrel = dpct
     endif
         dqmax = dqrel * qmax / 100
         say 'Setting CINT using DQREL: 'dqrel'%, DQMAX: 'dqmax

   if( dqmax > 0 )
      'd log10('dqmax')'
       n = subwrd(result,4)
   else
       n = 0
   endif
   say '    Log Factor: 'n
   if( n<0 ) ; n = n-2 ; endif
   'getint 'n
            n = result
   if( n>0 )
       if( n<=2 )
           n = 0
        else
           n = n+2
        endif
   endif

   if( fixpltfact != NULL )
       'd 'fixpltfact
        n =subwrd(result,4)
   endif

   say 'Diff Scaling Factor: 'n

   if( fixpltcint != NULL )
      'd  'fixpltcint
           fixpltcint =subwrd(result,4)
       cint = fixpltcint
   else
      'd 'dqmax'/1e'n
       cint = subwrd(result,4)
   endif

      'shades 'cint
      'define qdif = (qmod-qobs)/1e'n
      'd qdif'
      'cbarn -snum 0.55 -xmid 4.25 -ymid 0.4'

      'set gxout stat'
      'd qdif'
      qdifminmax = sublin(result,8)
      qdifmin    = subwrd(qdifminmax,4)
      qdifmax    = subwrd(qdifminmax,5)
      say 'QDIF_Max Value: 'qdifmax
      say 'QDIF_Min Value: 'qdifmin
      'set gxout shaded'
      'd abs('qdifmin')'
             aqdifmin = subwrd(result,4)
      'd abs('qdifmax')'
             aqdifmax = subwrd(result,4)
      if( aqdifmin > aqdifmax ) ; aqdifmax = aqdifmin ; endif
      say 'Absolute QDIF_MAX: ' aqdifmax

      'set gxout stat'
      'd maskout(qdif,abs(qobs))'
         dqminmax = sublin(result,8)
         dqmin    = subwrd(qminmax,4)
         dqmax    = subwrd(qminmax,5)
      say 'DQMin Value: 'dqmin
      say 'DQMax Value: 'dqmax
      'set gxout shaded'

'stats maskout(qdif,abs(qobs))'
 avgdif = subwrd(result,1)
 stddif = subwrd(result,2)

k = m + n

'set vpage off'
'set string 1 l 4'
'set strsiz 0.065'
'draw string 0.05 0.08 ( EXPID:  'expid' )'

'set string 1 c 6'
'set strsiz .125'
'draw string 4.25 10.85 'title

'set strsiz .11'
if( level = 0 )
    if( m != 0 )
   'draw string 4.25 10.62 'expdsc'  'season' ('nmod') (x 10**'m')'
    else
   'draw string 4.25 10.62 'expdsc'  'season' ('nmod')'
    endif
   'draw string 4.25 7.22 'obsname'  'season' ('nobs')  ('climate')'
   if( k != 0 )
   'draw string 4.25 3.80 Difference (Top-Middle) (x 10**'k')'
   else
   'draw string 4.25 3.80 Difference (Top-Middle)'
   endif
else
    if( m != 0 )
   'draw string 4.25 10.62 'expdsc'  'level'-mb  'season' ('nmod') (x 10**'m')'
    else
   'draw string 4.25 10.62 'expdsc'  'level'-mb  'season' ('nmod')'
    endif
   'draw string 4.25 7.22 'obsname'  'season' ('nobs')  ('climate')'
   if( k != 0 )
   'draw string 4.25 3.80 'level'-mb  Difference (Top-Middle) (x 10**'k')'
   else
   'draw string 4.25 3.80 'level'-mb  Difference (Top-Middle)'
   endif
endif

                date = getdate (begdate)
bmnthm = subwrd(date,1)
byearm = subwrd(date,2)
                date = getdate (enddate)
emnthm = subwrd(date,1)
eyearm = subwrd(date,2)
                date = getdate (begdateo)
bmntho = subwrd(date,1)
byearo = subwrd(date,2)
                date = getdate (enddateo)
emntho = subwrd(date,1)
eyearo = subwrd(date,2)

'set string 1 l 4'
'set strsiz .08'

'draw string 0.050 10.25 Beg: 'bmnthm' 'byearm
'draw string 0.050 10.10 End: 'emnthm' 'eyearm
'draw string 0.050 9.85  Max: 'qmodmax
'draw string 0.050 9.70  Min: 'qmodmin
'draw string 0.050 9.40 Mean: 'avgmod
'draw string 0.050 9.25  Std: 'stdmod

'draw string 0.050 6.85 Beg: 'bmntho' 'byearo
'draw string 0.050 6.70 End: 'emntho' 'eyearo
'draw string 0.050 6.45  Max: 'qobsmax
'draw string 0.050 6.30  Min: 'qobsmin
'draw string 0.050 6.00 Mean: 'avgobs
'draw string 0.050 5.85  Std: 'stdobs

'draw string 0.050 3.45 Beg: 'bmntho' 'byearo
'draw string 0.050 3.30 End: 'emntho' 'eyearo
'draw string 0.050 3.05  Max: 'qdifmax
'draw string 0.050 2.90  Min: 'qdifmin
'draw string 0.050 2.60 Mean: 'avgdif
'draw string 0.050 2.45  Std: 'stddif

if( CINTDIFF != 'NULL' )
   'set strsiz .07'
   'draw string 0.050 1.77 Plot represents'
   'draw string 0.050 1.62 values > 'dqrel' %'
   'draw string 0.050 1.47 Relative Difference'
   'draw string 0.050 1.32 ( DQ/QMax )'
endif

'myprint -name 'output'/hdiag_'PRFX''anal'_'EXPORT'.'GC'_'level'.'season
'set clab on'

'set mproj latlon'
return

function getdate (date,month,year)
       num = 1
       bit = substr(date,num,1)
while( bit != '' )
       num = num+1
       bit = substr(date,num,1)
endwhile
       loc = num-7
     month = substr(date,loc  ,3)
      year = substr(date,loc+3,4)
return month' 'year

function getlength (string)
tb = ""
i = 1
while (i<=80)
blank = substr(string,i,1)
if( blank = tb )
length = i-1
i = 81
else
i = i + 1
endif
endwhile
return length

