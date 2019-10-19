from rsf.proj import *
import math

suspike  = os.path.join(os.getenv('CWPROOT'),'bin/suspike')
sugain   = os.path.join(os.getenv('CWPROOT'),'bin/sugain')
sufilter = os.path.join(os.getenv('CWPROOT'),'bin/sufilter')


#--------------------------#
#Basic Utils - Raanan Dafni#
#--------------------------#

#*************************************************************************************************#
#******************************************** VELOCITY *******************************************#
#*************************************************************************************************#

#-------------------------------------------------------------------------------
# Velocity Model Building
def MakeVel (fileOut='velModel', **params) :

	Flow (fileOut, None,
		'''
		math n1=%(nx)d o1=0 d1=%(dx)f output="%(Z0)f+%(dzdx)f*x1" |
		unif2 n1=%(nz)d d1=%(dz)f n2=%(nx)d d2=%(dx)f v00=%(V1)f,%(V2)f |
		put n3=%(ny)d d3=%(dy)f o3=%(origy)f
		''' % params )

#-------------------------------------------------------------------------------
# Velocity Model Building (Multi-Layer)
def MakeVel_MultiLayer (fileOut='velModel', **params) :

        Flow (fileOut, None,
                '''
                unif2 n1=%(nz)d d1=%(dz)f n2=%(nx)d d2=%(dx)f v00=%(V1)f,%(V2)f,%(V3)f,%(V4)f z0=%(Z0)f,%(Z1)f,%(Z2)f x0=0,0,0 |
                put n3=%(ny)d d3=%(dy)f o3=%(origy)f
                ''' % params )

#-------------------------------------------------------------------------------
# Square Velocity
def SquareVel (fileOut='csqModel', fileIn='velModel') :

	Flow (fileOut, fileIn,
		'''
		mul ${SOURCES[0]} | put data_type=csq
		''')

#-------------------------------------------------------------------------------
# Smooth Velocity
def SmoothVel (fileOut='csqModelSm', fileIn='csqModel', **params) :

	Flow (fileOut, fileIn,
		'''
		smooth rect1=%(rect1)d rect2=%(rect2)d repeat=%(repeat)d 
		''' % params )

#-------------------------------------------------------------------------------
# Pertub Velocity
def PerturbVel (fileOut='csqPerturb', fileIn1='csqModel', fileIn2='csqModelSm') :

	Flow (fileOut, [fileIn1, fileIn2], 
		'''
		add scale=1,-1 ${SOURCES[1]}
		''')

#-------------------------------------------------------------------------------
# Inverse Velocity
def InvVel (fileOut='SlModel', fileIn='velModel') :

	Flow(fileOut, fileIn,
	     	'''
	     	math output=1/input
	     	''')

#-------------------------------------------------------------------------------
# Transpose Velocity
def TranspVel (fileOut='transpModel', fileIn='Model') :

	Flow(fileOut, fileIn,
	     	'''
	     	transp plane12 |
		put id1=1 id2=0
	     	''')

#-------------------------------------------------------------------------------
# Transpose Extended Velocity
def TranspExtVel (fileOut='transpExtModel', fileIn='Model') :

	Flow(fileOut, fileIn,
	     	'''
	     	transp plane=12 |
		put id1=1 id2=0 id3=101 dim=2 gdim=3
	     	''')

#-------------------------------------------------------------------------------
# Extend Velocity
def ExtVel (fileOut='ExtModel', fileIn='Model', **params) :

	Flow (fileOut, fileIn,
		'''
		pad beg3=%(nh)d end3=%(nh)d |
		put d3=%(dh)f o3=%(oh)f dim=2 gdim=3 id1=0 id2=1 id3=101
		''' % params )


#-------------------------------------------------------------------------------
# Bulk Modulus model Building
def MakeBulkMod (fileOut='bulkmodModel', fileIn1='csqModel', fileIn2='denModel') :

	Flow (fileOut, [fileIn1, fileIn2], 
		'''
		mul ${SOURCES[1]}
		''')

#*************************************************************************************************#
#******************************************** WAVELET ********************************************#
#*************************************************************************************************#

#Create band-pass Base Wavelet
def CreateWaveletBP (fileOut='wavelet_base.su', **params) :

        #Band Pass wavelet
        Flow (fileOut, None, suspike + ' nt=%(nt)d dt=%(dt)f it1=%(it)d ntr=1 offset=0 ix1=1 nspk=1 | ' % params 
	                   + sugain + ' scale=%(scale)f | ' % params + sufilter + ' f=%(f1)f,%(f2)f,%(f3)f,%(f4)f ' % params )

#-------------------------------------------------------------------------------
# Create Base Wavelet
def CreateWavelet (fileOut1='wavelet', fileOut2='wavelethead', fileOut3='wavelet_base.su', **params) :

	#Gaussian derivative wavelet
	Flow (fileOut1, None,
		'''
		math n1=%(nt)d d1=%(dt)f output="%(f_s1)g*(x1-%(t0)g)*exp(-%(f_s2)g*(x1-%(t0)g)^2)" |
		put o1=%(orig)f
		''' % params )

	#create header file
	Flow (fileOut2, fileOut1,
		'''
		segyheader
		''')

	#convert to SU format
	Flow (fileOut3, [fileOut1, fileOut2],
		'''
		suwrite endian=0 tfile=${SOURCES[1]}
		''')

#-------------------------------------------------------------------------------
# Create Base Wavelet (asg to acd)
def CreateWavelet_asg2acd (fileOut1='wavelet', fileOut2='wavelethead', fileOut3='wavelet_base.su', **params) :

        #second Gaussian derivative wavelet (Ricker)
        Flow (fileOut1, None,
                '''
                math n1=%(nt)d d1=%(dt)f output="%(f_s1)g*(1-2*%(f_s2)g*(x1-%(t0)g)^2)*exp(-%(f_s2)g*(x1-%(t0)g)^2)" |
                put o1=%(orig)f |
		add scale=%(den)g
                ''' % params )

        #create header file
        Flow (fileOut2, fileOut1,
                '''
                segyheader
                ''')

        #convert to SU format
        Flow (fileOut3, [fileOut1, fileOut2],
                '''
                suwrite endian=0 tfile=${SOURCES[1]}
                ''')

#-------------------------------------------------------------------------------
# Design Velocity Geometry
def DesignVelGeom (fileOut1='geomVelGrid', fileOut2='VelMap', **params) :

	Flow (fileOut1, None, 
		'''
		spike n1=%(nz)d d1=%(dz)f n2=%(nx)d o2=0 d2=1 mag=0
		''' % params)

	Flow (fileOut2, fileOut1,'window n1=1')

#-------------------------------------------------------------------------------
# Create Velocity Headers
def CreateVelHeaders (fileOut1='Velhead', fileOut2='Velhead.su', fileIn1='geomVelGrid', fileIn2='VelMap', **params) :

	Flow ('traclvel', fileIn2, 'math output=x1                           | dd type=int' % params )
	Flow ('sxvel',    fileIn2, 'math output="%(ox)f+%(dx)f*(x1-1)"       | dd type=int' % params )
	Flow ('gxvel',    fileIn2, 'math output="%(ox)f+%(dx)f*(x1-1)"       | dd type=int' % params )
	Flow ('ofstvel',  'sxvel gxvel', 'add scale=-1,1 ${SOURCES[1]}')

	#create header file
	Flow (fileOut1, [fileIn1, 'traclvel', 'sxvel', 'gxvel', 'ofstvel'],
		'''
		segyheader tracl=${SOURCES[1]} sx=${SOURCES[2]} gx=${SOURCES[3]} offset=${SOURCES[4]}
		''')

	#convert to SU format
	Flow (fileOut2, [fileIn1, fileOut1],
		'''
		suwrite tfile=${SOURCES[1]} endian=0
		''')


#*************************************************************************************************#
#************************************** ACQUSITION GEOMETRY **************************************#
#*************************************************************************************************#

#-------------------------------------------------------------------------------
# Design Shot-Geophone Geometry
def DesignShotGeoGeom (fileOut1='geomGrid', fileOut2='SGMap', **params) :

	Flow (fileOut1, None, 
		'''
		spike n1=%(nt)d d1=%(dt)f n2=%(ng)d o2=%(origg)d d2=%(dg)d
		      n3=%(ns)d o3=%(origs)d d3=%(ds)d mag=0
		''' % params)

	Flow (fileOut2, fileOut1,'window n1=1')

#-------------------------------------------------------------------------------
# Create Shot-Gathers Headers
def CreateSGHeaders (fileOut1='SGhead', fileOut2='SGhead.su', fileIn1='geomGrid', fileIn2='SGMap', **params) :

	#create the header fields
	if params['ns'] == 1:
		Flow ('tracl', fileIn2, 'math output=x1                           | dd type=int' % params )
		Flow ('selev', fileIn2, 'math output="-1*%(dz)f"                  | dd type=int' % params )
		Flow ('gelev', fileIn2, 'math output="-1*%(dz)f"                  | dd type=int' % params )
		Flow ('sx',    fileIn2, 'math output=%(sht0)f                     | dd type=int' % params )
		Flow ('gx',    fileIn2, 'math output="%(geo0)f+%(geoint)f*(x1-1)" | dd type=int' % params )
		Flow ('ofst',  'sx gx', 'add scale=-1,1 ${SOURCES[1]}')
                Flow ('cdp',   fileIn2, 'math output="%(fcdp)d+x1"                | dd type=int' % params )
	else:
		Flow ('tracl', fileIn2, 'math output="x1+%(ng)d*(x2-1)"                             | dd type=int' % params )
		Flow ('selev', fileIn2, 'math output="-1*%(dz)f"                                    | dd type=int' % params )
		Flow ('gelev', fileIn2, 'math output="-1*%(dz)f"                                    | dd type=int' % params )
		Flow ('sx',    fileIn2, 'math output="%(sht0)f+%(shtint)f*(x2-1)"                   | dd type=int' % params )
		Flow ('gx',    fileIn2, 'math output="%(geo0)f+%(geoint)f*(x1-1)+%(shtint)f*(x2-1)" | dd type=int' % params )
		Flow ('ofst',  'sx gx', 'add scale=-1,1 ${SOURCES[1]}')
		Flow ('cdp',   fileIn2, 'math output="%(fcdp)d+x1+%(acdp)d*(x2-1)"                  | dd type=int' % params )

	#create header file
	Flow (fileOut1, [fileIn1, 'tracl', 'selev', 'gelev', 'sx', 'gx', 'ofst', 'cdp'],
		'''
		segyheader tracl=${SOURCES[1]} selev=${SOURCES[2]}  gelev=${SOURCES[3]}
			      sx=${SOURCES[4]}    gx=${SOURCES[5]} offset=${SOURCES[6]} cdp=${SOURCES[7]}
		''')

	#convert to SU format
	Flow (fileOut2, [fileIn1, fileOut1],
		'''
		suwrite tfile=${SOURCES[1]} endian=0
		''')

#-------------------------------------------------------------------------------
# Create Shot-Gathers Headers (FIXED SPREAD)
def CreateSGHeaders_FIXED (fileOut1='SGhead', fileOut2='SGhead.su', fileIn1='geomGrid', fileIn2='SGMap', **params) :

	#create the header fields
	if params['ns'] == 1:
		Flow ('tracl', fileIn2, 'math output=x1                           | dd type=int' % params )
		Flow ('selev', fileIn2, 'math output="-1*%(dz)f"                  | dd type=int' % params )
		Flow ('gelev', fileIn2, 'math output="-1*%(dz)f"                  | dd type=int' % params )
		Flow ('sx',    fileIn2, 'math output=%(sht0)f                     | dd type=int' % params )
		Flow ('gx',    fileIn2, 'math output="%(geo0)f+%(geoint)f*(x1-1)" | dd type=int' % params )
		Flow ('ofst',  'sx gx', 'add scale=-1,1 ${SOURCES[1]}')
	else:
		Flow ('tracl', fileIn2, 'math output="x1+%(ng)d*(x2-1)"                             | dd type=int' % params )
		Flow ('selev', fileIn2, 'math output="-1*%(dz)f"                                    | dd type=int' % params )
		Flow ('gelev', fileIn2, 'math output="-1*%(dz)f"                                    | dd type=int' % params )
		Flow ('sx',    fileIn2, 'math output="%(sht0)f+%(shtint)f*(x2-1)"                   | dd type=int' % params )
		Flow ('gx',    fileIn2, 'math output="%(geo0)f+%(geoint)f*(x1-1)" 					| dd type=int' % params )
		Flow ('ofst',  'sx gx', 'add scale=-1,1 ${SOURCES[1]}')

	#create header file
	Flow (fileOut1, [fileIn1, 'tracl', 'selev', 'gelev', 'sx', 'gx', 'ofst'],
		'''
		segyheader tracl=${SOURCES[1]} selev=${SOURCES[2]}  gelev=${SOURCES[3]}
			      sx=${SOURCES[4]}    gx=${SOURCES[5]} offset=${SOURCES[6]}
		''')

	#convert to SU format
	Flow (fileOut2, [fileIn1, fileOut1],
		'''
		suwrite tfile=${SOURCES[1]} endian=0
		''')

#-------------------------------------------------------------------------------
# Create Towed Streamer Source Trces
def CreatTowedArray (fileOut='waveletSG.su', fileIn1='wavelet_base.su', fileIn2='SGhead.su') :

	Flow (fileOut, [fileIn1, fileIn2],
		'''
		towed_array src=${SOURCES[0]} data=${SOURCES[1]} towed=${TARGETS[0]}
		''', stdin=0, stdout=0)


#*************************************************************************************************#
#******************************************* MODELING ********************************************#
#*************************************************************************************************#

# Non-Linearized Modeling
def Modeling_acd (fileOut='bornSG.su', fileIn1='waveletSG.su', fileIn2='csqModel', fileIn3='SGhead.su', **params) :

	Flow (fileOut,[fileIn1, fileIn2, fileIn3],
		'''
		/bin/cp ${SOURCES[2]} $TARGET &&
		acd deriv=%(deriv)d adjoint=%(adj)d order=%(ord)d cmin=%(cmin)f cmax=%(cmax)f
		source=${SOURCES[0]} csq=${SOURCES[1]} data=$TARGET
		''' % params, 
		stdin=0,stdout=-1,workdir='bornSG.work')

#-------------------------------------------------------------------------------
# Linearized Modeling (Born)
def ModelingBorn_acd (fileOut='bornSG.su', fileIn1='waveletSG.su', fileIn2='csqModelSm', fileIn3='csqPerturb', fileIn4='SGhead.su', **params) :

	Flow (fileOut,[fileIn1, fileIn2, fileIn3, fileIn4],
		'''
	  	/bin/cp ${SOURCES[3]} $TARGET &&
		acd deriv=%(deriv)d adjoint=%(adj)d order=%(ord)d cmin=%(cmin)f cmax=%(cmax)f
		source=${SOURCES[0]} csq=${SOURCES[1]} csq_d1=${SOURCES[2]} data=$TARGET
		dump_lda=%(dump_lda)d dump_term=%(dump_term)d''' % params, 
		stdin=0,stdout=-1,workdir='bornSG.work')


#*************************************************************************************************#
#******************************************* MIGRATION *******************************************#
#*************************************************************************************************#

# RTM
def RTM_acd (fileOut='Mig', fileIn1='bornSG.su', fileIn2='csqModelSm', fileIn3='waveletSG.su' , fileIn4='csqPerturbExtTransp', **params) :

	Flow(fileOut,[fileIn1, fileIn2, fileIn3, fileIn4],
     		'''
     		scale < ${SOURCES[3]} > $TARGET dscale=0.0 &&
     		acd deriv=%(deriv)d adjoint=%(adj)d nsnaps=%(nsnaps)d order=%(ord)d cmin=%(cmin)f cmax=%(cmax)f 
     		data=${SOURCES[0]} csq=${SOURCES[1]} source=${SOURCES[2]} csq_b1=$TARGET sampord=1
     		''' % params,
		stdin=0,stdout=-1,workdir='migSG.work')


