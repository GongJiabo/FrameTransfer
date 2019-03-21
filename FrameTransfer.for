	Program FrameTransfer
      implicit none

	  integer  ::  year, month, day, hour, minute, types
	  real*8   ::  second
 	  real*8  ::  x_input, y_input, z_input
     .              ,x_output,y_output,z_output,tx,ty,tz,ox,oy,oz
     
C          ����Ҫת�������������Ϣ   
	    x_input =4243518.86441863d0      
	    y_input =-10079880.92486593d0
	    z_input =-4943382.36494401d0  
	    
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                        �ο���׼                                             C
C                                                                                             C    
C               epoch1                  2008.0                 7.0                 1.0        C
C               epoch2                     0.0                 0.0          0.00000000        C
C          J2000  pos       4243518.86441863d0-10079880.92486593d0 -4943382.36494401d0        C
C                 vel       2869.15380529916d0  3125.81235293791d0 -3956.00353758515d0        C
C          ECF    pos      10633969.92409196d0  2561753.18400690d0 -4940099.72532129d0        C
C                 vel      -2433.43320886780d0  2564.67135150358d0 -3953.46890766155d0        C
C                                                                                             C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	  
	  write(*,*)  x_input,y_input,z_input
	  
  	  call frametrans(2008,7,1,0,0,0.0, 1 ,x_input,y_input,
     .          z_input,x_output,y_output,z_output)
         
 	  tx = x_output
 	  ty = y_output
 	  tz = z_output
 
 	  write(*,*)  tx,ty,tz    
C 	  call frametrans(2008,7,1,12,0,0.0, 0 ,tx,ty,tz
C    .                ,ox,oy,oz) 
C  	  write(*,*) ox,oy,oz
      end program

	subroutine frametrans(
     .    year,month,day             !�۲����ڣ��������
     .   ,hour,minute,second         !�۲�ʱ�䣬�������
     .   ,types                      !ת�����ͣ�0Ϊ�ӵع�ϵת��������ϵ��1Ϊ�ӹ���ϵת�����ع�ϵ
     .	 ,x_input,y_input,z_input    !��ת��������ֵ�� �������
     .   ,x_output,y_output,z_output !����ת���������ֵ���������
     .                       )
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                         yang.z    2008.03.21                         cc
cc      ������ѡת�����Ͷ�����������ڵع�ϵ�͹���ϵ֮�����ת����      cc
cc      ����ϵ�ο���ԪΪJ2000.0��2000��1��1��12ʱ����ת��ʱ��ԪΪ       cc
cc      ����Ĺ۲�ʱ�䣬�㷨ͨ�����ֵ㡢�������������߼��������        cc
cc      ���¶�����ʵ������ת�����μ������ĵ�����ѧԭ���߲�����      cc
cc      ��GPS��Զ�λ����ѧģ�͡� κ���䣬��ï��                        cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc	
	implicit none

	integer   ::  year, month, day, hour, minute, types
	real*8    ::  second
	real*8    ::  x_input,y_input,z_input,posin(3)
     .              ,x_output,y_output,z_output,posout(3)

      integer ::  yr,mo
      real*8  ::  JD, MJD            !����ʱ��������ա�Լ��������
	real*8  ::  t_TDB              !�۲�ʱ�̼��ο���Ԫ��J2000.0��������������

	integer ::  i,j,k              !ѭ������

	real*8  ::  kesiA,thitaA,zA    !����
	real*8  ::  epsilonA           !�Ƴཻ�ǵı仯  

	real*8  ::  deltaphi
     .             ,deltepthilon

	integer ::  KL(106)            !  L��L'��F��D��OMG��ϵ��
     .            ,KLP(106)
     .            ,KF(106)
     .            ,KD(106)
     .            ,KOMG(106)

	real*8  ::  AI(106),AAI(106)   !�ƾ��¶����仯��
     .            ,BI(106),BBI(106)  !�����¶����仯��

	real*8  ::  m_L, m_LP, m_F, m_D, m_OMG
	real*8  ::  P(3,3), N(3,3), R(3,3), W(3,3)
	real*8  ::  PRZZ(3,3), PRYTHITA(3,3), PRZKESI(3,3)
     .            ,NRXEPTILON(3,3), NRZPHI(3,3), NRXEPTILON2(3,3)
     .            ,BUFFERP(3,3), BUFFERN(3,3),WRXYP(3,3),WRYXP(3,3)
     .            ,BUFPN(3,3),BUFPNR(3,3),CTSTOCIS(3,3),CISTOCTS(3,3)
     .            ,BUFPN1(3,3),PUSAIN(3,3)
     
	real*8  ::  GST        !tʱ�̸��������Ӻ���ʱ��SOFA��
     .           , GAST       !tʱ�̸��������Ӻ���ʱ(�̲�)
     .           , GMST       !��������ƽ����ʱ
     .           , GMST0      !����ʱΪ0h�ĸ�������ƽ����ʱ
     .           , gama       !����ʱ������ʱ֮��
     .           , deltaUT    !UT1-UTC
     .           , UTC        !Э������ʱ
     .           , modangl
     
    
      double precision :: pi
     .                    ,  sec2rad
     
	integer ::  control
	integer ::  mm, dd
	real*8  ::  XPOLE, YPOLE

	character(len = 4)  :: t_year
	character(len = 80) :: buffer

	pi = 3.141592653589793238462643
      sec2rad=4.84813681109536e-6
        
	  data KL /0,0,0,0,0,1,0,0,1,0,1,0,-1,1,0,-1,-1,1,2,-2,0,2,2,1,0
     .          ,0,-1,0,0,-1,0,1,0,2,-1,1,0,0,1,0,-2,0,2,1,1,0,0,2,1,1
     .          ,0,0,1,2,0,1,1,-1,0,1,3,-2,1,-1,1,-1,0,-2,2,3,1,0,1,1,1
     .          ,0,0,0,1,1,1,1,2,0,0,-2,2,0,0,0,0,1,3,-2,-1,0,0,-1,2,2
     .          ,2,2,1,-1,-1,0/

	  data KLP /0,0,0,0,1,0,1,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
     .           ,0,0,2,2,0,1,0,-1,0,0,0,-1,0,1,1,0,0,0,0,0,0,-1,0,-1,0
     .           ,0,1,0,0,1,1,-1,-1,-1,-1,0,0,0,0,0,0,-2,0,0,0,1,0,0,0,1
     .           ,1,1,1,0,0,0,0,0,0,0,0,0,-1,0,0,1,1,0,0,0,0,1,0,1,0
     .           ,0,0,-1,0,-1,0/

	  data KF /0,2,2,0,0,0,2,2,2,2,0,2,2,0,0,2,0,2,0,2,2,2,0,2,2
     .          ,2,2,0,2,0,0,0,0,-2,2,2,2,2,0,2,0,0,2,0,2,0,2,2,0,0
     .          ,0,0,-2,0,2,0,0,2,2,2,2,2,2,2,0,2,2,0,0,0,2,2,0,2,0
     .          ,0,2,-2,-2,-2,2,0,0,2,2,2,2,2,-2,4,0,2,2,2,0,-2,2,4,0,0
     .          ,2,-2,0,0,0,0/

	  data KD /0,-2,0,0,0,0,-2,0,0,-2,-2,-2,0,0,2,2,0,0,-2,0
     .          ,2,0,0,-2,0,-2,0,0,-2,2,0,-2,0,0,2,2,0,2,-2,0
     .          ,2,2,-2,2,-2,-2,-2,0,0,-1,1,-2,0,-2,-2,0,-1,2,2,0
     .          ,0,0,0,4,0,-2,-2,0,0,0,0,1,2,2,-2,2,-2,2,2,-2
     .          ,-2,-4,-4,4,-1,4,2,0,0,-2,0,-2,-2,2,0,2,0,0,-2,2
     .          ,-2,0,-2,1,2,1/

	  data KOMG /1,2,2,2,0,0,2,1,2,2,0,1,2,1,0,2,1,1,0,1,2,2,0,2,0
     .            ,0,1,0,2,1,1,1,1,0,1,2,2,1,0,2,1,1,2,0,1,1,1,1,0,0
     .            ,0,0,0,1,1,0,0,2,2,2,2,2,0,2,2,1,1,1,1,0,2,2,1,1,1
     .            ,0,0,0,0,0,0,0,0,2,2,2,2,1,1,2,2,2,2,2,2,1,1,2,0,0
     .            ,1,1,0,1,1,0/

	  data AI/-17.1996,-1.3187,-0.2274, 0.2062, 0.1426, 0.0712,-0.0517
     .          ,-0.0386,-0.0301, 0.0217,-0.0158, 0.0129, 0.0123, 0.0063
     .          , 0.0063,-0.0059,-0.0058,-0.0051, 0.0048, 0.0046,-0.0038
     .          ,-0.0031, 0.0029, 0.0029, 0.0026,-0.0022, 0.0021, 0.0017
     .          ,-0.0016, 0.0016,-0.0015,-0.0013,-0.0012, 0.0011,-0.0010
     .          ,-0.0008,-0.0007,-0.0007,-0.0007, 0.0007,-0.0006,-0.0006
     .          , 0.0006, 0.0006, 0.0006,-0.0005,-0.0005,-0.0005, 0.0005
     .          ,-0.0004,-0.0004,-0.0004, 0.0004, 0.0004, 0.0004,-0.0003
     .          ,-0.0003,-0.0003,-0.0003,-0.0003,-0.0003,-0.0003, 0.0003
     .          ,-0.0002,-0.0002,-0.0002,-0.0002,-0.0002, 0.0002, 0.0002
     .          , 0.0002, 0.0002,-0.0001,-0.0001,-0.0001,-0.0001,-0.0001
     .          ,-0.0001,-0.0001,-0.0001,-0.0001,-0.0001,-0.0001,-0.0001
     .          ,-0.0001,-0.0001,-0.0001,-0.0001,-0.0001, 0.0001, 0.0001
     .          , 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001
     .          , 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001
     .          , 0.0001/

	  data AAI /-0.01742,-0.00016,-0.00002, 0.00002,-0.00034,0.00001
     .           , 0.00012,-0.00004,     0.0,-0.00005,     0.0,0.00001
     .           ,     0.0, 0.00001,     0.0,     0.0,-0.00001,    0.0
     .           ,     0.0,     0.0,     0.0,     0.0,     0.0,    0.0
     .           ,     0.0,     0.0,     0.0,-0.00001, 0.00001,    0.0
     .           ,     0.0,     0.0,     0.0,     0.0,     0.0,    0.0
     .           ,     0.0,     0.0,     0.0,     0.0,     0.0,    0.0
     .           ,     0.0,     0.0,     0.0,     0.0,     0.0,    0.0
     .           ,     0.0,     0.0,     0.0,     0.0,     0.0,    0.0
     .           ,     0.0,     0.0,     0.0,     0.0,     0.0,    0.0
     .           ,     0.0,     0.0,     0.0,     0.0,     0.0,    0.0
     .           ,     0.0,     0.0,     0.0,     0.0,     0.0,    0.0
     .           ,     0.0,     0.0,     0.0,     0.0,     0.0,    0.0
     .           ,     0.0,     0.0,     0.0,     0.0,     0.0,    0.0
     .           ,     0.0,     0.0,     0.0,     0.0,     0.0,    0.0
     .           ,     0.0,     0.0,	   0.0,		0.0,	 0.0,	 0.0
     .		   ,	 0.0,	  0.0,	   0.0,		0.0,	 0.0,	 0.0
     .		   ,	 0.0,	  0.0,	   0.0,     0.0/

	  data BI /9.2025,0.5736,0.0977,-0.0895,0.0054,-0.0007,0.0224,0.02
     .          ,0.0129,-0.0095,-0.0001,-0.007,-0.0053,-0.0033,-0.0002
     .          ,0.0026,0.0032,0.0027,0.0001,-0.0024,0.0016,0.0013
     .          ,-0.0001,-0.0012,-0.0001,0.0,-0.001,0.0,0.0007,-0.0008
     .          ,0.0009,0.0007,0.0006,0.0,0.0005,0.0003,0.0003,0.0003
     .          ,0.0,-0.0003,0.0003,0.0003,-0.0003,0.0,-0.0003,0.0003
     .          ,0.0003,0.0003,0.0,0.0,0.0,0.0,0.0,-0.0002,-0.0002,0.0
     .          ,0.0,0.0001,0.0001,0.0001,0.0001,0.0001,0.0,0.0001
     .          ,0.0001,0.0001,0.0001,0.0001,-0.0001,0.0,-0.0001,-0.0001
     .          ,0.0,0.0001,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0
     .          ,0.0001,0.0,0.0,0.0,0.0,0.0,-0.0001,0.0,-0.0001,-0.0001
     .          ,0.0,0.0,0.0,0.0,0.0,-0.0001,0.0,0.0,0.0,0.0,0.0/

	  data BBI /0.00089,-0.00031,-0.00005,0.00005,-0.00001,0.0
     .           ,-0.00006,0.0,-0.00001,0.00003,0.0,0.0,0.0,0.0,0.0,0.0
     .           ,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0
     .           ,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0
     .           ,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0
     .           ,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0
     .           ,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0
     .           ,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0
     .           ,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0/

cc    ��������ʱ����������ա�����������
        if(month.LE.2) then
          yr = year - 1
	    mo = month + 12
	  else
          yr = year
          mo = month
        end if
 
  	  JD = int(365.25*yr) + int(30.6001*(mo+1)) + day + hour/24.d0
     .     + minute/1440.d0 + second/86400.d0 + 1720981.5
	  MJD = JD - 2400000.5
	  t_TDB = (JD - 2451545.0)/36525.0d0
cc    ��������ƽ�����m_L��̫��ƽ�����m_LP������ƽ�����m_F�������뿪̫����ƽ����m_D��
cc        ����������ƽ�ƾ�m_OMG����λΪ�ȣ�
C  	  m_L = (134.d0+57/60.d0+46.733/3600.d0)
C     .      + (1325*360.d0 + 198.d0+52/60.d0+2.633/3600.d0)*t_TDB
C     .      + (31.31/3600.d0)*t_TDB**2 + (0.064/3600.d0)*t_TDB**3
C
C         m_LP = (357.d0+31/60.d0+39.804/3600.d0)
C     .       + (99*360.d0 + 359.d0+3/60.d0+1.244/3600.d0)*t_TDB
C     .       - (0.577/3600.d0)*t_TDB**2 - (0.012/3600.d0)*t_TDB**3
ccc
C 	  m_F = (93.d0+16/60.d0+18.877/3600.d0)
C     .      + (1342*360.d0 + 82.d0+1/60.d0+3.137/3600.d0)*t_TDB
C     .      - (13.257/3600.d0)*t_TDB**2 + (0.011/3600.d0)*t_TDB**3

C	  m_D = (297.d0+51/60.d0+1.307/3600.d0)
C     .      + (1236*360.d0 + 307.d0+6/60.d0+41.328/3600.d0)*t_TDB
C     .      - (6.891/3600.d0)*t_TDB**2 + (0.019/3600.d0)*t_TDB**3
C 
C 	  m_OMG = (125.d0+2/60.d0+40.28/3600.d0)
C     .        - (5*360.d0 + 134.d0+8/60.d0+10.539/3600.d0)*t_TDB
C     .        + (7.455/3600.d0)*t_TDB**2 + (0.008/3600.d0)*t_TDB**3

cc    �����ĵ�����ѧԭ���и�¼F��ʽ
	
 	  m_L   = 134.96d0 + 13.064993d0*(MJD - 51544.5)
 	  m_LP  = 357.53d0 + 0.985600d0 *(MJD - 51544.5)
 	  m_F   = 93.27d0  + 13.22935d0 *(MJD - 51544.5)
 	  m_D   = 297.85d0 + 12.190749d0*(MJD - 51544.5)
 	  m_OMG = 125.04d0 - 0.052954d0 *(MJD - 51544.5)

cc    ���ƾ��¶��ͽ����¶�����λΪ���룩
	  deltaphi = 0.d0
	  deltepthilon = 0.d0
	  do i = 1,106
	    modangl = mod((KL(i)*m_L + KLP(i)*m_LP + KF(i)*m_F
     .                    + KD(i)*m_D + KOMG(i)*m_OMG), 360.d0)
		deltaphi = deltaphi + (AI(i) + AAI(i)*t_TDB)
     .             * dsind(modangl)
		deltepthilon = deltepthilon + (BI(i) + BBI(i)*t_TDB)
     .                 * dcosd(modangl)		
	  end do
	  
c	deltaphi = (-6.857 - 0.007*t_TDB)*dsind(mod(m_OMG,360.d0)) 
c     .          + 0.083 * dsind(mod(2*m_OMG,360.d0)) 
c     .          - 0.506 * dsind(mod(2*m_LP,360.d0))
c     .          - 0.081 * dsind(mod(2*m_L,360.d0))
	  
cc    ���������㣨��λΪ���룩
	  kesiA    = 2306.2181 * t_TDB + 0.30188 * t_TDB**2 
     .           + 0.017998 * t_TDB**3
	  thitaA   = 2004.3109 * t_TDB - 0.42665 * t_TDB**2 
     .           - 0.041833 * t_TDB**3
	  zA       = 2306.2181 * t_TDB + 1.09468 * t_TDB**2
     .           + 0.018203 * t_TDB**3
cc    �۲�˲���ƽ�Ƴཻ�ǣ���λΪ���룩
	  epsilonA = 84381.448 - 46.815 * t_TDB - 0.00059 * t_TDB**2
     .           + 0.001813 * t_TDB**3
     
cc	�����¶�����

	  NRXEPTILON(1,1) = 1.d0
	  NRXEPTILON(2,1) = 0.d0
	  NRXEPTILON(3,1) = 0.d0
	  NRXEPTILON(1,2) = 0.d0
	  NRXEPTILON(2,2) = dcosd((epsilonA)/3600.d0)
	  NRXEPTILON(3,2) = -dsind((epsilonA)/3600.d0)
	  NRXEPTILON(1,3) = 0.d0
	  NRXEPTILON(2,3) = dsind((epsilonA)/3600.d0)
	  NRXEPTILON(3,3) = dcosd((epsilonA)/3600.d0)

	  NRZPHI(1,1) = dcosd(-deltaphi/3600.d0)
	  NRZPHI(2,1) = -dsind(-deltaphi/3600.d0)
	  NRZPHI(3,1) = 0.d0
	  NRZPHI(1,2) = dsind(-deltaphi/3600.d0)
	  NRZPHI(2,2) = dcosd(-deltaphi/3600.d0)
	  NRZPHI(3,2) = 0.d0
	  NRZPHI(1,3) = 0.d0
	  NRZPHI(2,3) = 0.d0
	  NRZPHI(3,3) = 1.d0

	  NRXEPTILON2(1,1) = 1.d0
	  NRXEPTILON2(2,1) = 0.d0
	  NRXEPTILON2(3,1) = 0.d0
	  NRXEPTILON2(1,2) = 0.d0
	  NRXEPTILON2(2,2) = dcosd((-epsilonA - deltepthilon)/3600.d0)
	  NRXEPTILON2(3,2) = -dsind((-epsilonA - deltepthilon)/3600.d0)
	  NRXEPTILON2(1,3) = 0.d0
	  NRXEPTILON2(2,3) = dsind((-epsilonA - deltepthilon)/3600.d0)
	  NRXEPTILON2(3,3) = dcosd((-epsilonA - deltepthilon)/3600.d0)

	  BUFFERN = matmul(NRXEPTILON2, NRZPHI)
	  N       = matmul(BUFFERN, NRXEPTILON)

cc    ����������
        
	  PRZZ(1,1) = dcosd(-zA/3600.d0)
	  PRZZ(2,1) = -dsind(-zA/3600.d0)
	  PRZZ(3,1) = 0.d0
	  PRZZ(1,2) = dsind(-zA/3600.d0)
	  PRZZ(2,2) = dcosd(-zA/3600.d0)
	  PRZZ(3,2) = 0.d0
	  PRZZ(1,3) = 0.d0
	  PRZZ(2,3) = 0.d0
	  PRZZ(3,3) = 1.d0

	  PRYTHITA(1,1) = dcosd(thitaA/3600.d0)
	  PRYTHITA(2,1) = 0.d0
	  PRYTHITA(3,1) = dsind(thitaA/3600.d0)
	  PRYTHITA(1,2) = 0.d0
	  PRYTHITA(2,2) = 1.d0
	  PRYTHITA(3,2) = 0.d0
	  PRYTHITA(1,3) = -dsind(thitaA/3600.d0)
	  PRYTHITA(2,3) = 0.d0
	  PRYTHITA(3,3) = dcosd(thitaA/3600.d0)

	  PRZKESI(1,1) = dcosd(-kesiA/3600.d0)
	  PRZKESI(2,1) = -dsind(-kesiA/3600.d0)
	  PRZKESI(3,1) = 0.d0
	  PRZKESI(1,2) = dsind(-kesiA/3600.d0)
	  PRZKESI(2,2) = dcosd(-kesiA/3600.d0)
	  PRZKESI(3,2) = 0.d0
	  PRZKESI(1,3) = 0.d0
	  PRZKESI(2,3) = 0.d0
	  PRZKESI(3,3) = 1.d0

	  BUFFERP = matmul(PRZZ, PRYTHITA)
	  P       = matmul(BUFFERP, PRZKESI)
	  
cc    �������ڶ�ȡ��ӦERP�ļ�����1986.1.1--2008.3.20������ȡXPOLE, YPOLE, UT1-UTC
	  open(unit = 32, file = "temp")
	  write(32,"(I4.4)") year
	  close(32)
	  open(unit = 32, file = "temp")
        read(32,"(A4)") t_year
	  close(32,status = 'delete')
	  open(unit = 101, file = "C04_"//t_year//".ERP")
	  read(101,"(//////A80)") buffer
	  control = 0
	  do while(control.EQ.0)
	    read(101,"(5X,I2,X,I2,7X,D8.5,X,D8.5,X,D9.6)") 
     .        mm, dd, XPOLE, YPOLE, deltaUT
          if((mm.EQ.month).AND.(dd.EQ.day)) then
		   control = -1
		end if
	  end do 
	  close(101)
	  
c       ������ת����

	  WRXYP(1,1) = 1.d0
	  WRXYP(2,1) = 0.d0
	  WRXYP(3,1) = 0.d0
	  WRXYP(1,2) = 0.d0
	  WRXYP(2,2) = dcosd(-YPOLE/3600.d0)
	  WRXYP(3,2) = -dsind(-YPOLE/3600.d0)
	  WRXYP(1,3) = 0.d0
	  WRXYP(2,3) = dsind(-YPOLE/3600.d0)
	  WRXYP(3,3) = dcosd(-YPOLE/3600.d0)

	  WRYXP(1,1) = dcosd(-XPOLE/3600.d0)
	  WRYXP(2,1) = 0.d0
	  WRYXP(3,1) = dsind(-XPOLE/3600.d0)
	  WRYXP(1,2) = 0.d0
	  WRYXP(2,2) = 1.d0
	  WRYXP(3,2) = 0.d0
	  WRYXP(1,3) = -dsind(-XPOLE/3600.d0)
	  WRYXP(2,3) = 0.d0
	  WRYXP(3,3) = dcosd(-XPOLE/3600.d0)

	  W = matmul(WRYXP,WRXYP)

	  GMST0 = 6*3600.d0 + 41*60.d0 + 50.54841 + 8640184.812866*t_TDB
     .        + 0.093104*t_TDB**2 - 6.2d-6*t_TDB**3

	  gama = 1.002737909350795 + 5.9006d-11*t_TDB - 5.9d-15*t_TDB**2

	  GMST = GMST0 + gama*(deltaUT + hour*3600 + minute*60 + second)

        GST =2*pi*GMST/86400.d0 + 
     .       pi*((deltaphi*dcosd(epsilonA/3600.d0)+0.00264*dsind(m_OMG)
     .              + 0.000063*dsind(2*m_OMG))/3600.d0)/180.d0
        GAST=2*pi*GMST/86400.d0+
     .       pi*((deltaphi*dcosd((epsilonA + deltepthilon)/3600.d0))/
     .              3600.d0)/180.d0
                    
c       GAST��z����ת����,����ԭ������GAST,����SOFA���ǰ���GST�Ĺ�ʽת�ģ����Ǹ��������Ӻ���ʱ��
  
	  R(1,1) = dcos(GAST)
	  R(2,1) = -dsin(GAST)
	  R(3,1) = 0.d0
	  R(1,2) = dsin(GAST)
	  R(2,2) = dcos(GAST)
	  R(3,2) = 0.d0
	  R(1,3) = 0.d0
	  R(2,3) = 0.d0
	  R(3,3) = 1.d0

	  BUFPN = matmul(W,R)
	  BUFPNR = matmul(BUFPN,N)
	  
C       ��Z����תdeltaphi*dcosd(epsilonA+deltepthilon)�� ת��TLE(Two-Line Elements)ʱ�õ��˾���,���ļ����� ����ת����ʽ.docx
        
        PUSAIN(1,1)=dcosd((deltaphi*dcosd((epsilonA+deltepthilon)
     .              /3600))/3600)
        PUSAIN(2,1)=-dsind((deltaphi*dcosd((epsilonA+deltepthilon)
     .              /3600))/3600)
        PUSAIN(3,1)=0
        PUSAIN(1,2)=dsind((deltaphi*dcosd((epsilonA+deltepthilon)
     .              /3600))/3600)
        PUSAIN(2,2)=dcosd((deltaphi*dcosd((epsilonA+deltepthilon)
     .              /3600))/3600)
        PUSAIN(1,3)=0
        PUSAIN(2,3)=0
        PUSAIN(3,3)=1
        
        BUFPN1=matmul(N,P)
        
	  CISTOCTS = matmul(BUFPNR,P)
	  CTSTOCIS = transpose(CISTOCTS)

        posin(1) = x_input 
	  posin(2) = y_input
	  posin(3) = z_input

	  if(types.EQ.1) then
	    posout = matmul(CISTOCTS,posin)
	  else
	    posout = matmul(CTSTOCIS,posin)
	  end if

        x_output = posout(1)
        y_output = posout(2)
        z_output = posout(3)
	end subroutine