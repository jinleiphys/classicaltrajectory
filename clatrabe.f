        	IMPLICIT REAL*8(A-H,O-Z)

	DIMENSION Y(4),WS(115)
c	REAL Y(4), WS(115)
	COMMON V0,R0,A,ap,at,zp,zt
	EXTERNAL EXTERN
	OPEN(UNIT=8,FILE='TRAJ2.tex',STATUS='NEW',ACCESS='SEQUENTIAL')
	OPEN(UNIT=9,FILE='coulbarr2.tex',STATUS='NEW',ACCESS='SEQUENTIAL')
	print *,' en minL maxL dl ap at zp zt'
	READ *,IE,iamin,IA, dl,ap,at,zp,zt	

	NT=115*2
	NT1=NT*2
	A=.55
	R0=1.25*(ap**.333333333333+at**.333333333)
	RZ=20.
	RI=1.275*(ap**.333333333333+at**.333333333)
	AI=.3
	amn=1.0432
	amu=amn*ap*at/(ap+at)
       
c  Coulomb length parameter
		DO 3 J=ie,ie
	EF=REAL(J)  ! lab energy 
	 ECM=AT*Ef/(AP+AT)
         AC=ZP*ZT*1.44/2./ECM

      

	IF(J-10)4,4,5   ! if J-10 > 0 then 5, J-10<=0 , then 4 
  4     V0=-35.8-(0.3*EF-3.)
        GO TO 6
  5	V0=-35.8
  6	CONTINUE
	

c	WRITE(8,*)'SET LIMITS Y -15 15'

         ep=-2.39
        be=6.8
        ro=4.85
        GAMMA=8.18
        vm=35.4
c        wv=-gamma*(ef-ep)**4/((ef-ep)**4+vm)
c        WSS=-be*(ef-ep)**2/((ef-ep)**2+ro*ro)+.8313*wv! be/gamma=.8313
c	print *,v0,wv,wss
	DO 2 IJ=iamin,ia,dl!impact parameter
	B1=REAL(IJ)!(SQRT(REAL(IJ)**2+eta**2)+eta)/(SQRT(2*amu*EF)/6.5821) 
       	dd=sqrt(b1**2+ac**2)+ac  ! distance of closest approach
	vc=zp*zt*1.44/sqrt(rz**2+b1**2) ! position of Coulomb barrier 

	 V=SQRT(2*(EF)/(amn*ap)) ! velocity 
	eta=zp*zt*1.44/(6.5821*v)
	TI=RZ/V    ! starting time 
	DELTA=.1
C
C       TRAJECTORY
C
	XT=-TI
	Y(1)=-RZ
	Y(2)=V
	Y(3)=B1
	Y(4)=0.0
        AL=b1*SQRT(2*amu*(EF-vc))/6.5821    ! angular momentum , hbar = 6.5821
	PRINT *,'b D  AL', B1,dd,al
	DO 1 I=2,NT
        CX=ABS(Y(1))
c	IF(CX.GT.20)GO TO 1

	
	CALL INTSTEP(4,DELTA,XT,Y,EXTERN) ! 4 for number of equations, Delta: integration step 
! EXTERN is subroutine 
	XI=Y(1)
	XJ=Y(3)
	R=XI**2+XJ**2
	R1=SQRT(R)
	WRITE(8,*)XI,XJ
c        IF(R1.GT.15)GO TO 1
c	EX=EXP((R1-RI)/AI)
c	WS(I)=WV/(1+EX)+4.*WSS*EX/(1.+EX)**2
  1	CONTINUE
c	print *,ws
	WRITE(8,*)'         '
c	CALL SIMP(WS,SUM,NT,DELTA)
c	T=SUM
c	ts=xt+ti
c	tt=sum/ts
	teta=atan((xj-b1)/xi)*180/3.141593
	thcm=2*atan(ac/b1)*180/3.141593
	tetalab=atan(sin(thcm)/(cos(thcm)+ap/at))*180/3.141593
	print *,'teta  tetalab tetacm',teta, tetalab,thcm
	PRINT *,EF,IJ,v,Y(2)
	write (9,*) b1,teta
  2     CONTINUE
	
  3     CONTINUE
	ir0=2*r0/.1+2
	do 22 i=1,ir0
	ax=-r0+real(i-1)*.1
        ay=sqrt(r0**2-ax**2)
	write(9,*),ax,ay
   22   continue

	STOP
	END

	SUBROUTINE SIMP(W,SUM,NT,DELTA)
	IMPLICIT REAL*8(A-H,O-Z)
!	REAL*8(A-H,O-Z)
	DIMENSION W(201)
!	REAL W(201)
	SOM=0.0
	NTT=(NT-1)/2
	DO 1 IX=1,NTT
	II=IX*2
	IJ=(IX*2)+1
	SOM=SOM+(2.*W(II)+W(IJ))
   1    CONTINUE
	SUM=(2.*SOM+W(1)-W(NT))*DELTA/3.
	RETURN
	END


	SUBROUTINE INTSTEP(NN,HH,XX,Y,EXTERN)
 	IMPLICIT REAL*8(A-H,O-Z)
!	REAL*8(A-H,O-Z)
	DIMENSION Y(4),G(4),O(4),Z(4)
!	REAL Y(4),G(4),O(4),Z(4)	

	N=NN
	H=HH
	H1=0.5*H
	H2=H/6.
	X=XX
	XHH=X+H1
	XH=X+H
	CALL EXTERN(X,Y,G)
	DO 2 J=1,N
    2   Z(J)=Y(J)+H1*G(J)
	CALL EXTERN(XHH,Z,O)
	DO 3 J=1,N
	G(J)=G(J)+2.*O(J)
    3   Z(J)=Y(J)+H1*O(J)
        CALL EXTERN(XHH,Z,O)
	DO 4 J=1,N
	G(J)=G(J)+2.*O(J)
    4   Z(J)=Y(J)+H*O(J)
        CALL EXTERN(XH,Z,O)
	DO 5 J=1,N
    5   Y(J)=Y(J)+H2*(O(J)+G(J))
	XX=XH
		RETURN
	END

	SUBROUTINE EXTERN(X,Y,F)
	IMPLICIT REAL*8(A-H,O-Z)
!	REAL*8(A-H,O-Z)

	DIMENSION Y(4),F(4)
!	REAL Y(4), F(4)
	COMMON V0,R0,A,ap,at,zp,zt
	AMU=1.0432*ap*at/(ap+at)	!REDUCED MASS
	XX=Y(1)
	XY=Y(3)
	AR=XX**2+XY**2
	BR=SQRT(AR)
	RE=(BR-R0)
	IF(RE)1,1,2
  1     vcx=zp*zt*1.44*xx/(r0**3)  ! not the potential but the derivative 
	vcy=zp*zt*1.44*xy/(r0**3)
	go to 4
  2      vcx=zp*zt*1.44*xx/(br**3)
	vcy=zp*zt*1.44*xy/(br**3)
	continue
  4	V0A=V0/A/BR
	EX=EXP((BR-R0)/A)	
	F(1)=Y(2)
	F(3)=Y(4)
	F(2)=V0A*XX*EX/(1.+EX)**2+vcx  ! the derivate of the potential. with nuclear potential 
      !$ +(V0A*XX*EX/(1.+EX)**2)*xx/(br**3)
      !1 +
        F(4)=V0A*XY*EX/(1.+EX)**2+vcy
      !$ +(V0A*XX*EX/(1.+EX)**2)*xy/(b**3)
      !1
	RETURN
	END
