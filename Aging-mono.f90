program Aging


!implicit none
Integer, Parameter      :: nkmax = 2**12, ntc=10, ntmax=1000, nt2max=1000000, nrmax=250
Real(8), Parameter      ::    R=1.d0, alfa=1.305d0 
Real(8), Parameter      ::  D=1.d0
Real(8)			::pi, fact, a, b, etha, rho, k,gama, kc, suma, tol, tr1, ethavw,lamdak, rhovw, du, tn, dk, dr, jj, u
Real(8)			:: ethamax, errorF, tolF,  h, t, suma1, suma2, Dl,large, error, ethamin, Friccion, ymax, ymin, suma7
Real(8), Dimension(nkmax) :: Sb,landa, landaI, SbI, kap, kaps, gamalist
Real(8), Dimension(ntmax) :: Dz, bdt, wm2
Real(8), Dimension(nrmax) ::rd,gdr
Real(8), Dimension(nt2max)	:: t2, Fs1, bdt2
Real(8)			:: Dt, sumacri,sumagr
Real(8), Dimension(nkmax,ntmax)	:: Fn, Fs
Complex(8), Parameter	:: I = Cmplx(0.d0,1.d0)
Integer		        :: flag, nk, it, nt, l, decim, ml, tim1, mlm, npi, my1, my2, my3, my4, my5, my6, my7,ii
Real(8), Dimension(nkmax) :: rk,Ski, Skf, Sk0, alpha, kd, rk2
 character*4 sub


!npi=32


pi=4.d0*atan(1.)

fact =1.d0/(6.d0*(pi**2))



etha=0.10d0! + (lx-1)*dx0

rho = (6.0d0*etha)/(pi*R**3)


!stop

!open(5,file='Gdraging'//subx//'.dat',status="Unknown")
open(9,file='Skaging.dat',status="Unknown")
open(10,file="sdk_ini.dat")			!%%%% phi0=0.30, T0=10.d0 %%%%
open(11,file='sdk_fin.dat')

do nk=1,nkmax
read(10,*) rk(nk),Ski(nk)
read(11,*) rk2(nk),Skf(nk)

kd(nk)=rk(nk)


alpha(nk)=2.d0*(kd(nk)**2)*D/Skf(nk)					!definición de la alfa(k)


enddo
close(10)
close(11)


dk = kd(2)-kd(1)


!stop


du=0.00001d0!*(9.3999915E-03)

tn=0.0d0
u=0.d0

 nu=6
 !np=9	

! do ii=1,nu
! do jj=1.d0,9.8d0,0.2d0

do u=42.3d0, 1000.d0,du

it=it+1
write(9,*) '  '





  do nk=1,nkmax
					
Sk0(nk)=Ski(nk)*exp(-alpha(nk)*u) + Skf(nk)*(1.d0 - exp(-alpha(nk)*u))			!definicion del factor de estructura dependiente del tiempo

 Sb(nk)=Sk0(nk)			!factores de estructura que le entran al programa dinámico
 SbI(nk)=1.d0/Sb(nk)

 call first_minimum( ymax, ymin) 



!write(9,*) kd(nk), Sb(nk)!(Sb(nk)-Skf(nk))/(Ski(nk) - Skf(nk))


 enddo


 
 call evalua_flag

	print*, "flag=",flag, u

		
!!!!! aqui indicamos los valores de k donde se evaluara el factor de estructura dinámico

	do nk=1, nkmax

	if(kd(nk) >= ymax) then 
	exit
	end if
	
	end do
	
	my1=nk

	do nk=1, nkmax

	if(kd(nk) >= 5.d0) then 
	exit
	end if
	
	end do
	
	my2=nk


	do nk=1, nkmax

	if(kd(nk) >= 6.5d0) then 
	exit
	end if
	
	end do
	
	my3=nk

	do nk=1, nkmax

	if(kd(nk) >= 7.3d0) then 
	exit
	end if
	
	end do
	
	my4=nk

	do nk=1, nkmax

	if(kd(nk) >= 8.d0) then 
	exit
	end if
	
	end do
	
	my5=nk

	do nk=1, nkmax

	if(kd(nk) >= 9.d0) then 
	exit
	end if
	
	end do
	
	my6=nk

	do nk=1, nkmax

	if(kd(nk) >= 10.d0) then 
	exit
	end if
	
	end do
	
	my7=nk



!!!!!! fin de los k's en los cuales se evaluara el factor de estructura dinámico

!stop
	do nk=1, nkmax

	if(kd(nk) >= ymax) then 
	exit
	end if
	
	end do
	
	ml=nk

print*,ml, ymax		!indica el valor de nk en el cual se evaluara la Fs



!kc= 1.118d0*ymax

kc=2.d0*pi*alfa


!   definicion de la funcion interpoladora lamda %%%%%%
do nk=1,nkmax

landa(nk)= 1.d0/(1.d0 + (kd(nk)/kc)**2)
landaI(nk)= 1.d0 + (kd(nk)/kc)**2

!write(12,*) kd(nk),landa(nk)

!write(9,*) kd(nk),  Sb(nk) !, alpha(nk)
enddo


!stop





 if( flag==-1) then
  write(*,*)'la u de arresto es',sngl(u), "gamma=",sngl(gama)

 do nk=1,nkmax
					
Sk0(nk)=Ski(nk)*exp(-alpha(nk)*u) + Skf(nk)*(1.d0 - exp(-alpha(nk)*u))			!definicion del factor de estructura dependiente del tiempo

 Sb(nk)=Sk0(nk)			!factores de estructura que le entran al programa dinámico
 SbI(nk)=1.d0/Sb(nk)

 call first_minimum( ymax, ymin) 



write(9,*) kd(nk), Sb(nk)!(Sb(nk)-Skf(nk))/(Ski(nk) - Skf(nk))


 enddo


    stop
    else
 continue   

 end if

goto 121




 !mlm=etha*1000
 
 !call Indice(mlm,sub) 



     Open(Unit=2, File = 'Fns.dat', Status = "Unknown")
     Open(Unit=3, File = 'Fself.dat', Status = "Unknown")
     Open(Unit=4, File = 'meansquared.dat', Status = "Unknown")
     
write(2,*) '  '
write(3,*) '  '
write(4,*) '  '

     
     

h=1.d-7
friccion=0.0d0
suma6=0.d0

do nt=1, ntc

	t= nt*h
	suma=0.0d0
	do nk=1, nkmax
	
	
	Fn(nk,nt)=(exp(-((kd(nk))**2)*D*(SbI(nk))*t))*Sb(nk)
	Fs(nk,nt)=exp(-((kd(nk))**2)*D*t)
	
	suma= suma + ((kd(nk))**4)*(((Sb(nk)- 1.d0)/Sb(nk))**2)*Fs(nk,nt)*Fn(nk,nt)

	end do

	Dz(nt)= D*(fact/rho)*suma*dk

	write(2,*) t,Fs(ml,nt)	!, Fn(3,nt),  Fs(3,nt)

	write(3,*) t,Fs(ml,nt)


	
!!!!! empezamos a calcular el desplazmiento a tiempos cortos va lineal


bdt(nt)=D-Dz(1)*t


!!!!! empezamos a calcular el desplazmiento a tiempos cortos por el metodo 2

suma6= suma6 + bdt(nt)*h

wm2(nt)=suma6

write(4,*) t,wm2(nt)


end do



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!time medios

tolF=1.d-6

do nt=(ntc + 1), ntmax

t= nt*h

errorF=1.d0

Dz(nt)=Dz(nt-1)
bdt(nt)=bdt(nt-1)


	do while(errorF >tolF)

	suma=0.0d0
	
	
		!$OMP PARALLEL	NUM_THREADS(npi) DEFAULT(NONE)	&
		!$OMP	SHARED(kd,dk,nk,nt,Dz,Fn,Fs,h,Sb,SbI, suma, landa)	&
		!$OMP 	PRIVATE(t,l,suma1,suma2, it, kap, kaps,k)
		!$OMP DO REDUCTION(+:suma)
		
		
		
	do nk=1,nkmax
	!k=k0 + (nk-1)*dk

	suma1=0.0d0
	suma2=0.0d0
	
	do l=2, nt/2
	suma1= suma1 + (Dz(nt - l + 1) - Dz(nt -l))*Fn(nk,l) + (Fn(nk,nt - l + 1) - Fn(nk,nt - l))*Dz(l)
	suma2= suma2 + (Dz(nt - l + 1) - Dz(nt -l))*Fs(nk,l) + (Fs(nk,nt - l + 1) - Fs(nk,nt - l))*Dz(l)
	end do
	
	kap(nk)= 1.d0/h + ((kd(nk))**2)*D*SbI(nk) + Dz(1)*landa(nk)
	kaps(nk)= 1.d0/h + ((kd(nk))**2)*D + landa(nk)*Dz(1)
	


Fn(nk,nt)=(1.d0/kap(nk))*(landa(nk)*Dz(nt)*(Sb(nk) - Fn(nk,1)) - landa(nk)*Dz(nt/2)*Fn(nk,nt/2) + (1.d0/h +&
landa(nk)*Dz(1))*Fn(nk,nt-1) + landa(nk)*Dz(nt-1)*Fn(nk,1)) - (landa(nk)/kap(nk))*suma1

Fs(nk,nt)=(1.d0/kaps(nk))*(landa(nk)*Dz(nt)*(1.d0 - Fs(nk,1)) - landa(nk)*Dz(nt/2)*Fs(nk,nt/2) + (1.d0/h +&
landa(nk)*Dz(1))*Fs(nk,nt-1) + landa(nk)*Dz(nt-1)*Fs(nk,1)) - (landa(nk)/kaps(nk))*suma2

suma= suma + ((kd(nk))**4)*(((Sb(nk)- 1.d0)/Sb(nk))**2)*Fs(nk,nt)*Fn(nk,nt)


	end do!fin de k
	
	!$OMP END DO
	!$OMP END PARALLEL



	Dt=D*(fact/rho)*suma*dk
	errorF=abs((Dz(nt) - Dt)/Dt)

	Dz(nt)= Dt

	end do !fin do while


	write(2,*) t,Fs(ml,nt)
	write(3,*) t,Fs(ml,nt)

!!!!!!!!!!!!!!!!!!!!!!  b(t)





suma7=0.0d0
do l=2, nt/2
suma7= suma7 + (Dz(nt - l + 1) - Dz(nt -l))*bdt(l) + (bdt(nt - l + 1) - bdt(nt - l))*Dz(l)			!mean square displacement
enddo

bdt(nt)= (1.d0/((1.d0/h) + Dz(1)))*(-Dz(nt/2)*bdt(nt/2) - (Dz(nt)- Dz(nt-1))*bdt(1) + ((1.d0/h) + Dz(1))*bdt(nt-1) - suma7)



!!!!! empezamos a calcular el desplazmiento a tiempos cortos por el metodo 2

suma6= suma6 + bdt(nt)*h

wm2(nt)=suma6

write(4,*) t,wm2(nt)


end do !find de nt

!close (unit=2)

!stop


Friccion=Friccion + sum(Dz)*h			!empezamos a calcular el coeficiente de friccion=friccion


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!tiempos largos


it=0.d0
  do decim=1,35
it=it + 1
      
	do nt=1,ntmax/2
	
	Fn(:,nt)=Fn(:,2*nt)
	Fs(:,nt)=Fs(:,2*nt)
	Dz(nt)=Dz(2*nt)
	bdt(nt)=bdt(2*nt)
	end do


	 h=2*h

do nt= (ntmax/2) + 1,ntmax

	t=nt*h

	errorF=1.d0

	Dz(nt)=Dz(nt-1)
	bdt(nt)=bdt(nt-1)

 	do while(errorF >tolF)

	suma=0.0d0



		!$OMP PARALLEL	NUM_THREADS(npi) DEFAULT(NONE)	&
		!$OMP	SHARED(kd,dk,nk,nt,Dz,Fn,Fs,h,Sb,SbI, suma, landa)	&
		!$OMP 	PRIVATE(t,l,suma1,suma2, it, kap, kaps,k)
		!$OMP DO REDUCTION(+:suma)



	do nk=1,nkmax
	!k=k0 + (nk-1)*dk

	suma1=0.0d0
	suma2=0.0d0
 	
	do l=2, nt/2
	suma1= suma1 + (Dz(nt - l + 1) - Dz(nt -l))*Fn(nk,l) + (Fn(nk,nt - l + 1) - Fn(nk,nt - l))*Dz(l)
	suma2= suma2 + (Dz(nt - l + 1) - Dz(nt -l))*Fs(nk,l) + (Fs(nk,nt - l + 1) - Fs(nk,nt - l))*Dz(l)
	
	end do
	



	kap(nk)= 1.d0/h + ((kd(nk))**2)*D*SbI(nk) + Dz(1)*landa(nk)
	kaps(nk)= 1.d0/h + ((kd(nk))**2)*D + landa(nk)*Dz(1)


Fn(nk,nt)=(1.d0/kap(nk))*(landa(nk)*Dz(nt)*(Sb(nk) - Fn(nk,1)) - landa(nk)*Dz(nt/2)*Fn(nk,nt/2) + (1.d0/h +&
landa(nk)*Dz(1))*Fn(nk,nt-1) + landa(nk)*Dz(nt-1)*Fn(nk,1)) - (landa(nk)/kap(nk))*suma1

Fs(nk,nt)=(1.d0/kaps(nk))*(landa(nk)*Dz(nt)*(1.d0 - Fs(nk,1)) - landa(nk)*Dz(nt/2)*Fs(nk,nt/2) + (1.d0/h +&
landa(nk)*Dz(1))*Fs(nk,nt-1) + landa(nk)*Dz(nt-1)*Fs(nk,1)) - (landa(nk)/kaps(nk))*suma2

suma= suma + ((kd(nk))**4)*(((Sb(nk)- 1.d0)/Sb(nk))**2)*Fs(nk,nt)*Fn(nk,nt)

	end do!fin de k


	!$OMP END DO
	!$OMP END PARALLEL



	Dt=D*(fact/rho)*suma*dk
	errorF=abs((Dz(nt) - Dt)/Dt)

	Dz(nt)= Dt


	end do !fin do while

  write(2,*) t,Fs(ml,nt)
  write(3,*) t,Fs(ml,nt)

!!!!!!!!!!!!!!!!!!!!!!  b(t)



suma7=0.0d0
do l=2, nt/2
suma7= suma7 + (Dz(nt - l + 1) - Dz(nt -l))*bdt(l) + (bdt(nt - l + 1) - bdt(nt - l))*Dz(l)			!mean square displacement
enddo

bdt(nt)= (1.d0/((1.d0/h) + Dz(1)))*(-Dz(nt/2)*bdt(nt/2) - (Dz(nt)- Dz(nt-1))*bdt(1) + ((1.d0/h) + Dz(1))*bdt(nt-1) - suma7)



!!!!! empezamos a calcular el desplazmiento a tiempos cortos por el metodo 2

suma6= suma6 + bdt(nt)*h

wm2(nt)=suma6

write(4,*) t,wm2(nt)


end do !find de nt



Friccion=Friccion + h*sum(Dz((ntmax/2) + 1:ntmax)) 


print*, it

end do   !fin de decimaciones
close (unit=2)


Dl=1.d0/(1.d0 + Friccion)

print*, "coficiente de difusión",sngl(Dl)			!me da por salida el coeficiente de difusión a tiempos largos


 Open(Unit=80, File = 'diffusionCoef_un.dat', Status = "Unknown")
 Open(Unit=90, File = 'diffusionCoef_tn.dat', Status = "Unknown")
 Open(Unit=25, File = 'relaxtime_un.dat', Status = "Unknown")
 Open(Unit=26, File = 'relaxtime_tn.dat', Status = "Unknown")
 Open(Unit=27, File = 'sdkvariasktime_un.dat', Status = "Unknown")
 


write(80,320) sngl(u), sngl(Dl)					! coeficiente de difusión como funcion de u
write(90,320) sngl(tn), sngl(Dl)				! coeficiente de difusión como funcion de tn


!!!!!!!!!!!! Aquí se calcula el tiempo de relajación del sistema para cada Fs(ml,nt)


if( flag==-1) then
  write(*,*)'la u de arresto es',u
    stop
    else
 continue   

 end if



 open(7,File='Fns.dat')

tim1=0.d0
do nt2=1,nt2max

read(7,*) t2(nt2), Fs1(nt2)
	
tim1=nt2

	
   if(Fs1(nt2).le. 1.d0/exp(1.d0))then 
   exit
   endif
   
enddo

tr1=t2(tim1)

print*, "tiempo de relajacion=",sngl(tr1)			!me da por salida el tiempo de relajación

write(25,320) sngl(u), sngl(tr1)				! tiempo de relajación como funcion de u

write(26,320) sngl(tn), sngl(tr1)				! tiempo de relajación como funcion de tn

320 format(2f13.7)
close (unit=7)


!!!!!!!!!!! me escribe el archivo de sb(nk,t)

write(27,310) u, tn, Sb(my1), Sb(my2), Sb(my3), Sb(my4), Sb(my5), Sb(my6), Sb(my7)
310 format(10f13.7)

  

tn=tn + (jj*du-u)/Dl			!!!!!secuencia para sacar t_{n+1}=t_{n} + du/Dl, siendo Dl el coeficientte de difusion



! !%%%%
! u= jj*du 
! !%%%% 
! 
! !go to 121
!   
! end do! ! 
! 
! du=du*10.0d0  

121 continue
! 
 enddo



contains




Subroutine Indice(ml,sub)      

     Implicit Integer (i-n)
     Implicit Real (A-H,O-Z)
     Character*10 Doschr
     Character*4 sub
     Character*1 Str1
     Character*2 Str2
     Character*3 Str3
     Character*4 Str4
     Character*5 Str5
     
     If (ml.Lt.10) Then
        Write(Str1,933)ml
        sub=Str1
     Elseif (ml.Lt.100) Then  
        Write(Str2,935)ml
        sub=Str2
     Elseif (ml.Lt.1000) Then
        Write(Str3,937)ml
        sub=Str3
     Elseif (ml.Lt.10000) Then
        Write(Str4,939)ml
        Sub=Str4
     Elseif (ml.Lt.100000) Then
        Write(Str5,940)ml
        sub=Str5
     Endif

933  Format(I1)
935  Format(I2)
937  Format(I3)
939  Format(I4)
940  Format(I5)
     Return
End subroutine




! Subroutine Indice(mlx,subx)      
! 
!      Implicit Integer (i-n)
!      Implicit Real (A-H,O-Z)
!      Character*10 Doschr
!      Character*4 subx
!      Character*1 Str1
!      Character*2 Str2
!      Character*3 Str3
!      Character*4 Str4
!      Character*5 Str5
!      
!      If (mlx.Lt.10) Then
!         Write(Str1,933)mlx
!         subx=Str1
!      Elseif (mlx.Lt.100) Then  
!         Write(Str2,935)mlx
!         subx=Str2
!      Elseif (mlx.Lt.1000) Then
!         Write(Str3,937)mlx
!         subx=Str3
!      Elseif (mlx.Lt.10000) Then
!         Write(Str4,939)mlx
!         Subx=Str4
!      Elseif (mlx.Lt.100000) Then
!         Write(Str5,940)mlx
!         subx=Str5
!      Endif
! 
! 933  Format(I1)
! 935  Format(I2)
! 937  Format(I3)
! 939  Format(I4)
! 940  Format(I5)
!      Return
! End subroutine





subroutine first_minimum(ymax, ymin)		!calcula el maximo y el primer minimo despues del maximo del factor de estructura	lyr
	Implicit none
	
	Double precision  ymax, ymin
	integer nk, imax, imin
	Double precision Smax, Smin
        
       
	Smax= sb(1)						!aqui se calcula el maximo del factor de estructura
	do nk=1, nkmax
	If (sb(nk).ge.Smax) then 
	Smax=Sb(nk)
	imax=nk

	end if
	end do
   	ymax= kd(imax)

	Smin= Sb(imax)
	do nk=imax, nkmax					!aqui calculamos el minimo del factor de estructura :p
	If (Sb(nk).le.Smin) then 
	Smin=Sb(nk)
	imin=nk
	end if
	end do
   	ymin= kd(imin)

	return 
end  subroutine first_minimum





subroutine evalua_flag				
integer	:: nk, it
real(8)	:: k, suma

tol=1.d-6
large=1.d+6
gama=1.d-8
it=0
error=1.d0

do while ( error > tol  .and. gama < large)
	it= it + 1

	sumacri= 0.0d0

	do nk=1, nkmax

				
		sumacri= sumacri + ((kd(nk))**4)*(((Sb(nk) - 1.d0)**2)*((landa(nk))**2))/(((landa(nk))*(Sb(nk)) + &
		((kd(nk))**2)*gama)*(landa(nk) + ((kd(nk))**2)*gama))

	end do

	sumacri=(fact/rho)*sumacri*dk

	error= abs(((1.d0/sumacri) - gama)/gama)

	gama=1.d0/sumacri

	gamalist(it)= gama

!print*, "it, gama =", it, gama

end do


	flag =sign(1.d0, gamalist(it) - 2.d0*gamalist(it - 1) + gamalist(it - 2))

end subroutine evalua_flag


! 
! 
! subroutine Fstruc
! integer ::nk
! real(8) ::k
! 
! 
! a=(1.d0 + 2.d0*ethavw)/((1.d0 - ethavw)**2)
! b=-(3.d0/2.d0)*(R*ethavw)/((1.d0- ethavw)**2)
! 
! 
! 
! Open(Unit=1, File = "bar2.dat", Status = "Unknown")
! 
! Do nk=1, nkmax
! 
! 	k = (k0 + (nk-1)*dk)*lamdak
! 
! 	Sb(nk)= 1.d0/(Q(k)*Q(-k))
! 
! 	SbI(nk)=Q(k)*Q(-k)
! 
! 	landa(nk)= 1.d0/((1.d0 + ((k0 + (nk-1)*dk)/kc)**2))
! 	landaI(nk)= 1.d0 + ((k0 + (nk-1)*dk)/kc)**2
! 
! 	write(1,*) (k0 + (nk-1)*dk), Sb(nk) , landa(nk)
! 
! End Do 
! 
! close(Unit=1)
! 
! end subroutine Fstruc
! 
! 
! complex(8) function Q(k)
! real(8)	:: k
! 
! Q=1.d0 - ((pi*rhovw)/(k**3))*(2.d0*exp(I*k*R)*(b*k + a*(I + k*R)) - I*((2.d0*b*k*(-I + k*R)) + a*(2.d0 + (k**2)*(R**2))))
! 
! end function Q


end program Aging
