!!!!programa para calcular el factor de estructura de 2 yukawas usando la solucion de wertheim para hs
program hard
Implicit None
Integer, Parameter      :: nkmax = 2**12
Real(8)		        :: etha, k, Ts, zeta, rho, errorF,Tmin, Tmax,tolF, kappa,pi, ethavw, lamdak, rhovw
Real(8)			:: gama, kc, suma, tol, large, error, fact, kspin, A, zeta2
Real(8), Parameter      :: sigma= 1.0d0, dk = .01d0, k0 = .01d0, alfa=1.305d0
Real(8), Dimension(nkmax) :: sdkArray,landa, gamalist
Complex(8), Parameter	:: I = Cmplx(0.d0,1.d0)
Integer		        :: nk, nt,nx,flag,l,flagc


pi=4.d0*atan(1.)

gamalist=0.d0


zeta=1.d0
zeta2=0.5d0
A=0.5d0



 etha=0.15d0


!$$$$$$$$$$$$$$$$$$$$$ correccion verlet $$$$$$$$$$$$4

ethavw=etha-(etha**2)/16.0d0	
				!!
lamdak = (ethavw/etha)**(1.0d0/3.0d0)		

rhovw = (6.0d0*ethavw)/(pi*(sigma**3))

!$$$$$$$$$$$$$$$$$$$$$ end $$$$$$$$$$$$4


rho=(6.0d0*etha)/(pi*(sigma**3))


! 
  tolF=1.d-6
  errorF=1.d0

 
 Ts=0.68d0! 
! 	
print*, "etha=",etha, "temperatura=",Ts



call sdkbaxter





Contains

Real(8) Function C0(k)
Real(8)	:: k

!!!!calculo de la C(K) de esfera dura, aqui solo esta la transformada de fourier de la c(r) previamente calculada.



 C0=-(4.d0*pi/(1.d0 - ethavw)**4)*(((1.d0 + 2.d0*ethavw)**2)*(sin(k*sigma)- k*sigma*cos(k*sigma))/(k**3) -&
(6.d0*ethavw/sigma)*((1.d0 + ethavw/2.d0)**2)*((2.d0 - (k*sigma)**2)*cos(k*sigma) - 2.d0 + 2.d0*(k*sigma)*sin(k*sigma))/k**4 +&
(ethavw*(1.d0 +2.d0*ethavw)**2)*(-(24.d0 - 12.d0*(k*sigma)**2 + (k*sigma)**4)*cos(k*sigma) + 4.d0*(6.d0 +&
 k*sigma*((k*sigma)**2 - 6.d0)*sin(k*sigma)))/(2.d0*(sigma**3)*k**6))


End Function 




Real(8) Function bu(k)
Real(8)	:: k

!!!!calculo de la C(K) de yukawa, aqui solo esta la transformada de fourier de la c(r) previamente calculada.


 bu=-(4.d0*pi*sigma/(Ts*k))*(1.d0/((k*sigma)**2 + zeta**2))*((k*(sigma**2))*cos(k*sigma) + (zeta*sigma)*sin(k*sigma)) +&
 (4.d0*pi*sigma*A/k)*(1.d0/((k*sigma)**2 + zeta2**2))*((k*(sigma**2))*cos(k*sigma) + (zeta2*sigma)*sin(k*sigma)) 

End Function 




subroutine sdkbaxter


Open(Unit=1, File = "sdk_Tf.dat", Status = "Unknown")

Do nk=1, nkmax

	k=(k0 + (nk-1)*dk)*lamdak

	sdkArray(nk)= 1.d0/(1.d0 - rhovw*C0(k) + rho*bu(k))
	
	landa(nk)= 1.d0/((1.d0 + ((k0 + (nk-1)*dk)/kc)**2))
 	
    	write(1,*) sngl((k0 + (nk-1)*dk)),sngl(sdkArray(nk))

End Do 


close(Unit=1)




end subroutine sdkbaxter



end program hard
