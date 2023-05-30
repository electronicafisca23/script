      program ConductividadElectricaGaAs
      implicit none
      real*8 volt, a, W, Acu, C, q, m0, mefp, mefn, k, Temp, Temps
      real*8 Ecmax, Ecmin, Ef, d, E, vp, vn, vs, tn,tp, mup, mun, Ek,n
      real*8 j,js,Nd,Ndp,Ed,ns
      integer i,l
      real*8 Egrid1, DOS, ge, integral
      real*8 f, fs, Egrid2,sum, sums, Ecult,Ev, t


      !DECLARAR VARIABLES INICIALES

      !parametros de los pozos de potencial en amstrong
      a = 10.0
      W = 20.0
      !Reescalar a metros
      Acu = W*10.0**(-10.0)
      !Magnitudes varias
      d =0.10 !longitud material en m
      C = 100.0 !n§ de capas que atraviesa
      q = -1.6E-19  !carga electr¢n
      m0 = 9.109E-31 !masa electr¢n en kg
      mefn = 1.52158*m0  !masa efectiva electrones en bc
      mefp=-2.90796*m0  !masa efectiva huecos en bv
      k = 8.617332385E-5 !cte boltzmann en eV/K
      Temp = 300.0 !K
      Ecmax = 0.739 !eV   Conduccion
      Ecmin = 0.699 !eV   Conduccion
      Ev=0.196 !eV Valencia
      Ecult = 6.0171 !Esta es la utlima energ¡a de mi archivo DOS para poder hacer el bucle eV
      Ed=Ecmin-0.01 !Nivel donante eV

      !Calcular mi nivel de Fermi CON
      Call Bisectriz (Ef,Ev,Ecmin,k,Temp)
      !SinDopaje
      !Ef=(Ev+Ecmin)/2.0+3.0/2.0*k*Temp*Log(-mefp/mefn)

      !calcular el numero de impurezas donantes ionizadas
      Nd= 10.0**(3.0)
      Ndp=Nd/(1.0+0.5*exp((Ef-Ed)/(k*Temp)))

      !Abrir archivos
      open(11,file="resultados.dat")
      open(12,file="superconductor.dat")
      open(14,file="parametros.dat")
      write(14,*) "Parametros iniciales"
      write(14,*) "A A","W A","T K", "Ev eV","Ec eV","Ecmax eV","Ef eV"
      write(14,*) A, W, Temp, Ev, Ecmin, Ecmax, Ef
      write(14,*) "mefn m0", "mefp m0"
      write(14,*) mefn, mefp
      !datos del bucle0
      write(14,*) "n m^-1 ","ndp m^-1 ","E V/m ","tn s","tp s","vn","vp"


      ! C lculos
       do i=0,10
       !abrir archivos DOS

        open(10,file="1DDOSplot.txt")

        !calcular todo lo necesario para cada voltaje
        volt = (i)*0.02+0.5
        E = volt/d !CAMPO ELECTRICO
        vn = sqrt(-(Acu*C*q*E)/mefn) !VELOCIDADES
        vp= sqrt((Acu*C*q*E)/mefp)
        vs= sqrt(-q*Volt*2.0/m0) !SUPERCONDUCTOR
        tp = C*Acu/vp  !TIEMPO DE RELAJACIOO EN s
        tn = C*Acu/vn
        !mun = (q*t)/mefn
        !mup = (-q*t)/mefp

        !declarar variables para la integral de n
        sum=0.0
        Egrid1=Ecmin
        Egrid2=0.0

        !Leer encabezado

        read(10,*)
        read(10,*)
        read(10,*)
        read(10,*)

        ! Calcular n

        do while (Egrid2.lt.Ecmax)

           read(10,*) DOS, Egrid2

           if ((Egrid2.gt.Ecmin).and.(Egrid2.le.Ecult))then

              f = 1/(1+exp((Egrid2-Ef)/(k*Temp))) !Estadisitica de Fermi
              ge = DOS*4.72686e+10 !Densidad de estados del archivo
              sum = sum+f*ge*(Egrid2-Egrid1)
              Egrid1=Egrid2

           endif

        enddo

        integral=sum
        n = integral
        !CON
        j = ((n+ndp)*q**2.0*E)*(tn/mefn-tp/mefp)
        !SIN
      !  j = (n*q**2.0*E)*(tn/mefn-tp/mefp)

        !guardar datos en archivos
        write(14,*) n, ndp, E, tn, tp, vn, vp
        write(11,*) Volt, j
        !Ahora para el estado de superconductividad
        ns=n
        js=-q*(ns+ndp)*vs

        write(12,*) Volt, js
        close(10)
      end do

      pause
      stop
      end
      
      !Algoritmo de m‚todos numericos que sirve para resolver ecuaciones
      !no lineales. En este caso para determinar la energ¡a de fermi
      !usando las ecuaciones de neutralidad de carga
      
      Subroutine Bisectriz (Ef,Ev,Ec,k,T)
      implicit none
      integer n
      real*8 a,b,c,cmn,e,f,xpn,x
      real*8 Ev,Ec,fa,fb,fc,k,T, Ef
      
      !error
      e=0.000001

      !intervalo inicial
      a=Ev+0.01
      b=Ec-0.01
      n=0.d0
      
      !llamar a las funciones
      fb=f(b,Ec,Ev,k,T)
200   fa=f(a,Ec,Ev,k,T)
      c=(a+b)/2.0
      fc=f(c,Ec,Ev,k,T)

      !El algoritmo se podr¡a hacer con un do while
      if ((b-a).lt.e) then
        write(*,*) "La Enegr¡a de Fermi es:  ",c ," eV"
      else if(fa*fc.lt.0.d0)then
          b=c
          goto 200
      else if (fa*fc.gt.0.d0)then
          a=c
          goto 200
      endif

      Ef=c

      return
      end

      !Declaro la funcion que tengo que igualar a 0

      real*8 function f(Ef,Ec,Ev,k,T)
      implicit none
      real*8 Ef
      real*8 n, p, Na, Nd, Nc, Nv, Nam, Ndp, T,k
      real*8 Ec, Ev, Ea, Ed

      !Datos de mi material GaAs
      Na= 5.0*10**(1.0)
      Nd= 10.0*10**(3.0)
      Nc= 4.7*10**(4.0)
      Nv= 7.0*10**(5.0)
      Ea=Ev+0.01
      Ed=Ec-0.01

      
      !Calculos
      n=Nc*(T/300.0)**(3.0/2.0)*exp(-(Ec-Ef)/(k*T))
      p=Nv*(T/300.0)**(3.0/2.0)*exp(-(Ef-Ev)/(k*T))

      Nam=Na/(1.0+0.5*exp((Ea-Ef)/(k*T)))
      Ndp=Nd/(1.0+0.5*exp((Ef-Ed)/(k*T)))

      !condicion de neutralidad de carga
      f=p+ndp-n-nam

      return
      end

