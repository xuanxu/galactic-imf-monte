c        federico-andrea b.- giacomo g. (viareggio-monterspertoli)
c        programma di evoluzione chimica di un modello galattico
c          cÁlculo delle matrici di produzione degli elementi
c                            qfm(k,i,j)
c                     versione 1.0 - 18.03.96
c
c	mercedes mollá - marta gavilán 
c	versión 2.0 - noviembre 2000 
c
c     mercedes mollá - juanjo bazán
c     versión fi montecarlo 2003
c
c	según ferrini,1992,apj,387,138
c	      portinari,1998,aa,334,505
c
c     versione renzini voli (1980), ar=1.5, sni,e=1
c	comprende la contribución de supernovas del tipo ia y ib
c  ---------------------------------------------------------------
c     Índices de la matriz q(i,j)
c       1 - h     4 - he4    7 - n14    10 - ne    13 - s
c       2 - d     5 - c12    8 - c13    11 - mg    14 - ca
c       3 - he3   6 - o16    9 - n.r.   12 - si    15 - fe
c ----------------------------------------------------------------
c
       real m,mmax,mmin,msep,msn2,minf,msup,cs,feh,f
       integer g, elflag
       dimension vm(800,2),vna(800),vnb(800),et(800)
       dimension tm(80),t(19,80),w(10)
       dimension qm(15,15),qmn(15,15),qmsa(15,15),qmsb(15,15)
       dimension x(6),xx(6)
       character*4 mod
       character*30 fimonte

       common/const/nmas,imax,irid,imax1,jmax,ic
       common/masse/mmax,mmin,msep,msn2,bmin,bmax,rm
       common/sndat/alf,umalf,gamma,costfmu,cs

       common/intrp/beta,erre
       common/abund/x1,x2,x3,x4,x12,x16,x14,x13,x,feh
       common/tab/tm,t
       common/intgr/nw,w
       common/tasas/elflag

c ----declaro las variables necesarias para el montecarlo
       integer ya,seed,totstar
       real montm, montmi, montmf, montns, fimax, montx, monty
       real h
       dimension vmontfi(800),vmontfi2(800)
       common/aaaa/vmontfi,vm,vmontfi2,lmt

c  abundancias si [Fe/H]=no solar
c
c   hacemos la z=0.02
        real z
        z=0.02


        x1=-3.4479*z+0.7740
        x4=2.4479*z+0.2260
        feh=alog10(z/0.02)
      
       x12=x12*(10**feh)
       x16=x16*(10**feh)
       x13=x13*(10**feh)
       x14=x14*(10**feh)       
       do k=1,6
       x(k)=x(k)*(10**feh)       
        end do

     

c  ---------------------------------------------------------------
c   lectura de los parametros del programa (paramt)
       open(10,file='pargra',status='old')
2       format (f6.2)
4       format (i6)
       read(10,4) iw,ic,lm2,lblk
       read(10,2) msep,msn2,bmin,bmax
       read(10,2) alf,gamma,ttot
       close(10)

c  ---------------------------------------------------------------
c   cálculo de las variables auxiliares
       imax1=imax
       if (ic.eq.0) imax1=irid
       tsep=tau(msep)
       delt=tsep/lm2
       delt1=lblk*delt
       lm1=int(1+(lm2*ttot)/(tsep*lblk))
c        if (iw.eq.1) write (6,5) delt,lm2,lm1,lm1*delt1
5       format(' delta t = ',e9.4/' memorie 2-1 = ',2i5/
     1    ' t.tot = ',f7.3)
       bmaxm=bmax/2.
       umalf=1.-alf
       costfmu=(1.+gamma)*2.**(1.+gamma)
       sw=0.0
       do 11 ip=1,nw
        sw=sw+w(ip)
11     end do
        sw=(nw-1)/sw
       do 12 ip=1,nw
        w(ip)=w(ip)*sw
12     end do


c ----------------------------------------------------------------
c   	lectura de datos de las masas eyectadas
      open(10,file='exp',status='old')
c ----------------------------------------------------------------
c     se lee  en el mismo orden que la matriz q. como no se eyecta
c	deuterio, aparece con valor 0 
c       1 - h     4 - he4    7 - n14    10 - ne    13 - s
c       2 - d     5 - c12    8 - c13    11 - mg    14 - ca
c       3 - he3   6 - o16    9 - n.r.   12 - si    15 - fe
c      16 - remanentes   17 - core de helio (no se lee, se calcula)
c	 18 - c13s	19 - n14s
c ----------------------------------------------------------------
c 50 	format(f6.2,f9.5,f8.5,2f8.5,10f8.5,2f9.5)
 50    format(f6.2,16(e11.0))
          j=1
 81    read(10,*,end=80) tm(j),t(1,j),t(3,j),t(4,j),t(5,j),
     1    t(6,j),t(7,j),t(8,j),t(10,j),t(11,j),t(12,j),
     1    t(13,j),t(14,j),t(15,j),t(16,j),t(18,j),t(19,j)
          if(tm(j).ge.100)goto 80
            j=j+1
          go to 81
80     continue
           nmas=j-1
       close(10)
c---------------------------------------------------------


c -------cálculo de la nucleosíntesis explosiva
c
c ----- generación del fichero de intervalos de masas
       open(15,file='intmas')
91     format (3i5/e15.8,f7.2)
       write (15,91) lm2,lblk,lm1,delt,lm1*delt1
93     format(1p,2e14.7,i6)

       minf=mmax
       tsup=0.0
       do 95 i=1,lm2        
	msup=minf
        minf=emme(delt*i)
        if (minf.gt.mmax) minf=mmax
        write(15,93) msup,minf,i
        vm(i,1)=minf
        vm(i,2)=msup
        tinf=tsup
        tsup=delt*i
c        vna(i)=(ennea(tsup)-ennea(tinf))*eta
c        vnb(i)=(enneb(tsup)-enneb(tinf))*eta
c        et(i)=etout(tsup)-etout(tinf)
95     end do
        minf=msep
       do 96 i=1,lm1
        msup=minf
        tinf=tsup
        minf=emme(delt1*i)
        if (minf.gt.msep) minf=msep
        write(15,93) msup,minf,i
        ii=i+lm2
        vm(ii,1)=minf
        vm(ii,2)=msup
        tsup=delt1*i
        if(tsup.le.tinf)tsup=tinf
c        vna(ii)=(ennea(tsup)-ennea(tinf))*eta
c        vnb(ii)=(enneb(tsup)-enneb(tinf))*eta
c        print*,minf,msup,vna(ii),vnb(ii),eta
96     end do
       close(15)
        lmt=lm1+lm2

       totstar=0.21*1e6
       open(unit=1,file='nest.dat')
       write(1,1)totstar
 1     format(1x,I10)
       print*,totstar
       close(1)
     


c ---------- loop de cálculo
       open(20,file='qmsnr')
       open(21,file='fisnr')
       open(13,file='intermedio')
c       open(30,file='sn')     
c       open(30,file='fim')
       write (21,91) lm2,lblk,lm1,delt,lm1*delt1

c ---------- cálculo de qm por snei
       call compq(rm,qmsa,1)
       call compq(rm,qmsb,2)
       g=0
       do 100 im=1,lm2+lm1
        minf=vm(im,1)
        msup=vm(im,2)
        stm=(msup-minf)/(nw-1)
        fik=0.
        fisik=0
        fisiik=0.


c 101    format(I4,1x,6(E9.3,1x))

        do 120 i=1,imax1
          do 120 j=1,jmax
            qm(i,j)=0.
120     end do
c        vnak=1e6*vna(im)
c        vnbk=1e6*vnb(im)
c        etk=1e6*et(im)
c        print*,vnak,vnbk,etk,minf,msup
        if (msup.lt.mmin.or.stm.eq.0.) goto 300

c  ------------- loop de integración a nw puntos
        do 200 ip=1,nw

c---------- busco un flag para el fichero de control
	elflag=ip

c -------------------
          m=msup-stm*(ip-1)
c           m=minf+stm*(ip-1)
          call compq(m,qmn,0)

c----------- cálculo de las funciones fi
          f=totstar*w(ip)*stm
          fm1=f*fi1(m)
          fm12=fm1+(f*fi2(m,0))
          fm2s=f*fi2(m,1)
          fmr=f*fir(m)
          fik=fik+fm12
          fisik=fisik+fm2s/m

          if (m.gt.msn2) fisiik=fisiik+fm1/m
c          if(m.le.msn2)print*,m,fisik,fisiik,fir(m),fmr,fm1,f,montns
c----------- cálculo de los qm
          do 130 i=1,imax1
         
            do 130 j=1,jmax
              qm(i,j)=qm(i,j)+fm12*qmn(i,j)
c quito la parte de las snIb, haciendo vnbk = 0
          if (m.lt.bmaxm) qm(i,j)=qm(i,j)+fmr*qmsa(i,j)

c          if (m.lt.bmaxm) qm(i,j)=qm(i,j)+1.4*(qmsa(i,j)*vnak+
c     1                                qmsb(i,j)*vnbk)
c          if (m.lt.bmaxm) qm(i,j)=qm(i,j)+1.4*(qmsa(i,j)*vnak)
     
130        end do
200       end do
c       if(m.gt.80) print*,'despues',m,fik,fm12,fm1,fi2(m,0)
300     continue

c ------------ fin del loop de integración a nw puntos

c ------------ escritura de datos en los ficheros de salida

 53    format(4f11.3,f11.3)
 52    format(9(E12.6,1x))
        if (m.gt.1) g=1
        write(20,52)((qm(i,j),j=1,jmax),i=1,imax1)
        write(21,53) fik,fisiik,fisik

100     continue

c ------------ fin del loop de lectura y cálculo
c
        close(20)
        close(21)
        close(13)
        close(30)

        end

c =================================================================

c ------------------------------------
      subroutine compq(m,qm,isn)
c -------------------------------------
c        isn=0 normal stars;  isn=1 supernovae ia ; isn=2 supernovae ib
       real m,mmax,mmin,msep,msn2,bmin,bmax,rm,x1x4
       real q4,qc,wc,qns,qc13s,xn,xc13,xc,xo
       integer elflag
       dimension qm(15,15),el(19),csniz02(2,11),csniz0(2,11),
     1 csni(2,11),x(6),xx(6)
       common/const/nmas,imax,irid,imax1,jmax,ic      
       common/masse/mmax,mmin,msep,msn2,bmin,bmax,rm
       common/intrp/beta,erre
       common/abund/x1,x2,x3,x4,x12,x16,x14,x13,x,feh
       common/tasas/elflag
    
       data csniz02/0.,1.1,0.0483,0.,0.143,0.,0.,0.,0.,0.,
     1 0.00202,0.,0.0085,0.,0.154,0.,0.0846,0.,0.0119,0.,
     1 0.626,0.3/

       data csniz0/0.,1.1,0.051,0.,0.133,0.,0.,0.,0.,0.,0.00229,
     1 0.,0.0158,0.,0.142,0.,0.0914,0.,0.0181,0.,0.68,0.3/      
       do 1 j=1,jmax
        do 1 i=1,imax1
          qm(i,j)=0.
1      end do
	if(m.lt.mmin) return
c ----------------------------------------------------
c separo los procedimientos de estrellas de los de SN, y mando la
c parte de SN al final de la subrutina
c -----------------------------------------------------
	if(isn.ne.0) goto 1000

        
c ---------------------------
        call interp (m,el)

c ***********************************************************************
c cambio todo el procedimiento
c ***************************************************************************
c ---- elementos de la matriz q. primero calculo las magnitudes
c      d: remanente	
c      q4: core de helio
c      qc: core de co
c      wc: metales sintetizados nuevos y eyectados
c      qns: core de nitrógeno secundario
c      qc13s: core de carbono13 secundario
c	 xn: abundancia fraccional de nitrógeno en wc
c   	 xc13: abundancia fraccional de carbono13 en wc
c 	 xc: abundancia fraccional de carbono en wc
c	 xo: abundancia fraccional de oxígeno en wc
c	 xx(i): abundancia fraccional de ne, mg, si, s ca fe (i:1..6) en wc
c --------------------------------------------------------------------
	x1x4 = x1+x4
	if(isn.eq.0)d=el(16)
	q4=1-(el(1)/x1)
	qc=(((q4*x1)+x4-el(4))/x1x4)
	wc=qc-d

c --- controlo que wc no sea negativo
	if (wc.lt.0) wc=0

c --- el n y el c secundarios se calculan de forma distinta para
c estrellas  masivas que para intermedias o pequeñas. (hasta 8 masas
c solares).  Como no tengo nada clara la producción de ns y c13s, meto
c una columna más a la tabla de entrada, que será el valor de eyección
c de ns y c13s para estrellas no masivas (renzini) y cero para el
c resto. se correponden con  los elementos el(18) el carbono 13 y el(19)
c el nitrógeno 14
c controlo que la producción secundaria masiva es cero.
c Para que no interpole estre estrellas masivas y no masivas, pongo un control,
c de forma que si es primario, el(18) y el(19) sean cero, y si es secundario
c el(7) y el(8) sean cero, pero solamente entre estrellas masivas.

         if (m.gt.8) then
                 el(7) = 0
                 el(8) = 0
             end if


	qns=(el(19)/(x12+x13+x16))-((1-qc)*x14/(x12+x13+x16))
     1  +qc
	qc13s = (el(18)/x12)-((1-qns)*(x13/x12))+qns
	
	
	if (wc.ne.0) then
           xn=el(7)/(x1x4*wc)
	   xc13=el(8)/(x1x4*wc)
	   xc=el(5)/(x1x4*wc)-(1-qc13s)*x12/(x1x4*wc)	
	   xo=el(6)/(x1x4*wc)-(1-qns)*x16/(x1x4*wc)
	   do 20 i=10,15
		xx(i-9)=(el(i)-((1-d)*x(i-9)))/(wc*x1x4)
20	   end do
	else
	   xn = 0
	   xc13 = 0
	   xc = 0
	   xo = 0
	   do 21 i=1,6
		xx(i) = 0
21         end do
        end if
c------- controlo que no haya ninguna producción negativa
	if (xc.lt.0) xc =0
	if (xo.lt.0) xo =0
	if (xc13.lt.0) xc13 =0
	if (xn.lt.0) xn =0
	do 30 i=1,6
		if (xx(i).lt.0) xx(i)=0
30	end do

c       
c--------------------------------------------------------------
c	q3, w3 y w2 son magnitudes que no aparecen en el artículo de
c	l.portinari, pero sí en el de ferrini, y se usan para el 
c	cálculo con deuterio y helio3. las fórmulas no son las
c	del artículo, sino las del programa antiguo, que no se
c	parecen en nada.
c---------------------------------------------------------------
c  ---------------  cálculo de q3 para todo el rango de masas
        if (m.le.3.) q3 =q4
        if (m.gt.3. .and. m.le. 8.)q3 = 0.282+0.026*m
        if (m.gt.8. .and. m.le.15.) q3 = 0.33+0.02*m
        if (m.gt.15. .and. m.le. 25.) q3 = 0.525+0.007*m
        if (m.gt.25. .and. m.le. 50.) q3 = 0.63+0.00288*m
        if (m.gt.50.) q3 = 0.73+ 0.0008*m

c ---------------  cálculo de w2
        w2=0.0
c ---------------  cálculo de w3
c     schatzman:       
        w32=1.063e-3*(m**(-2))*(1.-d)/x1
c       el(4)=0.0
	if(m.lt.5) then
		if (m.ge.0.8 .and. m.lt.2.) w3=-3.47e-4*m+7.79e-4
		if (m.ge.2. .and. m.lt.3.) w3=-4.43e-5*m+1.74e-4
		if (m.ge.3. .and. m.lt.5.) w3=-1.15e-5*m+7.53e-5
		w3 =w3*(1-d)/x1
	else
		w3=0
	end if
c       print*, m,w3,w2
c---Cambio la parte correspondiente al H y al He.
c          q3=1-(el(2)/x2)   <
c          w3=(el(3)/x1)-((1-q3)*x3/x1)
c Fin del cambio

c------------------------------------------------------------------
c	control de los datos intermedio
c	escribo en el fichero intermedio los datos intermedios
c------------------------------------------------------------------
 	
62     format(f6.2,18e12.3)
       if (elflag.ne.1) goto 65       
        write(13,62) m,w2,w3,q4,qc,d,wc,qns,qc13s,xn,xc13,
     1                    xc,xo,(xx(i),i=1,6)
65     continue

c---------------------------------------------------------------
c definición de los elementos de la matriz q(i,j)
c la parte correspondiente al c12, o16, n14 y c13 depende de
c la masa. para estrellas de más de 8 mo no hay producción secun-
c daria, y los elementos de la matriz entre 5 y 8 solo son no
c nulos en la diagonal. no se corresponde con la matriz de 
c portinari
c---------------------------------------------------------------
	do 22 i=1,15
		do 23 j=1,15
			qm(i,j) = 0
23		end do
22	end do

	qm(1,1)=1-q4-w3
	qm(1,2)=-0.5*(1-d)
        qm(1,3)=0
        qm(1,4)=0
        qm(2,2)=0
	qm(3,1)=w3
	qm(3,2)=1.5*(1-q3)
	qm(3,3)=1-q3
	qm(4,1)=q4-qc
	qm(4,2)=1.5*(q3-qc)
	qm(4,3)=q3-qc      
	qm(4,4)=1-qc         

	qm(5,1)=xc*wc
	qm(5,2)=1.5*xc*wc
	qm(5,3)=xc*wc
	qm(5,4)=xc*wc
	qm(6,1)=xo*wc
	qm(6,2)=1.5*xo*wc
	qm(6,3)=xo*wc
	qm(6,4)=xo*wc
	qm(7,1)=xn*wc
	qm(7,2)=1.5*xn*wc
	qm(7,3)=xn*wc
	qm(7,4)=xn*wc
	qm(8,1)=xc13*wc
	qm(8,2)=1.5*xc13*wc
	qm(8,3)=xc13*wc
	qm(8,4)=xc13*wc
	do 25 i=10,15
	   do 26 j=1,4
	      qm(i,j)=xx(i-9)*wc
		if (j.eq.2) qm(i,j) = qm(i,j)*1.5 
26	   end do
25	end do

	qm(9,5)=wc
	qm(9,6)=wc
	qm(9,7)=wc
	qm(9,8)=wc
	do 28 i=9,15
	   qm(i,i)=1-d
28	end do
c----------------------------------------------------------------------
c	el cálculo de los elementos centrales: c12, n14, c13, o16
c	depende de los elementos secundarios. Cuando se supone
c	que no hay elementos secundarios (todo primario)
c	los valores qns, qc13s son 0, y  sólo queda 1-qc en la diagonal
c       pero eso no ocurre nunca!. 
c ---------------------------------------------------------------------

           qm(5,5)=1-qc13s
	   qm(6,6)=1-qns
	   qm(7,7)=1-qc
	   qm(8,8)=1-qns
	   qm(7,5)=qns-qc
	   qm(7,6)=qns-qc
	   qm(7,8)=qns-qc
	   qm(8,5)=qc13s-qns

c Impongo la condición de que no haya ningún elemento negativo
	do 2020 i=1,15
	    do 2030 j=1,15
	        if(qm(i,j).lt.0) qm(i,j)=0
2030	    end do	
2020	end do
 	return 
c ---------------------------------
c Parte correspondiente a las SN  -
c ---------------------------------

1000    continue

c --------------------------------------------------------------
c  leo las eyecciones de SN en función de la metalicidad. Hay dos
c  DATA, uno para baja metalicidad (Z=0.004, Z=0.0004) y otro 
c  para metalicidad más ata (Z=0.008, Z=0.02 y Z= 0.0317)
c --------------------------------------------------------------
       	if(feh.ge.-0.3)then
	do 200 i=1,2
	   do 200 j=1,11
	      csni(i,j) = csniz02(i,j)
200	end do
	end if

	if(feh.lt.-1.3)then
	do 210 i=1,2
	   do 210 j=1,11
	      csni(i,j) = csniz0(i,j)
210	end do
	end if

c --------------------------------------------------------------
c       aquí se incorporan las eyecciones por sn que se han 
c       obtenido de la lectura del data. 
c --------------------------------------------------------------
       do 70 i=1,15
	  el(i)=0
70     end do
       den=0.99*m
       do 71 i=4,8
         el(i)=el(i)+ csni(isn,i-3)/den
71     end do
       do 72 i=10,15
         el(i)=csni(isn,i-4)/den
72     end do
c ------------------------------------------------
c calculo el valor de la remanente en función de las eyecciones
c -------------------------------------------------
	d=1.4
	do 73 i=1,11
	   d=d-csni(isn,i)
73	end do
c ----------------------------------------------------------------
c hago un procedimiento exclusivo para calcular la matriz q en el
c caso de las sn
c ----------------------------------------------------------------
	do 74 i=5,8
   	   do 75 j=1,4
	      qm(i,j) = el(i)
75	   end do
74	end do
	do 76 i=10,15
   	   do 77 j=1,4
	      qm(i,j) = el(i)
77	   end do
76      end do
	do 79 i=5,8
	   qm(9,i)=0
79      end do
	do 80 i=9,15
	   qm(i,i)=1-d
80      end do

c------------------------------------------------------------------------------
c Normalizo las columnas de la matriz de SN para que la suma de 1-d
c------------------------------------------------------------------------------

c	columna = 0
c	do 84 j=1,4
c	   do 82 i = 1,15
c	      columna = columna+qm(i,j)
c82	   continue
c	   do 84 i=1,15
c	      qm(i,j) = qm(i,j)*(1-d)/columna
c84	   continue
c
 	return
 	end
c
c *****************************************************************

c --- interpolación ------------
      subroutine interp (m,el)
c ------------------------------
       real m
       dimension tm(80),t(19,80),el(19)
       common/tab/tm,t

       k=0
1      k=k+1
       if(m.gt.tm(k)) goto 1
       km1=k-1
       p=(tm(k)-m)/(tm(k)-tm(km1))
       do i=1,19
        d=t(i,k)-t(i,km1)
        el(i)=(t(i,k)-p*d)/m
       end do
       return
       end

c ***************************************************
c       	definición de funciones             *
c ------------------------------------------- 
 
c ----------------------------------------------- 
c  función fi originaria (di ferrini-palla-penco) 
       function fiteorica(m) 
       real m 
       a=alog10(m) 
       fiteorica=2.0865*10**(-sqrt(0.73+a*(1.92+a*2.07)))/m**0.52 
       return 
       end 
c ------------------------------------------- 
c  función fim1 primaria de binaria
       function fim1(m)
       real m,mmax,mmin,msep,msn2
       common/masse/mmax,mmin,msep,msn2,bmin,bmax,rm
       common/sndat/alf,umalf,gamma,costfmu,cs
       common/intgr/nw,w
       dimension w(10)
       fim1=0.
       binf=amax1(bmin,m)
       bsup=amin1(bmax,2.*m)
       stm=(bsup-binf)/(nw-1)
       if (stm.le.0.0) return
c --- integral 
       do 1 i=1,nw
        bm=binf+(i-1)*stm
        fim1=fim1+w(i)*effe(1.-m/bm)*fiteorica(bm)*m/(bm*bm)
1      end do
       fim1=alf*stm*fim1
       return
       end
c -------------------------------------------
c  función fi1 para estrellas normales y primarias de binarias
       function fi1(m)
       real m,mmax,mmin,msep,msn2
       common/masse/mmax,mmin,msep,msn2,bmin,bmax,rm
       common/sndat/alf,umalf,gamma,costfmu,cs
       fi1=fiteorica(m)
       if (m.ge.bmin.and.m.le.bmax) fi1=umalf*fi1
       fi1=fi1+fim1(m)
       return
       end

c -------------------------------------------
c  función fi2 secundaria de binarias (is=1 se dan snei)
       function fi2(m,is)
       real m,mmax,mmin,msep,msn2
       common/masse/mmax,mmin,msep,msn2,bmin,bmax,rm
       common/sndat/alf,umalf,gamma,costfmu,cs
       common/intgr/nw,w
       dimension w(10)
       fi2=0.
       bsup=bmax
       if (is.eq.1) bsup=amin1(bmax,msn2+m)
       binf=amax1(bmin,2.*m)
       stm=(bsup-binf)/(nw-1)
       if (stm.le.0.0) return
c   integrale
       do 1 i=1,nw
        bm=binf+(i-1)*stm
        fi2=fi2+w(i)*effe(m/bm)*fiteorica(bm)*m/(bm*bm)
   1   end do
       fi2=alf*stm*fi2
       return
       end
c --------------------------------------------------------
c funcion fir remanente de primarias en función de la masa de secundarias
        function fir(m)
         real m,mmax,mmin,msep,msn2
       common/masse/mmax,mmin,msep,msn2,bmin,bmax,rm
       common/sndat/alf,umalf,gamma,costfmu,cs
       common/intgr/nw,w
       dimension w(10), el(19)
       fir = 0
       bsup=amin1(bmax, msn2+m)
       binf=amax1(bmin,2.*m)
       stm=(bsup-binf)/(nw-1)
       if (stm.le.0.0) return
c   integrale
       do 1 i=1,nw
       bm=binf+(i-1)*stm
c contribución de las remanentes primarias con el(1)=d(bm-m)
        if (m.lt.8.) then
        resta = 0
        end if
        resta = bm-m
        call interp (resta,el)
        fir=fir+w(i)*effe(m/bm)*fiteorica(bm)*el(16)*(bm-m)/(bm*bm)
1           continue
        fir = alf*stm*fir
       return
       end     
c --------------------------------------------------------
       function effe(rmu)
       common/sndat/alf,umalf,gamma,costfmu,cs
       effe=costfmu*rmu**gamma
       return
       end
c -------------------------------------------
       function tau(emme)
       if (emme.gt.10.) goto 10
       if (emme.le.0.1) emme=0.1
       tau=8.*emme**(-2.8)
       return
10     tau=0.05*emme**(-0.6)
       return
       end
c -------------------------------------------
       function emme(tau)
       if (tau.le.1.e-6) tau=1.e-6
       if (tau.lt.0.0127) goto 10
       emme=(tau/8.)**(-1/2.8)
       return
10     emme=(tau/0.05)**(-1/0.6)
       return
       end
c -------------------------------------------
       function ennea(t)
       real t
       ennea=0.
       if (t.eq.0.) return
       a=alog10(t)
       b=-1.4
       if(a.gt.b)  ennea=0.003252*(a-b)
       return
       end
c --------------------------------------------
       function enneb(t)
       real t
       b=-1.2
       c=-0.1
       enneb=0.
       if (t.eq.0.) return
       a=alog10(t)
       a=amin1(a,c)
       if(a.gt.b) enneb=0.02497*(a-b)
       return
       end
c----------------------------------------------------------------------
       function etout(t)
       real t,tc,rt
       etout=0.
       if (t.eq.0.)return
       tc=5.3e-5
       rt=(tc/t)**0.4
       if (t.gt.tc) goto 21
       etout=8.67e3*t
       return
21     etout=(1-0.44*(rt**2)*(1-0.41*rt)-0.22*(rt**2))
       return
       end
c ---------------------------------------------
       block data
       real mmax,mmin,msep,msn2
       common/const/nmas,imax,irid,imax1,jmax,ic
       common/masse/mmax,mmin,msep,msn2,bmin,bmax,rm
       common/sndat/alf,umalf,gamma,costfmu,cs
       common/intrp/beta,erre
       common/abund/x1,x2,x3,x4,x12,x16,x14,x13,x,feh
       common/intgr/nw,w
       dimension w(10)
       dimension x(6)
       data nmas,imax,irid,jmax,ic/80,15,10,9,0/
       data mmax,mmin,rm/100.,0.8,1.4/
       data beta,erre/0.6,8./
       data x1,x2,x3,x4,x12,x16,x14,x13,x/
     1 .705,.48e-4,.29e-4,.275,.28e-2,.76e-2,.082e-2,
     1 .37e-4,.17e-2,.64e-3,.70e-3,.48e-3,.65e-4,.12e-2/
       data nw,w/7,41.,216.,27.,272.,27.,216.,41.,3*0./
c      data cs/1/
       data cs/0.3/	
       end


