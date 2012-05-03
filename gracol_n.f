c ----------------------------------------------------------------
c               f.ferrini - c.pardi - u.penco (pisa)
c        programma di evoluzione chimica di un modello galattico
c                   versione 3.0 - 27.11.90
c
c        routine di integrazione: lsode
c        le equazioni dei metalli tengono conto dei ritardi nella
c        produzione e utilizzano la matrice q(i,j) di matteucci
c        e comprendono il contributo di supernovae di tipo i.
c        si calcolano i rates numerici di snei e sneii
c -----------------------------------------------------------------
c        e' necessario predisporre i files da leggere:
c     pargra :
c         iw      flag di scrittura video
c         ic      flag di matrice q completa
c         lm2     dim. memoria circolare stelle massicce     max=350
c         lblk    num. passi per mem. stelle piccole
c         msep    valore separat. tra stelle massicce e piccola m
c         msn2    limite inferiore snii
c         bmin    limite inf. binarie per sni
c         bmax    limite sup. binarie      "
c         alf     frazione di binarie
c         gamma   param. distr. binarie
c         ttot    tempo totale integrazione (gyr)
c         nsez      sezioni di output
c         t1, dt1    tempo e intervallo sezione 1
c         t2, dt2                idem   sezione 2 ...
c         rtol  tolleranza relativa
c         atol  tolleranza assoluta
c      ----------------------------------------------------------
c         rgal   raggio totale galassia (in kpc)
c         rring  raggio interno anjllo
c         delr   larghezza anjllo 
c         tcoll  tempo collasso alone (in 10e7 yr) reciproco ff
c         epsk   efficienza per k1+k2
c         epsmu       "         mu
c         epsh        "     per h1+h2+hh
c         eps1h       "     per per h1+h2 (da h)
c         epsa        "     per a1+a2+aa
c         eps1a       "     per a1+a2 (da a)   
c        -------------------------------- una volta ......
c         sa12   coeff. formaz. per inter. nubi-stelle di
c         aa            distr.     "           "
c         sh12          formaz. stelle per coll. nube-nube di
c         hh            distr.  nube-nube        di
c         sak12         formaz. stelle framm.sp. ha
c         en     esponente legge framm. spont.   ha e di
c         um     coeff. formaz. spont. nubi      di
c         ff     coeff. di infall gas ha su di
c  --------------------------------------------------------------
c     fisnr :  fi1(k) + fi2(k)
c                media di imf su interv.di massa k.mo
c                per (fi1) stelle singole, primarie di binarie e
c                (fi2) secondarie di binarie
c               fisii(k), fisi(k)
c                frazione numerica dei sistemi che danno
c                luogo a snii e sni.
c     qmsnr :  qfm(k,i,j) fraz.di elem. j che diventa i
c                ed espulsa mediata su intervallo k.mo
c                con pesi fi(m) per stelle ad evol. normale
c                e fir(m) per remn. di primarie di binarie
c     incipit : valori iniziali delle varie fasi e
c                delle abbondanze
c -------------------------------------------------------------
       external f,jac
       character*6 num
       double precision atol,rtol,rwork,t,tout,delt,y
       double precision a1,a2,aa,h1,h2,hh,ak1,ak2,en,um,ff
       double precision sum,sa12,sh12,sak12
       double precision rgal, rring, delr, tcoll, epsk
       double precision epsmu, epsh, eps1h, epsa, eps1a
       double precision htot, atot, valo,vdisc,denom
       double precision fnq,qfm,q1,qma,umd,totstar
       double precision fi,fik,fisii,fisiik,fisi,fisik
       double precision stnad
       dimension y(42),rwork(4600),iwork(90)
       dimension tmp(500),dtmp(500),npas(500),nrec(500)
       dimension fi(800),fisii(800),fisi(800),qfm(800,15,15),q1(15)
       dimension il(800)
       character*80 msf,disk,dione,dihc,dimgfe,cstars,
     1 rasne,halo,hahc,haone,hamgfe,report
c ------- halo -----
       double precision hamas,psiha,psiham,psihp,psihpm,tsh,tshm,tsh1
       double precision xhp,xhpm,remh,wh,w1h,w2h,proh,proha
       double precision prohl,prohl0,prohl1,prohal,prohal1
       double precision rh,rh0,rh1,r1h,r2h,r1h0,r1h1
       double precision rsnhii,rsnhi
       dimension psihp(350),psihpm(350),xhp(15,350),xhpm(15,350)
       dimension proh(15),prohl(15),prohl0(15),prohl1(15)
       common/dati1h/r1h,r1h0,r1h1,r2h,w1h,w2h,wh,tsh,tshm,remh,
     1           rsnhi,rsnhii
       common/dati2h/proh
c ------- disk -----
       double precision dimas,psidi,psidim,psidp,psidpm,tsd,tsdm,tsd1
       double precision xdp,xdpm,remd,wd,w1d,w2d,prod,prodi
       double precision prodl,prodl0,prodl1,prodil,prodil1
       double precision rd,rd0,rd1,r1d,r2d,r1d0,r1d1
       double precision rsndii,rsndi
       dimension psidp(350),psidpm(350),xdp(15,350),xdpm(15,350)
       dimension prod(15),prodl(15),prodl0(15),prodl1(15)
       common/dati1d/r1d,r1d0,r1d1,r2d,w1d,w2d,wd,tsd,tsdm,remd,
     1           rsndi,rsndii
       common/dati2d/prod,iod
c
       common/plsode/itol,itask,istate,iopt,lrw,liw,mf
       common/pmod/nmt,nmtr,nmtp,imas
       common/dati/a1,a2,aa,h1,h2,hh,ak1,ak2,en,um,ff
c --------------------------------------------------------------
       data neq /42/
       data itol,itask,istate,iopt,lrw,liw,mf/1,4,1,0,4600,90,21/
       data nmt,nmtr,nmtp,imas/15,10,8,0/
       data r1h,r1h0,r1h1,r2h,w1h,w2h,wh,tsh,tshm,remh,rsnhi,
     1           rsnhii/12*0.d0/
       data r1d,r1d0,r1d1,r2d,w1d,w2d,wd,tsd,tsdm,remd,rsndi,
     1           rsndii/12*0.d0/
c --------------------------------------------------------------
c              variabili
c
c     offset per    alone: ioh=0  - disco: iod=6+nmt
c         y(off+1) = gas           y(off+2) = nubi
c         y(off+3) = st.massicce   y(off+4) = altre stelle
c         y(off+5) = remnants      y(off+6) = tot. st. form.
c      da y(off+7) a y(off+nmt+6)=abbondanze relative di elementi
c
c        psihp,psidp,xhp,xdp valori di sfr media e abbondanze
c           conservate nella memoria circolare (2)
c        psihpm,psidpm,xhpm,xdpm valori di sfr media e abbondanze
c           conservate nella memoria (1)
c
c            inizializzazioni
c
c ---------------------------------------------------------------
        open(unit=90,file='modelo')
       read(90,*)num
        close(90)
            report='report_'//num
      open (25,file=report)
c
c ---------------------------------------------------------------
c     lettura dei parametri
c ---------------------------------------------------------------
       open (90,file='pargra',status='old')
       open (91,file='qmsnr',status='old')
       open (92,file='fisnr',status='old')
c ----------------------------------------------
2001   format  (3x,'lettura dei parametri e coefficienti...')
       write  (25,2001)
2003   format (3i5/e15.8,f7.2)
       read (92,2003) lm2f,lblkf,lm1,delt,ttot
2005   format (i6)
2020   format (i3)
2006   format (f6.1)
2021   format (1x,2f12.6)
       read (90,2005) iw,ic,lm2,lblk
          if (iw.eq.1) write (6,2001)
       read (90,2006) emsep,emsn2,bmin,bmax
       read (90,2006) alf,gamma,ttotn
c       print *,iw,ic,lm2,lblk
c       print *,alf,gamma,ttotn
       open (unit=1,file='tiemposcal.dat',status='old')
       read (1,2020) nsez
          tt=0.0
       do 210 i=1,nsez
2007   format (2f6.1)
        read (1,2021) tmp(i),dtmp(i)
           npas(i)=int(0.5+(dtmp(i)/delt))
                    
          if(npas(i).eq.0)npas(i)=1
          nrec(i)=int(0.5+(tmp(i)/(npas(i)*delt)))
          if(nrec(i).eq.0)nrec(i)=1
          tt=tt+delt*nrec(i)*npas(i)
         print *,tmp(i),dtmp(i),npas(i),nrec(i),tt
210    end do
      close(1)
213   format (d14.4)
       read (90,213) rtol,atol
c    ---------------------------
       read (90,213) rgal,rring,rsol,delr,tcollsol,esc
       read (90,213) epsk,epsmu,epsh,eps1h,epsa,eps1a
       close(90)
c    -------------------------------------

       tcolzero=tcollsol*exp(-rsol/esc)
       tcoll=tcolzero*exp(rring/esc)
       denom=rgal*rring*delr
       valo=denom*(1.-(rring/rgal)**2.)**0.5
       vdisc=denom*alt(rring)/rgal
         sak12=0.19*epsk/(valo**0.5)
         um=0.19*epsmu/(vdisc**0.5)
         ff=1./tcoll
c         htot=2.97*eps1h/(vdisc**0.5)
         htot = 2.97*eps1h/vdisc
         sh12=epsh*htot
         hh=(1.-epsh)*htot
         atot= eps1a*0.33d8
         sa12=epsa*atot
         aa=(1.-epsa)*atot
         en=1.5
c ---------------------------------------------------------------
         write (25,3001) valo, vdisc
 3001    format (3x,'valo, vdisc = ',2d14.4)
         write (25,3002) sak12
 3002    format (3x,'k1+k2 = ',d14.4)
         write (25,3003) um
 3003    format (3x,'mu = ',d14.4)
         write (25,3004) ff
 3004    format (3x,'ff = ',d14.4)
         write (25,3005) sh12,hh
 3005    format (3x,'h1+h2, hh = ',2d14.4)
         write (25,3006) sa12,aa
 3006    format (3x,'a1+a2, aa = ',2d14.4)
c  ---------------------------------------------------
c   controlli e calcolo di variabili ausiliarie
       if (lm1.gt.200) goto 805
       if (lm2f.ne.lm2.or.lblkf.ne.lblk) goto 815
       if (tt.gt.ttot) goto 825
c
       if (ic.eq.0) nmt=nmtr
       iod=nmt+6
       iod6=iod+6
       nmq=lm1+lm2
       imas=nmq
       khi=lm2
c ---------------------------------------------------------------
c          conversione in unita' interne di tempo (1.e7 yr)
       delt=100.*delt
c ---------------------------------------------------------------
c          azzeramento delle memorie  (per prudenza...)
101    format (3x,'azzeramento memorie...')
       write (25,101)
        if (iw.eq.1) write (6,101)
       do 110 j=1,lm2
        psihp(j)=0.d0
        psidp(j)=0.d0
        do 110 i=1,nmt
          xhp(i,j)=0.d0
          xdp(i,j)=0.d0
110    end do
       do 120 j=1,lm1
        psihpm(j)=0.d0
        psidpm(j)=0.d0
        do 120 i= 1,nmt
          xhpm(i,j)=0.d0
          xdpm(i,j)=0.d0
120       end do
c ---------------------------------------------------------------
c       lettura delle matrici qfm(k,i,j) e delle fix(k)
       do 220 k=1,nmq
        do 220 i=1,nmt
          do 220 j=nmtp+2,nmt
            qfm(k,i,j)=0.d0
220     end do
       do 230 k=1,imas
        fi(k)=0.d0
        fisi(k)=0.d0
        fisii(k)=0.d0
230    end do
       open(unit=1,file='nest.dat')
       read(1,*)totstar
c      fnq=1.d-6
       fnq=1/totstar
       close(1)
       print*,totstar,fnq
c ----------------------------------------------
202   format (3x,'lettura delle qfm(k,i,j)...')
      write (25,202)
        if (iw.eq.1) write (6,202)
       do 240 k=1,nmq
        do 241 i=1,nmt
          read(91,290) (q1(j),j=1,nmtp+1)
          do 241 j=1,nmtp+1
            qfm(k,i,j)=fnq*q1(j)
241     end do
        umd=qfm(k,nmtp+1,nmtp+1)
        do 240 j=nmtp+2,nmt
          qfm(k,j,j)=umd
240     end do
c ----------------------------------------------
203    format (3x,'lettura delle fix(k)...')
       write (25,203)
        if (iw.eq.1) write (6,203)
       sum=0.d0
       ipnz=0
       do 250 k=1,imas
        read(92,291)fik,fisiik,fisik
        fi(k)=fnq*fik
        if (k.le.khi) sum=sum+fi(k)
        fisii(k)=fnq*fisiik
        if (fisik.ne.0.0.or.ipnz.ne.0) goto 245
          ksn=k
          goto 250
245     continue
          fisi(k)=fnq*fisik
          ipnz=1
250     end do

 291   format(4f11.3,f11.3)
 290   format(9(E12.6,1x))
       close(91)
       close(92)
c -------- determinazione del rapporto dei coefficienti --------
204    format (3x,'somma delle masse grandi...',f10.6)
       write(25,204) sum
        if (iw.eq.1) write(6,204) sum
       a2=sa12*sum
       a1=sa12-a2
       ak2=sak12*sum
       ak1=sak12-ak2
       h2=sh12*sum
       h1=sh12-h2
c ---------------------------------------------------------------
c         lettura delle condizioni iniziali
c    ------------------------------------------------------------
       open(93,file='incipit',status='old')
301    format (3x,'lettura delle condizioni iniziali...')
       write (25,301)
        if (iw.eq.1) write (6,301)
       read(93,213) (y(i),i=1,neq)
       close (93)
c      y(1)=1.d0-y(iod+1)
c --------------------------------
       do 310 i=1,nmt
        ip6=i+6
        xhp(i,1)=y(ip6)
        xhpm(i,1)=y(ip6)
        ip6=ip6+iod
        xdp(i,1)=y(ip6)
        xdpm(i,1)=y(ip6)
310    end do
c ---------------------------------------------------------------
c        apertura delle files

            msf='msf_'//num
            disk='disk_'//num
            dihc='dihc_'//num
            dione='dione_'//num
            dimgfe='dimgfe_'//num
            halo='halo_'//num
            hahc='hahc_'//num
            haone='haone_'//num
            hamgfe='hamgfe_'//num
            rasne='rasne_'//num
            cstars='cstars_'//num


        open (10,file=msf    )
        open (11,file=halo   )
        open (12,file=hahc   )
        open (13,file=haone  )
        open (14,file=hamgfe )
        open (15,file=disk   )
        open (16,file=dihc   )
        open (17,file=dione  )
        open (18,file=dimgfe )
        open (19,file=cstars )
        open (20,file=rasne  )

c
 
c --------- inizializzazione dell'integrazione  -------------
        t = 0.d0
        tout = delt
        rwork(1)=tout
        it  = 1
c ---------------------------------------------------------------
c           scrittura nei files dei valori iniziali
        hamas=y(1)
        dimas=y(iod+1)
        t100=t/100
        psiha3=0.d0
        psidi3=0.d0
        sum=1.d0
          if (iw.eq.1) write (6,910) t100,hamas,dimas,
     1 sum,psiha3,psidi3
        write (10,910) t100, hamas,dimas, sum, psiha3, psidi3
        write (11,911) y(1),y(2),y(3),y(4),y(5)
        write (12,912) (y(6+i),i=1,5)
        write (13,912) (y(6+5+i),i=1,5)
        write (14,912) (y(6+10+i),i=1,5)
        write (15,911) y(iod+1),y(iod+2),y(iod+3),y(iod+4),y(iod+5)
        write (16,912) (y(iod6+i),i=1,5)
        write (17,912) (y(iod6+5+i),i=1,5)
        write (18,912) (y(iod6+10+i),i=1,5)
        write (19,911) y(6), y(iod+6)
        write (20,912) rsnhii,rsndii,rsnhi,rsndi
910     format (1x,f10.6,5(1x,f13.9))
911     format (1x,5f15.11)
912     format (1x,5d14.6)
c ---------------------------------------------------------------
c
c          inizio loop di integrazione
c
c ---------------------------------------------------------------
501    format(1x,'inizio della integrazione ...')
       write(25,501)
        if (iw.eq.1) write(6,501)
c
        do 510 isez = 1,nsez
         do 510 iout = 1,nrec(isez)
          do 520 jout = 1,npas(isez)
         call lsode (f,neq,y,t,tout,itol,rtol,atol,itask,istate,
     1                  iopt,rwork,lrw,iwork,liw,jac,mf)
            if (istate.lt.0) go to 845
c ------- calcolo di indici -------------------------------
        it1=it/lblk
        ir=mod(it,lblk)
        it2=1+mod(it-1,lm2)
c ------- verifica della conservazione della massa --------
        hamas = y(1)+y(2)+y(3)+y(4)+y(5)
        dimas = y(iod+1)+y(iod+2)+y(iod+3)+y(iod+4)+y(iod+5)
        sum = hamas+ dimas
c       if(dabs(mtot-sum).le.1.d-4) goto 21
c22     format (//' warning... massa non conservata...:  sum = ',d15.6)
c        write (25,22) sum
c          if (iw.eq.1) write (6,22) sum
c21     hamas=hamas/sum
c       dimas=dimas/sum
c       do 521 i=1,5
c        y(i)=y(i)/sum
c        y(iod+i)=y(iod+i)/sum
c521    continue
c ----- calcolo della s.f.r. -----------------------------
       tsh1=y(6)
       psiha=(tsh1-tsh)/delt
       tsh=tsh1
       tsd1=y(iod+6)
       psidi=(tsd1-tsd)/delt
       tsd=tsd1
       if (ir.ne.0) goto 531
        psiham=(tsh1-tshm)/(lblk*delt)
        tshm=tsh1
        psidim=(tsd1-tsdm)/(lblk*delt)
        tsdm=tsd1
531    continue
c ------ trasferimento nella mem. circolare dei valori calcolati
       psihp(it2)=psiha
       psidp(it2)=psidi
       do 540 i=1,nmt
        xhp(i,it2)=y(6+i)
        xdp(i,it2)=y(iod6+i)
540    end do
c ------ trasferim. nella memoria masse piccole ----------------
       if (ir.ne.0) goto 551
        psihpm(it1)=psiham
        psidpm(it1)=psidim
        do 550 i=1,nmt
          xhpm(i,it1)=y(6+i)
          xdpm(i,it1)=y(iod6+i)
550      end do
551    continue
c ---------------------------------------------------------------
c     calcolo dei termini con ritardo da utilizzare in f
c     r1h,r2h rate di massa che esce dalla fase stelle 1,2
c     r1d,r2d rate di massa che esce dalla fase stelle 1,2
c     prohl(i),proh(i) rate di prod. di elem. i in ha da stelle 1+2
c     prodl(i),prod(i) rate di prod. di elem. i in di da stelle 1+2
c     w1h,w2h rate di restit.totale al gas in ha da stelle 1,2
c     w1d,w2d rate di restit.totale al gas in di da stelle 1,2
c ---------------------------------------------------------------
c        contributo stelle massicce (2)
       r2h=0.d0
       r2d=0.d0
       kmax=min(khi,it)
       do 560 k=1,kmax
        ilk=it2+1-k
        if (ilk.le.0) ilk=lm2+ilk
        il(k)=ilk
        fik=fi(k)
        r2h=r2h+psihp(ilk)*fik
        r2d=r2d+psidp(ilk)*fik
560    end do
       w2h=0.d0
       w2d=0.d0
       do 570 i=1,nmt
        proha=0.d0
        prodi=0.d0
        do 571 k=1,kmax
          rh=0.d0
          rd=0.d0
          ilk=il(k)
          do 572 j=1,nmt
            qma=qfm(k,i,j)
            if(qma.eq.0.0) goto 572
              rh=rh+xhp(j,ilk)*qma
              rd=rd+xdp(j,ilk)*qma
572       end do
          proha=proha+psihp(ilk)*rh
          prodi=prodi+psidp(ilk)*rd
571     end do
        w2h=w2h+proha
        w2d=w2d+prodi
        proh(i)=proha
        prod(i)=prodi
570    end do
c ---------------------------------------------------------------
c          contributo stelle piccole (1)
c            calcolo dei termini che variano ogni lblk passi
c     if (ir.ne.0) goto 581
        r1h0=r1h1
        r1h1=0.d0
        r1d0=r1d1
        r1d1=0.d0
        kmax=min(it1+1,imas-khi)
        do 580 k=2,kmax
          fik=fi(khi+k)/lblk
          ips=it1+2-k
          r1h1=r1h1+psihpm(ips)*fik
          r1d1=r1d1+psidpm(ips)*fik
580     end do
        do 590 i=1,nmt
          prohal1=0.d0
          prodil1=0.d0
          kmax=min(it1+1,nmq-khi)
          do 591 k=2,kmax
            rh1=0.d0
            rd1=0.d0
            ips=it1+2-k
            do 592 j=1,nmt
              qma=qfm(khi+k,i,j)
              if(qma.eq.0.0) goto 592
                qma=qma/lblk
                rh1=rh1+xhpm(j,ips)*qma
                rd1=rd1+xdpm(j,ips)*qma
592         end do
            prohal1=prohal1+psihpm(ips)*rh1
            prodil1=prodil1+psidpm(ips)*rd1
591      end do
          prohl0(i)=prohl1(i)
          prohl1(i)=prohal1
          prodl0(i)=prodl1(i)
          prodl1(i)=prodil1
590     end do
581     continue
c   ------ calcolo dei termini ad ogni passo ---------
       r1h=(lblk-ir)*r1h0+ir*r1h1
       r1d=(lblk-ir)*r1d0+ir*r1d1
       w1h=0.d0
       w1d=0.d0
       do 595 i=1,nmt
        prohal=(lblk-ir)*prohl0(i)+ir*prohl1(i)
        prodil=(lblk-ir)*prodl0(i)+ir*prodl1(i)
        prohl(i)=prohal
        prodl(i)=prodil
        w1h=w1h+prohal
        w1d=w1d+prodil
595    end do
c --------------------------------------------------------------
c        si sommano le produzioni di elementi da stelle 1 e 2
       do 596 i=1,nmt
        proh(i)=prohl(i)+proh(i)
        prod(i)=prodl(i)+prod(i)
596    end do
       wh=w1h+w2h
       wd=w1d+w2d
c ---------------------------------------------------------------
c        incremento dei contatori
       it=it+1
       tout=tout+delt
       rwork(1)=tout
c
520    end do
c         fine loop interno di integrazione
c ---------------------------------------------------------------
c         calcolo dei rates di supernovae
c ---------------------------------------------------------------
c     supernovae ii :
c        contributo da stelle massicce (2)
       rsnhii=0.d0
       rsndii=0.d0
       kmax=min(ksn,it)
       do 610 k=1,kmax
        ilk=it2+1-k
        if (ilk.le.0) ilk=lm2+ilk
        fisiik=fisii(k)
        rsnhii=rsnhii+psihp(ilk)*fisiik
        rsndii=rsndii+psidp(ilk)*fisiik
610    end do
c ---------------------------------------------------------------
c     supernovae i :
c        contributo da stelle tra 8 e 4 (massicce)
       rsnhi=0.d0
       rsndi=0.d0
       if (it.le.ksn) goto 690
       kmax=min(khi,it)
       do 620 k=ksn,kmax
        ilk=it2+1-k
        if (ilk.le.0) ilk=lm2+ilk
        fisik=fisi(k)
        rsnhi=rsnhi+psihp(ilk)*fisik
        rsndi=rsndi+psidp(ilk)*fisik
620     end do
c -------------------------------------------
c        contributo da stelle tra  4 e mmin (piccole)
        if (it.lt.khi) goto 690
        rh0=0.d0
        rh1=0.d0
        rd0=0.d0
        rd1=0.d0
        kmax=min(it1+1,imas-khi)
        do 630 k=2,kmax
        fisik=fisi(khi+k)/lblk
        ips0=it1+1-k
        ips1=ips0+1
        rh0=rh0+psihpm(ips0)*fisik
        rd0=rd0+psidpm(ips0)*fisik
        rh1=rh1+psihpm(ips1)*fisik
        rd1=rd1+psidpm(ips1)*fisik
630     end do
        if (khi+kmax.eq.imas) goto 631
        fisik=fisi(khi+kmax+1)/lblk
        rh1=rh1+psihpm(ips0)*fisik
        rd1=rd1+psidpm(ips0)*fisik
631     rsnhi=rsnhi+(lblk-ir)*rh0+ir*rh1
        rsndi=rsndi+(lblk-ir)*rd0+ir*rd1
690     continue
c ---------------------------------------------------------------
c     scrittura dei risultati
c ---------------------------------------------------------------
        t100=t/100
        psiha3=1000*psiha
        psidi3=1000*psidi
          if(iw.eq.1) write(6,910)t100,hamas,dimas,sum,psiha3,psidi3
        write (10,910) t100, hamas,dimas,sum, psiha3, psidi3
        write (11,911) y(1), y(2), y(3), y(4), y(5)
        write (12,912) (y(6+i),i=1,5)
        write (13,912) (y(6+5+i),i=1,5)
        write (14,912) (y(6+10+i),i=1,5)
        write (15,911) y(iod+1),y(iod+2),y(iod+3),y(iod+4),y(iod+5)
        write (16,912) (y(iod6+i),i=1,5)
        write (17,912) (y(iod6+5+i),i=1,5)
        write (18,912) (y(iod6+10+i),i=1,5)
        write (19,911) y(6), y(iod+6)
        write (20,912) rsnhii,rsndii,rsnhi,rsndi
510     end do
c          fine loop di scrittura risultati e di sezioni
c ---------------------------------------------------------------
c          chiusura delle files e del programma
c ---------------------------------------------------------------
9995      continue
        close (10)
        close (11)
        close (12)
        close (13)
        close (14)
        close (15)
        close (16)
        close (17)
        close (18)
        close (19)
        close (20)
        close (25)
      stop
c ---------------------------------------------------------------
c         situazioni di errore
c ---------------------------------------------------------------
800    format (//' memoria insufficiente per stelle p.m. ...')
805    write (25,800)
        if (iw.eq.1) write (6,800)
       goto 9991
810    format (//' parametri incoerenti con files esterne...')
815    write (25,810)
         if (iw.eq.1) write (6,810)
       goto 9991
820    format (//' tempo richiesto maggiore del possibile...')
825    write (25,820)
        if (iw.eq.1) write (6,820)
9991   close (92)
       close (93)
       close (25)
       stop
840    format (//' programma interrotto in lsode... istate = ',i3)
845    write (25,840) istate
        if (iw.eq.1) write (6,840) istate
       goto 9995
       end
c ---------------------------------------------------------------
c    altezza del disco galattico
c -----------------------------------------------------
       function alt(r)
       double precision r
       alt = 0.2
       return
       end
c ---------------------------------------------------------------
c     sistema di equazioni del modello
c ---------------------------------------------------------------
       subroutine  f (neq, t, y, ydot)
       double precision t, y, ydot
       dimension y(42), ydot(42)
       dimension proh(15),prod(15)
       double precision a1,a2,aa,h1,h2,hh,ak1,ak2,en,um,ff
       double precision r1h,r1h0,r1h1,r2h,w1h,w2h,wh
       double precision tsh,tshm,remh,rsnhi,rsnhii,proh
       double precision r1d,r1d0,r1d1,r2d,w1d,w2d,wd
       double precision tsd,tsdm,remd,rsndi,rsndii,prod
       double precision hgen,a1hgen,a2hgen,ffhg,dcc,dcs,udgen,dgc
       double precision a1dcs,a2dcs,aadcs,h1dcc,h2dcc,hhdcc
       common/pmod/nmt,nmtr,nmtp,imas
       common/dati/a1,a2,aa,h1,h2,hh,ak1,ak2,en,um,ff
       common/dati1h/r1h,r1h0,r1h1,r2h,w1h,w2h,wh,tsh,tshm,remh,
     1           rsnhi,rsnhii
       common/dati2h/proh
       common/dati1d/r1d,r1d0,r1d1,r2d,w1d,w2d,wd,tsd,tsdm,remd,
     1           rsndi,rsndii
       common/dati2d/prod,iod
       hg=y(1)
       dg=y(iod+1)
       dc=y(iod+2)
       ds=y(iod+3)
       hgen=hg**en
       a1hgen=ak1*hgen
       a2hgen=ak2*hgen
       ffhg=ff*hg
       dcc=dc*dc
       dcs=dc*ds
       udgen=um*dg**en
       a1dcs=a1*dcs
       a2dcs=a2*dcs
       aadcs=aa*dcs
       h1dcc=h1*dcc
       h2dcc=h2*dcc
       hhdcc=hh*dcc
       dgc=dg+dc
       ydot(1) = wh-(a1hgen+a2hgen)-ffhg
       ydot(2) = 0.d0
       ydot(3) = a2hgen-r2h
       ydot(4) = a1hgen-r1h
       ydot(5) = r1h+r2h-wh
       ydot(6) = a1hgen+a2hgen
       ydot(iod+1) = wd-udgen+ffhg+aadcs+hhdcc
       ydot(iod+2) = udgen-(a1dcs+a2dcs+aadcs)-(h1dcc+h2dcc+hhdcc)
       ydot(iod+3) = h2dcc+a2dcs-r2d
       ydot(iod+4) = h1dcc+a1dcs-r1d
       ydot(iod+5) = r1d+r2d-wd
       ydot(iod+6) = a1dcs+a2dcs+h1dcc+h2dcc
       do 1 i=1,nmt
        im=6+i
        ydot(im)=(-y(im)*wh+proh(i))/hg
        imd=iod+im
        ydot(imd)=(-y(imd)*wd-(y(imd)-y(im))*ffhg+prod(i))/dgc
1      end do
       return
       end
c ---------------------------------------------------------------
c     matrice jacobiana
c ---------------------------------------------------------------
       subroutine jac (neq, t, y, ml, mu, pd, nrpd)
       double precision pd, t, y
       dimension y(42), pd(nrpd,42)
       dimension proh(15),prod(15)
       double precision a1,a2,aa,h1,h2,hh,ak1,ak2,en,um,ff
       double precision r1h,r1h0,r1h1,r2h,w1h,w2h,wh
       double precision tsh,tshm,remh,rsnhi,rsnhii,proh
       double precision r1d,r1d0,r1d1,r2d,w1d,w2d,wd
       double precision tsd,tsdm,remd,rsndi,rsndii,prod
       double precision en1,hgen1,ddc,whhg,hgq,dgc,wdff,dgcq,
     1 ffdgc
       double precision ffhg,exp1,exp2,dify
       common/pmod/nmt,nmtr,nmtp,imas
       common/dati/a1,a2,aa,h1,h2,hh,ak1,ak2,en,um,ff
       common/dati1h/r1h,r1h0,r1h1,r2h,w1h,w2h,wh,tsh,tshm,remh,
     1           rsnhi,rsnhii
       common/dati2h/proh
       common/dati1d/r1d,r1d0,r1d1,r2d,w1d,w2d,wd,tsd,tsdm,remd,
     1           rsndi,rsndii
       common/dati2d/prod,iod
       en1=en-1.d0
       hg=y(1)
       dg=y(iod+1)
       dc=y(iod+2)
       ds=y(iod+3)
       hgen1=hg**en1
       ddc=2.d0*dc
       whhg=wh/hg
       hgq=hg*hg
       dgc=dg+dc
       wdff=wd*ff
       dgcq=dgc*dgc
       ffdgc=ff/dgc
       ffhg=ff*hg
       exp1=ffhg/dgc
       exp2=-(ffhg+wd)/dgc
c
       pd(3,1) = en*ak2*hgen1
       pd(4,1) = en*ak1*hgen1
       pd(6,1) = pd(3,1)+pd(4,1)
       pd(1,1) = -ff-pd(6,1)
       pd(iod+1,1) = ff
       pd(iod+1,iod+1) = -en*um*dg**en1
       pd(iod+2,iod+1) = -pd(iod+1,iod+1)
       pd(iod+1,iod+2) = hh*ddc+aa*ds
       pd(iod+3,iod+2) = a2*ds+h2*ddc
       pd(iod+4,iod+2) = a1*ds+h1*ddc
       pd(iod+6,iod+2) = pd(iod+3,iod+2)+pd(iod+4,iod+2)
       pd(iod+2,iod+2) = -pd(iod+1,iod+2)-pd(iod+6,iod+2)
       pd(iod+1,iod+3) = aa*dc
       pd(iod+3,iod+3) = a2*dc
       pd(iod+4,iod+3) = a1*dc
       pd(iod+6,iod+3) = pd(iod+3,iod+3)+pd(iod+4,iod+3)
       pd(iod+2,iod+3) = -pd(iod+1,iod+3)-pd(iod+6,iod+3)
       do 1 i=1,nmt
        im=6+i
        pd(im,1) =(wh*y(im)-proh(i))/hgq
        pd(im,im) = -whhg
        imd=iod+im
        dify=y(im)-y(imd)
        pd(imd,1) = dify*ffdgc
        pd(imd,iod+1) = (y(imd)*wd-dify*ffhg-prod(i))/dgcq
        pd(imd,iod+2) = pd(imd,iod+1)
        pd(imd,im) = exp1
        pd(imd,imd) = exp2
1      end do
       return
       end
c ---------------------------------------------------------------
       double precision function dfloat(i)
c
       integer i
       dfloat = i
       return
       end
c ---------------------------------------------------------------
