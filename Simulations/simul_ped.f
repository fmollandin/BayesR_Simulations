      program simulation
      implicit none

      integer(kind=2),dimension(:,:),allocatable::phe
      integer(kind=2),dimension(:,:,:),allocatable::alq
      integer(kind=4),dimension(:),allocatable::chr,all,imc,imt,nmc
      integer(kind=4),dimension(:),allocatable::listm,listmc,listq
      integer(kind=4),dimension(:),allocatable::chrom,pos,listpot

      integer(kind=4),dimension(:),allocatable::nmcum,listqc,ncuq,pere
      integer(kind=4),dimension(:),allocatable::mere
      real(kind=4), dimension(:), allocatable ::vq,effet,y,g
      real(kind=8), dimension(:), allocatable ::freq
      logical, dimension(:), allocatable ::  typ
      integer(kind=4),dimension(:),allocatable:: iani
      
      character*100 fichqtl,fichtypq,ficped,namechr,namechrs,fichperf
      character*100 fichfreq,fichchr,fichchrs
      character*3500000 ligne
      character*128 jour
      logical finfich,tecrit      
      character*1   separateur,optmq,fi,fo
      integer nqg,nqm,nqp,nq,ntyp,nm,nm1,nmax,nt,ntt,np,ipos,nc1,nchrom1
      integer inmt,i,ia,j,jj,k,l,ng,irc,prc,mrc,nchrom,ichr,nqc,im,iq,nc
      integer io
      real vart,h2,pcqtl,ve,vg,rgauss,sdg,sdg1,sdg2,x,xx
      real xmin,xmax,moy,st,minmaf,seed2
      real v1g,v1m,v1p
      intrinsic random_seed,random_number     
      call fdate(jour)
      print *,jour
      print *

      print '(a)','****************************************************'
      print '(a)','***                                  '
      print '(a)','***  Simulation de performances a partir de typages '
      print '(a)','***                v1 - 11/02/2020                  '
      print '(a)','***                                                 '
      print '(a)','****************************************************'
      print *


c *** lecture des parametres ***    
      vart=100.       
      separateur=' '
      
      print *,'Nb de gros, moyens, petits QTL '
      read (5,*) nqg,nqm,nqp
      print *,'Gros, moyens, petits : ',nqg,nqm,nqp

      print *,'h2, part de variance genetique due aux QTL'
      read (5,*) h2, pcqtl
      print *,'h2 et part de variance des qtl : ',h2, pcqtl 

      print *,'Variances relatives des qtl gros moyens petits qtl'
      read (5,*) v1g,v1m,v1p
      print *,'Var relatives des Gros, moyens, petits :',v1g,v1m,v1p

      print *,'Nombre de chromosomes (0 si un seul fichier)'
      read (5,*) nchrom
      print *,'Nombre de chromosomes : ',nchrom  
      nchrom1=nchrom 
      if (nchrom.eq.0) then
        nchrom1=0; nchrom=1
        end if   

      print *,'Fichier pedigree'
      read (5,*) ficped
      print *,'Fichier pedigree : ',trim(ficped)

      print *,'Fichiers marqueurs en entree'
      read (5,*) namechr
      print *,'Fichier marqueur lu  : ',trim(namechr)

      print *,'Fichiers marqueurs en sortie (no si pas de sortie)'
      read (5,*) namechrs
      print *,'Fichier marqueur ecrit: ',trim(namechrs)
    
      print *,'Fichier info sur les QTL'
      read (5,*) fichqtl
      print *,'Fichier QTL : ',trim(fichqtl)

      print *,'Fichier genotypes QTL'
      read (5,*) fichtypq
      print *,'Fichier genot QTL : ',trim(fichtypq)

      print *,'Fichier performances'
      read (5,*) fichperf
      print *,'Fichier Perf : ',trim(fichperf)

      print *,'Fichier de frequences'
      read (5,*) fichfreq
      print *,'Fichier de frequences : ',trim(fichfreq)

      minmaf=0.05
      print *,'Limite de MAF pour les QTL ?'
      read (5,*,iostat=io) minmaf
      print *,'Limite de MAF pour les QTL : ',minmaf
      if (minmaf.lt.0..or.minmaf.gt.0.5) stop 'Maf sup 0.5'
      
      print *,'On garde les QTL dans le fichier marqueur ? (o/n)'
      optmq='o'
      read (5,*,iostat=io) optmq
      if (optmq.ne.'n'.and.optmq.ne.'o') stop 'o ou n'
      if (optmq.eq.'n') print *,'Les qtl sont elimines des marqueurs'
      if (optmq.eq.'o') print *,'Les marqueurs incluent les qtl'

      print *,'Format des fichiers de typages/phases input (P ou T) ?'
      fi='P'
      read (5,*,iostat=io) fi
      if (fi.eq.'p') fi='P'; if (fi.eq.'t') fi='T'; 
      print *,'Format en entree : ',fi
      if (fi.ne.'P'.and.fi.ne.'T') stop 'format different de P et T'

      print *,'Format des fichiers de typages/phases output (P ou T) ?'
      fo='T'
      read (5,*,iostat=io) fo
      if (fo.eq.'p') fo='P'; if (fo.eq.'t') fo='T'; 
      print *,'Format en entree : ',fo
      if (fo.ne.'P'.and.fo.ne.'T') stop 'format different de P et T'

      print *,'Seed ?'
      read (5,*) seed2
      print *,'Seed : ',seed2


c initialisation de la racine de tirage aleatoire
      call init_random_seed(seed2)

c Faut il d'autres formats ?           
      if (namechrs.eq.'no'.and.fi.ne.fo) then
        print *,'Pas de fichiers marqueurs demandes '
        print *,'alors que les formats E/S sont differents ?'
        stop
        end if
      if (namechrs.eq.'no'.and.optmq.eq.'n') then
        print *,'Pas de fichiers marqueurs demandes '
        print *,'alors que les qtl sont supprimes des marqueurs ?'
        stop
        end if

      tecrit=.true.
      if (namechrs.eq.'no') tecrit=.false.
      if (fi.eq.fo.and.namechrs.ne.'no'.and.optmq.eq.'o') then
        print *,'La reecriture des fichiers marqueurs est inutile'
        print *,'Elle n est pas realisee'
        tecrit=.false.
        end if

      
c lecture du fichier pedigree pour verifier qu'il est correct (irc de 1 a n, parent<produit) et compter
      ng=0
      open (2,file=ficped,form='formatted')      
      do 
       read (2,*,iostat=io) irc,prc,mrc
       if (io.ne.0) exit
       ng=ng+1
       if (ng.ne.irc) stop 'les individus non de 1 a n dans le pedigree'
       if (prc.ge.irc) stop 'pere > individu'
       if (mrc.ge.irc) stop 'mere > individu'
       end do
      close (2)
      print *,'Nombre d individus dans le pedigree : ',ng
      
c y=phenotype, g=valeur genetique, typ=indicateur typé ou non
      allocate (y(ng),g(ng),iani(ng),pere(ng),mere(ng),typ(ng))
      y(:)=0. ; typ(:)=.false.

c 2eme lecture du fichier pedigree et stockage
      open (2,file=ficped,form='formatted')      
      do i=1,ng
       read (2,*) irc,prc,mrc
       if (prc.lt.0) prc=0; if (mrc.lt.0) mrc=0
       pere(i)=prc ; mere(i)=mrc
       end do
      close (2)
      

c variance des qtl
      nq=nqg+nqm+nqp
      print *,'Nombre de qtl : ', nq

      if (v1p.le.0.) stop 'variance des petits qtl <= 0.'
c on rapporte les variances de qtl à celle des petits qtl
      v1g=v1g/v1p
      v1m=v1m/v1p

c xx variance expliquee par les qtl, en equivalent petits qtl
      xx=nqg*v1g + nqm*v1m + nqp
      x=vart*h2*pcqtl/xx  ! part des petits qtl
      
c vq=variance de chaque qtl, chr=chrom du qtl, imc=position intra chrom, imt=position globale, effet=effet, all=allele a effet positif
      allocate (vq(nq),chr(nq),imc(nq),imt(nq),effet(nq),all(nq))
      do i=1, nqg; vq(i)=v1g*x;      end do
      do i=1, nqm; vq(nqp+i)=v1m*x;  end do
      do i=1, nqp; vq(nqp+nqm+i)=x; end do


c on calcule le nombre de marqueurs par chrom et les freq alleliques
c nmc = nb marq / chrom, nmcum = nb de marq cumules sur les chrom precedents, ncuq = ?
      allocate (nmc(nchrom),nmcum(nchrom+1),ncuq(nchrom+1))
      do i=1, nchrom; nmc(i)=0; end do
      nm=0 ; ntyp=0 ; nmax=0

c lecture de la premiere ligne des typages pour compter le nb de marqueurs  
      DO ICHR = 1, nchrom
       if (nchrom1.eq.0) then
         fichchr=namechr
       else
         write(fichchr,'(a,i0)') trim(namechr),ichr ; ! cette ligne construit le nom du fichier pour le chrom ichr
         end if
       open(3,file=fichchr,form='formatted')
       call lire_ligne(3,ligne,separateur,nc,finfich)                                        
       if (finfich) stop 'Le fichier de phenotypes est vide'
       close(3)
       call lecligne(ligne,nc)    ; ! cette ligne compte le nombre de champs nc sur la ligne
       if (fi.eq.'P') then
         nc=nc-2
       else
         nc=(nc-1)/2
         end if
       if (nc.gt.nmax) nmax=nc
       nmc(ichr)=nc
       nmcum(ichr)=nm
       nm=nm+nc
       end do
      print *,'Nombre total de marqueurs : ',nm
      print *,'Nombre maxi de marqueurs par chromosome: ',nmax
      nmcum(nchrom+1)=nm

c freq = frequences de l allele 2 par marqueur, phe = typages pour un animal et un chrom
      allocate (freq(nm),phe(2,nmax))
      freq(:)=0.

c exploration des fichiers de typages pour compter le nb d individus, remplir typ, et verifier qu ils sont coherents entre eux
      DO ICHR = 1, nchrom
       if (nchrom1.eq.0) then
         fichchr=namechr
       else
         write(fichchr,'(a,i0)') trim(namechr),ichr ; ! cette ligne construit le nom du fichier pour le chrom ichr
         end if
      open(3,file=fichchr,form='formatted')  
      io=0 ; nt=0 ; ntyp=0
      do while(io.eq.0)
        read(3,*,iostat=io) irc
        if (io.ne.0) exit
        nt=nt+1
        if (irc.gt.ng) stop 'Individu type absent de pedigree'
        typ(irc)=.true. 
        end do
      close (3)
      if (ichr.eq.1) then
        ntt=nt
      else if (nt.ne.ntt) then
        stop 'fichiers de typages non coherents'
        end if
      end do
      print *, 'Nombre de lignes de typages : ',ntt
      ntyp=nt
      if (fi.eq.'P') then
       ntyp=ntyp/2
       print *, 'Nombre d individus types : ',ntyp
       end if


c lecture fichiers de typages pour calculer les frequences
      DO ICHR = 1, nchrom
       if (nchrom1.eq.0) then
         fichchr=namechr
       else
         write(fichchr,'(a,i0)') trim(namechr),ichr ; ! cette ligne construit le nom du fichier pour le chrom ichr
         end if
      open(3,file=fichchr,form='formatted')  
      io=0  ; nc=nmc(ichr) ; nm1=nmcum(ichr)
      do ia=1, ntyp
       if (fi.eq.'P') then   
        read(3,*) irc,i,(phe(1,j),j=1,nc)
        if (i.ne.1) stop 'typages non phases'
        read(3,*) irc,i,(phe(2,j),j=1,nc)
        if (i.ne.2) stop 'typages non phases'
       else
        read(3,*) irc,((phe(i,j),i=1,2),j=1,nc)
        end if
       j=nm1
       do i=1, nc
        j=j + 1
        if (phe(1,i).eq.2) freq(j)=freq(j)+1.
        if (phe(2,i).eq.2) freq(j)=freq(j)+1.
	end do
       if (ichr.eq.1) iani(ia)=irc
       end do	

       close(3)
       END DO
	
c fin de calcul des frequences et sauvegarde dans fichfreq
      do i=1, nm
        freq(i)=freq(i)/dfloat(2*ntyp)
	end do

      open(9,file=fichfreq)
      do i=1, nm
        write(9,'(i6,f8.3)') i,freq(i)
        end do
      close(9)

c on determine les qtl parmi les marqueurs, puis ses caracteristiques
      allocate(listq(nm),listm(nm))
      listq(:)=0 ; listm(:)=0

c on determine les np marqueurs potentiels pour etre qtl, sur la base de la maf
      allocate(listpot(nm),chrom(nm),pos(nm))
      np=0 ; i=0
      do ichr=1, nchrom
       do ipos=1, nmc(ichr) 
        i=i+1
        chrom(i)=ichr; pos(i)=ipos        
        if (freq(i).ge.minmaf.and.freq(i).le.(1.-minmaf)) then  ; 
          np=np+1
          listpot(np)=i  ; ! listpot = liste des np marqueurs potentiels
          end if
        end do
       end do
      if (np.le.nq) stop 'moins de marqueurs potentiels que de qtl'
      print *,'Nombre de marqueurs a MAF non extreme : ',np

! caracteristique des nq qtl : effet (> ou <0), all=allele a effet positif, chr=chromosome, imc=position intr chrom, imt=numero d ordre global du marq, vq=variance
! listq(num_marq)=num_qtl

c on choisit les qtl par tirage dans np possibilites, on stocke l info, en retire le marqueur de la liste, et on passe au suivant
      do i=1, nq
c        j=1+np*rand(x)
        call random_number(x)
        j=1 + np*x
	if (j.lt.1.or.j.gt.np) stop 'erreur dans le tirage des qtl'  
	listq(listpot(j))=i
        imt(i)=listpot(j)
        iq=listpot(j)

        do k=j+1, np
          listpot(k-1)=listpot(k)
          end do
        np=np-1

	chr(i)=chrom(iq)
	imc(i)=pos(iq)

c tirage de l effet (positif)
	effet(i)=.5*sqrt(vq(i)/(2.*freq(iq)*(1.-freq(iq))))
c tirage de l allele a effet positif
	all(i)=1
        call random_number(x)
	if (x.gt.0.5)then
          all(i)=2 
          effet(i)=-effet(i)
          end if

	end do
	
c sauvegarde des caracteristiques des QTL dans fichier fichqtl
      open (4,file=fichqtl,form='formatted')
      do i=1, nq
        write (4,200) i,chr(i),imc(i),imt(i),vq(i),freq(imt(i)),effet(i)
  200 format (I6,i3,2i7,3f12.5)	
	end do
      close(4)	

c mini maxi des effets	
      xmin=effet(1); xmax=xmin
      do i=2, nq
	  if (effet(i).lt.xmin) xmin=effet(i)
	  if (effet(i).gt.xmax) xmax=effet(i)
	  end do
      print *,'Effets mini / maxi : ',xmin, xmax
	

c liste des marqueurs non qtl - sert a preparer les fichiers de marqueurs sans qtl
      if (optmq.eq.'n') then
       j=0
       do i=1, nm
        if (listq(i).eq.0) then
	  j=j+1
	  listm(i)=j ;  ! listm(1:nm) contient 0 pour les qtl, 1:nm-nq pour les marqueurs non qtl ?
	  end if
	end do  
      else
       do i=1, nm
         listm(i)=i
	 end do
       end if

c alq=alleles aux qtl, listqc contient le numero de qtl pour les qtl d un chromosome donné, listmc contient les marqueurs non qtl
      allocate (alq(2,nq,ntyp),listqc(nq),listmc(nm))

c relecture des chromosomes pour traitement des typages   
      DO ICHR = 1, nchrom
       if (nchrom1.eq.0) then
         fichchr=namechr
         fichchrs=namechrs
       else
         write(fichchr,'(a,i0)') trim(namechr),ichr ; ! cette ligne construit le nom du fichier pour le chrom ichr
         write(fichchrs,'(a,i0)') trim(namechrs),ichr ; ! cette ligne construit le nom du fichier pour le chrom ichr
         end if
       print *,'lu    : ', fichchr (1:len_trim(fichchr))	
       print *,'ecrit : ', fichchrs(1:len_trim(fichchrs))

c recherche de la liste des qtl sur le chromosome, stock de leur numéro dans listqc
       nqc=0	
       do iq=1, nq
	  if (chr(iq).eq.ichr) then
	    nqc=nqc+1
	    listqc(nqc)=iq
	    end if
	  end do  

c liste des marqueurs non qtl intra chromosome, pour svg
       nc1=0
       do i=nmcum(ichr)+1, nmcum(ichr+1)
        if (listm(i).ne.0) then
	  nc1=nc1+1
	  listmc(nc1)=i-nmcum(ichr)
	  end if
	end do  
      print *,'Chr',ichr,': ',nmc(ichr),' locus dont ', nqc,' qtl'
	
      open(3,file=fichchr,form='formatted')
      open(4,file=fichchrs,form='formatted')

c lecture proprement dite pour les ntyp animaux
      nc=nmc(ichr) ; nm1=nmcum(ichr)
      do ia=1, ntyp
       if (fi.eq.'P') then   
        read(3,*) irc,i,(phe(1,j),j=1,nc)
        read(3,*) irc,i,(phe(2,j),j=1,nc)
       else
        read(3,*) irc,((phe(i,j),i=1,2),j=1,nc)
        end if

c traitement des qtl : ajout des effets a y
       do im=1, nqc
         iq=listqc(im)	
c effets des qtl sur les performances
	 if (phe(1,imc(iq)).eq.1) then
            y(irc)=y(irc)+effet(iq)
         else
            y(irc)=y(irc)-effet(iq)
            end if
	 if (phe(2,imc(iq)).eq.1) then
            y(irc)=y(irc)+effet(iq)
         else
            y(irc)=y(irc)-effet(iq)
            end if

c stock des alleles qtl, dans l ordre des qtl. On les sauvegardera a la fin, apres trt de tous les chromosomes
         alq(1,iq,ia)=phe(1,imc(iq))    
         alq(2,iq,ia)=phe(2,imc(iq))
	 end do

c reecriture des marqueurs si demande. 
        if (tecrit) then
         if (fo.eq.'P') then
          write (4,100) irc,1,(phe(1,listmc(k)),k=1,nc1)
          write (4,100) irc,2,(phe(2,listmc(k)),k=1,nc1)
        else
          write (4,100) irc,((phe(i,listmc(k)),i=1,2),k=1,nc1)
          end if
        end if

        end do
  100 format (I7,65000i2)	

      close(3)
      close(4)
      
      END DO	
      
c sauvegarde des typages qtl      
      open (4,file=fichtypq,form='formatted')
      do i=1, ntyp
        write (4,100) iani(i),((alq(k,j,i),k=1,2),j=1,nq)
	end do
      close(4)	      
	 
c ajout de la valeur polygenique	 
      vg=h2*vart*(1.-pcqtl)  ! variance poly
      sdg=sqrt(vg)           ! ecart type poly 
      sdg1=sqrt(0.75*vg)     ! ecart type poly si un parent connu
      sdg2=sqrt(0.5*vg)      ! ecart type poly si deux parents connus
      ve=sqrt(vart*(1.-h2))  ! ecart type residuel
      
      moy=0.; st=0.
      	 	 
      do i=1, ng
        call normal(x)
        if (pere(i).eq.0.and.mere(i).eq.0) then
	   g(i)=x*sdg
        else if (mere(i).eq.0) then
	   g(i)=0.5*g(pere(i)) + x*sdg1 
        else if (pere(i).eq.0) then
	   g(i)=0.5*g(mere(i)) + x*sdg1 
        else 
	   g(i)=0.5*(g(pere(i))+g(mere(i))) + x*sdg2
	   end if
	   
        call normal(x)
	y(i)=y(i)+g(i)+ x*ve  
	end do
	
c sauvegarde des performances	
      open (4,file=fichperf,form='formatted')
      moy=0.; st=0.;
      do i=1, ng
        if (typ(i)) then 
	  write (4,*) i,y(i)
          moy=moy+y(i)
	  st=st+y(i)*y(i)
	  end if
	end do
      close(4)	
      if (ntyp.gt.0) moy=moy/ntyp
      if (ntyp.gt.1) then
        st=(st-ntyp*moy*moy)/(ntyp-1)
      else
        st=0.
	end if
      st=sqrt(st)
      print *,'Moyenne des performances : ',moy
      print *,'Ecart type               : ',st	  	

      stop
      end	
 	   

c lecture de ligne ****************************************
      SUBROUTINE lire_ligne(iunit,ligne,separ,nc,fin)                                        
      IMPLICIT none      
      integer nc,i,j,k,l,iunit,ll
      character*3500000 ligne 
      character*1 separ
      logical fin

      fin=.false.
      read(iunit,'(a3500000)',end=5) ligne
      l=len_trim(ligne)
      if (separ.eq.' ') return

      if (ligne(l:l).eq.separ) then
           ligne(l:l)=' '; l=l-1
           end if
      ll=l
      do i=l-1,1,-1
           if (ligne(i:i).eq.separ) then
             if (ligne(i+1:i+1).eq.separ) then
               ligne(i+2:ll+1)=ligne(i+1:ll)
               ligne(i+1:i+1)='_'
               ll=ll+1              
               end if
             end if
           end do
       l=ll
         
       nc=1

       do i=1, l
           if (ligne(i:i).eq.' ') then
c             ligne(i:i)='_'
           else if (ligne(i:i).eq.separ) then
             ligne(i:i)=' '
             nc=nc+1
             end if
           end do

      return 
                                                   
    5 fin=.true.
      return

      end                        
 
c **************************************************************************
      subroutine lecligne(ligne,nc)
      implicit none
      character*3500000 ligne
      integer i,j,k,l,nc
      i=1 
      j=len_trim(ligne) 
      do while(ligne(i:i).eq.' '.and.i.lt.j)
        i=i+1
        end do

      do while(ligne(j:j).eq.' '.and.j.gt.1)
        j=j-1
        end do

      nc=0
      do k=i,j-1
        if (ligne(k:k).ne.' '.and.ligne(k+1:k+1).eq.' ')
     *   nc=nc+1
        end do
       nc=nc+1
       return
       end
       
c **************************************************************************  
      subroutine normal(x)
      implicit none
      real x,y
      integer i,n
      intrinsic random_number
      n=12
      x=0.
      do i=1, n
           call random_number(y)
           x=x+y
           end do
      x=((x/float(n))-0.5)*12.
      return
      end

     
      subroutine init_random_seed(seed1)
!     one integer for seed in parameter file
!     seed1 >0: initialized with seed1
      integer :: a, b, clock
      integer,dimension(:), allocatable :: seed

      call random_seed(size = b)
      allocate(seed(b))
      do a=1,n
         seed(i)=abs(seed1)+(a-1)
      enddo
      call random_seed(put = seed)
      deallocate(seed)
      end subroutine init_random_seed
