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
      print '(a)','***  Simulation of performances from genotypes      '
      print '(a)','***                v1 - 11/02/2020                  '
      print '(a)','***                                                 '
      print '(a)','****************************************************'
      print *


! *** Parameters screaning ***    
      vart=100.       
      separateur=' '
      
      print *,'Nb de big, medium, small QTL '
      read (5,*) nqg,nqm,nqp
      print *,'Big, medium, small : ',nqg,nqm,nqp

      print *,'h2, share of genetic variance due to QTLs'
      read (5,*) h2, pcqtl
      print *,'h2 and QTL share of genetic variance : ',h2, pcqtl 

      print *,'Relatives variances of big, medium and small QTL'
      read (5,*) v1g,v1m,v1p
      print *,'Relatives variance of big, medium, small :',v1g,v1m,v1p

      print *,'Number of chromosomes (0 if one file only)'
      read (5,*) nchrom
      print *,'Number chromosomes : ',nchrom  
      nchrom1=nchrom 
      if (nchrom.eq.0) then
        nchrom1=0; nchrom=1
        end if   

      print *,'Pedigree file'
      read (5,*) ficped
      print *,'Pedigree file : ',trim(ficped)

      print *,'Markers file input'
      read (5,*) namechr
      print *,'Read markers file  : ',trim(namechr)

      print *,'Markers file output (no if no output)'
      read (5,*) namechrs
      print *,'Markers file wrote: ',trim(namechrs)
    
      print *,'QTL info file'
      read (5,*) fichqtl
      print *,'QTL info file : ',trim(fichqtl)

      print *,'QTL genotype files'
      read (5,*) fichtypq
      print *,'QTL genotype files : ',trim(fichtypq)

      print *,'Simulated performances file'
      read (5,*) fichperf
      print *,'Simulated performances file : ',trim(fichperf)

      print *,'Frequency file'
      read (5,*) fichfreq
      print *,'Frequency file : ',trim(fichfreq)

      minmaf=0.05
      print *,'Minimum MAF for the QTLs ?'
      read (5,*,iostat=io) minmaf
      print *,'Minimum MAF for the QTL : ',minmaf
      if (minmaf.lt.0..or.minmaf.gt.0.5) stop 'Maf sup 0.5'
      
      print *,'Do we keep QTL in the markers file ? (o/n)'
      optmq='o'
      read (5,*,iostat=io) optmq
      if (optmq.ne.'n'.and.optmq.ne.'o') stop 'o ou n'
      if (optmq.eq.'n') print *,'QTLs are removed from markers'
      if (optmq.eq.'o') print *,'QTLs are included in markers'

      print *,'Format of typing/phases files input (P ou T) ?'
      fi='P'
      read (5,*,iostat=io) fi
      if (fi.eq.'p') fi='P'; if (fi.eq.'t') fi='T'; 
      print *,'Input file : ',fi
      if (fi.ne.'P'.and.fi.ne.'T') stop 'file format dif from P and T'

      print *,'Format of typing/phases files output (P ou T) ?'
      fo='T'
      read (5,*,iostat=io) fo
      if (fo.eq.'p') fo='P'; if (fo.eq.'t') fo='T'; 
      print *,'Output file : ',fo
      if (fo.ne.'P'.and.fo.ne.'T') stop 'file format dif from P and T'

      print *,'Seed ?'
      read (5,*) seed2
      print *,'Seed : ',seed2


! Initialization of the seed
      call init_random_seed(seed2)

! Do we need other formats ?           
      if (namechrs.eq.'no'.and.fi.ne.fo) then
        print *,'No markers file asked '
        print *,'while the I/O formats are different ?'
        stop
        end if
      if (namechrs.eq.'no'.and.optmq.eq.'n') then
        print *,'No markers file asked '
        print *,'while the QTLs are removed from the markers?'
        stop
        end if

      tecrit=.true.
      if (namechrs.eq.'no') tecrit=.false.
      if (fi.eq.fo.and.namechrs.ne.'no'.and.optmq.eq.'o') then
        print *,'No need to rewrite the markers file '
        print *,'it is not achieved'
        tecrit=.false.
        end if

      
! Read the pedigree file to check that it is correct (irc from 1 to n, parent<product) and count
      ng=0
      open (2,file=ficped,form='formatted')      
      do 
       read (2,*,iostat=io) irc,prc,mrc
       if (io.ne.0) exit
       ng=ng+1
       if (ng.ne.irc) stop 'individuals not from 1 to n in the pedigree'
       if (prc.ge.irc) stop 'father > individuals'
       if (mrc.ge.irc) stop 'mother > individuals'
       end do
      close (2)
      print *,'Number of individuals in the pedigree : ',ng
      
! y=phenotype, g=genetic value, typ=typed or untyped indicator
      allocate (y(ng),g(ng),iani(ng),pere(ng),mere(ng),typ(ng))
      y(:)=0. ; typ(:)=.false.

! 2nd reading of pedigree and stockage files
      open (2,file=ficped,form='formatted')      
      do i=1,ng
       read (2,*) irc,prc,mrc
       if (prc.lt.0) prc=0; if (mrc.lt.0) mrc=0
       pere(i)=prc ; mere(i)=mrc
       end do
      close (2)
      

! QTL variance
      nq=nqg+nqm+nqp
      print *,'QTL numbers : ', nq

      if (v1p.le.0.) stop 'Variance of small QTLS <= 0.'
      
! We report the variances of QTL to that of the small QTL
      v1g=v1g/v1p
      v1m=v1m/v1p

! xx variance explained by QTL, in small QTL equivalent
      xx=nqg*v1g + nqm*v1m + nqp
      x=vart*h2*pcqtl/xx  ! part des petits qtl
      
! vq=variance of each QTL, chr=QTL chromosome, imc=intra chromosome position, imt=overall position, effet=effect, all=positive effect allele
      allocate (vq(nq),chr(nq),imc(nq),imt(nq),effet(nq),all(nq))
      do i=1, nqg; vq(i)=v1g*x;      end do
      do i=1, nqm; vq(nqg+i)=v1m*x;  end do
      do i=1, nqp; vq(nqg+nqm+i)=x; end do


! we calculate the number of markers per chromosome and the allelic frequency
! nmc = nb marq / chrom, nmcum = numbers of cumulated markers in previous chromosomes, ncuq = ?
      allocate (nmc(nchrom),nmcum(nchrom+1),ncuq(nchrom+1))
      do i=1, nchrom; nmc(i)=0; end do
      nm=0 ; ntyp=0 ; nmax=0

! Count the numbers of markers
      DO ICHR = 1, nchrom
       if (nchrom1.eq.0) then
         fichchr=namechr
       else
         write(fichchr,'(a,i0)') trim(namechr),ichr ; 
         end if
       open(3,file=fichchr,form='formatted')
       call lire_ligne(3,ligne,separateur,nc,finfich)                                        
       if (finfich) stop 'Le fichier de phenotypes est vide'
       close(3)
       call lecligne(ligne,nc)    ; 
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
      print *,'Notal number of markers : ',nm
      print *,'Maximum number of markers per chromosome: ',nmax
      nmcum(nchrom+1)=nm

! freq = frequencies of the 2nd allele per marker, phe = genotype for one individual and one chromosome
      allocate (freq(nm),phe(2,nmax))
      freq(:)=0.
      
!c exploring genotype files to count individuals,fill typ and see if there are coherent
      DO ICHR = 1, nchrom
       if (nchrom1.eq.0) then
         fichchr=namechr
       else
         write(fichchr,'(a,i0)') trim(namechr),ichr ; 
         end if
      open(3,file=fichchr,form='formatted')  
      io=0 ; nt=0 ; ntyp=0
      do while(io.eq.0)
        read(3,*,iostat=io) irc
        if (io.ne.0) exit
        nt=nt+1
        if (irc.gt.ng) stop 'Individual missing from pedigree'
        typ(irc)=.true. 
        end do
      close (3)
      if (ichr.eq.1) then
        ntt=nt
      else if (nt.ne.ntt) then
        stop 'Non coherent genotype file'
        end if
      end do
      print *, 'Number of genotype lines : ',ntt
      ntyp=nt
      if (fi.eq.'P') then
       ntyp=ntyp/2
       print *, 'Number of indivuals : ',ntyp
       end if


! Read genotypes to compute frequencies
      DO ICHR = 1, nchrom
       if (nchrom1.eq.0) then
         fichchr=namechr
       else
         write(fichchr,'(a,i0)') trim(namechr),ichr ; 
         end if
      open(3,file=fichchr,form='formatted')  
      io=0  ; nc=nmc(ichr) ; nm1=nmcum(ichr)
      do ia=1, ntyp
       if (fi.eq.'P') then   
        read(3,*) irc,i,(phe(1,j),j=1,nc)
        if (i.ne.1) stop 'typages no phases'
        read(3,*) irc,i,(phe(2,j),j=1,nc)
        if (i.ne.2) stop 'typages no phases'
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
	
! End of frequencies computation, saving 
      do i=1, nm
        freq(i)=freq(i)/dfloat(2*ntyp)
	end do

      open(9,file=fichfreq)
      do i=1, nm
        write(9,'(i6,f8.3)') i,freq(i)
        end do
      close(9)

! we determine the QTLs among the markers, then its characteristics
      allocate(listq(nm),listm(nm))
      listq(:)=0 ; listm(:)=0

! we determine the np potential marker to be QTL, based on the MAF
      allocate(listpot(nm),chrom(nm),pos(nm))
      np=0 ; i=0
      do ichr=1, nchrom
       do ipos=1, nmc(ichr) 
        i=i+1
        chrom(i)=ichr; pos(i)=ipos        
        if (freq(i).ge.minmaf.and.freq(i).le.(1.-minmaf)) then  ; 
          np=np+1
          listpot(np)=i  ; ! listpot = listof potential markers
          end if
        end do
       end do
      if (np.le.nq) stop 'moins de marqueurs potentiels que de qtl'
      print *,'Number of markers with non extrem MAF : ',np

! informations on the nq QTLs : effect (> or <0), all=positive allele effect, chr=chromosome, imc=intra chromosome position, imt=overall position, vq=variance
! listq(num_marq)=num_qtl

! we randomly chose the QTLs in the np possibilites, we keep the info, remove it from the list and continue with the next marker
      do i=1, nq
!        j=1+np*rand(x)
        call random_number(x)
        j=1 + np*x
	if (j.lt.1.or.j.gt.np) stop 'Error in the QTL drawing'  
	listq(listpot(j))=i
        imt(i)=listpot(j)
        iq=listpot(j)

        do k=j+1, np
          listpot(k-1)=listpot(k)
          end do
        np=np-1

	chr(i)=chrom(iq)
	imc(i)=pos(iq)

! draw of the positive allele effect
	effet(i)=.5*sqrt(vq(i)/(2.*freq(iq)*(1.-freq(iq))))
! draw if positive or negative effect
	all(i)=1
        call random_number(x)
	if (x.gt.0.5)then
          all(i)=2 
          effet(i)=-effet(i)
          end if

	end do
	
! saving of QTLs info in fichqtl file
      open (4,file=fichqtl,form='formatted')
      do i=1, nq
        write (4,200) i,chr(i),imc(i),imt(i),vq(i),freq(imt(i)),effet(i)
  200 format (I6,i3,2i7,3f12.5)	
	end do
      close(4)	

! effects mini maxi 	
      xmin=effet(1); xmax=xmin
      do i=2, nq
	  if (effet(i).lt.xmin) xmin=effet(i)
	  if (effet(i).gt.xmax) xmax=effet(i)
	  end do
      print *,'Effects mini / maxi : ',xmin, xmax
	

! without QTL markers file 
      if (optmq.eq.'n') then
       j=0
       do i=1, nm
        if (listq(i).eq.0) then
	  j=j+1
	  listm(i)=j ;  ! listm(1:nm) 0 for the QTLs, 1:nm-nq for non QTL markers?
	  end if
	end do  
      else
       do i=1, nm
         listm(i)=i
	 end do
       end if

! alq=QTLs allele, listqc includes QTLs index for a given chromosome, listmc includes non QTL markers
      allocate (alq(2,nq,ntyp),listqc(nq),listmc(nm))

! rereading chromosomes for genotypes use   
      DO ICHR = 1, nchrom
       if (nchrom1.eq.0) then
         fichchr=namechr
         fichchrs=namechrs
       else
         write(fichchr,'(a,i0)') trim(namechr),ichr ; 
         write(fichchrs,'(a,i0)') trim(namechrs),ichr ; 
         end if
       print *,'lu    : ', fichchr (1:len_trim(fichchr))	
       print *,'ecrit : ', fichchrs(1:len_trim(fichchrs))

! research of the QTL list on the chromosome, saving of their index in listqc 
       nqc=0	
       do iq=1, nq
	  if (chr(iq).eq.ichr) then
	    nqc=nqc+1
	    listqc(nqc)=iq
	    end if
	  end do  

! list of no QTL markers in the chromosome for svg
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

! reading of individuals genotypes
      nc=nmc(ichr) ; nm1=nmcum(ichr)
      do ia=1, ntyp
       if (fi.eq.'P') then   
        read(3,*) irc,i,(phe(1,j),j=1,nc)
        read(3,*) irc,i,(phe(2,j),j=1,nc)
       else
        read(3,*) irc,((phe(i,j),i=1,2),j=1,nc)
        end if

! QTLs using, effects add to y 
       do im=1, nqc
         iq=listqc(im)	
! QTL effects on phenotypes
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

! saving of QTL alleles, in QTL order (we save it at the end after processing all the chromosomes)
         alq(1,iq,ia)=phe(1,imc(iq))    
         alq(2,iq,ia)=phe(2,imc(iq))
	 end do

! rewritting of markers if asked 
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
      
! saving of QTL genotypes      
      open (4,file=fichtypq,form='formatted')
      do i=1, ntyp
        write (4,100) iani(i),((alq(k,j,i),k=1,2),j=1,nq)
	end do
      close(4)	      
	 
! add polygenic value	 
      vg=h2*vart*(1.-pcqtl)  ! polygenic variance
      sdg=sqrt(vg)           ! polygenic standard deviation
      sdg1=sqrt(0.75*vg)     ! polygenic standard deviation if one known parent
      sdg2=sqrt(0.5*vg)      ! polygenic standard deviation if two known parents
      ve=sqrt(vart*(1.-h2))  ! residuals standard deviation
      
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
	
! phenotypes saving	
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
 	   

! lines reading ****************************************
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
             ligne(i:i)='_'
           else if (ligne(i:i).eq.separ) then
             ligne(i:i)=' '
             nc=nc+1
             end if
           end do

      return 
                                                   
    5 fin=.true.
      return

      end                        
 
! **************************************************************************
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
       
! **************************************************************************  
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
