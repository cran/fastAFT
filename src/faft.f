      subroutine faft(x,dlt,z,wt,size,npred,
     &     epl,epu,
     &     para,intp,lzr,rsdl,rfrac,zbar,proj,uppm,drct,lp,atrisk,
     &     lrgh,iew,bt,va,chis,ci,ynci,
     &     cef,cva,zsum,zssq,vec1,vec2,inde,mat)
      integer size,npred,dlt(size)
      double precision x(size),z(size,npred),wt(size),
     &     epl,epu,para(0:npred)
      integer intp(npred+1),lzr(size)
      double precision rsdl(size),rfrac(npred+1),
     &     zbar(0:npred),proj(0:npred,npred+1),
     &     uppm(npred+1,npred+1),drct(0:npred),lp(size),
     &     atrisk(0:npred,0:npred)

      integer lrgh,iew,ynci
      double precision bt(npred,3),va(npred,npred),chis(3),
     &     ci(npred,2)
      double precision cef(npred),cva(npred,npred),
     &     zsum(npred,size,2),zssq(npred,npred,size,2),
     &     vec1(npred),vec2(npred),
     &     inde(npred,npred),mat(npred,npred)

      logical ini,cont,succ
      integer i,j,k,p,q,tmp,cnt,ifr
      double precision chi,dtmp,dcmp

      iew=0
      ini=.true.
      do i=1,size
         wt(i)=1.0d0
      enddo
      call aqm(x,dlt,z,wt,size,npred,
     &     epl,epu,bt(1,3),va,succ,
     &     para,intp,lzr,rsdl,rfrac,zbar,proj,uppm,drct,lp,atrisk,
     &     zssq)
      if(.not. succ) then
         iew=1
         return
      endif
      call efmmt(x,dlt,z,size,npred,bt(1,3),lrgh,ini,
     &     cef,cva,chis(3),lzr,zsum,zssq,rsdl,atrisk,proj,uppm)
      do p=1,npred
         vec1(p)=cef(p)/dble(size)
      enddo

      do q=1,npred
         do p=1,npred
            bt(p,2)=bt(p,3)+va(p,q)
         enddo

         call efmmt(x,dlt,z,size,npred,bt(1,2),lrgh,ini,
     &        cef,cva,chi,lzr,zsum,zssq,rsdl,atrisk,proj,uppm)
         do p=1,npred
            mat(p,q)=cef(p)/dble(size)-vec1(p)
         enddo
      enddo

      call inverse(mat,npred,ifr,proj,uppm)
      if(ifr .eq. 0) iew=2

      do p=1,npred
         do q=1,npred
            inde(p,q)=0.0d0
            do j=1,npred
               inde(p,q)=inde(p,q)+va(p,j)*mat(j,q)
            enddo
         enddo
      enddo

      do p=1,npred
         bt(p,2)=bt(p,3)
      enddo
      call newton(x,dlt,z,size,npred,bt(1,2),lrgh,ini,cef,
     &     cva,chis(2),lzr,zsum,zssq,rsdl,va,proj,uppm,
     &     0,inde,1,vec1,vec2)

      do p=1,npred
         bt(p,1)=bt(p,2)
      enddo
      call newton(x,dlt,z,size,npred,bt(1,1),lrgh,ini,cef,
     &     cva,chis(1),lzr,zsum,zssq,rsdl,va,proj,uppm,
     &     0,inde,0,vec1,vec2)

      call efmmt(x,dlt,z,size,npred,bt(1,1),lrgh,ini,
     &     cef,cva,chi,lzr,zsum,zssq,rsdl,va,proj,uppm)
      call sandwich(npred,size,inde,cva,va,mat)

      if(ynci .ne. 1 .or. iew .ne. 0) return
      tmp=1
      do i=1,npred
         do k=1,2
            tmp=-tmp

            cnt=0
            do while(cnt .eq. 0 .or. (dtmp .lt. -1.0d-10 .and.
     &           cnt .le. 20))
               dtmp=dble(tmp)*1.959964d0*dsqrt(va(i,i))
     &              *1.1d0**dble(cnt)
               do p=1,npred
                  if(p .eq. i) then
                     rfrac(p)=bt(p,1)+dtmp
                  else
                     rfrac(p)=bt(p,1)+dtmp*va(i,p)/va(i,i)
                  endif
               enddo

               call newton(x,dlt,z,size,npred,rfrac,lrgh,ini,cef,
     &              cva,chi,lzr,zsum,zssq,rsdl,mat,proj,uppm,
     &              i,inde,0,zbar,drct)
               dtmp=dsqrt(chi)-1.959964d0
               cnt=cnt+1
            enddo

            if(dtmp .lt. -1.0d-10) iew=3

            cont=.true.
            if(dtmp .lt. 1.0d-10) cont=.false.

            do while(cont)
               do p=1,npred
                  vec1(p)=(rfrac(p)-bt(p,1))
     &                 *(1.0d0-1.959964d0/dsqrt(chi))
               enddo

               cnt=0
               do while(cnt .eq. 0 .or. (dcmp .lt. -1.0d-10 .and.
     &              cnt .le. 20))
                  do p=1,npred
                     vec2(p)=rfrac(p)-vec1(p)/2.0d0**dble(cnt)
                  enddo

                  call newton(x,dlt,z,size,npred,vec2,lrgh,ini,cef,
     &                 cva,chi,lzr,zsum,zssq,rsdl,mat,proj,uppm,
     &                 i,inde,0,zbar,drct)
                  dcmp=dsqrt(chi)-1.959964d0
                  cnt=cnt+1
               enddo
               cnt=cnt-1

               if(cnt .ge. 20) then
                  cont=.false.
               else
                  do p=1,npred
                     rfrac(p)=vec2(p)
                  enddo
                  if(dcmp .lt. 1.0d-10 .or.
     &                 dabs(vec1(i))/2.0d0**dble(cnt) .lt. 1.0d-6)
     &                 cont=.false.
               endif
            enddo

            ci(i,k)=rfrac(i)

         enddo
      enddo

      end

      subroutine inverse(a,dim,ifr,b,c)
      integer dim,ifr
      double precision a(dim,dim),b(dim,dim),
     &     c(dim,dim)

      integer p,q,n

      ifr=1
      do p=1,dim
         do q=1,p-1
            b(q,p)=0.0d0
            do n=1,dim
               b(q,p)=b(q,p)+a(n,p)*a(n,q)
            enddo
            do n=1,dim
               a(n,p)=a(n,p)-b(q,p)*a(n,q)
            enddo
         enddo
         b(p,p)=0.0d0
         do n=1,dim
            b(p,p)=b(p,p)+a(n,p)**2
         enddo
         b(p,p)=dsqrt(b(p,p))
         if(b(p,p) .gt. 1.0d-10) then
            do n=1,dim
               a(n,p)=a(n,p)/b(p,p)
            enddo
         else
            ifr=0
            b(p,p)=1.0d0
         endif
      enddo

      do p=dim,1,-1
         do q=1,p-1
            c(p,q)=0.0d0
         enddo
         c(p,p)=1.0d0/b(p,p)
         do q=p+1,dim
            c(p,q)=0.0d0
            do n=p+1,q
               c(p,q)=c(p,q)-b(p,n)*c(n,q)
            enddo
            c(p,q)=c(p,q)/b(p,p)
         enddo
      enddo

      do p=1,dim
         do q=1,dim
            b(p,q)=0.0d0
            do n=1,dim
               b(p,q)=b(p,q)+c(p,n)*a(q,n)
            enddo
         enddo
      enddo

      do p=1,dim
         do q=1,dim
            a(p,q)=b(p,q)
         enddo
      enddo

      end

      subroutine newton(x,dlt,z,size,npred,b,lrgh,ini,cef,cva,
     &     chi,ord,zsum,zssq,rsdl,mat1,mat2,mat3,
     &     elm,inde,once,vec1,vec2)
      logical ini
      integer size,npred,dlt(size),lrgh,ord(size)
      double precision x(size),z(size,npred),b(npred),
     &     cef(npred),cva(npred,npred),chi,zsum(npred,size,2),
     &     zssq(npred,npred,size,2),
     &     rsdl(size),mat1(npred,npred),mat2(npred,npred),
     &     mat3(npred,npred)
      integer elm,once
      double precision inde(npred,npred),vec1(npred),vec2(npred)

      logical cont
      integer p,q,cnt
      double precision dcmp,dtmp

      call efmmt(x,dlt,z,size,npred,b,lrgh,ini,cef,cva,dcmp,
     &     ord,zsum,zssq,rsdl,mat1,mat2,mat3)
      cont=.true.
      do while(cont)
         do p=1,npred
            vec1(p)=0.0d0
            do q=1,npred
               vec1(p)=vec1(p)+inde(p,q)*(cef(q)/dble(size))
            enddo
         enddo
         if(elm .ne. 0) then
            call sandwich(npred,size,inde,cva,mat1,mat2)
            do p=1,npred
               if(p .ne. elm) vec1(p)=vec1(p)-
     &              vec1(elm)*mat1(p,elm)/mat1(elm,elm)
            enddo
            vec1(elm)=0.0d0
         endif

         cnt=0
         do while(cnt .eq. 0 .or.
     &        (chi .ge. dcmp .and. cnt .le. 20))
            do p=1,npred
               vec2(p)=b(p)-vec1(p)/2.0d0**dble(cnt)
            enddo
            call efmmt(x,dlt,z,size,npred,vec2,lrgh,ini,cef,cva,chi,
     &           ord,zsum,zssq,rsdl,mat1,mat2,mat3)
            cnt=cnt+1
         enddo
         cnt=cnt-1

         dtmp=0.0d0
         do p=1,npred
            dtmp=dtmp+vec1(p)**2.0d0
         enddo
         dtmp=dsqrt(dtmp)/2.0d0**dble(cnt)

         if(cnt .eq. 20) then
            cont=.false.
            chi=dcmp
         else
            do p=1,npred
               b(p)=vec2(p)
            enddo
            if(once .eq. 1 .or. dtmp .lt. 1.0d-6) cont=.false.
            dcmp=chi
         endif
      enddo

      end

      subroutine sandwich(npred,size,inde,cva,va,a)
      integer npred,size
      double precision inde(npred,npred),cva(npred,npred),
     &     va(npred,npred),a(npred,npred)

      integer i,p,q

      do p=1,npred
         do q=1,p
            va(p,q)=cva(p,q)/dble(size)**2.0d0
         enddo
      enddo
      do p=1,npred
         do q=p+1,npred
            va(p,q)=va(q,p)
         enddo
      enddo

      do p=1,npred
         do q=1,npred
            a(p,q)=0.0d0
            do i=1,npred
               a(p,q)=a(p,q)+inde(p,i)*va(i,q)
            enddo
         enddo
      enddo
      do p=1,npred
         do q=1,p
            va(p,q)=0.0d0
            do i=1,npred
               va(p,q)=va(p,q)+a(p,i)*inde(q,i)
            enddo
         enddo
      enddo
      do p=1,npred
         do q=p+1,npred
            va(p,q)=va(q,p)
         enddo
      enddo

      end

      subroutine aqm(x,dlt,z,wt,size,npred,
     &     epl,epu,ibt,ibd,succ,
     &     para,intp,lzr,rsdl,rfrac,zbar,proj,uppm,drct,lp,atrisk,
     &     work)
      logical succ
      integer size,npred,dlt(size)
      double precision x(size),z(size,npred),wt(size),
     &     epl,epu,ibt(npred),ibd(npred,npred)
      integer intp(npred+1),lzr(size)
      double precision para(0:npred),rsdl(size),
     &     rfrac(npred+1),zbar(0:npred),proj(0:npred,npred+1),
     &     uppm(npred+1,npred+1),drct(0:npred),lp(size),
     &     atrisk(0:npred,0:npred),work(0:npred,0:npred)

      double precision logdet

      integer nintp,eintp,bd,ninf

      logical cont
      integer i,p,q
      double precision taul,cdf,inc,irisk,dtmp

      para(0)=x(1)
      do i=2,size
         para(0)=dmin1(para(0),x(i))
      enddo
      do p=1,npred
         para(p)=0.0d0
      enddo
      do i=1,size
         rsdl(i)=x(i)-para(0)
         if(rsdl(i) .gt. 1.0d-10) then
            lzr(i)=2
         else
            lzr(i)=1
         endif
      enddo
      do p=0,npred
         zbar(p)=0.0d0
      enddo
      do i=1,size
         zbar(0)=zbar(0)+wt(i)
         do p=1,npred
            zbar(p)=zbar(p)+z(i,p)*wt(i)
         enddo
      enddo
      nintp=0

      cdf=0.0d0
      do p=0,npred
         do q=0,p
            atrisk(p,q)=0.0d0
         enddo
      enddo
      do i=1,size
         if(dlt(i) .eq. 1) then
            p=0
            q=0
            atrisk(p,q)=atrisk(p,q)+wt(i)
            do p=1,npred
               q=0
               atrisk(p,q)=atrisk(p,q)+z(i,p)*wt(i)
               do q=1,p
                  atrisk(p,q)=atrisk(p,q)+z(i,p)*z(i,q)*wt(i)
               enddo
            enddo
         endif
      enddo
      irisk=-1.0d0

      do p=1,npred
         ibt(p)=0.0d0
         do q=1,npred
            ibd(p,q)=0.0d0
         enddo
      enddo

      taul=-1.0d0
      succ=.true.

      cont=.true.
      do while(cont)
         call minstep(rsdl,dlt,z,wt,size,npred,zbar,
     &        intp,rfrac,nintp,lzr,para,bd,
     &        proj,uppm,drct,lp)

         if(irisk .lt. -5.0d-1) then
            irisk=logdet(atrisk,npred+1,work,ninf)
            if(ninf .eq. 1) then
               succ=.false.
               return
            endif
            dtmp=1.0d0
         else
            dtmp=logdet(atrisk,npred+1,work,ninf)-irisk
         endif
         if(taul .lt. -5.0d-1 .and. (ninf .eq. 1 .or.
     &        dtmp .lt. dlog(epl)*dble(npred+1))) taul=cdf

         if(ninf .eq. 1 .or.
     &        dtmp .le. dlog(epu)*dble(npred+1)) cont=.false.

         if(cont) then
            eintp=0
            do i=1,nintp
               eintp=eintp+dlt(intp(i))
            enddo

            if(eintp .eq. 0) then
               inc=1.0d0
            else
               call incstep(dlt,z,wt,size,npred,zbar,
     &              intp,rfrac,eintp,nintp,lzr,
     &              inc,proj,uppm,atrisk,drct)
            endif

            if(taul .gt. -5.0d-1) then
               dtmp=(1.0d0-cdf)*inc
               do p=1,npred
                  ibt(p)=ibt(p)+para(p)*dtmp
                  do q=1,p
                     ibd(p,q)=ibd(p,q)+para(p)*para(q)*dtmp
                  enddo
               enddo
            endif

            cdf=1.0d0-(1.0d0-cdf)*(1.0d0-inc)

            if(cdf .gt. 1.0d0-1.0d-10)
     &           cont=.false.
         endif
      enddo

      if(taul .lt. -5.0d-1 .or. cdf-taul .le. 0.0d0) then
         succ=.false.
         return
      endif

      do p=1,npred
         ibt(p)=ibt(p)/(cdf-taul)
      enddo
      do p=1,npred
         do q=1,p
            ibd(p,q)=ibd(p,q)/(cdf-taul)-ibt(p)*ibt(q)
         enddo
      enddo
      call cholesky(ibd,npred,ibd)

      end

      double precision function logdet(atrisk,dim,mat,ninf)
      integer dim,ninf
      double precision atrisk(dim,dim),mat(dim,dim)

      integer i,j,k

      do i=1,dim
         do j=1,i
            mat(i,j)=atrisk(i,j)
         enddo
      enddo
      do i=1,dim
         do j=i+1,dim
            mat(i,j)=mat(j,i)
         enddo
      enddo

      logdet=0.0d0
      ninf=0
      i=1
      do while(i .le. dim)
         if(i .ge. 2) then
            do j=1,dim
               do k=1,min0(i,j)-1
                  mat(i,j)=mat(i,j)-mat(i,k)*mat(k,j)
               enddo
               if(j .lt. i) mat(i,j)=mat(i,j)/mat(j,j)
            enddo
         endif
         if(mat(i,i) .lt. 1.0d-10) then
            ninf=1
            return
         endif
         i=i+1
      enddo

      logdet=dlog(mat(1,1))
      do i=2,dim
         logdet=logdet+dlog(mat(i,i))
      enddo

      end

      subroutine incstep(dlt,z,wt,size,npred,zbar,
     &     intp,rfrac,eintp,nintp,lzr,inc,
     &     proj,uppm,atrisk,drct)
      integer size,npred,dlt(size),intp(npred+1),eintp,
     &     nintp,lzr(size)
      double precision z(size,npred),wt(size),rfrac(npred+1),
     &     zbar(0:npred),inc
      double precision proj(0:npred,npred+1),
     &     uppm(npred+1,npred+1),atrisk(0:npred,0:npred),
     &     drct(npred+1)

      logical first
      integer i,j,p,q,k,l
      double precision dtmp

      do i=1,nintp
         drct(i)=0.0d0
         do p=0,npred
            drct(i)=drct(i)+zbar(p)*proj(p,i)
         enddo
      enddo
      k=0
      i=nintp
      do while(k .lt. eintp)
         do j=i+1,nintp
            drct(i)=drct(i)-uppm(i,j)*drct(j)
         enddo
         drct(i)=drct(i)/uppm(i,i)
         k=k+dlt(intp(i))
         i=i-1
      enddo

      inc=1.0d0
      i=1
      do while(i .le. nintp)
         if(dlt(intp(i)) .eq. 1 .and. dabs(drct(i)) .gt. 1.0d-10) then
            if(drct(i) .gt. 0.0d0) then
               dtmp=rfrac(i)*wt(intp(i))/drct(i)
            else
               dtmp=(rfrac(i)-1.0d0)*wt(intp(i))/drct(i)
            endif
            if(dtmp .gt. 1.0d-50) then
               inc=dmin1(inc,dtmp)
            else
               do j=i+1,nintp
                  intp(j-1)=intp(j)
                  rfrac(j-1)=rfrac(j)
                  drct(j-1)=drct(j)
               enddo
               nintp=nintp-1
               i=i-1
            endif
         endif
         i=i+1
      enddo

      do p=0,npred
         zbar(p)=zbar(p)*(1.0d0-inc)
      enddo
      first=.true.
      i=1
      l=nintp+1
      do while(i .le. nintp)
         if(dlt(intp(i)) .eq. 1) then
            rfrac(i)=rfrac(i)-drct(i)/wt(intp(i))*inc
            p=0
            q=0
            atrisk(p,q)=atrisk(p,q)-drct(i)*inc
            do p=1,npred
               q=0
               atrisk(p,q)=atrisk(p,q)-z(intp(i),p)
     &              *drct(i)*inc
               do q=1,p
                  atrisk(p,q)=atrisk(p,q)-z(intp(i),p)
     &                 *z(intp(i),q)*drct(i)*inc
               enddo
            enddo

            if(rfrac(i) .lt. 1.0d-10 .or.
     &           rfrac(i) .gt. 1.0d0-1.0d-10) then
               if(first) l=i
               first=.false.
               lzr(intp(i))=1
               if(rfrac(i) .lt. 1.0d-10) lzr(intp(i))=-1
               do j=i,nintp-1
                  drct(j)=drct(j+1)
                  intp(j)=intp(j+1)
                  rfrac(j)=rfrac(j+1)
               enddo
               nintp=nintp-1
               i=i-1
            endif
         endif
         i=i+1
      enddo

      i=1
      j=nintp
      first=.true.
      do while(i .lt. j)
         if(dlt(intp(i)) .eq. 0) then
            do while(j .gt. i .and. dlt(intp(j)) .eq. 0)
               j=j-1
            enddo
            if(i .lt. j) then
               if(first) then
                  first=.false.
                  if(i .lt. l) l=i
               endif
               k=intp(i)
               intp(i)=intp(j)
               intp(j)=k
               dtmp=rfrac(i)
               rfrac(i)=rfrac(j)
               rfrac(j)=dtmp
               i=i+1
               j=j-1
            endif
         else
            i=i+1
         endif
      enddo

      call orth(z,size,npred,intp,proj,uppm,l,nintp)

      end

      subroutine minstep(rsdl,dlt,z,wt,size,npred,zbar,
     &     intp,rfrac,nintp,lzr,beta,bd,proj,uppm,
     &     drct,lp)
      integer size,npred,dlt(size),intp(npred+1),nintp,
     &     lzr(size),bd
      double precision rsdl(size),z(size,npred),wt(size),
     &     zbar(0:npred),rfrac(npred+1),beta(0:npred)
      double precision proj(0:npred,npred+1),
     &     uppm(npred+1,npred+1),drct(0:npred),lp(size)

      integer ocode(2),i,p,mark,cyc,nconst,drop

      bd=0

      mark=nintp
      do while(mark .gt. 0 .and. dlt(intp(mark)) .eq. 0)
         mark=mark-1
      enddo

      cyc=npred+1
      nconst=0
      do while(nconst .lt. cyc-mark)
         drop=-1
         if(nintp .eq. cyc) then
            lzr(intp(mark+1))=1
            if(dlt(intp(mark+1)) .eq. 0) then
               zbar(0)=zbar(0)+wt(intp(mark+1))
     &              *(1-rfrac(mark+1))
               do p=1,npred
                  zbar(p)=zbar(p)+z(intp(mark+1),p)
     &                 *wt(intp(mark+1))*(1-rfrac(mark+1))
               enddo
            else
               if(rfrac(mark+1) .lt. 1.0d-10)
     &              lzr(intp(mark+1))=-1
            endif
            drop=intp(mark+1)
            do i=mark+1,nintp-1
               intp(i)=intp(i+1)
               rfrac(i)=rfrac(i+1)
            enddo
            call orth(z,size,npred,intp,proj,uppm,mark+1,cyc-1)
            nintp=nintp-1
         endif

         if(nintp .lt. npred+1) intp(nintp+1)=drop

         call line(rsdl,dlt,z,wt,size,npred,proj,zbar,
     &        intp,rfrac,nintp,lzr,beta,ocode,drct,lp)
         if(ocode(1) .eq. 1) bd=1
         if(ocode(2) .eq. 0) then
            cyc=nintp
         else
            call orth(z,size,npred,intp,proj,uppm,nintp,nintp)
         endif

         if(ocode(2) .eq. 1 .and. drop .eq. intp(nintp)) then
            nconst=nconst+1
         else
            nconst=0
            if(ocode(2) .eq. 1 .and. cyc .eq. nintp) nconst=1
         endif

         if(cyc .lt. npred+1 .and. cyc .gt. mark .and.
     &        nconst .eq. cyc-mark) then
            intp(nintp+1)=-1
            call line(rsdl,dlt,z,wt,size,npred,proj,zbar,
     &           intp,rfrac,nintp,lzr,beta,ocode,drct,lp)
            if(ocode(1) .eq. 1) then
               bd=1
               if(ocode(2) .eq. 1) then
                  cyc=npred+1
                  call orth(z,size,npred,intp,proj,uppm,nintp,nintp)
                  nconst=0
                  if(cyc .eq. nintp) nconst=1
               endif
            else
               if(ocode(2) .eq. 1) nintp=nintp-1
            endif
         endif

      enddo

      end

      subroutine orth(z,size,npred,intp,proj,uppm,m,n)
      integer size,npred,intp(npred+1),m,n
      double precision z(size,npred),proj(0:npred,npred+1),
     &     uppm(npred+1,npred+1)
      integer i,j,p

      do i=m,n
         proj(0,i)=1.0d0
         do p=1,npred
            proj(p,i)=z(intp(i),p)
         enddo
         do j=1,i-1
            uppm(j,i)=0.0d0
            do p=0,npred
               uppm(j,i)=uppm(j,i)+proj(p,j)*proj(p,i)
            enddo
            do p=0,npred
               proj(p,i)=proj(p,i)-uppm(j,i)*proj(p,j)
            enddo
         enddo
         uppm(i,i)=0.0d0
         do p=0,npred
            uppm(i,i)=uppm(i,i)+proj(p,i)**2
         enddo
         uppm(i,i)=dsqrt(uppm(i,i))
         do p=0,npred
            proj(p,i)=proj(p,i)/uppm(i,i)
         enddo
      enddo

      end

      subroutine line(rsdl,dlt,z,wt,size,npred,proj,zbar,
     &     intp,rfrac,nintp,lzr,beta,ocode,drct,lp)
      integer size,npred,dlt(size),intp(npred+1),nintp,lzr(size),
     &     ocode(2)
      double precision rsdl(size),z(size,npred),wt(size),
     &     proj(0:npred,npred+1),zbar(0:npred),rfrac(npred+1),
     &     beta(0:npred)
      double precision drct(0:npred),lp(size)

      logical first,ftry
      integer i,j,p,act
      double precision dtmp,dcmp

      ftry=.false.
      if(nintp .lt. npred+1 .and. intp(nintp+1) .gt. 0) ftry=.true.

      ocode(1)=0
      ocode(2)=0

 10   continue
      do p=0,npred
         drct(p)=zbar(p)
      enddo
      do i=1,nintp
         dtmp=0.0d0
         do p=0,npred
            dtmp=dtmp+proj(p,i)*drct(p)
         enddo
         do p=0,npred
            drct(p)=drct(p)-dtmp*proj(p,i)
         enddo
      enddo

      if(ftry) then
         ftry=.false.
         i=intp(nintp+1)
         call caseone(dlt,z,wt,size,npred,proj,zbar,
     &        intp,rfrac,nintp,lzr,ocode,drct,lp,i,act)
         if(act .eq. 1) return
         if(act .eq. 2) goto 10
      endif

      dtmp=0.0d0
      do p=0,npred
         dtmp=dtmp+drct(p)**2
      enddo
      if(dtmp .lt. 1.0d-20) return

 20   continue
      do i=1,size
         if(iabs(lzr(i)) .eq. 1 .and. dlt(i) .eq. 1) then
            call caseone(dlt,z,wt,size,npred,proj,zbar,
     &           intp,rfrac,nintp,lzr,ocode,drct,lp,i,act)
            if(act .eq. 1) return
         endif
      enddo
      do i=1,size
         if(iabs(lzr(i)) .eq. 1 .and. dlt(i) .eq. 0) then
            call caseone(dlt,z,wt,size,npred,proj,zbar,
     &           intp,rfrac,nintp,lzr,ocode,drct,lp,i,act)
            if(act .eq. 1) return
            if(act .eq. 2) goto 10
         endif
      enddo

      first=.true.
      do i=1,size
         if(iabs(lzr(i)) .eq. 2) then
            lp(i)=drct(0)
            do p=1,npred
               lp(i)=lp(i)+z(i,p)*drct(p)
            enddo
            if(lp(i)*dble(lzr(i)) .gt. 2.0d-10) then
               dcmp=rsdl(i)/lp(i)
               if(first .or. dcmp .lt. dtmp) then
                  first=.false.
                  dtmp=dcmp
               endif
            endif
         endif
      enddo
      if(.not. first) then
         do i=1,size
            if(lzr(i) .ne. 0) then
               rsdl(i)=rsdl(i)-dtmp*lp(i)
               j=1
               if(lzr(i) .lt. 0) j=-1
               if(dabs(rsdl(i)) .lt. 1.0d-10) then
                  lzr(i)=1
               else
                  lzr(i)=2
               endif
               lzr(i)=j*lzr(i)
            endif
         enddo
         do p=0,npred
            beta(p)=beta(p)+dtmp*drct(p)
         enddo
         ocode(1)=1

         goto 20
      endif

      end

      subroutine caseone(dlt,z,wt,size,npred,proj,zbar,
     &     intp,rfrac,nintp,lzr,ocode,drct,lp,i,act)
      integer size,npred,dlt(size),intp(npred+1),nintp,lzr(size),
     &     ocode(2),i,act
      double precision z(size,npred),wt(size),
     &     proj(0:npred,npred+1),zbar(0:npred),rfrac(npred+1)
      double precision drct(0:npred),lp(size)

      integer j,p
      double precision dtmp,dcmp

      act=0
      lp(i)=drct(0)
      do p=1,npred
         lp(i)=lp(i)+z(i,p)*drct(p)
      enddo
      dcmp=lp(i)*dble(lzr(i))
      if(dcmp .gt. 1.0d-10) then
         if(dlt(i) .eq. 1) then
            act=1
            nintp=nintp+1
            rfrac(nintp)=dble(lzr(i)+1)/2.0d0
            intp(nintp)=i
            lzr(i)=0
            ocode(2)=1
         else
            proj(0,npred+1)=1.0d0
            do p=1,npred
               proj(p,npred+1)=z(i,p)
            enddo
            do j=1,nintp
               dtmp=0.0d0
               do p=0,npred
                  dtmp=dtmp+proj(p,j)*proj(p,npred+1)
               enddo
               do p=0,npred
                  proj(p,npred+1)=proj(p,npred+1)
     &                 -dtmp*proj(p,j)
               enddo
            enddo
            dtmp=proj(0,npred+1)
            do p=1,npred
               dtmp=dtmp+proj(p,npred+1)*z(i,p)
            enddo
            dtmp=dtmp*wt(i)
            if(dcmp-dtmp .lt. -1.0d-10) then
               act=1
               nintp=nintp+1
               intp(nintp)=i
               dtmp=dcmp/dtmp
               zbar(0)=zbar(0)-dble(lzr(i))*wt(i)*dtmp
               do p=1,npred
                  zbar(p)=zbar(p)-dble(lzr(i))*z(i,p)
     &                 *wt(i)*dtmp
               enddo
               if(lzr(i) .eq. 1) then
                  rfrac(nintp)=1.0d0-dtmp
               else
                  rfrac(nintp)=dtmp
               endif
               lzr(i)=0
               ocode(2)=1
            else
               act=2
               zbar(0)=zbar(0)-wt(i)*dble(lzr(i))
               do p=1,npred
                  zbar(p)=zbar(p)-z(i,p)*wt(i)*dble(lzr(i))
               enddo
               lzr(i)=-lzr(i)
            endif
         endif
      endif

      end

      subroutine efmmt(x,dlt,z,size,npred,b,lrgh,ini,cef,cva,chi,
     &     ord,zsum,zssq,rsdl,mat1,mat2,mat3)
      logical ini
      integer size,npred,dlt(size),lrgh,ord(size)
      double precision x(size),z(size,npred),b(npred),
     &     cef(npred),cva(npred,npred),chi,zsum(npred,size,2),
     &     zssq(npred,npred,size,2),
     &     rsdl(size),mat1(npred,npred),mat2(npred,npred),
     &     mat3(npred,npred)

      integer i,j,tmp,p,q,ifr

      do i=1,size
         rsdl(i)=x(i)
         do p=1,npred
            rsdl(i)=rsdl(i)-b(p)*z(i,p)
         enddo
      enddo

      if(ini) then
         do i=1,size
            ord(i)=i
         enddo
      endif

      do i=2,size
         j=i
         do while(j .gt. 1 .and.
     &        (rsdl(ord(j-1)) .gt. rsdl(ord(j)) .or.
     &        (rsdl(ord(j-1)) .eq. rsdl(ord(j)) .and.
     &        (dlt(ord(j-1)) .lt. dlt(ord(j))))))
            if(.not. ini) then
               if(dlt(ord(j-1)) .eq. 1) then
                  call adsu(z,size,npred,lrgh,cef,cva,ord,
     &                 zsum,zssq,j-1,-1)
                  do p=1,npred
                     zsum(p,ord(j-1),1)=zsum(p,ord(j-1),1)-z(ord(j),p)
                     do q=1,p
                        zssq(p,q,ord(j-1),1)=zssq(p,q,ord(j-1),1)
     &                       -z(ord(j),p)*z(ord(j),q)
                     enddo
                  enddo
               endif
               if(dlt(ord(j)) .eq. 1) then
                  call adsu(z,size,npred,lrgh,cef,cva,ord,
     &                 zsum,zssq,j,-1)
                  do p=1,npred
                     zsum(p,ord(j),1)=zsum(p,ord(j),1)+z(ord(j-1),p)
                     do q=1,p
                        zssq(p,q,ord(j),1)=zssq(p,q,ord(j),1)
     &                       +z(ord(j-1),p)*z(ord(j-1),q)
                     enddo
                  enddo
               endif
            endif

            tmp=ord(j)
            ord(j)=ord(j-1)
            ord(j-1)=tmp

            if(.not. ini) then
               if(dlt(ord(j-1)) .eq. 1)
     &              call adsu(z,size,npred,lrgh,cef,cva,ord,
     &              zsum,zssq,j-1,1)
               if(dlt(ord(j)) .eq. 1)
     &              call adsu(z,size,npred,lrgh,cef,cva,ord,
     &              zsum,zssq,j,1)
            endif

            j=j-1
         enddo
      enddo

      if(ini) then
         do p=1,npred
            cef(p)=0.0d0
            zsum(p,ord(size),1)=z(ord(size),p)
            zsum(p,ord(size),2)=0.0d0
            do q=1,p
               cva(p,q)=0.0d0
               zssq(p,q,ord(size),1)=z(ord(size),p)*z(ord(size),q)
               zssq(p,q,ord(size),2)=0.0d0
            enddo
         enddo
         do i=size-1,1,-1
            do p=1,npred
               zsum(p,ord(i),1)=zsum(p,ord(i+1),1)+z(ord(i),p)
               do q=1,p
                  zssq(p,q,ord(i),1)=zssq(p,q,ord(i+1),1)
     &                 +z(ord(i),p)*z(ord(i),q)
               enddo
            enddo
            if(dlt(ord(i)) .eq. 1)
     &           call adsu(z,size,npred,lrgh,cef,cva,ord,zsum,zssq,
     &           i,1)
         enddo
         ini=.false.
      endif

      do p=1,npred
         do q=1,p
            mat1(p,q)=cva(p,q)
         enddo
      enddo
      do p=1,npred
         do q=p+1,npred
            mat1(p,q)=cva(q,p)
         enddo
      enddo
      call inverse(mat1,npred,ifr,mat2,mat3)
      chi=0.0d0
      do p=1,npred
         chi=chi+cef(p)**2.0d0*mat1(p,p)
         do q=1,p-1
            chi=chi+2.0d0*cef(p)*cef(q)*mat1(p,q)
         enddo
      enddo

      end

      subroutine adsu(z,size,npred,lrgh,cef,cva,ord,zsum,zssq,
     &     ind,as)
      integer size,npred,lrgh,ord(size),ind,as
      double precision z(size,npred),
     &     cef(npred),cva(npred,npred),zsum(npred,size,2),
     &     zssq(npred,npred,size,2)

      integer p,q

      if(as .eq. 1) then
         do p=1,npred
            zsum(p,ord(ind),2)=z(ord(ind),p)
     &           -zsum(p,ord(ind),1)/dble(size-ind+1)
            if(lrgh .eq. 1) zsum(p,ord(ind),2)=zsum(p,ord(ind),2)
     &           *dble(size-ind+1)/dble(size)
            cef(p)=cef(p)+zsum(p,ord(ind),2)
            do q=1,p
               zssq(p,q,ord(ind),2)=zssq(p,q,ord(ind),1)
     &              /dble(size-ind+1)
     &              -zsum(p,ord(ind),1)*zsum(q,ord(ind),1)
     &              /dble(size-ind+1)**2.0d0
               if(lrgh .eq. 1) zssq(p,q,ord(ind),2)
     &              =zssq(p,q,ord(ind),2)
     &              *(dble(size-ind+1)/dble(size))**2.0d0
               cva(p,q)=cva(p,q)+zssq(p,q,ord(ind),2)
            enddo
         enddo
      endif

      if(as .eq. -1) then
         do p=1,npred
            cef(p)=cef(p)-zsum(p,ord(ind),2)
            do q=1,p
               cva(p,q)=cva(p,q)-zssq(p,q,ord(ind),2)
            enddo
         enddo
      endif

      end

      subroutine cholesky(a,dim,l)
      integer dim
      double precision a(dim,dim),l(dim,dim)

      integer i,j,k

      do i=1,dim
         l(i,i)=a(i,i)
         do j=1,i-1
            l(i,i)=l(i,i)-l(i,j)**2.0d0
         enddo
         l(i,i)=dsqrt(l(i,i))
         do j=i+1,dim
            l(j,i)=a(j,i)
            do k=1,i-1
               l(j,i)=l(j,i)-l(j,k)*l(i,k)
            enddo
            l(j,i)=l(j,i)/l(i,i)
         enddo
      enddo

      end
