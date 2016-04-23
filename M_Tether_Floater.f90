Module M_Tether_Floater

  use M_General,              only: pi,Ly,Bp,dt,RoDrop,gz,ZFree,fileplace,Froude              
  use M_Solid,                only: Solid_Bot,Solid_ResB_Anal,Solid_Mas_Anal,Solid_S_Ini,Solid_H_Ini,Solid_TetAng_Ini,  &
                                  & R_Cyl_1,Point_Cyl_1,Xbar,YBar,ZBar,omegax,omegay,omegaz,ubar,vbar,wbar
  use M_Tower,                only: MNac,MTow,MRot
  use M_Math,                 only: LENGTHF,CROSS
  integer,parameter                               ::   tetn=1
  real(8),dimension(1:tetn,3)                     ::   point_Teth,point_Teth_old,point_Teth_old2
  real(8),dimension(1:tetn)                       ::   Ox,Oy,Oz,LZero,LNew     
  real(8)                                         ::   KTeth,ETeth,DTethO,DTethI,TickTeth,leg   
  real(8)                                         ::   L_Rig_Teth,M_Bal_Teth,Buoy_Bal_Teth   
    
    
  contains


    subroutine Tether_Constant_Ini()
      implicit none 

      
      Include "Par_Constant_Tether_Floater.txt"
       
       
      Ox(1)=Bp(1,1)   
      Oy(1)=Bp(1,2)     
      Oz(1)=Bp(1,3)-R_Cyl_1-L_Rig_Teth+(Solid_S_Ini**2)/(2*L_Rig_Teth**2)+Solid_H_Ini

      return
    end subroutine 


    subroutine Tether_Dynamic_Ini()
      implicit none
      integer :: i 
      Include "Par_Dynamic_Tether_Floater.txt"

      do i=1,tetn
        print*,"Tether",i, "length in the first/current  time step configuration is:", &
        &  LENGTHF( ox(i)-point_Teth(i,1),oy(i)-point_Teth(i,2),oz(1)-point_Teth(i,3) ) 
      end do 

      print*,"Initial tethers values are ended"

      return
    end subroutine



    subroutine Tether(FTether,MomTether)
      implicit none     
      real(8),dimension(3)                      ,intent(out)   :: FTether,MomTether
      real(8),dimension(3)                                     :: FSTether,MomSTether,rTmp
      integer i
     
      oz(1) =Point_Teth(1,3)-L_Rig_Teth &
      &     +abs(point_Teth(1,1)-ox(1))**2/(2*L_Rig_Teth**2)   

      FTether(:)=0
      MomTether(:)=0
      do i=1,tetn   !! number of tethers !!
    
        !FSTether(1)= (M_Bal_Teth*gz+Buoy_Bal_Teth*abs(gz))*( Ox(i)-point_Teth(i,1) )/( Oz(i)-point_Teth(i,3) )
        FSTether(1)= (M_Bal_Teth*gz+Buoy_Bal_Teth*abs(gz))*( point_Teth(i,1)-Ox(i) )/L_Rig_Teth
        FSTether(2)=0  !KTeth* ( Lnew(i)-Lzero(i) )* ( Oy(i)-point_Teth(i,2) )/Lnew(i)
        FSTether(3)= M_Bal_Teth*gz+Buoy_Bal_Teth*abs(gz)  !KTeth* ( Lnew(i)-Lzero(i) )* ( Oz(i)-point_Teth(i,3) )/Lnew(i)

      
        rTmp(1)=point_Teth(i,1)-XBar
        rTmp(2)=point_Teth(i,2)-YBar 
        rTmp(3)=point_Teth(i,3)-ZBar
        MomSTether(:)=CROSS(rTmp,FSTether)
        FTether(:)  =  FTether(:)+  FSTether(:)
        MomTether(:)=MomTether(:)+MomSTether(:)
      end do 

      MomTether(1)=0
      MomTether(3)=0

      PRINT*,'Tether vertical   force',FTether(3)/1000, 'KN'
      PRINT*,'Tether horizental force',FTether(1)/1000, 'KN'

      return
    end subroutine 



    subroutine Teth_Sec_Sav()
      implicit none 

      point_Teth_old(:,:)=point_Teth(:,:)

      return           
    end subroutine
 
    subroutine Teth_It_Sav()
      implicit none   
 
      point_Teth_old2(:,:)=point_Teth(:,:)

      return           
    end subroutine



    subroutine Teth_Sec_Cor()
      implicit none 

      point_Teth(:,:)=0.5*(point_Teth(:,:)+point_Teth_old(:,:))

      return    
    end subroutine 



    subroutine Teth_update()
      implicit none 

      call Teth_Update_Posi()
         
      return 
    end subroutine 


    subroutine Teth_Update_Posi()
      implicit none 
 
      real(8),dimension(tetn,3)                       ::   point_Teth_n
      integer                                         ::   i

      do i=1,tetn
        point_Teth_n(i,1)=point_Teth_old2(i,1) +dt*( ubar + omegay*(point_Teth(i,3)-zbar)  -omegaz*(point_Teth(i,2)-ybar) )
        point_Teth_n(i,2)=point_Teth_old2(i,2) +dt*( vbar - omegax*(point_Teth(i,3)-zbar)  +omegaz*(point_Teth(i,1)-xbar) )
        point_Teth_n(i,3)=point_Teth_old2(i,3) +dt*( wbar + omegax*(point_Teth(i,2)-ybar)  -omegay*(point_Teth(i,1)-xbar) )
      end do 

      point_Teth(:,:)=point_Teth_n(:,:)
      
      return      
    end subroutine 


    Subroutine Tether_Dynamic_Write()   
      implicit none 
      integer i
  
511   format (3(1x,e23.15))

      OPEN(unit=3045,file=fileplace//"Teth_Save.dat",STATUS='REPLACE')   
     
      do i=1,tetn
        write (3045,511) point_Teth(i,1),point_Teth(i,2),point_Teth(i,3)
      end do     
   

      close(3045)

      return
    end subroutine Tether_Dynamic_Write



    Subroutine Tether_Dynamic_Read()   
      implicit none 
      integer i
  
511   format (3(1x,e23.15))

      OPEN(unit=3045,file=fileplace//"Teth_Save.dat")   
     
      do i=1,tetn
        read (3045,511) point_Teth(i,1),point_Teth(i,2),point_Teth(i,3)
      end do     
   

      close(3045)

      return
    end subroutine Tether_Dynamic_Read


end module 






