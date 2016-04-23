  Module M_Tower

   use M_General,             only: pi,gz,dt,Froude,fileplace
   use M_Solid,               only: Point_Cyl_1,asolidx,AnAcY,xbar,ybar,zbar
   use M_math,                only: CROSS,LENGTHF

   real(8),save :: MRot,MNac,MTow,DiB,DiT,Di,TickB,TickT,Tick,TowEH,ETow,TurbH
   real(8),save :: OffLen,DYawRotor,DYawNacelle,RSpinV,Le,ICTow,ITow,IRotY,IRotX,RRot,PlCgH
   real(8)      :: gama(3),CorTowSp

   contains

    subroutine Tower(FTower,MomTower,NumImp,tp,Iterate,torder)
    implicit none    
    real(8),intent(out)  ,dimension(3)   :: FTower,MomTower
    integer,intent(in)                   :: tp,Numimp,Iterate,torder
         
    real(8),save         :: SolidTeta
    real(8),dimension(3) :: RM
     
    SolidTeta=ATAN( (point_Cyl_1(1,1)-point_Cyl_1(2,1))/(point_Cyl_1(1,3)-point_Cyl_1(2,3)) )
   
    if(Numimp.eq.1) then 
      Gama(3)=Gama(2) 
      Gama(2)=Gama(1)
    end if
 
    Gama(1)=2*Gama(2)-Gama(3)                                             &  
                 & +( dt*dt/(IRotY+ITow+(MRot+MNac)*TurbH**2) )            &
                 &*(                                                       &
                 & -  MTow             *TowEH          *aSolidX            &
                 & - (MRot+MNac)       *TurbH          *aSolidX            & 
                 & - (MTow)     * PlCgH*TowEH          *AnAcY              &
                 & - (MRot+MNac)* PlCgH*TurbH          *AnAcY              &
                                                                       
                 & +   MTow            *TowEH          *abs(gz)*Gama(2)    &
                 & +  (MRot+MNac)      *TurbH          *abs(gz)*Gama(2)    &
                  
                 & -  CorTowSp*3*ETow*ICTow/(TurbH)*(Gama(2)-SolidTeta)    &
                 & )
   
    FTower(1)= -(MNac+MRot)      *(TurbH*(Gama(1)-2*Gama(2)+Gama(3))/(dt**2)) &
            &  -MTow            *(TowEH*(Gama(1)-2*Gama(2)+Gama(3))/(dt**2)) &
            & -(MNac+MRot+MTow) *(ASolidX+PlCgH*AnAcY) 
    FTower(2)=0
    FTower(3)=-(MTow+MRot+MNac)*abs(gz)   !! assuming the centifugal acceleration which is a nonlinear term is negligible
    
    MomTower(1)=0
    MomTower(2)=3*ETow*ICTow/(TurbH)*(Gama(1)-SolidTeta) 
    MomTower(3)=0

    RM(1)=Point_Cyl_1(1,1)-xbar 
    RM(2)=Point_Cyl_1(1,2)-ybar
    RM(3)=Point_Cyl_1(1,3)-zbar
   
   MomTower(:)=MomTower(:)+CROSS(RM,FTower)

 
   if (torder.eq.2.AND.(numimp+1).eq.Iterate) then
   write(2025,2014)  (tp+1)*dt,Gama(1)*180/pi,(Gama(1)-2*Gama(2)+Gama(3))/(dt*dt)*180/pi,SolidTeta*180/pi,asolidx,AnAcy*180/pi, &
                    & 3*ETow*ICTow/(TurbH)*(Gama(1)-SolidTeta),MomTower(2),FTower(1),FTower(3)
                     
   2014  format (10(1x,e15.7))

      !! in correction step  (prediction correction method)  and the last iteration in whcih the value is finilized
      Gama(1)=0.5*(Gama(1)+Gama(3))
      Gama(2)=Gama(3)


      end if

  end subroutine
 




  
subroutine Tower_Constant_Ini()
 implicit none 

include "Par_Constant_Tower.txt"


Di=0.5*(DiB+DiT)
Tick=0.5*(TickT+TickB)
ICTow=pi/64.0*((Di+2*Tick)**4-Di**4)
ITow=1.0/12.0*MTow*(  3*( (Di/2)**2 +((Di+2*Tick)/2)**2 )+TowEH**2  )
IRotY=1.0/4.0*MRot*RRot**2


 end subroutine 

subroutine Tower_Dynamic_Ini()
 implicit none 


OPEN(2025,file='TowerMotion.plt')
write(2025,2015) 'variables="t","Gama","GamaA","Teta","aSolidX","AnAcy","mom","mom Cg","Forcex","Forcez"'  
2015  format (A150)


PlCgH=LENGTHF(point_Cyl_1(1,1)-xbar,point_Cyl_1(1,2)-ybar,point_Cyl_1(1,3)-zbar)  !Towerb+Len-ZBar
!! This Gama intial position is true if we do not wish independent inital tilting for the tower 
Gama(:)=ATAN( (point_Cyl_1(1,1)-point_Cyl_1(2,1))/(point_Cyl_1(1,3)-point_Cyl_1(2,3)) )
print*, "Zbar=",ZBar,"PlCgH=",PlCgH

end subroutine Tower_Dynamic_Ini






Subroutine Tower_Dynamic_Write()

   
   implicit none
   
   OPEN(unit=3025,file=fileplace//"Tower_Dynamic_Save.dat",STATUS='REPLACE') 

522  format (1x,e23.15)

   write(3025,522) gama(1)
   write(3025,522) gama(2)
   write(3025,522) gama(3)
   write(3025,522) PlCgH

   call flush (3025)
   close(3025)

return
end subroutine Tower_Dynamic_Write


Subroutine Tower_Dynamic_Read()

   implicit none

   OPEN(unit=3025,file=fileplace//"Tower_Dynamic_Save.dat")
522  format (1x,e23.15)

   read(3025,522) gama(1)
   read(3025,522) gama(2)
   read(3025,522) gama(3)
   read(3025,522) PlCgH

   close(3025)

   OPEN(2025,file='TowerMotion.plt',POSITION='APPEND',STATUS='OLD')

return
end subroutine Tower_Dynamic_Read



end module  
