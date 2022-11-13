subroutine latgen ()                                           
	use cl
	use bl1
	use bl2
	implicit real*8 (a-h,o-z)

!---------------------------------------------------------------------  
      sqrt3=sqrt(3.d0)
      alpha(1)=alpha(1)*pi/180.0d0                                             
      alpha(2)=alpha(2)*pi/180.0d0                                             
      alpha(3)=alpha(3)*pi/180.0d0                                             
      pia(1)=2.d0*pi/aa                                                 
      pia(2)=2.d0*pi/bb                                                 
      pia(3)=2.d0*pi/cc                                                 
      if(lattic(1:1).eq.'H') goto 10                                    
      if(lattic(1:1).eq.'S') goto 20                                    
      if(lattic(1:1).eq.'P') goto 20                                    
      if(lattic(1:1).eq.'F') goto 30                                    
      if(lattic(1:1).eq.'B') goto 40                                    
      if(lattic(1:1).eq.'C') goto 50                                    
      if(lattic(1:1).eq.'R') goto 60                                    
!
!                                                                       
!.....hexagonal lattice                                                 
 10   continue                                                          
      br1(1,1)=sqrt3/2.d0*aa                                        
      br1(1,2)=-1.d0/2.d0*aa                                        
      br1(1,3)=0.0d0                                                    
      br1(2,1)=0.0d0                                                    
      br1(2,2)=bb                                                   
      br1(2,3)=0.0d0                                                    
      br1(3,1)=0.0d0                                                    
      br1(3,2)=0.0d0                                                    
      br1(3,3)=cc                                                   
                                                                       
      rvfac=2.d0/sqrt(3.d0)                                             
      ortho=.false.                                             
      goto 100                                                          
!                                                                       
!.....rhombohedral case                                                    
 60   br1(1,1)=1.d0/2.d0/sqrt(3.d0)*aa                                   
      br1(1,2)=1.d0/2.d0/sqrt(3.d0)*aa                                   
      br1(1,3)=-1.d0/sqrt(3.d0)*aa                                         
      br1(2,1)=-1.0d0/2.d0*bb                                                  
      br1(2,2)=1.0d0/2.d0*bb                                                    
      br1(2,3)=0.0d0*bb                                                    
      br1(3,1)=1.0d0/3.0d0*cc                                                 
      br1(3,2)=1.0d0/3.0d0*cc                                               
      br1(3,3)=1.0d0/3.0d0*cc                                                 
                                                   
      rvfac=6.d0/sqrt(3.d0)
      ortho=.false.                                             
      goto 100                                                          
!                                                                       
!.....primitive lattice                                                 
!                                                                       
  20  continue

      sinbc=sin(alpha(1))
      cosab=cos(alpha(3))
      cosac=cos(alpha(2))
      cosbc=cos(alpha(1))
      br1(1,1)= sqrt(1-((cosab-cosac*cosbc)/sinbc)**2-cosac**2)*aa
      br1(1,2)= aa*(cosab-cosac*cosbc)/sinbc
      br1(1,3)= aa*cosac
      br1(2,1)= 0.0
      br1(2,2)= bb*sinbc
      br1(2,3)= bb*cosbc
      br1(3,1)= 0.0
      br1(3,2)= 0.0
      br1(3,3)= cc

      rvfac= 1.d0/wurzel
      ortho=.true.
      if(abs(alpha(1)-pi/2.d0).gt.0.0001) ortho=.false.
      if(abs(alpha(2)-pi/2.d0).gt.0.0001) ortho=.false.
      if(abs(alpha(3)-pi/2.d0).gt.0.0001) ortho=.false.
!
      goto 100
!                                                                       
!.....fc lattice                                                        
 30   continue                                                          
      br1(1,1)=aa                                                   
      br1(1,2)=0.0d0                                                    
      br1(1,3)=0.0d0                                                    
      br1(2,1)=0.0d0                                                    
      br1(2,2)=bb                                                   
      br1(2,3)=0.0d0                                                    
      br1(3,2)=0.0d0                                                    
      br1(3,1)=0.0d0                                                    
      br1(3,3)=cc                                                   
                                                                    
      rvfac=4.d0                                                        
      ortho=.true.                                             
      goto 100                                                          
!                                                                       
!.....bc lattice                                                        
 40   continue                                                          
      br1(1,1)=aa                                                   
      br1(1,2)=0.0d0                                                    
      br1(1,3)=0.0d0                                                    
      br1(2,1)=0.0d0                                                    
      br1(2,2)=bb                                                   
      br1(2,3)=0.0d0                                                    
      br1(3,1)=0.0d0                                                    
      br1(3,2)=0.0d0                                                    
      br1(3,3)=cc                                                   

      rvfac=2.d0
      ortho=.true.                                             
      goto 100                                                    
!             
 50   continue                                                          
      if(lattic(2:3).eq.'XZ') goto 51                                    
      if(lattic(2:3).eq.'YZ') goto 52                                    
!.....cxy lattice                                                          
      br1(1,1)=aa                                                   
      br1(1,2)=0.0d0                                                    
      br1(1,3)=0.0d0                                                    
      br1(2,1)=0.0d0                                                    
      br1(2,2)=bb                                                   
      br1(2,3)=0.0d0                                                    
      br1(3,1)=0.0d0                                                    
      br1(3,2)=0.0d0                                                    
      br1(3,3)=cc                                                   
                                                                      
      rvfac=2.d0                                                        
      ortho=.true.                                             
      goto 100                                                          
!                                                                       
!.....cxz case (cxz lattice build up)                                     
 51   continue                                     
!.....cxz orthorombic case 
      if(abs(alpha(3)-pi/2.0d0).lt.0.0001) then
         br1(1,1)=aa                                                   
         br1(1,2)=0.0d0                                                    
         br1(1,3)=0.0d0                                                    
         br1(2,1)=0.0d0                                                    
         br1(2,2)=bb                                                   
         br1(2,3)=0.0d0                                                    
         br1(3,1)=0.0d0                                                    
         br1(3,2)=0.0d0                                                    
         br1(3,3)=cc                                                   
                                                                      
         rvfac=2.0                                                         
         ortho=.true.                                             
         goto 100                                                          
      else
!.....cxz monoclinic case 
!         write(*,*) '  gamma not equal 90'
         sinab=sin(alpha(3))
         cosab=cos(alpha(3))
!                                                                       
         br1(1,1)= aa*sinab 
         br1(1,2)= aa*cosab
         br1(1,3)= 0.0                                                   
         br1(2,1)= 0.0                                                      
         br1(2,2)= bb                                                     
         br1(2,3)= 0.0                                                     
         br1(3,1)= 0.0                                                     
         br1(3,2)= 0.0                                                     
         br1(3,3)= cc                                                     
                                                                      
         rvfac=2.0/sinab                                                   
         ortho=.false.                                             
         goto 100                                                          
      endif
!                                                                       
!.....cyz case (cyz lattice build up)                                     
 52   continue                                     
      br1(1,1)=aa                                                   
      br1(1,2)=0.0d0                                                    
      br1(1,3)=0.0d0                                                    
      br1(2,1)=0.0d0                                                    
      br1(2,2)=bb                                                   
      br1(2,3)=0.0d0                                                    
      br1(3,1)=0.0d0                                                    
      br1(3,2)=0.0d0                                                    
      br1(3,3)=cc                                                   
                                                                       
      rvfac=2.0                                                         
      ortho=.true.                                             
      goto 100                                                          
!                                                                       
!.....define volume of unit cell                                        
 100  continue                                                          
      vol=aa*bb*cc/rvfac                                                
!
      return
      end                                                               

